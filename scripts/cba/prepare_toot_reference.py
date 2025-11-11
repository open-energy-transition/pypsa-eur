# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

"""
Prepare TOOT reference network by ensuring all CBA projects are included.

For TOOT methodology, the reference network should contain ALL CBA projects.
The simple_2030.nc network already contains projects where in_reference2030=True.
This script adds the missing projects (where in_reference2030=False) to create
the complete TOOT reference network with all projects included.

Handles reference year to planning horizon mapping:
- reference2030 → horizon 2030
- reference2035 → horizon 2040 (if 2035 not available)
"""

import logging

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_toot_reference")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load simplified network
    n = pypsa.Network(snakemake.input.network)

    # Load project definitions
    transmission_projects = pd.read_csv(snakemake.input.transmission_projects)

    # Get planning horizons from config
    planning_horizons = int(snakemake.wildcards.planning_horizons)

    logger.info(f"\n{'=' * 80}")
    logger.info("PREPARING TOOT REFERENCE NETWORK")
    logger.info(f"{'=' * 80}")
    logger.info(f"Current planning horizon: {planning_horizons}")

    # For TOOT: Add projects that are NOT in the reference (in_reference=False)
    # the reference network should have ALL projects
    projects_added = []
    projects_created = []
    projects_in_base = []
    projects_missing_buses = []

    # Add transmission projects for this horizon
    for _, project in transmission_projects.iterrows():
        in_reference: bool = project[f"in_reference{planning_horizons}"]

        bus0 = project["bus0"]
        bus1 = project["bus1"]
        link_id = f"{bus0}-{bus1}-DC"
        reverse_link_id = f"{bus1}-{bus0}-DC"

        if not in_reference:
            # Project is NOT in the base network - ADD it to create TOOT reference
            capacity_0to1 = project["p_nom 0->1"]
            capacity_1to0 = project["p_nom 1->0"]

            links_exist = link_id in n.links.index and reverse_link_id in n.links.index

            if links_exist:
                # Links exist - add capacity to existing links
                original_capacity = n.links.loc[link_id, "p_nom"]
                original_capacity_reverse = n.links.loc[reverse_link_id, "p_nom"]

                n.links.loc[link_id, "p_nom"] += capacity_0to1
                n.links.loc[reverse_link_id, "p_nom"] += capacity_1to0

                projects_added.append(project["project_id"])
                logger.info(
                    f"Project {project['project_id']} ({project['project_name']}) ADDED to existing links:"
                )
                logger.info(
                    f"    Link {link_id}: {original_capacity:.0f} → {n.links.loc[link_id, 'p_nom']:.0f} MW (+{capacity_0to1:.0f} MW)"
                )
                logger.info(
                    f"    Link {reverse_link_id}: {original_capacity_reverse:.0f} → {n.links.loc[reverse_link_id, 'p_nom']:.0f} MW (+{capacity_1to0:.0f} MW)"
                )
            else:
                # Links don't exist but buses do - CREATE new links
                # Get a template DC link for attributes
                dc_links = n.links[n.links.carrier == "DC"]
                if len(dc_links) == 0:
                    logger.error(
                        "Cannot create new links: No DC links found in network as template"
                    )
                    projects_missing_buses.append(project["project_id"])
                    continue

                template = dc_links.iloc[0]

                # Create forward link
                n.add(
                    "Link",
                    link_id,
                    bus0=bus0,
                    bus1=bus1,
                    carrier="DC",
                    p_nom=capacity_0to1,
                    p_nom_extendable=False,
                    length=0,  # Will be calculated if needed
                    marginal_cost=template.marginal_cost,  # TODO: adjust if needed
                    capital_cost=template.capital_cost,  # Not relevant in market simulation, will be set in CBA B1
                )

                # Create reverse link
                n.add(
                    "Link",
                    reverse_link_id,
                    bus0=bus1,
                    bus1=bus0,
                    carrier="DC",
                    p_nom=capacity_1to0,
                    p_nom_extendable=False,
                    length=0,
                    marginal_cost=template.marginal_cost,
                    capital_cost=template.capital_cost,
                )

                projects_created.append(project["project_id"])
                logger.info(
                    f"Project {project['project_id']} ({project['project_name']}) CREATED new links:"
                )
                logger.info(f"    Link {link_id}: 0 → {capacity_0to1:.0f} MW (new)")
                logger.info(
                    f"    Link {reverse_link_id}: 0 → {capacity_1to0:.0f} MW (new)"
                )

        else:
            # Project IS already in the base network - log it for transparency
            if link_id in n.links.index and reverse_link_id in n.links.index:
                capacity_forward = n.links.loc[link_id, "p_nom"]
                capacity_reverse = n.links.loc[reverse_link_id, "p_nom"]

                projects_in_base.append(project["project_id"])
                logger.info(
                    f"Project {project['project_id']} ({project['project_name']}) already in base network:"
                )
                logger.info(
                    f"    Link {link_id}: {capacity_forward:.0f} MW (no change)"
                )
                logger.info(
                    f"    Link {reverse_link_id}: {capacity_reverse:.0f} MW (no change)"
                )
            else:
                projects_missing_buses.append(project["project_id"])
                logger.warning(
                    f" Project {project['project_id']} ({project['project_name']}): Marked in_reference but links not found"
                )
                logger.warning(f"    Missing: {link_id} and/or {reverse_link_id}")

    # Analyze border aggregation (multiple projects per border)
    border_projects = {}
    new_links_list = []
    for _, project in transmission_projects.iterrows():
        in_reference: bool = project[f"in_reference{planning_horizons}"]
        bus0, bus1 = project["bus0"], project["bus1"]

        if not in_reference:
            border = f"{bus0}-{bus1}"
            if border not in border_projects:
                border_projects[border] = []
            border_projects[border].append(project["project_id"])

            # Track new links with project details
            if project["project_id"] in projects_created:
                new_links_list.append(
                    {
                        "border": border,
                        "id": project["project_id"],
                        "name": project["project_name"],
                    }
                )

    borders_with_multiple = {
        b: pids for b, pids in border_projects.items() if len(pids) > 1
    }

    # Summary
    logger.info(f"\n{'=' * 80}")
    logger.info("TOOT REFERENCE PREPARATION SUMMARY")
    logger.info(f"{'=' * 80}")

    if projects_added:
        logger.info(f" Added capacity to {len(projects_added)} existing links")

    if projects_created:
        logger.info(f" Created {len(projects_created)} new links")
        for link_info in new_links_list:
            logger.info(
                f"    {link_info['border']}: project {link_info['id']} ({link_info['name']})"
            )

    if projects_in_base:
        logger.info(f" {len(projects_in_base)} projects already in base network")

    total_in_reference = (
        len(projects_added) + len(projects_created) + len(projects_in_base)
    )
    logger.info(f" Total projects in TOOT reference: {total_in_reference}")

    if borders_with_multiple:
        logger.info(
            f" {len(borders_with_multiple)} borders with multiple projects (capacity aggregated):"
        )
        for border, pids in sorted(borders_with_multiple.items()):
            logger.info(f"    {border}: projects {pids}")

    if projects_missing_buses:
        logger.info(
            f" {len(projects_missing_buses)} projects excluded (buses not in network): {projects_missing_buses}"
        )

    logger.info(f"{'=' * 80}\n")

    # Save reference network with all projects
    n.export_to_netcdf(snakemake.output.network)
    logger.info(
        f"TOOT reference network saved (horizon {planning_horizons}, {len(projects_added)} capacity added, {len(projects_created)} links created)"
    )
