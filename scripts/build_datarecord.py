# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Exports the network handed to the solver as a datarecord.

Reads the composed network of one planning horizon, applies the preparation
that `solve_network` applies before it builds the model, and writes the
network as a standalone datarecord with `datarecord.tools.pypsa.PyPSA`.
The shapes are set here only, since not every composed network carries
them: the clustered onshore and offshore regions on the buses, and the
country and offshore shapes on a `country` dim that `bus.country` refers to
(`set_shapes`).
Custom constraints of `solve_network` are not exported: a solve of the record
uses PyPSA's own math.

The script needs the `datarecord` pixi environment:
`pixi run -e datarecord snakemake -call build_datarecord`.
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
from datarecord import Workspace
from datarecord.tools.pypsa import PyPSA

from scripts._helpers import configure_logging, set_scenario_config
from scripts.solve_network import prepare_network

logger = logging.getLogger(__name__)

STUDY_REGION = "EU"


def set_shapes(n: pypsa.Network, inputs) -> None:
    """
    Replace the network's shapes with the bus regions and the country shapes.

    Bus regions are `Bus` shapes of type `onshore` and `offshore`. Country and
    offshore shapes are shapes of the component `country`, which the record
    holds as a `country` dim, with the same two types. Where a bus names the
    study region `EU`, it gets the dissolve of every country's shape of each
    type. A `bus.country` naming no country shape is dropped, since the record
    refers `bus.country` to the `country` dim.
    """
    n.remove("Shape", n.shapes.index)
    sources = {
        ("Bus", "onshore"): inputs.onshore_regions,
        ("Bus", "offshore"): inputs.offshore_regions,
        ("country", "onshore"): inputs.country_shapes,
        ("country", "offshore"): inputs.offshore_shapes,
    }
    for (component, which), path in sources.items():
        shapes = gpd.read_file(path).to_crs(n.crs)
        if component == "country" and (n.buses.country == STUDY_REGION).any():
            region = gpd.GeoDataFrame(
                {"name": [STUDY_REGION]}, geometry=[shapes.union_all()], crs=n.crs
            )
            shapes = pd.concat([shapes, region], ignore_index=True)
        n.add(
            "Shape",
            f"{component} " + shapes.name + f" {which}",
            geometry=shapes.geometry.values,
            idx=shapes.name.values,
            component=component,
            type=which,
        )
    named = n.c.Shape.static.query("component == 'country'").idx
    unknown = n.buses.country.notna() & ~n.buses.country.isin(named)
    if unknown.any():
        logger.warning(
            f"Dropping bus countries without a country shape: "
            f"{sorted(n.buses.country[unknown].unique())}."
        )
    n.buses.loc[unknown, "country"] = np.nan


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_datarecord",
            configfiles="config/test/config.preaggregated.yaml",
            horizon=2050,
        )
    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)

    solve_opts = snakemake.params.solving["options"]
    np.random.seed(solve_opts.get("seed", 123))

    n = pypsa.Network(snakemake.input.network)
    prepare_network(
        n,
        config=snakemake.config,
        solve_opts=solve_opts,
        foresight=snakemake.params.foresight,
        planning_horizons=snakemake.wildcards.horizon,
        co2_sequestration_potential=snakemake.params.co2_sequestration_potential,
        limit_max_growth=snakemake.params.sector["limit_max_growth"],
        rolling_horizon=solve_opts["rolling_horizon"],
    )

    for c in n.components[["Bus", "Line", "Link"]]:
        flags = c.static.columns.intersection(["under_construction", "underground"])
        c.static[flags] = c.static[flags].fillna(False).astype(bool)

    set_shapes(n, snakemake.input)

    ws = Workspace()
    ws.write(PyPSA.to_datarecord(n, ws.con), snakemake.output.directory)
    logger.info(f"Wrote the datarecord to {snakemake.output.directory}.")
