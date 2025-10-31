# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT

"""
EV demand profile processor.

This script prepares regional EV demand profiles.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import (
    configure_logging,
    generate_periodic_profiles,
    get_snapshots,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def prepare_ev_demand_profiles(
    ev_demand_path: str,
    traffic_fn: str,
    snapshots: pd.DatetimeIndex,
) -> pd.DataFrame:
    """
    Parse and prepare EV demand profiles.

    This function processes regional hydrogen electrolysis capacity data, calculates
    regional distributions, and redistributes unmapped capacities based on existing
    regional patterns.

    Args:
        regional_gb_data_path (str): Path to the regional GB data CSV file containing
                                   hydrogen electrolysis capacity data by region and year.

    Returns:
        pd.DataFrame: Processed grid-connected electrolysis capacities with MultiIndex
                     ['bus', 'year'] and 'data' values in original units in MW.
                     Includes both originally mapped capacities and redistributed
                     unmapped capacities.

    Processing steps:
        1. Filter regional data for hydrogen electrolysis technology
        2. Group mapped data by bus (region) and year
        3. Calculate regional distribution proportions for each year
        4. Identify and aggregate unmapped capacities (NaN bus values)
        5. Redistribute unmapped capacities based on regional distribution patterns
        6. Combine mapped and redistributed capacities into final dataset
    """
    # Read regional EV demand data
    regional_ev_demand = pd.read_csv(ev_demand_path)

    # Read averaged weekly traffic counts from the year 2010-2015
    traffic = pd.read_csv(traffic_fn, skiprows=2, usecols=["count"]).squeeze("columns")

    # Determine nodes (regions) from the regional EV demand data
    nodes = regional_ev_demand["bus"].unique()

    # Create annual profile take account time zone + summer time
    transport_shape = generate_periodic_profiles(
        dt_index=snapshots,
        nodes=nodes,
        weekly_profile=traffic.values,
    )
    transport_shape = transport_shape / transport_shape.sum()

    return transport_shape


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the file path
    ev_demand_path = snakemake.input.ev_demand

    # Define snapshots and nodes
    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day, tz="UTC"
    )

    # Prepare EV demand profile shape
    ev_demand_profiles_shape = prepare_ev_demand_profiles(
        ev_demand_path,
        traffic_fn=snakemake.input.traffic_data_KFZ,
        snapshots=snapshots,
    )

    # Save the EV demand profiles
    ev_demand_profiles_shape.to_csv(snakemake.output.ev_demand_profile_shape)
    logger.info(
        f"Grid-connected electrolysis capacities saved to {snakemake.output.grid_electrolysis_capacities}"
    )
