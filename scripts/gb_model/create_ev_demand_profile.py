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

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import get_regional_distribution

logger = logging.getLogger(__name__)


def prepare_ev_demand_profiles(
    ev_demand_path: str,
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
    # Read regional GB data
    regional_gb_data = pd.read_csv(regional_gb_data_path)

    # Select regional hydrogen electrolysis data
    electrolysis_data = regional_gb_data[
        (regional_gb_data["Technology Detail"].str.lower() == "hydrogen electrolysis")
    ]

    # Calculate regional grid-connected electrolysis capacities
    regional_grid_electrolysis_capacities = electrolysis_data.groupby(["bus", "year"])[
        "data"
    ].sum()

    # Calculate regional distribution of grid-connected electrolysis capacities
    electrolysis_distribution = get_regional_distribution(
        regional_grid_electrolysis_capacities
    )

    # Map unmapped grid-connected electrolysis data
    unmapped_grid_electrolysis_capacities = (
        electrolysis_data[electrolysis_data["bus"].isnull()]
        .groupby("year")["data"]
        .sum()
    )
    unmapped_grid_electrolysis_capacities = (
        unmapped_grid_electrolysis_capacities * electrolysis_distribution
    )

    # Combine mapped and unmapped grid-connected electrolysis data
    grid_electrolysis_capacities = (
        regional_grid_electrolysis_capacities + unmapped_grid_electrolysis_capacities
    )

    # Rename series to 'p_nom'
    grid_electrolysis_capacities.name = "p_nom"

    return grid_electrolysis_capacities


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the file path
    ev_demand_path = snakemake.input.ev_demand

    # Prepare EV demand profiles
    ev_demand_profiles = prepare_ev_demand_profiles(ev_demand_path)

    # Save the EV demand profiles
    ev_demand_profiles.to_csv(snakemake.output.ev_demand_profiles)
    logger.info(
        f"Grid-connected electrolysis capacities saved to {snakemake.output.grid_electrolysis_capacities}"
    )
