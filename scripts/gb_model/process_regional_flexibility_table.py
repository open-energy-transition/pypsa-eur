# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Regional Flexibility data processor.

This script splits flexibility data into regionsfrom the FES workbook.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import get_regional_distribution

logger = logging.getLogger(__name__)


def parse_regional_flexibility_data(
    flexibility_data_path: str,
    regional_gb_data_path: str,
    regional_distribution_reference: list,
) -> pd.DataFrame:
    """
    Parse the flexibility data to obtain required flexibility capacity in the required format.

    Args:
        flexibility_data_path (str): Filepath to the flexibility data CSV file containing
                                   flexibility capacity data by technology and year
        technology_detail (list): List of technology details relevant to a particular
                                flexibility type (e.g., EV DSM technologies)
        fes_scenario (str): FES scenario name to filter by (e.g., "leading the way")
        year_range (list): Two-element list [start_year, end_year] defining the year range to include

    Returns:
        pd.Series: Series containing aggregated flexibility capacity indexed by year.
                  Values represent total flexibility capacity in MW for the specified
                  technology types and scenario.

    Processing steps:
        1. Load and pre-format flexibility data from CSV file
        2. Filter by technology type, scenario, and year range
        3. Convert units from GW to MW and aggregate by year
    """
    # Load annual flexibility data
    flexibility_data = pd.read_csv(flexibility_data_path, index_col=0)

    # Load regional GB data
    regional_gb_data = pd.read_csv(regional_gb_data_path)

    # Obtain regional reference for distribution
    regional_reference = regional_gb_data[
        regional_gb_data["Technology Detail"].str.lower().isin(regional_distribution_reference)
    ]

    # Group by bus and year
    regional_reference = regional_reference.groupby(["bus", "year"])["data"].sum()

    # Get regional distribution
    regional_dist = get_regional_distribution(regional_reference)

    # Fill backward for NaN values (there are 0 instances in early years for regional data, but not in annual data)
    regional_dist = regional_dist.bfill()

    # Distribute flexibility data regionally
    flexibility_regional = regional_dist * flexibility_data['p_nom']

    # Set name as p_nom
    flexibility_regional.name = "p_nom"

    return flexibility_regional


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem, flexibility_type="fes_ev_dsm")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the inputs
    flexibility_data_path = snakemake.input.flexibility
    regional_gb_data_path = snakemake.input.regional_gb_data

    # Parse input data
    flexibility_type = snakemake.params.flexibility_type
    regional_distribution_reference = snakemake.params.regional_distribution_reference[flexibility_type]

    df_regional_flexibility = parse_regional_flexibility_data(
        flexibility_data_path,
        regional_gb_data_path,
        regional_distribution_reference,
    )

    # Write regional flexibility dataframe to csv file
    df_regional_flexibility.to_csv(snakemake.output.regional_flexibility)
    logger.info(f"Regional flexibility data saved to {snakemake.output.regional_flexibility}")
