# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT

"""
Baseline electricty demand data processor.

This script processes baseline electricity demand data from the FES workbook.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def parse_demand_data(
    regional_gb_data_path: str,
) -> pd.DataFrame:
    """
    Parse the regional gb data to obtain baseline electricity demand in the required format

    Args:
        regional_gb_data_path(str): Filepath to the regional data file containing the electricity baseline demand

    Returns:
        pd.DataFrame : MultiIndex dataframe containing the baseline electricity demand indexed by PyPSA bus and year
    """
    df_regional_data = pd.read_csv(regional_gb_data_path)

    # Filter the data for baseline demand
    df_demand = df_regional_data.query("`Technology Detail` == 'Baseline Demand'")

    # Convert from GWh to MWh
    df_demand.loc[:, "data"] *= 1000

    df_demand_grouped = df_demand.groupby(["bus", "year"])["data"].sum()

    df_demand_grouped.name = "p_set"

    return df_demand_grouped


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the regional gb data file path
    regional_gb_data_path = snakemake.input.regional_gb_data

    # Parse input data
    df_baseline_electricity_demand = parse_demand_data(regional_gb_data_path)

    # Write demand dataframe to csv file
    df_baseline_electricity_demand.to_csv(
        snakemake.output["baseline_electricity_demand"]
    )
    logger.info(
        f"Baseline electricity demand data saved to {snakemake.output['baseline_electricity_demand']}"
    )
