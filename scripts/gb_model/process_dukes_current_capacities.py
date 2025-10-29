# SPDX-FileCopyrightText:  gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT


"""
GSP-level data table generator.

This is a script to combine the BB1 sheet with the BB2 (metadata) sheet of the FES workbook.
"""

import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import map_points_to_regions

logger = logging.getLogger(__name__)


def _add_region_cols(df: pd.DataFrame, region_df: pd.DataFrame) -> pd.DataFrame:
    cols_to_add = region_df[["name", "TO_region"]].rename(columns={"name": "bus"})
    df = df.assign(
        **{col: float("nan") for col in cols_to_add.columns if col not in df.columns}
    )
    return df.fillna(cols_to_add)


def process_dukes_data(
    df_dukes: pd.DataFrame,
    gdf_regions: gpd.GeoDataFrame,
) -> pd.DataFrame:
    region_data = map_points_to_regions(
        df_dukes,
        gdf_regions,
        "Y-Coordinate",
        "X-Coordinate",
        "EPSG:27700",
        snakemake.params.target_crs,
    )
    df_dukes = _add_region_cols(df_dukes, region_data)

    initial_drop = df_dukes["X-Coordinate"].isnull().sum()
    bus_drop = (df_dukes["X-Coordinate"].notnull() & df_dukes["bus"].isnull()).sum()

    logger.info(
        f"DUKES capacity table | Dropping {initial_drop} rows with missing coordinate data plus {bus_drop} rows that could not be mapped to a bus."
    )

    return df_dukes.dropna(subset=["bus"])


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the file paths
    gdf_regions = gpd.read_file(snakemake.input.regions)
    gdf_regions_offshore = gpd.read_file(snakemake.input.regions_offshore)
    gdf_regions_all = pd.concat([gdf_regions, gdf_regions_offshore]).query(
        "country == 'GB'"
    )
    sheet_config = snakemake.params.sheet_config
    sheet_name = sheet_config.pop("sheet_name")
    df_dukes = pd.read_excel(
        snakemake.input.dukes_data, sheet_name=sheet_name, **sheet_config
    )
    df_dukes_cleaned = process_dukes_data(df_dukes, gdf_regions_all)

    df_dukes_cleaned.to_csv(snakemake.output.csv, index=False)
