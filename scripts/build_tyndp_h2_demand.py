# SPDX-FileCopyrightText: Contributors to Open-TYNDP <https://github.com/open-energy-transition/open-tyndp>
#
# SPDX-License-Identifier: MIT

import logging
import multiprocessing as mp
from functools import partial
from pathlib import Path

import pandas as pd
from _helpers import configure_logging, get_snapshots, set_scenario_config
from tqdm import tqdm

logger = logging.getLogger(__name__)


def extract_from_sheet(
    h2_demand_fn: str, sheet_name: str, weather_year: int
) -> pd.Series:
    """
    Extract weather year data from a single sheet in the H2 demand Excel file.

    Parameters
    ----------
    h2_demand_fn : str
        Path to the TYNDP hydrogen demand Excel file.
    sheet_name : str
        Name of the sheet to read.
    weather_year : int
        The weather year to extract.

    Returns
    -------
    pd.Series
        Series with Date-Hour multi-index containing the weather year data.
        Returns empty Series if weather year not found or if reading fails.
    """
    try:
        # Read sheet starting from row 10 (0-indexed), so header is row 11
        df_sheet = pd.read_excel(
            h2_demand_fn,
            sheet_name=sheet_name,
            engine="openpyxl",
            header=10  # Row 11 becomes header (0-indexed)
        )

        # Set first two columns as index (Date and Hour)
        df_sheet = df_sheet.set_index([df_sheet.columns[0], df_sheet.columns[1]])

        # Find column with specified weather year (columns are integers)
        if weather_year in df_sheet.columns:
            # Extract weather year data as a Series
            series_weather_year = df_sheet[weather_year]
            logger.info(
                f"Sheet '{sheet_name}' - Extracted {weather_year} data, length: {len(series_weather_year)}"
            )
            return series_weather_year
        else:
            logger.warning(
                f"Sheet '{sheet_name}' - No {weather_year} weather year column found. "
                f"Available years: {df_sheet.columns.tolist()}"
            )
            return pd.Series()

    except Exception as e:
        logger.warning(
            f"Sheet '{sheet_name}' - Failed to extract data: {type(e).__name__}: {e}"
        )
        return pd.Series()


def read_h2_demand_zone(
    h2_demand_fn: str, zone_name: str, weather_year: int, snapshots: pd.DatetimeIndex
) -> pd.DataFrame:
    """
    Read hydrogen demand data from TYNDP Excel file for a single zone.

    Each Excel file contains multiple sheets (one per node). Each sheet has:
    - Row 11 (index 10) as header with weather years as column names
    - First two columns: Date and Hour (used as multi-index)
    - Remaining columns: different weather years (1982-2019)

    This function extracts the specified weather year data from each sheet and
    combines them into a single DataFrame, reindexed to match the snapshots.

    Parameters
    ----------
    h2_demand_fn : str
        Path to the TYNDP hydrogen demand Excel file.
    zone_name : str
        Name of the zone (e.g., "ZONE_1", "ZONE_2").
    weather_year : int
        The weather year to extract from the data.
    snapshots : pd.DatetimeIndex
        Target snapshots to use as index.

    Returns
    -------
    pd.DataFrame
        DataFrame with snapshots as index and columns named by sheet names,
        containing the specified weather year data.
    """
    try:
        # Read all sheet names
        excel_file = pd.ExcelFile(h2_demand_fn, engine="openpyxl")
        sheet_names = excel_file.sheet_names

        logger.info(f"Read H2 demand data for {zone_name} from {h2_demand_fn}")
        logger.info(f"Found {len(sheet_names)} sheets: {sheet_names}")
        logger.info(f"Extracting weather year: {weather_year}")

        # Dictionary to store series from each sheet
        sheet_series = {}

        for sheet_name in sheet_names:
            # Extract weather year data from this sheet
            series_weather_year = extract_from_sheet(h2_demand_fn, sheet_name, weather_year)

            # Add to dictionary if data was successfully extracted
            if not series_weather_year.empty:
                sheet_series[sheet_name] = series_weather_year

        # Combine all series into a single DataFrame
        if sheet_series:
            df_combined = pd.DataFrame(sheet_series)

            # Reindex to match snapshots
            df_combined.index = df_combined.index.map(
                lambda t: t.replace(year=snapshots[0].year)
            )
            df_combined = df_combined.reindex(snapshots)

            logger.info(f"{zone_name} - Combined shape: {df_combined.shape}")
            logger.info(f"{zone_name} - Columns: {df_combined.columns.tolist()}")
            return df_combined
        else:
            logger.warning(f"{zone_name} - No data extracted, returning empty DataFrame with snapshots index")
            return pd.DataFrame(index=snapshots)

    except Exception as e:
        logger.warning(
            f"Failed to read {zone_name} from {h2_demand_fn}: {type(e).__name__}: {e}"
        )
        logger.warning(f"Creating empty DataFrame for {zone_name} with snapshots index")
        return pd.DataFrame(index=snapshots)


def load_h2_demand_zone(
    zone_number: int,
    base_path: str,
    scenario: str,
    pyear: int,
    cyear: int,
    snapshots: pd.DatetimeIndex,
) -> pd.DataFrame:
    """
    Load H2 demand data for a single zone file.

    This function is designed to be called in parallel for multiple zones.

    Parameters
    ----------
    zone_number : int
        Zone number (1 or 2).
    base_path : str
        Base path to the demand profiles directory.
    scenario : str
        Scenario name (e.g., "DE", "GA", "NT").
    pyear : int
        Planning year.
    cyear : int
        Climatic/weather year.
    snapshots : pd.DatetimeIndex
        Target snapshots to use as index.

    Returns
    -------
    pd.DataFrame
        DataFrame with snapshots as index and columns named by sheet names.
    """
    zone_name = f"ZONE_{zone_number}"
    h2_demand_path = f"{base_path}/{scenario}/{pyear}/H2_{zone_name}.xlsx"

    logger.info(f"Processing {zone_name} from {h2_demand_path}")

    return read_h2_demand_zone(h2_demand_path, zone_name, cyear, snapshots)


def combine_h2_zones(df_z1: pd.DataFrame, df_z2: pd.DataFrame) -> pd.DataFrame:
    """
    Combine hydrogen demand data from multiple zones.

    Each zone DataFrame has nodes as columns. This function combines them
    horizontally, keeping the Date-Hour multi-index.

    Parameters
    ----------
    df_z1 : pd.DataFrame
        Hydrogen demand data for Zone 1 (columns are node names).
    df_z2 : pd.DataFrame
        Hydrogen demand data for Zone 2 (columns are node names).

    Returns
    -------
    pd.DataFrame
        Combined hydrogen demand data with all nodes from both zones.
    """
    # Combine both zones horizontally (add columns from both)
    df_combined = pd.concat([df_z1, df_z2], axis=1)

    logger.info(f"Combined H2 demand data - Shape: {df_combined.shape}")
    logger.info(f"Combined H2 demand data - Columns: {df_combined.columns.tolist()}")

    return df_combined


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_tyndp_h2_demand", planning_horizons="2030")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = snakemake.params["scenario"]
    pyear = int(snakemake.wildcards.planning_horizons)
    cyear = get_snapshots(snakemake.params.snapshots)[0].year
    base_path = snakemake.input.h2_demand

    logger.info(
        f"Processing H2 demand for scenario: {scenario}, planning year: {pyear}, weather year: {cyear}"
    )

    # Load H2 demand for both zones in parallel
    zones = [1, 2]
    tqdm_kwargs = {
        "ascii": False,
        "unit": " zone",
        "total": len(zones),
        "desc": "Loading TYNDP H2 demand data",
    }

    func = partial(
        load_h2_demand_zone,
        base_path=base_path,
        scenario=scenario,
        pyear=pyear,
        cyear=cyear,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        h2_demand_zones = list(tqdm(pool.imap(func, zones), **tqdm_kwargs))

    # Combine zones horizontally (concatenate columns)
    h2_demand = pd.concat(h2_demand_zones, axis=1)

    logger.info(f"Combined H2 demand data - Shape: {h2_demand.shape}")
    logger.info(f"Combined H2 demand data - Columns: {h2_demand.columns.tolist()}")

    # Export to CSV
    h2_demand.to_csv(snakemake.output.h2_demand, index=True)

    logger.info(f"H2 demand data exported to {snakemake.output.h2_demand}")
