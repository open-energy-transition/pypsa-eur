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


def multiindex_to_datetimeindex(df, year):
    """
    Convert MultiIndex with ('Date', 'Hour') to DatetimeIndex.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with MultiIndex of ('Date', 'Hour')
    year : int
        Year to use for the DatetimeIndex
        
    Returns
    -------
    pd.DataFrame
        DataFrame with DatetimeIndex
    """
    # Reset index to work with the data
    df_reset = df.reset_index()
    
    # Parse the date string and combine with year
    # Date format is 'DD.MM.' and Hour is 1-24
    df_reset['datetime'] = pd.to_datetime(
        df_reset['Date'].str.strip('.') + f'.{year} ' + 
        (df_reset['Hour'] - 1).astype(str) + ':00',
        format='%d.%m.%Y %H:%M'
    )
    
    # Set as index and drop the old columns
    df_new = df_reset.set_index('datetime').drop(columns=['Date', 'Hour'])
    
    return df_new
    
    
def read_h2_excel(demand_fn: str, scenario: str, pyear: int, cyear: int, h2zone: int) -> pd.DataFrame:

    try:
        data = pd.read_excel(
                demand_fn,
                header=10,
                index_col=[0, 1],
                sheet_name=None,
                usecols=lambda name: name == "Date" or name == "Hour" or name == int(cyear),
                )
        
        
        demand = pd.concat(data, axis=1).droplevel(1, axis=1)
        # Reindex to match snapshots
        demand = multiindex_to_datetimeindex(demand, year=cyear)
        # rename UK in GB
        demand.columns = demand.columns.str.replace("UK", "GB")
        
        demand.columns = [f"{col[:2]} H2 Z{h2zone}" for col in demand.columns]
        demand.columns.name = "Bus"
        
    except Exception as e:
        logger.warning(
            f"Failed to read H2 demand for scenario {scenario} and pyear {pyear}: {type(e).__name__}: {e}"
        )
        demand = pd.DataFrame()

    
    return demand


def load_h2_demand(
    fn: str, scenario: str, pyear: int, cyear: int
):
    """
    Load hydrogen demand files into dictionary of dataframes. Filter for specific climatic year and format data.
    """
    
    # NT scenario H2 demand is not split into 2 zones
    if scenario == "NT":
        if pyear == 2050:
            logger.warning(
                "2050 hydrogen demand data are not defined for NT in 2024 TYNDP cycle. Falling back to 2040."
            )
            pyear = 2040
            
        demand_fn = Path(
            fn,
            scenario,
            "H2 demand profiles",
            f"H2 {pyear}",
            f"{scenario}_{pyear}.xlsx",
        )

        # check if climate year is available
        if int(cyear) < 1982 or int(cyear) > 2019:
            logger.warning(
                "Snapshot year doesn't match available TYNDP data. Falling back to 2009."
            )
            cyear = 2009
            
        demand = read_h2_excel(demand_fn, scenario, pyear, cyear, h2zone=2)
    # DE and GA scenario H2 demand is split into 2 zones
    if scenario in ["DE", "GA"]:
        
        if int(cyear) not in [1995, 2008, 2009]:
            logger.warning(
                "Snapshot year doesn't match available TYNDP data. Falling back to 2009."
            )
            cyear = 2009
        
        demand = {}
        for h2zone in [1, 2]:
            
            demand_fn = Path(
                fn,
                scenario,
                str(pyear),
                f"H2_ZONE_{h2zone}.xlsx",
            )
            
            demand[h2zone] = read_h2_excel(demand_fn, scenario, pyear, cyear, h2zone=h2zone)
        demand = pd.concat(demand, axis=1).droplevel(0, axis=1)
        
    
    return demand

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_tyndp_h2_demand",
                                   planning_horizons="2040",
                                   configfiles="config/test/config.tyndp.yaml")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    scenario = "NT" # snakemake.params["scenario"]
    planning_horizons = snakemake.params["planning_horizons"]
    snapshots = get_snapshots(snakemake.params.snapshots)
    cyear = 2009 # snapshots[0].year
    fn = snakemake.input.h2_demand

    logger.info(
        f"Processing H2 demand for scenario: {scenario}, planning year: {pyear}, weather year: {cyear}"
    )

    # Load and prep hydrogen demand
    tqdm_kwargs = {
        "ascii": False,
        "unit": " pyear",
        "total": len(planning_horizons),
        "desc": "Loading TYNDP hydrogen demand data",
    }

    func = partial(
        load_h2_demand,
        fn=fn,
        scenario=scenario,
        pyear=pyear,
        cyear=cyear,
        snapshots=snapshots,
    )

    with mp.Pool(processes=snakemake.threads) as pool:
        demand = list(tqdm(pool.imap(func, planning_horizons), **tqdm_kwargs))

    # Combine zones horizontally (concatenate columns)
    h2_demand = pd.concat(h2_demand_zones, axis=1)

    # Export to CSV
    h2_demand.to_csv(snakemake.output.h2_demand, index=True)


