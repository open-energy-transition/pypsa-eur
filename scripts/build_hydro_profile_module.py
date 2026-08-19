# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Aggregate the per-plant inflow of `module_hydropower` onto network buses.

The module resolves inflow at plant level, so plants are assigned to buses with
the same spatial join `build_powerplants` uses for the rest of the fleet, and
summed per bus and carrier.

Outputs
-------

- `resources/profile_{clusters}_hydro.nc`:

| Field | Dimensions | Description |
| --- | --- | --- |
| inflow | time, bus, carrier | Inflow to the state of charge (in MW), e.g. due to river inflow in hydro reservoir. |
"""

import logging

import geopandas as gpd
import pandas as pd
import xarray as xr

from scripts._helpers import configure_logging, get_snapshots, set_scenario_config
from scripts.build_powerplants import fill_unoccupied_holes, map_to_country_bus

logger = logging.getLogger(__name__)

CARRIER = {"Run-Of-River": "ror", "Reservoir": "hydro"}


def build_hydro_profile_module(
    inflow_mwh_fn: str,
    powerplants_fn: str,
    regions_fn: str,
    output_fn: str,
    snapshots: pd.DatetimeIndex,
) -> None:
    regions = gpd.read_file(regions_fn).dissolve("name")
    regions["geometry"] = fill_unoccupied_holes(regions)

    ppl = gpd.read_parquet(powerplants_fn).rename(columns={"country": "Country"})
    ppl["carrier"] = ppl["technology"].map(CARRIER)

    inflow = pd.read_parquet(inflow_mwh_fn)
    ppl = map_to_country_bus(ppl, regions).set_index("powerplant_id")
    ppl = ppl.reindex(inflow.columns)

    dropped = ppl[["bus", "carrier"]].isna().any(axis="columns")
    if dropped.any():
        stats = ppl[dropped].groupby("Country")["output_capacity_mw"].sum().div(1e3)
        logger.warning(
            f"Dropping {dropped.sum()} plants outside the network regions "
            f"({stats.sum():.1f} GW) [GW per country]\n{stats.round(2).to_string()}"
        )
    ppl = ppl[~dropped]

    per_bus = (
        inflow[ppl.index]
        .T.groupby([ppl["bus"], ppl["carrier"]])
        .sum()
        .T.rename_axis(index="time")
    )
    profile = (
        xr.DataArray(per_bus, name="inflow")
        .unstack("dim_1")
        .fillna(0.0)
        .transpose("time", "bus", "carrier")
    )
    profile.attrs = {"long_name": "Hydropower inflow", "units": "MW"}

    profile = profile.sel(time=snapshots)

    profile.to_netcdf(output_fn)
    logger.info(
        f"Wrote inflow for {profile.sizes['bus']} buses "
        f"({ppl['output_capacity_mw'].sum() / 1e3:.1f} GW) to {output_fn}"
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_hydro_profile_module", clusters=50)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    snapshots = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )

    build_hydro_profile_module(
        snakemake.input.inflow_mwh,
        snakemake.input.powerplants,
        snakemake.input.regions_onshore,
        snakemake.output.profile,
        snapshots,
    )
