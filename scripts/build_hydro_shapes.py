# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Expose the pypsa-eur clustered onshore regions to `module_hydropower`.

A module shape is a network bus. The module partitions its work by shape and
aggregates inflow per shape, so handing it the clustered regions makes its
`aggregated/{plant_type}_inflow_pu.parquet` output already resolved per bus.
Region holes are filled the same way `build_powerplants` fills them, so plants
sitting in a gap between regions are not lost.

`country_id` carries the code EIA publishes rather than a plain alpha-3
conversion. The module builds its generation series id as
`INTL.33-12-{country_id}-BKWH.A`; EIA publishes Kosovo as `XKS`, while
`country_converter` maps `XK` to `XKX`. Without the override the module's
`prepare_statistics` raises for Kosovo, a country the per-country path models.

Coverage is asserted here rather than filtered away: a shape is a bus, so
dropping shapes would delete buses from the network. A country whose series is
genuinely missing therefore has to fail loudly at this point instead of losing
its buses silently inside the module.

Outputs
-------

- `resources/modules/hydropower/base_s_{clusters}/shapes.parquet`
"""

import logging
import re
import zipfile

import country_converter as coco
import geopandas as gpd
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.build_powerplants import fill_unoccupied_holes

logger = logging.getLogger(__name__)

# EIA International bulk series for annual hydropower generation, e.g.
# "INTL.33-12-DEU-BKWH.A" (33 = hydropower, BKWH = billion kWh, A = annual).
EIA_HYDRO_SERIES = re.compile(rb'"series_id":"INTL\.33-12-([A-Z]{3})-BKWH\.A"')

# Alpha-3 codes EIA publishes under a different code than `country_converter`
# returns. Verified against the module's bulk archive (225 hydropower series):
# `XKX` is absent, `XKS` is present.
EIA_COUNTRY_CODES = {"XKX": "XKS"}


def eia_hydro_countries(eia_bulk_fn: str) -> set[str]:
    """Scan the EIA bulk archive for countries with a hydropower series."""
    with zipfile.ZipFile(eia_bulk_fn) as zf, zf.open("INTL.txt") as f:
        return {
            m.group(1).decode() for line in f if (m := EIA_HYDRO_SERIES.search(line))
        }


def build_hydro_shapes(regions_fn: str, eia_bulk_fn: str, output_fn: str) -> None:
    cc = coco.CountryConverter()
    regions = gpd.read_file(regions_fn).dissolve("name")
    regions["geometry"] = fill_unoccupied_holes(regions)

    # Bus names are prefixed with the alpha-2 country code, the same convention
    # `map_to_country_bus` relies on to keep plants within their own country.
    country_id = cc.pandas_convert(
        pd.Series(regions.index.str[:2], index=regions.index), to="ISO3"
    ).replace(EIA_COUNTRY_CODES)

    shapes = gpd.GeoDataFrame(
        {
            "shape_id": regions.index,
            "country_id": country_id,
            "shape_class": "land",
        },
        geometry=regions.geometry,
        crs=regions.crs,
    )

    uncovered = sorted(set(shapes["country_id"]) - eia_hydro_countries(eia_bulk_fn))
    if uncovered:
        raise ValueError(
            f"No EIA hydropower generation series for {uncovered}, which "
            "module_hydropower needs to normalise inflow. Every modelled "
            "country must be covered, because dropping its shapes would remove "
            f"its buses from the network. Check {EIA_COUNTRY_CODES} for a code "
            "that EIA publishes under a different name."
        )

    shapes.reset_index(drop=True).to_parquet(output_fn)
    logger.info(
        f"Wrote {len(shapes)} onshore regions across "
        f"{shapes['country_id'].nunique()} countries to {output_fn}"
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_hydro_shapes", clusters=50)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    build_hydro_shapes(
        snakemake.input.regions_onshore,
        snakemake.input.eia_bulk,
        snakemake.output.shapes,
    )
