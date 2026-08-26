# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Expose the geo_boundaries land shapes to `module_hydropower`.

Countries without an EIA hydropower generation series are dropped: the module
normalises inflow against it and its `prepare_statistics` raises on the first
missing series.

Outputs
-------

- `resources/modules/hydropower/{shapes}_shapes.parquet`
"""

import logging
import re
import zipfile

import geopandas as gpd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

# EIA International bulk series for annual hydropower generation, e.g.
# "INTL.33-12-DEU-BKWH.A" (33 = hydropower, BKWH = billion kWh, A = annual).
EIA_HYDRO_SERIES = re.compile(rb'"series_id":"INTL\.33-12-([A-Z]{3})-BKWH\.A"')


def eia_hydro_countries(eia_bulk_fn: str) -> set[str]:
    """Scan the EIA bulk archive for countries with a hydropower series."""
    with zipfile.ZipFile(eia_bulk_fn) as zf, zf.open("INTL.txt") as f:
        return {
            m.group(1).decode() for line in f if (m := EIA_HYDRO_SERIES.search(line))
        }


def prepare_hydro_shapes(geo_shapes_fn: str, eia_bulk_fn: str, output_fn: str) -> None:
    shapes = gpd.read_parquet(geo_shapes_fn).query("shape_class == 'land'")
    covered = eia_hydro_countries(eia_bulk_fn)

    dropped = sorted(set(shapes["country_id"].dropna()) - covered)
    if dropped:
        logger.warning(
            f"Dropping {len(dropped)} country(ies) without EIA hydropower "
            f"generation, which module_hydropower cannot normalise: {dropped}"
        )
    shapes = shapes[shapes["country_id"].isin(covered)]

    shapes.to_parquet(output_fn)
    logger.info(
        f"Wrote {len(shapes)} land shapes across "
        f"{shapes['country_id'].nunique()} countries to {output_fn}"
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_hydro_shapes")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    prepare_hydro_shapes(
        snakemake.input.geo_shapes,
        snakemake.input.eia_bulk,
        snakemake.output.shapes,
    )
