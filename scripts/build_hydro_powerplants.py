# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Convert the powerplantmatching hydro fleet into a `module_hydropower` input.

Pumped storage is excluded (closed loop, no natural inflow). Plants with an
unknown technology default to run-of-river: in ppm 0.8.1 all 197 such plants
have `Set=PP` with no storage evidence, and 91% of labelled non-PHS hydro PP
plants are run-of-river. Plants outside the shape countries are dropped here so
they do not consume the module's `max_dropped` budget.

The output follows the module's `PowerplantSchema`, plus a `country` column
(alpha-2) that the schema filters out but `build_hydro_profile_module` uses to
map plants onto buses.

Outputs
-------

- `resources/modules/hydropower/{shapes}/powerplants.parquet`
"""

import logging

import country_converter as coco
import geopandas as gpd
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

DEFAULT_TECHNOLOGY = "Run-Of-River"
EXCLUDED_TECHNOLOGY = "Pumped Storage"


def build_hydro_powerplants(
    powerplants_fn: str, shapes_fn: str, output_fn: str
) -> None:
    cc = coco.CountryConverter()
    countries = set(pd.read_parquet(shapes_fn, columns=["country_id"])["country_id"])

    ppl = pd.read_csv(powerplants_fn, index_col=0)
    ppl = ppl[ppl["Fueltype"].eq("Hydro") & ppl["Technology"].ne(EXCLUDED_TECHNOLOGY)]
    ppl = ppl[cc.pandas_convert(ppl["Country"], to="ISO3").isin(countries)]
    if ppl.empty:
        raise ValueError(
            f"No hydro powerplants in {powerplants_fn} for countries {sorted(countries)}."
        )

    out = gpd.GeoDataFrame(
        {
            "powerplant_id": ppl.index.astype(str),
            "output_capacity_mw": ppl["Capacity"],
            "technology": ppl["Technology"].fillna(DEFAULT_TECHNOLOGY),
            "start_year": ppl["DateIn"].fillna(0),
            "end_year": ppl["DateOut"].fillna(9999),
            "country": cc.pandas_convert(ppl["Country"], to="ISO2"),
        },
        geometry=gpd.points_from_xy(ppl["lon"], ppl["lat"]),
        crs="EPSG:4326",
    ).dropna(subset=["output_capacity_mw", "geometry"])

    logger.info(
        f"Writing {len(out)} of {len(ppl)} hydro powerplants "
        f"({out['output_capacity_mw'].sum() / 1e3:.1f} GW) to {output_fn}\n"
        f"{out['technology'].value_counts().to_string()}"
    )
    out.reset_index(drop=True).to_parquet(output_fn)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_hydro_powerplants", shapes="default")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    build_hydro_powerplants(
        snakemake.input.powerplants,
        snakemake.input.shapes,
        snakemake.output.powerplants,
    )
