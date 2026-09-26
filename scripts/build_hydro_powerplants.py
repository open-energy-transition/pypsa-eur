# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Convert the pypsa-eur hydro fleet into a `module_hydropower` input.

The fleet is `powerplants_s_{clusters}.csv`, the table `add_electricity` builds
the network from, so the module and the network see the same plants.
Technologies are renamed to carriers with `renewable.hydro.technology_mapping`,
the mapping `add_electricity` uses. Hydro plants without a mapped technology
(missing in powerplantmatching) are left out here, as they are left out of the
network; the module would reject them anyway, since `technology` may not be
empty. Pumped storage is passed on: the module keeps only the plant types in
its own `technology_mapping` and ignores the rest.

The output follows the module's `PowerplantSchema`. Extra columns are filtered
out by the schema, and the module assigns each plant to a shape itself.

Outputs
-------

- `resources/modules/hydropower/base_s_{clusters}/powerplants.parquet`
"""

import logging

import geopandas as gpd
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def build_hydro_powerplants(
    powerplants_fn: str, technology_mapping: dict[str, str], output_fn: str
) -> None:
    ppl = pd.read_csv(powerplants_fn, index_col=0)
    ppl = ppl[ppl["Fueltype"].eq("Hydro")]
    carrier = ppl["Technology"].map(technology_mapping)

    unmapped = carrier.isna()
    if unmapped.any():
        logger.info(
            f"Leaving out {unmapped.sum()} hydro plants "
            f"({ppl.loc[unmapped, 'Capacity'].sum() / 1e3:.2f} GW) without a "
            "technology in `renewable.hydro.technology_mapping`; they are not "
            "attached to the network either."
        )
    ppl, carrier = ppl[~unmapped], carrier[~unmapped]
    if ppl.empty:
        raise ValueError(f"No mapped hydro powerplants in {powerplants_fn}.")

    out = gpd.GeoDataFrame(
        {
            "powerplant_id": ppl.index.astype(str),
            "output_capacity_mw": ppl["Capacity"],
            "technology": carrier,
            "start_year": ppl["DateIn"].fillna(0),
            "end_year": ppl["DateOut"].fillna(9999),
        },
        geometry=gpd.points_from_xy(ppl["lon"], ppl["lat"]),
        crs="EPSG:4326",
    ).dropna(subset=["output_capacity_mw"])

    logger.info(
        f"Writing {len(out)} hydro powerplants "
        f"({out['output_capacity_mw'].sum() / 1e3:.1f} GW) to {output_fn}\n"
        f"{out['technology'].value_counts().to_string()}"
    )
    out.reset_index(drop=True).to_parquet(output_fn)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_hydro_powerplants", clusters=50)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    build_hydro_powerplants(
        powerplants_fn=snakemake.input.powerplants,
        technology_mapping=snakemake.params.technology_mapping,
        output_fn=snakemake.output.powerplants,
    )
