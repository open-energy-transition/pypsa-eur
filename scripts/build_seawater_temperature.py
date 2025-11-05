# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build seawater temperature data from Copernicus Marine Service.

This script downloads historical seawater temperature data for use in sea water
heat pump calculations. It retrieves potential temperature (thetao) data from
the global ocean physics reanalysis dataset at daily resolution.

The data covers European coastal areas at a spatial resolution of 0.083° and
includes near-surface depths (5-15m) suitable for heat pump applications.

Relevant Settings
-----------------

.. code:: yaml

    # No specific configuration required
    # Uses year wildcard from Snakemake rule

Inputs
------
- None (downloads from Copernicus Marine Service)

Outputs
-------
- `data/seawater_temperature_{year}.nc`: NetCDF file containing seawater temperature data

Notes
-----
Requires Copernicus Marine Service credentials configured via copernicusmarine package.
See https://marine.copernicus.eu/ for account setup and API access.
"""

import logging
import os

import copernicusmarine
import requests

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_seawater_temperature",
            clusters="39",
            opts="",
            ll="vopt",
            sector_opts="",
            planning_horizons=2050,
        )

    # Configure logging and scenario
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    if snakemake.input.data != []:
        if snakemake.params.default_cutout == "be-03-2013-era5":
            logger.info("Retrieving test-cutout seawater temperature data.")

            # url = "https://zenodo.org/records/15828866/files/seawater_temperature.nc"
            url = snakemake.input.seawater_temperature

            response = requests.get(url, stream=True)
            response.raise_for_status()

            with open(snakemake.output.data, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            logger.info(
                f"Successfully downloaded test-cutout seawater temperature data to {snakemake.output.seawater_temperature}"
            )
        else:
            raise ValueError (
                f"Archive data does not exist for cutout {snakemake.params.default_cutout}. Use primary source to download cutout"
            )
    else:
        # Download seawater temperature data from Copernicus Marine Service
        # Dataset: Global Ocean Physics Reanalysis (daily, 0.083° resolution)
        # Variable: thetao (potential temperature in °C)
        # Spatial coverage: European waters (-12°W to 42°E, 33°N to 72°N)
        # Depth range: 5-15m (suitable for heat pump intake depths)
        logger.info(
            f"Downloading seawater temperature data for year {snakemake.wildcards.year}"
        )

        copernicusmarine_cutout=copernicusmarine.subset(
            dataset_id=snakemake.params.dataset_id,
            start_datetime=f"{snakemake.wildcards.year}-01-01",
            end_datetime=f"{int(snakemake.wildcards.year)}-12-31",
            minimum_longitude=snakemake.params.longitude[0],  # Western European boundary
            maximum_longitude=snakemake.params.longitude[1],  # Eastern European boundary
            minimum_latitude=snakemake.params.latitude[0],  # Southern European boundary
            maximum_latitude=snakemake.params.latitude[1],  # Northern European boundary
            variables=snakemake.params.variables,  # Potential temperature [°C]
            minimum_depth=snakemake.params.depth[0],  # Near-surface depth for heat pumps [m]
            maximum_depth=snakemake.params.depth[1],  # Near-surface depth for heat pumps [m]
            output_filename=snakemake.output.seawater_temperature,
        )

        # Verify successful download
        if not os.path.exists(snakemake.output.seawater_temperature):
            raise FileNotFoundError(
                f"Failed to retrieve seawater temperature data and save to {snakemake.output.seawater_temperature}. "
                f"One reason might be missing Copernicus Marine login info. "
                f"See the copernicusmarine package documentation for details."
            )

        logger.info(
            f"Successfully downloaded seawater temperature data to {snakemake.output.seawater_temperature}"
        )
