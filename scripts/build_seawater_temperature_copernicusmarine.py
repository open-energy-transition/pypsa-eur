# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build seawater temperature data from Copernicus Marine Service.

This script downloads historical seawater temperature data for use in sea water
heat pump calculations. It retrieves potential temperature (thetao) data from
the global ocean physics reanalysis dataset at daily resolution.

The data covers European coastal areas at a spatial resolution of 0.083° and
includes near-surface depths (depth range set in config) suitable for heat pump applications.

Relevant Settings
-----------------

.. code:: yaml

    copernicusmarine:
        dataset_id: ID of the dataset to download from copernicusmarine package
        variables: Potential temperature [°C]
        depth: Near-surface depth for heat pumps [m]

    # Uses cutout wildcard from Snakemake rule

Inputs
------
- None (downloads from Copernicus Marine Service)

Outputs
-------
- `data/seawater_temperature_copernicusmarine_{cutout}.nc`: NetCDF file containing seawater temperature data

Notes
-----
Requires Copernicus Marine Service credentials configured via copernicusmarine package.
For account setup and API access see this [documentation article](https://help.marine.copernicus.eu/en/articles/8185007-copernicus-marine-toolbox-credentials-configuration#h_f60e2fce95) on the [Copernicus Marine website](https://marine.copernicus.eu/).
"""

import logging
import os

import copernicusmarine

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)


def build_cutout(
    dataset_id: (str),
    latitude: (list),
    longitude: (list),
    variables: (list),
    depth: (list),
    year_range: (list),
    output_path: (str),
) -> None:
    """
    Download seawater temperature data from Copernicus Marine Service.

    Parameters
    ----------
    dataset_id: str
        ID of the Global Ocean Physics Reanalysis dataset (daily, 0.083° resolution)
    latitude: list
        Min. and max. latitude
    longitude: list
        Min. and max. longitude
    variables: list
        Variables to download. "thetao" refers to temperature
    depth: list
    Depth range
    year_range: list
    Years for which to download data from Jan through Dec
    output_path: str
        Output path to store temeperature data

    Notes
    -----
    Requires login for Copernicusmarine API.
    """
    _ = copernicusmarine.subset(
        dataset_id=dataset_id,
        start_datetime=f"{int(year_range[0])}-01-01",
        end_datetime=f"{int(year_range[1])}-12-31",
        minimum_longitude=longitude[0],
        maximum_longitude=longitude[1],
        minimum_latitude=latitude[0],
        maximum_latitude=latitude[1],
        variables=variables,
        minimum_depth=depth[0],
        maximum_depth=depth[1],
        output_filename=output_path,
    )

    # Verify successful download
    if not os.path.exists(output_path):
        raise FileNotFoundError(
            f"Failed to retrieve seawater temperature data and save to {output_path}. "
            f"One reason might be missing Copernicus Marine login info. "
            f"See the copernicusmarine package documentation for details."
        )

    logger.info(f"Successfully downloaded seawater temperature data to {output_path}")


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

    # Load snakemake parameters
    cutout = snakemake.wildcards.cutout
    cutout_dict = snakemake.params.cutout_dict

    # Build cutout
    build_cutout(
        dataset_id=snakemake.params.dataset_id,
        latitude=cutout_dict[cutout]["y"],
        longitude=cutout_dict[cutout]["x"],
        variables=snakemake.params.variables,
        depth=snakemake.params.depth,
        year_range=cutout_dict[cutout]["time"],
        output_path=snakemake.output.seawater_temperature,
    )

    logger.info(
        f"Downloading seawater temperature data for cutout '{cutout}' completed."
    )
