# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Coordinate reference systems configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#crs
"""

from pydantic import Field

from scripts.lib.validation.config._base import ConfigModel


class CrsConfig(ConfigModel):
    """Configuration for `crs` settings."""

    projected: str = Field(
        "epsg:3035",
        description="Projected coordinate reference system, used for area and distance calculations. Defaults to ETRS89-extended / LAEA Europe, whose units are metres.",
    )
    geographic: str = Field(
        "epsg:4326",
        description="Geographic coordinate reference system, used for storing and exchanging coordinates. Defaults to WGS 84, whose units are degrees.",
    )
