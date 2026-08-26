# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Modules configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#modules
"""

from pathlib import Path

from pydantic import Field, FilePath

from scripts.lib.validation.config._base import ConfigModel


class _CrsConfig(ConfigModel):
    """Configuration for `modules.crs` settings."""

    projected: str = Field(
        "epsg:3035",
        description="Projected coordinate reference system handed to the modules, used for area and distance calculations. Defaults to ETRS89-extended / LAEA Europe, whose units are metres.",
    )
    geographic: str = Field(
        "epsg:4326",
        description="Geographic coordinate reference system handed to the modules, used for storing and exchanging coordinates. Defaults to WGS 84, whose units are degrees.",
    )


class _ModuleConfig(ConfigModel):
    """Configuration for module settings"""

    config_path: FilePath = Field(
        description="Path to module default configuration file. Can be relative to the project directory or absolute.",
    )
    version: str = Field(
        description="Module version to use.",
    )


class _GeoBoundariesModuleConfig(_ModuleConfig):
    """Configuration for module settings"""

    scenario: str = Field(
        default="default",
        description="Scenario to use for the geo_boundaries module.",
    )


class _HydropowerModuleConfig(_ModuleConfig):
    """Configuration for module settings"""

    scenario: str = Field(
        default="default",
        description="Scenario to use for the hydropower module.",
    )


class ModulesConfig(ConfigModel):
    """Configuration for modules."""

    crs: _CrsConfig = Field(
        default_factory=_CrsConfig,
        description="Coordinate reference systems shared by all composed modules, overriding the module defaults.",
    )
    geo_boundaries: _GeoBoundariesModuleConfig = Field(
        default=_GeoBoundariesModuleConfig(
            config_path=Path("config/modules/geo_boundaries.yaml"), version="v1.0.1"
        ),
        description="Configuration for the geo_boundaries module.",
    )
    hydropower: _HydropowerModuleConfig = Field(
        default=_HydropowerModuleConfig(
            config_path=Path("config/modules/hydropower.yaml"), version="v0.2.2"
        ),
        description="Configuration for the hydropower module, used when `renewable.hydro.source` is `module`.",
    )
