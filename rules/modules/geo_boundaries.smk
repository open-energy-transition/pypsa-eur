# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

MODULE_NAME = "geo_boundaries"


def geo_boundaries_module_config() -> dict:
    """Module configuration with the keys left `null` in the module config filled in.

    `crs` is taken from `modules.crs`, so it is defined once for every composed
    module rather than repeated in each module config.
    """
    module_cfg = module_config(MODULE_NAME)
    module_cfg["crs"] = dict(config["modules"]["crs"])
    return module_cfg


module geo_boundaries:
    pathvars:
        shapes=f"resources/modules/{MODULE_NAME}/{{scenario}}.parquet",
        logs=f"logs/modules/{MODULE_NAME}",
        resources=f"data/modules/{MODULE_NAME}",
        results=f"resources/modules/{MODULE_NAME}/results",
    snakefile:
        github(
            "modelblocks-org/module_geo_boundaries",
            path="workflow/Snakefile",
            tag=config["modules"][MODULE_NAME]["version"],
        )
    config:
        geo_boundaries_module_config()


use rule * from geo_boundaries exclude all as geo_boundaries_*
