# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

# Module repository: https://github.com/modelblocks-org/module_hydropower
# (see also https://www.modelblocks.org/modules/)

from scripts._helpers import get_snapshots

MODULE_NAME = "hydropower"
HYDRO_DIR = f"resources/modules/{MODULE_NAME}"
HYDRO_SCENARIO = config["modules"][MODULE_NAME]["scenario"]
HYDRO_SHAPES = f"{HYDRO_DIR}/{HYDRO_SCENARIO}_shapes.parquet"


def hydropower_module_config() -> dict:
    """Module configuration with the keys left `null` in the module config filled in.

    `years` is derived from the pypsa-eur snapshots, and `crs` is taken from
    `modules.crs`, so it is defined once for every composed module rather than
    repeated in each module config.
    """
    snapshot_years = get_snapshots(
        config["snapshots"], config["enable"]["drop_leap_day"]
    ).year.unique()
    module_cfg = module_config(MODULE_NAME)
    module_cfg["years"] = dict(
        start=int(snapshot_years.min()),
        end=int(snapshot_years.max()) + 1,
    )
    module_cfg["crs"] = dict(config["modules"]["crs"])
    return module_cfg


module hydropower:
    pathvars:
        shapes=HYDRO_SHAPES,
        powerplants=f"{HYDRO_DIR}/{{shapes}}/powerplants.parquet",
        disaggregated_inflow=f"{HYDRO_DIR}/{{shapes}}/disaggregated/inflow_mwh.parquet",
        logs=f"logs/modules/{MODULE_NAME}",
        resources=f"data/modules/{MODULE_NAME}",
        results=f"{HYDRO_DIR}/results",
    snakefile:
        github(
            "modelblocks-org/module_hydropower",
            path="workflow/Snakefile",
            tag=config["modules"][MODULE_NAME]["version"],
        )
    config:
        hydropower_module_config()


use rule * from hydropower exclude all as hydropower_*
