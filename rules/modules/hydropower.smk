# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from scripts._helpers import get_snapshots

MODULE_NAME = "hydropower"
HYDRO_DIR = f"resources/modules/{MODULE_NAME}"
HYDRO_SCENARIO = config["modules"][MODULE_NAME]["scenario"]
HYDRO_SHAPES = f"{HYDRO_DIR}/{HYDRO_SCENARIO}_shapes.parquet"


def validate_hydropower_years():
    if config["renewable"]["hydro"]["source"] != "module":
        return
    years = module_config(MODULE_NAME)["years"]
    snapshot_years = get_snapshots(
        config["snapshots"], config["enable"]["drop_leap_day"]
    ).year.unique()
    missing = sorted(set(snapshot_years) - set(range(years["start"], years["end"])))
    if missing:
        raise ValueError(
            f"module_hydropower covers {years['start']}-{years['end'] - 1} but the "
            f"snapshots need {missing}. Adjust `years` in "
            f"{config['modules'][MODULE_NAME]['config_path']}."
        )


validate_hydropower_years()


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
        module_config(MODULE_NAME)


use rule * from hydropower exclude all as hydropower_*
