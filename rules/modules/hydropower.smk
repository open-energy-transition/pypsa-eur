# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

# Module repository: https://github.com/modelblocks-org/module_hydropower
# (see also https://www.modelblocks.org/modules/)

from scripts._helpers import get_snapshots

MODULE_NAME = "hydropower"
HYDRO_DIR = f"resources/modules/{MODULE_NAME}"
HYDRO_PARTITION = "base_s_{clusters}"
HYDRO_PLANT_TYPES = {"ror": "run_of_river", "hydro": "reservoir"}


def hydro_module_inflow_pu(plant_type: str) -> str:
    """Per-bus, per-unit inflow aggregated by the module for one plant type."""
    return f"{HYDRO_DIR}/{HYDRO_PARTITION}/aggregated/{plant_type}_inflow_pu.parquet"


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
        shapes=f"{HYDRO_DIR}/{{shapes}}/shapes.parquet",
        powerplants=f"{HYDRO_DIR}/{{shapes}}/powerplants.parquet",
        disaggregated_inflow=f"{HYDRO_DIR}/{{shapes}}/disaggregated/inflow_mwh.parquet",
        aggregated_inflow_pu=f"{HYDRO_DIR}/{{shapes}}/aggregated/{{plant_type}}_inflow_pu.parquet",
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


rule build_hydro_shapes:
    input:
        regions_onshore=resources(f"regions_onshore_{HYDRO_PARTITION}.geojson"),
        eia_bulk=rules.hydropower_download_eia.output["zipfile"],
    output:
        shapes=f"{HYDRO_DIR}/{HYDRO_PARTITION}/shapes.parquet",
    log:
        logs("build_hydro_shapes_{clusters}.log"),
    benchmark:
        benchmarks("build_hydro_shapes_{clusters}")
    threads: 1
    resources:
        mem_mb=2000,
    message:
        "Exposing {wildcards.clusters} onshore regions to module_hydropower as shapes"
    script:
        scripts("build_hydro_shapes.py")


rule build_hydro_powerplants:
    input:
        powerplants=rules.retrieve_powerplants.output["powerplants"],
        shapes=rules.build_hydro_shapes.output["shapes"],
    output:
        powerplants=f"{HYDRO_DIR}/{HYDRO_PARTITION}/powerplants.parquet",
    log:
        logs("build_hydro_powerplants_{clusters}.log"),
    benchmark:
        benchmarks("build_hydro_powerplants_{clusters}")
    threads: 1
    resources:
        mem_mb=2000,
    message:
        "Preparing hydro powerplants for module_hydropower"
    script:
        scripts("build_hydro_powerplants.py")
