import pandas as pd


wildcard_constraints:
    cba_project=r"(stor|trans)\d+",


rule retrieve_tyndp_cba_projects:
    params:
        # TODO Integrate into Zenodo tyndp data bundle
        url="https://storage.googleapis.com/open-tyndp-data-store/CBA_projects.zip",
        source="CBA project explorer",
    input:
        "data/tyndp_2024_bundle",
    output:
        dir=directory("data/tyndp_2024_bundle/cba_projects"),
    log:
        "logs/retrieve_tyndp_cba_projects",
    retries: 2
    script:
        "../../scripts/retrieve_additional_tyndp_data.py"


# read in transmission and storage projects from excel sheets they should get a
# project_name column with trans{num} or stor{num} i'd start with transmission projects
# only and then if we see we need to extract different information for storage projects
# we split those out into another file instead and change the workflow accordingly
checkpoint read_projects:
    input:
        dir="data/tyndp_2024_bundle/cba_projects",
    output:
        projects=resources("cba/projects.csv"),
    script:
        "../../scripts/cba/read_projects.py"


def input_sb_network(w):
    scenario = config_provider("scenario")(w)
    expanded_wildcards = {
        "clusters": scenario["clusters"],
        "opts": scenario["opts"],
        "sector_opts": scenario["sector_opts"],
    }
    match config_provider("foresight")(w):
        case "perfect":
            expanded_wildcards["planning_horizons"] = "all"
        case "myopic":
            pass
        case _:
            raise ValueError('config["foresight"] must be one of "perfect" or "myopic"')

    return expand(
        RESULTS
        + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        **expanded_wildcards,
        ignore_missing=True,
    )


# extract a planning horizon from the SB optimized network and apply the simplifications
# necessary to get to the general CBA reference network
rule simplify_sb_network:
    input:
        network=input_sb_network,
        # TODO whatever additional data is needed
    output:
        network=resources("cba/networks/simple_{planning_horizons}.nc"),
    script:
        "../../scripts/cba/simplify_sb_network.py"


# build the reference network for toot with all projects included
# maybe worth to merge with pint rule, if they turn out to be very similar
rule prepare_toot_reference:
    input:
        network=rules.simplify_sb_network.output.network,
        projects=rules.read_projects.output.projects,
    output:
        network=resources("cba/toot/networks/reference_{planning_horizons}.nc"),
    script:
        "../../scripts/cba/prepare_toot_reference.py"


# build the reference network for pint with all projects removed
# maybe worth to merge with toot rule, if they turn out to be very similar
rule prepare_pint_reference:
    input:
        network=rules.simplify_sb_network.output.network,
        projects=rules.read_projects.output.projects,
    output:
        network=resources("cba/pint/networks/reference_{planning_horizons}.nc"),
    script:
        "../../scripts/cba/prepare_pint_reference.py"


# remove the single project {cba_project} from the toot reference network
# currently this can be either a trans123 or a stor123 project
rule prepare_toot_project:
    input:
        network=rules.prepare_toot_reference.output.network,
        projects=rules.read_projects.output.projects,
    output:
        network=resources(
            "cba/toot/networks/project_{cba_project}_{planning_horizons}.nc"
        ),
    script:
        "../../scripts/cba/prepare_toot_project.py"


# add the single project {cba_project} to the pint reference network
# currently this can be either a trans123 or a stor123 project
rule prepare_pint_project:
    input:
        network=rules.prepare_pint_reference.output.network,
        projects=rules.read_projects.output.projects,
    output:
        network=resources(
            "cba/pint/networks/project_{cba_project}_{planning_horizons}.nc"
        ),
    script:
        "../../scripts/cba/prepare_pint_project.py"


# solve any of the prepared networks, ie a reference or a project network
# should reuse/import functions from solve_network.py
rule solve_cba_network:
    input:
        network=resources("cba/{cba_method}/networks/{name}_{planning_horizons}.nc"),
    output:
        network=resources("cba/{cba_method}/postnetworks/{name}_{planning_horizons}.nc"),
    script:
        "../../scripts/cba/solve_cba_network.py"


# compute all metrics for a single pint or toot project comparing reference and project solution
rule compute_metrics:
    input:
        reference=resources(
            "cba/{cba_method}/postnetworks/reference_{planning_horizons}.nc"
        ),
        project=resources(
            "cba/{cba_method}/postnetworks/project_{cba_project}_{planning_horizons}.nc"
        ),
    output:
        metrics=RESULTS
        + "cba/{cba_method}/project_{cba_project}_{planning_horizons}.csv",
    script:
        "../../scripts/cba/compute_metrics.py"


def input_metrics(w):
    """
    List all metric csv
    """
    projects = pd.read_csv(checkpoints.read_projects.get(**w).output.projects)
    planning_horizons = config_provider("scenario", "planning_horizons")(w)

    # assumes a one row per project csv file, maybe needs to be changed if we use
    # a tidy format for the projects file.
    cba_projects = projects.loc[projects["type"] == "transmission"].itertuples()

    return expand(
        rules.compute_metrics.output.metrics,
        cba_project=cba_projects,
        cba_method=["toot", "pint"],  # maybe promote to config
        planning_horizons=planning_horizons,
    )


# assemble the metrics for all projects into a single overview csv (ideally just concatenating it)
rule assemble_metrics:
    input:
        metrics=input_metrics,
    output:
        metrics=RESULTS + f"cba/metrics.csv",
    script:
        "../../scripts/cba/assemble_metrics.py"


# pseudo-rule, to run enable running cba with snakemake cba --configfile config/config.tyndp.yaml
rule cba:
    input:
        rules.assemble_metrics.output.metrics,
