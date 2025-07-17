<!--
SPDX-FileCopyrightText: Open Energy Transition gGmbH, Ember, and contributors to the Ember Flexibility Study
SPDX-License-Identifier: CC-BY-4.0
-->

# Ember Flexibility Study

<img src="https://raw.githubusercontent.com/open-energy-transition/Ember-Flexibility-Study/refs/heads/master/doc/img/ember_oet.png" alt="Open Energy Transition Logo" width="200" height="234" align="right">

This repository contains the code and analysis for the **Ember Flexibility Study**, conducted in collaboration with [Ember](https://ember-climate.org/) and [Open Energy Transition (OET)](https://openenergytransition.org/). The study investigates clean flexibility options for Europe's energy system, building on the PyPSA-Eur framework. All results are computed from raw data and code to ensure full reproducibility.

This repository is a soft-fork of [OET-PyPSA-Eur](https://github.com/open-energy-transition/pypsa-eur) and contains the entire project **Clean Flexibility for Europe's Energy System** supported by [Open Energy Transition (OET)](https://openenergytransition.org/)<sup>*</sup>, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

This repository is maintained using [OET's soft-fork strategy](https://open-energy-transition.github.io/handbook/docs/Engineering/SoftForkStrategy). OET's primary aim is to contribute as much as possible to the open source (OS) upstream repositories. For long-term changes that cannot be directly merged upstream, the strategy organizes and maintains OET forks, ensuring they remain up-to-date and compatible with upstream, while also supporting future contributions back to the OS repositories.

---

# Repository structure

* `benchmarks`: will store `snakemake` benchmarks (does not exist initially)
* `config`: configurations used in the study
* `cutouts`: will store raw weather data cutouts from `atlite` (does not exist initially)
* `data`: includes input data that is not produced by any `snakemake` rule
* `doc`: includes all files necessary to build the `readthedocs` documentation of PyPSA-Eur
* `envs`: includes all the `mamba` environment specifications to run the workflow
* `logs`: will store log files (does not exist initially)
* `notebooks`: includes all the `notebooks` used for ad-hoc analysis
* `report`: contains all files necessary to build the report; plots and result files are generated automatically
* `rules`: includes all the `snakemake`rules loaded in the `Snakefile`
* `resources`: will store intermediate results of the workflow which can be picked up again by subsequent rules (does not exist initially)
* `results`: will store the solved PyPSA network data, summary files and plots (does not exist initially)
* `scripts`: includes all the Python scripts executed by the `snakemake` rules to build the model

# Installation and Usage

## 1. Fork and clone the repository

First fork the repository [Ember-Flexibility-Study](https://github.com/open-energy-transition/Ember-Flexibility-Study/) on GitHub to your own account. Please make sure to check the box `Copy the master branch only`. Then, clone your fork locally:

    git clone https://github.com/<your-username>/Ember-Flexibility-Study.git

Once you have cloned your fork, you should add the following upstream remotes to keep your repository up to date with the main projects:

- Add the main [Ember-Flexibility-Study](https://github.com/open-energy-transition/Ember-Flexibility-Study/) repository as `upstream`:

      git remote add upstream https://github.com/open-energy-transition/Ember-Flexibility-Study.git

- Add the main OET soft-fork of PyPSA-Eur as `upstream_pypsa_eur_oet`:

      git remote add upstream_pypsa_eur_oet https://github.com/open-energy-transition/pypsa-eur.git

This setup allows you to fetch and integrate changes from both the main study repository and the OET soft-fork of PyPSA-Eur.

## 2. Merging changes from upstream repositories

To keep your fork up to date, you can merge changes from the master branch of either `upstream_pypsa_eur_oet` or `upstream` as follows:

- To merge changes from the OET soft-fork of PyPSA-Eur:

      git fetch upstream_pypsa_eur_oet
      git merge upstream_pypsa_eur_oet/master

- To merge changes from the main Ember-Flexibility-Study repository:

      git fetch upstream
      git merge upstream/master

Resolve any conflicts if they arise, then push the updates to your fork if needed.

## 3. Installation

Clone the repository:

    git clone https://github.com/open-energy-transition/{{repository}}

You need [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Users may also prefer to use [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) or [conda](https://docs.conda.io/projects/conda/en/latest/index.html). Using `mamba`, you can create an environment from within you can run it:

    mamba env create -f environment.yaml

Activate the newly created `{{project_short_name}}` environment:

    mamba activate {{project_short_name}}

## 4. Run the analysis

    snakemake -call

This will run all analysis steps to reproduce results and build the report.

To generate a PDF of the dependency graph of all steps `resources/dag.pdf` run:

    snakemake -c1 dag

---

# What is Snakemake?

[Snakemake](https://snakemake.readthedocs.io/) is a workflow management system that enables reproducible and scalable data analyses. It allows you to define complex pipelines in a readable Python-based language, automatically handling dependencies, job execution, and resource management. Snakemake is widely used in scientific computing for automating data processing, analysis, and reporting.

## Defining Rules in Snakemake

Snakemake workflows are built from modular units called **rules**. Each rule specifies how to create output files from input files, using scripts or shell commands. Rules define the steps of your workflow and their dependencies, making it easy to manage complex pipelines.

## Main Snakemake Command-Line Keys

Here are some of the most important command-line options (keys) to control the workflow:

- `-j`, `--jobs [N]`: Set the maximum number of jobs to run in parallel (e.g., `-j 4`).
- `-c`, `--cores [N]`: Specify the number of CPU cores to use (e.g., `-c 1`).
- `-n`, `--dryrun`: Show what would be executed without actually running the workflow.
- `-s`, `--snakefile [FILE]`: Specify a custom Snakefile (default is `Snakefile`).
- `-R`, `--rerun-incomplete`: Re-run jobs with incomplete output files.
- `--unlock`: Unlock the working directory if a previous run was interrupted.
- `--dag`: Print the directed acyclic graph (DAG) of jobs in the workflow.
- `--forceall`: Force the execution of all rules, regardless of output file timestamps.
- `-k`, `--keep-going`: Continue as much as possible after an error.
- `--config [KEY=VALUE,...]`: Override config file values from the command line.

For a full list of options, see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) or run `snakemake --help`.

## Important Rules in PyPSA-Eur

In [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur), some of the most important rules that structure the workflow include:

- **retrieve**: Downloads and prepares all required input data.
- **build_network**: Constructs the base energy system network from input data.
- **prepare_sector**: Prepares sector-coupling data (e.g., heating, transport).
- **solve_network**: Runs the optimization to solve the energy system model.
- **postprocess**: Processes and analyzes the results after solving.
- **plot_network**: Generates plots and visualizations from the results.
- **report**: Builds the final report or documentation from the results.

These rules are typically defined in separate `.smk` files (e.g., `rules/retrieve.smk`, `rules/build_electricity.smk`) and are orchestrated by the main `Snakefile`.

## How to run the workflow
This section builds upon the detailed description of the available PyPSA-Eur [configurations](https://pypsa-eur.readthedocs.io/en/latest/configuration.html) and [wildcards](https://pypsa-eur.readthedocs.io/en/latest/wildcards.html), and the tutorials for the [electricy-only](https://pypsa-eur.readthedocs.io/en/latest/tutorial.html) and for the [sector-coupled](https://pypsa-eur.readthedocs.io/en/latest/tutorial_sector.html) models.

You can use the `config/validation_config_2023.yaml` to run the workflow with different configurations. The `validate_ember_networks` rule will run the validation and trigger the plotting routines for the Ember Flexibility Study. The command is as follows:

`snakemake validate_ember_networks --configfile config/validation_config_2023.yaml`

---

# Contributing and Support

We strongly welcome anyone interested in contributing to this project. If you have any ideas, suggestions or encounter problems, feel invited to file issues or make pull requests on GitHub.

## Issue a pull request and merging it
To issue a pull request to the `master` branch of the upstream repository [Ember-Flexibility-Study](https://github.com/open-energy-transition/Ember-Flexibility-Study/), please follow the [instructions](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) and follow the instructions from the pull request [template](https://github.com/open-energy-transition/Ember-Flexibility-Study/blob/master/.github/pull_request_template.md).

## Raise issues, bugs or feature requests
For **issues, bugs and feature requests**, please use the [GitHub Issues page](https://github.com/open-energy-transition/{{repository}}/issues).

# Licence

The code in this repository is released as free software under the [MIT License](https://opensource.org/licenses/MIT), see [`doc/licenses.rst`](doc/licenses.rst). However, different licenses and terms of use may apply to the various input data, see [`doc/data_sources.rst`](doc/data_sources.rst).