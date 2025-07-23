# SPDX-FileCopyrightText: Open Energy Transition gGmbH, Ember, and contributors to the Ember Flexibility Study
#
# SPDX-License-Identifier: MIT
import time
# rules/emberdata.smk

from pathlib import Path

import requests

DOWNLOADS = {
    Path("validation", "ember_data", "yearly_full_release_long_format.csv"):
        "https://storage.googleapis.com/emb-prod-bkt-publicdata/public-downloads/yearly_full_release_long_format.csv",
    Path("validation", "ember_data", "europe_monthly_full_release_long_format.csv"):
        "https://storage.googleapis.com/emb-prod-bkt-publicdata/public-downloads/europe_monthly_full_release_long_format.csv",
    Path("validation", "entsoe_data", "physical_energy_power_flows_2023.csv"):
        "https://www.entsoe.eu/publications/data/power-stats/2023/physical_energy_power_flows_2023.csv"
}

rule download_ember_data:
    output:
        [str(path) for path in DOWNLOADS.keys()]
    run:
        import urllib.request
        import yaml

        config_path = Path("config", "validation_2023.yaml")

        # Load config if it exists
        cfg = {}
        if config_path.exists():
            with config_path.open("r") as f:
                cfg = yaml.safe_load(f)

        for filepath, url in DOWNLOADS.items():
            # Create directory if it doesn't exist
            filepath.parent.mkdir(parents=True, exist_ok=True)

            # Download file if it doesn't already exist
            if not filepath.exists():
                logger.info(f"Downloading {url} -> {filepath}")
                response = requests.get(url)
                response.raise_for_status()  # Raise an error for non-200 responses
                with open(filepath,"wb") as f:
                    f.write(response.content)

                # Confirm file creation
                while not filepath.exists():
                    logger.info(f"Waiting for {filepath} to appear...")
                    time.sleep(1)

                else:
                    logger.info(f"Already exists: {filepath}")

rule download_ember_NTC_data:
    output:
        file="validation/ember_data/Reg_NTC"
    shell:
        """
        gdown https://drive.google.com/uc?id=1GTo4UrI_X9ZCsgtM4KobO_pw-TTrEoDy -O {output.file}
        """

rule download_eurostat:
    output:
        "validation/eurostatdata/eurostat_nrg_bal_c_2023.csv"
    shell:
        """
        curl -L "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/nrg_bal_c?format=SDMX-CSV&startPeriod=2023&endPeriod=2023&lang=en&geo=EU27_2020&unit=KTOE" -o {output}
        """

rule download_jrc_idees:
    output:
        "validation/eurostatdata/JRC-IDEES-2021_EU27.zip"
    shell:
        """
        curl -L "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/JRC-IDEES/JRC-IDEES-2021_v1/JRC-IDEES-2021_EU27.zip" -o {output}
        """

rule download_hotmaps:
    output:
        "validation/eurostatdata/Industrial_Database.csv"
    shell:
        """
        curl -L "https://gitlab.com/hotmaps/industrial_sites/industrial_sites_Industrial_Database/-/raw/master/data/Industrial_Database.csv" -o {output}
        """

rule extract_jrc_idees:
    input:
        "validation/eurostatdata/JRC-IDEES-2021_EU27.zip"
    output:
        directory("validation/eurostatdata/jrc_idees/")
    shell:
        """
        unzip {input} -d {output}
        """