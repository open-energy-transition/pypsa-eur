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
                with open(filepath,"wb") as f:
                    f.write(response.content)

                # Confirm file creation
                while not filepath.exists():
                    logger.info(f"Waiting for {filepath} to appear...")
                    time.sleep(1)

                else:
                    logger.info(f"Already exists: {filepath}")
