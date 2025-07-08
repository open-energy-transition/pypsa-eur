# rules/emberdata.smk

DOWNLOADS = {
    "validation/ember_data/yearly_full_release_long_format.csv":
        "https://storage.googleapis.com/emb-prod-bkt-publicdata/public-downloads/yearly_full_release_long_format.csv",
    "validation/ember_data/europe_monthly_full_release_long_format.csv":
        "https://storage.googleapis.com/emb-prod-bkt-publicdata/public-downloads/europe_monthly_full_release_long_format.csv"
}

rule download_ember_data:
    output:
        list(DOWNLOADS.keys())
    run:
        import os
        import urllib.request
        import yaml

        config_path = "config/validation_2023.yaml"
        if os.path.exists(config_path):
            with open(config_path, "r") as f:
                cfg = yaml.safe_load(f)
        else:
            cfg = {}

        if cfg.get("download_ember_data", False):
            for filepath, url in DOWNLOADS.items():
                os.makedirs(os.path.dirname(filepath), exist_ok=True)
                if not os.path.exists(filepath):
                    print(f"Downloading {url} -> {filepath}")
                    urllib.request.urlretrieve(url, filepath)
                else:
                    print(f"Already exists: {filepath}")
        else:
            print("Skipping ember data download (flag is false or not set).")
