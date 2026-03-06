import argparse
import os
import pathlib
import re
import subprocess

import ruamel.yaml
from ruamel.yaml.scalarstring import DoubleQuotedScalarString

DEFAULT_BENCHMARKS = [
    "pypsa-eur-elec",
    "pypsa-eur-sec"
]


def load_yaml(path: pathlib.Path) -> dict:
    y = ruamel.yaml.YAML()
    with open(path, "r") as f:
        return y.load(f)


def save_yaml(path: pathlib.Path, data: dict):
    y = ruamel.yaml.YAML()
    y.preserve_quotes = True
    y.width = 4096
    with open(path, "w") as f:
        y.dump(data, f)


def recursive_merge(base: dict, override: dict) -> dict:
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            recursive_merge(base[k], v)
        else:
            base[k] = v
    return base


def select_yaml_files(benchmark_name: str, yaml_dir: pathlib.Path):
    benchmark_yaml_map = {
        "pypsa-eur-elec": [
            "pypsa-eur-elec.yaml",
#            "pypsa-eur-elec-uc.yaml",
#            "pypsa-eur-elec-dfp.yaml",
#            "pypsa-eur-elec-trex_copt.yaml",
#            "pypsa-eur-elec-trex_copt-dfp.yaml",
#            "pypsa-eur-elec-trex_vopt.yaml",
#            "pypsa-eur-elec-trex_vopt-dfp.yaml",
#            "pypsa-eur-elec-uc-dfp.yaml",
#            "pypsa-eur-elec-trex_copt-uc.yaml",
#            "pypsa-eur-elec-trex_copt-uc-dfp.yaml",
#            "pypsa-eur-elec-trex_vopt-uc.yaml",
#            "pypsa-eur-elec-trex_vopt-uc-dfp.yaml",
        ],
        "pypsa-eur-sec": [
#            "pypsa-eur-sec.yaml",
#            "pypsa-eur-sec-trex_copt.yaml",
#            "pypsa-eur-sec-trex_vopt.yaml",
        ],
    }

    matches = []
    for name in benchmark_yaml_map.get(benchmark_name, []):
        p = yaml_dir / name
        if p.exists():
            matches.append(p)

    if not matches:
        raise ValueError(f"No YAML found for {benchmark_name}")

    return matches


def validate_time_resolution(res: str) -> str:
    m = re.match(r"^(\d+)[hH]$", res)
    if not m:
        raise argparse.ArgumentTypeError(f"Invalid time resolution: {res}")
    return f"{m.group(1)}H"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--benchmark_name", choices=DEFAULT_BENCHMARKS, default="pypsa-eur-elec")
    p.add_argument("--file_extension", default=".lp")
    p.add_argument("--output_dir", default="/tmp/")
    p.add_argument("-n", "--dry_run", action="store_true")
    p.add_argument("-c", "--clusters", nargs="+", type=int, default=[41])
    p.add_argument("-r", "--time_resolutions", nargs="+", type=validate_time_resolution,
                   default=["1h"])
    return p.parse_args()

def run(cmd: str):
    proc = subprocess.Popen(
        cmd.split(),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    for line in proc.stdout:
        print(line, end="")
    proc.wait()

    if proc.returncode != 0:
        raise RuntimeError(f"Snakemake failed (rc={proc.returncode}): {cmd}")


def generate_benchmark(name: str,
                       config_path: pathlib.Path,
                       dry: bool,
                       clusters: int,
                       t_res: str):

    dry_flag = "-n" if dry else ""

    # ONLY FIRST-STAGE SOLVES ARE KEPT
    if "elec" in name:
        run(
            f"snakemake --snakefile Snakefile solve_elec_networks "
            f"--configfile {config_path} {dry_flag} --cores all"
        )
    else:
        run(
            f"snakemake --snakefile Snakefile solve_sector_networks "
            f"--configfile {config_path} {dry_flag} --cores all"
        )


if __name__ == "__main__":

    args = parse_args()

    base_dir = pathlib.Path(__file__).parent
    default_cfg = load_yaml(base_dir / "config" / "config.default.yaml")
    benchmark_files = select_yaml_files(args.benchmark_name, base_dir / "config")

    bench_type = "elec" if "elec" in args.benchmark_name else "sec"
    horizon = "2050"

    for bench_file in benchmark_files:
        bench_cfg = load_yaml(bench_file)

        for clusters in args.clusters:
            for t_res in args.time_resolutions:

                # Build merged config
                final_cfg = recursive_merge(default_cfg.copy(), bench_cfg.copy())

                final_cfg.setdefault("scenario", {})
                final_cfg["scenario"]["clusters"] = [clusters]
                final_cfg["scenario"]["planning_horizons"] = [horizon]

                final_cfg.setdefault("clustering", {})
                final_cfg["clustering"].setdefault("temporal", {})
                if bench_type == "elec":
                    final_cfg["clustering"]["temporal"]["resolution_elec"] = t_res
                else:
                    final_cfg["clustering"]["temporal"]["resolution_sector"] = t_res

                # Fix NO/False issue
                if "countries" in final_cfg:
                    cleaned = []
                    for c in final_cfg["countries"]:
                        if c is False or c == "NO":
                            cleaned.append(DoubleQuotedScalarString("NO"))
                        else:
                            cleaned.append(DoubleQuotedScalarString(str(c)))
                    final_cfg["countries"] = cleaned

                out_cfg = bench_file.with_stem(
                    f"{bench_file.stem}_{clusters}_{t_res}_{horizon}"
                )
                save_yaml(out_cfg, final_cfg)

                # We still generate MPS/LP for the **first stage only**
                mps_path = (
                    pathlib.Path(args.output_dir)
                    / out_cfg.with_suffix(f".{args.file_extension.lstrip('.')}")
                        .name
                )

                os.environ["ONLY_GENERATE_PROBLEM_FILE"] = str(mps_path)

                generate_benchmark(
                    args.benchmark_name,
                    out_cfg,
                    args.dry_run,
                    clusters,
                    t_res,
                )

                out_cfg.unlink(missing_ok=True)

