# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Exports the network handed to the solver as a datarecord with its mathspec.

Reads the composed network of the last planning horizon (for perfect foresight
the stacked multi-period network) and applies the same preparation that
`solve_network` applies before it builds the model. It writes the network as a
datarecord with `datarecord.tools.pypsa.PyPSA`, copies the base mathspec spec
that datarecord ships, and writes two fragments next to it: the global
constraints of the network and the custom constraints of `solve_network`.
The script merges the three files, loads them with mathspec and checks that
every name the merged spec binds to is in the written record.

The script needs the `datarecord` pixi environment:
`pixi run -e datarecord snakemake -call build_datarecord`.
"""

import logging
import shutil
from pathlib import Path

import mathspec
import numpy as np
import pypsa
import yaml
from datarecord import Record, connect, write_record
from datarecord.tools.pypsa import MATHSPEC, PyPSA

from scripts._helpers import configure_logging, set_scenario_config
from scripts.lib.mathspec_fragments import extra_functionality, global_constraints
from scripts.solve_network import prepare_network

logger = logging.getLogger(__name__)


def drop_bus_country_duplicates(n: pypsa.Network) -> None:
    """
    Drop a component `country` column that repeats the country of its bus.

    The bus country is the source of truth, and a record holds one attribute
    `country`. A component whose `country` differs from the country of its
    bus (`bus`, or `bus0` for a branch) raises.
    """
    for c in n.components:
        if c.name == "Bus" or "country" not in c.static:
            continue
        port = "bus" if "bus" in c.static else "bus0"
        bus_country = c.static[port].map(n.buses.country)
        clash = c.static.index[c.static.country != bus_country]
        if not clash.empty:
            raise ValueError(
                f"{c.name} components {list(clash)} carry a country that differs "
                "from the country of their bus."
            )
        c.static.drop(columns="country", inplace=True)


def prefix_shared_bus_columns(n: pypsa.Network) -> None:
    """
    Prefix a custom bus column with ``bus_`` where a component has it too.

    A record has one flat attribute namespace, and datarecord refuses a custom
    bus column whose name a custom component column also takes (``tags``,
    ``under_construction`` from OSM). This applies the prefix that datarecord
    applies to registry clashes such as ``bus_x``.
    """
    registry = set(n.components.buses.defaults.index)
    entity_columns = {
        column
        for c in n.components
        if c.name not in ("Bus", "Carrier", "GlobalConstraint", "Shape")
        for column in c.static.columns
    }
    shared = [
        column
        for column in n.buses.columns
        if column not in registry and column in entity_columns
    ]
    n.buses.rename(columns={column: f"bus_{column}" for column in shared}, inplace=True)


def rename_shared_one_port_names(n: pypsa.Network) -> None:
    """
    Give a one-port component a name of its own where a branch shares it.

    A record names an entity once across all component types; PyPSA scopes
    names per type, and PyPSA-Eur gives a load and the link that supplies it
    one name (``EU shipping oil``). The one-port takes the suffix
    `` (<type>)``; branch names stay, as the custom constraints read them.
    PyPSA's ``rename_component_names`` fails on a generator here, as it looks
    for a ``generator`` column in every component type.
    """
    branches = {
        name for c in n.components[list(n.branch_components)] for name in c.static.index
    }
    for c in n.components[["Generator", "Load", "Store", "StorageUnit"]]:
        shared = c.static.index.intersection(list(branches))
        mapping = {name: f"{name} ({c.name})" for name in shared}
        c.static.rename(index=mapping, inplace=True)
        for frame in c.dynamic.values():
            frame.rename(columns=mapping, inplace=True)
        if c.name == "Generator":
            n.buses["generator"] = n.buses.generator.replace(mapping)


def export_datarecord(
    n: pypsa.Network,
    config: dict,
    custom_extra_functionality: str | None,
    output: Path,
) -> None:
    """Write the record, the base spec and both fragments, then check them."""
    con = connect()
    write_record(None, PyPSA.to_datarecord(n), con, uri=str(output / "record"))
    math = output / "math"
    math.mkdir(parents=True)
    shutil.copy(MATHSPEC, math / MATHSPEC.name)
    base = yaml.safe_load(MATHSPEC.read_text())
    fragments = {
        "global_constraints.yaml": global_constraints(n, base),
        "extra_functionality.yaml": extra_functionality(
            n, config, custom_extra_functionality, base
        ),
    }
    for name, fragment in fragments.items():
        (math / name).write_text(yaml.safe_dump(fragment, sort_keys=False))
    spec = mathspec.merge([math / MATHSPEC.name, *(math / name for name in fragments)])
    PyPSA.check_spec(spec, Record.at(str(output / "record"), con))
    logger.info(
        f"Wrote the datarecord and a spec of {len(spec.constraints)} constraints "
        f"to {output}."
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_datarecord",
            configfiles="config/test/config.keep_assets.yaml",
        )
    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)

    solve_opts = snakemake.params.solving["options"]
    np.random.seed(solve_opts.get("seed", 123))

    n = pypsa.Network(snakemake.input.network)
    prepare_network(
        n,
        config=snakemake.config,
        solve_opts=solve_opts,
        foresight=snakemake.params.foresight,
        planning_horizons=str(snakemake.params.planning_horizons[-1]),
        co2_sequestration_potential=snakemake.params.co2_sequestration_potential,
        limit_max_growth=snakemake.params.sector["limit_max_growth"],
        rolling_horizon=solve_opts["rolling_horizon"],
    )
    drop_bus_country_duplicates(n)
    rename_shared_one_port_names(n)
    prefix_shared_bus_columns(n)
    export_datarecord(
        n,
        snakemake.config,
        snakemake.params.custom_extra_functionality or None,
        Path(snakemake.output.directory),
    )
