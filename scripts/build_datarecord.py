# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Exports the network handed to the solver as a datarecord with its mathspec.

Reads the composed network of the last planning horizon (for perfect foresight
the stacked multi-period network) and applies the same preparation that
`solve_network` applies before it builds the model. It writes the network as a
datarecord with `datarecord.tools.pypsa.PyPSA`, together with the data of the
constraints that PyPSA's model does not hold: the global constraints of the
network and the custom constraints of `solve_network`
(`scripts.lib.mathspec_constraints`). It copies the base mathspec spec that
datarecord ships and the static math of those constraints
(`scripts/lib/mathspec/`) into `math/`, merges the three files, loads them
with mathspec and checks that every name the merged spec binds to is in the
written record.

The script needs the `datarecord` pixi environment:
`pixi run -e datarecord snakemake -call build_datarecord`.
"""

import logging
import shutil
from pathlib import Path

import mathspec
import numpy as np
import pypsa
from datarecord import Record, combine, connect, write_record
from datarecord.tools.pypsa import MATHSPEC, PyPSA

from scripts._helpers import configure_logging, set_scenario_config
from scripts.lib.mathspec_constraints import SPECS, constraint_data, extension
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


def drop_port_voltages(n: pypsa.Network) -> None:
    """
    Drop the ``voltage_bus0`` and ``voltage_bus1`` columns of a branch.

    They repeat ``v_nom`` of the bus at each port, which the record carries on
    the bus. A branch whose port voltage differs from ``v_nom`` of its bus raises.
    """
    for c in n.components[list(n.branch_components)]:
        for port in c.ports:
            column = f"voltage_bus{port}"
            if column not in c.static:
                continue
            bus_v_nom = c.static[f"bus{port}"].map(n.buses.v_nom)
            clash = c.static.index[c.static[column] != bus_v_nom]
            if not clash.empty:
                raise ValueError(
                    f"{c.name} components {list(clash)} carry a {column} that differs "
                    "from v_nom of their bus."
                )
            c.static.drop(columns=column, inplace=True)


def prefix_shared_entity_columns(n: pypsa.Network) -> None:
    """
    Prefix a custom component column with ``entity_`` where a bus column has its name.

    A record has one flat attribute namespace, so a custom component column
    whose name a bus column also takes (``location``, ``unit``, ``v_nom``, or
    ``tags`` from OSM) cannot be told apart from it. The bus column keeps its
    name, as the bus registry fixes most of them.
    """
    bus_columns = set(n.buses.columns) | set(n.components.buses.defaults.index)
    for c in n.components:
        if c.name in ("Bus", "Carrier", "GlobalConstraint", "Shape"):
            continue
        shared = [
            column
            for column in c.static.columns
            if column in bus_columns and column not in c.defaults.index
        ]
        c.static.rename(
            columns={column: f"entity_{column}" for column in shared}, inplace=True
        )


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


def add_zero_gas_usage(n: pypsa.Network) -> None:
    """
    Give every carrier a `gas_usage` of zero where no gas limit set one.

    The static math of primary energy limits reads `gas_usage` beside
    `co2_emissions`, so the record carries it for every network.
    """
    if "gas_usage" not in n.carriers:
        n.carriers["gas_usage"] = 0.0


def export_datarecord(
    n: pypsa.Network,
    config: dict,
    custom_extra_functionality: str | None,
    output: Path,
) -> None:
    """Write the record with the constraint data, copy the static math, then check them."""
    base = PyPSA.to_datarecord(n)
    data = constraint_data(n, config, custom_extra_functionality)
    record = combine(base, extension(data, base.schema))
    con = connect()
    write_record(None, record, con, uri=str(output / "record"))
    math = output / "math"
    math.mkdir(parents=True)
    for source in (MATHSPEC, *SPECS):
        shutil.copy(source, math / source.name)
    spec = mathspec.merge([math / source.name for source in (MATHSPEC, *SPECS)])
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
    drop_port_voltages(n)
    rename_shared_one_port_names(n)
    prefix_shared_entity_columns(n)
    add_zero_gas_usage(n)
    export_datarecord(
        n,
        snakemake.config,
        snakemake.params.custom_extra_functionality or None,
        Path(snakemake.output.directory),
    )
