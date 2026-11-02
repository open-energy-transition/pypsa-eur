# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Solves a datarecord and an input scenario patch of it with PyPSA's own math.

Reads the multi-period record written by `build_datarecord` and stores it as
the root of a layered datarecord. It builds the network with
`datarecord.tools.pypsa.PyPSA`, optimises it without any extra functionality
and adds the results as a layer over the inputs. A second branch holds a
scenario patch of the inputs (`PATCH`) and the results of solving it. The
revision ids are written to a JSON file beside the store:

    inputs ─┬─ results
            └─ patch ── patch_results

The `solve` wildcard picks how the investment periods are solved:

- ``perfect``: all periods in one optimisation.
- ``myopic``: one period after the other, each starting from the capacities
  the previous one built.
- ``overnight<year>``: the one period ``<year>``, from the existing assets
  alone.

Myopic and overnight results are written from a `pypsa.NetworkCollection` of
the period networks, so every results layer has the shape of a
perfect-foresight one.

The script needs the `datarecord` pixi environment:
`pixi run -e datarecord snakemake -call solve_datarecord`.
"""

import json
import logging
from pathlib import Path

import pandas as pd
import pypsa
from datarecord import RecordLike, WorkingRecord, Workspace
from datarecord.store import Revision
from datarecord.tools.pypsa import PyPSA
from duckdb import SQLExpression

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

PATCH: list[tuple[str, list[str] | None, float]] = [
    ("capital_cost", ["offwind-ac", "offwind-dc", "offwind-float"], 1.3),
    ("p_nom_max", ["onwind"], 0.5),
    ("nom_max_onwind", None, 0.5),
    ("p_max_pu", ["solar", "solar-hsat"], 0.9),
]


def scenario_patch(record: RecordLike) -> WorkingRecord:
    """
    Scale each attribute in `PATCH`.

    An entry with carriers scales the attribute on the entities of those
    carriers; one without scales every per-period bus limit
    ``<attribute>_<period>`` on every bus.
    """
    w = WorkingRecord(record)
    entity = record.properties["entity"].select("entity", "carrier")
    for attribute, carriers, factor in PATCH:
        names = (
            [attribute]
            if carriers
            else [a for a in record.attributes if a.startswith(f"{attribute}_")]
        )
        if not names:
            raise KeyError(f"The record has no attribute {attribute}_<period>.")
        for name in names:
            rows = record.attributes[name]
            if carriers:
                rows = rows.join(entity, "entity").filter(
                    f"carrier IN ({', '.join(repr(c) for c in carriers)})"
                )
            keys = [
                c
                for c in rows.columns
                if c not in ("attribute", "breakpoint", "value", "carrier")
            ]
            w.set(
                name,
                rows.select(*keys, SQLExpression(f"value * {factor}").alias("value")),
            )
    return w


def period_network(
    n: pypsa.Network, period: int, previous: pypsa.Network | None
) -> pypsa.Network:
    """
    Cut the network of one investment period out of the multi-period `n`.

    It keeps the assets active in `period` and the global constraints of
    `period` or of no period. An extendable asset built before `period` takes
    its capacity from `previous`, the solved network of the period before,
    and is no longer extendable; without `previous` (overnight) a build-year
    cohort ``<name>-<build_year>`` built before `period` is dropped, while an
    asset of no cohort, such as a line, stays extendable.
    """
    snapshots = n.snapshots[n.snapshots.get_level_values("period") == period]
    m = n.copy(snapshots=snapshots, investment_periods=[period])
    for c in m.components:
        if c.static.empty or "build_year" not in c.static:
            continue
        attrs = c._operational_attrs
        static = c.static
        earlier = static[attrs["nom_extendable"]] & (static.build_year < period)
        if previous is None:
            cohort = pd.Series(
                [
                    name.endswith(f"-{year}")
                    for name, year in zip(static.index, static.build_year)
                ],
                index=static.index,
            )
            drop = ~m.get_active_assets(c.name, period) | (earlier & cohort)
            m.remove(c.name, static.index[drop])
            continue
        m.remove(c.name, static.index[~m.get_active_assets(c.name, period)])
        fixed = c.static.index[earlier.reindex(c.static.index)]
        built = previous.components[c.name].static.loc[fixed, f"{attrs['nom']}_opt"]
        c.static.loc[fixed, attrs["nom"]] = built
        c.static.loc[fixed, attrs["nom_extendable"]] = False
    gc = m.global_constraints
    other = gc.investment_period.notna() & (gc.investment_period != period)
    m.remove("GlobalConstraint", gc.index[other])
    return m


def optimize(n: pypsa.Network, solving: dict, log_fn: str) -> None:
    """Optimise `n` with PyPSA's own math and the configured solver."""
    solver_options = solving["solver"]["options"]
    status, condition = n.optimize(
        solver_name=solving["solver"]["name"],
        solver_options=solving["solver_options"][solver_options]
        if solver_options
        else {},
        transmission_losses=solving["options"]["transmission_losses"],
        multi_investment_periods=True,
        log_fn=log_fn,
    )
    if status != "ok":
        raise RuntimeError(f"Solving failed with status {status} ({condition}).")


def solve(inputs: Revision, mode: str, solving: dict, log_fn: str) -> Revision:
    """Solve the record of `inputs` by `mode` and add its results as a child."""
    missing = PyPSA.verify(inputs.record)
    if missing:
        raise RuntimeError(missing.describe())
    n = PyPSA.build(inputs.record)
    if mode == "perfect":
        optimize(n, solving, log_fn)
        solved: pypsa.Network | pypsa.NetworkCollection = n
        collection_dim = None
        objective = n.objective
    else:
        periods = list(n.investment_periods) if mode == "myopic" else [int(mode[9:])]
        members: list[pypsa.Network] = []
        for period in periods:
            previous = members[-1] if mode == "myopic" and members else None
            m = period_network(n, period, previous)
            optimize(m, solving, log_fn)
            members.append(m)
        solved = pypsa.NetworkCollection(
            members, index=pd.Index(periods, name="period")
        )
        collection_dim = "period"
        objective = sum(m.objective for m in members)
    results = PyPSA.to_datarecord(
        solved,
        inputs.store.workspace.con,
        results_only=True,
        base=inputs.record,
        collection_dim=collection_dim,
    )
    logger.info(f"{mode}: objective {objective:.6e}.")
    return inputs.child(results)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_datarecord",
            configfiles="config/test/config.preaggregated.yaml",
            horizon=2050,
            solve="perfect",
        )
    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)

    mode = snakemake.wildcards.solve
    solving = snakemake.params.solving
    ws = Workspace()
    store = ws.store(str(Path(snakemake.output.store).resolve()))
    inputs = store.root(ws.at(snakemake.input.record))

    results = solve(inputs, mode, solving, snakemake.log.solver)
    patch = inputs.child(scenario_patch(inputs.record).edits)
    patch_results = solve(patch, mode, solving, snakemake.log.solver_patch)

    revisions = {
        "inputs": inputs,
        "results": results,
        "patch": patch,
        "patch_results": patch_results,
    }
    Path(snakemake.output.revisions).write_text(
        json.dumps({name: str(rev.id) for name, rev in revisions.items()}, indent=2)
    )
    logger.info(f"Wrote {snakemake.output.store}.")
