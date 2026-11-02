# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Combines the layered datarecords of several solve modes into one tree.

Each store written by `solve_datarecord` holds the inputs, a scenario patch
and the results of one solve mode over each. This script copies those layers
into one store, the inputs and the patch once, and the results of every mode
as a branch of each:

    inputs ─┬─ <mode> ...
            └─ patch ── <mode> ...

The revision ids are written to a JSON file beside the store, the results
under ``<mode>`` and ``patch/<mode>``.

The script needs the `datarecord` pixi environment.
"""

import json
import logging
from pathlib import Path
from uuid import UUID

from datarecord import Workspace

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "combine_datarecords",
            configfiles="config/test/config.preaggregated.yaml",
            horizon=2050,
        )
    configure_logging(snakemake)  # pylint: disable=E0606
    set_scenario_config(snakemake)

    ws = Workspace()
    store = ws.store(str(Path(snakemake.output.store).resolve()))
    inputs = store.root(ws.at(snakemake.input.record))
    revisions = {"inputs": inputs}

    sources = {
        Path(revisions_fn).parent.name: (
            ws.store(str(Path(store_dir).resolve())),
            json.loads(Path(revisions_fn).read_text()),
        )
        for store_dir, revisions_fn in zip(
            snakemake.input.stores, snakemake.input.revisions
        )
    }
    first, first_ids = next(iter(sources.values()))
    patch = inputs.child(first.layer(UUID(first_ids["patch"])))
    revisions["patch"] = patch
    for mode, (source, ids) in sources.items():
        revisions[mode] = inputs.child(source.layer(UUID(ids["results"])))
        revisions[f"patch/{mode}"] = patch.child(
            source.layer(UUID(ids["patch_results"]))
        )

    Path(snakemake.output.revisions).write_text(
        json.dumps({name: str(rev.id) for name, rev in revisions.items()}, indent=2)
    )
    logger.info(f"Combined {', '.join(sources)} into {snakemake.output.store}.")
