<!-- SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur> -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Composing and Solving Networks

One rule assembles the network for each planning horizon from the clustered
electricity network and the prepared sector data, and one rule solves it. Both
run the same way for electricity-only and sector-coupled studies and for all
[foresight](../design/foresight.md) modes; the foresight mode decides which
inputs the composition step needs. Dispatch-only analyses on an already solved
network fix the expanded capacities and re-solve operation, optionally in a
rolling horizon.

## Temporal aggregation

{{ rules("time_aggregation") }}

## Composing

{{ rules("compose_network") }}

## Solving

{{ rules("solve_network", "solve_operations_network") }}

## Exporting and solving a datarecord

`build_datarecord` writes the network that the solver gets as a datarecord.
`solve_datarecord` builds a network from that record, solves it with PyPSA's
own math, without the custom constraints of `solve_network`, and writes the
solved network as a datarecord and as netCDF. Both need the pixi environment
`datarecord`. `config/test/config.preaggregated.yaml` runs them for five
Belgian regions over 2013 in 365 time segments, with perfect foresight over
2030, 2040 and 2050 and every power plant aggregated per region and carrier:
`pixi run -e datarecord snakemake -c4 --configfile config/test/config.preaggregated.yaml -- results/test-elec-preaggregated-perfect/networks/solved_datarecord_2050.nc`.

{{ rules("build_datarecord", "solve_datarecord") }}

## Library modules

These modules are not rules of their own. Their functions are imported and
called by [compose_network][] in the order the foresight mode requires.

### `add_electricity`

::: add_electricity

### `prepare_network`

::: prepare_network

### `prepare_sector_network`

::: prepare_sector_network

### `add_existing_baseyear`

::: add_existing_baseyear

### `add_brownfield`

::: add_brownfield

### `prepare_perfect_foresight`

::: prepare_perfect_foresight
