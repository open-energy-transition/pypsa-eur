# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Solves optimal operation and capacity for a network with the option to
iteratively optimize while updating line reactances.

This script is used for optimizing the electrical network as well as the
sector coupled network.

Description
-----------

Total annual system costs are minimised with PyPSA. The full formulation of the
linear optimal power flow (plus investment planning
is provided in the
`documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#linear-optimal-power-flow>`_.

The optimization is based on the :func:`network.optimize` function.
Additionally, some extra constraints specified in :mod:`solve_network` are added.

.. note::

    The rules ``solve_elec_networks`` and ``solve_sector_networks`` run
    the workflow for all scenarios in the configuration file (``scenario:``)
    based on the rule :mod:`solve_network`.
"""

import importlib
import logging
import os
import re
import sys
from functools import partial
from typing import Any

import country_converter as coco
import linopy
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
import yaml
from _benchmark import memory_logger
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
    #    create_tuples,
)
from add_electricity import add_missing_carriers, load_costs
from prepare_sector_network import get
from pypsa.descriptors import get_activity_mask
from pypsa.descriptors import get_switchable_as_dense as get_as_dense

cc = coco.CountryConverter()

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


class ObjectiveValueError(Exception):
    pass


def add_land_use_constraint_perfect(n: pypsa.Network) -> None:
    """
    Add global constraints for tech capacity limit.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance

    Returns
    -------
    pypsa.Network
        Network with added land use constraints
    """
    logger.info("Add land-use constraint for perfect foresight")

    def compress_series(s):
        def process_group(group):
            if group.nunique() == 1:
                return pd.Series(group.iloc[0], index=[None])
            else:
                return group

        return s.groupby(level=[0, 1]).apply(process_group)

    def new_index_name(t):
        # Convert all elements to string and filter out None values
        parts = [str(x) for x in t if x is not None]
        # Join with space, but use a dash for the last item if not None
        return " ".join(parts[:2]) + (f"-{parts[-1]}" if len(parts) > 2 else "")

    def check_p_min_p_max(p_nom_max):
        p_nom_min = n.generators[ext_i].groupby(grouper).sum().p_nom_min
        p_nom_min = p_nom_min.reindex(p_nom_max.index)
        check = (
            p_nom_min.groupby(level=[0, 1]).sum()
            > p_nom_max.groupby(level=[0, 1]).min()
        )
        if check.sum():
            logger.warning(
                f"summed p_min_pu values at node larger than technical potential {check[check].index}"
            )

    grouper = [n.generators.carrier, n.generators.bus, n.generators.build_year]
    ext_i = n.generators.p_nom_extendable
    # get technical limit per node and investment period
    p_nom_max = n.generators[ext_i].groupby(grouper).min().p_nom_max
    # drop carriers without tech limit
    p_nom_max = p_nom_max[~p_nom_max.isin([np.inf, np.nan])]
    # carrier
    carriers = p_nom_max.index.get_level_values(0).unique()
    gen_i = n.generators[(n.generators.carrier.isin(carriers)) & (ext_i)].index
    n.generators.loc[gen_i, "p_nom_min"] = 0
    # check minimum capacities
    check_p_min_p_max(p_nom_max)
    # drop multi entries in case p_nom_max stays constant in different periods
    # p_nom_max = compress_series(p_nom_max)
    # adjust name to fit syntax of nominal constraint per bus
    df = p_nom_max.reset_index()
    df["name"] = df.apply(
        lambda row: f"nom_max_{row['carrier']}"
        + (f"_{row['build_year']}" if row["build_year"] is not None else ""),
        axis=1,
    )

    for name in df.name.unique():
        df_carrier = df[df.name == name]
        bus = df_carrier.bus
        n.buses.loc[bus, name] = df_carrier.p_nom_max.values


def add_land_use_constraint(n: pypsa.Network, planning_horizons: str) -> None:
    """
    Add land use constraints for renewable energy potential.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    planning_horizons : str
        The planning horizon year as string

    Returns
    -------
    pypsa.Network
        Modified PyPSA network with constraints added
    """
    # warning: this will miss existing offwind which is not classed AC-DC and has carrier 'offwind'

    for carrier in [
        "solar",
        "solar rooftop",
        "solar-hsat",
        "onwind",
        "offwind-ac",
        "offwind-dc",
        "offwind-float",
    ]:
        ext_i = (n.generators.carrier == carrier) & ~n.generators.p_nom_extendable
        existing = (
            n.generators.loc[ext_i, "p_nom"]
            .groupby(n.generators.bus.map(n.buses.location))
            .sum()
        )
        existing.index += f" {carrier}-{planning_horizons}"
        n.generators.loc[existing.index, "p_nom_max"] -= existing

    # check if existing capacities are larger than technical potential
    existing_large = n.generators[
        n.generators["p_nom_min"] > n.generators["p_nom_max"]
    ].index
    if len(existing_large):
        logger.warning(
            f"Existing capacities larger than technical potential for {existing_large},\
                        adjust technical potential to existing capacities"
        )
        n.generators.loc[existing_large, "p_nom_max"] = n.generators.loc[
            existing_large, "p_nom_min"
        ]

    n.generators["p_nom_max"] = n.generators["p_nom_max"].clip(lower=0)


def add_solar_potential_constraints(n: pypsa.Network, config: dict) -> None:
    """
    Add constraint to make sure the sum capacity of all solar technologies (fixed, tracking, ets. ) is below the region potential.

    Example:
    ES1 0: total solar potential is 10 GW, meaning:
           solar potential : 10 GW
           solar-hsat potential : 8 GW (solar with single axis tracking is assumed to have higher land use)
    The constraint ensures that:
           solar_p_nom + solar_hsat_p_nom * 1.13 <= 10 GW
    """
    land_use_factors = {
        "solar-hsat": config["renewable"]["solar"]["capacity_per_sqkm"]
        / config["renewable"]["solar-hsat"]["capacity_per_sqkm"],
    }
    rename = {"Generator-ext": "Generator"}

    solar_carriers = ["solar", "solar-hsat"]
    solar = n.generators[
        n.generators.carrier.isin(solar_carriers) & n.generators.p_nom_extendable
    ].index

    solar_today = n.generators[
        (n.generators.carrier == "solar") & (n.generators.p_nom_extendable)
    ].index
    solar_hsat = n.generators[(n.generators.carrier == "solar-hsat")].index

    if solar.empty:
        return

    land_use = pd.DataFrame(1, index=solar, columns=["land_use_factor"])
    for carrier, factor in land_use_factors.items():
        land_use = land_use.apply(
            lambda x: (x * factor) if carrier in x.name else x, axis=1
        )

    location = pd.Series(n.buses.index, index=n.buses.index)
    ggrouper = n.generators.loc[solar].bus
    rhs = (
        n.generators.loc[solar_today, "p_nom_max"]
        .groupby(n.generators.loc[solar_today].bus.map(location))
        .sum()
        - n.generators.loc[solar_hsat, "p_nom"]
        .groupby(n.generators.loc[solar_hsat].bus.map(location))
        .sum()
        * land_use_factors["solar-hsat"]
    ).clip(lower=0)

    lhs = (
        (n.model["Generator-p_nom"].rename(rename).loc[solar] * land_use.squeeze())
        .groupby(ggrouper)
        .sum()
    )

    logger.info("Adding solar potential constraint.")
    n.model.add_constraints(lhs <= rhs, name="solar_potential")


def add_co2_sequestration_limit(
    n: pypsa.Network,
    limit_dict: dict[str, float],
    planning_horizons: str | None,
) -> None:
    """
    Add a global constraint on the amount of Mt CO2 that can be sequestered.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    limit_dict : dict[str, float]
        CO2 sequestration potential limit constraints by year.
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight
    """

    if not n.investment_periods.empty:
        periods = n.investment_periods
        limit = pd.Series(
            {
                f"co2_sequestration_limit-{period}": limit_dict.get(period, 200)
                for period in periods
            }
        )
        names = limit.index
    else:
        limit = get(limit_dict, int(planning_horizons))
        periods = [np.nan]
        names = pd.Index(["co2_sequestration_limit"])

    n.add(
        "GlobalConstraint",
        names,
        sense=">=",
        constant=-limit * 1e6,
        type="operational_limit",
        carrier_attribute="co2 sequestered",
        investment_period=periods,
    )


def add_carbon_constraint(n: pypsa.Network, snapshots: pd.DatetimeIndex) -> None:
    glcs = n.global_constraints.query('type == "co2_atmosphere"')
    if glcs.empty:
        return
    for name, glc in glcs.iterrows():
        carattr = glc.carrier_attribute
        emissions = n.carriers.query(f"{carattr} != 0")[carattr]

        if emissions.empty:
            continue

        # stores
        bus_carrier = n.stores.bus.map(n.buses.carrier)
        stores = n.stores[bus_carrier.isin(emissions.index) & ~n.stores.e_cyclic]
        if not stores.empty:
            last = n.snapshot_weightings.reset_index().groupby("period").last()
            last_i = last.set_index([last.index, last.timestep]).index
            final_e = n.model["Store-e"].loc[last_i, stores.index]
            time_valid = int(glc.loc["investment_period"])
            time_i = pd.IndexSlice[time_valid, :]
            lhs = final_e.loc[time_i, :] - final_e.shift(snapshot=1).loc[time_i, :]

            rhs = glc.constant
            n.model.add_constraints(lhs <= rhs, name=f"GlobalConstraint-{name}")


def add_carbon_budget_constraint(n: pypsa.Network, snapshots: pd.DatetimeIndex) -> None:
    glcs = n.global_constraints.query('type == "Co2Budget"')
    if glcs.empty:
        return
    for name, glc in glcs.iterrows():
        carattr = glc.carrier_attribute
        emissions = n.carriers.query(f"{carattr} != 0")[carattr]

        if emissions.empty:
            continue

        # stores
        bus_carrier = n.stores.bus.map(n.buses.carrier)
        stores = n.stores[bus_carrier.isin(emissions.index) & ~n.stores.e_cyclic]
        if not stores.empty:
            last = n.snapshot_weightings.reset_index().groupby("period").last()
            last_i = last.set_index([last.index, last.timestep]).index
            final_e = n.model["Store-e"].loc[last_i, stores.index]
            time_valid = int(glc.loc["investment_period"])
            time_i = pd.IndexSlice[time_valid, :]
            weighting = n.investment_period_weightings.loc[time_valid, "years"]
            lhs = final_e.loc[time_i, :] * weighting

            rhs = glc.constant
            n.model.add_constraints(lhs <= rhs, name=f"GlobalConstraint-{name}")


def add_max_growth(n: pypsa.Network, opts: dict) -> None:
    """
    Add maximum growth rates for different carriers.
    """

    # take maximum yearly difference between investment periods since historic growth is per year
    factor = n.investment_period_weightings.years.max() * opts["factor"]
    for carrier in opts["max_growth"].keys():
        max_per_period = opts["max_growth"][carrier] * factor
        logger.info(
            f"set maximum growth rate per investment period of {carrier} to {max_per_period} GW."
        )
        n.carriers.loc[carrier, "max_growth"] = max_per_period * 1e3

    for carrier in opts["max_relative_growth"].keys():
        max_r_per_period = opts["max_relative_growth"][carrier]
        logger.info(
            f"set maximum relative growth per investment period of {carrier} to {max_r_per_period}."
        )
        n.carriers.loc[carrier, "max_relative_growth"] = max_r_per_period


def add_retrofit_gas_boiler_constraint(
    n: pypsa.Network, snapshots: pd.DatetimeIndex
) -> None:
    """
    Allow retrofitting of existing gas boilers to H2 boilers and impose load-following must-run condition on existing gas boilers.
    Modifies the network in place, no return value.

    n : pypsa.Network
        The PyPSA network to be modified
    snapshots : pd.DatetimeIndex
        The snapshots of the network
    """
    c = "Link"
    logger.info("Add constraint for retrofitting gas boilers to H2 boilers.")
    # existing gas boilers
    mask = n.links.carrier.str.contains("gas boiler") & ~n.links.p_nom_extendable
    gas_i = n.links[mask].index
    mask = n.links.carrier.str.contains("retrofitted H2 boiler")
    h2_i = n.links[mask].index

    n.links.loc[gas_i, "p_nom_extendable"] = True
    p_nom = n.links.loc[gas_i, "p_nom"]
    n.links.loc[gas_i, "p_nom"] = 0

    # heat profile
    cols = n.loads_t.p_set.columns[
        n.loads_t.p_set.columns.str.contains("heat")
        & ~n.loads_t.p_set.columns.str.contains("industry")
        & ~n.loads_t.p_set.columns.str.contains("agriculture")
    ]
    profile = n.loads_t.p_set[cols].div(
        n.loads_t.p_set[cols].groupby(level=0).max(), level=0
    )
    # to deal if max value is zero
    profile.fillna(0, inplace=True)
    profile.rename(columns=n.loads.bus.to_dict(), inplace=True)
    profile = profile.reindex(columns=n.links.loc[gas_i, "bus1"])
    profile.columns = gas_i

    rhs = profile.mul(p_nom)

    dispatch = n.model["Link-p"]
    active = get_activity_mask(n, c, snapshots, gas_i)
    rhs = rhs[active]
    p_gas = dispatch.sel(Link=gas_i)
    p_h2 = dispatch.sel(Link=h2_i)

    lhs = p_gas + p_h2

    n.model.add_constraints(lhs == rhs, name="gas_retrofit")


def prepare_network(
    n: pypsa.Network,
    solve_opts: dict,
    foresight: str,
    planning_horizons: str | None,
    co2_sequestration_potential: dict[str, float],
    limit_max_growth: dict[str, Any] | None = None,
) -> None:
    """
    Prepare network with various constraints and modifications.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    solve_opts : Dict
        Dictionary of solving options containing clip_p_max_pu, load_shedding etc.
    foresight : str
        Planning foresight type ('myopic' or 'perfect')
    planning_horizons : str or None
        The current planning horizon year or None for perfect foresight
    co2_sequestration_potential : Dict[str, float]
        CO2 sequestration potential constraints by year

    Returns
    -------
    pypsa.Network
        Modified PyPSA network with added constraints
    """
    if "clip_p_max_pu" in solve_opts:
        for df in (
            n.generators_t.p_max_pu,
            n.generators_t.p_min_pu,
            n.links_t.p_max_pu,
            n.links_t.p_min_pu,
            n.storage_units_t.inflow,
        ):
            df.where(df > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    if load_shedding := solve_opts.get("load_shedding"):
        # intersect between macroeconomic and surveybased willingness to pay
        # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
        # TODO: retrieve color and nice name from config
        n.add("Carrier", "load", color="#dd2e23", nice_name="Load shedding")
        buses_i = n.buses.index
        if not np.isscalar(load_shedding):
            # TODO: do not scale via sign attribute (use Eur/MWh instead of Eur/kWh)
            load_shedding = 1e2  # Eur/kWh

        n.add(
            "Generator",
            buses_i,
            " load",
            bus=buses_i,
            carrier="load",
            sign=1e-3,  # Adjust sign to measure p and p_nom in kW instead of MW
            marginal_cost=load_shedding,  # Eur/kWh
            p_nom=1e9,  # kW
        )

    if solve_opts.get("curtailment_mode"):
        n.add("Carrier", "curtailment", color="#fedfed", nice_name="Curtailment")
        n.generators_t.p_min_pu = n.generators_t.p_max_pu
        buses_i = n.buses.query("carrier == 'AC'").index
        n.add(
            "Generator",
            buses_i,
            suffix=" curtailment",
            bus=buses_i,
            p_min_pu=-1,
            p_max_pu=0,
            marginal_cost=-0.1,
            carrier="curtailment",
            p_nom=1e6,
        )

    if solve_opts.get("noisy_costs"):
        for t in n.iterate_components():
            # if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if "marginal_cost" in t.df:
                t.df["marginal_cost"] += 1e-2 + 2e-3 * (
                    np.random.random(len(t.df)) - 0.5
                )

        for t in n.iterate_components(["Line", "Link"]):
            t.df["capital_cost"] += (
                1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)
            ) * t.df["length"]

    if solve_opts.get("nhours"):
        nhours = solve_opts["nhours"]
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760.0 / nhours

    if foresight == "myopic":
        add_land_use_constraint(n, planning_horizons)

    if foresight == "perfect":
        add_land_use_constraint_perfect(n)
        if limit_max_growth is not None and limit_max_growth["enable"]:
            add_max_growth(n, limit_max_growth)

    if n.stores.carrier.eq("co2 sequestered").any():
        limit_dict = co2_sequestration_potential
        add_co2_sequestration_limit(
            n, limit_dict=limit_dict, planning_horizons=planning_horizons
        )


def calculate_grid_score(
    n: pypsa.Network, include_techs: list, name: str, include_ci=False
) -> None:
    """
    Calculates the time-series grid supply score for each nodes, based on the share of energy generated by the technologies specified in include_techs.

    NOTE: This calculation reflects the generation-based score, not the consumption-based score.
    If the goal is to assess consumption, the score must be weighted by the score of imported electricity.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    include_techs : list
        Configuration dictionary containing solver settings
    name: str
        Name of the score (e.g. cfe, res)
    include_ci: bool
        Set if CI related generators and links are included in the score or not.

    Returns
    -------
    pypsa.Network
        Modified PyPSA network with added attribute of {name}_score in n.buses and n.buses_t
    """

    def get_values(n, df, df_t, bus_col, include_techs, include_ci=False):
        # Map low-voltage bus to main grid bus
        grid_buses = n.buses[n.buses.carrier == "AC"].index
        low_voltage_map = (
            n.links[
                (n.links.carrier == "electricity distribution grid")
                & ~n.links.bus0.isin(grid_buses)
            ]
            .set_index("bus0")["bus1"]
            .to_dict()
        )

        # Prepare and annotate the time series data
        df_t = df_t.T.copy()
        df_t = df_t.join(df[[bus_col, "carrier"]])
        df_t["bus"] = df_t[bus_col].map(low_voltage_map).fillna(df_t[bus_col])

        # Filter out grid specific carriers
        exclude_carriers = {"electricity distribution grid", "AC", "DC"}
        df_t = df_t[df_t["bus"].isin(grid_buses) & ~df_t.carrier.isin(exclude_carriers)]

        # Remove CI if include_ci is False
        if not include_ci:
            remain_index = df[df["ci"].isin([np.NaN, ""])].index
            df_t = df_t[df_t.index.isin(remain_index)]

        # Aggregate values for included technologies and all carriers
        df_t_clean = (
            df_t[df_t.carrier.isin(include_techs)].groupby("bus")[n.snapshots].sum().T
        )
        df_t_all = df_t.groupby("bus")[n.snapshots].sum().T

        return df_t_clean, df_t_all

    df_gen_clean, df_gen_total = get_values(
        n, n.generators, n.generators_t.p, "bus", include_techs, include_ci=include_ci
    )
    df_link_clean, df_link_total = get_values(
        n, n.links, n.links_t.p1, "bus1", include_techs, include_ci=include_ci
    )

    n.buses_t[f"{name}_p"] = (
        pd.concat([df_gen_clean, -df_link_clean], axis=1).T.groupby(level=0).sum().T
    )
    all_p = pd.concat([df_gen_total, -df_link_total], axis=1).T.groupby(level=0).sum().T

    n.buses_t[f"{name}_score"] = n.buses_t[f"{name}_p"] / all_p
    n.buses[f"{name}_score"] = n.buses_t[f"{name}_p"].sum() / all_p.sum()

    if n.buses_t[f"{name}_score"].empty:
        grid_buses = n.buses[n.buses.carrier == "AC"].index
        n.buses_t[f"{name}_score"] = pd.DataFrame(
            0, index=n.snapshots, columns=grid_buses
        )
        logger.info(f"{name}_score currently is empty")
    else:
        global_score = round(
            n.buses_t[f"{name}_p"].sum().sum() / all_p.sum().sum() * 100, 2
        )
        logger.info(f"The average {name}_score is: {global_score}%")


def add_CCL_constraints(
    n: pypsa.Network, config: dict, planning_horizons: str | None
) -> None:
    """
    Add CCL (country & carrier limit) constraint to the network.

    Add minimum and maximum levels of generator nominal capacity per carrier
    for individual countries. Opts and path for agg_p_nom_minmax.csv must be defined
    in config.yaml. Default file is available at data/agg_p_nom_minmax.csv.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    config : dict
        Configuration dictionary
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight

    Example
    -------
    scenario:
        opts: [Co2L-CCL-24h]
    electricity:
        agg_p_nom_limits: data/agg_p_nom_minmax.csv
    """

    assert planning_horizons is not None, (
        "add_CCL_constraints are not implemented for perfect foresight, yet"
    )

    agg_p_nom_minmax = pd.read_csv(
        config["solving"]["agg_p_nom_limits"]["file"], index_col=[0, 1], header=[0, 1]
    )[planning_horizons]
    logger.info("Adding generation capacity constraints per carrier and country")
    p_nom = n.model["Generator-p_nom"]

    gens = n.generators.query("p_nom_extendable").rename_axis(index="Generator-ext")
    if config["solving"]["agg_p_nom_limits"]["agg_offwind"]:
        rename_offwind = {
            "offwind-ac": "offwind-all",
            "offwind-dc": "offwind-all",
            "offwind": "offwind-all",
        }
        gens = gens.replace(rename_offwind)
    grouper = pd.concat([gens.bus.map(n.buses.country), gens.carrier], axis=1)
    lhs = p_nom.groupby(grouper).sum().rename(bus="country")

    if config["solving"]["agg_p_nom_limits"]["include_existing"]:
        gens_cst = n.generators.query("~p_nom_extendable").rename_axis(
            index="Generator-cst"
        )
        gens_cst = gens_cst[
            (gens_cst["build_year"] + gens_cst["lifetime"]) >= int(planning_horizons)
        ]
        if config["solving"]["agg_p_nom_limits"]["agg_offwind"]:
            gens_cst = gens_cst.replace(rename_offwind)
        rhs_cst = (
            pd.concat(
                [gens_cst.bus.map(n.buses.country), gens_cst[["carrier", "p_nom"]]],
                axis=1,
            )
            .groupby(["bus", "carrier"])
            .sum()
        )
        rhs_cst.index = rhs_cst.index.rename({"bus": "country"})
        rhs_min = agg_p_nom_minmax["min"].dropna()
        idx_min = rhs_min.index.join(rhs_cst.index, how="left")
        rhs_min = rhs_min.reindex(idx_min).fillna(0)
        rhs = (rhs_min - rhs_cst.reindex(idx_min).fillna(0).p_nom).dropna()
        rhs[rhs < 0] = 0
        minimum = xr.DataArray(rhs).rename(dim_0="group")
    else:
        minimum = xr.DataArray(agg_p_nom_minmax["min"].dropna()).rename(dim_0="group")

    index = minimum.indexes["group"].intersection(lhs.indexes["group"])
    if not index.empty:
        n.model.add_constraints(
            lhs.sel(group=index) >= minimum.loc[index], name="agg_p_nom_min"
        )

    if config["solving"]["agg_p_nom_limits"]["include_existing"]:
        rhs_max = agg_p_nom_minmax["max"].dropna()
        idx_max = rhs_max.index.join(rhs_cst.index, how="left")
        rhs_max = rhs_max.reindex(idx_max).fillna(0)
        rhs = (rhs_max - rhs_cst.reindex(idx_max).fillna(0).p_nom).dropna()
        rhs[rhs < 0] = 0
        maximum = xr.DataArray(rhs).rename(dim_0="group")
    else:
        maximum = xr.DataArray(agg_p_nom_minmax["max"].dropna()).rename(dim_0="group")

    index = maximum.indexes["group"].intersection(lhs.indexes["group"])
    if not index.empty:
        n.model.add_constraints(
            lhs.sel(group=index) <= maximum.loc[index], name="agg_p_nom_max"
        )


def add_EQ_constraints(n, o, scaling=1e-1):
    """
    Add equity constraints to the network.

    Currently this is only implemented for the electricity sector only.

    Opts must be specified in the config.yaml.

    Parameters
    ----------
    n : pypsa.Network
    o : str

    Example
    -------
    scenario:
        opts: [Co2L-EQ0.7-24h]

    Require each country or node to on average produce a minimal share
    of its total electricity consumption itself. Example: EQ0.7c demands each country
    to produce on average at least 70% of its consumption; EQ0.7 demands
    each node to produce on average at least 70% of its consumption.
    """
    # TODO: Generalize to cover myopic and other sectors?
    float_regex = r"[0-9]*\.?[0-9]+"
    level = float(re.findall(float_regex, o)[0])
    if o[-1] == "c":
        ggrouper = n.generators.bus.map(n.buses.country)
        lgrouper = n.loads.bus.map(n.buses.country)
        sgrouper = n.storage_units.bus.map(n.buses.country)
    else:
        ggrouper = n.generators.bus
        lgrouper = n.loads.bus
        sgrouper = n.storage_units.bus
    load = (
        n.snapshot_weightings.generators
        @ n.loads_t.p_set.groupby(lgrouper, axis=1).sum()
    )
    inflow = (
        n.snapshot_weightings.stores
        @ n.storage_units_t.inflow.groupby(sgrouper, axis=1).sum()
    )
    inflow = inflow.reindex(load.index).fillna(0.0)
    rhs = scaling * (level * load - inflow)
    p = n.model["Generator-p"]
    lhs_gen = (
        (p * (n.snapshot_weightings.generators * scaling))
        .groupby(ggrouper.to_xarray())
        .sum()
        .sum("snapshot")
    )
    # TODO: double check that this is really needed, why do have to subtract the spillage
    if not n.storage_units_t.inflow.empty:
        spillage = n.model["StorageUnit-spill"]
        lhs_spill = (
            (spillage * (-n.snapshot_weightings.stores * scaling))
            .groupby(sgrouper.to_xarray())
            .sum()
            .sum("snapshot")
        )
        lhs = lhs_gen + lhs_spill
    else:
        lhs = lhs_gen
    n.model.add_constraints(lhs >= rhs, name="equity_min")


def add_BAU_constraints(n: pypsa.Network, config: dict) -> None:
    """
    Add business-as-usual (BAU) constraints for minimum capacities.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network instance
    config : dict
        Configuration dictionary containing BAU minimum capacities
    """
    mincaps = pd.Series(config["electricity"]["BAU_mincapacities"])
    p_nom = n.model["Generator-p_nom"]
    ext_i = n.generators.query("p_nom_extendable")
    ext_carrier_i = xr.DataArray(ext_i.carrier.rename_axis("Generator-ext"))
    lhs = p_nom.groupby(ext_carrier_i).sum()
    rhs = mincaps[lhs.indexes["carrier"]].rename_axis("carrier")
    n.model.add_constraints(lhs >= rhs, name="bau_mincaps")


# TODO: think about removing or make per country
def add_SAFE_constraints(n, config):
    """
    Add a capacity reserve margin of a certain fraction above the peak demand.
    Renewable generators and storage do not contribute. Ignores network.

    Parameters
    ----------
        n : pypsa.Network
        config : dict

    Example
    -------
    config.yaml requires to specify opts:

    scenario:
        opts: [Co2L-SAFE-24h]
    electricity:
        SAFE_reservemargin: 0.1
    Which sets a reserve margin of 10% above the peak demand.
    """
    peakdemand = n.loads_t.p_set.sum(axis=1).max()
    margin = 1.0 + config["electricity"]["SAFE_reservemargin"]
    reserve_margin = peakdemand * margin
    conventional_carriers = config["electricity"]["conventional_carriers"]  # noqa: F841
    ext_gens_i = n.generators.query(
        "carrier in @conventional_carriers & p_nom_extendable"
    ).index
    p_nom = n.model["Generator-p_nom"].loc[ext_gens_i]
    lhs = p_nom.sum()
    exist_conv_caps = n.generators.query(
        "~p_nom_extendable & carrier in @conventional_carriers"
    ).p_nom.sum()
    rhs = reserve_margin - exist_conv_caps
    n.model.add_constraints(lhs >= rhs, name="safe_mintotalcap")


def add_operational_reserve_margin(n, sns, config):
    """
    Build reserve margin constraints based on the formulation given in
    https://genxproject.github.io/GenX/dev/core/#Reserves.

    Parameters
    ----------
        n : pypsa.Network
        sns: pd.DatetimeIndex
        config : dict

    Example:
    --------
    config.yaml requires to specify operational_reserve:
    operational_reserve: # like https://genxproject.github.io/GenX/dev/core/#Reserves
        activate: true
        epsilon_load: 0.02 # percentage of load at each snapshot
        epsilon_vres: 0.02 # percentage of VRES at each snapshot
        contingency: 400000 # MW
    """
    reserve_config = config["electricity"]["operational_reserve"]
    EPSILON_LOAD = reserve_config["epsilon_load"]
    EPSILON_VRES = reserve_config["epsilon_vres"]
    CONTINGENCY = reserve_config["contingency"]

    # Reserve Variables
    n.model.add_variables(
        0, np.inf, coords=[sns, n.generators.index], name="Generator-r"
    )
    reserve = n.model["Generator-r"]
    summed_reserve = reserve.sum("Generator")

    # Share of extendable renewable capacities
    ext_i = n.generators.query("p_nom_extendable").index
    vres_i = n.generators_t.p_max_pu.columns
    if not ext_i.empty and not vres_i.empty:
        capacity_factor = n.generators_t.p_max_pu[vres_i.intersection(ext_i)]
        p_nom_vres = (
            n.model["Generator-p_nom"]
            .loc[vres_i.intersection(ext_i)]
            .rename({"Generator-ext": "Generator"})
        )
        lhs = summed_reserve + (
            p_nom_vres * (-EPSILON_VRES * xr.DataArray(capacity_factor))
        ).sum("Generator")

        # Total demand per t
        demand = get_as_dense(n, "Load", "p_set").sum(axis=1)

        # VRES potential of non extendable generators
        capacity_factor = n.generators_t.p_max_pu[vres_i.difference(ext_i)]
        renewable_capacity = n.generators.p_nom[vres_i.difference(ext_i)]
        potential = (capacity_factor * renewable_capacity).sum(axis=1)

        # Right-hand-side
        rhs = EPSILON_LOAD * demand + EPSILON_VRES * potential + CONTINGENCY

        n.model.add_constraints(lhs >= rhs, name="reserve_margin")

    # additional constraint that capacity is not exceeded
    gen_i = n.generators.index
    ext_i = n.generators.query("p_nom_extendable").index
    fix_i = n.generators.query("not p_nom_extendable").index

    dispatch = n.model["Generator-p"]
    reserve = n.model["Generator-r"]

    capacity_variable = n.model["Generator-p_nom"].rename(
        {"Generator-ext": "Generator"}
    )
    capacity_fixed = n.generators.p_nom[fix_i]

    p_max_pu = get_as_dense(n, "Generator", "p_max_pu")

    lhs = dispatch + reserve - capacity_variable * xr.DataArray(p_max_pu[ext_i])

    rhs = (p_max_pu[fix_i] * capacity_fixed).reindex(columns=gen_i, fill_value=0)

    n.model.add_constraints(lhs <= rhs, name="Generator-p-reserve-upper")


def add_TES_energy_to_power_ratio_constraints(n: pypsa.Network) -> None:
    """
    Add TES constraints to the network.

    For each TES storage unit, enforce:
        Store-e_nom - etpr * Link-p_nom == 0

    Parameters
    ----------
    n : pypsa.Network
        A PyPSA network with TES and heating sectors enabled.

    Raises
    ------
    ValueError
        If no valid TES storage or charger links are found.
    RuntimeError
        If the TES storage and charger indices do not align.
    """
    indices_charger_p_nom_extendable = n.links.index[
        n.links.index.str.contains("water tanks charger|water pits charger")
        & n.links.p_nom_extendable
    ]
    indices_stores_e_nom_extendable = n.stores.index[
        n.stores.index.str.contains("water tanks|water pits")
        & n.stores.e_nom_extendable
    ]

    if indices_charger_p_nom_extendable.empty or indices_stores_e_nom_extendable.empty:
        raise ValueError(
            "No valid extendable charger links or stores found for TES energy to power constraints."
        )

    energy_to_power_ratio_values = n.links.loc[
        indices_charger_p_nom_extendable, "energy to power ratio"
    ].values

    linear_expr_list = []
    for charger, tes, energy_to_power_value in zip(
        indices_charger_p_nom_extendable,
        indices_stores_e_nom_extendable,
        energy_to_power_ratio_values,
    ):
        charger_var = n.model["Link-p_nom"].loc[charger]
        if not tes == charger.replace(" charger", ""):
            # e.g. "DE0 0 urban central water tanks charger-2050" -> "DE0 0 urban central water tanks-2050"
            raise RuntimeError(
                f"Charger {charger} and TES {tes} do not match. "
                "Ensure that the charger and TES are in the same location and refer to the same technology."
            )
        store_var = n.model["Store-e_nom"].loc[tes]
        linear_expr = store_var - energy_to_power_value * charger_var
        linear_expr_list.append(linear_expr)

    # Merge the individual expressions
    merged_expr = linopy.expressions.merge(
        linear_expr_list, dim="Store-ext, Link-ext", cls=type(linear_expr_list[0])
    )

    n.model.add_constraints(merged_expr == 0, name="TES_energy_to_power_ratio")


def add_TES_charger_ratio_constraints(n: pypsa.Network) -> None:
    """
    Add TES charger ratio constraints.

    For each TES unit, enforce:
        Link-p_nom(charger) - efficiency * Link-p_nom(discharger) == 0

    Parameters
    ----------
    n : pypsa.Network
        A PyPSA network with TES and heating sectors enabled.

    Raises
    ------
    ValueError
        If no valid TES discharger or charger links are found.
    RuntimeError
        If the charger and discharger indices do not align.
    """
    indices_charger_p_nom_extendable = n.links.index[
        n.links.index.str.contains("water tanks charger|water pits charger")
        & n.links.p_nom_extendable
    ]
    indices_discharger_p_nom_extendable = n.links.index[
        n.links.index.str.contains("water tanks discharger|water pits discharger")
        & n.links.p_nom_extendable
    ]

    if (
        indices_charger_p_nom_extendable.empty
        or indices_discharger_p_nom_extendable.empty
    ):
        raise ValueError(
            "No valid extendable TES discharger or charger links found for TES charger ratio constraints."
        )

    for charger, discharger in zip(
        indices_charger_p_nom_extendable, indices_discharger_p_nom_extendable
    ):
        if not charger.replace(" charger", " ") == discharger.replace(
            " discharger", " "
        ):
            # e.g. "DE0 0 urban central water tanks charger-2050" -> "DE0 0 urban central water tanks-2050"
            raise RuntimeError(
                f"Charger {charger} and discharger {discharger} do not match. "
                "Ensure that the charger and discharger are in the same location and refer to the same technology."
            )

    eff_discharger = n.links.efficiency[indices_discharger_p_nom_extendable].values
    lhs = (
        n.model["Link-p_nom"].loc[indices_charger_p_nom_extendable]
        - n.model["Link-p_nom"].loc[indices_discharger_p_nom_extendable]
        * eff_discharger
    )

    n.model.add_constraints(lhs == 0, name="TES_charger_ratio")


def add_battery_constraints(n):
    """
    Add constraint ensuring that charger = discharger, i.e.
    1 * charger_size - efficiency * discharger_size = 0
    """
    if not n.links.p_nom_extendable.any():
        return

    discharger_bool = n.links.index.str.contains("battery discharger")
    charger_bool = n.links.index.str.contains("battery charger")

    dischargers_ext = n.links[discharger_bool].query("p_nom_extendable").index
    chargers_ext = n.links[charger_bool].query("p_nom_extendable").index

    eff = n.links.efficiency[dischargers_ext].values
    lhs = (
        n.model["Link-p_nom"].loc[chargers_ext]
        - n.model["Link-p_nom"].loc[dischargers_ext] * eff
    )

    n.model.add_constraints(lhs == 0, name="Link-charger_ratio")


def add_lossy_bidirectional_link_constraints(n):
    if not n.links.p_nom_extendable.any() or not any(n.links.get("reversed", [])):
        return

    carriers = n.links.loc[n.links.reversed, "carrier"].unique()  # noqa: F841
    backwards = n.links.query(
        "carrier in @carriers and p_nom_extendable and reversed"
    ).index
    forwards = backwards.str.replace("-reversed", "")
    lhs = n.model["Link-p_nom"].loc[backwards]
    rhs = n.model["Link-p_nom"].loc[forwards]
    n.model.add_constraints(lhs == rhs, name="Link-bidirectional_sync")


def add_chp_constraints(n):
    electric = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("electric")
    )
    heat = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("heat")
    )

    electric_ext = n.links[electric].query("p_nom_extendable").index
    heat_ext = n.links[heat].query("p_nom_extendable").index

    electric_fix = n.links[electric].query("~p_nom_extendable").index
    heat_fix = n.links[heat].query("~p_nom_extendable").index

    p = n.model["Link-p"]  # dimension: [time, link]

    # output ratio between heat and electricity and top_iso_fuel_line for extendable
    if not electric_ext.empty:
        p_nom = n.model["Link-p_nom"]

        lhs = (
            p_nom.loc[electric_ext]
            * (n.links.p_nom_ratio * n.links.efficiency)[electric_ext].values
            - p_nom.loc[heat_ext] * n.links.efficiency[heat_ext].values
        )
        n.model.add_constraints(lhs == 0, name="chplink-fix_p_nom_ratio")

        rename = {"Link-ext": "Link"}
        lhs = (
            p.loc[:, electric_ext]
            + p.loc[:, heat_ext]
            - p_nom.rename(rename).loc[electric_ext]
        )
        n.model.add_constraints(lhs <= 0, name="chplink-top_iso_fuel_line_ext")

    # top_iso_fuel_line for fixed
    if not electric_fix.empty:
        lhs = p.loc[:, electric_fix] + p.loc[:, heat_fix]
        rhs = n.links.p_nom[electric_fix]
        n.model.add_constraints(lhs <= rhs, name="chplink-top_iso_fuel_line_fix")

    # back-pressure
    if not electric.empty:
        lhs = (
            p.loc[:, heat] * (n.links.efficiency[heat] * n.links.c_b[electric].values)
            - p.loc[:, electric] * n.links.efficiency[electric]
        )
        n.model.add_constraints(lhs <= rhs, name="chplink-backpressure")


def add_pipe_retrofit_constraint(n):
    """
    Add constraint for retrofitting existing CH4 pipelines to H2 pipelines.
    """
    if "reversed" not in n.links.columns:
        n.links["reversed"] = False
    gas_pipes_i = n.links.query(
        "carrier == 'gas pipeline' and p_nom_extendable and ~reversed"
    ).index
    h2_retrofitted_i = n.links.query(
        "carrier == 'H2 pipeline retrofitted' and p_nom_extendable and ~reversed"
    ).index

    if h2_retrofitted_i.empty or gas_pipes_i.empty:
        return

    p_nom = n.model["Link-p_nom"]

    CH4_per_H2 = 1 / n.config["sector"]["H2_retrofit_capacity_per_CH4"]
    lhs = p_nom.loc[gas_pipes_i] + CH4_per_H2 * p_nom.loc[h2_retrofitted_i]
    rhs = n.links.p_nom[gas_pipes_i].rename_axis("Link-ext")

    n.model.add_constraints(lhs == rhs, name="Link-pipe_retrofit")


def add_flexible_egs_constraint(n):
    """
    Upper bounds the charging capacity of the geothermal reservoir according to
    the well capacity.
    """
    well_index = n.links.loc[n.links.carrier == "geothermal heat"].index
    storage_index = n.storage_units.loc[
        n.storage_units.carrier == "geothermal heat"
    ].index

    p_nom_rhs = n.model["Link-p_nom"].loc[well_index]
    p_nom_lhs = n.model["StorageUnit-p_nom"].loc[storage_index]

    n.model.add_constraints(
        p_nom_lhs <= p_nom_rhs,
        name="upper_bound_charging_capacity_of_geothermal_reservoir",
    )


def add_import_limit_constraint(n: pypsa.Network, sns: pd.DatetimeIndex):
    """
    Add constraint for limiting green energy imports (synthetic and biomass).
    Does not include fossil fuel imports.
    """

    import_links = n.links.loc[n.links.carrier.str.contains("import")].index
    import_gens = n.generators.loc[n.generators.carrier.str.contains("import")].index

    limit = n.config["sector"]["imports"]["limit"]
    limit_sense = n.config["sector"]["imports"]["limit_sense"]

    if (import_links.empty and import_gens.empty) or not np.isfinite(limit):
        return

    weightings = n.snapshot_weightings.loc[sns, "generators"]

    # everything needs to be in MWh_fuel
    eff = n.links.loc[import_links, "efficiency"]

    p_gens = n.model["Generator-p"].loc[sns, import_gens]
    p_links = n.model["Link-p"].loc[sns, import_links]

    lhs = (p_gens * weightings).sum() + (p_links * eff * weightings).sum()

    rhs = limit * 1e6

    n.model.add_constraints(lhs, limit_sense, rhs, name="import_limit")


def add_co2_atmosphere_constraint(n, snapshots):
    glcs = n.global_constraints[n.global_constraints.type == "co2_atmosphere"]

    if glcs.empty:
        return
    for name, glc in glcs.iterrows():
        carattr = glc.carrier_attribute
        emissions = n.carriers.query(f"{carattr} != 0")[carattr]

        if emissions.empty:
            continue

        # stores
        bus_carrier = n.stores.bus.map(n.buses.carrier)
        stores = n.stores[bus_carrier.isin(emissions.index) & ~n.stores.e_cyclic]
        if not stores.empty:
            last_i = snapshots[-1]
            lhs = n.model["Store-e"].loc[last_i, stores.index]
            rhs = glc.constant

            n.model.add_constraints(lhs <= rhs, name=f"GlobalConstraint-{name}")


def res_capacity_constraints(n):
    """
    Restrict the deployment of renewable capacities for the same carrier within the same buses.
    """
    rename = {"Generator-ext": "Generator"}

    for carrier in ["solar", "onwind"]:
        ext_carrier = n.generators[
            (n.generators.carrier == carrier) & n.generators.p_nom_extendable
        ]
        p_nom_max = (
            ext_carrier[ext_carrier.p_nom_max != np.inf].groupby("bus").p_nom_max.sum()
        )
        gen = (
            n.model["Generator-p_nom"]
            .rename(rename)
            .loc[ext_carrier.index]
            .groupby(ext_carrier.bus)
            .sum()
        )

        n.model.add_constraints(gen <= p_nom_max, name=f"RES_capacity-{carrier}")


def ember_res_target(n):
    """
    Set a system-wide national RES constraints based on NECPs.

    In comparison to Iegor's 247-cfe paper, this RES target is based on energy generated based on EMBER 2030 Global Renewable Target Tracker.
    CI related generators and links are excluded in this constraint to avoid big overshoot of national RES targets due to CI-procured portfolio.
    Note that EU RE directive counts corporate PPA within NECPs.
    """
    # --- Load and prepare RES targets ---
    df_ember = pd.read_excel(
        "https://storage.googleapis.com/emb-prod-bkt-publicdata/public-downloads/res_tracker/outputs/targets_download.xlsx",
        sheet_name="share_target_wide",
    )

    # Convert ISO3 to ISO2, keeping "EU" unchanged
    df_ember["country"] = df_ember["country_code"].apply(
        lambda code: code
        if code == "EU"
        else cc.convert(names=code, src="ISO3", to="ISO2")
    )

    # --- Define technologies and weights ---
    procurement = n.config["procurement"]
    res_target = procurement["res_target"]
    res_tech = procurement["grid_policy"]["renewable_carriers"]
    weights = n.snapshot_weightings["generators"]

    # --- Helper function to filter and assign country ---
    def get_carriers(df, bus_col):
        bus_list = n.buses[n.buses.carrier == "AC"].index
        grid_carriers = ["electricity distribution grid", "AC", "DC"]

        return (
            df[
                df[bus_col].isin(bus_list)
                & ~df["carrier"].isin(grid_carriers)
                & df["ci"].isin([np.NaN, ""])
            ]
            .copy()
            .assign(country=lambda d: d[bus_col].map(n.buses["country"]))
        )

    # --- Helper function to factor in powerplant efficiencies ---
    def get_link_model(n, df, weights):
        return (
            n.model["Link-p"].loc[:, df.index]
            * df.loc[df.index, "efficiency"]
            * weights
        )

    # Find EU and national targets
    eu_target = df_ember.loc[df_ember.country == "EU", "res_share_target"].values[0]
    countries = list(filter(None, n.buses.country.unique()))
    df_country = df_ember[df_ember.country.isin(countries)]
    country_target = df_country.groupby("country").res_share_target.sum()

    if res_target == "both" and len(countries) == len(country_target):
        logger.info(
            f"All {str(len(countries))} countries have national targets, disable EU-wide RES share target to prevent overconstraints."
        )
        res_target = "country"

    # --- EU-wide RES target constraint ---
    if res_target in ["EU", "both"]:
        logger.info(f"Set EU-wide RES share target to {eu_target}%")

        # --- Apply for generators and links ---
        all_gen_carrier = get_carriers(n.generators, "bus")
        all_link_carrier = get_carriers(n.links, "bus1")

        # Separate RES carriers
        res_gen_carrier = all_gen_carrier.query("carrier in @res_tech")
        res_link_carrier = all_link_carrier.query("carrier in @res_tech")

        all_gen = n.model["Generator-p"].loc[:, all_gen_carrier.index] * weights
        all_link = get_link_model(n, all_link_carrier, weights)

        all_eu = all_gen.sum() + all_link.sum()

        res_gen = n.model["Generator-p"].loc[:, res_gen_carrier.index] * weights
        res_link = get_link_model(n, res_link_carrier, weights)

        res_eu = res_gen.sum() + res_link.sum()

        n.model.add_constraints(
            res_eu == (eu_target / 100) * all_eu, name="EU_res_constraint"
        )

    # --- Country-level RES target constraint ---
    if res_target in ["country", "both"]:
        logger.info(f"Set national RES share targets to {country_target}")

        # Filter carrier dataframes to relevant countries
        all_gen_carrier = get_carriers(n.generators, "bus").query(
            "country in @country_target.index"
        )
        all_link_carrier = get_carriers(n.links, "bus1").query(
            "country in @country_target.index"
        )

        res_gen_carrier = all_gen_carrier.query("carrier in @res_tech")
        res_link_carrier = all_link_carrier.query("carrier in @res_tech")

        # Compute RES and total by country
        all_gen = n.model["Generator-p"].loc[:, all_gen_carrier.index] * weights
        all_link = get_link_model(n, all_link_carrier, weights)

        all_country = (
            all_gen.sum(dim="snapshot").groupby(all_gen_carrier.country).sum()
            + all_link.sum(dim="snapshot").groupby(all_link_carrier.country).sum()
        )

        res_gen = n.model["Generator-p"].loc[:, res_gen_carrier.index] * weights
        res_link = get_link_model(n, res_link_carrier, weights)

        res_country = (
            res_gen.sum(dim="snapshot").groupby(res_gen_carrier.country).sum()
            + res_link.sum(dim="snapshot").groupby(res_link_carrier.country).sum()
        )

        n.model.add_constraints(
            res_country == (country_target / 100) * all_country,
            name="country_res_constraint",
        )


def res_annual_matching_constraints(n):
    """
    Implement strategies for annual renewable procurement matching.

    The total generation from all CI-related generators (renewable carriers) and links (conventional/clean carriers) must equal to its own load consumption.
    """
    weights = n.snapshot_weightings["generators"]
    energy_matching = n.config["procurement"]["energy_matching"] / 100

    for name in n.config["procurement"]["ci"]:
        gen_ci = list(n.generators.query("ci == @name").index)
        links_ci = list(n.links.query("ci == @name").index)

        gen_sum = (n.model["Generator-p"].loc[:, gen_ci] * weights).sum()
        link_sum = (
            n.model["Link-p"].loc[:, links_ci]
            * n.links.loc[links_ci].efficiency
            * weights
        ).sum()
        lhs = gen_sum + link_sum

        total_load = (n.loads_t.p_set[name + " load"] * weights).sum()

        # Note equality sign
        n.model.add_constraints(
            lhs == energy_matching * total_load, name=f"RES_annual_matching_{name}"
        )


def excess_constraints(n):
    """
    Each CI bus must meet its own load consumption before exporting any energy back to the grid.
    """
    weights = n.snapshot_weightings["generators"]

    for name in n.config["procurement"]["ci"]:
        ci_export = n.model["Link-p"].loc[:, [name + " export"]]
        excess = (ci_export * weights).sum()
        total_load = (n.loads_t.p_set[name + " load"] * weights).sum()
        share = n.config["procurement"][
            "excess_share"
        ]  # 'sliding': max(0., energy_matching - 0.8)

        n.model.add_constraints(
            excess <= share * total_load, name=f"Excess_constraint_{name}"
        )


def extra_functionality(
    n: pypsa.Network, snapshots: pd.DatetimeIndex, planning_horizons: str | None = None
) -> None:
    """
    Add custom constraints and functionality.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance with config and params attributes
    snapshots : pd.DatetimeIndex
        Simulation timesteps
    planning_horizons : str, optional
        The current planning horizon year or None in perfect foresight

    Collects supplementary constraints which will be passed to
    ``pypsa.optimization.optimize``.

    If you want to enforce additional custom constraints, this is a good
    location to add them. The arguments ``opts`` and
    ``snakemake.config`` are expected to be attached to the network.
    """
    config = n.config
    constraints = config["solving"].get("constraints", {})
    if constraints["BAU"] and n.generators.p_nom_extendable.any():
        add_BAU_constraints(n, config)
    if constraints["SAFE"] and n.generators.p_nom_extendable.any():
        add_SAFE_constraints(n, config)
    if constraints["CCL"] and n.generators.p_nom_extendable.any():
        add_CCL_constraints(n, config, planning_horizons)

    reserve = config["electricity"].get("operational_reserve", {})
    if reserve.get("activate"):
        add_operational_reserve_margin(n, snapshots, config)

    if EQ_o := constraints["EQ"]:
        add_EQ_constraints(n, EQ_o.replace("EQ", ""))

    if {"solar-hsat", "solar"}.issubset(
        config["electricity"]["renewable_carriers"]
    ) and {"solar-hsat", "solar"}.issubset(
        config["electricity"]["extendable_carriers"]["Generator"]
    ):
        add_solar_potential_constraints(n, config)

    if n.config.get("sector", {}).get("tes", False):
        if n.buses.index.str.contains(
            r"urban central heat|urban decentral heat|rural heat",
            case=False,
            na=False,
        ).any():
            add_TES_energy_to_power_ratio_constraints(n)
            add_TES_charger_ratio_constraints(n)

    add_battery_constraints(n)
    add_lossy_bidirectional_link_constraints(n)
    add_pipe_retrofit_constraint(n)
    if n._multi_invest:
        add_carbon_constraint(n, snapshots)
        add_carbon_budget_constraint(n, snapshots)
        add_retrofit_gas_boiler_constraint(n, snapshots)
    else:
        add_co2_atmosphere_constraint(n, snapshots)

    if config["sector"]["enhanced_geothermal"]["enable"]:
        add_flexible_egs_constraint(n)

    if config["sector"]["imports"]["enable"]:
        add_import_limit_constraint(n, snapshots)

    if n.params.custom_extra_functionality:
        source_path = n.params.custom_extra_functionality
        assert os.path.exists(source_path), f"{source_path} does not exist"
        sys.path.append(os.path.dirname(source_path))
        module_name = os.path.splitext(os.path.basename(source_path))[0]
        module = importlib.import_module(module_name)
        custom_extra_functionality = getattr(module, module_name)
        custom_extra_functionality(n, snapshots, snakemake)  # pylint: disable=E0601

    if (
        n.params.procurement_enable
        and str(n.params.procurement["year"]) == planning_horizons
    ):
        procurement = config["procurement"]
        strategy = procurement["strategy"]
        energy_matching = procurement["energy_matching"]
        res_capacity_constraints(n)

        if procurement["res_target"]:
            ember_res_target(n)

        if strategy == "vol-match":
            logger.info(f"Setting annual volume matching of {energy_matching}%")
            res_annual_matching_constraints(n)
            excess_constraints(n)
        else:
            logger.info("no target set")


def check_objective_value(n: pypsa.Network, solving: dict) -> None:
    """
    Check if objective value matches expected value within tolerance.

    Parameters
    ----------
    n : pypsa.Network
        Network with solved objective
    solving : Dict
        Dictionary containing objective checking parameters

    Raises
    ------
    ObjectiveValueError
        If objective value differs from expected value beyond tolerance
    """
    check_objective = solving["check_objective"]
    if check_objective["enable"]:
        atol = check_objective["atol"]
        rtol = check_objective["rtol"]
        expected_value = check_objective["expected_value"]
        if not np.isclose(n.objective, expected_value, atol=atol, rtol=rtol):
            raise ObjectiveValueError(
                f"Objective value {n.objective} differs from expected value "
                f"{expected_value} by more than {atol}."
            )


def solve_network(
    n: pypsa.Network,
    config: dict,
    params: dict,
    solving: dict,
    rule_name: str | None = None,
    planning_horizons: str | None = None,
    **kwargs,
) -> None:
    """
    Solve network optimization problem.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network instance
    config : Dict
        Configuration dictionary containing solver settings
    params : Dict
        Dictionary of solving parameters
    solving : Dict
        Dictionary of solving options and configuration
    rule_name : str, optional
        Name of the snakemake rule being executed
    planning_horizons : str, optional
            The current planning horizon year or None in perfect foresight
    **kwargs
        Additional keyword arguments passed to the solver

    Returns
    -------
    n : pypsa.Network
        Solved network instance
    status : str
        Solution status
    condition : str
        Termination condition

    Raises
    ------
    RuntimeError
        If solving status is infeasible or warning
    ObjectiveValueError
        If objective value differs from expected value
    """
    set_of_options = solving["solver"]["options"]
    cf_solving = solving["options"]

    kwargs["multi_investment_periods"] = config["foresight"] == "perfect"
    kwargs["solver_options"] = (
        solving["solver_options"][set_of_options] if set_of_options else {}
    )
    kwargs["solver_name"] = solving["solver"]["name"]
    kwargs["extra_functionality"] = partial(
        extra_functionality, planning_horizons=planning_horizons
    )
    kwargs["transmission_losses"] = cf_solving.get("transmission_losses", False)
    kwargs["linearized_unit_commitment"] = cf_solving.get(
        "linearized_unit_commitment", False
    )
    kwargs["assign_all_duals"] = cf_solving.get("assign_all_duals", False)
    kwargs["io_api"] = cf_solving.get("io_api", None)

    if kwargs["solver_name"] == "gurobi":
        logging.getLogger("gurobipy").setLevel(logging.CRITICAL)

    rolling_horizon = cf_solving.pop("rolling_horizon", False)
    skip_iterations = cf_solving.pop("skip_iterations", False)
    if not n.lines.s_nom_extendable.any():
        skip_iterations = True
        logger.info("No expandable lines found. Skipping iterative solving.")

    # add to network for extra_functionality
    n.config = config
    n.params = params

    if rolling_horizon and rule_name == "solve_operations_network":
        kwargs["horizon"] = cf_solving.get("horizon", 365)
        kwargs["overlap"] = cf_solving.get("overlap", 0)
        n.optimize.optimize_with_rolling_horizon(**kwargs)
        status, condition = "", ""
    # elif (
    #     n.params.procurement_enable
    #     and str(n.params.procurement["year"]) == planning_horizons
    #     and n.params.procurement["strategy"] == "247-cfe"
    # ):
    #     status, condition = optimize_cfe_iteratively(n, config, **kwargs)
    elif skip_iterations:
        status, condition = n.optimize(**kwargs)
    else:
        kwargs["track_iterations"] = cf_solving["track_iterations"]
        kwargs["min_iterations"] = cf_solving["min_iterations"]
        kwargs["max_iterations"] = cf_solving["max_iterations"]
        if cf_solving["post_discretization"].pop("enable"):
            logger.info("Add post-discretization parameters.")
            kwargs.update(cf_solving["post_discretization"])
        status, condition = n.optimize.optimize_transmission_expansion_iteratively(
            **kwargs
        )

    if not rolling_horizon:
        if status != "ok":
            logger.warning(
                f"Solving status '{status}' with termination condition '{condition}'"
            )
        check_objective_value(n, solving)

    if "warning" in condition:
        raise RuntimeError("Solving status 'warning'. Discarding solution.")

    if "infeasible" in condition:
        labels = n.model.compute_infeasibilities()
        logger.info(f"Labels:\n{labels}")
        n.model.print_infeasibilities()
        raise RuntimeError("Solving status 'infeasible'. Infeasibilities computed.")


def strip_network(n: pypsa.Network, config: dict) -> None:
    """
    Removes unnecessary components from a pypsa network.

    Args:
    - n (pypsa.Network): The network object to be stripped.

    Returns:
    - None
    """
    ci_names = config["ci"].keys()
    ci_locations = [config["ci"][ci_name]["location"] for ci_name in ci_names]
    zone = set(n.buses.country[bus] for bus in ci_locations)

    # Perform queries and combine results into a single set
    bus_core = n.buses[n.buses["country"].isin(zone)].index.unique()
    combined_lines = n.lines[n.lines.bus1.isin(bus_core) | n.lines.bus0.isin(bus_core)]
    combined_links = n.links[n.links.bus1.isin(bus_core) | n.links.bus0.isin(bus_core)]

    # Combine the results of bus0 and bus1 in lines and links
    bus_connect = (
        set(combined_lines.bus0.unique())
        | set(combined_lines.bus1.unique())
        | set(combined_links.bus0.unique())
        | set(combined_links.bus1.unique())
    )

    zone_all = set(n.buses.country[bus] for bus in bus_connect)
    nodes_to_keep = n.buses[n.buses["country"].isin(zone_all)].index.unique()

    n.remove("Bus", n.buses.index.symmetric_difference(nodes_to_keep))

    # make sure lines are kept
    n.lines.carrier = "AC"

    for c in n.iterate_components(
        ["Generator", "Link", "Line", "Store", "StorageUnit", "Load"]
    ):
        if c.name in ["Link", "Line"]:
            location_boolean = c.df.bus0.isin(nodes_to_keep) & c.df.bus1.isin(
                nodes_to_keep
            )
        else:
            location_boolean = c.df.bus.isin(nodes_to_keep)
        to_keep = c.df.index[location_boolean]
        to_drop = c.df.index.symmetric_difference(to_keep)
        n.remove(c.name, to_drop)


def load_profile(
    n: pypsa.Network,
    location: str,
    profile_shape: str,
    config,
) -> pd.Series:
    """
    Create daily load profile for C&I buyers based on config setting.

    Args:
    - n (object): object
    - profile_shape (str): shape of the load profile, must be one of 'baseload' or 'industry'
    - config (dict): config settings

    Returns:
    - pd.Series: annual load profile for C&I buyers
    """

    procurement = config["procurement"]
    scaling = n.snapshot_weightings.objective.sum() / len(
        n.snapshot_weightings.objective
    )  # e.g., 3 for 3H time resolution

    shapes = {
        "baseload": [1 / 24] * 24,
        "industry": [0.009] * 5
        + [0.016, 0.031, 0.07, 0.072, 0.073, 0.072, 0.07]
        + [0.052, 0.054, 0.066, 0.07, 0.068, 0.063]
        + [0.035] * 2
        + [0.045] * 2
        + [0.009] * 2,
    }

    try:
        shape = shapes[profile_shape]
    except KeyError:
        print(
            f"'profile_shape' option must be one of 'baseload' or 'industry'. Now is {profile_shape}."
        )
        sys.exit()

    if procurement["strategy"] == "ref":
        load = 0.0
    else:
        country = n.buses.country[location]
        load_year = (
            pd.read_csv(procurement["load"])
            .groupby("Country Code")["2023"]
            .sum()[country]
        )  # GWh
        load = load_year / 8760 * 1000 * procurement["participation"] / 100  # MW

    load_day = load * 24  # daily load in MWh
    load_profile_day = pd.Series(shape) * load_day
    load_profile_year = pd.concat([load_profile_day] * 365)

    if scaling != 1.0:
        load_profile_year.index = pd.date_range(
            start="2013-01-01", periods=len(load_profile_year), freq="h"
        )
        profile = (
            load_profile_year.resample(f"{int(scaling)}h")
            .mean()
            .reindex(n.snapshots, method="nearest")
        )
    else:
        profile = load_profile_year.set_axis(n.snapshots)

    return profile


def add_ci(n: pypsa.Network, year: str, config: dict, costs: pd.DataFrame) -> None:
    """
    Add C&I buyer(s) to the network.

    Args:
    - n: pypsa.Network to which the C&I buyer(s) will be added.
    - year: the year of optimisation based on config setting.

    Returns:
    - None
    """
    # tech_palette options
    procurement = config["procurement"]
    clean_techs = procurement["technology"]["generation_tech"]
    storage_techs = procurement["technology"]["storage_tech"]
    ci = procurement["ci"]
    strategy = procurement["strategy"]
    scope = procurement["scope"]

    for name in ci.keys():
        location = ci[name]["location"]
        profile = procurement["profile"]

        n.add("Bus", name, country="")

        n.add(
            "Link",
            f"{name}" + " export",
            bus0=name,
            bus1=location,
            marginal_cost=0.1,  # large enough to avoid optimization artifacts, small enough not to influence PPA portfolio
            p_nom=1e6,
            reversed=False,
        )

        n.add(
            "Link",
            f"{name}" + " import",
            bus0=location,
            bus1=name,
            marginal_cost=0.001,  # large enough to avoid optimization artifacts, small enough not to influence PPA portfolio
            p_nom=1e6,
            reversed=False,
        )

        n.add(
            "Load",
            f"{name}" + " load",
            carrier="electricity",
            bus=name,
            p_set=load_profile(n, location, profile, config),
            ci=name,  # C&I markers used in constraints
        )

        # C&I following voluntary clean energy procurement is a share of C&I load -> subtract it from node's profile
        n.loads_t.p_set[location] -= n.loads_t.p_set[f"{name}" + " load"]

        # ===================== Adding Dispatchable Technologies =====================
        # ============================================================================

        gen_implemented = {
            "nuclear": {
                "carrier": "uranium",
                "carrier_nodes": "EU uranium",
                "unit": "MWh_th",
            },
            "allam": {
                "carrier": "gas",
                "carrier_nodes": "EU gas",
                "unit": "MWh_LHV",
            },
            "geothermal": {
                "carrier": "geothermal",
                "carrier_nodes": "EU enhanced geothermal systems",
                "unit": "MWh_th",
            },
        }
        gen_not_implemented = list(
            set(clean_techs).difference(
                list(gen_implemented.keys()) + ["onwind", "solar"]
            )
        )
        gen_available_carriers = list(
            set(clean_techs).intersection(gen_implemented.keys())
        )
        if len(gen_not_implemented) > 0:
            logger.warning(
                f"{gen_not_implemented} are not yet implemented as Clean technologies for CI in PyPSA-Eur"
            )

        for generator in gen_available_carriers:
            carrier = gen_implemented[generator]["carrier"]
            carrier_nodes = gen_implemented[generator]["carrier_nodes"]

            if carrier_nodes not in n.buses.index:
                logger.info(f"Missing buses: {carrier_nodes}. Adding them now for CI")
                n.add(
                    "Bus",
                    carrier_nodes,
                    carrier=carrier,
                    location="EU",
                    unit=gen_implemented[generator]["unit"],
                )

            n.add(
                "Link",
                name + " " + generator,
                bus0=carrier_nodes,
                bus1=name,
                bus2="co2 atmosphere",
                marginal_cost=costs.at[generator, "efficiency"]
                * costs.at[generator, "VOM"],  # NB: VOM is per MWel
                capital_cost=costs.at[generator, "efficiency"]
                * costs.at[generator, "capital_cost"],  # NB: fixed cost is per MWel
                p_nom_extendable=True if strategy else False,
                p_max_pu=0.7
                if carrier == "uranium"
                else 1,  # be conservative for nuclear (maintenance or unplanned shut downs)
                carrier=generator,
                efficiency=costs.at[generator, "efficiency"],
                efficiency2=costs.at[carrier, "CO2 intensity"],
                lifetime=costs.at[generator, "lifetime"],
                reversed=False,
                ci=name,  # C&I markers used in constraints
            )

        add_missing_carriers(n, gen_available_carriers)

        logger.info(f"Include {gen_available_carriers} for the CI: {name}.")

        # ===================== Adding Variable Renewable Technologies =====================
        # ==================================================================================

        res_available_carriers = list(
            set(clean_techs).intersection(["onwind", "solar"])
        )

        for carrier in res_available_carriers:
            if scope == "node" or strategy == "247-cfe":
                res_df = pd.DataFrame(index=[name])
                res_df["gen_name"] = name + " " + carrier
                res_df["gen_template"] = location + " " + carrier + f"-{year}"
            else:
                if scope == "country":
                    zone = n.buses.loc[location, "country"]
                    index = n.buses[
                        (n.buses.carrier == "AC") & (n.buses.country == zone)
                    ].index
                else:  # scope == "all" is the default
                    index = n.buses[
                        (n.buses.carrier == "AC") & (n.buses.country != "")
                    ].index

                res_df = pd.DataFrame(index=index)
                res_df["gen_name"] = [f"{i} {name} {carrier}" for i in res_df.index]
                res_df["gen_template"] = [f"{i} {carrier}-{year}" for i in res_df.index]

            p_max_pu_df = n.generators_t.p_max_pu[res_df["gen_template"]]
            p_max_pu_df = p_max_pu_df.rename(
                columns=res_df.set_index("gen_template")["gen_name"].to_dict()
            )

            n.add(
                "Generator",
                res_df["gen_name"],
                carrier=carrier,
                bus=res_df.index,
                p_nom_extendable=True if strategy else False,
                p_max_pu=p_max_pu_df,
                capital_cost=costs.at[carrier, "capital_cost"],
                marginal_cost=costs.at[carrier, "marginal_cost"],
                ci=name,  # C&I markers used in constraints
            )

        logger.info(
            f"Include {res_available_carriers} for the CI: {name} with the scope: {scope}."
        )

        # ===================== Adding Storage Technologies =====================
        # =======================================================================

        max_hours = config["max_hours"]

        # check for not implemented storage technologies
        storage_implemented = [
            "H2",
            "li-ion battery",
            "iron-air battery",
            "lfp",
            "vanadium",
            "lair",
            "pair",
        ]
        storage_not_implemented = list(
            set(storage_techs).difference(storage_implemented)
        )
        storage_available_carriers = list(
            set(storage_techs).intersection(storage_implemented)
        )
        if len(storage_not_implemented) > 0:
            logger.warning(
                f"{storage_not_implemented} are not yet implemented as Storage technologies in PyPSA-Eur"
            )
        available_carriers_max_hours = [
            f"{carrier} {max_hour}h"
            for carrier in storage_available_carriers
            if carrier in max_hours
            for max_hour in max_hours[carrier]
        ]
        missing_carriers = list(
            set(available_carriers_max_hours).difference(n.carriers.index)
        )
        n.add("Carrier", missing_carriers)

        lookup_store = {
            "H2": "electrolysis",
            "li-ion battery": "battery inverter",
            "iron-air battery": "iron-air battery charge",
            "lfp": "Lithium-Ion-LFP-bicharger",
            "vanadium": "Vanadium-Redox-Flow-bicharger",
            "lair": "Liquid-Air-charger",
            "pair": "Compressed-Air-Adiabatic-bicharger",
        }
        lookup_dispatch = {
            "H2": "fuel cell",
            "li-ion battery": "battery inverter",
            "iron-air battery": "iron-air battery discharge",
            "lfp": "Lithium-Ion-LFP-bicharger",
            "vanadium": "Vanadium-Redox-Flow-bicharger",
            "lair": "Liquid-Air-discharger",
            "pair": "Compressed-Air-Adiabatic-bicharger",
        }

        for carrier in storage_available_carriers:
            for max_hour in max_hours[carrier]:
                roundtrip_correction = 0.5 if carrier == "li-ion battery" else 1
                cost_carrier = "H2 tank" if carrier == "H2" else carrier
                n.add(
                    "StorageUnit",
                    name,
                    suffix=f" {carrier} {max_hour}h",
                    bus=name,
                    carrier=f"{carrier} {max_hour}h",
                    p_nom_extendable=True if strategy else False,
                    capital_cost=costs.at[
                        f"{cost_carrier} {max_hour}h", "capital_cost"
                    ],
                    marginal_cost=0.0,
                    efficiency_store=costs.at[lookup_store[carrier], "efficiency"]
                    ** roundtrip_correction,
                    efficiency_dispatch=costs.at[lookup_dispatch[carrier], "efficiency"]
                    ** roundtrip_correction,
                    max_hours=max_hour,
                    cyclic_state_of_charge=True,
                    lifetime=costs.at[f"{cost_carrier} {max_hour}h", "lifetime"],
                    ci=name,  # C&I markers used in constraints
                )

        logger.info(f"Include {storage_available_carriers} for the CI: {name}.")


# %%
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_sector_network_myopic",
            run="baseline-3H",
            opts="",
            clusters="39",
            configfiles="config/config.meta.yaml",
            ll="v1.0",
            sector_opts="",
            planning_horizons="2030",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    solve_opts = snakemake.params.solving["options"]

    np.random.seed(solve_opts.get("seed", 123))

    n = pypsa.Network(snakemake.input.network)
    planning_horizons = snakemake.wildcards.get("planning_horizons", None)

    prepare_network(
        n,
        solve_opts=snakemake.params.solving["options"],
        foresight=snakemake.params.foresight,
        planning_horizons=planning_horizons,
        co2_sequestration_potential=snakemake.params["co2_sequestration_potential"],
        limit_max_growth=snakemake.params.get("sector", {}).get("limit_max_growth"),
    )

    logging_frequency = snakemake.config.get("solving", {}).get(
        "mem_logging_frequency", 30
    )
    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=logging_frequency
    ) as mem:
        if (
            snakemake.params.procurement_enable
            and str(snakemake.params.procurement["year"]) == planning_horizons
        ):
            print("procurement_enable is activated")
            procurement = snakemake.params.procurement

            if procurement["strip_network"]:
                print("stript_network is activated")
                strip_network(n, procurement)

            Nyears = n.snapshot_weightings.objective.sum() / 8760.0
            costs = load_costs(
                snakemake.input.costs,
                snakemake.params.costs,
                snakemake.params.max_hours,
                Nyears,
            )
            add_ci(n, snakemake.wildcards.planning_horizons, snakemake.params, costs)

        solve_network(
            n,
            config=snakemake.config,
            params=snakemake.params,
            solving=snakemake.params.solving,
            planning_horizons=planning_horizons,
            rule_name=snakemake.rule,
        )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    grid_policy = snakemake.params.procurement.get("grid_policy", False)
    if grid_policy:
        res_techs = grid_policy["renewable_carriers"]
        clean_techs = grid_policy["clean_carriers"]
        calculate_grid_score(n, res_techs, "res")
        calculate_grid_score(n, clean_techs, "cfe")
        calculate_grid_score(n, res_techs, "res_w_ci", include_ci=True)
        calculate_grid_score(n, clean_techs, "cfe_w_ci", include_ci=True)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)

    with open(snakemake.output.config, "w") as file:
        yaml.dump(
            n.meta,
            file,
            default_flow_style=False,
            allow_unicode=True,
            sort_keys=False,
        )
