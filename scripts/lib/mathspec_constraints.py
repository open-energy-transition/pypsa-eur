# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
The data of the constraints a datarecord of a PyPSA-Eur network keeps beside PyPSA's model.

The math is static: `mathspec/global_constraints.yaml` holds the global
constraints of the network, collapsed by family, and
`mathspec/extra_functionality.yaml` holds what `solve_network` adds. Both
compose onto the base spec that datarecord ships. This module writes their
data as tables of the record: a `datarecord.FrameRecord` whose schema is the
dims, relations and parameters the two files declare and the network record
does not. A constraint has data only where the network and the config
activate it, and an active constraint that the math cannot state raises
:class:`UnsupportedConstraintError`.
"""

import ast
import re
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import Any

import narwhals as nw
import numpy as np
import pandas as pd
import pypsa
import yaml
from datarecord import AttributeSpec, Dimension, FrameRecord, Group, Schema
from datarecord.combine import combine_schemas

from scripts.solve_network import emission_stores, get_glc_sense

MATH = Path(__file__).with_name("mathspec")
SPECS = (MATH / "global_constraints.yaml", MATH / "extra_functionality.yaml")
DTYPES = {"float": nw.Float64, "int": nw.Int64, "bool": nw.Boolean, "str": nw.String}
BOUNDS = {"<=": "upper", ">=": "lower", "==": "equal"}
LAND_USE = "tech_capacity_expansion_limit"
KINDS = {
    ("primary_energy", "co2_emissions"): "co2_limit",
    ("primary_energy", "gas_usage"): "gas_limit",
    ("co2_atmosphere", "co2_emissions"): "co2_atmosphere",
    ("Co2Budget", "co2_emissions"): "co2_budget",
    ("operational_limit", "co2 sequestered"): "co2_sequestration_limit",
    ("operational_limit", "solid biomass"): "biomass_limit",
    ("operational_limit", "unsustainable solid biomass"): "unsustainable_biomass_limit",
    ("operational_limit", "municipal solid waste"): "msw_limit",
}
TRANSMISSION_KINDS = {
    "transmission_volume_expansion_limit": "transmission_volume_limit",
    "transmission_expansion_cost_limit": "transmission_cost_limit",
}
CAPACITY_TYPES = ["Generator", "Line", "Transformer", "Link", "Store", "StorageUnit"]
ONE_PORT_TYPES = ["Generator", "StorageUnit", "Store"]

Data = dict[str, Any]


class UnsupportedConstraintError(NotImplementedError):
    pass


def _declarations() -> dict[str, dict[str, Any]]:
    specs = [yaml.safe_load(path.read_text()) for path in SPECS]
    return {
        kind: {name: entry for spec in specs for name, entry in spec[kind].items()}
        for kind in ("dimensions", "relations", "parameters")
    }


def _group(relation: Mapping[str, Any]) -> Group:
    key, values = relation["key"], relation.get("values")
    if values is None:
        return Group(over=list(key))
    if isinstance(values, Mapping):
        return Group(over={key: key, **values})
    return Group(over=[key], into=values)


def _with_parents(dims: Iterable[str], base: Schema) -> set[str]:
    closed = set(dims)
    for dim in list(closed):
        closed |= (
            _with_parents(base.dimensions[dim].within, base)
            if dim in base.dimensions
            else set()
        )
    return closed


def extension_schema(base: Schema) -> Schema:
    """The dims, relations and parameters of the static math that ``base`` lacks."""
    declared = _declarations()
    groups = {
        name: _group(relation)
        for name, relation in declared["relations"].items()
        if name not in base.groups and name != "snapshot_period"
    }
    attributes = {
        name: AttributeSpec(
            dtype=DTYPES[spec.get("dtype", "float")](),
            dims=frozenset(spec["dims"]),
            description=spec.get("description"),
        )
        for name, spec in declared["parameters"].items()
    }
    new = {
        dim: Dimension(dtype=nw.String(), description=spec.get("description"))
        for dim, spec in declared["dimensions"].items()
        if dim not in base.dimensions
    }
    used = {d for spec in attributes.values() for d in spec.dims}
    for group in groups.values():
        used |= set(group.over.values())
    shared = _with_parents(used - set(new), base)
    return Schema(
        dimensions={d: base.dimensions[d] for d in sorted(shared)} | new,
        attributes=attributes,
        groups=groups,
    )


def _long(
    attribute: str, values: pd.Series | float | str, columns: tuple[str, ...]
) -> pd.DataFrame:
    if isinstance(values, pd.Series):
        frame = values.rename("value").reset_index()
    else:
        frame = pd.DataFrame({"value": [values]})
    frame["attribute"] = attribute
    frame["breakpoint"] = np.nan
    return frame[list(columns)]


def extension(data: Data, base: Schema) -> FrameRecord:
    """
    Pack ``data``, keyed by the name the math reads, into tables of the record.

    A dim maps to its labels, an attribute to a Series over its dims (a scalar
    for one over no dims), and a relation to a frame of its coordinates.
    """
    schema = extension_schema(base)
    combined = combine_schemas(base, schema)
    unknown = sorted(
        set(data) - set(schema.attributes) - set(schema.groups) - set(schema.dimensions)
    )
    if unknown:
        raise ValueError(f"No static math declares {unknown}.")
    dim_frames: dict[str, pd.DataFrame] = {
        dim: pd.DataFrame({dim: pd.Index(data.get(dim, []), dtype=str)})
        for dim in schema.dimensions
        if dim not in base.dimensions
    }
    attribute_frames = {}
    for name, spec in schema.attributes.items():
        if name not in data:
            continue
        if len(spec.dims) == 1:
            (dim,) = spec.dims
            column = data[name].rename(name).rename_axis(dim).reset_index()
            frame = dim_frames.get(dim)
            dim_frames[dim] = (
                column if frame is None else frame.merge(column, on=dim, how="outer")
            )
        else:
            attribute_frames[name] = _long(
                name, data[name], combined.long_columns_for(name)
            )
    group_frames = {
        name: data[name][list(group.coordinates)]
        for name, group in schema.groups.items()
        if name in data
    }
    return FrameRecord(
        schema,
        dim_frames=dim_frames,
        group_frames=group_frames,
        attribute_frames=attribute_frames,
    )


def _ones(index: pd.Index | pd.MultiIndex) -> pd.Series:
    return pd.Series(1, index=index)


def _standing(n: pypsa.Network) -> pd.DataFrame:
    """Whether each entity with a capacity stands in each period, by PyPSA's activity rule."""
    return pd.concat(
        [
            pd.concat({p: c.get_active_assets(p) for p in n.investment_periods}, axis=1)
            for c in n.components[CAPACITY_TYPES]
        ]
    )


def _extendable(n: pypsa.Network) -> pd.Index:
    return pd.Index(
        [name for c in n.components[CAPACITY_TYPES] for name in c.extendables]
    )


def _entity_carriers(n: pypsa.Network) -> pd.Series:
    return pd.concat(
        [
            c.static.carrier
            for c in n.components[CAPACITY_TYPES]
            if "carrier" in c.static
        ]
    )


def _indicator(mask: pd.DataFrame) -> pd.Series:
    stacked = mask.rename_axis(index="entity", columns="period").stack()
    return _ones(stacked.index[stacked])


def _at_bus(values: pd.Series, buses: pd.Series) -> pd.Series:
    """``values``, keyed by entity first, with the bus of each entity after it."""
    frame = values.rename("value").reset_index()
    frame.insert(1, "bus", frame.entity.map(buses).to_numpy())
    return frame.set_index(list(frame.columns[:-1])).value


def _one_port_buses(n: pypsa.Network) -> pd.Series:
    return pd.concat([c.static.bus for c in n.components[ONE_PORT_TYPES]])


def _period_ends(n: pypsa.Network, entities: pd.Index, values: pd.Series) -> pd.Series:
    """``values`` of each period at its last snapshot, for each of ``entities``."""
    ends = n.snapshots.to_frame(index=False).groupby("period").tail(1)
    ends.columns = ["period", "snapshot"]
    frame = ends.merge(pd.DataFrame({"entity": entities}), how="cross")
    frame["value"] = frame.period.map(values).to_numpy()
    return frame.set_index(["period", "snapshot", "entity"]).value


def _require_multi_period(n: pypsa.Network) -> None:
    if n.investment_periods.empty:
        raise ValueError("The constraint math needs a network with investment periods.")
    if n.has_scenarios:
        raise UnsupportedConstraintError(
            "The constraint math does not cover stochastic networks."
        )


def _storage_check(n: pypsa.Network, carriers: Iterable[str]) -> None:
    for c, cyclic, per_period in (
        (n.c.stores, "e_cyclic", "e_initial_per_period"),
        (
            n.c.storage_units,
            "cyclic_state_of_charge",
            "state_of_charge_initial_per_period",
        ),
    ):
        static = c.static
        continuous = static.index[
            static.carrier.isin(list(carriers))
            & ~static[cyclic]
            & ~static[per_period]
            & static.active
        ]
        if not continuous.empty:
            raise UnsupportedConstraintError(
                f"Non-cyclic {c.name} components {list(continuous)} deplete over all "
                "investment periods; PyPSA refuses them where period years differ "
                "from 1, and the math covers only storage that resets each period."
            )


def _transmission_carriers(attribute: str) -> list[str]:
    return [re.sub(r"[\[\]\(\)]", "", c.strip()) for c in attribute.split(",")]


def _glc_kind(name: str, glc: pd.Series) -> tuple[str, list[str]]:
    """The constraint a row of ``n.global_constraints`` is, and the carriers it selects."""
    if glc.type in TRANSMISSION_KINDS:
        return TRANSMISSION_KINDS[glc.type], _transmission_carriers(
            glc.carrier_attribute
        )
    kind = KINDS.get((glc.type, glc.carrier_attribute))
    if kind is None:
        raise UnsupportedConstraintError(
            f"Global constraint {name} of type {glc.type} on {glc.carrier_attribute!r} "
            "has no mathspec translation."
        )
    selects = glc.type == "operational_limit"
    return kind, [glc.carrier_attribute] if selects else []


def global_constraint_data(n: pypsa.Network) -> Data:
    """
    The limits of ``n.global_constraints``, one parameter per constraint and sense.

    A row binding per investment period sets ``<kind>_<bound>[period]``, one over
    the horizon ``<kind>_horizon_<bound>``, the bound being the row's sense.
    A constraint that selects carriers reads them from ``<kind>_carrier``.
    """
    _require_multi_period(n)
    glcs = n.global_constraints[n.global_constraints.type != LAND_USE]
    if glcs.empty:
        return {}
    declared = _declarations()["parameters"]
    limits: dict[str, dict[int | None, float]] = {}
    selections: dict[str, list[str]] = {}
    for name, glc in glcs.iterrows():
        kind, carriers = _glc_kind(name, glc)
        bound = BOUNDS[get_glc_sense(glc)]
        per_period = pd.notna(glc.investment_period)
        parameter = f"{kind}_{bound}" if per_period else f"{kind}_horizon_{bound}"
        if parameter not in declared:
            raise UnsupportedConstraintError(
                f"Global constraint {name} would set '{parameter}', which no static "
                "math declares."
            )
        period = int(glc.investment_period) if per_period else None
        if per_period and period not in n.investment_periods:
            raise UnsupportedConstraintError(
                f"Global constraint {name} names period {period}, which the network does not have."
            )
        if period in limits.setdefault(parameter, {}):
            raise UnsupportedConstraintError(
                f"Global constraints set '{parameter}' twice for {period}."
            )
        limits[parameter][period] = float(glc.constant)
        if carriers and selections.setdefault(kind, carriers) != carriers:
            raise UnsupportedConstraintError(
                f"Global constraints of '{kind}' select different carriers."
            )
    kinds = {_glc_kind(name, glc)[0] for name, glc in glcs.iterrows()}
    if kinds & {"co2_atmosphere", "co2_budget"}:
        emitting = n.carriers.co2_emissions[n.carriers.co2_emissions != 0].index
        emission_stores(n, emitting)
    for kind, attribute in (("co2_limit", "co2_emissions"), ("gas_limit", "gas_usage")):
        if kind in kinds:
            weights = n.carriers[attribute]
            _storage_check(n, weights[weights.fillna(0) != 0].index)
    for carriers in selections.values():
        _storage_check(n, carriers)
    data: Data = {
        parameter: values[None]
        if None in values
        else pd.Series(values).rename_axis("period")
        for parameter, values in limits.items()
    }
    data |= {
        f"{kind}_carrier": _ones(
            pd.Index([c for c in carriers if c in n.carriers.index], name="carrier")
        )
        for kind, carriers in selections.items()
    }
    return data | global_constraint_weights(n)


def global_constraint_weights(n: pypsa.Network) -> Data:
    """The entities each family of global constraints sums over, as weights per entity."""
    if not n.generators_t.efficiency.empty:
        raise UnsupportedConstraintError(
            f"Generators {list(n.generators_t.efficiency.columns)} have a time-varying "
            "efficiency; the math reads one fuel weight per generator."
        )
    extendable = _standing(n).loc[_extendable(n)]
    branches = extendable[extendable.index.isin(n.lines.index.union(n.links.index))]
    years = n.investment_period_weightings.years
    stores, units = n.stores, n.storage_units
    final_stores = stores.index[
        stores.active & ~stores.e_cyclic & stores.e_initial_per_period
    ]
    final_units = units.index[
        units.active
        & ~units.cyclic_state_of_charge
        & units.state_of_charge_initial_per_period
    ]
    emitting = stores.bus.map(n.buses.carrier).map(n.carriers.co2_emissions)
    co2_stores = stores.index[~stores.e_cyclic & emitting.fillna(0).ne(0)]
    return {
        "glc_generator_fuel": 1 / n.generators.efficiency,
        "glc_store_final": _period_ends(n, final_stores, years),
        "glc_storage_unit_final": _period_ends(n, final_units, years),
        "glc_co2_store": _period_ends(n, co2_stores, _ones(years.index)),
        "glc_branch_horizon": _ones(branches.index[branches.any(axis=1)]),
    }


def _pairs(name: str, pairs: Iterable[tuple[str, str]]) -> Data:
    frame = pd.DataFrame(list(pairs), columns=["entity", "partner"])
    return {name: frame} if not frame.empty else {}


def _positional(left: pd.Index, right: pd.Index, source: str) -> list[tuple[str, str]]:
    if len(left) != len(right):
        raise UnsupportedConstraintError(
            f"{source} pairs {len(left)} with {len(right)} entities; linopy aligns "
            "them by an outer join then, which the math does not state."
        )
    return list(zip(left, right))


def bau_constraints(n: pypsa.Network, config: Mapping[str, Any]) -> Data:
    """``solve_network.add_BAU_constraints``."""
    mincaps = pd.Series(config["electricity"]["BAU_mincapacities"])
    extendable = n.generators.query("p_nom_extendable")
    return {"bau_mincapacity": mincaps[extendable.carrier.unique()].astype(float)}


def safe_constraints(n: pypsa.Network, config: Mapping[str, Any]) -> Data:
    """``solve_network.add_SAFE_constraints``."""
    peak = n.loads_t.p_set.sum(axis=1).max()
    gens = n.generators[
        n.generators.carrier.isin(config["electricity"]["conventional_carriers"])
    ]
    return {
        "safe_peak_capacity": float(
            peak * (1.0 + config["electricity"]["SAFE_reservemargin"])
        ),
        "safe_extendable": _ones(gens.index[gens.p_nom_extendable]),
        "safe_fixed": _ones(gens.index[~gens.p_nom_extendable]),
    }


def ccl_constraints(n: pypsa.Network, config: Mapping[str, Any]) -> Data:
    raise UnsupportedConstraintError(
        "solve_network.add_CCL_constraints is not implemented for perfect foresight "
        "(solve_network asserts a planning horizon), so it has no mathspec translation."
    )


def operational_reserve_margin(n: pypsa.Network, config: Mapping[str, Any]) -> Data:
    """``solve_network.add_operational_reserve_margin``."""
    reserve = config["electricity"]["operational_reserve"]
    vres = [
        "ror" if c == "hydro" else c
        for c in config["electricity"]["renewable_carriers"]
    ]
    gens = n.generators.query("carrier != 'load'")
    no_profile = gens.index[gens.carrier.isin(vres)].difference(
        n.generators_t.p_max_pu.columns
    )
    if not no_profile.empty:
        raise KeyError(
            f"Renewable generators {list(no_profile)} have no p_max_pu profile."
        )
    vres_gens = gens[gens.carrier.isin(vres)]
    return {
        "reserve_generator": pd.Series(True, index=gens.index),
        "reserve_vres_extendable": _ones(vres_gens.index[vres_gens.p_nom_extendable]),
        "reserve_vres_fixed": _ones(vres_gens.index[~vres_gens.p_nom_extendable]),
        "reserve_epsilon_load": float(reserve["epsilon_load"]),
        "reserve_epsilon_vres": float(reserve["epsilon_vres"]),
        "reserve_contingency": float(reserve["contingency"]),
    }


def eq_constraints(n: pypsa.Network, option: str) -> Data:
    """``solve_network.add_EQ_constraints``, per node."""
    if option[-1] == "c":
        raise UnsupportedConstraintError(
            "add_EQ_constraints per country needs a country dimension; the record "
            "carries country only as a bus attribute."
        )
    static_loads = n.loads.index.difference(n.loads_t.p_set.columns)
    if n.loads.p_set[static_loads].ne(0).any():
        raise UnsupportedConstraintError(
            "add_EQ_constraints reads loads_t.p_set only; the record merges static "
            f"p_set of loads {list(static_loads)} into it."
        )
    return {"equity_min_share": float(re.findall(r"[0-9]*\.?[0-9]+", option)[0])}


def solar_potential_constraints(n: pypsa.Network, config: Mapping[str, Any]) -> Data:
    """``solve_network.add_solar_potential_constraints``."""
    factor = (
        config["renewable"]["solar"]["capacity_per_sqkm"]
        / config["renewable"]["solar-hsat"]["capacity_per_sqkm"]
    )
    solar = n.generators[
        n.generators.carrier.isin(["solar", "solar-hsat"])
        & n.generators.p_nom_extendable
    ]
    if solar.empty:
        return {}
    if (
        solar.index.str.contains("solar-hsat") != (solar.carrier == "solar-hsat")
    ).any():
        raise UnsupportedConstraintError(
            "add_solar_potential_constraints scales generators by name; names and "
            "carriers of the solar generators disagree."
        )
    factors = pd.Series({"solar": 1.0, "solar-hsat": factor})
    factors = factors[factors > 0]
    weighted = solar[solar.carrier.isin(factors.index)]
    reference = weighted[weighted.carrier == "solar"]
    tracking = n.generators[
        n.generators.carrier.isin(factors.index.drop("solar", errors="ignore"))
    ]
    left = (
        reference.groupby("bus")
        .p_nom_max.sum()
        .sub(
            (tracking.p_nom * tracking.carrier.map(factors)).groupby(tracking.bus).sum()
        )
    )
    limit = left.dropna().clip(lower=0)
    return {
        "solar_potential_generator": _at_bus(
            weighted.carrier.map(factors).rename_axis("entity"), weighted.bus
        ),
        "solar_potential_limit": limit[np.isfinite(limit)],
    }


def tes_energy_to_power_ratio_constraints(n: pypsa.Network) -> Data:
    """``solve_network.add_TES_energy_to_power_ratio_constraints``."""
    chargers = n.links.index[
        n.links.index.str.contains("water tanks charger|water pits charger")
        & n.links.p_nom_extendable
    ]
    stores = n.stores.index[
        n.stores.index.str.contains("water tanks|water pits")
        & n.stores.e_nom_extendable
    ]
    if chargers.empty or stores.empty:
        return {}
    pairs = list(zip(stores, chargers))
    for tes, charger in pairs:
        if tes != charger.replace(" charger", ""):
            raise RuntimeError(f"Charger {charger} and TES {tes} do not match.")
    ratio = n.links.loc[[c for _, c in pairs], "energy to power ratio"].set_axis(
        [s for s, _ in pairs]
    )
    return _pairs("tes_store_charger", pairs) | {"tes_energy_to_power_ratio": ratio}


def tes_charger_ratio_constraints(n: pypsa.Network) -> Data:
    """``solve_network.add_TES_charger_ratio_constraints``."""
    kinds = "water tanks {0}|water pits {0}|aquifer thermal energy storage {0}"
    ext = n.links.p_nom_extendable
    chargers = n.links.index[n.links.index.str.contains(kinds.format("charger")) & ext]
    dischargers = n.links.index[
        n.links.index.str.contains(kinds.format("discharger")) & ext
    ]
    if chargers.empty or dischargers.empty:
        return {}
    pairs = list(zip(chargers, dischargers))
    for charger, discharger in pairs:
        if charger.replace(" charger", " ") != discharger.replace(" discharger", " "):
            raise RuntimeError(
                f"Charger {charger} and discharger {discharger} do not match."
            )
    return _pairs("tes_charger_discharger", pairs) | {
        "discharger_efficiency": _partner_efficiency(n, pairs)
    }


def _partner_efficiency(n: pypsa.Network, pairs: list[tuple[str, str]]) -> pd.Series:
    return n.links.efficiency[[d for _, d in pairs]].set_axis([c for c, _ in pairs])


def battery_constraints(n: pypsa.Network) -> Data:
    """``solve_network.add_battery_constraints``."""
    ext = n.links.query("p_nom_extendable")
    if ext.empty:
        return {}
    dischargers = ext.index[ext.index.str.contains("battery discharger")]
    chargers = ext.index[ext.index.str.contains("battery charger")]
    pairs = _positional(chargers, dischargers, "solve_network.add_battery_constraints")
    return _pairs("battery_charger_discharger", pairs) | {
        "discharger_efficiency": _partner_efficiency(n, pairs)
    }


def lossy_bidirectional_link_constraints(n: pypsa.Network) -> Data:
    """``solve_network.add_lossy_bidirectional_link_constraints``."""
    links = n.links
    if (
        not links.p_nom_extendable.any()
        or "reversed" not in links
        or not links.reversed.any()
    ):
        return {}
    carriers = links.loc[links.reversed, "carrier"].unique()
    backwards = links.index[
        links.carrier.isin(carriers)
        & links.p_nom_extendable
        & links.reversed
        & links.active
    ]
    forwards = backwards.str.replace("-reversed", "")
    missing = forwards.difference(links.index[links.p_nom_extendable])
    if not missing.empty:
        raise KeyError(
            f"Forward links {list(missing)} of reversed links are not extendable."
        )
    return _pairs("reversed_forward", zip(backwards, forwards))


def pipe_retrofit_constraint(n: pypsa.Network, config: Mapping[str, Any]) -> Data:
    """``solve_network.add_pipe_retrofit_constraint``."""
    links = n.links.assign(
        reversed=n.links["reversed"] if "reversed" in n.links else False
    )
    candidates = links.p_nom_extendable & ~links.reversed.astype(bool) & links.active
    gas = links.index[candidates & (links.carrier == "gas pipeline")]
    h2 = links.index[candidates & (links.carrier == "H2 pipeline retrofitted")]
    if gas.empty or h2.empty:
        return {}
    pairs = _positional(gas, h2, "solve_network.add_pipe_retrofit_constraint")
    return _pairs("pipe_retrofit", pairs) | {
        "pipe_retrofit_ch4_per_h2": 1 / config["sector"]["H2_retrofit_capacity_per_CH4"]
    }


def flexible_egs_constraint(n: pypsa.Network) -> Data:
    """``solve_network.add_flexible_egs_constraint``."""
    wells = n.links.index[n.links.carrier == "geothermal heat"]
    storage = n.storage_units.index[n.storage_units.carrier == "geothermal heat"]
    return _pairs(
        "egs_storage_well",
        _positional(storage, wells, "solve_network.add_flexible_egs_constraint"),
    )


def import_limit_constraint(n: pypsa.Network, config: Mapping[str, Any]) -> Data:
    """``solve_network.add_import_limit_constraint``."""
    links = n.links.index[n.links.carrier.str.contains("import")]
    gens = n.generators.index[n.generators.carrier.str.contains("import")]
    limit = config["sector"]["imports"]["limit"]
    if (links.empty and gens.empty) or not np.isfinite(limit):
        return {}
    nyears = n.snapshot_weightings.generators.sum() / 8760
    sense = config["sector"]["imports"]["limit_sense"]
    return {
        "import_weight": pd.concat(
            [pd.Series(1.0, index=gens), n.links.efficiency[links]]
        ),
        f"import_limit_{BOUNDS['==' if sense == '=' else sense]}": float(
            limit * 1e6 * nyears
        ),
    }


def retrofit_gas_boiler_constraint(n: pypsa.Network) -> Data:
    """
    ``solve_network.add_retrofit_gas_boiler_constraint``.

    Each existing gas boiler and the H2 boiler at its position (linopy aligns
    the two by position) follow the heat load at the gas boiler's ``bus1``,
    over its peak in each period, times the gas boiler's ``p_nom``.
    """
    gas = n.links.index[
        n.links.carrier.str.contains("gas boiler") & ~n.links.p_nom_extendable
    ]
    h2 = n.links.index[n.links.carrier.str.contains("retrofitted H2 boiler")]
    if gas.empty or h2.empty:
        return {}
    source = "solve_network.add_retrofit_gas_boiler_constraint"
    pairs = _positional(gas, h2, source)
    columns = n.loads_t.p_set.columns
    heat = columns[
        columns.str.contains("heat")
        & ~columns.str.contains("industry")
        & ~columns.str.contains("agriculture")
    ]
    loads = pd.Series(heat, index=n.loads.bus[heat])
    if loads.index.duplicated().any():
        raise UnsupportedConstraintError(f"{source}: a bus carries several heat loads.")
    following = [
        (boiler, loads[n.links.bus1[boiler]])
        for boiler, _ in pairs
        if n.links.bus1[boiler] in loads.index
    ]
    p_set = n.loads_t.p_set[sorted({load for _, load in following})]
    profile = (
        (p_set / p_set.groupby(level=0).transform("max")).replace(np.inf, 0).fillna(0)
    )
    profile = profile.rename_axis(
        index=["period", "snapshot"], columns="entity"
    ).stack()
    return (
        _pairs("gas_retrofit", pairs)
        | _pairs("gas_retrofit_heat_load", following)
        | {"heat_load_profile": profile}
    )


def land_use_data(n: pypsa.Network) -> Data:
    """The limits of the ``tech_capacity_expansion_limit`` global constraints of ``add_land_use_constraint``."""
    gens = n.generators[
        n.generators.p_nom_extendable
        & np.isfinite(n.generators.p_nom_max)
        & n.generators.build_year.isin(n.investment_periods)
    ]
    if gens.empty:
        return {}
    one_ports = _one_port_buses(n)
    standing = _standing(n).loc[_extendable(n).intersection(one_ports.index)]
    keys = [gens.index, gens.bus, gens.build_year.astype(int)]
    return {
        "land_use_limit": gens.p_nom_max.set_axis(
            pd.MultiIndex.from_arrays(keys, names=["entity", "bus", "period"])
        ),
        "land_use_weight": _at_bus(_indicator(standing), one_ports),
    }


def growth_data(n: pypsa.Network) -> Data:
    """PyPSA's ``define_growth_limit``: an extendable entity counts in the first period it stands in."""
    limited = n.carriers.index[np.isfinite(n.carriers.max_growth)]
    extendable = _standing(n).loc[_extendable(n)]
    carriers = _entity_carriers(n).reindex(extendable.index)
    extendable = extendable[carriers.isin(limited)]
    if extendable.empty:
        return {}
    return {"growth_first": _indicator(extendable.cumsum(axis=1) == 1)}


def land_use_check(n: pypsa.Network) -> None:
    """
    Fail unless the ``tech_capacity_expansion_limit`` global constraints are the limits ``add_land_use_constraint`` sets.

    The math states the land-use limit from ``p_nom_max`` itself, so a
    constraint of this type with any other limit would be one it misses.
    """
    gens = n.generators[n.generators.p_nom_extendable]
    limits = gens.groupby(["carrier", "bus", "build_year"]).p_nom_max.min()
    limits = limits[np.isfinite(limits)].rename(None)
    glcs = n.global_constraints[n.global_constraints.type == LAND_USE]
    keys = [glcs.carrier_attribute, glcs.bus, glcs.investment_period.astype(int)]
    stated = glcs.constant.set_axis(
        pd.MultiIndex.from_arrays(keys, names=limits.index.names)
    ).rename(None)
    if (glcs.sense != "<=").any() or not stated.sort_index().equals(
        limits.sort_index()
    ):
        raise UnsupportedConstraintError(
            "The tech_capacity_expansion_limit global constraints are not the "
            "land-use limits that add_land_use_constraint derives from p_nom_max."
        )


def growth_check(n: pypsa.Network) -> None:
    """Fail where PyPSA clips a negative relative growth that the math reads as is."""
    negative = n.carriers.index[n.carriers.max_relative_growth < 0]
    if not negative.empty:
        raise UnsupportedConstraintError(
            f"Carriers {list(negative)} have a negative max_relative_growth, which PyPSA clips to 0."
        )


def _custom_extra_functionality_is_empty(path: str) -> bool:
    """Whether the custom function only holds a docstring and ``pass``, as the shipped one."""
    tree = ast.parse(Path(path).read_text())
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == Path(path).stem
    )
    return all(
        isinstance(statement, ast.Pass)
        or (
            isinstance(statement, ast.Expr)
            and isinstance(statement.value, ast.Constant)
        )
        for statement in function.body
    )


def _merge(parts: Iterable[Data]) -> Data:
    merged: Data = {}
    for part in parts:
        for name, value in part.items():
            if name not in merged:
                merged[name] = value
            elif isinstance(value, pd.Series):
                merged[name] = pd.concat([merged[name], value], verify_integrity=True)
            else:
                raise ValueError(f"Two constraints set '{name}'.")
    return merged


def extra_functionality_data(
    n: pypsa.Network, config: Mapping[str, Any], custom_extra_functionality: str | None
) -> Data:
    """The data of what ``solve_network`` adds, under the same config switches."""
    _require_multi_period(n)
    if custom_extra_functionality and not _custom_extra_functionality_is_empty(
        custom_extra_functionality
    ):
        raise UnsupportedConstraintError(
            f"Custom extra functionality {custom_extra_functionality} has no math to export."
        )
    land_use_check(n)
    growth_check(n)
    constraints = config["solving"]["constraints"]
    extendable_generators = n.generators.p_nom_extendable.any()
    switched: list[tuple[bool, Callable[[], Data]]] = [
        (
            constraints["BAU"] and extendable_generators,
            lambda: bau_constraints(n, config),
        ),
        (
            constraints["SAFE"] and extendable_generators,
            lambda: safe_constraints(n, config),
        ),
        (
            constraints["CCL"] and extendable_generators,
            lambda: ccl_constraints(n, config),
        ),
        (
            config["electricity"]["operational_reserve"]["activate"],
            lambda: operational_reserve_margin(n, config),
        ),
        (
            bool(constraints["EQ"]),
            lambda: eq_constraints(n, constraints["EQ"].replace("EQ", "")),
        ),
        (
            {"solar-hsat", "solar"}.issubset(
                config["electricity"]["renewable_carriers"]
            )
            and {"solar-hsat", "solar"}.issubset(
                config["electricity"]["extendable_carriers"]["Generator"]
            ),
            lambda: solar_potential_constraints(n, config),
        ),
        (True, lambda: land_use_data(n)),
        (True, lambda: growth_data(n)),
        (True, lambda: tes_energy_to_power_ratio_constraints(n)),
        (True, lambda: tes_charger_ratio_constraints(n)),
        (True, lambda: battery_constraints(n)),
        (True, lambda: lossy_bidirectional_link_constraints(n)),
        (True, lambda: pipe_retrofit_constraint(n, config)),
        (True, lambda: retrofit_gas_boiler_constraint(n)),
        (
            config["sector"]["enhanced_geothermal"]["enable"],
            lambda: flexible_egs_constraint(n),
        ),
        (
            config["sector"]["imports"]["enable"],
            lambda: import_limit_constraint(n, config),
        ),
    ]
    return _merge(build() for active, build in switched if active)


def constraint_data(
    n: pypsa.Network, config: Mapping[str, Any], custom_extra_functionality: str | None
) -> Data:
    """The data of every constraint the static math states, for this network and config."""
    return _merge(
        [
            global_constraint_data(n),
            extra_functionality_data(n, config, custom_extra_functionality),
        ]
    )
