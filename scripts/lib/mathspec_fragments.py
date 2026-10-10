# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Mathspec fragments for the constraints that a datarecord does not carry.

The base spec shipped with datarecord (``datarecord.tools.pypsa.MATHSPEC``)
states PyPSA's linear optimal power flow over the names of a record. It has no
global constraints and no custom constraints. This module writes them as two
mathspec fragments that compose onto the base spec with ``mathspec.merge``:

- :func:`global_constraints` translates ``n.global_constraints`` row by row.
- :func:`extra_functionality` translates the constraints that
  ``solve_network.extra_functionality`` adds, under the same config switches,
  and the PyPSA constraints that ``solve_network.prepare_network`` sets up
  through data (land-use limits per bus and carrier, growth limits per carrier).

Each fragment inlines its data (constants, carrier and entity selections) and
reads only names of the base spec or of the record. A constraint that is active
and that this module cannot translate raises :class:`UnsupportedConstraintError`.
"""

import ast
import re
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pypsa
from pypsa.descriptors import nominal_attrs

from scripts.solve_network import emission_stores, get_glc_sense

Piece = dict[str, dict[str, Any]]

SECTIONS = ("variables", "expressions", "constraints")
EXTRA_PARAMETERS = {
    "length": {"dims": ["entity"]},
    "period_weighting_years": {"dims": ["period"]},
}
BUS_CARRIER = {
    "key": "bus",
    "values": "carrier",
    "missing": "absent",
    "description": "the carrier of each bus",
}
EXTENDABLE = "(p_nom_extendable OR s_nom_extendable OR e_nom_extendable)"
ONE_PORTS = ("Generator", "StorageUnit", "Store")
LAST_IN_PERIOD = "position(snapshot, within=snapshot_period[period]) == -1"
PERIOD_YEARS = "at(period_weighting_years, by=snapshot_period[period])"
STANDING_IN_PERIOD = (
    "active AND build_year <= period_year AND build_year + lifetime > period_year"
)
QUOTED = re.compile(r"'[^']*'|\"[^\"]*\"")
IDENTIFIER = re.compile(r"[A-Za-z_]\w*")


class UnsupportedConstraintError(NotImplementedError):
    pass


def _ident(*parts: object) -> str:
    return re.sub(r"\W+", "_", "_".join(str(p) for p in parts)).strip("_")


def _label(value: object) -> str:
    text = str(value)
    if "'" not in text:
        return f"'{text}'"
    if '"' not in text:
        return f'"{text}"'
    raise ValueError(f"Label {text!r} holds both quote characters.")


def _num(value: float) -> str:
    if not np.isfinite(value):
        raise ValueError(f"Cannot inline the non-finite value {value}.")
    return repr(float(value))


def _any_of(name: str, labels: Iterable[object]) -> str:
    labels = list(dict.fromkeys(labels))
    if not labels:
        raise ValueError(f"No labels to select '{name}' by.")
    return "(" + " OR ".join(f"{name} == {_label(v)}" for v in labels) + ")"


def _is(*types: str) -> str:
    return _any_of("entity_entity_type", types)


def _in_period(period: int) -> str:
    return f"snapshot_period[period] == {int(period)}"


def _standing(period: int) -> str:
    return (
        f"active AND build_year <= {int(period)} "
        f"AND build_year + lifetime > {int(period)}"
    )


def _standing_in(periods: Iterable[int]) -> str:
    return "(" + " OR ".join(f"({_standing(p)})" for p in periods) + ")"


def _sense(sense: str) -> str:
    return {"=": "==", "==": "==", "<=": "<=", ">=": ">="}[sense]


def _single(dims: list[str], when: str, expression: object) -> dict[str, Any]:
    return {
        "dims": dims,
        "cases": {"selected": {"when": when, "expression": str(expression)}},
        "otherwise": "0",
    }


def _cases(dims: list[str], cases: Mapping[str, tuple[str, object]]) -> dict[str, Any]:
    return {
        "dims": dims,
        "cases": {
            key: {"when": when, "expression": str(expression)}
            for key, (when, expression) in cases.items()
        },
        "otherwise": "0",
    }


def _constraint(dims: list[str], expression: str, description: str, **kw: str) -> dict:
    return {"dims": dims, "expression": expression, "description": description} | kw


def _merge(pieces: Iterable[Piece]) -> Piece:
    merged: Piece = {section: {} for section in SECTIONS}
    for piece in pieces:
        for section, entries in piece.items():
            for key, entry in entries.items():
                if merged[section].get(key, entry) != entry:
                    raise ValueError(f"Two pieces declare the {section[:-1]} '{key}'.")
                merged[section][key] = entry
    return {section: entries for section, entries in merged.items() if entries}


def _strings(value: Any) -> Iterable[str]:
    if isinstance(value, str):
        yield value
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from _strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from _strings(item)


def _given(base: Mapping[str, Any], body: Piece) -> dict[str, dict]:
    names = {
        token
        for text in _strings(body)
        for token in IDENTIFIER.findall(QUOTED.sub("", text))
    }
    declared = {key for entries in body.values() for key in entries}
    parameters = base["parameters"] | EXTRA_PARAMETERS
    given = {
        "parameters": {
            name: {k: v for k, v in parameters[name].items() if k in ("dims", "dtype")}
            for name in sorted(names & parameters.keys() - declared)
        },
        "variables": {
            name: {"dims": base["variables"][name]["dims"]}
            for name in sorted(names & base["variables"].keys() - declared)
        },
        "expressions": {
            name: {"dims": base["expressions"][name]["dims"]}
            for name in sorted(names & base["expressions"].keys() - declared)
        },
    }
    return {kind: entries for kind, entries in given.items() if entries}


def fragment(
    base: Mapping[str, Any], description: str, pieces: Iterable[Piece]
) -> dict[str, Any]:
    """Assemble pieces into a fragment that loads on its own and merges onto ``base``."""
    body = _merge(pieces)
    return {
        "description": description,
        "given": _given(base, body),
        "dimensions": base["dimensions"],
        "relations": base["relations"] | {"bus_carrier": BUS_CARRIER},
        **body,
    }


def _require_multi_period(n: pypsa.Network) -> None:
    if n.investment_periods.empty:
        raise ValueError("Mathspec fragments need a network with investment periods.")
    if n.has_scenarios:
        raise UnsupportedConstraintError(
            "Mathspec fragments do not cover stochastic networks."
        )


def _period_of(glc: pd.Series) -> int | None:
    period = glc["investment_period"]
    return None if pd.isna(period) else int(period)


def _require_period(name: str, glc: pd.Series) -> int:
    period = _period_of(glc)
    if period is None:
        raise ValueError(
            f"Global constraint '{name}' of type {glc.type} needs a period."
        )
    return period


def _active_in(c: pypsa.Components, index: pd.Index, period: int | None) -> pd.Index:
    return c.filter_by_active_assets(index, period)


GENERATOR_EFFICIENCY = {
    "generator_port": _single(
        ["entity", "bus"], f"{_is('Generator')} AND role == 'attached'", 1
    ),
    "generator_efficiency": {
        "dims": ["scenario", "snapshot", "entity"],
        "expression": "sum(efficiency * generator_port, over=bus)",
        "description": "the efficiency of a generator, read at its bus",
    },
}
STORAGE_ATTRIBUTES = {
    "StorageUnit": (
        "cyclic_state_of_charge",
        "state_of_charge_initial_per_period",
        "state_of_charge_initial",
    ),
    "Store": ("e_cyclic", "e_initial_per_period", "e_initial"),
}


def _storage_terms(
    n: pypsa.Network,
    c: pypsa.Components,
    weights: pd.Series,
    period: int | None,
    primary: bool,
) -> tuple[dict[str, tuple[str, str]], float]:
    """
    The final-level terms of non-cyclic storage in a PyPSA energy limit.

    Mirrors ``define_primary_energy_limit`` and ``define_operational_limit``
    for a multi-period network. Returns the cases of the weight on the level
    ``e`` and the constant that PyPSA adds to the left-hand side: the initial
    level, weighted as the final one.
    """
    cyclic, per_period, initial = STORAGE_ATTRIBUTES[c.name]
    static = c.static
    selected = static.index[static.carrier.isin(weights.index) & ~static[cyclic]]
    selected = _active_in(c, selected, period if primary else None)
    if selected.empty:
        return {}, 0.0
    static = static.loc[selected]
    years = n.investment_period_weightings.years
    if period is not None:
        years = years.loc[[period]]
    continuous = static.index[~static[per_period]]
    if not continuous.empty and years.ne(1).any():
        raise NotImplementedError(
            f"Found non-cyclic {c.name} components {list(continuous)} with continuous "
            "depletion over multiple investment periods combined with investment "
            "period year weightings != 1, which PyPSA refuses."
        )
    initial_levels = static.carrier.map(weights) * static[initial]
    constant = float(initial_levels[~static[per_period]].sum())
    constant += float(initial_levels[static[per_period]].sum() * years.sum())
    final = "position(snapshot) == -1" if period is None else LAST_IN_PERIOD
    selection = "" if period is None else f" AND {_in_period(period)}"
    cases = {}
    for i, (carrier, value) in enumerate(weights.items()):
        kind = f"{_is(c.name)} AND entity_carrier == {_label(carrier)} AND NOT {cyclic}"
        cases[f"continuous_{i}"] = (
            f"{kind} AND NOT {per_period} AND {final}{selection}",
            f"-{_num(value)}",
        )
        cases[f"per_period_{i}"] = (
            f"{kind} AND {per_period} AND {LAST_IN_PERIOD}{selection}",
            f"-{_num(value)} * {PERIOD_YEARS}",
        )
    return cases, constant


def _energy_limit(
    n: pypsa.Network, name: str, glc: pd.Series, weights: pd.Series, primary: bool
) -> Piece:
    """``define_primary_energy_limit`` (``primary``) or ``define_operational_limit``."""
    period = _period_of(glc)
    if period is not None and period not in n.investment_periods:
        return {}
    key = _ident("GlobalConstraint", name)
    selection = "" if period is None else f" AND {_in_period(period)}"
    divisor = " / generator_efficiency" if primary else ""
    gens = n.generators.index[n.generators.carrier.isin(weights.index)]
    gens = gens.intersection(n.c.generators.active_assets)
    expressions: dict[str, Any] = {}
    lhs, constant = [], 0.0
    if not gens.empty:
        expressions[f"{key}_generator"] = _cases(
            ["scenario", "snapshot", "entity"],
            {
                f"carrier_{i}": (
                    f"{_is('Generator')} AND entity_carrier == {_label(carrier)}{selection}",
                    f"{_num(value)} * snapshot_weighting_generators * {PERIOD_YEARS}{divisor}",
                )
                for i, (carrier, value) in enumerate(weights.items())
            },
        )
        lhs.append(f"sum(p * {key}_generator, over=[snapshot, entity])")
        if primary:
            expressions |= GENERATOR_EFFICIENCY
    for c in (n.c.storage_units, n.c.stores):
        cases, initial = _storage_terms(n, c, weights, period, primary)
        if not cases:
            continue
        storage_key = f"{key}_{_ident(c.name)}"
        expressions[storage_key] = _cases(["scenario", "snapshot", "entity"], cases)
        lhs.append(f"sum(e * {storage_key}, over=[snapshot, entity])")
        constant += initial
    if not lhs:
        return {}
    return {
        "expressions": expressions,
        "constraints": {
            key: _constraint(
                ["scenario"],
                f"{' + '.join(lhs)} {_sense(get_glc_sense(glc))} "
                f"{_num(glc.constant - constant)}",
                f"GlobalConstraint '{name}' of type {glc.type}, as PyPSA builds it.",
            )
        },
    }


def primary_energy(n: pypsa.Network, name: str, glc: pd.Series) -> Piece:
    weights = n.carriers[glc.carrier_attribute].dropna()
    weights = weights[weights != 0]
    if weights.empty:
        return {}
    return _energy_limit(n, name, glc, weights, primary=True)


def operational_limit(n: pypsa.Network, name: str, glc: pd.Series) -> Piece:
    weights = pd.Series(1.0, index=[glc.carrier_attribute])
    return _energy_limit(n, name, glc, weights, primary=False)


def _co2_stores(n: pypsa.Network, glc: pd.Series) -> pd.Index:
    emissions = n.carriers[glc.carrier_attribute]
    emissions = emissions[emissions != 0]
    if emissions.empty:
        return pd.Index([])
    return emission_stores(n, emissions.index).index


def co2_atmosphere(n: pypsa.Network, name: str, glc: pd.Series) -> Piece:
    """``solve_network.add_carbon_constraint``: the emissions of one period."""
    stores = _co2_stores(n, glc)
    if stores.empty:
        return {}
    period = _require_period(name, glc)
    key = _ident("GlobalConstraint", name)
    store = f"{_any_of('entity', stores)} AND {LAST_IN_PERIOD}"
    cases = {"final": (f"{store} AND {_in_period(period)}", 1)}
    periods = list(n.investment_periods)
    if periods.index(period) > 0:
        previous = periods[periods.index(period) - 1]
        cases["previous"] = (f"{store} AND {_in_period(previous)}", -1)
    return {
        "expressions": {f"{key}_level": _cases(["snapshot", "entity"], cases)},
        "constraints": {
            key: _constraint(
                ["scenario"],
                f"sum(e * {key}_level, over=[snapshot, entity]) "
                f"{_sense(get_glc_sense(glc))} {_num(glc.constant)}",
                f"GlobalConstraint '{name}' (co2_atmosphere), as "
                "solve_network.add_carbon_constraint builds it: the rise of the "
                f"atmospheric CO2 level over the period {period}.",
            )
        },
    }


def co2_budget(n: pypsa.Network, name: str, glc: pd.Series) -> Piece:
    """``solve_network.add_carbon_budget_constraint``."""
    stores = _co2_stores(n, glc)
    if stores.empty:
        return {}
    period = _require_period(name, glc)
    key = _ident("GlobalConstraint", name)
    weighting = n.investment_period_weightings.loc[period, "years"]
    when = f"{_any_of('entity', stores)} AND {LAST_IN_PERIOD} AND {_in_period(period)}"
    return {
        "expressions": {
            f"{key}_level": _single(["snapshot", "entity"], when, _num(weighting))
        },
        "constraints": {
            key: _constraint(
                ["scenario"],
                f"sum(e * {key}_level, over=[snapshot, entity]) "
                f"{_sense(get_glc_sense(glc))} {_num(glc.constant)}",
                f"GlobalConstraint '{name}' (Co2Budget), as "
                "solve_network.add_carbon_budget_constraint builds it.",
            )
        },
    }


def _transmission_carriers(glc: pd.Series) -> list[str]:
    return [
        re.sub(r"[\[\]\(\)]", "", c.strip()) for c in glc.carrier_attribute.split(",")
    ]


def _transmission(n: pypsa.Network, name: str, glc: pd.Series, cost: bool) -> Piece:
    """``define_transmission_volume_expansion_limit`` or ``..._cost_limit``."""
    carriers = _transmission_carriers(glc)
    period = _period_of(glc)
    periods = list(n.investment_periods) if period is None else [period]
    selected = [
        _active_in(
            c,
            c.extendables.intersection(c.static.index[c.static.carrier.isin(carriers)]),
            periods,
        )
        for c in n.components[["Line", "Link"]]
    ]
    if all(s.empty for s in selected):
        return {}
    key = _ident("GlobalConstraint", name)
    if not cost:
        weight = "length"
    elif period is None:
        weight = "capital_cost * sum(capital_weight, over=period)"
    else:
        weight = "capital_cost"
    when = (
        f"{_is('Line', 'Link')} AND {EXTENDABLE} "
        f"AND {_any_of('entity_carrier', carriers)} AND {_standing_in(periods)}"
    )
    return {
        "expressions": {f"{key}_weight": _single(["entity"], when, weight)},
        "constraints": {
            key: _constraint(
                [],
                f"sum(nom * {key}_weight, over=entity) "
                f"{_sense(get_glc_sense(glc))} {_num(glc.constant)}",
                f"GlobalConstraint '{name}' of type {glc.type}, as PyPSA builds it.",
            )
        },
    }


def transmission_volume_expansion_limit(
    n: pypsa.Network, name: str, glc: pd.Series
) -> Piece:
    return _transmission(n, name, glc, cost=False)


def transmission_expansion_cost_limit(
    n: pypsa.Network, name: str, glc: pd.Series
) -> Piece:
    return _transmission(n, name, glc, cost=True)


def tech_capacity_expansion_limit(n: pypsa.Network, name: str, glc: pd.Series) -> Piece:
    """``define_tech_capacity_expansion_limit``."""
    carrier = glc.carrier_attribute
    period = _period_of(glc)
    standing = "active" if period is None else _standing(period)
    selected = [
        _active_in(
            c,
            c.extendables.intersection(c.static.index[c.static.carrier == carrier]),
            period,
        )
        for c in n.components
        if c.name in nominal_attrs and "carrier" in c.static
    ]
    if all(s.empty for s in selected):
        return {}
    key = _ident("GlobalConstraint", name)
    when = f"{EXTENDABLE} AND entity_carrier == {_label(carrier)} AND {standing}"
    bus = glc["bus"] if "bus" in glc.index else ""
    if pd.isna(bus) or bus == "":
        weight = _single(["entity"], when, 1)
        lhs = f"sum(nom * {key}_weight, over=entity)"
    else:
        port = "(role == 'input' OR role == 'attached')"
        weight = _single(
            ["entity", "bus"], f"{when} AND {port} AND bus == {_label(bus)}", 1
        )
        lhs = f"sum(nom * {key}_weight, over=[entity, bus])"
    return {
        "expressions": {f"{key}_weight": weight},
        "constraints": {
            key: _constraint(
                [],
                f"{lhs} {_sense(get_glc_sense(glc))} {_num(glc.constant)}",
                f"GlobalConstraint '{name}' of type tech_capacity_expansion_limit.",
            )
        },
    }


GLOBAL_CONSTRAINTS: dict[str, Callable[[pypsa.Network, str, pd.Series], Piece]] = {
    "primary_energy": primary_energy,
    "operational_limit": operational_limit,
    "transmission_volume_expansion_limit": transmission_volume_expansion_limit,
    "transmission_expansion_cost_limit": transmission_expansion_cost_limit,
    "tech_capacity_expansion_limit": tech_capacity_expansion_limit,
    "co2_atmosphere": co2_atmosphere,
    "Co2Budget": co2_budget,
}


def global_constraints(n: pypsa.Network, base: Mapping[str, Any]) -> dict[str, Any]:
    """Translate every row of ``n.global_constraints`` into one fragment."""
    _require_multi_period(n)
    unknown = set(n.global_constraints.type) - GLOBAL_CONSTRAINTS.keys()
    if unknown:
        raise UnsupportedConstraintError(
            f"No mathspec translation for global constraint types {sorted(unknown)}."
        )
    pieces = [
        GLOBAL_CONSTRAINTS[glc.type](n, name, glc)
        for name, glc in n.global_constraints.iterrows()
    ]
    return fragment(
        base,
        "PyPSA-Eur global constraints (n.global_constraints), with data inlined.",
        pieces,
    )


def _pairs(
    name: str, rows: list[tuple[dict[str, float], str, float]], source: str
) -> Piece:
    """One scalar constraint per pair of entities: sum of coefficient * nom."""
    expressions, constraints = {}, {}
    for i, (coefficients, sense, constant) in enumerate(rows):
        key = f"{name}_{i}"
        expressions[f"{key}_coefficient"] = _cases(
            ["entity"],
            {
                f"entity_{j}": (f"entity == {_label(entity)}", _num(value))
                for j, (entity, value) in enumerate(coefficients.items())
            },
        )
        constraints[key] = _constraint(
            [],
            f"sum(nom * {key}_coefficient, over=entity) {sense} {_num(constant)}",
            f"{source} for {', '.join(coefficients)}.",
        )
    return {"expressions": expressions, "constraints": constraints}


def bau_constraints(n: pypsa.Network, config: Mapping[str, Any]) -> Piece:
    """``solve_network.add_BAU_constraints``."""
    mincaps = pd.Series(config["electricity"]["BAU_mincapacities"])
    carriers = n.generators.query("p_nom_extendable").carrier.unique()
    rhs = mincaps[carriers]
    return {
        "expressions": {
            "bau_mincaps_generator": _single(
                ["entity"], f"{_is('Generator')} AND p_nom_extendable", 1
            ),
            "bau_mincaps_rhs": _cases(
                ["carrier"],
                {
                    f"carrier_{i}": (f"carrier == {_label(c)}", _num(v))
                    for i, (c, v) in enumerate(rhs.items())
                },
            ),
        },
        "constraints": {
            "bau_mincaps": _constraint(
                ["carrier"],
                "sum(nom * bau_mincaps_generator, over=entity, by=entity_carrier[carrier]) >= bau_mincaps_rhs",
                "solve_network.add_BAU_constraints: minimum extendable generator capacity per carrier.",
                where=_any_of("carrier", rhs.index),
            )
        },
    }


def safe_constraints(n: pypsa.Network, config: Mapping[str, Any]) -> Piece:
    """``solve_network.add_SAFE_constraints``."""
    peakdemand = n.loads_t.p_set.sum(axis=1).max()
    reserve_margin = peakdemand * (1.0 + config["electricity"]["SAFE_reservemargin"])
    conventional = config["electricity"]["conventional_carriers"]
    existing = n.generators.query(
        "~p_nom_extendable & carrier in @conventional"
    ).p_nom.sum()
    when = f"{_is('Generator')} AND p_nom_extendable AND {_any_of('entity_carrier', conventional)}"
    return {
        "expressions": {"safe_mintotalcap_generator": _single(["entity"], when, 1)},
        "constraints": {
            "safe_mintotalcap": _constraint(
                [],
                f"sum(nom * safe_mintotalcap_generator, over=entity) >= {_num(reserve_margin - existing)}",
                "solve_network.add_SAFE_constraints: conventional capacity above the peak demand.",
            )
        },
    }


def ccl_constraints(n: pypsa.Network, config: Mapping[str, Any]) -> Piece:
    raise UnsupportedConstraintError(
        "solve_network.add_CCL_constraints is not implemented for perfect foresight "
        "(solve_network asserts a planning horizon), so it has no mathspec translation."
    )


def operational_reserve_margin(n: pypsa.Network, config: Mapping[str, Any]) -> Piece:
    """``solve_network.add_operational_reserve_margin``."""
    reserve = config["electricity"]["operational_reserve"]
    vres = [
        "ror" if c == "hydro" else c
        for c in config["electricity"]["renewable_carriers"]
    ]
    gens = n.generators.query("carrier != 'load'")
    vres_i = gens.index[gens.carrier.isin(vres)]
    no_profile = vres_i.difference(n.generators_t.p_max_pu.columns)
    if not no_profile.empty:
        raise KeyError(
            f"Renewable generators {list(no_profile)} have no p_max_pu profile."
        )
    generator = f"{_is('Generator')} AND NOT entity_carrier == 'load'"
    vres_gen = f"{generator} AND {_any_of('entity_carrier', vres)}"
    rhs = (
        f"{_num(reserve['epsilon_load'])} * sum(p_set * reserve_margin_load, over=entity) "
        f"+ {_num(reserve['epsilon_vres'])} * sum(p_max_pu * p_nom * reserve_margin_vres_fixed, over=entity) "
        f"+ {_num(reserve['contingency'])}"
    )
    lhs = (
        "sum(Generator_r, over=entity) "
        f"- {_num(reserve['epsilon_vres'])} * sum(nom * p_max_pu * reserve_margin_vres_extendable, over=entity)"
    )
    description = "solve_network.add_operational_reserve_margin"
    return {
        "variables": {
            "Generator_r": {
                "dims": ["scenario", "snapshot", "entity"],
                "where": generator,
                "bounds": {"lower": 0},
                "description": "reserve a generator holds",
            }
        },
        "expressions": {
            "reserve_margin_load": _single(["entity"], _is("Load"), 1),
            "reserve_margin_vres_fixed": _single(
                ["entity"], f"{vres_gen} AND NOT p_nom_extendable", 1
            ),
            "reserve_margin_vres_extendable": _single(
                ["entity"], f"{vres_gen} AND p_nom_extendable", 1
            ),
        },
        "constraints": {
            "reserve_margin": _constraint(
                ["scenario", "snapshot"],
                f"{lhs} >= {rhs}",
                f"{description}: the reserve requirement.",
            ),
            "Generator_p_reserve_upper_extendable": _constraint(
                ["scenario", "snapshot", "entity"],
                "p + Generator_r - p_max_pu * nom <= 0",
                f"{description}: dispatch and reserve within extendable capacity.",
                where="Generator_r AND p AND p_nom_extendable",
            ),
            "Generator_p_reserve_upper_extendable_idle": _constraint(
                ["scenario", "snapshot", "entity"],
                "Generator_r - p_max_pu * nom <= 0",
                f"{description}: reserve within extendable capacity where the unit does not dispatch.",
                where="Generator_r AND NOT p AND p_nom_extendable",
            ),
            "Generator_p_reserve_upper_fixed": _constraint(
                ["scenario", "snapshot", "entity"],
                "p + Generator_r <= p_max_pu * p_nom",
                f"{description}: dispatch and reserve within fixed capacity.",
                where="Generator_r AND p AND NOT p_nom_extendable",
            ),
            "Generator_p_reserve_upper_fixed_idle": _constraint(
                ["scenario", "snapshot", "entity"],
                "Generator_r <= p_max_pu * p_nom",
                f"{description}: reserve within fixed capacity where the unit does not dispatch.",
                where="Generator_r AND NOT p AND NOT p_nom_extendable",
            ),
        },
    }


def eq_constraints(n: pypsa.Network, option: str, scaling: float = 1e-1) -> Piece:
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
    level = float(re.findall(r"[0-9]*\.?[0-9]+", option)[0])
    lhs = "sum(p * equity_min_generator, over=[snapshot, entity])"
    inflow = "0"
    expressions = {
        "equity_min_generator": _single(
            ["snapshot", "entity", "bus"],
            f"{_is('Generator')} AND role == 'attached'",
            f"snapshot_weighting_generators * {_num(scaling)}",
        ),
        "equity_min_load": _single(
            ["snapshot", "entity", "bus"],
            f"{_is('Load')} AND role == 'attached'",
            "snapshot_weighting_generators",
        ),
    }
    if not n.storage_units_t.inflow.empty:
        expressions["equity_min_storage_unit"] = _single(
            ["snapshot", "entity", "bus"],
            f"{_is('StorageUnit')} AND role == 'attached'",
            "snapshot_weighting_stores",
        )
        lhs += f" - {_num(scaling)} * sum(spill * equity_min_storage_unit, over=[snapshot, entity])"
        inflow = "sum(inflow * equity_min_storage_unit, over=[snapshot, entity])"
    load = "sum(p_set * equity_min_load, over=[snapshot, entity])"
    return {
        "expressions": expressions,
        "constraints": {
            "equity_min": _constraint(
                ["scenario", "bus"],
                f"{lhs} >= {_num(scaling)} * ({_num(level)} * {load} - {inflow})",
                "solve_network.add_EQ_constraints: each node produces a share of its demand.",
                where=f"count({_is('Load')} AND role == 'attached', over=entity) >= 1",
            )
        },
    }


def solar_potential_constraints(n: pypsa.Network, config: Mapping[str, Any]) -> Piece:
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
    attached = "role == 'attached'"
    extendable = f"{_is('Generator')} AND p_nom_extendable AND {attached}"
    today = f"{extendable} AND entity_carrier == 'solar'"
    hsat = f"{_is('Generator')} AND {attached} AND entity_carrier == 'solar-hsat'"
    expressions = {
        "solar_potential_generator": _cases(
            ["entity", "bus"],
            {
                "solar": (f"{extendable} AND entity_carrier == 'solar'", 1),
                "solar_hsat": (
                    f"{extendable} AND entity_carrier == 'solar-hsat'",
                    _num(factor),
                ),
            },
        ),
        "solar_potential_today": _single(["entity", "bus"], today, 1),
        "solar_potential_hsat": _single(["entity", "bus"], hsat, 1),
        "solar_potential_rhs": {
            "dims": ["bus"],
            "expression": (
                "sum(p_nom_max * solar_potential_today, over=entity) "
                f"- {_num(factor)} * sum(p_nom * solar_potential_hsat, over=entity)"
            ),
            "description": "the solar potential of a bus left after the existing tracking solar",
        },
    }
    rows = f"count({today}, over=entity) >= 1 AND count({hsat}, over=entity) >= 1"
    lhs = "sum(nom * solar_potential_generator, over=entity)"
    description = "solve_network.add_solar_potential_constraints"
    return {
        "expressions": expressions,
        "constraints": {
            "solar_potential": _constraint(
                ["bus"],
                f"{lhs} <= solar_potential_rhs",
                f"{description}.",
                where=f"{rows} AND solar_potential_rhs >= 0",
            ),
            "solar_potential_clipped": _constraint(
                ["bus"],
                f"{lhs} <= 0",
                f"{description}, where the potential left is clipped to zero.",
                where=f"{rows} AND solar_potential_rhs < 0",
            ),
        },
    }


def tes_energy_to_power_ratio_constraints(n: pypsa.Network) -> Piece:
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
    ratios = n.links.loc[chargers, "energy to power ratio"]
    rows = []
    for charger, tes, ratio in zip(chargers, stores, ratios):
        if not tes == charger.replace(" charger", ""):
            raise RuntimeError(f"Charger {charger} and TES {tes} do not match.")
        rows.append(({tes: 1.0, charger: -ratio}, "==", 0.0))
    return _pairs(
        "TES_energy_to_power_ratio",
        rows,
        "solve_network.add_TES_energy_to_power_ratio_constraints",
    )


def tes_charger_ratio_constraints(n: pypsa.Network) -> Piece:
    """``solve_network.add_TES_charger_ratio_constraints``."""
    kinds = "water tanks {0}|water pits {0}|aquifer thermal energy storage {0}"
    chargers = n.links.index[
        n.links.index.str.contains(kinds.format("charger")) & n.links.p_nom_extendable
    ]
    dischargers = n.links.index[
        n.links.index.str.contains(kinds.format("discharger"))
        & n.links.p_nom_extendable
    ]
    if chargers.empty or dischargers.empty:
        return {}
    rows = []
    for charger, discharger in zip(chargers, dischargers):
        if charger.replace(" charger", " ") != discharger.replace(" discharger", " "):
            raise RuntimeError(
                f"Charger {charger} and discharger {discharger} do not match."
            )
        rows.append(
            ({charger: 1.0, discharger: -n.links.efficiency[discharger]}, "==", 0.0)
        )
    return _pairs(
        "TES_charger_ratio", rows, "solve_network.add_TES_charger_ratio_constraints"
    )


def _positional(left: pd.Index, right: pd.Index, source: str) -> list[tuple[str, str]]:
    if len(left) != len(right):
        raise UnsupportedConstraintError(
            f"{source} pairs {len(left)} with {len(right)} entities; linopy aligns "
            "them by an outer join then, which this module does not translate."
        )
    return list(zip(left, right))


def battery_constraints(n: pypsa.Network) -> Piece:
    """``solve_network.add_battery_constraints``."""
    if not n.links.p_nom_extendable.any():
        return {}
    ext = n.links.query("p_nom_extendable")
    dischargers = ext.index[ext.index.str.contains("battery discharger")]
    chargers = ext.index[ext.index.str.contains("battery charger")]
    source = "solve_network.add_battery_constraints"
    rows = [
        ({charger: 1.0, discharger: -n.links.efficiency[discharger]}, "==", 0.0)
        for charger, discharger in _positional(chargers, dischargers, source)
    ]
    return _pairs("Link_charger_ratio", rows, source)


def lossy_bidirectional_link_constraints(n: pypsa.Network) -> Piece:
    """``solve_network.add_lossy_bidirectional_link_constraints``."""
    if (
        not n.links.p_nom_extendable.any()
        or "reversed" not in n.links
        or not n.links.reversed.any()
    ):
        return {}
    carriers = n.links.loc[n.links.reversed, "carrier"].unique()
    links = n.links
    backwards = links.index[
        links.carrier.isin(carriers)
        & links.p_nom_extendable
        & links.reversed
        & links.active
    ]
    forwards = backwards.str.replace("-reversed", "")
    missing = forwards.difference(n.links.index[n.links.p_nom_extendable])
    if not missing.empty:
        raise KeyError(
            f"Forward links {list(missing)} of reversed links are not extendable."
        )
    rows = [({b: 1.0, f: -1.0}, "==", 0.0) for b, f in zip(backwards, forwards)]
    return _pairs(
        "Link_bidirectional_sync",
        rows,
        "solve_network.add_lossy_bidirectional_link_constraints",
    )


def pipe_retrofit_constraint(n: pypsa.Network, config: Mapping[str, Any]) -> Piece:
    """``solve_network.add_pipe_retrofit_constraint``."""
    links = n.links.assign(
        reversed=n.links["reversed"] if "reversed" in n.links else False
    )
    gas = links.query(
        "carrier == 'gas pipeline' and p_nom_extendable and ~reversed and active"
    ).index
    h2 = links.query(
        "carrier == 'H2 pipeline retrofitted' and p_nom_extendable and ~reversed and active"
    ).index
    if gas.empty or h2.empty:
        return {}
    ch4_per_h2 = 1 / config["sector"]["H2_retrofit_capacity_per_CH4"]
    source = "solve_network.add_pipe_retrofit_constraint"
    rows = [
        ({g: 1.0, h: ch4_per_h2}, "==", n.links.p_nom[g])
        for g, h in _positional(gas, h2, source)
    ]
    return _pairs("Link_pipe_retrofit", rows, source)


def flexible_egs_constraint(n: pypsa.Network) -> Piece:
    """``solve_network.add_flexible_egs_constraint``."""
    wells = n.links.index[n.links.carrier == "geothermal heat"]
    storage = n.storage_units.index[n.storage_units.carrier == "geothermal heat"]
    source = "solve_network.add_flexible_egs_constraint"
    rows = [
        ({s: 1.0, w: -1.0}, "<=", 0.0) for s, w in _positional(storage, wells, source)
    ]
    return _pairs("upper_bound_charging_capacity_of_geothermal_reservoir", rows, source)


def import_limit_constraint(n: pypsa.Network, config: Mapping[str, Any]) -> Piece:
    """``solve_network.add_import_limit_constraint``."""
    nyears = n.snapshot_weightings.generators.sum() / 8760
    links = n.links[n.links.carrier.str.contains("import")]
    gens = n.generators[n.generators.carrier.str.contains("import")]
    limit = config["sector"]["imports"]["limit"]
    sense = _sense(config["sector"]["imports"]["limit_sense"])
    if (links.empty and gens.empty) or not np.isfinite(limit):
        return {}
    efficiency = links.groupby(["carrier", "build_year"]).efficiency
    if efficiency.nunique().gt(1).any():
        raise UnsupportedConstraintError(
            "add_import_limit_constraint weights import links by their static "
            "efficiency, which differs within one carrier and build year."
        )
    cases = {
        f"generator_{i}": (
            f"{_is('Generator')} AND entity_carrier == {_label(c)}",
            "snapshot_weighting_generators",
        )
        for i, c in enumerate(gens.carrier.unique())
    }
    for i, ((carrier, year), eff) in enumerate(efficiency.first().items()):
        cases[f"link_{i}"] = (
            f"{_is('Link')} AND entity_carrier == {_label(carrier)} AND build_year == {int(year)}",
            f"{_num(eff)} * snapshot_weighting_generators",
        )
    return {
        "expressions": {"import_limit_weight": _cases(["snapshot", "entity"], cases)},
        "constraints": {
            "import_limit": _constraint(
                ["scenario"],
                f"sum(p * import_limit_weight, over=[snapshot, entity]) {sense} {_num(limit * 1e6 * nyears)}",
                "solve_network.add_import_limit_constraint: green energy imports.",
            )
        },
    }


def retrofit_gas_boiler_constraint(n: pypsa.Network) -> Piece:
    """
    ``solve_network.add_retrofit_gas_boiler_constraint``.

    Each existing gas boiler and the H2 boiler at its position (linopy aligns
    the two by position) follow the heat load at the gas boiler's ``bus1``,
    scaled by its peak in each period, times the gas boiler's ``p_nom``. The
    peak is data, so it is inlined per load and period.
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
    peak = n.loads_t.p_set[heat].groupby(level=0).max()
    expressions, constraints = {}, {}
    for i, (boiler, retrofit) in enumerate(pairs):
        bus = n.links.bus1[boiler]
        if bus not in loads.index:
            continue
        load = loads[bus]
        static = n.links.loc[boiler]
        periods = [
            p
            for p in n.investment_periods
            if static.active
            and static.build_year <= p < static.build_year + static.lifetime
        ]
        if not periods:
            continue
        key = f"gas_retrofit_{i}"
        expressions[f"{key}_coefficient"] = _cases(
            ["entity"],
            {
                "gas": (f"entity == {_label(boiler)}", 1),
                "h2": (f"entity == {_label(retrofit)}", 1),
            },
        )
        scale = (static.p_nom / peak[load]).replace(np.inf, 0).fillna(0)
        expressions[f"{key}_profile"] = _cases(
            ["snapshot", "entity"],
            {
                f"period_{p}": (
                    f"entity == {_label(load)} AND {_in_period(p)}",
                    _num(scale[p]),
                )
                for p in periods
            },
        )
        constraints[key] = _constraint(
            ["scenario", "snapshot"],
            f"sum(p * {key}_coefficient, over=entity) == sum(p_set * {key}_profile, over=entity)",
            f"{source} for {boiler} and {retrofit}, following {load}.",
            where="(" + " OR ".join(_in_period(p) for p in periods) + ")",
        )
    return {"expressions": expressions, "constraints": constraints}


def land_use_constraint(n: pypsa.Network) -> Piece:
    """
    ``solve_network.add_land_use_constraint``, through the PyPSA bus columns it sets.

    PyPSA reads a bus column ``nom_max_{carrier}_{period}`` as an upper limit on
    the extendable capacity of that carrier at that bus, among the assets that
    stand in that period. The land-use constraint sets it to the smallest
    ``p_nom_max`` of the extendable generators of that carrier, bus and build
    year. "At most the smallest" is "at most each", which is how this writes it.
    """
    columns = n.buses.columns[n.buses.columns.str.startswith("nom_")]
    gens = n.generators[n.generators.p_nom_extendable]
    limits = gens.groupby(["carrier", "bus", "build_year"]).p_nom_max.min()
    limits = limits[np.isfinite(limits)]
    pieces = []
    for column in columns:
        if not re.fullmatch(r"nom_max_.+_\d{4}", column):
            raise UnsupportedConstraintError(
                f"Bus column '{column}' is a PyPSA capacity limit per bus and carrier "
                "that add_land_use_constraint does not set."
            )
        carrier, period = column.removeprefix("nom_max_").rsplit("_", 1)
        actual = n.buses[column].dropna().sort_index()
        expected = limits.xs((carrier, int(period)), level=["carrier", "build_year"])
        if not actual.equals(
            expected.sort_index().rename(column).rename_axis(actual.index.name)
        ):
            raise UnsupportedConstraintError(
                f"Bus column '{column}' is not the land-use limit that "
                "add_land_use_constraint derives from p_nom_max."
            )
        key = _ident("land_use_constraint", carrier, period)
        selected = f"{EXTENDABLE} AND entity_carrier == {_label(carrier)} AND role == 'attached'"
        pieces.append(
            {
                "expressions": {
                    f"{key}_weight": _single(
                        ["entity", "bus"],
                        f"{_is(*ONE_PORTS)} AND {selected} AND {_standing(int(period))}",
                        1,
                    ),
                    f"{key}_capacity": {
                        "dims": ["bus"],
                        "expression": f"sum(nom * {key}_weight, over=entity)",
                    },
                },
                "constraints": {
                    key: _constraint(
                        ["entity", "bus"],
                        f"{key}_capacity <= p_nom_max",
                        f"solve_network.add_land_use_constraint (Bus-{column}): extendable "
                        f"{carrier} at a bus in {period}, at most the potential of each "
                        f"generator built in {period}.",
                        where=(
                            f"{_is('Generator')} AND {selected} AND build_year == {int(period)} "
                            "AND p_nom_max < inf"
                        ),
                    )
                },
            }
        )
    return _merge(pieces) if pieces else {}


def max_growth(n: pypsa.Network) -> Piece:
    """PyPSA's ``define_growth_limit``, from the carrier columns ``add_max_growth`` sets."""
    growth = n.carriers.max_growth[n.carriers.max_growth != np.inf]
    relative = n.carriers.max_relative_growth.clip(lower=0)
    expressions, constraints = {}, {}
    first = f"({STANDING_IN_PERIOD}) AND NOT shift({STANDING_IN_PERIOD}, along=period, offset=1)"
    for carrier, limit in growth.items():
        key = _ident("Carrier_growth_limit", carrier)
        when = f"{EXTENDABLE} AND entity_carrier == {_label(carrier)} AND {first}"
        expressions[f"{key}_weight"] = _single(["entity", "period"], when, 1)
        built = f"sum(nom * {key}_weight, over=entity)"
        lhs = built
        if relative[carrier] > 0:
            lhs += f" - {_num(relative[carrier])} * shift({built}, along=period, offset=1, edge=0)"
        constraints[key] = _constraint(
            ["period"],
            f"{lhs} <= {_num(limit)}",
            f"PyPSA growth limit of {carrier} (solve_network.add_max_growth): capacity "
            "first standing in a period.",
        )
    return (
        {"expressions": expressions, "constraints": constraints} if constraints else {}
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


def extra_functionality(
    n: pypsa.Network,
    config: Mapping[str, Any],
    custom_extra_functionality: str | None,
    base: Mapping[str, Any],
) -> dict[str, Any]:
    """Translate what ``solve_network`` adds beyond PyPSA, under its config switches."""
    _require_multi_period(n)
    if custom_extra_functionality and not _custom_extra_functionality_is_empty(
        custom_extra_functionality
    ):
        raise UnsupportedConstraintError(
            f"Custom extra functionality {custom_extra_functionality} has no math to export."
        )
    constraints = config["solving"]["constraints"]
    extendable_generators = n.generators.p_nom_extendable.any()
    pieces = [land_use_constraint(n), max_growth(n)]
    if constraints["BAU"] and extendable_generators:
        pieces.append(bau_constraints(n, config))
    if constraints["SAFE"] and extendable_generators:
        pieces.append(safe_constraints(n, config))
    if constraints["CCL"] and extendable_generators:
        pieces.append(ccl_constraints(n, config))
    if config["electricity"]["operational_reserve"]["activate"]:
        pieces.append(operational_reserve_margin(n, config))
    if constraints["EQ"]:
        pieces.append(eq_constraints(n, constraints["EQ"].replace("EQ", "")))
    solar = {"solar-hsat", "solar"}
    if solar.issubset(config["electricity"]["renewable_carriers"]) and solar.issubset(
        config["electricity"]["extendable_carriers"]["Generator"]
    ):
        pieces.append(solar_potential_constraints(n, config))
    pieces += [
        tes_energy_to_power_ratio_constraints(n),
        tes_charger_ratio_constraints(n),
        battery_constraints(n),
        lossy_bidirectional_link_constraints(n),
        pipe_retrofit_constraint(n, config),
        retrofit_gas_boiler_constraint(n),
    ]
    if config["sector"]["enhanced_geothermal"]["enable"]:
        pieces.append(flexible_egs_constraint(n))
    if config["sector"]["imports"]["enable"]:
        pieces.append(import_limit_constraint(n, config))
    return fragment(
        base,
        "PyPSA-Eur custom constraints (solve_network.extra_functionality and "
        "prepare_network), with data inlined.",
        pieces,
    )
