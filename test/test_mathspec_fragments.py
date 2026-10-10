# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import copy
from pathlib import Path

import pandas as pd
import pypsa
import pytest
import yaml

from scripts.lib import mathspec_fragments as mf
from scripts.solve_network import add_land_use_constraint, add_max_growth

pytestmark = pytest.mark.datarecord

PERIODS = [2030, 2040]


@pytest.fixture(scope="module")
def base() -> dict:
    from datarecord.tools.pypsa import MATHSPEC

    return yaml.safe_load(MATHSPEC.read_text())


@pytest.fixture(scope="function")
def network() -> pypsa.Network:
    n = pypsa.Network()
    snapshots = pd.date_range("2013-01-01", periods=3, freq="h")
    n.set_snapshots(pd.MultiIndex.from_product([PERIODS, snapshots]))
    n.investment_periods = PERIODS
    n.investment_period_weightings.loc[:, "years"] = 10.0
    n.investment_period_weightings.loc[:, "objective"] = [10.0, 7.0]
    n.add("Carrier", ["AC", "gas", "solar", "solar-hsat", "co2", "battery", "hydro"])
    n.add("Carrier", ["battery charger", "battery discharger", "H2 pipeline", "H2"])
    n.carriers.loc["gas", ["co2_emissions", "gas_usage"]] = [0.2, 1.0]
    n.carriers.loc["co2", "co2_emissions"] = -1.0
    n.add("Bus", ["b0", "b1"], carrier="AC", country="DE")
    n.add("Bus", ["b0 battery", "b0 H2", "b1 H2"], carrier=["battery", "H2", "H2"])
    n.add("Bus", "co2 atmosphere", carrier="co2")
    loads = pd.DataFrame(10.0, n.snapshots, ["b0 load", "b1 load"])
    n.add("Load", loads.columns, bus=["b0", "b1"], carrier="AC", p_set=loads)
    n.add("Generator", "b0 gas", bus="b0", carrier="gas", p_nom=20, efficiency=0.5)
    for year in PERIODS:
        built = dict(p_nom_extendable=True, build_year=year, lifetime=25)
        n.add("Generator", f"b0 gas-{year}", bus="b0", carrier="gas", **built)
        for bus in ("b0", "b1"):
            for carrier, p_nom_max in (("solar", 50), ("solar-hsat", 40)):
                n.add(
                    "Generator",
                    f"{bus} 0 {carrier}-{year}",
                    bus=bus,
                    carrier=carrier,
                    p_nom_max=p_nom_max,
                    p_max_pu=pd.Series(0.5, n.snapshots),
                    **built,
                )
        battery = dict(carrier="battery charger", bus0="b0", bus1="b0 battery")
        n.add("Link", f"b0 battery charger-{year}", efficiency=0.95, **battery, **built)
        battery = dict(carrier="battery discharger", bus0="b0 battery", bus1="b0")
        n.add(
            "Link", f"b0 battery discharger-{year}", efficiency=0.9, **battery, **built
        )
        n.add(
            "Store",
            f"b0 battery-{year}",
            bus="b0 battery",
            carrier="battery",
            e_nom_extendable=True,
            e_cyclic_per_period=True,
            build_year=year,
            lifetime=25,
        )
        pipe = dict(carrier="H2 pipeline", length=100, **built)
        n.add("Link", f"H2 pipeline-{year}", bus0="b0 H2", bus1="b1 H2", **pipe)
        n.add(
            "Link",
            f"H2 pipeline-{year}-reversed",
            bus0="b1 H2",
            bus1="b0 H2",
            reversed=True,
            **pipe,
        )
    n.links["reversed"] = n.links.reversed.fillna(False).astype(bool)
    n.add(
        "Line",
        "b0-b1",
        bus0="b0",
        bus1="b1",
        carrier="AC",
        x=0.1,
        s_nom_extendable=True,
        length=200,
        capital_cost=3,
    )
    n.add(
        "Store",
        "co2 atmosphere",
        bus="co2 atmosphere",
        carrier="co2",
        e_nom=1e6,
        e_min_pu=-1,
    )
    n.add(
        "StorageUnit",
        "b1 hydro",
        bus="b1",
        carrier="hydro",
        p_nom=5,
        max_hours=6,
        state_of_charge_initial=2.0,
        state_of_charge_initial_per_period=True,
        inflow=pd.Series(1.0, n.snapshots),
    )
    for c in n.components[
        ["Generator", "Link", "Line", "Load", "Store", "StorageUnit"]
    ]:
        c.static["technology_group"] = c.static.carrier
    glc = dict(carrier_attribute="co2_emissions", sense="<=")
    n.add(
        "GlobalConstraint",
        "CO2Limit-upper-2040",
        type="co2_atmosphere",
        investment_period=2040,
        constant=500.0,
        **glc,
    )
    n.add(
        "GlobalConstraint",
        "CO2Budget",
        type="Co2Budget",
        investment_period=2040,
        constant=1e4,
        **glc,
    )
    n.add(
        "GlobalConstraint",
        "GasLimit",
        type="primary_energy",
        carrier_attribute="gas_usage",
        sense="<=",
        constant=1e5,
    )
    n.add(
        "GlobalConstraint",
        "lv_limit",
        type="transmission_volume_expansion_limit",
        carrier_attribute="AC, DC",
        sense="<=",
        constant=1e4,
    )
    n.add(
        "GlobalConstraint",
        "lc_limit",
        type="transmission_expansion_cost_limit",
        carrier_attribute="AC, DC",
        sense="<=",
        constant=1e6,
    )
    n.add(
        "GlobalConstraint",
        "hydro limit",
        type="operational_limit",
        carrier_attribute="hydro",
        sense="<=",
        constant=1e3,
        investment_period=2040,
    )
    n.add(
        "GlobalConstraint",
        "solar cap",
        type="tech_capacity_expansion_limit",
        carrier_attribute="solar",
        sense="<=",
        constant=1e3,
        investment_period=2040,
    )
    add_land_use_constraint(n)
    add_max_growth(
        n,
        {
            "factor": 1.3,
            "max_growth": {"solar": 2},
            "max_relative_growth": {"solar": 1.5},
        },
    )
    return n


@pytest.fixture(scope="function")
def solve_config(config: dict) -> dict:
    solve_config = copy.deepcopy(config)
    solve_config["electricity"].update(
        renewable_carriers=["solar", "solar-hsat"],
        conventional_carriers=["gas"],
        BAU_mincapacities={"gas": 1.0, "solar": 0.0, "solar-hsat": 0.0},
        SAFE_reservemargin=0.1,
    )
    solve_config["electricity"]["extendable_carriers"]["Generator"] = [
        "solar",
        "solar-hsat",
        "gas",
    ]
    return solve_config


SWITCHES = {
    "defaults": (
        {},
        {
            "land_use_constraint_solar_2040",
            "Link_charger_ratio_1",
            "Link_bidirectional_sync_0",
            "Carrier_growth_limit_solar",
            "solar_potential",
        },
    ),
    "BAU": ({"BAU": True}, {"bau_mincaps"}),
    "SAFE": ({"SAFE": True}, {"safe_mintotalcap"}),
    "EQ": ({"EQ": "EQ0.7"}, {"equity_min"}),
    "reserve": (
        {"operational_reserve": True},
        {"reserve_margin", "Generator_p_reserve_upper_fixed_idle"},
    ),
}


@pytest.mark.parametrize(
    ("switches", "expected"), SWITCHES.values(), ids=SWITCHES.keys()
)
def test_fragments_compose_onto_base_spec(
    network: pypsa.Network,
    solve_config: dict,
    base: dict,
    switches: dict,
    expected: set[str],
) -> None:
    import mathspec
    from datarecord.tools.pypsa import MATHSPEC, PyPSA

    solve_config["electricity"]["operational_reserve"]["activate"] = switches.pop(
        "operational_reserve", False
    )
    solve_config["solving"]["constraints"].update(switches)
    fragments = [
        mf.global_constraints(network, base),
        mf.extra_functionality(network, solve_config, None, base),
    ]
    for fragment in fragments:
        mathspec.to_spec(fragment)
    spec = mathspec.merge([MATHSPEC, *fragments])
    PyPSA.check_spec(spec, PyPSA.to_datarecord(network))
    assert expected <= set(spec.constraints)


@pytest.mark.parametrize(
    ("name", "rhs"),
    [
        ("CO2Limit_upper_2040", "<= 500.0"),
        ("CO2Budget", "<= 10000.0"),
        ("GasLimit", "<= 100000.0"),
        ("lv_limit", "<= 10000.0"),
        ("lc_limit", "<= 1000000.0"),
        ("hydro_limit", "<= 980.0"),
        ("solar_cap", "<= 1000.0"),
    ],
)
def test_global_constraint_constants(
    network: pypsa.Network, base: dict, name: str, rhs: str
) -> None:
    constraint = mf.global_constraints(network, base)["constraints"][
        f"GlobalConstraint_{name}"
    ]
    assert constraint["expression"].endswith(rhs)


def test_co2_atmosphere_reads_the_rise_over_its_period(
    network: pypsa.Network, base: dict
) -> None:
    level = mf.global_constraints(network, base)["expressions"][
        "GlobalConstraint_CO2Limit_upper_2040_level"
    ]
    assert {key: case["expression"] for key, case in level["cases"].items()} == {
        "final": "1",
        "previous": "-1",
    }
    assert "snapshot_period[period] == 2030" in level["cases"]["previous"]["when"]


def test_battery_pairs_inline_the_discharger_efficiency(
    network: pypsa.Network, solve_config: dict, base: dict
) -> None:
    fragment = mf.extra_functionality(network, solve_config, None, base)
    cases = fragment["expressions"]["Link_charger_ratio_0_coefficient"]["cases"]
    assert [(case["when"], case["expression"]) for case in cases.values()] == [
        ("entity == 'b0 battery charger-2030'", "1.0"),
        ("entity == 'b0 battery discharger-2030'", "-0.9"),
    ]


def test_unknown_global_constraint_type_raises(
    network: pypsa.Network, base: dict
) -> None:
    network.add("GlobalConstraint", "odd", type="odd_type", constant=1.0)
    with pytest.raises(mf.UnsupportedConstraintError, match="odd_type"):
        mf.global_constraints(network, base)


@pytest.mark.parametrize(
    "switches",
    [{"CCL": True}, {"EQ": "EQ0.7c"}],
    ids=["CCL", "EQ per country"],
)
def test_untranslatable_constraint_raises(
    network: pypsa.Network, solve_config: dict, base: dict, switches: dict
) -> None:
    solve_config["solving"]["constraints"].update(switches)
    with pytest.raises(mf.UnsupportedConstraintError):
        mf.extra_functionality(network, solve_config, None, base)


@pytest.mark.parametrize(
    ("body", "raises"),
    [('    """Nothing."""\n    pass\n', False), ("    n.add('Bus', 'x')\n", True)],
    ids=["empty", "custom"],
)
def test_custom_extra_functionality(
    network: pypsa.Network,
    solve_config: dict,
    base: dict,
    tmp_path: Path,
    body: str,
    raises: bool,
) -> None:
    path = tmp_path / "custom.py"
    path.write_text(f"def custom(n, snapshots, snakemake):\n{body}")
    if raises:
        with pytest.raises(mf.UnsupportedConstraintError, match="custom.py"):
            mf.extra_functionality(network, solve_config, str(path), base)
    else:
        mf.extra_functionality(network, solve_config, str(path), base)


@pytest.mark.parametrize(("country", "raises"), [("DE", False), ("FR", True)])
def test_component_country_must_match_its_bus(
    network: pypsa.Network, country: str, raises: bool
) -> None:
    from scripts.build_datarecord import drop_bus_country_duplicates

    network.loads["country"] = country
    if raises:
        with pytest.raises(ValueError, match="b0 load"):
            drop_bus_country_duplicates(network)
    else:
        drop_bus_country_duplicates(network)
        assert "country" not in network.loads


def test_export_datarecord(
    network: pypsa.Network, solve_config: dict, tmp_path: Path
) -> None:
    from scripts.build_datarecord import (
        export_datarecord,
        prefix_shared_bus_columns,
        rename_shared_one_port_names,
    )

    network.add("Load", "H2 pipeline-2030", bus="b0 H2", p_set=1.0)
    network.buses["tags"] = "bus"
    network.lines["tags"] = "line"
    rename_shared_one_port_names(network)
    prefix_shared_bus_columns(network)
    export_datarecord(network, solve_config, None, tmp_path / "composed")
    assert "H2 pipeline-2030 (Load)" in network.loads.index
    assert "bus_tags" in network.buses
    assert sorted(p.name for p in (tmp_path / "composed" / "math").iterdir()) == [
        "extra_functionality.yaml",
        "global_constraints.yaml",
        "pypsa.mathspec.yaml",
    ]
