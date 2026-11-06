# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for selected helpers in scripts/add_electricity.py."""

import numpy as np
import pandas as pd
import pypsa
import pytest
import xarray as xr

from scripts.add_electricity import (
    assign_unit_dates,
    attach_load,
    attach_renewable_estimates,
    attach_storageunits,
    attach_stores,
    cohort_build_years,
    conventional_candidates,
    estimate_efficiency,
    group_small_units,
    last_active_horizon,
    load_and_aggregate_powerplants,
    reduce_candidate_potentials,
)
from scripts.lib.validation.config.conventional import _EstimateEfficienciesConfig


def test_attach_load(tmp_path):
    """Clustered demand is attached per bus with scaling applied."""

    times = pd.date_range("2000-01-01", periods=2, freq="h")
    buses = ["zone_1", "zone_2"]
    values = np.array([[1.0, 2.0], [3.0, 4.0]])

    data = xr.DataArray(
        values,
        coords={"time": times, "bus": buses},
        dims=["time", "bus"],
        name="electricity demand (MW)",
    )
    load_path = tmp_path / "electricity_demand.nc"
    data.to_netcdf(load_path)

    n = pypsa.Network()
    n.set_snapshots(times)
    n.add("Bus", buses)

    attach_load(n, load_path.as_posix(), scaling=2.0)

    assert sorted(n.loads.index) == buses
    assert sorted(n.loads_t.p_set.columns) == buses
    np.testing.assert_allclose(n.loads_t.p_set[buses].values, 2.0 * values)


def test_attach_storageunits_energy_basis():
    """A dispatched-basis `max_hours` sizes the store to sustain it, at no extra cost."""
    costs = pd.DataFrame(
        {
            "capital_cost": {"iron-air": 157959.0, "battery": 60630.0},
            "marginal_cost": {"iron-air": 0.0, "battery": 0.0},
            "lifetime": {"iron-air": 17.5, "battery": 17.5},
            "efficiency": {
                "iron-air battery charge": 0.74,
                "iron-air battery discharge": 0.63,
                "battery inverter": 0.96,
            },
        }
    )
    max_hours = {"iron-air": 100, "battery": 6}

    n = pypsa.Network()
    n.add("Bus", ["bus_1", "bus_2"])
    attach_storageunits(n, costs, n.buses.index, ["iron-air", "battery"], max_hours)

    su = n.storage_units.set_index("carrier")

    iron_air = su.loc["iron-air"]
    np.testing.assert_allclose(iron_air.max_hours * iron_air.efficiency_dispatch, 100)
    np.testing.assert_allclose(iron_air.capital_cost, 157959.0)

    # Stored-basis carriers must be left alone, despite efficiency_dispatch < 1.
    np.testing.assert_allclose(su.loc["battery"].max_hours, 6)


def test_attach_stores_energy_basis():
    """A dispatched-basis store cost is converted to cost per MWh stored."""
    costs = pd.DataFrame(
        {
            "capital_cost": {
                "iron-air battery": 1000.0,
                "iron-air battery charge": 0.0,
                "iron-air battery discharge": 0.0,
            },
            "marginal_cost": 0.0,
            "lifetime": 17.5,
            "efficiency": {
                "iron-air battery": 1.0,
                "iron-air battery charge": 0.74,
                "iron-air battery discharge": 0.63,
            },
        }
    )

    n = pypsa.Network()
    n.add("Bus", ["bus_1"])
    attach_stores(n, costs, n.buses.index, ["iron-air"])

    np.testing.assert_allclose(n.stores.capital_cost, 630.0)


def test_estimate_efficiency():
    """Efficiency rises with build year (retrofit first), is clipped and degrades with age."""
    ppl = pd.DataFrame(
        {
            "carrier": ["CCGT", "coal", "coal", "biomass", "CCGT"],
            "datein": [2000, 1950, 1970, 2000, np.nan],
            "dateretrofit": [np.nan, np.nan, 2010, np.nan, np.nan],
        }
    )
    config = _EstimateEfficienciesConfig().model_dump()

    eta = estimate_efficiency(ppl, config)

    expected = [
        0.48 * (1 - 0.015),
        0.28 * (1 - 0.065),
        0.43 * (1 - 0.005),
        np.nan,
        np.nan,
    ]
    np.testing.assert_allclose(eta, expected)


@pytest.mark.parametrize(
    ("earliest_retirement", "expected_lifetime"),
    [(None, [40.0, 50.0, 45.0]), (2045, [45.0, 50.0, 45.0])],
)
def test_assign_unit_dates(earliest_retirement, expected_lifetime):
    units = pd.DataFrame(
        {
            "carrier": ["CCGT", "CCGT", "CCGT", "coal"],
            "datein": [2000.0, np.nan, 1990.0, np.nan],
            "dateout": [np.nan, 2045.0, 2025.0, 2040.0],
        }
    )
    lifetime = pd.Series(40.0, index=units.index)

    result = assign_unit_dates(units, lifetime, 2030, earliest_retirement)

    assert result.index.tolist() == [0, 1, 3]
    assert result.build_year.tolist() == [2000, 1995, 1995]
    assert result.lifetime.tolist() == expected_lifetime


def test_last_active_horizon():
    units = pd.DataFrame(
        {"build_year": [1990, 2005, 2016, 2035], "lifetime": [30, 30, 25, 30]}
    )

    result = last_active_horizon(units, [2030, 2040, 2050])

    np.testing.assert_array_equal(result, [np.nan, 2030, 2040, 2050])


@pytest.mark.parametrize(
    ("min_capacity", "expected"),
    [
        (0.0, {"a": 0.5, "b": 0.7, "c": 0.4, "d": 3.0}),
        (1.0, {"small-2010-25": 1.2, "small-2015-25": 0.4, "d": 3.0}),
    ],
)
def test_group_small_units(min_capacity, expected):
    """Small units of one candidate, build year and lifetime merge; others stay."""
    units = pd.DataFrame(
        {
            "capacity": [0.5, 0.7, 0.4, 3.0],
            "build_year": [2010, 2010, 2015, 2010],
            "lifetime": [25.0, 25.0, 25.0, 25.0],
            "lon": [4.0, 4.1, 4.2, 4.3],
            "lat": [50.0, 50.1, 50.2, 50.3],
        },
        index=["a", "b", "c", "d"],
    )
    candidates = pd.Series("BE0 0 solar", index=units.index)

    result = group_small_units(units, candidates, min_capacity)

    assert result.capacity.to_dict() == pytest.approx(expected)
    assert (result.candidate == "BE0 0 solar").all()
    if min_capacity:
        assert result.geometry["small-2010-25"].startswith("MULTIPOINT")


def test_conventional_candidates():
    n = pypsa.Network()
    n.add("Bus", ["ac", "dc"], carrier=["AC", "DC"])
    costs = pd.DataFrame(
        {
            "efficiency": [0.4],
            "VOM": [5.0],
            "fuel": [20.0],
            "capital_cost": [1e4],
            "lifetime": [25.0],
        },
        index=["OCGT"],
    )

    candidates = conventional_candidates(n, costs, ["OCGT"])

    assert candidates.index.tolist() == ["ac OCGT"]
    assert candidates.loc[
        "ac OCGT", ["p_nom", "build_year", "marginal_cost"]
    ].tolist() == [
        0.0,
        0,
        55.0,
    ]


def test_reduce_candidate_potentials():
    n = pypsa.Network()
    n.add("Bus", "b")
    n.add(
        "Generator",
        ["b 0 onwind", "b 1 onwind", "b 0 onwind 7", "b 0 onwind estimate"],
        bus="b",
        carrier="onwind",
        p_nom=[0.0, 0.0, 30.0, 80.0],
        p_nom_max=[100.0, 100.0, np.inf, np.inf],
        asset_status=["candidate", "candidate", "existing", "estimate"],
        resource_class=[0, 1, 0, 0],
    )

    reduce_candidate_potentials(n, ["onwind"])

    assert n.generators.p_nom_max[["b 0 onwind", "b 1 onwind"]].tolist() == [0.0, 100.0]


@pytest.mark.parametrize(
    ("keep_assets", "horizons", "expected"),
    [
        (True, [2030, 2040], {"DE0 1 CCGT 0": 100.0, "DE0 1 CCGT 1": 300.0}),
        (False, [2030, 2040], {"DE0 1 CCGT 2030": 100.0, "DE0 1 CCGT 2040": 300.0}),
        (False, None, {"DE0 1 CCGT": 400.0}),
    ],
)
def test_load_and_aggregate_powerplants_keeps_units_at_cluster_bus(
    tmp_path, keep_assets, horizons, expected
):
    ppl_fn = tmp_path / "powerplants.csv"
    pd.DataFrame(
        {
            "Name": ["a", "b"],
            "Fueltype": ["Natural Gas", "Natural Gas"],
            "Technology": ["CCGT", "CCGT"],
            "Set": "PP",
            "Country": "DE",
            "Capacity": [100.0, 300.0],
            "Efficiency": [0.5, 0.6],
            "DateIn": [2000.0, 2010.0],
            "DateOut": [2040.0, 2050.0],
            "lat": [50.0, 51.0],
            "lon": [10.0, 11.0],
            "bus": "DE0 1",
        }
    ).to_csv(ppl_fn)
    costs = pd.DataFrame(
        {
            "VOM": [4.0],
            "FOM": [3.0],
            "efficiency": [0.55],
            "capital_cost": [1e4],
            "marginal_cost": [50.0],
            "fuel": [20.0],
            "lifetime": [30.0],
        },
        index=["CCGT"],
    )

    ppl = load_and_aggregate_powerplants(
        ppl_fn, costs, keep_assets=keep_assets, horizons=horizons
    )

    assert ppl.p_nom.to_dict() == expected, "units keep their cluster bus"
    assert (ppl.bus == "DE0 1").all()


@pytest.mark.parametrize(
    ("horizons", "expected"),
    [
        ([2030, 2040, 2050], [2000, 2010, 2020, 2024, np.nan]),
        ([2030], [2000, 2024, 2024, 2024, np.nan]),
    ],
)
def test_cohort_build_years(horizons, expected):
    build_year = pd.Series([1995, 2005, 2015, 2022, 2025])
    result = cohort_build_years(build_year, horizons, 30, 2024)
    pd.testing.assert_series_equal(result, pd.Series(expected), check_names=False)


def test_attach_renewable_estimates(monkeypatch):
    """Missing capacity per cohort becomes one estimate per cohort, retired cohorts are dropped."""
    additions = pd.DataFrame(
        [[20.0, 30.0, 40.0, 50.0]],
        index=pd.MultiIndex.from_tuples(
            [("Onshore", "BE")], names=["Technology", "Country"]
        ),
        columns=pd.Index([1995, 2005, 2015, 2022], name="Year"),
    )
    monkeypatch.setattr(
        "scripts.add_electricity.irena_additions", lambda *args: additions
    )
    n = pypsa.Network()
    n.set_snapshots(range(2))
    n.add("Bus", "b", country="BE")
    n.add(
        "Generator",
        "b 0 onwind",
        bus="b",
        carrier="onwind",
        p_nom_extendable=True,
        p_nom_max=100.0,
        p_max_pu=pd.Series([0.4, 0.6]),
        lifetime=30.0,
        resource_class=0,
        asset_status="candidate",
    )
    n.add(
        "Generator",
        "b 0 onwind unit",
        bus="b",
        carrier="onwind",
        p_nom=10.0,
        build_year=2015,
        lifetime=30.0,
        asset_status="existing",
    )

    attach_renewable_estimates(
        n, 2024, {"Onshore": "onwind"}, ["BE"], [2030, 2040, 2050]
    )

    estimates = n.generators.query("asset_status == 'estimate'")
    assert estimates[["p_nom", "build_year"]].to_dict("index") == {
        "b 0 onwind estimate-2010": {"p_nom": 30.0, "build_year": 2010},
        "b 0 onwind estimate-2020": {"p_nom": 30.0, "build_year": 2020},
        "b 0 onwind estimate-2024": {"p_nom": 50.0, "build_year": 2024},
    }
