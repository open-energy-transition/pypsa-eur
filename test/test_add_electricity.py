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
    attach_storageunits,
    attach_stores,
    conventional_candidates,
    estimate_efficiency,
    group_small_units,
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


def test_assign_unit_dates():
    units = pd.DataFrame(
        {
            "carrier": ["CCGT", "CCGT", "CCGT", "coal"],
            "datein": [2000.0, np.nan, 1990.0, np.nan],
            "dateout": [np.nan, 2045.0, 2025.0, 2040.0],
        }
    )
    lifetime = pd.Series(40.0, index=units.index)

    result = assign_unit_dates(units, lifetime, first_horizon=2030)

    assert result.index.tolist() == [0, 1]
    assert result.build_year.tolist() == [2000, 1995]
    assert result.lifetime.tolist() == [40.0, 50.0]


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
    ("keep_assets", "expected"),
    [
        (True, {"DE0 1 CCGT 0": 100.0, "DE0 1 CCGT 1": 300.0}),
        (False, {"DE0 1 CCGT": 400.0}),
    ],
)
def test_load_and_aggregate_powerplants_keeps_units_at_cluster_bus(
    tmp_path, keep_assets, expected
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
        ppl_fn, costs, keep_assets=keep_assets, first_horizon=2030
    )

    assert ppl.p_nom.to_dict() == expected, "units keep their cluster bus"
    assert (ppl.bus == "DE0 1").all()
