# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for selected helpers in scripts/prepare_sector_network.py."""

from types import SimpleNamespace

import pandas as pd
import pypsa
import pytest

from scripts.prepare_sector_network import (
    connect_existing_chp_units,
    convert_existing_units_to_links,
    existing_units,
)


@pytest.fixture
def unit_network() -> pypsa.Network:
    n = pypsa.Network()
    n.add("Bus", ["b1", "b2"], carrier="AC")
    n.add(
        "Generator",
        ["b1 CCGT 1", "b2 nuclear 2", "b2 biomass 3", "b1 OCGT", "b1 solar"],
        bus=["b1", "b2", "b2", "b1", "b1"],
        carrier=["CCGT", "nuclear", "biomass", "OCGT", "solar"],
        p_nom=[100.0, 500.0, 30.0, 0.0, 10.0],
        p_nom_extendable=[False, False, False, True, False],
        p_max_pu=[1.0, 0.8, 1.0, 1.0, 1.0],
        efficiency=[0.5, 0.33, 0.3, 0.4, 1.0],
        build_year=[2000, 1980, 2010, 0, 2015],
        lifetime=[40.0, 50.0, 30.0, 30.0, 25.0],
    )
    n.generators["geometry"] = "POINT (0 0)"
    return n


@pytest.fixture
def unit_costs() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "VOM": [4.0, 3.0, 2.0],
            "capital_cost": [1e4, 1e5, 4e4],
            "CO2 intensity": [0.2, 0.0, 0.0],
            "efficiency": [0.4, 0.33, 0.3],
            "efficiency-heat": [0.0, 0.0, 0.6],
        },
        index=["CCGT", "nuclear", "central solid biomass CHP"],
    ).reindex(["CCGT", "nuclear", "central solid biomass CHP", "gas", "uranium"])


def test_convert_existing_units_to_links(unit_network, unit_costs):
    n = unit_network
    spatial = SimpleNamespace(
        gas=SimpleNamespace(nodes=["b1 gas", "b2 gas"], locations=["b1", "b2"]),
        uranium=SimpleNamespace(nodes=["EU uranium"], locations=["EU"]),
        biomass=SimpleNamespace(nodes=["EU solid biomass"], locations=["EU"]),
    )
    unit_costs.loc["gas", "CO2 intensity"] = 0.2
    unit_costs.loc["uranium", "CO2 intensity"] = 0.0

    fuels = convert_existing_units_to_links(n, existing_units(n), unit_costs, spatial)

    assert fuels == {"gas", "uranium", "biomass"}
    assert n.generators.index.tolist() == ["b1 OCGT", "b1 solar"]
    links = n.links
    assert links.loc[
        ["b1 CCGT 1", "b2 nuclear 2", "b2 biomass 3"], "bus0"
    ].tolist() == [
        "b1 gas",
        "EU uranium",
        "EU solid biomass",
    ]
    assert links.loc["b1 CCGT 1", ["p_nom", "efficiency", "efficiency2"]].tolist() == [
        200.0,
        0.5,
        0.2,
    ]
    assert links.at["b1 CCGT 1", "marginal_cost"] == 2.0
    assert links.at["b2 nuclear 2", "p_max_pu"] == 0.8
    assert links.at["b2 biomass 3", "carrier"] == "urban central solid biomass CHP"
    assert links.at["b2 biomass 3", "p_nom"] == pytest.approx(100.0)
    assert (links.build_year.tolist(), links.asset_status.unique().tolist()) == (
        [2000, 1980, 2010],
        ["existing"],
    )
    assert not links.p_nom_extendable.any()


def test_connect_existing_chp_units():
    n = pypsa.Network()
    n.add("Bus", ["b1", "b2", "b1 urban central heat"])
    n.add(
        "Link",
        ["b1 biomass 1", "b2 biomass 2"],
        bus1=["b1", "b2"],
        carrier="urban central solid biomass CHP",
    )

    connect_existing_chp_units(n)

    assert n.links.bus2.tolist() == ["b1 urban central heat", ""]
