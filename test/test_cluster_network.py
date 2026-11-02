# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for helper utilities in scripts/cluster_network.py."""

import geopandas as gpd
import pandas as pd
import pypsa
import pytest
from shapely.geometry import Polygon

from scripts.cluster_network import (
    busmap_from_shapes,
    cluster_reference_busmaps,
    sanitize_busmap,
)


def test_sanitize_busmap_strips_values():
    """Whitespace in busmap values is removed while indices stay untouched."""

    raw = pd.Series([" zone_a ", "zone_b  "], index=[" bus1 ", "bus2  "], name="busmap")

    cleaned = sanitize_busmap(raw)

    assert cleaned.index.tolist() == [" bus1 ", "bus2  "]
    assert cleaned.tolist() == ["zone_a", "zone_b"]
    assert cleaned.name == "busmap"


def test_busmap_from_shapes_returns_trimmed_labels():
    """Busmap inferred from shapes produces sanitized labels."""

    n = pypsa.Network()
    n.add("Bus", "bus1", x=0.0, y=0.0, carrier="AC", country="AA")
    n.add("Bus", "bus2", x=2.0, y=0.0, carrier="AC", country="AA")

    shapes = gpd.GeoDataFrame(
        {
            "name": [" region_1 ", "region_2  "],
            "geometry": [
                Polygon([(-1, -1), (1, -1), (1, 1), (-1, 1)]),
                Polygon([(1, -1), (3, -1), (3, 1), (1, 1)]),
            ],
        },
        crs="EPSG:4326",
    )

    busmap = busmap_from_shapes(n, shapes)

    assert busmap.index.tolist() == ["bus1", "bus2"]
    assert busmap.tolist() == ["region_1", "region_2"]


@pytest.mark.parametrize(
    ("busmap", "expected"),
    [
        (["b1", "b2", "b3", "b4"], ["r1", "r1", "r2", "r2"]),
        (["c1", "c1", "c1", "c2"], ["r2", "r2"]),
    ],
    ids=["unclustered", "clustered"],
)
def test_cluster_reference_busmaps_takes_reference_bus_with_most_load(busmap, expected):
    buses = pd.Index(["b1", "b2", "b3", "b4"])
    busmaps = pd.DataFrame({"busmap_2": ["r1", "r1", "r2", "r2"]}, index=buses)
    load = pd.Series([1.0, 2.0, 5.0, 0.0], index=buses)

    result = cluster_reference_busmaps(busmaps, pd.Series(busmap, index=buses), load)

    assert result.busmap_2.tolist() == expected
    assert result.index.tolist() == sorted(set(busmap))
