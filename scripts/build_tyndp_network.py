# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import geopandas as gpd
import pandas as pd
from geopandas import GeoDataFrame
from shapely.geometry import LineString

from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

GEO_CRS = "EPSG:4326"
DISTANCE_CRS = "EPSG:3035"
BUSES_COLUMNS = [
    "station_id",
    "voltage",
    "dc",
    "symbol",
    "under_construction",
    "tags",
    "x",
    "y",
    "country",
    "geometry",
]
LINES_COLUMNS = [
    "bus0",
    "bus1",
    "voltage",
    "circuits",
    "length",
    "underground",
    "under_construction",
    "tags",
    "geometry",
]
LINKS_COLUMNS = [
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "length",
    "underground",
    "under_construction",
    "tags",
    "geometry",
]
TRANSFORMERS_COLUMNS = [
    "bus0",
    "bus1",
    "voltage_bus0",
    "voltage_bus1",
    "s_nom",
    "station_id",
    "geometry",
]
CONVERTERS_COLUMNS = [
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "geometry",
]


def format_bz_names(s: str):
    s = (s
         .replace("DK_1", "DKW1")
         .replace("DK_2", "DKE1")
         .replace("GB", "UK00")
         .replace("IT_NORD", "ITN1")
         .replace("IT_SUD", "ITS1")
         .replace("NO_1", "NOS0")
         .replace("NO_3", "NOM1")
         .replace("NO_4", "NON1")
         .replace('SE_', 'SE0')
         .replace('_', '')
         .ljust(4, "0")
         )[:4]
    return s


def build_shapes(
        bz_fn,
        geo_crs: str = GEO_CRS,
        distance_crs: str = DISTANCE_CRS
):
    """
    Process bidding zones from the shape file and calculate representative point. Deduce the country shapes and their representative point.

    Parameters
    ----------
        - bz_fn (str | Path): Path to bidding zone shape file.
        - geo_crs (CRS, optional): Coordinate reference system for geographic calculations. Defaults to GEO_CRS.
        - distance_crs (CRS, optional): Coordinate reference system to use for distance calculations. Defaults to DISTANCE_CRS.

    Returns
    -------
        - bidding_shapes: A GeoDataFrame including bidding zone geometry, representative point and id.
        - country_shapes: A GeoDataFrame including country geometry and representative point.
    """
    # Bidding zone shapes
    bidding_shapes = (
        gpd.read_file(bz_fn)
        .assign(
            bz_id=lambda df: df["id"].apply(format_bz_names),
            node=lambda df: df.geometry.to_crs(distance_crs).representative_point().to_crs(geo_crs),
            x=lambda df: df["node"].x,
            y=lambda df: df["node"].y
        )
        .set_index("bz_id")
    )

    # Country shapes
    country_shapes = (
        bidding_shapes
        .dissolve(by="country")[["geometry"]]
        .assign(
            node=lambda df: df.geometry.to_crs(distance_crs).representative_point().to_crs(geo_crs),
            x=lambda df: df["node"].x,
            y=lambda df: df["node"].y
        )
    )

    return bidding_shapes, country_shapes


def build_buses(
        buses_fn,
        bidding_shapes: GeoDataFrame,
        country_shapes: GeoDataFrame,
        geo_crs: str = GEO_CRS,
        distance_crs: str = DISTANCE_CRS
):
    """
    Extend the node list for both electricity and hydrogen with attributes, incl. country and coordinates.

    Parameters
    ----------
        - buses_fn (str | Path): Path to bidding zone shape file.
        - bidding_shapes (GeoDataFrame): A GeoDataFrame including bidding zone geometry, representative point and id.
        - country_shapes (GeoDataFrame): A GeoDataFrame including country geometry and representative point.
        - geo_crs (CRS, optional): Coordinate reference system for geographic calculations. Defaults to GEO_CRS.
        - distance_crs (CRS, optional): Coordinate reference system to use for distance calculations. Defaults to DISTANCE_CRS.


    Returns
    -------
        - buses: A GeoDataFrame of electrical buses including country and coordinates.
        - buses_h2: A GeoDataFrame of hydrogen buses including country and coordinates.
    """
    buses = (
        pd.read_excel(buses_fn)
        .merge(bidding_shapes[["country", "node", "x", "y"]], how="left", left_on="NODE", right_index=True)
        .rename({"NODE": "bus_id", "node": "geometry"}, axis=1)
        .assign(
            station_id=lambda df: df["bus_id"],
            voltage=380,  # TODO Improve assumption
            dc=None,
            symbol="Substation",
            under_construction="f",
            tags=lambda df: df["bus_id"],
        )
        .set_index("bus_id")
        [BUSES_COLUMNS]
    )
    buses = gpd.GeoDataFrame(buses, geometry="geometry", crs=geo_crs)

    buses_h2 = (
        country_shapes[["node", "x", "y"]]
        .reset_index()
        .rename({"node": "geometry"}, axis=1)
        .assign(
            bus_id=lambda df: df[["country"]] + " H2",
            station_id=lambda df: df["bus_id"],
            voltage=None,
            dc="f",
            symbol="Substation",
            under_construction="f",
            tags=lambda df: df["bus_id"],
        )
        .set_index("bus_id")
        [BUSES_COLUMNS]
    )
    buses_h2 = gpd.GeoDataFrame(buses_h2, geometry="geometry", crs=geo_crs)

    return buses, buses_h2


def build_links(
        grid_fn,
        buses: GeoDataFrame,
        geo_crs: str = GEO_CRS,
        distance_crs: str = DISTANCE_CRS
):
    """
    Process reference grid information to produce link data. p_nom are NTC values.

    Parameters
    ----------
        - grid_fn (str | Path): Path to bidding zone shape file.
        - geo_crs (CRS, optional): Coordinate reference system for geographic calculations. Defaults to GEO_CRS.
        - distance_crs (CRS, optional): Coordinate reference system to use for distance calculations. Defaults to DISTANCE_CRS.

    Returns
    -------
        - links: A GeoDataFrame including NTC from the reference grid.

    """
    links = pd.read_excel(grid_fn)
    links[["bus0", "bus1"]] = links.Border.str.split("-", expand=True)

    # Create forward and reverse direction dataframes
    forward_links = (
        links[["bus0", "bus1", "Summary Direction 1"]]
        .rename(columns={"Summary Direction 1": "p_nom"})
    )

    reverse_links = (
        links[["bus1", "bus0", "Summary Direction 2"]]
        .rename(columns={
            "bus1": "bus0",
            "bus0": "bus1",
            "Summary Direction 2": "p_nom"
        })
    )

    # Combine into unidirectional links
    links = pd.concat([forward_links, reverse_links])

    # Add missing attributes
    links = (
        links
        .merge(buses["geometry"], how="left", left_on="bus0", right_index=True)
        .merge(buses["geometry"], how="left", left_on="bus1", right_index=True, suffixes=("0", "1"))
        .dropna()  # TODO Remove this when all nodes are known
    )
    links["geometry"] = gpd.GeoSeries([LineString([p0, p1]) for p0, p1 in zip(links["geometry0"], links["geometry1"])])
    links = gpd.GeoDataFrame(links, geometry="geometry", crs=geo_crs)

    links = (
        links
        .assign(
            link_id=lambda df: df["bus0"] + "-" + df["bus1"] + "-DC",
            voltage=380,  # TODO Improve assumption
            length=lambda df: df.geometry.to_crs(distance_crs).length,
            underground="t",
            under_construction="f",
            tags=lambda df: df["bus0"] + " > " + df["bus1"],
        )
        .set_index("link_id")
        [LINKS_COLUMNS]
    )

    return links


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_tyndp_network")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Build node coordinates
    bidding_shapes, country_shapes = build_shapes(snakemake.input.bidding_shapes)
    buses, buses_h2 = build_buses(snakemake.input.buses, bidding_shapes, country_shapes)

    # Build links
    links = build_links(snakemake.input.reference_grid, buses)

    lines = gpd.GeoDataFrame(columns=LINES_COLUMNS, geometry="geometry").set_index(pd.Index([], name="line_id"))
    converters = gpd.GeoDataFrame(columns=CONVERTERS_COLUMNS, geometry="geometry").set_index(pd.Index([], name="converter_id"))
    transformers = gpd.GeoDataFrame(columns=TRANSFORMERS_COLUMNS, geometry="geometry").set_index(pd.Index([], name="transformer_id"))

    # Export to csv for base_network
    buses.to_csv(snakemake.output["substations"], quotechar="'")
    buses_h2.to_csv(snakemake.output["substations_h2"], quotechar="'")
    lines.to_csv(snakemake.output["lines"], quotechar="'")
    links.to_csv(snakemake.output["links"], quotechar="'")
    converters.to_csv(snakemake.output["converters"], quotechar="'")
    transformers.to_csv(snakemake.output["transformers"], quotechar="'")

    # Export to GeoJSON for quick validations
    buses.to_file(snakemake.output["substations_geojson"])
    buses_h2.to_file(snakemake.output["substations_h2_geojson"])
    lines.to_file(snakemake.output["lines_geojson"])
    links.to_file(snakemake.output["links_geojson"])
    converters.to_file(snakemake.output["converters_geojson"])
    transformers.to_file(snakemake.output["transformers_geojson"])
