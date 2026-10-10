# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Add electricity components to the clustered network: load, existing conventional and hydro power plants, renewable generators and extendable storage.

Load time series are attached from the clustered electricity demand. Power
plants from [build_powerplants][] are aggregated per bus and carrier
(optionally by efficiency class) and attached as conventional generators with
fuel costs, efficiencies and emissions from the cost assumptions; hydro plants
become run-of-river generators, reservoir and pumped hydro storage units with
inflow from [build_hydro_profile][]. Wind and solar generators are added with
the profiles from [build_renewable_profiles][] and zero existing capacity,
unless existing capacities are estimated from IRENASTAT or powerplantmatching.
Transmission costs of lines and links are set from the cost assumptions.
Carriers in `electricity.extendable_carriers` are added with zero initial
capacity. Storage is either a `StorageUnit` whose energy-to-power ratio is
fixed by `electricity.max_hours` (one investment variable) or a `Store` on an
extra bus with separate charging and discharging links (three investment
variables).

!!! note "Called from compose_network"
    This module has no Snakemake rule of its own. Its `main` function is
    called by [compose_network][].
"""

import logging
from collections.abc import Iterable
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
import xarray as xr
from pypsa.clustering.spatial import DEFAULT_ONE_PORT_STRATEGIES, normed_or_uniform

from scripts._helpers import (
    PYPSA_V1,
    rename_techs,
    update_p_nom_max,
)

if PYPSA_V1:
    pypsa.options.params.add.return_names = True

STORE_LOOKUP = {
    "battery": {
        "store": "battery storage",
        "bicharger": "battery inverter",
        "roundtrip_correction": 0.5,
    },
    "home battery": {
        "store": "home battery storage",
        "bicharger": "home battery inverter",
        "roundtrip_correction": 0.5,
    },
    "li-ion": {
        "store": "battery storage",
        "bicharger": "battery inverter",
        "roundtrip_correction": 0.5,
    },
    "lfp": {
        "store": "Lithium-Ion-LFP-store",
        "bicharger": "Lithium-Ion-LFP-bicharger",
    },
    "vanadium": {
        "store": "Vanadium-Redox-Flow-store",
        "bicharger": "Vanadium-Redox-Flow-bicharger",
    },
    "lair": {
        "store": "Liquid-Air-store",
        "charger": "Liquid-Air-charger",
        "discharger": "Liquid-Air-discharger",
    },
    "pair": {
        "store": "Compressed-Air-Adiabatic-store",
        "bicharger": "Compressed-Air-Adiabatic-bicharger",
    },
    "iron-air": {
        "store": "iron-air battery",
        "charger": "iron-air battery charge",
        "discharger": "iron-air battery discharge",
        # duration and EUR/kWh are quoted on dispatched, not stored, energy
        "energy_basis": "dispatched",
    },
    "H2": {
        "store": "hydrogen storage underground",
        "charger": "electrolysis",
        "discharger": "fuel cell",
    },
    "H2 tank": {
        "store": "hydrogen storage tank type 1 including compressor",
        "charger": "electrolysis",
        "discharger": "fuel cell",
    },
}

idx = pd.IndexSlice

logger = logging.getLogger(__name__)


def normed(s: pd.Series) -> pd.Series:
    """
    Normalize a pandas Series by dividing each element by the sum of all elements.

    Parameters
    ----------
    s : pd.Series
        Input series to normalize

    Returns
    -------
    pd.Series
        Normalized series where all elements sum to 1
    """
    return s / s.sum()


def flatten(t: Iterable[Any]) -> str:
    return " ".join(map(str, t))


def calculate_annuity(n: float, r: float | pd.Series) -> float | pd.Series:
    """
    Calculate the annuity factor for an asset with lifetime n years and discount rate r.

    The annuity factor is used to calculate the annual payment required to pay off a loan
    over n years at interest rate r. For example, annuity(20, 0.05) * 20 = 1.6.

    Parameters
    ----------
    n : float
        Lifetime of the asset in years
    r : float | pd.Series
        Discount rate (interest rate). Can be a single float or a pandas Series of rates.

    Returns
    -------
    float | pd.Series
        Annuity factor. Returns a float if r is float, or pd.Series if r is pd.Series.

    Examples
    --------
    >>> calculate_annuity(20, 0.05)
    0.08024258718774728
    """
    if isinstance(r, pd.Series):
        return pd.Series(1 / n, index=r.index).where(
            r == 0, r / (1.0 - 1.0 / (1.0 + r) ** n)
        )
    elif r > 0:
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    else:
        return 1 / n


def add_missing_carriers(n, carriers):
    """
    Function to add missing carriers to the network without raising errors.
    """
    missing_carriers = set(carriers) - set(n.carriers.index)
    if len(missing_carriers) > 0:
        n.add("Carrier", missing_carriers)


def sanitize_carriers(n, config):
    """
    Sanitize the carrier information in a PyPSA Network object.

    The function ensures that all unique carrier names are present in the network's
    carriers attribute, and adds nice names and colors for each carrier according
    to the provided configuration dictionary.

    Parameters
    ----------
    n : pypsa.Network
        A PyPSA Network object that represents an electrical power system.
    config : dict
        A dictionary containing configuration information, specifically the
        "plotting" key with "nice_names" and "tech_colors" keys for carriers.

    Returns
    -------
    None
        The function modifies the 'n' PyPSA Network object in-place, updating the
        carriers attribute with nice names and colors.

    Warnings
    --------
    Raises a warning if any carrier's "tech_colors" are not defined in the config dictionary.
    """

    for c in n.components:
        if "carrier" in c.static:
            add_missing_carriers(n, c.static.carrier)

    carrier_i = n.carriers.index
    nice_names = (
        pd.Series(config["plotting"]["nice_names"])
        .reindex(carrier_i)
        .fillna(carrier_i.to_series())
    )
    n.carriers["nice_name"] = n.carriers.nice_name.where(
        n.carriers.nice_name != "", nice_names
    )

    tech_colors = config["plotting"]["tech_colors"]
    colors = pd.Series(tech_colors).reindex(carrier_i)
    # try to fill missing colors with tech_colors after renaming
    missing_colors_i = colors[colors.isna()].index
    colors[missing_colors_i] = missing_colors_i.map(rename_techs).map(tech_colors)
    if colors.isna().any():
        missing_i = list(colors.index[colors.isna()])
        logger.warning(f"tech_colors for carriers {missing_i} not defined in config.")
    n.carriers["color"] = n.carriers.color.where(n.carriers.color != "", colors)


def sanitize_locations(n):
    if "location" in n.buses.columns:
        n.buses["x"] = n.buses.x.where(n.buses.x != 0, n.buses.location.map(n.buses.x))
        n.buses["y"] = n.buses.y.where(n.buses.y != 0, n.buses.location.map(n.buses.y))
        n.buses["country"] = n.buses.country.where(
            n.buses.country.ne("") & n.buses.country.notnull(),
            n.buses.location.map(n.buses.country),
        )


def estimate_efficiency(df: pd.DataFrame, config: dict) -> pd.Series:
    """
    Estimate power plant efficiencies from carrier and build year.

    Parameters
    ----------
    df : pd.DataFrame
        Power plants with columns `carrier`, `datein` and `dateretrofit`.
    config : dict
        Settings of `conventional: estimate_efficiencies`.

    Returns
    -------
    pd.Series
        Estimated efficiencies, NaN where carrier parameters or build year are missing.
    """
    p = pd.DataFrame(config["parameters"]).T.reindex(df.carrier).set_index(df.index)
    year = df.dateretrofit.combine_first(df.datein)

    eta = (p.efficiency + p.slope * (year - p.year)).clip(
        lower=p.efficiency, upper=p["max"]
    )

    age = config["reference_year"] - year - config["degradation_start"]
    return eta * (1 - age.clip(lower=0) * config["degradation_rate"])


def add_co2_emissions(n, costs, carriers):
    """
    Add CO2 emissions to the network's carriers attribute.
    """
    suptechs = n.carriers.loc[carriers].index.str.split("-").str[0]
    n.carriers.loc[carriers, "co2_emissions"] = costs.loc[
        suptechs, "CO2 intensity"
    ].values


def assign_unit_dates(
    units: pd.DataFrame, default_lifetime: pd.Series, first_horizon: int
) -> pd.DataFrame:
    """
    Set build year and lifetime of individual units from their commissioning and decommissioning years.

    Missing commissioning years are filled with the mean of the carrier, missing
    decommissioning years with the commissioning year plus the default lifetime
    (30 years where unknown). Units without commissioning year and units retired
    before the first planning horizon are dropped.

    Parameters
    ----------
    units : pd.DataFrame
        Units with columns `carrier`, `datein` and `dateout`.
    default_lifetime : pd.Series
        Lifetime per unit used where `dateout` is missing.
    first_horizon : int
        The first planning horizon.

    Returns
    -------
    pd.DataFrame
        Units with columns `build_year` and `lifetime`.
    """
    datein = units.datein.fillna(units.groupby("carrier").datein.transform("mean") // 1)
    if datein.isna().any():
        logger.warning(
            f"Dropping {datein.isna().sum()} units without commissioning year "
            f"of carriers {units.carrier[datein.isna()].unique().tolist()}."
        )
    units = units.assign(datein=datein).dropna(subset="datein")
    dateout = units.dateout.combine_first(
        units.datein + default_lifetime.reindex(units.index).fillna(30)
    )
    units = units.assign(
        build_year=units.datein.astype(int), lifetime=dateout - units.datein
    )
    retired = units.build_year + units.lifetime <= first_horizon
    logger.info(f"Dropping {retired.sum()} units retired by {first_horizon}.")
    return units.loc[~retired]


def points_to_wkt(df: pd.DataFrame) -> pd.Series:
    """Return WKT points in EPSG:4326 from the `lon` and `lat` columns."""
    points = gpd.GeoSeries(gpd.points_from_xy(df.lon, df.lat), index=df.index)
    return points.to_wkt()


def load_and_aggregate_powerplants(
    ppl_fn: str,
    costs: pd.DataFrame,
    consider_efficiency_classes: bool | list[float] = False,
    aggregation_strategies: dict = None,
    exclude_carriers: list = None,
    estimate_efficiencies: dict | None = None,
    keep_assets: bool = False,
    first_horizon: int | None = None,
) -> pd.DataFrame:
    if not aggregation_strategies:
        aggregation_strategies = {}

    if not exclude_carriers:
        exclude_carriers = []

    carrier_dict = {
        "ocgt": "OCGT",
        "ccgt": "CCGT",
        "bioenergy": "biomass",
        "ccgt, thermal": "CCGT",
        "hard coal": "coal",
    }
    tech_dict = {
        "Run-Of-River": "ror",
        "Reservoir": "hydro",
        "Pumped Storage": "PHS",
    }
    ppl = (
        pd.read_csv(ppl_fn, index_col=0, dtype={"bus": "str"})
        .powerplant.to_pypsa_names()
        .rename(columns=str.lower)
        .replace({"carrier": carrier_dict, "technology": tech_dict})
    )

    # Replace carriers "natural gas" and "hydro" with the respective technology;
    # OCGT or CCGT and hydro, PHS, or ror)
    ppl["carrier"] = ppl.carrier.where(
        ~ppl.carrier.isin(["hydro", "natural gas"]), ppl.technology
    )

    cost_columns = [
        "VOM",
        "FOM",
        "efficiency",
        "capital_cost",
        "marginal_cost",
        "fuel",
        "lifetime",
    ]
    ppl = ppl.join(costs[cost_columns], on="carrier", rsuffix="_r")

    efficiency = ppl.efficiency
    if estimate_efficiencies and estimate_efficiencies["enable"]:
        efficiency = efficiency.combine_first(
            estimate_efficiency(ppl, estimate_efficiencies)
        )
    ppl["efficiency"] = efficiency.combine_first(ppl.efficiency_r)
    ppl["marginal_cost"] = (
        ppl.carrier.map(costs.VOM) + ppl.carrier.map(costs.fuel) / ppl.efficiency
    )

    if keep_assets:
        units = assign_unit_dates(ppl, ppl.lifetime, first_horizon)
        units["geometry"] = points_to_wkt(units)
        units.index = units.bus + " " + units.carrier + " " + units.index.astype(str)
        return units

    ppl["lifetime"] = (ppl.dateout - ppl.datein).fillna(np.inf)
    ppl["build_year"] = ppl.datein.fillna(0).astype(int)

    strategies = {
        **DEFAULT_ONE_PORT_STRATEGIES,
        **{"country": "first"},
        **aggregation_strategies.get("generators", {}),
    }
    strategies = {k: v for k, v in strategies.items() if k in ppl.columns}

    to_aggregate = ~ppl.carrier.isin(exclude_carriers)
    df = ppl[to_aggregate].copy()

    if consider_efficiency_classes:
        quantiles = (
            [0.1, 0.9]
            if consider_efficiency_classes is True
            else consider_efficiency_classes
        )
        for c in df.carrier.unique():
            df_c = df.query("carrier == @c")
            thresholds = df_c.efficiency.quantile(quantiles).tolist()
            if thresholds[0] < thresholds[-1]:
                unique_thresholds, indices = np.unique(thresholds, return_index=True)
                labels = ["Q0"] + [
                    f"Q{int(quantiles[i] * 100)}" for i in sorted(indices)
                ]
                bins = [0.0] + unique_thresholds.tolist() + [1.0]
                suffix = pd.cut(df_c.efficiency, bins=bins, labels=labels).astype(str)
                df.update({"carrier": df_c.carrier + " " + suffix + " efficiency"})

    grouper = ["bus", "carrier"]
    weights = df.groupby(grouper).p_nom.transform(normed_or_uniform)

    for k, v in strategies.items():
        if v == "capacity_weighted_average":
            df[k] = df[k] * weights
            strategies[k] = pd.Series.sum

    aggregated = df.groupby(grouper, as_index=False).agg(strategies)
    aggregated.index = aggregated.bus + " " + aggregated.carrier
    aggregated.build_year = aggregated.build_year.astype(int)
    aggregated.carrier = aggregated.carrier.str.replace(
        r" Q\d+ efficiency", "", regex=True
    )

    disaggregated = ppl[~to_aggregate].copy()
    disaggregated.index = (
        disaggregated.bus
        + " "
        + disaggregated.carrier
        + " "
        + disaggregated.index.astype(str)
        + " "
        + disaggregated.name
    )
    disaggregated = disaggregated[aggregated.columns]

    return pd.concat([aggregated, disaggregated])


def attach_load(
    n: pypsa.Network,
    load_fn: str,
    scaling: float = 1.0,
) -> None:
    """
    Attach clustered load data to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the load data to.
    load_fn : str
        Path to the clustered load data file.
    scaling : float, optional
        Scaling factor for the load data, by default 1.0.
    """
    load = (
        xr.open_dataarray(load_fn)
        .to_dataframe()
        .squeeze(axis=1)
        .unstack(level="time")
        .T
    )

    logger.info(f"Load data scaled by factor {scaling}.")
    load *= scaling

    n.add("Load", load.columns, bus=load.columns, p_set=load, carrier="electricity")


def set_transmission_costs(
    n: pypsa.Network,
    costs: pd.DataFrame,
    line_length_factor: float = 1.0,
    link_length_factor: float = 1.0,
) -> None:
    """
    Set the transmission costs for lines and links in the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to set the transmission costs for.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    line_length_factor : float, optional
        Factor to scale the line length, by default 1.0.
    link_length_factor : float, optional
        Factor to scale the link length, by default 1.0.
    """
    n.lines["capital_cost"] = (
        n.lines["length"]
        * line_length_factor
        * costs.at["HVAC overhead", "capital_cost"]
    )

    if n.links.empty:
        return

    dc_b = n.links.carrier == "DC"

    # If there are no dc links, then the 'underwater_fraction' column
    # may be missing. Therefore we have to return here.
    if n.links.loc[dc_b].empty:
        return

    costs = (
        n.links.loc[dc_b, "length"]
        * link_length_factor
        * (
            (1.0 - n.links.loc[dc_b, "underwater_fraction"])
            * costs.at["HVDC overhead", "capital_cost"]
            + n.links.loc[dc_b, "underwater_fraction"]
            * costs.at["HVDC submarine", "capital_cost"]
        )
        + costs.at["HVDC inverter pair", "capital_cost"]
    )
    n.links.loc[dc_b, "capital_cost"] = costs


def attach_wind_and_solar(
    n: pypsa.Network,
    costs: pd.DataFrame,
    profile_filenames: dict,
    carriers: list | set,
    extendable_carriers: list | set,
    line_length_factor: float = 1.0,
    landfall_lengths: dict = None,
    keep_assets: bool = False,
) -> None:
    """
    Attach wind and solar generators to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the generators to.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    profile_filenames : dict
        Dictionary containing the paths to the wind and solar profiles.
    carriers : list | set
        List of renewable energy carriers to attach.
    extendable_carriers : list | set
        List of extendable renewable energy carriers.
    line_length_factor : float, optional
        Factor to scale the line length, by default 1.0.
    landfall_lengths : dict, optional
        Dictionary containing the landfall lengths for offshore wind, by default None.
    keep_assets : bool, optional
        If True, add the resource class and the class region as WKT geometry.
    """
    add_missing_carriers(n, carriers)

    if landfall_lengths is None:
        landfall_lengths = {}

    for car in carriers:
        if car == "hydro":
            continue

        landfall_length = landfall_lengths.get(car, 0.0)

        with xr.open_dataset(profile_filenames["profile_" + car]) as ds:
            if ds.indexes["bus"].empty:
                continue

            # if-statement for compatibility with old profiles
            if "year" in ds.indexes:
                ds = ds.sel(year=ds.year.min(), drop=True)

            ds = ds.stack(bus_bin=["bus", "bin"])

            supcar = car.split("-", 2)[0]
            if supcar == "offwind":
                distance = ds["average_distance"].to_pandas()
                distance.index = distance.index.map(flatten)
                submarine_cost = costs.at[car + "-connection-submarine", "capital_cost"]
                underground_cost = costs.at[
                    car + "-connection-underground", "capital_cost"
                ]
                connection_cost = line_length_factor * (
                    distance * submarine_cost + landfall_length * underground_cost
                )

                # Take 'offwind-float' capital cost for 'float', and 'offwind' capital cost for the rest ('ac' and 'dc')
                midcar = car.split("-", 2)[1]
                if midcar == "float":
                    capital_cost = (
                        costs.at[car, "capital_cost"]
                        + costs.at[car + "-station", "capital_cost"]
                        + connection_cost
                    )
                else:
                    capital_cost = (
                        costs.at["offwind", "capital_cost"]
                        + costs.at[car + "-station", "capital_cost"]
                        + connection_cost
                    )
                logger.info(
                    f"Added connection cost of {connection_cost.min():0.0f}-{connection_cost.max():0.0f} Eur/MW/a to {car}"
                )
            else:
                capital_cost = costs.at[car, "capital_cost"]

            buses = ds.indexes["bus_bin"].get_level_values("bus")
            bus_bins = ds.indexes["bus_bin"].map(flatten)

            p_nom_max = ds["p_nom_max"].to_pandas()
            p_nom_max.index = p_nom_max.index.map(flatten)

            p_max_pu = ds["profile"].to_pandas()
            p_max_pu.columns = p_max_pu.columns.map(flatten)

            n.add(
                "Generator",
                bus_bins,
                suffix=" " + car,
                bus=buses,
                carrier=car,
                p_nom_extendable=car in extendable_carriers["Generator"],
                p_nom_max=p_nom_max,
                marginal_cost=costs.at[supcar, "marginal_cost"],
                capital_cost=capital_cost,
                efficiency=costs.at[supcar, "efficiency"],
                p_max_pu=p_max_pu,
                lifetime=costs.at[supcar, "lifetime"],
            )

            if keep_assets:
                names = bus_bins + " " + car
                class_regions = gpd.read_file(profile_filenames[f"class_regions_{car}"])
                geometry = class_regions.set_index(["bus", "bin"]).geometry.to_wkt()
                geometry.index = geometry.index.map(flatten) + " " + car
                n.generators.loc[names, "resource_class"] = ds.indexes[
                    "bus_bin"
                ].get_level_values("bin")
                n.generators.loc[names, "geometry"] = geometry.reindex(names)
                n.generators.loc[names, "asset_status"] = "candidate"


def conventional_candidates(
    n: pypsa.Network, costs: pd.DataFrame, carriers: list
) -> pd.DataFrame:
    """
    Build one extendable candidate per AC bus and carrier with zero capacity.

    Parameters
    ----------
    n : pypsa.Network
        The network with the candidate buses.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    carriers : list
        The carriers of the candidates.

    Returns
    -------
    pd.DataFrame
        Candidates in the format of the power plant data, indexed by `{bus} {carrier}`.
    """
    buses = n.buses.index[n.buses.carrier == "AC"]
    index = pd.MultiIndex.from_product([buses, carriers], names=["bus", "carrier"])
    candidates = index.to_frame(index=False)
    candidates.index = candidates.bus + " " + candidates.carrier
    tech = costs.loc[candidates.carrier].set_axis(candidates.index)
    return candidates.assign(
        p_nom=0.0,
        efficiency=tech.efficiency,
        marginal_cost=tech.VOM + tech.fuel / tech.efficiency,
        capital_cost=tech.capital_cost,
        build_year=0,
        lifetime=tech.lifetime,
    )


def attach_conventional_generators(
    n: pypsa.Network,
    costs: pd.DataFrame,
    ppl: pd.DataFrame,
    conventional_carriers: list,
    extendable_carriers: dict,
    renewable_carriers: set,
    conventional_params: dict,
    conventional_inputs: dict,
    unit_commitment: pd.DataFrame = None,
    fuel_price: pd.DataFrame = None,
    keep_assets: bool = False,
):
    """
    Attach conventional generators to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the generators to.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    ppl : pd.DataFrame
        DataFrame containing the power plant data.
    conventional_carriers : list
        List of conventional energy carriers.
    extendable_carriers : dict
        Dictionary of extendable energy carriers.
    renewable_carriers : set
        Set of carriers added separately in attach_wind_and_solar().
    conventional_params : dict
        Dictionary of conventional generator parameters.
    conventional_inputs : dict
        Dictionary of conventional generator inputs.
    unit_commitment : pd.DataFrame, optional
        DataFrame containing unit commitment data, by default None.
    fuel_price : pd.DataFrame, optional
        DataFrame containing fuel price data, by default None.
    keep_assets : bool, optional
        If True, existing units of conventional carriers are fixed and one
        extendable candidate per bus is added for each extendable carrier.
    """
    if ppl.empty:
        return

    carriers = list(
        set(conventional_carriers)
        | set(extendable_carriers["Generator"]) - set(renewable_carriers)
    )

    candidates_i = pd.Index([])
    if keep_assets:
        candidates = conventional_candidates(
            n, costs, sorted(set(extendable_carriers["Generator"]) & set(carriers))
        )
        candidates_i = candidates.index
        ppl = pd.concat([ppl.query("carrier in @conventional_carriers"), candidates])

    ppl = ppl.query("carrier in @carriers").copy()

    # reduce carriers to those in power plant dataset
    carriers = list(set(carriers) & set(ppl.carrier.unique()))
    add_missing_carriers(n, carriers)
    add_co2_emissions(n, costs, carriers)

    if unit_commitment is not None:
        committable_attrs = ppl.carrier.isin(unit_commitment).to_frame("committable")
        for attr in unit_commitment.index:
            default = n.components["Generator"].defaults.loc[attr, "default"]
            committable_attrs[attr] = ppl.carrier.map(unit_commitment.loc[attr]).fillna(
                default
            )
            # These values are given per MW, but PyPSA expects them in total.
            if attr in {"start_up_cost", "shut_down_cost", "stand_by_cost"}:
                committable_attrs[attr] *= ppl.p_nom
    else:
        committable_attrs = {}

    if fuel_price is not None:
        fuel_price = fuel_price.assign(
            OCGT=fuel_price["gas"], CCGT=fuel_price["gas"]
        ).drop("gas", axis=1)
        missing_carriers = list(set(carriers) - set(fuel_price))
        fuel_price = fuel_price.assign(**costs.fuel[missing_carriers])
        fuel_price = fuel_price.reindex(ppl.carrier, axis=1)
        fuel_price.columns = ppl.index
        marginal_cost = fuel_price.div(ppl.efficiency).add(ppl.carrier.map(costs.VOM))
    else:
        marginal_cost = ppl.marginal_cost

    # Define generators using modified ppl DataFrame
    caps = ppl.groupby("carrier").p_nom.sum().div(1e3).round(2)
    logger.info(f"Adding {len(ppl)} generators with capacities [GW]pp \n{caps}")

    n.add(
        "Generator",
        ppl.index,
        carrier=ppl.carrier,
        bus=ppl.bus,
        p_nom_min=ppl.p_nom.where(ppl.carrier.isin(conventional_carriers), 0),
        p_nom=ppl.p_nom.where(ppl.carrier.isin(conventional_carriers), 0),
        p_nom_extendable=(
            ppl.index.isin(candidates_i)
            if keep_assets
            else ppl.carrier.isin(extendable_carriers["Generator"])
        ),
        efficiency=ppl.efficiency,
        marginal_cost=marginal_cost,
        capital_cost=ppl.capital_cost,
        build_year=ppl.build_year,
        lifetime=ppl.lifetime,
        **committable_attrs,
    )

    for carrier in set(conventional_params) & set(carriers):
        # Generators with technology affected
        idx = n.generators.query("carrier == @carrier").index

        for attr in list(set(conventional_params[carrier]) & set(n.generators)):
            values = conventional_params[carrier][attr]

            if f"conventional_{carrier}_{attr}" in conventional_inputs:
                # Values affecting generators of technology k country-specific
                # First map generator buses to countries; then map countries to p_max_pu
                df = pd.read_csv(
                    conventional_inputs[f"conventional_{carrier}_{attr}"], index_col=0
                )
                try:
                    df.columns = df.columns.astype(int)
                    year = n.snapshots[0].year
                    values = df[year]
                except (ValueError, TypeError, KeyError):
                    values = df.iloc[:, -1]  # take last column if year selection fails
                bus_values = n.buses.country.map(values)
                n.generators.update(
                    {attr: n.generators.loc[idx].bus.map(bus_values).dropna()}
                )
            else:
                # Single value affecting all generators of technology k indiscriminantely of country
                n.generators.loc[idx, attr] = values


def attach_existing_batteries(
    n: pypsa.Network,
    costs: pd.DataFrame,
    ppl: pd.DataFrame,
) -> None:
    """Attach existing battery storage units from the power plant dataset."""
    batt = ppl.query('carrier == "battery"')
    if batt.empty:
        return

    add_missing_carriers(n, ["battery"])
    efficiency = np.sqrt(costs.at["battery inverter", "efficiency"])
    batt["max_hours"] = batt.max_hours.fillna(batt.max_hours.median())

    n.add(
        "StorageUnit",
        batt.index,
        carrier="battery",
        bus=batt.bus,
        p_nom=batt.p_nom,
        capital_cost=batt.capital_cost,
        max_hours=batt.max_hours,
        efficiency_store=efficiency,
        efficiency_dispatch=efficiency,
        cyclic_state_of_charge=True,
    )

    stats = (
        batt.groupby("country")
        .p_nom.sum()
        .div(1e3)
        .round(3)
        .sort_values(ascending=False)
    )
    logger.info(f"Added {len(batt)} existing battery storage units\n({stats} MW)")


def attach_hydro(
    n: pypsa.Network,
    costs: pd.DataFrame,
    ppl: pd.DataFrame,
    profile_hydro: str,
    hydro_capacities: str,
    carriers: list,
    **params,
):
    """
    Attach hydro generators and storage units to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the hydro units to.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    ppl : pd.DataFrame
        DataFrame containing the power plant data.
    profile_hydro : str
        Path to the hydro profile data.
    hydro_capacities : str
        Path to the hydro capacities data.
    carriers : list
        List of hydro energy carriers.
    **params :
        Additional parameters for hydro units.
    """
    add_missing_carriers(n, carriers)
    add_co2_emissions(n, costs, carriers)

    ror = ppl.query('carrier == "ror"')
    phs = ppl.query('carrier == "PHS"')
    hydro = ppl.query('carrier == "hydro"')

    country = ppl["bus"].map(n.buses.country).rename("country")

    inflow_idx = ror.index.union(hydro.index)
    inflow_t = pd.DataFrame()
    if not inflow_idx.empty:
        dist_key = ppl.loc[inflow_idx, "p_nom"].groupby(country).transform(normed)

        with xr.open_dataarray(profile_hydro) as inflow:
            inflow_countries = pd.Index(country[inflow_idx])
            missing_c = inflow_countries.unique().difference(
                inflow.indexes["countries"]
            )
            assert missing_c.empty, (
                f"'{profile_hydro}' is missing "
                f"inflow time-series for at least one country: {', '.join(missing_c)}"
            )

            inflow_t = (
                inflow.sel(countries=inflow_countries)
                .rename({"countries": "name"})
                .assign_coords(name=inflow_idx)
                .transpose("time", "name")
                .to_pandas()
                .multiply(dist_key, axis=1)
            )

    if "ror" in carriers and not ror.empty:
        n.add(
            "Generator",
            ror.index,
            carrier="ror",
            bus=ror["bus"],
            p_nom=ror["p_nom"],
            efficiency=costs.at["ror", "efficiency"],
            capital_cost=costs.at["ror", "capital_cost"],
            weight=ror["p_nom"],
            p_max_pu=(
                inflow_t[ror.index]
                .divide(ror["p_nom"], axis=1)
                .where(lambda df: df <= 1.0, other=1.0)
            ),
        )

    if "PHS" in carriers and not phs.empty:
        # fill missing max hours to params value and
        # assume no natural inflow due to lack of data
        max_hours = params.get("PHS_max_hours", 6)
        phs = phs.replace({"max_hours": {0: max_hours, np.nan: max_hours}})
        n.add(
            "StorageUnit",
            phs.index,
            carrier="PHS",
            bus=phs["bus"],
            p_nom=phs["p_nom"],
            capital_cost=costs.at["PHS", "capital_cost"],
            max_hours=phs["max_hours"],
            efficiency_store=np.sqrt(costs.at["PHS", "efficiency"]),
            efficiency_dispatch=np.sqrt(costs.at["PHS", "efficiency"]),
            cyclic_state_of_charge=True,
        )

    if "hydro" in carriers and not hydro.empty:
        hydro_max_hours = params.get("hydro_max_hours")

        assert hydro_capacities is not None, "No path for hydro capacities given."

        hydro_stats = pd.read_csv(
            hydro_capacities, comment="#", na_values="-", index_col=0
        )
        e_target = hydro_stats["E_store[TWh]"].clip(lower=0.2) * 1e6
        e_installed = hydro.eval("p_nom * max_hours").groupby(hydro.country).sum()
        e_missing = e_target - e_installed
        missing_mh_i = hydro.query("max_hours.isnull() or max_hours == 0").index
        # some countries may have missing storage capacity but only one plant
        # which needs to be scaled to the target storage capacity
        missing_mh_single_i = hydro.index[
            ~hydro.country.duplicated() & hydro.country.isin(e_missing.dropna().index)
        ]
        missing_mh_i = missing_mh_i.union(missing_mh_single_i)

        if hydro_max_hours == "energy_capacity_totals_by_country":
            # watch out some p_nom values like IE's are totally underrepresented
            max_hours_country = (
                e_missing / hydro.loc[missing_mh_i].groupby("country").p_nom.sum()
            )

        elif hydro_max_hours == "estimate_by_large_installations":
            max_hours_country = (
                hydro_stats["E_store[TWh]"] * 1e3 / hydro_stats["p_nom_discharge[GW]"]
            )
        else:
            raise ValueError(f"Unknown hydro_max_hours method: {hydro_max_hours}")

        max_hours_country.clip(0, inplace=True)

        missing_countries = pd.Index(hydro["country"].unique()).difference(
            max_hours_country.dropna().index
        )
        if not missing_countries.empty:
            logger.warning(
                f"Assuming max_hours=6 for hydro reservoirs in the countries: {', '.join(missing_countries)}"
            )
        hydro_max_hours = hydro.max_hours.where(
            (hydro.max_hours > 0) & ~hydro.index.isin(missing_mh_single_i),
            hydro.country.map(max_hours_country),
        ).fillna(6)

        if params.get("flatten_dispatch", False):
            buffer = params.get("flatten_dispatch_buffer", 0.2)
            average_capacity_factor = inflow_t[hydro.index].mean() / hydro["p_nom"]
            p_max_pu = (average_capacity_factor + buffer).clip(upper=1)
        else:
            p_max_pu = 1

        n.add(
            "StorageUnit",
            hydro.index,
            carrier="hydro",
            bus=hydro["bus"],
            p_nom=hydro["p_nom"],
            max_hours=hydro_max_hours,
            capital_cost=costs.at["hydro", "capital_cost"],
            marginal_cost=costs.at["hydro", "marginal_cost"],
            p_max_pu=p_max_pu,  # dispatch
            p_min_pu=0.0,  # store
            efficiency_dispatch=costs.at["hydro", "efficiency"],
            efficiency_store=0.0,
            cyclic_state_of_charge=True,
            inflow=inflow_t.loc[:, hydro.index],
        )


def renewable_powerplants(
    smk_inputs: dict, tech_map: dict[str, str]
) -> dict[str, gpd.GeoDataFrame]:
    """
    Read renewable power plants and join them to the resource class regions.

    Parameters
    ----------
    smk_inputs : dict
        The rule inputs with `powerplants` and `class_regions_{carrier}`.
    tech_map : dict[str, str]
        A dictionary mapping fuel types to carrier names.

    Returns
    -------
    dict[str, gpd.GeoDataFrame]
        Power plants with columns `bus` and `bin` of their class region, by carrier.
    """
    df = pd.read_csv(smk_inputs["powerplants"], index_col=0).drop("bus", axis=1)
    technology_b = ~df.Technology.isin(["Onshore", "Offshore"])
    df["Fueltype"] = df.Fueltype.where(technology_b, df.Technology).replace(
        {"Solar": "PV"}
    )

    plants = {}
    for fueltype, carrier in tech_map.items():
        class_regions = gpd.read_file(smk_inputs.get(f"class_regions_{carrier}"))
        df_fueltype = df.query("Fueltype == @fueltype")
        geometry = gpd.points_from_xy(df_fueltype.lon, df_fueltype.lat)
        caps = gpd.GeoDataFrame(df_fueltype, geometry=geometry, crs=4326)
        plants[carrier] = caps.sjoin(class_regions)
    return plants


def add_fixed_copies_of_candidates(
    n: pypsa.Network, candidates: pd.Index, names: pd.Index, **attrs
) -> None:
    """
    Add non-extendable generators with the attributes and profiles of renewable candidates.

    Parameters
    ----------
    n : pypsa.Network
        The network with the candidates.
    candidates : pd.Index
        The candidate of each new generator.
    names : pd.Index
        The names of the new generators.
    **attrs
        Further attributes of the new generators, aligned with `names`.
    """
    static = n.generators.loc[candidates].set_axis(names)
    p_max_pu = n.generators_t.p_max_pu[candidates].set_axis(names, axis=1)
    n.add(
        "Generator",
        names,
        bus=static.bus,
        carrier=static.carrier,
        marginal_cost=static.marginal_cost,
        capital_cost=static.capital_cost,
        efficiency=static.efficiency,
        resource_class=static.resource_class,
        p_max_pu=p_max_pu,
        p_nom_extendable=False,
        asset_status="existing",
        **attrs,
    )


def attach_renewable_units(
    n: pypsa.Network,
    tech_map: dict[str, str],
    smk_inputs: dict,
    first_horizon: int,
) -> None:
    """
    Attach each renewable power plant as a fixed generator with the profile of its resource class.

    Parameters
    ----------
    n : pypsa.Network
        The network with the renewable candidates.
    tech_map : dict[str, str]
        A dictionary mapping fuel types to carrier names.
    smk_inputs : dict
        The rule inputs with `powerplants` and `class_regions_{carrier}`.
    first_horizon : int
        The first planning horizon.
    """
    for carrier, plants in renewable_powerplants(smk_inputs, tech_map).items():
        plants = plants.rename(columns=str.lower).assign(carrier=carrier)
        candidates = plants.bus + " " + plants.bin.astype(str) + " " + carrier
        lifetime = candidates.map(n.generators.lifetime)
        units = assign_unit_dates(plants, lifetime, first_horizon)
        candidates = candidates[units.index]
        names = candidates + " " + units.index.astype(str)
        logger.info(f"Adding {len(units)} {carrier} units.")
        add_fixed_copies_of_candidates(
            n,
            pd.Index(candidates),
            pd.Index(names),
            p_nom=units.capacity.set_axis(names),
            build_year=units.build_year.set_axis(names),
            lifetime=units.lifetime.set_axis(names),
            geometry=points_to_wkt(units).set_axis(names),
        )


def attach_renewable_powerplants(
    n: pypsa.Network, tech_map: dict[str, list[str]], smk_inputs: list[str]
) -> None:
    """
    Attach renewable capacities from the GEM dataset to the network.

    Args:
    - n: The PyPSA network to attach the capacities to.
    - tech_map: A dictionary mapping fuel types to carrier names.

    Returns:
    - None
    """
    tech_string = ", ".join(tech_map.values())
    logger.info(
        f"Using powerplantmatching renewable capacities for carriers {tech_string}."
    )

    for carrier, caps in renewable_powerplants(smk_inputs, tech_map).items():
        caps = caps.groupby(["bus", "bin"]).Capacity.sum()
        caps.index = caps.index.map(flatten) + " " + carrier

        n.generators.update({"p_nom": caps.dropna()})
        n.generators.update({"p_nom_min": caps.dropna()})


def irena_capacities(year: int, tech_map: dict, countries: list) -> pd.Series:
    """Return IRENASTAT capacities by technology and country for a year."""
    capacities = pm.data.IRENASTAT().powerplant.convert_country_to_alpha2()
    capacities = capacities.query(
        "Year == @year and Technology in @tech_map and Country in @countries"
    )
    capacities = capacities.groupby(["Technology", "Country"]).Capacity.sum()

    logger.info(
        f"Heuristics applied to distribute renewable capacities [GW]: "
        f"\n{capacities.groupby('Technology').sum().div(1e3).round(2)}"
    )
    return capacities


def country_capacities(
    capacities: pd.Series, ppm_technology: str, countries: list
) -> pd.Series:
    """Return the capacities of one technology for all countries, zero where missing."""
    if ppm_technology in capacities.index.get_level_values("Technology"):
        return capacities.loc[ppm_technology].reindex(countries, fill_value=0.0)
    return pd.Series(0.0, index=countries)


def attach_renewable_estimates(
    n: pypsa.Network, year: int, tech_map: dict, countries: list
) -> None:
    """
    Attach the IRENASTAT capacities not covered by renewable units as fixed estimate generators.

    The remainder per country is distributed to the candidates of the country
    with the heuristic `n.generators_t.p_max_pu.mean() * n.generators.p_nom_max`.
    Each estimate is named `{candidate} estimate`, is built in `year` and has
    the lifetime of its candidate.

    Parameters
    ----------
    n : pypsa.Network
        The network with renewable units and candidates.
    year : int
        The year of the IRENASTAT capacities.
    tech_map : dict
        A dictionary mapping fuel types to carrier names.
    countries : list
        The country codes to estimate capacities for.
    """
    if not len(countries) or not len(tech_map):
        return

    capacities = irena_capacities(year, tech_map, countries)
    country = n.generators.bus.map(n.buses.country)
    dist = n.generators_t.p_max_pu.mean() * n.generators.p_nom_max

    for ppm_technology, tech in tech_map.items():
        is_tech = n.generators.carrier == tech
        is_candidate = n.generators.asset_status == "candidate"
        candidates_i = n.generators.index[is_tech & is_candidate]
        units_i = n.generators.index[is_tech & ~is_candidate]
        stats = country_capacities(capacities, ppm_technology, countries)
        existent = n.generators.p_nom[units_i].groupby(country[units_i]).sum()
        missing = stats - existent.reindex(countries, fill_value=0.0)
        p_nom = (
            dist[candidates_i]
            .groupby(country[candidates_i])
            .transform(lambda s: normed(s) * missing[s.name])
        )
        p_nom = p_nom[p_nom > 0.1]
        logger.info(f"Adding {len(p_nom)} {tech} estimates.")
        if p_nom.empty:
            continue
        names = p_nom.index + " estimate"
        add_fixed_copies_of_candidates(
            n,
            p_nom.index,
            names,
            p_nom=p_nom.set_axis(names),
            build_year=year,
            lifetime=n.generators.lifetime[p_nom.index].set_axis(names),
        )
        n.generators.loc[names, "asset_status"] = "estimate"


def reduce_candidate_potentials(n: pypsa.Network, carriers: list) -> None:
    """
    Reduce the potential of renewable candidates by the capacity of fixed generators in their class.

    Parameters
    ----------
    n : pypsa.Network
        The network with renewable units, estimates and candidates.
    carriers : list
        The renewable carriers.
    """
    gens = n.generators.query("carrier in @carriers")
    candidates = gens.query("asset_status == 'candidate'")
    fixed = gens.query("asset_status != 'candidate'")
    stem = fixed.bus + " " + fixed.resource_class.astype(int).astype(str)
    used = fixed.p_nom.groupby(stem + " " + fixed.carrier).sum()
    n.generators.loc[candidates.index, "p_nom_max"] = (
        candidates.p_nom_max - used.reindex(candidates.index, fill_value=0.0)
    ).clip(lower=0.0)


def estimate_renewable_capacities(
    n: pypsa.Network,
    year: int,
    tech_map: dict,
    expansion_limit: bool,
    countries: list,
):
    """
    Estimate a different between renewable capacities in the network and
    reported country totals from IRENASTAT dataset. Distribute the difference
    with a heuristic.

    Heuristic: n.generators_t.p_max_pu.mean() * n.generators.p_nom_max

    Args:
    - n: The PyPSA network.
    - year: The year of optimisation.
    - tech_map: A dictionary mapping fuel types to carrier names.
    - expansion_limit: Boolean value from config file
    - countries: A list of country codes to estimate capacities for.

    Returns:
    - None
    """
    if not len(countries) or not len(tech_map):
        return

    capacities = irena_capacities(year, tech_map, countries)

    for ppm_technology, tech in tech_map.items():
        tech_i = n.generators.query("carrier == @tech").index
        stats = country_capacities(capacities, ppm_technology, countries)
        country = n.generators.bus[tech_i].map(n.buses.country)
        existent = n.generators.p_nom[tech_i].groupby(country).sum()
        missing = stats - existent
        dist = n.generators_t.p_max_pu.mean() * n.generators.p_nom_max

        n.generators.loc[tech_i, "p_nom"] += (
            dist[tech_i]
            .groupby(country)
            .transform(lambda s: normed(s) * missing[s.name])
            .where(lambda s: s > 0.1, 0.0)  # only capacities above 100kW
        )
        n.generators.loc[tech_i, "p_nom_min"] = n.generators.loc[tech_i, "p_nom"]

        if expansion_limit:
            assert np.isscalar(expansion_limit)
            logger.info(
                f"Reducing capacity expansion limit to {expansion_limit * 100:.2f}% of installed capacity."
            )
            n.generators.loc[tech_i, "p_nom_max"] = (
                expansion_limit * n.generators.loc[tech_i, "p_nom_min"]
            )


def get_available_storage_carriers(carriers):
    """
    Filter and register available storage carriers from a given list.
    """
    implemented = set(STORE_LOOKUP.keys())
    input_carriers = set(carriers)

    not_implemented = input_carriers - implemented
    if not_implemented:
        logger.warning(
            "The following carriers are not implemented as storage technologies in PyPSA-Eur and will be skipped:\n - "
            + "\n - ".join(sorted(not_implemented))
        )

    available_carriers = sorted(input_carriers & implemented)

    return available_carriers


def attach_storageunits(
    n: pypsa.Network,
    costs: pd.DataFrame,
    buses_i: list,
    extendable_carriers: list,
    max_hours: dict,
):
    """
    Attach storage units to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the storage units to.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    buses_i : list
        List of high voltage electricity buses.
    extendable_carriers : list
        List of extendable storage units carrier names.
    max_hours : dict
        Dictionary of maximum hours for storage units.
    """
    available_carriers = get_available_storage_carriers(extendable_carriers)
    n.add("Carrier", available_carriers)

    for carrier in available_carriers:
        max_hour = max_hours.get(carrier)
        if max_hour is None:
            logger.warning(f"No max_hours defined for carrier '{carrier}'. Skipping.")
            continue

        lookup = STORE_LOOKUP[carrier]
        if "bicharger" in lookup:
            lookup_charge = lookup_discharge = lookup["bicharger"]
        else:
            lookup_charge = lookup["charger"]
            lookup_discharge = lookup["discharger"]

        roundtrip_correction = lookup.get("roundtrip_correction", 1)

        efficiency_store = costs.at[lookup_charge, "efficiency"] ** roundtrip_correction
        efficiency_dispatch = (
            costs.at[lookup_discharge, "efficiency"] ** roundtrip_correction
        )

        # A dispatched-basis `max_hours` counts hours at rated output, so the store must
        # be larger by 1/efficiency_dispatch
        if lookup.get("energy_basis") == "dispatched":
            max_hour /= efficiency_dispatch
            logger.info(
                f"'{carrier}' max_hours counts hours at rated output; sizing its store "
                f"to {max_hour:.2f} h at efficiency_dispatch={efficiency_dispatch:.2f}."
            )

        n.add(
            "StorageUnit",
            buses_i,
            " " + carrier,
            bus=buses_i,
            carrier=carrier,
            p_nom_extendable=True,
            capital_cost=costs.at[carrier, "capital_cost"],
            marginal_cost=costs.at[carrier, "marginal_cost"],
            efficiency_store=efficiency_store,
            efficiency_dispatch=efficiency_dispatch,
            max_hours=max_hour,
            cyclic_state_of_charge=True,
            lifetime=costs.at[carrier, "lifetime"],
        )

    logger.info(
        "Add the following technologies as storage units:\n - "
        + "\n - ".join(available_carriers)
    )


def attach_stores(
    n: pypsa.Network,
    costs: pd.DataFrame,
    buses_i: list,
    extendable_carriers: list,
):
    """
    Attach stores to the network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to attach the stores to.
    costs : pd.DataFrame
        DataFrame containing the cost data.
    buses_i : list
        List of high voltage electricity buses.
    extendable_carriers : list
        List of extendable storage carrier names.
    """
    available_carriers = get_available_storage_carriers(extendable_carriers)
    n.add("Carrier", available_carriers)

    for carrier in available_carriers:
        lookup = STORE_LOOKUP[carrier]
        lookup_store = lookup["store"]
        if "bicharger" in lookup:
            lookup_charge = lookup_discharge = lookup["bicharger"]
        else:
            lookup_charge = lookup["charger"]
            lookup_discharge = lookup["discharger"]

        roundtrip_correction = lookup.get("roundtrip_correction", 1)

        bus_names = buses_i + f" {carrier}"
        charge_name = "Electrolysis" if lookup_charge == "electrolysis" else "charger"
        discharge_name = (
            "Fuel Cell" if lookup_discharge == "fuel cell" else "discharger"
        )

        discharge_capital_cost = (
            0.0 if "bicharger" in lookup else costs.at[lookup_discharge, "capital_cost"]
        )
        if lookup_discharge == "fuel cell":
            # NB: fuel cell investment cost is per MWel
            discharge_capital_cost *= costs.at[lookup_discharge, "efficiency"]

        store_capital_cost = costs.at[lookup_store, "capital_cost"]
        if lookup.get("energy_basis") == "dispatched":
            # convert cost per MWh dispatched to cost per MWh stored
            store_capital_cost *= (
                costs.at[lookup_discharge, "efficiency"] ** roundtrip_correction
            )

        n.add(
            "Bus",
            bus_names,
            location=buses_i,
            carrier=carrier,
            x=n.buses.loc[list(buses_i)].x.values,
            y=n.buses.loc[list(buses_i)].y.values,
        )

        n.add(
            "Store",
            bus_names,
            bus=bus_names,
            e_cyclic=True,
            e_nom_extendable=True,
            carrier=carrier,
            capital_cost=store_capital_cost,
            lifetime=costs.at[lookup_store, "lifetime"],
        )

        n.add("Carrier", [f"{carrier} {charge_name}", f"{carrier} {discharge_name}"])

        n.add(
            "Link",
            bus_names,
            suffix=f" {charge_name}",
            bus0=buses_i,
            bus1=bus_names,
            carrier=f"{carrier} {charge_name}",
            efficiency=costs.at[lookup_charge, "efficiency"] ** roundtrip_correction,
            capital_cost=costs.at[lookup_charge, "capital_cost"],
            p_nom_extendable=True,
            marginal_cost=costs.at[lookup_charge, "marginal_cost"],
            lifetime=costs.at[lookup_charge, "lifetime"],
        )

        n.add(
            "Link",
            bus_names,
            suffix=f" {discharge_name}",
            bus0=bus_names,
            bus1=buses_i,
            carrier=f"{carrier} {discharge_name}",
            efficiency=costs.at[lookup_discharge, "efficiency"] ** roundtrip_correction,
            capital_cost=discharge_capital_cost,
            p_nom_extendable=True,
            marginal_cost=costs.at[lookup_discharge, "marginal_cost"],
            lifetime=costs.at[lookup_discharge, "lifetime"],
        )

    logger.info(
        "Add the following storage technologies as stores and links:\n - "
        + "\n - ".join(available_carriers)
    )


def main(
    n: pypsa.Network,
    inputs,
    params,
    costs: pd.DataFrame,
) -> None:
    logger.info("Adding electricity components")
    foresight = params.foresight
    keep_assets = params.clustering["keep_assets"]
    first_horizon = params.horizons[0]

    landfall_lengths = {
        tech: settings["landfall_length"]
        for tech, settings in params.renewable.items()
        if "landfall_length" in settings.keys()
    }

    if inputs.powerplants:
        ppl = load_and_aggregate_powerplants(
            inputs.powerplants,
            costs,
            consider_efficiency_classes=params.clustering[
                "consider_efficiency_classes"
            ],
            aggregation_strategies=params.clustering["aggregation_strategies"],
            exclude_carriers=params.clustering["exclude_carriers"],
            estimate_efficiencies=params.conventional["estimate_efficiencies"],
            keep_assets=keep_assets,
            first_horizon=first_horizon,
        )
    else:
        ppl = pd.DataFrame()

    attach_load(
        n,
        inputs["load"],
        params.load["scaling_factor"],
    )

    set_transmission_costs(
        n,
        costs,
        params.lines["length_factor"],
        params.links["length_factor"],
    )

    renewable_carriers = set(params.electricity["renewable_carriers"])
    extendable_carriers = params.electricity["extendable_carriers"]
    conventional_carriers = params.electricity["conventional_carriers"]
    conventional_inputs = {
        k: v for k, v in inputs.items() if k.startswith("conventional_")
    }

    if params.conventional["unit_commitment"]:
        unit_commitment = pd.read_csv(inputs.unit_commitment, index_col=0)
    else:
        unit_commitment = None

    if params.conventional["dynamic_fuel_price"]:
        fuel_price = pd.read_csv(
            inputs.fuel_price, index_col=0, header=0, parse_dates=True
        )
        fuel_price = fuel_price.reindex(n.snapshots).ffill()
    else:
        fuel_price = None

    attach_conventional_generators(
        n,
        costs,
        ppl,
        conventional_carriers,
        extendable_carriers,
        renewable_carriers,
        conventional_params=params.conventional,
        conventional_inputs=conventional_inputs,
        unit_commitment=unit_commitment,
        fuel_price=fuel_price,
        keep_assets=keep_assets,
    )

    attach_wind_and_solar(
        n,
        costs,
        inputs,
        renewable_carriers,
        extendable_carriers,
        params.lines["length_factor"],
        landfall_lengths,
        keep_assets=keep_assets,
    )

    if "hydro" in renewable_carriers:
        p = params.renewable["hydro"].copy()
        carriers = p.pop("carriers", [])
        attach_hydro(
            n,
            costs,
            ppl,
            inputs.profile_hydro,
            inputs.hydro_capacities,
            carriers,
            **p,
        )

    estimate_renewable_caps = params.electricity["estimate_renewable_capacities"]
    if keep_assets and estimate_renewable_caps["enable"]:
        tech_map = estimate_renewable_caps["technology_mapping"]
        if estimate_renewable_caps["from_powerplantmatching"]:
            attach_renewable_units(n, tech_map, inputs, first_horizon)
        if estimate_renewable_caps["from_irenastat"]:
            attach_renewable_estimates(
                n, estimate_renewable_caps["year"], tech_map, params.countries
            )
        reduce_candidate_potentials(n, sorted(renewable_carriers - {"hydro"}))
    elif estimate_renewable_caps["enable"]:
        if foresight != "overnight":
            logger.info(
                "Skipping renewable capacity estimation because they are added later "
                "in add_existing_baseyear with foresight mode 'myopic'."
            )
        else:
            tech_map = estimate_renewable_caps["technology_mapping"]
            expansion_limit = estimate_renewable_caps["expansion_limit"]
            year = estimate_renewable_caps["year"]

            if estimate_renewable_caps["from_powerplantmatching"]:
                attach_renewable_powerplants(n, tech_map, inputs)

            if estimate_renewable_caps["from_irenastat"]:
                estimate_renewable_capacities(
                    n, year, tech_map, expansion_limit, params.countries
                )

    update_p_nom_max(n)

    max_hours = params.electricity["max_hours"]
    attach_storageunits(
        n, costs, n.buses.index, extendable_carriers["StorageUnit"], max_hours
    )
    attach_stores(n, costs, n.buses.index, extendable_carriers["Store"])

    if params.electricity["estimate_battery_capacities"]:
        attach_existing_batteries(n, costs, ppl)

    if keep_assets and not ppl.empty:
        for c in n.components[["Generator", "StorageUnit"]]:
            units_i = c.static.index.intersection(ppl.index)
            c.static.loc[units_i, "geometry"] = ppl.geometry[units_i]

    logger.info("Completed electricity components")
