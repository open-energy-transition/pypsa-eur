import pandas as pd
import pypsa
import pycountry
import matplotlib.pyplot as plt
import seaborn as sns
import geopandas as gpd
import numpy as np
import os
from pathlib import Path
import cartopy.crs as ccrs
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches

# Configuration
def get_config():
    return snakemake.config if 'snakemake' in globals() else {'year': 2023}

config = get_config()
year = config.get('year', 2023)
network_path = snakemake.input.network if 'snakemake' in globals() else f"results/validation_{year}/networks/base_s_39_elec_.nc"
ember_capacity_data_path = snakemake.input.ember_capacity if 'snakemake' in globals() else f"validation/ember_data/yearly_full_release_long_format.csv"
ember_demand_data_path = snakemake.input.ember_demand if 'snakemake' in globals() else f"validation/ember_data/europe_monthly_full_release_long_format.csv"
regions_onshore_path = snakemake.input.regions_onshore if 'snakemake' in globals() else f"resources/validation_{year}/country_shapes.geojson"

# Plotting styles
plt.style.use("bmh")
sns.set_style("darkgrid")

# Countries
countries = ['AL', 'AT', 'BA', 'BE', 'BG', 'CH', 'CZ', 'DE', 'DK', 'EE', 'ES',
             'FI', 'FR', 'GB', 'GR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV',
             'ME', 'MK', 'NL', 'NO', 'PL', 'PT', 'RO', 'RS', 'SE', 'SI', 'SK']

# Load data
try:
    n = pypsa.Network(network_path)
    print(f"Loaded PyPSA network: {network_path}")
except Exception as e:
    print(f"Error loading PyPSA network: {e}")
    raise

try:
    ember_capacity = pd.read_csv(ember_capacity_data_path)
    print(f"Loaded Ember capacity data: {ember_capacity_data_path}")
except Exception as e:
    print(f"Error loading Ember capacity data: {e}")
    raise

try:
    ember_demand = pd.read_csv(ember_demand_data_path)
    print(f"Loaded Ember demand data: {ember_demand_data_path}")
except Exception as e:
    print(f"Error loading Ember demand data: {e}")
    raise

try:
    regions_onshore = gpd.read_file(regions_onshore_path)
    print(f"Loaded regions onshore: {regions_onshore_path}")
except Exception as e:
    print(f"Error loading regions onshore: {e}")
    raise

# Process Ember capacity data
def process_ember_capacity():
    print("Processing Ember capacity data")
    exclude_areas = ["Belarus", "Gibraltar", "Iceland", "Kosovo", "Moldova", "Russian Federation (the)", "Malta", "Cyprus"]
    ember_capacity_filtered = ember_capacity.query(
        f"Year == {year} and Continent == 'Europe' and Category == 'Capacity' and Subcategory == 'Fuel'"
    ).copy()
    ember_capacity_filtered = ember_capacity_filtered[~ember_capacity_filtered['Area'].isin(exclude_areas)]
    
    def iso3_to_iso2(iso3):
        try:
            return pycountry.countries.get(alpha_3=iso3).alpha_2
        except:
            return None
    
    ember_capacity_filtered["ISO"] = ember_capacity_filtered["ISO 3 code"].apply(iso3_to_iso2)
    ember_capacity_filtered = ember_capacity_filtered[["ISO", "Variable", "Value", "Unit"]]
    ember_capacity_pivot = ember_capacity_filtered.pivot_table(index="ISO", columns="Variable", values="Value").fillna(0)
    print(f"Processed Ember capacity data sample:\n{ember_capacity_pivot.head().to_string()}")
    return ember_capacity_pivot

# Process Ember demand data
def process_ember_demand():
    print("Processing Ember demand data")
    ember_demand_europe = ember_demand[ember_demand["Continent"] == "Europe"].copy()
    
    def iso3_to_iso2(iso3):
        try:
            return pycountry.countries.get(alpha_3=iso3).alpha_2
        except:
            return None
    
    ember_demand_europe["ISO"] = ember_demand_europe["ISO 3 code"].apply(iso3_to_iso2)
    ember_demand_europe = ember_demand_europe[ember_demand_europe["Date"].str.startswith(str(year))]
    ember_demand_europe = ember_demand_europe[ember_demand_europe["Unit"] == "TWh"]
    ember_demand_europe = ember_demand_europe[ember_demand_europe["Subcategory"] == "Demand"]
    ember_demand_europe = ember_demand_europe[["ISO", "Date", "Variable", "Value", "Unit"]]
    
    ember_demand_yearly = ember_demand_europe.groupby(["ISO", "Variable"], as_index=False)["Value"].sum()
    ember_demand_yearly["Unit"] = "TWh"
    ember_demand_yearly = ember_demand_yearly.set_index("ISO").drop(columns=["Unit", "Variable"])
    print(f"Processed Ember demand data sample:\n{ember_demand_yearly.head().to_string()}")
    return ember_demand_yearly

# Process PyPSA capacity data
def process_pypsa_capacity():
    print("Processing PyPSA capacity data")
    def merge_and_replace(df, new_col, cols_to_merge, drop_original=True):
        df[new_col] = df[cols_to_merge].sum(axis=1)
        if drop_original:
            df = df.drop(columns=cols_to_merge, errors="ignore")
        return df
    
    pypsa_country_tech = n.generators.groupby(['bus', 'carrier'])['p_nom'].sum().unstack(fill_value=0).reset_index()
    pypsa_country_tech['country'] = pypsa_country_tech['bus'].str[:2]
    
    sto_country_tech = n.storage_units.groupby(['bus', 'carrier'])['p_nom'].sum().unstack(fill_value=0).reset_index()
    sto_country_tech['country'] = sto_country_tech['bus'].str[:2]
    
    pypsa_country_tech = pd.concat([pypsa_country_tech, sto_country_tech], ignore_index=True)
    pypsa_country_tech['country'] = pypsa_country_tech['bus'].str[:2]
    
    pypsa_country_tech_merged = pypsa_country_tech.groupby('country').sum(numeric_only=True).reset_index()
    
    wind_cols = [col for col in pypsa_country_tech_merged.columns if col.startswith("onwind") or col.startswith("offwind")]
    solar_cols = [col for col in pypsa_country_tech_merged.columns if col.startswith("solar")]
    
    pypsa_country_tech_merged = merge_and_replace(pypsa_country_tech_merged, "Wind", wind_cols)
    pypsa_country_tech_merged = merge_and_replace(pypsa_country_tech_merged, "Solar", solar_cols)
    pypsa_country_tech_merged = merge_and_replace(pypsa_country_tech_merged, "Gas", ["CCGT", "OCGT"])
    pypsa_country_tech_merged = merge_and_replace(pypsa_country_tech_merged, "Coal", ["coal", "lignite"])
    pypsa_country_tech_merged = merge_and_replace(pypsa_country_tech_merged, "Hydro", ["hydro", "ror", "PHS"])
    pypsa_country_tech_merged = merge_and_replace(pypsa_country_tech_merged, "Other Fossil", ["oil"])
    pypsa_country_tech_merged = merge_and_replace(pypsa_country_tech_merged, "Nuclear", ["nuclear"])
    pypsa_country_tech_merged = merge_and_replace(pypsa_country_tech_merged, "Other Renewables", ["geothermal"])
    pypsa_country_tech_merged = merge_and_replace(pypsa_country_tech_merged, "Bioenergy", ["biomass"])
    
    pypsa_country_tech_merged = pypsa_country_tech_merged.set_index("country").div(1000).round(2)
    print(f"Processed PyPSA capacity data sample:\n{pypsa_country_tech_merged.head().to_string()}")
    return pypsa_country_tech_merged

# Process PyPSA demand data
def process_pypsa_demand():
    print("Processing PyPSA demand data")
    pypsa_loads = n.loads_t.p.groupby(lambda col: col[:2], axis=1).sum()
    pypsa_loads_yearly = pypsa_loads.sum(axis=0).to_frame(name="Value")
    pypsa_loads_yearly = (pypsa_loads_yearly / 1e6).round()
    print(f"Processed PyPSA demand data sample:\n{pypsa_loads_yearly.head().to_string()}")
    return pypsa_loads_yearly

# Plot network
def plot_network():
    print("Generating network plot")
    n.carriers['color'] = n.carriers['color'].fillna('gray')
    n.carriers.loc[n.carriers['color'] == '', 'color'] = 'gray'
    
    bus_scale = 5e4
    line_scale = 1e3
    bus_sizes = [10000, 50000]
    line_sizes = [5000, 10000]
    
    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={"projection": ccrs.EqualEarth(n.buses.x.mean())})
    
    gen = n.generators[n.generators.carrier != "load"].groupby(["bus", "carrier"]).p_nom.sum()
    sto = n.storage_units.groupby(["bus", "carrier"]).p_nom.sum()
    buses = pd.concat([gen, sto])
    
    with plt.rc_context({"patch.linewidth": 0.}):
        n.plot(
            bus_sizes=buses / bus_scale,
            bus_alpha=1.0,
            line_widths=n.lines.s_nom_opt / line_scale,
            link_widths=n.links.p_nom_opt / line_scale,
            line_colors="teal",
            ax=ax,
            margin=0.2,
            color_geomap=None,
        )
    
    regions_onshore.plot(
        ax=ax,
        facecolor="whitesmoke",
        edgecolor="white",
        aspect="equal",
        transform=ccrs.PlateCarree(),
        linewidth=0,
    )
    
    ax.set_extent(regions_onshore.total_bounds[[0, 2, 1, 3]])
    
    legend_kwargs = {"loc": "upper left", "frameon": False}
    legend_circles_dict = {"bbox_to_anchor": (1, 0.67), "labelspacing": 2.5, **legend_kwargs}
    
    add_legend_circles(
        ax,
        [s / bus_scale for s in bus_sizes],
        [f"{s / 1000} GW" for s in bus_sizes],
        legend_kw=legend_circles_dict,
    )
    
    add_legend_lines(
        ax,
        [s / line_scale for s in line_sizes],
        [f"{s / 1000} GW" for s in line_sizes],
        legend_kw={"bbox_to_anchor": (1, 0.8), **legend_kwargs},
    )
    
    existing_carriers = n.generators.carrier.unique()
    carriers_legend = n.carriers.loc[n.carriers.index.isin(existing_carriers)]
    
    add_legend_patches(
        ax,
        carriers_legend.color,
        carriers_legend.nice_name,
        legend_kw={"bbox_to_anchor": (1, 0), **legend_kwargs, "loc": "lower left"},
    )
    
    output_path = snakemake.output.network_plot if 'snakemake' in globals() else f"results/validation_{year}/plots/network_plot.png"
    print(f"Saving network plot: {output_path}")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

# Plot total capacity comparison
def plot_total_capacity_comparison(ember_capacity, pypsa_capacity):
    print("Generating total capacity comparison plot")
    ember_totals = ember_capacity.sum(axis=0)
    pypsa_totals = pypsa_capacity.sum(axis=0)
    
    all_techs = sorted(set(ember_totals.index) | set(pypsa_totals.index))
    ember_vals = [ember_totals.get(tech, 0) for tech in all_techs]
    pypsa_vals = [pypsa_totals.get(tech, 0) for tech in all_techs]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    width = 0.35
    x = np.arange(len(all_techs))
    
    ax.bar(x - width/2, ember_vals, width, label='Ember', color='#13ce74')
    ax.bar(x + width/2, pypsa_vals, width, label='PyPSA-Eur', color='#192238')
    
    ax.set_xticks(x)
    ax.set_xticklabels(all_techs, rotation=45, ha='right')
    ax.set_ylabel('Installed Capacity (GW)')
    ax.set_title('Comparison of Installed Capacity by Technology (Summed Across Europe)')
    ax.legend()
    plt.tight_layout()
    
    output_path = snakemake.output.total_capacity_plot if 'snakemake' in globals() else f"results/validation_{year}/plots/total_capacity_plot.png"
    print(f"Saving total capacity comparison plot: {output_path}")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

# Plot country-specific capacity comparison
def plot_capacity_comparison_horizontal(countries, ember_capacity, pypsa_capacity):
    print("Generating country-specific capacity comparison plot")
    n = min(len(countries), 6)
    fig, axes = plt.subplots(3, 2, figsize=(12, 12))
    axes = axes.flatten()
    
    for idx, country in enumerate(countries[:6]):
        ax = axes[idx]
        ember_row = ember_capacity.loc[country] if country in ember_capacity.index else None
        pypsa_row = pypsa_capacity.loc[country] if country in pypsa_capacity.index else None
        
        if ember_row is not None and pypsa_row is not None:
            techs = sorted(set(ember_row.index) | set(pypsa_row.index))
            techs = [t if t != "Other Renewables" else "Other RES" for t in techs]
            
            ember_vals = [ember_row.get(tech, 0) for tech in techs]
            pypsa_vals = [pypsa_row.get(tech, 0) for tech in techs]
            
            y = np.arange(len(techs))
            height = 0.35
            
            ax.barh(y - height/2, ember_vals, height, label='Ember', color='#13ce74')
            ax.barh(y + height/2, pypsa_vals, height, label='PyPSA-Eur', color='#192238')
            ax.set_yticks(y)
            ax.set_yticklabels(techs)
            ax.set_title(f"{country} Capacity Comparison")
            
            if idx in [4, 5]:
                ax.set_xlabel("Capacity (GW)")
            
            if idx == 0:
                ax.legend(loc="upper right")
    
    for j in range(len(countries[:6]), len(axes)):
        axes[j].axis('off')
    
    fig.suptitle(f"Installed Capacity Comparison {year}", weight="bold")
    plt.tight_layout()
    
    output_path = snakemake.output.country_capacity_plot if 'snakemake' in globals() else f"results/validation_{year}/plots/country_capacity_plot.png"
    print(f"Saving country-specific capacity comparison plot: {output_path}")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

# Plot demand comparison
def plot_demand_comparison(ember_demand, pypsa_demand):
    print("Generating demand comparison plot")
    pypsa_loads_yearly_flat = pypsa_demand["Value"]
    ember_demand_yearly_flat = ember_demand["Value"]
    
    common_countries = pypsa_loads_yearly_flat.index.intersection(ember_demand_yearly_flat.index)
    pypsa_vals = pypsa_loads_yearly_flat.loc[common_countries]
    ember_vals = ember_demand_yearly_flat.loc[common_countries]
    
    df_compare = pd.DataFrame({
        "PyPSA-Eur": pypsa_vals,
        "Ember": ember_vals
    }, index=common_countries)
    
    n_countries = len(df_compare)
    mid = n_countries // 2
    countries1 = df_compare.index[:mid]
    countries2 = df_compare.index[mid:]
    
    fig, axes = plt.subplots(1, 2, figsize=(9, 6.5), sharey=False)
    
    width = 0.35
    for ax, countries, title in zip(axes, [countries1, countries2], ["Countries 1", "Countries 2"]):
        y = np.arange(len(countries))
        ax.barh(y - width/2, df_compare.loc[countries, "PyPSA-Eur"], height=width, label='PyPSA-Eur', color='#192238')
        ax.barh(y + width/2, df_compare.loc[countries, "Ember"], height=width, label='Ember', color='#13ce74')
        ax.set_yticks(y)
        ax.set_yticklabels(countries)
        ax.set_xlabel('Yearly Electricity Demand (TWh)')
        ax.legend()
    
    fig.suptitle(f"Electricity Demand Comparison {year}", weight="bold")
    plt.tight_layout()
    
    output_path = snakemake.output.demand_plot if 'snakemake' in globals() else f"results/validation_{year}/plots/demand_plot.png"
    print(f"Saving demand comparison plot: {output_path}")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

# Main execution
if __name__ == "__main__":
    ember_capacity = process_ember_capacity()
    ember_demand = process_ember_demand()
    pypsa_capacity = process_pypsa_capacity()
    pypsa_demand = process_pypsa_demand()
    
    plot_network()
    plot_total_capacity_comparison(ember_capacity, pypsa_capacity)
    plot_capacity_comparison_horizontal(
        countries=["DE", "NL", "IT", "PL", "CZ", "GR"],
        ember_capacity=ember_capacity,
        pypsa_capacity=pypsa_capacity
    )
    plot_demand_comparison(ember_demand, pypsa_demand)