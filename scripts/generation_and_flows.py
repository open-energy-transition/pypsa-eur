import pandas as pd
import pypsa
import pycountry
import matplotlib.pyplot as plt
import numpy as np
import os
from itertools import combinations

# Configuration (updated for Snakemake for snakefile)
def get_config():
    return snakemake.config if 'snakemake' in globals() else {'year': 2023}

config = get_config()
year = config.get('year', 2023)
network_path = snakemake.input.network if 'snakemake' in globals() else "results/validation_2023/networks/base_s_39_elec_.nc"
ember_monthly_data_path = snakemake.input.ember_monthly if 'snakemake' in globals() else "validation/ember_data/europe_monthly_full_release_long_format.csv"
ember_yearly_data_path = snakemake.input.ember_yearly if 'snakemake' in globals() else "validation/ember_data/yearly_full_release_long_format.csv"
power_flows_data_path = snakemake.input.power_flows if 'snakemake' in globals() else "validation/entsoe_data/physical_energy_power_flows_2023.csv"

countries = ['AL', 'AT', 'BA', 'BE', 'BG', 'CH', 'CY', 'CZ', 'DE', 'DK', 'EE', 'ES',
             'FI', 'FR', 'GB', 'GR', 'HR', 'HU', 'IE', 'IT', 'LT', 'LU', 'LV', 'ME',
             'MK', 'MT', 'NL', 'NO', 'PL', 'PT', 'RO', 'RS', 'SE', 'SI', 'SK']

color_dict = {
    "Bioenergy": "#0f9970",
    "Gas": "#87530b",
    "Hard coal": "#0000006E",
    "Hydro": "#80b1d3",
    "Lignite": "#61340c",
    "Nuclear": "#f1dd02",
    "Offshore wind": "#0c0a48",
    "Onshore wind": "#2D9FAB",
    "Other fossil": "#545454",
    "Other renewables": "#9e1414",
    "Solar": "#ffaa00"
}

# Load data from ember
n = pypsa.Network(network_path)
ember_monthly = pd.read_csv(ember_monthly_data_path)

# Helper function to detect columns
def detect_column(columns, keywords):
    for keyword in keywords:
        for col in columns:
            if keyword.lower() in col.lower():
                return col
    return None

# Process Ember generation data (Yearly CSV)
def process_ember_generation_yearly():
    print(f"Processing yearly CSV: {ember_yearly_data_path}")
    if not os.path.exists(ember_yearly_data_path):
        print(f"Yearly CSV not found at {ember_yearly_data_path}. Falling back to monthly aggregation.")
        return None
    
    try:
        df = pd.read_csv(ember_yearly_data_path)
        print(f"Yearly CSV loaded. Shape: {df.shape}")
        print(f"Columns: {df.columns.tolist()}")
        
        iso_col = detect_column(df.columns, ['iso 3 code', 'iso3', 'country code'])
        variable_col = detect_column(df.columns, ['variable', 'fuel', 'technology'])
        value_col = detect_column(df.columns, ['value', 'generation', 'amount'])
        unit_col = detect_column(df.columns, ['unit'])
        subcategory_col = detect_column(df.columns, ['subcategory', 'category'])
        continent_col = detect_column(df.columns, ['continent', 'region'])
        year_col = detect_column(df.columns, ['year', 'date'])
        
        print(f"Detected columns: ISO={iso_col}, Variable={variable_col}, Value={value_col}, "
              f"Unit={unit_col}, Subcategory={subcategory_col}, Continent={continent_col}, Year={year_col}")
        
        def iso3_to_iso2(iso3):
            try:
                return pycountry.countries.get(alpha_3=iso3).alpha_2
            except:
                return None
        
        if iso_col:
            df["ISO"] = df[iso_col].apply(iso3_to_iso2)
        else:
            print("Warning: ISO column not found. Assuming 'ISO' column exists.")
            iso_col = "ISO"
        
        filters = []
        if continent_col:
            filters.append(df[continent_col] == "Europe")
        if iso_col:
            filters.append(df["ISO"].isin(countries))
        if unit_col:
            filters.append(df[unit_col] == "TWh")
        if subcategory_col:
            filters.append(df[subcategory_col] == "Fuel")
        if year_col:
            filters.append(df[year_col].astype(str).str.startswith(str(year)))
        
        if filters:
            df = df[np.logical_and.reduce(filters)]
        print(f"After filtering: Shape {df.shape}")
        
        required_cols = ["ISO", variable_col or "Variable", value_col or "Value"]
        available_cols = [col for col in required_cols if col in df.columns]
        df = df[available_cols]
        df = df.rename(columns={variable_col: "Variable", value_col: "Value"})
        
        df = df.groupby(["ISO", "Variable"], as_index=False)["Value"].sum()
        print(f"Processed yearly data sample:\n{df.head().to_string()}")
        return df.set_index(["ISO", "Variable"])
    except Exception as e:
        print(f"Error processing yearly CSV: {e}")
        return None

# Process Ember generation data (Monthly CSV) already downloaded by the workflow
def process_ember_generation_monthly():
    print(f"Processing monthly CSV: {ember_monthly_data_path}")
    df = ember_monthly[ember_monthly["Continent"] == "Europe"].copy()
    print(f"Monthly CSV loaded. Shape: {df.shape}")
    
    def iso3_to_iso2(iso3):
        try:
            return pycountry.countries.get(alpha_3=iso3).alpha_2
        except:
            return None
    
    df["ISO"] = df["ISO 3 code"].apply(iso3_to_iso2)
    df = df[df["ISO"].isin(countries)]
    df = df[df["Date"].str.startswith(str(year))]
    df = df[df["Unit"] == "TWh"]
    df = df[df["Subcategory"] == "Fuel"]
    df = df[["ISO", "Date", "Variable", "Value", "Unit"]]
    
    df_yearly = df.groupby(["ISO", "Variable"], as_index=False)["Value"].sum()
    df_yearly["Unit"] = "TWh"
    print(f"Processed monthly aggregated data sample:\n{df_yearly.head().to_string()}")
    return df_yearly.set_index(["ISO", "Variable"]).drop(["Unit"], axis=1)

# Wrapper function to choose processing method
def process_ember_generation(use_yearly=False):
    if use_yearly:
        result = process_ember_generation_yearly()
        if result is not None:
            return result
        print("Falling back to monthly aggregation due to yearly CSV issues.")
    return process_ember_generation_monthly()

# Process PyPSA-Eur model generation data
def process_pypsa_generation():
    pypsa_to_ember = {
        "biomass": "Bioenergy", "Bioenergy": "Bioenergy",
        "gas": "Gas", "Gas": "Gas", "CCGT": "Gas", "OCGT": "Gas",
        "coal": "Hard coal", "Hard coal": "Hard coal",
        "lignite": "Lignite", "Lignite": "Lignite",
        "hydro": "Hydro", "Hydro": "Hydro", "PHS": "Hydro", "ror": "Hydro",
        "Nuclear": "Nuclear", "nuclear": "Nuclear",
        "offwind-ac": "Offshore wind", "offwind-dc": "Offshore wind",
        "offwind-float": "Offshore wind", "Offshore wind": "Offshore wind",
        "onwind": "Onshore wind", "Onshore wind": "Onshore wind",
        "oil": "Other fossil", "Other fossil": "Other fossil",
        "geothermal": "Other renewables", "Other renewables": "Other renewables",
        "solar": "Solar", "solar-hsat": "Solar", "Solar": "Solar"
    }
    
    gen_meta = n.generators[["bus", "carrier"]].copy()
    gen_meta.loc[:, "country"] = gen_meta["bus"].str[:2]
    gen_energy = n.generators_t.p.sum(axis=0) / 1e6  # MWh to TWh
    gen_energy.index.name = "generator"
    gen_energy = gen_energy.reset_index().rename(columns={0: "Value"})
    gen_energy = gen_energy.merge(gen_meta, left_on="generator", right_index=True)
    gen_grouped = gen_energy.groupby(["country", "carrier"], as_index=False)["Value"].sum()
    gen_grouped["Ember_Variable"] = gen_grouped["carrier"].map(pypsa_to_ember).fillna(gen_grouped["carrier"])
    gen_renamed = gen_grouped.groupby(["country", "Ember_Variable"], as_index=False)["Value"].sum()
    gen_renamed = gen_renamed.rename(columns={"country": "ISO", "Ember_Variable": "Variable"}).round()
    
    sto_meta = n.storage_units[["bus", "carrier"]].copy()
    sto_meta.loc[:, "country"] = sto_meta["bus"].str[:2]
    sto_energy = n.storage_units_t.p.sum(axis=0) / 1e6  # MWh to TWh
    sto_energy.index.name = "storage_unit"
    sto_energy = sto_energy.reset_index().rename(columns={0: "Value"})
    sto_energy = sto_energy.merge(sto_meta, left_on="storage_unit", right_index=True)
    sto_grouped = sto_energy.groupby(["country", "carrier"], as_index=False)["Value"].sum()
    sto_grouped["Ember_Variable"] = sto_grouped["carrier"].map(pypsa_to_ember).fillna(sto_grouped["carrier"])
    sto_renamed = sto_grouped.groupby(["country", "Ember_Variable"], as_index=False)["Value"].sum()
    sto_renamed = sto_renamed.rename(columns={"country": "ISO", "Ember_Variable": "Variable"}).round()
    
    gen_and_sto = pd.concat([gen_renamed, sto_renamed], ignore_index=True)
    gen_and_sto = gen_and_sto.groupby(["ISO", "Variable"], as_index=False)["Value"].sum().round()
    print(f"PyPSA generation sample:\n{gen_and_sto.head().to_string()}")
    return gen_and_sto.set_index(["ISO", "Variable"])

# Compare Ember processing methods
def compare_ember_processing(df_yearly, df_monthly, country_iso="DE"):
    print(f"Comparing Ember processing for {country_iso}")
    
    def get_country_data(df, iso):
        if df is None:
            return pd.Series()
        if isinstance(df.index, pd.MultiIndex):
            data = df.unstack(level=1).fillna(0).loc[iso]
            data = data[data > 0]
        else:
            data = df.xs(iso, level="ISO")["Value"]
            data = data[data > 0]
        return data
    
    yearly_data = get_country_data(df_yearly, country_iso)
    monthly_data = get_country_data(df_monthly, country_iso)
    
    techs = list(set(yearly_data.index).union(monthly_data.index))
    yearly_values = [yearly_data.get(tech, 0) for tech in techs]
    monthly_values = [monthly_data.get(tech, 0) for tech in techs]
    
    print(f"Yearly data for {country_iso}:\n{yearly_data.to_string()}")
    print(f"Monthly aggregated data for {country_iso}:\n{monthly_data.to_string()}")
    
    plt.style.use('ggplot')
    fig, ax = plt.subplots(figsize=(10, 6))
    bar_width = 0.35
    x = np.arange(len(techs))
    ax.bar(x - bar_width/2, yearly_values, bar_width, label='Yearly CSV', color='skyblue')
    ax.bar(x + bar_width/2, monthly_values, bar_width, label='Monthly Aggregated', color='salmon')
    ax.set_xlabel('Technology')
    ax.set_ylabel('Generation (TWh)')
    ax.set_title(f'Ember Generation Comparison for {country_iso} (TWh, {year})')
    ax.set_xticks(x)
    ax.set_xticklabels(techs, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, axis='y')
    plt.tight_layout()
    
    output_path = snakemake.output.ember_comparison if 'snakemake' in globals() else f"results/validation_{year}/plots/ember_comparison_de.png"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

# Plotting functions
def plot_country_generation_mix_donut_subplots(ember_generation_yearly, country_isos, color_dict=None):
    n = len(country_isos)
    fig, axes = plt.subplots(3, 2, figsize=(7, 10))
    axes = axes.flatten()
    
    pivot_df = ember_generation_yearly.unstack(level=1).fillna(0)
    pivot_df.columns = pivot_df.columns.get_level_values(1)
    
    legend_handles = []
    legend_labels = []
    
    for idx, country_iso in enumerate(country_isos):
        ax = axes[idx]
        if country_iso not in pivot_df.index:
            ax.axis('off')
            ax.set_title(f"{country_iso} not found")
            continue
        
        data = pivot_df.loc[country_iso]
        data = data[data > 0]
        
        colors = [color_dict.get(tech, "#cccccc") for tech in data.index] if color_dict else plt.cm.Set2.colors[:len(data)]
        
        wedges, texts = ax.pie(
            data.values,
            labels=None,
            startangle=90,
            colors=colors,
            wedgeprops=dict(width=0.7),
            autopct=None
        )
        
        for i, wedge in enumerate(wedges):
            angle = (wedge.theta2 + wedge.theta1) / 2
            x = 0.7 * np.cos(np.deg2rad(angle))
            y = 0.7 * np.sin(np.deg2rad(angle))
            ax.text(x, y, f"{int(round(data.values[i]))}", ha='center', va='center', fontsize=10, color='white', fontweight='bold')
        
        centre_circle = plt.Circle((0, 0), 0.25, color='white', fc='white', linewidth=0)
        ax.add_artist(centre_circle)
        ax.text(0, 0, country_iso, ha='center', va='center', fontsize=18, fontweight='bold')
        
        if idx == 0:
            legend_handles = wedges
            legend_labels = data.index
    
    for j in range(len(country_isos), len(axes)):
        axes[j].axis('off')
    
    fig.legend(
        legend_handles, legend_labels,
        loc='upper center',
        bbox_to_anchor=(0.5, 1.05),
        ncol=4,
        fontsize=10,
        frameon=True
    )
    
    fig.suptitle(f"Yearly Electricity Generation by Technology\n(TWh, Ember {year})", fontsize=16, weight='bold', y=1.12)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    
    output_path = snakemake.output.donut_subplots if 'snakemake' in globals() else f"results/validation_{year}/plots/donut_subplots.png"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

def plot_country_generation_mix_donut_comparison(df1, df2, country_isos, color_dict=None, df1_label="Ember", df2_label="PyPSA"):
    n = len(country_isos)
    fig, axes = plt.subplots(n, 2, figsize=(6, 3 * n))
    plt.subplots_adjust(wspace=0.05)
    if n == 1:
        axes = np.array([axes])
    legend_handles = []
    legend_labels = []
    
    def pivot(df):
        if isinstance(df.index, pd.MultiIndex):
            p = df.unstack(level=1).fillna(0)
            p.columns = p.columns.get_level_values(1)
        else:
            p = df.copy()
        return p
    
    pivot1 = pivot(df1)
    pivot2 = pivot(df2)
    
    for idx, country_iso in enumerate(country_isos):
        for j, (pivot_df, label) in enumerate(zip([pivot1, pivot2], [df1_label, df2_label])):
            ax = axes[idx, j]
            if country_iso not in pivot_df.index:
                ax.axis('off')
                ax.set_title(f"{country_iso} not found")
                continue
            data = pivot_df.loc[country_iso]
            data = data[data > 0]
            colors = [color_dict.get(tech, "#cccccc") for tech in data.index] if color_dict else plt.cm.Set2.colors[:len(data)]
            wedges, _ = ax.pie(
                data.values,
                labels=None,
                startangle=90,
                colors=colors,
                wedgeprops=dict(width=0.7),
                autopct=None
            )
            for i, wedge in enumerate(wedges):
                angle = (wedge.theta2 + wedge.theta1) / 2
                x = 0.7 * np.cos(np.deg2rad(angle))
                y = 0.7 * np.sin(np.deg2rad(angle))
                ax.text(x, y, f"{int(round(data.values[i]))}", ha='center', va='center', fontsize=10, color='white', fontweight='bold')
            centre_circle = plt.Circle((0, 0), 0.25, color='white', fc='white', linewidth=0)
            ax.add_artist(centre_circle)
            total = int(round(data.sum()))
            ax.text(0, 0, f"{total}", ha='center', va='center', fontsize=14, fontweight='bold')
            ax.set_title(f"{label}\n{country_iso}" if idx == 0 else country_iso, fontsize=14, fontweight='bold')
            if idx == 0 and j == 0:
                legend_handles = wedges
                legend_labels = data.index
    
    for i in range(n, axes.shape[0]):
        for j in range(2):
            axes[i, j].axis('off')
    
    fig.legend(
        legend_handles, legend_labels,
        loc='upper center',
        bbox_to_anchor=(0.5, 1.02),
        ncol=4,
        fontsize=10,
        frameon=True
    )
    fig.suptitle("Yearly Electricity Generation [TWh]", fontsize=16, weight='bold', y=1.05)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    
    output_path = snakemake.output.donut_comparison if 'snakemake' in globals() else f"results/validation_{year}/plots/donut_comparison.png"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

# Power flows analysis
def analyze_power_flows():
    plt.style.use('ggplot')
    
    csv_file = power_flows_data_path
    
    for sep in [';', ',']:
        try:
            flows_2023 = pd.read_csv(csv_file, sep=sep)
            if len(flows_2023.columns) > 1:
                break
        except pd.errors.ParserError:
            continue
    else:
        raise ValueError("Unable to parse CSV file with ';' or ',' as separator. Check file format.")
    
    print("Columns in power flows CSV:", flows_2023.columns.tolist())
    print(f"Power flows data sample:\n{flows_2023.head().to_string()}")
    
    measure_time_col = detect_column(flows_2023.columns, ['measuretime'])
    from_area_col = detect_column(flows_2023.columns, ['fromareacode'])
    to_area_col = detect_column(flows_2023.columns, ['toareacode'])
    
    if not measure_time_col or not from_area_col or not to_area_col:
        raise KeyError("Required columns not found. Available columns: " + str(flows_2023.columns.tolist()))
    
    print(f"Detected 'MeasureTime' column: {measure_time_col}")
    print(f"Detected 'FromAreaCode' column: {from_area_col}")
    print(f"Detected 'ToAreaCode' column: {to_area_col}")
    
    flows_2023_monthly = flows_2023[flows_2023[measure_time_col] == 'Monthly Value']
    focus_countries = {'DE', 'NL', 'IT', 'PL', 'CZ', 'GR'}
    flows_2023_focus = flows_2023_monthly[
        flows_2023_monthly[from_area_col].isin(focus_countries) |
        flows_2023_monthly[to_area_col].isin(focus_countries)
    ]
    
    if flows_2023_focus.empty:
        raise ValueError("No data found for focus countries in ENTSO-E 2023 dataset")
    
    avg_flows_2023 = flows_2023_focus.groupby([from_area_col, to_area_col, 'Month'])['Value'].apply(
        lambda x: np.sum(np.abs(x))
    ).reset_index(name='Avg_Abs_Flow_2023_MW')
    avg_flows_2023['Connection'] = avg_flows_2023.apply(
        lambda row: f"{row[from_area_col]}-{row[to_area_col]}", axis=1
    )
    
    def get_country(bus):
        return bus[:2] if len(bus) >= 2 else None
    
    focus_lines = [
        line for line in n.lines.index
        if get_country(n.lines.loc[line, 'bus0']) in focus_countries or
           get_country(n.lines.loc[line, 'bus1']) in focus_countries
    ]
    
    connection_map = {
        line: "-".join(sorted([
            get_country(n.lines.loc[line, 'bus0']),
            get_country(n.lines.loc[line, 'bus1'])
        ])) for line in focus_lines
    }
    
    power_flows_focus = n.lines_t.p0[focus_lines].reset_index().melt(
        id_vars='snapshot', var_name='line', value_name='p0'
    )
    power_flows_focus['Connection'] = power_flows_focus['line'].map(connection_map)
    power_flows_focus['Month'] = power_flows_focus['snapshot'].dt.month
    power_flows_focus['abs_p0'] = np.abs(power_flows_focus['p0'])
    
    monthly_flows_pypsa = power_flows_focus.groupby(['Connection', 'Month'])['abs_p0'].mean().reset_index(
        name='Avg_Abs_Flow_PyPSA_MW'
    )
    
    comparison_df = pd.merge(
        avg_flows_2023[['Connection', 'Month', 'Avg_Abs_Flow_2023_MW']],
        monthly_flows_pypsa,
        on=['Connection', 'Month'],
        how='inner'
    )
    
    if comparison_df.empty:
        print("Warning: No common connections found between ENTSO-E 2023 and PyPSA-Eur data")
        return []
    
    valid_connections = comparison_df['Connection'].unique().tolist()
    print(f"Valid connections with data: {valid_connections}")
    
    power_flow_outputs = []
    for connection in valid_connections:
        output_path = snakemake.output.power_flows[valid_connections.index(connection)] if 'snakemake' in globals() else f"results/validation_{year}/plots/power_flow_{connection.replace('-', '_')}.png"
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        conn_data = comparison_df[comparison_df['Connection'] == connection]
        months = conn_data['Month']
        flow_2023 = conn_data['Avg_Abs_Flow_2023_MW']
        flow_pypsa = conn_data['Avg_Abs_Flow_PyPSA_MW']
        
        plt.figure(figsize=(12, 6))
        bar_width = 0.35
        x = np.arange(len(months))
        plt.bar(x - bar_width/2, flow_2023, bar_width, label='ENTSO-E 2023', color='skyblue')
        plt.bar(x + bar_width/2, flow_pypsa, bar_width, label='PyPSA-Eur 2023', color='salmon')
        plt.xlabel('Month')
        plt.ylabel('Average Absolute Power Flow (MW)')
        plt.title(f'Average Absolute Power Flow: {connection}')
        plt.xticks(x, months)
        plt.legend()
        plt.grid(True, axis='y')
        plt.tight_layout()
        
        plt.savefig(output_path, bbox_inches='tight', dpi=300)
        plt.close()
        power_flow_outputs.append(output_path)
    
    def get_net_export(n, country, snapshots):
        net_export = pd.Series(0.0, index=snapshots)
        for line in n.lines.index:
            bus0_country = get_country(n.lines.loc[line, 'bus0'])
            bus1_country = get_country(n.lines.loc[line, 'bus1'])
            if bus0_country == country and bus1_country != country:
                net_export += n.lines_t.p0[line]
            elif bus1_country == country and bus0_country != country:
                net_export += -n.lines_t.p0[line]
        return net_export
    
    for country in focus_countries:
        try:
            net_export = get_net_export(n, country, n.snapshots)
            avg_net_export = net_export.mean()
            print(f'Average net export for {country}: {avg_net_export:.2f} MW')
        except (KeyError, ValueError) as e:
            print(f"Error calculating net export for {country}: {e}")
    
    output_list_path = snakemake.output.power_flows_list if 'snakemake' in globals() else f"results/validation_{year}/plots/power_flows_list.txt"
    os.makedirs(os.path.dirname(output_list_path), exist_ok=True)
    with open(output_list_path, 'w') as f:
        for output_path in power_flow_outputs:
            f.write(output_path + '\n')
    
    return power_flow_outputs

# Main execution
if __name__ == "__main__":
    ember_generation_yearly = process_ember_generation(use_yearly=False)
    ember_generation_yearly_alt = process_ember_generation(use_yearly=True)
    
    compare_ember_processing(ember_generation_yearly_alt, ember_generation_yearly, country_iso="DE")
    
    pypsa_generation = process_pypsa_generation()
    
    plot_country_generation_mix_donut_subplots(
        ember_generation_yearly, 
        ["DE", "NL", "IT", "PL", "CZ", "GR"], 
        color_dict=color_dict
    )
    
    plot_country_generation_mix_donut_comparison(
        ember_generation_yearly, 
        pypsa_generation, 
        ["DE", "NL", "GR"], 
        color_dict=color_dict, 
        df1_label="Ember", 
        df2_label="PyPSA"
    )
    
    analyze_power_flows()