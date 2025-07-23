import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys


sns.set_style("whitegrid", {"grid.linestyle": "--", "grid.alpha": 0.7})

# Load PyPSA-Eur network from sector_demand rule input
if len(sys.argv) > 1:
    network_path = sys.argv[1]
else:
    network_path = "path/to/your/pypsa-eur/network.nc"  # Placeholder: replace with actual path if not using Snakemake

n = pypsa.Network(network_path)

# Output path for the plot from sector_demand rule output
if len(sys.argv) > 2:
    output_path = sys.argv[2]
else:
    output_path = "sector_demand.png"


# Function to extract the carrier name
def get_carrier(col):
    parts = col.split()
    if len(parts) == 2 and parts[1].isdigit():
        return "electricity"
    elif parts[1].isdigit():
        return ' '.join(parts[2:])
    else:
        return ' '.join(parts[1:])

# Map column names to carriers
carrier_map = {col: get_carrier(col) for col in n.loads_t.p_set.columns}

# Identify available carriers
unique_carriers = set(carrier_map.values())
print("Available carriers:")
for carrier in sorted(unique_carriers):
    print(f"- {carrier}")

# Calculate total demand per carrier in MWh
total_demand_per_carrier_MWh = n.loads_t.p_set.groupby(carrier_map, axis=1).sum().sum(axis=0)

# Convert to TWh
total_demand_per_carrier_TWh = total_demand_per_carrier_MWh / 1_000_000

# Print demand results for all carriers
print("\nTotal demand per carrier (TWh):")
for carrier, demand in total_demand_per_carrier_TWh.items():
    print(f"- {carrier}: {demand:.6f} TWh")

# Aggregate for comparison
pypsa_demand_by_carrier = {
    'Electricity': total_demand_per_carrier_TWh.get('electricity', 0),
    'Heat Demand': total_demand_per_carrier_TWh.get('rural heat', 0) + 
                   total_demand_per_carrier_TWh.get('urban central heat', 0) + 
                   total_demand_per_carrier_TWh.get('urban decentral heat', 0),
    'Land Transport EV': total_demand_per_carrier_TWh.get('land transport EV', 0),
    'Land Transport oil': total_demand_per_carrier_TWh.get('land transport oil', 0)
}

# Publicly available data (TWh)
public_data = {
    'Electricity': 2790,
    'Heat Demand': 2780, # From https://www.sciencedirect.com/science/article/pii/S0306261925009663
    'Land Transport EV': 180, # From https://www.iea.org/reports/global-ev-outlook-2025/outlook-for-energy-demand
    'Land Transport oil': 2749  # From https://www.iea.org/reports/global-ev-outlook-2025/outlook-for-energy-demand
}

# Create DataFrames for comparison
pypsa_df = pd.DataFrame({
    'Carrier': pypsa_demand_by_carrier.keys(),
    'PyPSA-Eur (TWh)': pypsa_demand_by_carrier.values()
})
public_df = pd.DataFrame({
    'Carrier': public_data.keys(),
    'Public Data (TWh)': public_data.values()
})

# Merge DataFrames, filling missing values with 0
comparison_df = pd.merge(pypsa_df, public_df, on='Carrier', how='outer').fillna(0)

# Calculate total demand for each source
total_pypsa = comparison_df['PyPSA-Eur (TWh)'].sum()
total_public = comparison_df['Public Data (TWh)'].sum()

# Add total rows
comparison_df = pd.concat([
    comparison_df,
    pd.DataFrame({
        'Carrier': ['Total PyPSA-Eur', 'Total Public'],
        'PyPSA-Eur (TWh)': [total_pypsa, 0],
        'Public Data (TWh)': [0, total_public]
    })
], ignore_index=True)

# Create the plot
plt.figure(figsize=(12, 6), dpi=100)

# Use a professional color palette
colors = sns.color_palette("Set2", n_colors=2)

# Plot with improved styling
x = range(len(comparison_df) - 2)  # Exclude total rows for individual carriers
width = 0.35

plt.bar([p - width/2 for p in x], comparison_df['PyPSA-Eur (TWh)'][:-2], width, label='PyPSA-Eur', color=colors[0], edgecolor='black')
plt.bar([p + width/2 for p in x], comparison_df['Public Data (TWh)'][:-2], width, label='Public Data', color=colors[1], edgecolor='black')

# Plot total bars
plt.bar(len(comparison_df) - 2, total_pypsa, width, color=colors[0], alpha=0.6, label='PyPSA-Eur Total', edgecolor='black')
plt.bar(len(comparison_df) - 1, total_public, width, color=colors[1], alpha=0.6, label='Public Data Total', edgecolor='black')

# Customize title and labels
plt.title("Comparison of Demand by Carrier (TWh)", fontsize=16, weight="bold", pad=15)
plt.xlabel("Carrier", fontsize=12, labelpad=10)
plt.ylabel("Demand (TWh)", fontsize=12, labelpad=10)

# Customize x-axis
plt.xticks([p for p in x] + [len(comparison_df) - 2, len(comparison_df) - 1],
           comparison_df['Carrier'][:-2].tolist() + ['Total PyPSA-Eur', 'Total Public'],
           rotation=45, ha='right', fontsize=10)

# Customize y-axis
plt.yticks(fontsize=10)
plt.grid(axis='y', linestyle='--', alpha=0.7, color="gray")

# Remove top and right spines
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_color('gray')
plt.gca().spines['bottom'].set_color('gray')

# Add value labels on top of bars
for i, (pypsa_val, public_val) in enumerate(zip(comparison_df['PyPSA-Eur (TWh)'][:-2], comparison_df['Public Data (TWh)'][:-2])):
    plt.text(i - width/2, pypsa_val + 10, f"{pypsa_val:.2f}", ha='center', va='bottom', fontsize=9, color='black')
    plt.text(i + width/2, public_val + 10, f"{public_val:.2f}", ha='center', va='bottom', fontsize=9, color='black')

plt.text(len(comparison_df) - 2, total_pypsa + 10, f"{total_pypsa:.2f}", ha='center', va='bottom', fontsize=9, color='black')
plt.text(len(comparison_df) - 1, total_public + 10, f"{total_public:.2f}", ha='center', va='bottom', fontsize=9, color='black')

# Add legend
plt.legend()

# Tight layout for better spacing
plt.tight_layout()

# Save the plot to the output path
plt.savefig(output_path)

# Show the plot
plt.show()

# Print comparison results
print("\nTotal demand per carrier (TWh):")
for carrier, pypsa_demand, public_demand in zip(comparison_df['Carrier'][:-2], comparison_df['PyPSA-Eur (TWh)'][:-2], comparison_df['Public Data (TWh)'][:-2]):
    print(f"- {carrier}: PyPSA-Eur = {pypsa_demand:.2f} TWh, Public Data = {public_demand:.2f} TWh")
print(f"- Total PyPSA-Eur: {total_pypsa:.2f} TWh")
print(f"- Total Public Data: {total_public:.2f} TWh")

print("\nDiscrepancy Note: Differences between PyPSA-Eur and public data may arise from broader regional coverage in PyPSA-Eur, varying data sources, or different time periods/assumptions.")