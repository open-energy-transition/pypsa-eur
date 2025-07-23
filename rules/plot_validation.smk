# SPDX-FileCopyrightText: Open Energy Transition gGmbH, Ember, and contributors to the Ember Flexibility Study
#
# SPDX-License-Identifier: MIT

# Snakemake rule to generate validation graphs for power generation and flows
rule plot_validation:
    input:
        network="results/validation_{year}/networks/base_s_{clusters}_elec_{opts}.nc",
        ember_monthly="validation/ember_data/europe_monthly_full_release_long_format.csv",
        ember_yearly="validation/ember_data/yearly_full_release_long_format.csv",
        power_flows="validation/entsoe_data/physical_energy_power_flows_2023.csv"
    output:
        donut_subplots="results/validation_{year}/plots/base_s_{clusters}_elec_{opts}_donut_subplots.png",
        donut_comparison="results/validation_{year}/plots/base_s_{clusters}_elec_{opts}_donut_comparison.png",
        ember_comparison="results/validation_{year}/plots/base_s_{clusters}_elec_{opts}_ember_comparison_de.png",
        power_flows_list="results/validation_{year}/plots/base_s_{clusters}_elec_{opts}_power_flows_list.txt"
    params:
        script="",
    script:
        "../scripts/generation_and_flows.py"


rule validate_ember_networks:
    input:
        expand(
            RESULTS + "plots/base_s_{clusters}_elec_{opts}_donut_subplots.png",
            year=2023,
            **config["scenario"],
            run=config["run"]["name"],
        ),


rule plot_capacity_demand:
    input:
        network="results/validation_{year}/networks/base_s_{clusters}_elec_{opts}.nc",
        ember_capacity="validation/ember_data/yearly_full_release_long_format.csv",
        ember_demand="validation/ember_data/europe_monthly_full_release_long_format.csv",
        regions_onshore=f"resources/validation_{{year}}/country_shapes.geojson"
    output:
        network_plot="results/validation_{year}/plots/network_plot_base_s_{clusters}_elec_{opts}.png",
        total_capacity_plot="results/validation_{year}/plots/total_capacity_plot_base_s_{clusters}_elec_{opts}.png",
        country_capacity_plot="results/validation_{year}/plots/country_capacity_plot_base_s_{clusters}_elec_{opts}.png",
        demand_plot="results/validation_{year}/plots/base_s_{clusters}_elec_{opts}_demand_plot.png"
    script:
        "../scripts/capacities_and_demand.py"

rule cross_country_capacity_comparison:
    input:
        network = "results/validation_2023/networks/base_s_39_elec_.nc",
        ember_csv = "validation/ember_data/REF_NTC.csv"
    output:
        comparison_csv = "results/validation_2023/plots/focus_countries_comparison.csv",
        bar_plot = "results/validation_2023/plots/focus_countries_comparison.png"
    script:
        "../scripts/cross_country_capacities.py"





rule plot_all_capacity_demand:
    input:
        expand(
            RESULTS + "plots/network_plot_base_s_{clusters}_elec_{opts}.png",
            year=2023,
            **config["scenario"],
            run=config["run"]["name"],
        ),
