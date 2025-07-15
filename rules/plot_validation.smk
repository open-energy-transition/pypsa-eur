# SPDX-FileCopyrightText: Open Energy Transition gGmbH, Ember, and contributors to the Ember Flexibility Study
#
# SPDX-License-Identifier: MIT

# Snakemake rule to generate validation graphs for power generation and flows
rule plot_validation:
    input:
        network="results/validation_{year}/networks/base_s_39_elec_.nc",
        ember_monthly="validation/ember_data/europe_monthly_full_release_long_format.csv",
        ember_yearly="validation/ember_data/yearly_full_release_long_format.csv",
        power_flows="validation/entsoe_data/physical_energy_power_flows_2023.csv"
    output:
        donut_subplots="results/validation_{year}/plots/donut_subplots.png",
        donut_comparison="results/validation_{year}/plots/donut_comparison.png",
        ember_comparison="results/validation_{year}/plots/ember_comparison_de.png",
        power_flows_list="results/validation_{year}/plots/power_flows_list.txt"
    params:
        script="scripts/generation_and_flows.py",
        plot_flag=config.get("Plot_Validation_generation_and_flows_graphs", True)
    run:
        from pathlib import Path
        import subprocess
        
        output_dir = f"results/validation_{wildcards.year}/plots"
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        if params.plot_flag:
            print(f"Running script: {params.script}")
            result = subprocess.run(["python", params.script], capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Script error: {result.stderr}")
                raise RuntimeError("Script execution failed")
            print(f"Script output: {result.stdout}")
        else:
            print("Plotting disabled via config. Creating placeholder outputs.")
            Path(output.donut_subplots).touch()
            Path(output.donut_comparison).touch()
            Path(output.ember_comparison).touch()
            Path(output.power_flows_list).touch()


rule plot_capacity_demand:
    input:
        network="results/validation_{year}/networks/base_s_39_elec_.nc",
        ember_capacity="validation/ember_data/yearly_full_release_long_format.csv",
        ember_demand="validation/ember_data/europe_monthly_full_release_long_format.csv",
        regions_onshore=f"resources/validation_{{year}}/country_shapes.geojson"
    output:
        network_plot="results/validation_{year}/plots/network_plot.png",
        total_capacity_plot="results/validation_{year}/plots/total_capacity_plot.png",
        country_capacity_plot="results/validation_{year}/plots/country_capacity_plot.png",
        demand_plot="results/validation_{year}/plots/demand_plot.png"
    params:
        script="scripts/capacities_and_demand.py",
        plot_flag=config.get("Plot_Capacity_Demand", True)
    run:
        from pathlib import Path
        import subprocess
        
        output_dir = f"results/validation_{wildcards.year}/plots"
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        if params.plot_flag:
            print(f"Running script: {params.script}")
            result = subprocess.run(["python", params.script], capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Script error: {result.stderr}")
                raise RuntimeError("Script execution failed")
            print(f"Script output: {result.stdout}")
        else:
            print("Plotting disabled via config. Creating placeholder outputs.")
            Path(output.network_plot).touch()
            Path(output.total_capacity_plot).touch()
            Path(output.country_capacity_plot).touch()
            Path(output.demand_plot).touch()