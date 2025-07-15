# Plotting Tools for PumpKin Analysis

This directory contains plotting tools for analyzing PumpKin chemical reaction network outputs.

## plotting_main.py - Unified Command-Line Interface

The `plotting_main.py` script provides a comprehensive command-line interface to run all plotting functions from the Outputs folder. It supports all the individual plotting modules and allows you to configure common parameters like number of cells, data directory, and output directory.

### Features

- **Species-specific plots**: Pathway analysis, net rates, heatmaps, and temperature correlations
- **Pathway rates plots**: Aggregate pathway analysis across all cells
- **Reaction rates plots**: Individual reaction rate analysis
- **Deleted pathways plots**: Analysis of deleted pathway contributions
- **Spaghetti plots**: Species density and mixing ratio analysis

### Usage

```bash
# Run all plotting functions with default settings
python plotting_main.py --all

# Run only species plots for specific species
python plotting_main.py --species-plots --species H H2 O OH --data-dir ./

# Run pathway rates plots with custom settings
python plotting_main.py --pathway-rates --num-cells 100 --top-n-pathways 15

# Run reaction rates plots with custom file paths
python plotting_main.py --reaction-rates --reactions-file ../qt_reactions_list.txt --rates-file ../qt_rates.txt

# Run deleted pathways analysis
python plotting_main.py --deleted-pathways --num-cells 233

# Run spaghetti plots
python plotting_main.py --spaghetti --species-file ../qt_species_list.txt --densities-file ../qt_densities.txt
```

### Common Parameters

- `--data-dir`: Directory containing data files (default: ./)
- `--output-dir`: Output directory for plots (default: ./Plots/)
- `--num-cells`: Number of cells to process (default: 233)
- `--file-pattern`: Pattern for cell files (default: pumpkin_output_cell_*.txt)
- `--species`: Species to plot (default: H H2 O OH)
- `--top-n`: Number of top pathways to show (default: 5)

### Dependencies

The tool requires the following Python packages:
- pandas
- numpy
- matplotlib
- seaborn
- re (built-in)
- os (built-in)
- traceback (built-in)

If dependencies are missing, the tool will gracefully skip unavailable modules and continue with available ones.

## Individual Plotting Modules

### species_specific_plots.py
Generates species-specific pathway analysis plots:
- Top production/consumption pathways
- Net rates per cell
- Pathway heatmaps
- Stacked production fractions
- Net percentage vs temperature
- Pathway-species matrix

### pathway_rates_plots.py
Analyzes pathway rates across all cells:
- Total pathway rates vs cell number
- Supports both single-line and multi-line pathway parsing
- Log-log plotting for rate visualization

### rates_of_reactions_plots.py
Plots individual reaction rates:
- Reaction rates vs cell number
- Supports thresholding for better visualization
- Both overview and zoomed-in plots

### deleted_pathways_plots.py
Analyzes contributions from deleted pathways:
- Heatmaps of production/consumption percentages
- Identifies species with high deleted pathway contributions
- Flags entries above threshold (>10%)

### spagehtti_plots.py
Plots species densities and mixing ratios:
- Number density vs cell number
- Mixing ratio vs cell number
- Log-log scale visualization

## Output Files

All plots are saved as PNG files in the specified output directory:
- `{species}_top_pathways.png`: Top pathways for each species
- `{species}_net_rate.png`: Net rates for each species
- `{species}_heatmap.png`: Pathway heatmaps for each species
- `{species}_stacked_fraction.png`: Stacked production fractions
- `{species}_net_pct_vs_temp.png`: Net percentage vs temperature
- `pathway_species_matrix.png`: Pathway-species occurrence matrix
- `TotalPathwayRates.png`: Total pathway rates across cells
- `reaction_rates.png`: Individual reaction rates
- `deleted_contributions.png`: Deleted pathway contributions heatmap
- `species_density_and_mixing_ratios.png`: Species densities and mixing ratios

## Getting Started

1. Ensure you have PumpKin output files in your data directory
2. Install required dependencies (pandas, numpy, matplotlib, seaborn)
3. Run the plotting tool with desired options
4. Check the output directory for generated plots

Example for a complete analysis:
```bash
python plotting_main.py --all --data-dir ./ --output-dir ./Plots/ --num-cells 233 --species H H2 O OH O2 H2O
```

This will generate all available plots for the specified species and save them in the Plots/ directory.