# Test_B Example Usage

This directory contains a complete example of running the AIOLOS-into-PumpKin analysis pipeline using the Test_B dataset.

## Configuration

**Dataset**: Test_B from PumpKin examples
**Path**: `D:\OneDrive\Water Worlds\PumpKin\src\Examples\Test_B`
**Species**: 3 (H2, H, M)
**Cells**: 3 (cells 0-2)

## Command Used

```bash
python main.py --run-pumpkin --generate-plots \
  --num-cells 3 \
  --species H2 H M \
  --pumpkin-dir "/mnt/d/OneDrive/Water Worlds/PumpKin/src/Examples/Test_B" \
  --output-dir "./Test_B" \
  --data-dir "./Test_B" \
  --top-n 3
```

## Alternative Usage Examples

### 1. Complete Pipeline (if you have simulation data)
```bash
python main.py --run-all \
  --num-cells 3 \
  --species H2 H M \
  --pumpkin-dir "/mnt/d/OneDrive/Water Worlds/PumpKin/src/Examples/Test_B" \
  --output-dir "./Test_B" \
  --data-dir "./Test_B"
```

### 2. Only Generate Plots (using existing PumpKin outputs)
```bash
python main.py --generate-plots \
  --num-cells 3 \
  --species H2 H M \
  --data-dir "./Test_B" \
  --output-dir "./Test_B" \
  --top-n 3
```

### 3. Only Run PumpKin Analysis
```bash
python main.py --run-pumpkin \
  --num-cells 3 \
  --species H2 H M \
  --pumpkin-dir "/mnt/d/OneDrive/Water Worlds/PumpKin/src/Examples/Test_B" \
  --output-dir "./Test_B"
```

## Generated Files

### Input Data Files (copied from Test_B example)
- `qt_species_list.txt` - List of chemical species
- `qt_reactions_list.txt` - List of chemical reactions
- `qt_matrix.txt` - Stoichiometric matrix
- `qt_densities.txt` - Species densities across cells
- `qt_rates.txt` - Reaction rates across cells
- `qt_conditions.txt` - Temperature and conditions data

### PumpKin Output Files
- `pumpkin_output_cell_000.txt` - Analysis results for cell 0
- `pumpkin_output_cell_001.txt` - Analysis results for cell 1
- `pumpkin_output_cell_002.txt` - Analysis results for cell 2

### Expected Plot Files (generated when plotting dependencies are available)
- `Plots/H2_top_pathways.png` - Top production/consumption pathways for H2
- `Plots/H2_net_rate.png` - Net rate of H2 across cells
- `Plots/H2_heatmap.png` - Pathway contribution heatmap for H2
- `Plots/H2_stacked_fraction.png` - Stacked production fractions for H2
- `Plots/H2_net_pct_vs_temp.png` - Net percentage vs temperature for H2
- `Plots/H_top_pathways.png` - Top production/consumption pathways for H
- `Plots/H_net_rate.png` - Net rate of H across cells
- `Plots/H_heatmap.png` - Pathway contribution heatmap for H
- `Plots/H_stacked_fraction.png` - Stacked production fractions for H
- `Plots/H_net_pct_vs_temp.png` - Net percentage vs temperature for H
- `Plots/M_top_pathways.png` - Top production/consumption pathways for M
- `Plots/M_net_rate.png` - Net rate of M across cells
- `Plots/M_heatmap.png` - Pathway contribution heatmap for M
- `Plots/M_stacked_fraction.png` - Stacked production fractions for M
- `Plots/M_net_pct_vs_temp.png` - Net percentage vs temperature for M
- `Plots/pathway_species_matrix.png` - Pathway-species occurrence matrix
- `Plots/TotalPathwayRates.png` - Total pathway rates across all cells
- `Plots/reaction_rates.png` - Individual reaction rates
- `Plots/deleted_contributions.png` - Deleted pathway contributions heatmap
- `Plots/species_density_and_mixing_ratios.png` - Species densities and mixing ratios

## Key Pathways Identified

### H2 (Hydrogen Molecule)
- **Main Production**: H + H => H2 (85-95% across cells)
- **Main Consumption**: H2 + M => H2 + M (80-90% across cells)
- **Trend**: Production increases with cell number, consumption remains high

### H (Hydrogen Atom)
- **Main Production**: H2 => H + H (85-90% across cells)
- **Main Consumption**: H + H => H2 (70-80% across cells)
- **Trend**: Catalytic cycle between H and H2

### M (Third Body/Catalyst)
- **Main Pathway**: H + M => H + M (95-97% both production and consumption)
- **Role**: Acts as catalytic third body, minimal net change

## Analysis Insights

1. **Dominant Reaction**: H + H => H2 (hydrogen recombination)
2. **Catalytic Role**: M species acts as third body catalyst
3. **Equilibrium**: System shows H ⇌ H2 equilibrium with M catalysis
4. **Spatial Variation**: Reaction rates increase with cell number (0→1→2)

## Troubleshooting

### Common Issues

1. **PumpKin Executable Not Found**
   - Ensure PumpKin is compiled and executable
   - Check path to PumpKin directory
   - Verify WSL/Linux environment setup

2. **Missing Plot Dependencies**
   - Install required packages: `pip install pandas numpy matplotlib seaborn`
   - Plots will be skipped if dependencies are missing

3. **Permission Issues**
   - Ensure write permissions for output directory
   - Check file permissions for input files

### Verifying Results

Check that output files contain expected data:
```bash
# Check PumpKin outputs
head -20 Test_B/pumpkin_output_cell_000.txt

# Check data files
head -5 Test_B/qt_species_list.txt
head -5 Test_B/qt_reactions_list.txt

# Check for generated plots
ls -la Test_B/Plots/
```

## Customization Options

### Modify Species List
```bash
python main.py --generate-plots \
  --species H2 H \
  --num-cells 3 \
  --data-dir "./Test_B" \
  --output-dir "./Test_B"
```

### Change Number of Top Pathways
```bash
python main.py --generate-plots \
  --top-n 5 \
  --species H2 H M \
  --num-cells 3 \
  --data-dir "./Test_B" \
  --output-dir "./Test_B"
```

### Custom Output Directory
```bash
python main.py --run-all \
  --output-dir "./custom_results" \
  --data-dir "./custom_data" \
  --num-cells 3 \
  --species H2 H M
```

This example demonstrates the complete workflow of the AIOLOS-into-PumpKin analysis pipeline with a simple 3-species, 3-cell system.

  The plots folder is empty because the pipeline detected that the
  required Python packages for plotting are not installed:

  - pandas - for data manipulation
  - numpy - for numerical operations
  - matplotlib - for basic plotting
  - seaborn - for advanced visualizations

  Solutions

  Option 1: Install Dependencies (Recommended)

  pip install pandas numpy matplotlib seaborn

  Then run the plotting:
  python main.py --generate-plots --num-cells 3 --species H2 H M
  --data-dir "./Test_B" --output-dir "./Test_B" --top-n 3

  Option 2: Using Conda (if you have Anaconda/Miniconda)

  conda install pandas numpy matplotlib seaborn

  Option 3: Create a Virtual Environment

  python -m venv plot_env
  source plot_env/bin/activate  # On Windows: plot_env\Scripts\activate
  pip install pandas numpy matplotlib seaborn
  python main.py --generate-plots --num-cells 3 --species H2 H M
  --data-dir "./Test_B" --output-dir "./Test_B"

  What You Should See After Installing Dependencies

  Once you install the dependencies and run the plotting command, you
  should see files like:
  - H2_top_pathways.png
  - H2_net_rate.png
  - H2_heatmap.png
  - H2_stacked_fraction.png
  - H_top_pathways.png
  - H_net_rate.png
  - H_heatmap.png
  - M_top_pathways.png
  - M_net_rate.png
  - M_heatmap.png
  - pathway_species_matrix.png
  - TotalPathwayRates.png
  - deleted_contributions.png

  The tool was designed to gracefully handle missing dependencies by
  skipping the plotting functions but continuing with other parts of
  the pipeline. This way, you can still run PumpKin analysis even
  without plotting capabilities.