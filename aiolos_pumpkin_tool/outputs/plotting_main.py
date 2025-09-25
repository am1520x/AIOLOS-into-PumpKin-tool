#!/usr/bin/env python3
"""
plotting_main.py - Comprehensive command-line tool for running all plotting functions

This script provides a unified interface to run all plotting functions in the \nOutputs folder:
- Species-specific plots (pathway analysis, rates, heatmaps)
- Pathway rates plots (aggregate pathway analysis)
- Reaction rates plots (individual reaction analysis)
- Deleted pathways plots (analysis of deleted pathway contributions)
- Spaghetti plots (species density and mixing ratios)
"""

import argparse
import os
import sys
import traceback
from pathlib import Path

# Import plotting modules with error handling
try:
    from species_specific_plots import (
        load_all_cells_as_dict,
        plot_top_species_pathways_by_rate,
        plot_species_net_rate,
        plot_pathway_heatmap,
        plot_stacked_production_fraction,
        plot_species_net_percentage_vs_temperature,
        compute_pathway_species_matrix,
        plot_pathway_species_heatmap,
    )

    SPECIES_PLOTS_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Species plots module not available: {e}")
    SPECIES_PLOTS_AVAILABLE = False

try:
    from pathway_rates_plots import (
        parse_multi_line_pathway_table,
        parse_single_line_pathway_table,
        main as pathway_rates_main,
    )

    PATHWAY_RATES_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Pathway rates plots module not available: {e}")
    PATHWAY_RATES_AVAILABLE = False

try:
    from rates_of_reactions_plots import (
        plot_reaction_rates,
        main as reaction_rates_main,
    )

    REACTION_RATES_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Reaction rates plots module not available: {e}")
    REACTION_RATES_AVAILABLE = False

try:
    from deleted_pathways_plots import (
        parse_file as parse_deleted_file,
        plot_deleted_contributions,
        main as deleted_pathways_main,
    )

    DELETED_PATHWAYS_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Deleted pathways plots module not available: {e}")
    DELETED_PATHWAYS_AVAILABLE = False

try:
    from spagehtti_plots import plot_densities_and_mixing_ratios, main as spaghetti_main

    SPAGHETTI_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Spaghetti plots module not available: {e}")
    SPAGHETTI_AVAILABLE = False


def run_species_plots(args):
    """Run all species-specific plotting functions."""
    if not SPECIES_PLOTS_AVAILABLE:
        print("Species plots module not available. Skipping...")
        return False

    print(f"Running species-specific plots for {len(args.species)} species...")

    # Load cell data
    cell_dict = load_all_cells_as_dict(args.data_dir, args.file_pattern)
    if not cell_dict:
        print(f"No cell data found in {args.data_dir} with pattern {args.file_pattern}")
        return False

    print(f"Loaded data for {len(cell_dict)} cells")

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Generate plots for each species
    for species in args.species:
        print(f"  Generating plots for species: {species}")

        # Top pathways plot
        if args.plot_top_pathways:
            plot_top_species_pathways_by_rate(
                cell_dict,
                species,
                top_n=args.top_n,
                save_path=os.path.join(args.output_dir, f"{species}_top_pathways.png"),
            )

        # Net rate plot
        if args.plot_net_rate:
            plot_species_net_rate(
                cell_dict,
                species,
                save_path=os.path.join(args.output_dir, f"{species}_net_rate.png"),
            )

        # Pathway heatmap
        if args.plot_heatmap:
            plot_pathway_heatmap(
                cell_dict,
                species,
                save_path=os.path.join(args.output_dir, f"{species}_heatmap.png"),
            )

        # Stacked production fraction
        if args.plot_stacked:
            plot_stacked_production_fraction(
                cell_dict,
                species,
                top_n=args.top_n,
                save_path=os.path.join(
                    args.output_dir, f"{species}_stacked_fraction.png"
                ),
            )

        # Net percentage vs temperature
        if args.plot_net_pct_temp and args.temp_file:
            plot_species_net_percentage_vs_temperature(
                cell_dict,
                args.temp_file,
                species,
                save_path=os.path.join(
                    args.output_dir, f"{species}_net_pct_vs_temp.png"
                ),
            )

    # Pathway-species matrix
    if args.plot_pathway_matrix:
        print("  Generating pathway-species matrix...")
        matrix = compute_pathway_species_matrix(cell_dict)
        plot_pathway_species_heatmap(
            matrix,
            save_path=os.path.join(args.output_dir, "pathway_species_matrix.png"),
        )

    print("Species plots completed successfully!")
    return True


def run_pathway_rates_plots(args):
    """Run pathway rates plotting functions."""
    if not PATHWAY_RATES_AVAILABLE:
        print("Pathway rates plots module not available. Skipping...")
        return False

    print("Running pathway rates plots...")

    try:
        pathway_rates_main(
            data_dir=args.data_dir,
            num_cells=args.num_cells,
            use_multi_line_parser=args.use_multi_line_parser,
            top_n_pathways=args.top_n_pathways,
            save_plot_path=os.path.join(args.output_dir, "TotalPathwayRates.png"),
        )
        print("Pathway rates plots completed successfully!")
        return True
    except Exception as e:
        print(f"Error in pathway rates plotting: {e}")
        traceback.print_exc()
        return False


def run_reaction_rates_plots(args):
    """Run reaction rates plotting functions."""
    if not REACTION_RATES_AVAILABLE:
        print("Reaction rates plots module not available. Skipping...")
        return False

    print("Running reaction rates plots...")

    try:
        reaction_rates_main(
            reactions_path=args.reactions_file,
            rates_path=args.rates_file,
            output_plot=os.path.join(args.output_dir, "reaction_rates.png"),
        )
        print("Reaction rates plots completed successfully!")
        return True
    except Exception as e:
        print(f"Error in reaction rates plotting: {e}")
        traceback.print_exc()
        return False


def run_deleted_pathways_plots(args):
    """Run deleted pathways plotting functions."""
    if not DELETED_PATHWAYS_AVAILABLE:
        print("Deleted pathways plots module not available. Skipping...")
        return False

    print("Running deleted pathways plots...")

    try:
        # Set global variables for deleted pathways module
        import deleted_pathways_plots

        deleted_pathways_plots.DATA_DIR = args.data_dir
        deleted_pathways_plots.N_CELLS = args.num_cells

        deleted_pathways_main()
        print("Deleted pathways plots completed successfully!")
        return True
    except Exception as e:
        print(f"Error in deleted pathways plotting: {e}")
        traceback.print_exc()
        return False


def run_spaghetti_plots(args):
    """Run spaghetti plots (species density and mixing ratios)."""
    if not SPAGHETTI_AVAILABLE:
        print("Spaghetti plots module not available. Skipping...")
        return False

    print("Running spaghetti plots...")

    try:
        spaghetti_main(
            species_file=args.species_file,
            densities_file=args.densities_file,
            output_plot=os.path.join(
                args.output_dir, "species_density_and_mixing_ratios.png"
            ),
        )
        print("Spaghetti plots completed successfully!")
        return True
    except Exception as e:
        print(f"Error in spaghetti plotting: {e}")
        traceback.print_exc()
        return False


def create_parser():
    """Create the argument parser with all options."""
    parser = argparse.ArgumentParser(
        description="Comprehensive plotting tool for PumpKin analysis results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all plots with default settings
  python plotting_main.py --all

  # Run only species plots for specific species
  python plotting_main.py --species-plots --species H H2 O OH --data-dir ./

  # Run pathway rates plots with custom settings
  python plotting_main.py --pathway-rates --num-cells 100 --top-n-pathways 15

  # Run reaction rates plots with custom file paths
  python plotting_main.py --reaction-rates \
      --reactions-file ../qt_reactions_list.txt --rates-file ../qt_rates.txt

  # Run deleted pathways analysis
  python plotting_main.py --deleted-pathways --num-cells 233

  # Run spaghetti plots
  python plotting_main.py --spaghetti \
      --species-file ../qt_species_list.txt --densities-file ../qt_densities.txt
        """,
    )

    # Main operation modes
    parser.add_argument("--all", action="store_true", help="Run all plotting functions")
    parser.add_argument(
        "--species-plots", action="store_true", help="Run species-specific plots"
    )
    parser.add_argument(
        "--pathway-rates", action="store_true", help="Run pathway rates plots"
    )
    parser.add_argument(
        "--reaction-rates", action="store_true", help="Run reaction rates plots"
    )
    parser.add_argument(
        "--deleted-pathways", action="store_true", help="Run deleted pathways plots"
    )
    parser.add_argument(
        "--spaghetti",
        action="store_true",
        help="Run spaghetti plots (densities/mixing ratios)",
    )

    # Common parameters
    parser.add_argument(
        "--data-dir",
        type=str,
        default="./",
        help="Directory containing data files (default: ./)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./Plots/",
        help="Output directory for plots (default: ./Plots/)",
    )
    parser.add_argument(
        "--num-cells",
        type=int,
        default=233,
        help="Number of cells to process (default: 233)",
    )
    parser.add_argument(
        "--file-pattern",
        type=str,
        default="pumpkin_output_cell_*.txt",
        help="Pattern for cell files (default: pumpkin_output_cell_*.txt)",
    )

    # Species-specific parameters
    parser.add_argument(
        "--species",
        nargs="+",
        default=["H", "H2", "O", "OH"],
        help="Species to plot (default: H H2 O OH)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Number of top pathways to show (default: 5)",
    )
    parser.add_argument(
        "--temp-file",
        type=str,
        default="qt_conditions.txt",
        help="Temperature conditions file (default: qt_conditions.txt)",
    )

    # Species plot type controls
    parser.add_argument(
        "--plot-top-pathways",
        action="store_true",
        default=True,
        help="Plot top pathways (default: True)",
    )
    parser.add_argument(
        "--plot-net-rate",
        action="store_true",
        default=True,
        help="Plot net rates (default: True)",
    )
    parser.add_argument(
        "--plot-heatmap",
        action="store_true",
        default=True,
        help="Plot pathway heatmaps (default: True)",
    )
    parser.add_argument(
        "--plot-stacked",
        action="store_true",
        default=True,
        help="Plot stacked production fractions (default: True)",
    )
    parser.add_argument(
        "--plot-net-pct-temp",
        action="store_true",
        default=True,
        help="Plot net percentage vs temperature (default: True)",
    )
    parser.add_argument(
        "--plot-pathway-matrix",
        action="store_true",
        default=True,
        help="Plot pathway-species matrix (default: True)",
    )

    # Pathway rates parameters
    parser.add_argument(
        "--use-multi-line-parser",
        action="store_true",
        default=True,
        help="Use multi-line pathway parser (default: True)",
    )
    parser.add_argument(
        "--top-n-pathways",
        type=int,
        default=25,
        help="Number of top pathways for pathway rates plot (default: 25)",
    )

    # Reaction rates parameters
    parser.add_argument(
        "--reactions-file",
        type=str,
        default="qt_reactions_list.txt",
        help="Reactions list file (default: qt_reactions_list.txt)",
    )
    parser.add_argument(
        "--rates-file",
        type=str,
        default="qt_rates.txt",
        help="Rates file (default: qt_rates.txt)",
    )

    # Spaghetti plot parameters
    parser.add_argument(
        "--species-file",
        type=str,
        default="qt_species_list.txt",
        help="Species list file (default: qt_species_list.txt)",
    )
    parser.add_argument(
        "--densities-file",
        type=str,
        default="qt_densities.txt",
        help="Densities file (default: qt_densities.txt)",
    )

    return parser


def main():
    """Main function to run the plotting tool."""
    parser = create_parser()
    args = parser.parse_args()

    # Check if no specific mode is selected
    if not any(
        [
            args.all,
            args.species_plots,
            args.pathway_rates,
            args.reaction_rates,
            args.deleted_pathways,
            args.spaghetti,
        ]
    ):
        print(
            "No plotting mode selected. Use --help for options or "
            "--all to run everything."
        )
        return 1

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    success_count = 0
    total_count = 0

    # Run selected plotting functions
    if args.all or args.species_plots:
        total_count += 1
        if run_species_plots(args):
            success_count += 1

    if args.all or args.pathway_rates:
        total_count += 1
        if run_pathway_rates_plots(args):
            success_count += 1

    if args.all or args.reaction_rates:
        total_count += 1
        if run_reaction_rates_plots(args):
            success_count += 1

    if args.all or args.deleted_pathways:
        total_count += 1
        if run_deleted_pathways_plots(args):
            success_count += 1

    if args.all or args.spaghetti:
        total_count += 1
        if run_spaghetti_plots(args):
            success_count += 1

    # Summary
    print(f"\n=== Plotting Summary ===")
    print(f"Successfully completed: {success_count}/{total_count} plotting functions")
    print(f"Output directory: {args.output_dir}")

    if success_count == total_count:
        print("All plotting functions completed successfully!")
        return 0
    else:
        print(f"Some plotting functions failed. Check output above for details.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
