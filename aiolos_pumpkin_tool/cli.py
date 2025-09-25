#!/usr/bin/env python3
"""
Command-line interface for the AIOLOS-into-PumpKin-tool package.

This module provides entry points for the package that can be used from the command line
after installation.
"""

import sys
import argparse
from .main import AIOLOSPumpKinPipeline
from .utils import validate_dependencies


def main():
    """Main command-line entry point."""
    # Check dependencies first
    deps = validate_dependencies()
    missing_deps = [dep for dep, available in deps.items() if not available]
    
    if missing_deps:
        print("Warning: Missing dependencies:")
        for dep in missing_deps:
            print(f"  - {dep}")
        if "numpy" in missing_deps or "pandas" in missing_deps:
            print("Critical dependencies missing. Install with: pip install numpy pandas")
            return 1
    
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="AIOLOS into PumpKin Analysis Pipeline",
        epilog="For more information, see the documentation."
    )
    
    # Add main arguments
    parser.add_argument(
        "--run-all", action="store_true",
        help="Run the complete analysis pipeline"
    )
    parser.add_argument(
        "--process-data", action="store_true",
        help="Process AIOLOS data only"
    )
    parser.add_argument(
        "--run-pumpkin", action="store_true",
        help="Run PumpKin analysis only"
    )
    parser.add_argument(
        "--generate-plots", action="store_true",
        help="Generate plots only"
    )
    
    # Configuration arguments
    parser.add_argument(
        "--num-cells", type=int, default=100,
        help="Number of cells to process (default: 100)"
    )
    parser.add_argument(
        "--data-dir", default="./data/",
        help="Directory containing simulation data (default: ./data/)"
    )
    parser.add_argument(
        "--output-dir", default="./results/",
        help="Output directory for results (default: ./results/)"
    )
    parser.add_argument(
        "--pumpkin-dir", 
        help="Path to PumpKin directory (auto-detected if not specified)"
    )
    
    # Data files
    parser.add_argument(
        "--reaction-file", default="thermo.reac",
        help="Reaction file to process (default: thermo.reac)"
    )
    parser.add_argument(
        "--species", nargs="+", default=["H", "H2", "O", "OH"],
        help="Species to analyze (default: H H2 O OH)"
    )
    
    # Advanced options
    parser.add_argument(
        "--timestep", type=int, default=85,
        help="Timestep to process (default: 85)"
    )
    parser.add_argument(
        "--top-n", type=int, default=5,
        help="Number of top pathways to show in plots (default: 5)"
    )
    parser.add_argument(
        "--processing-type", default="timesteps",
        choices=["timesteps", "radial"],
        help="Type of processing to perform (default: timesteps)"
    )
    parser.add_argument(
        "--output-folder", default="Testing",
        help="Output folder name (default: Testing)"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate arguments
    if not any([args.run_all, args.process_data, args.run_pumpkin, args.generate_plots]):
        print("Error: Must specify at least one action (--run-all, --process-data, --run-pumpkin, or --generate-plots)")
        parser.print_help()
        return 1
    
    try:
        # Create and run pipeline
        pipeline = AIOLOSPumpKinPipeline(args)
        
        if args.run_all:
            success = pipeline.run_complete_pipeline()
        else:
            success = True
            if args.process_data:
                success &= pipeline.process_aiolos_data()
            if args.run_pumpkin:
                success &= pipeline.run_pumpkin_analysis()
            if args.generate_plots:
                success &= pipeline.generate_visualizations()
        
        if success:
            print("Pipeline completed successfully!")
            return 0
        else:
            print("Pipeline completed with warnings or errors.")
            return 1
            
    except Exception as e:
        print(f"Error running pipeline: {e}")
        import traceback
        traceback.print_exc()
        return 1


def conversion_main():
    """Entry point for conversion tool only."""
    from .conversion_tool.main_command import main as conv_main
    return conv_main()


def plotting_main():
    """Entry point for plotting tools only."""
    try:
        from .outputs.plotting_main import main as plot_main
        return plot_main()
    except ImportError:
        print("Plotting dependencies not available. Install with: pip install matplotlib seaborn")
        return 1


if __name__ == "__main__":
    sys.exit(main())