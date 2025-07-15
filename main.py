#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py - Complete AIOLOS-to-PumpKin Analysis Pipeline

This script provides a comprehensive command-line interface for the complete analysis pipeline:
1. Process AIOLOS reaction files and simulation data
2. Run PumpKin analysis across multiple cells
3. Generate comprehensive plots and visualizations

The tool integrates all components of the AIOLOS-into-PumpKin-tool workflow into a single,
easy-to-use command-line interface suitable for both new users and advanced analysis.

Usage Examples:
    # Basic run with default settings
    python main.py --run-all

    # Custom configuration
    python main.py --run-all --num-cells 100 --data-dir ./data/ --output-dir ./results/

    # Run only specific steps
    python main.py --process-data --run-pumpkin --num-cells 50
    python main.py --generate-plots --species H H2 O OH

    # Advanced configuration
    python main.py --run-all --pumpkin-dir "/custom/path" --species H H2 O OH H2O --top-n 10

Requirements:
    - Python 3.7+
    - Required packages: pandas, numpy, matplotlib, seaborn
    - PumpKin executable accessible via WSL or local installation
    - AIOLOS simulation data files

For detailed documentation and examples, see the README.md file.
"""

import argparse
import os
import sys
import subprocess
import re
import traceback
from pathlib import Path

# Add Conversion_Tool to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'Conversion_Tool'))

class AIOLOSPumpKinPipeline:
    """Main pipeline class for AIOLOS to PumpKin analysis."""
    
    def __init__(self, args):
        """Initialize pipeline with command-line arguments."""
        self.args = args
        self.setup_directories()
        self.validate_configuration()
    
    def setup_directories(self):
        """Create necessary directories."""
        os.makedirs(self.args.output_dir, exist_ok=True)
        os.makedirs(self.args.data_dir, exist_ok=True)
    
    def validate_configuration(self):
        """Validate configuration and warn about potential issues."""
        if self.args.run_pumpkin and not os.path.exists(self.args.pumpkin_dir):
            print(f"Warning: PumpKin directory {self.args.pumpkin_dir} does not exist")
        
        if self.args.num_cells <= 0:
            raise ValueError("Number of cells must be positive")
    
    def process_aiolos_data(self):
        """Process AIOLOS reaction files and simulation data."""
        print("=" * 60)
        print("STEP 1: Processing AIOLOS Data")
        print("=" * 60)
        
        try:
            # Import conversion tools
            from processing_aiolos_reac_file import parse_reaction_data, process_reaction_line, transform_species_set
            from making_densities_file import process_timesteps
            from making_rates_file import make_rates
            
            # Process reaction files
            print("Processing reaction files...")
            if not self._process_reaction_files():
                print("Warning: Reaction file processing failed")
                return False
            
            # Process simulation data
            print("Processing simulation data...")
            if not self._process_simulation_data():
                print("Error: Simulation data processing failed")
                return False
            
            print("SUCCESS: AIOLOS data processing completed successfully")
            return True
            
        except ImportError as e:
            print(f"Error: Missing conversion tools - {e}")
            return False
        except Exception as e:
            print(f"Error in AIOLOS data processing: {e}")
            traceback.print_exc()
            return False
    
    def _process_reaction_files(self):
        """Process AIOLOS reaction files."""
        try:
            from processing_aiolos_reac_file import parse_reaction_data, process_reaction_line, transform_species_set
            
            input_file = self.args.reaction_file
            if not os.path.exists(input_file):
                print(f"Warning: Reaction file {input_file} not found")
                return False
            
            with open(input_file, 'r') as f:
                reac_text = f.read()
            
            reaction_list, stoich_list, species_list = parse_reaction_data(reac_text)
            species_list = transform_species_set(species_list)
            
            # Generate output files
            with open(os.path.join(self.args.data_dir, self.args.species_file), "w") as f:
                for i, sp in enumerate(species_list, start=1):
                    f.write(f"{i} {sp}\n")
            
            with open(os.path.join(self.args.data_dir, self.args.reactions_file), "w") as f:
                for idx, reaction in enumerate(reaction_list, start=1):
                    modified = process_reaction_line(reaction)
                    f.write(f"{idx} {modified}\n")
            
            with open(os.path.join(self.args.data_dir, "qt_matrix.txt"), "w") as f:
                for sp in species_list:
                    row = [str(int(coeff)) if isinstance(coeff, (int, float)) and coeff.is_integer() 
                          else str(coeff) for coeff in (r.get(sp, 0) for r in stoich_list)]
                    f.write("\t".join(row) + "\n")
            
            print(f"Generated: {self.args.species_file}, {self.args.reactions_file}, qt_matrix.txt")
            return True
            
        except Exception as e:
            print(f"Error processing reaction files: {e}")
            return False
    
    def _process_simulation_data(self):
        """Process simulation timesteps and generate density/rate files."""
        try:
            from making_densities_file import process_timesteps
            from making_rates_file import make_rates
            
            # Default species configuration (can be made configurable)
            species = [
                ["H", 1, 1, r"$\rm H$"],
                ["H-", 1, -1, r"$\rm H^-$"],
                ["H2", 2, 0, r"$\rm H_2$"],
                ["H2O", 18, 0, r"$\rm H_2O$"],
                ["H2O2", 34, 0, r"$\rm H_2O_2$"],
                ["H2Op", 18, +1, r"$\rm H_2O^+$"],
                ["H2p", 2, +1, r"$\rm H_2^+$"],
                ["H3Op", 19, +1, r"$\rm H_3O^+$"],
                ["H3p", 3, +1, r"$\rm H_3^+$"],
                ["O", 16, 0, r"$\rm O$"],
                ["O-", 16, -1, r"$\rm O^-$"],
                ["O2", 32, 0, r"$\rm O_2$"],
                ["O2-", 32, -1, r"$\rm O_2^-$"],
                ["O2H", 33, 0, r"$\rm O_2H$"],
                ["O2Hp", 33, +1, r"$\rm O_2H^+$"],
                ["O2p", 32, +1, r"$\rm O_2^+$"],
                ["OH", 17, 0, r"$\rm OH$"],
                ["OH-", 17, -1, r"$\rm OH^-$"],
                ["OHp", 17, +1, r"$\rm OH^+$"],
                ["Op", 16, +1, r"$\rm O^+$"],
                ["S1", 1, +1, r"$\rm H$"],
                ["e-", 5e-4, -1, r"$\rm e^-$"]
            ]
            
            timesteps = range(10, 19, 1)
            densities_file = os.path.join(self.args.data_dir, self.args.densities_file)
            rates_file = os.path.join(self.args.data_dir, self.args.rates_file)
            
            avg_T, num_den = process_timesteps(
                directory=self.args.simulation_dir,
                sim=self.args.simulation_name,
                timesteps=timesteps,
                species=species,
                rplanet=1.37e8,
                index=10,
                output_file=densities_file
            )
            
            make_rates(
                output_file=rates_file,
                avg_T_at_index=avg_T,
                number_densities=num_den,
                timesteps=timesteps
            )
            
            print(f"Generated: {self.args.densities_file}, {self.args.rates_file}, qt_conditions.txt")
            return True
            
        except Exception as e:
            print(f"Error processing simulation data: {e}")
            return False
    
    def run_pumpkin_analysis(self):
        """Run PumpKin analysis across all cells."""
        print("=" * 60)
        print("STEP 2: Running PumpKin Analysis")
        print("=" * 60)
        
        try:
            input_file_path = os.path.join(self.args.pumpkin_in_dir, 'input.txt')
            
            if not os.path.exists(input_file_path):
                print(f"Error: PumpKin input file not found at {input_file_path}")
                return False
            
            # Generate input for PumpKin (species selection)
            input_data = "\n".join(str(i) for i in range(1, len(self.args.species) + 1)) + "\n-1\n"
            
            success_count = 0
            
            for cell_index in range(self.args.num_cells):
                print(f"Processing cell {cell_index + 1}/{self.args.num_cells}...")
                
                # Update input file for current cell
                if not self._modify_pumpkin_input(input_file_path, cell_index):
                    print(f"Warning: Failed to modify input file for cell {cell_index}")
                    continue
                
                # Run PumpKin
                command = f'cd "{os.path.abspath(self.args.pumpkin_dir)}" && ./PumpKin {os.path.join(os.path.basename(os.path.dirname(self.args.pumpkin_in_dir)), os.path.basename(self.args.pumpkin_in_dir))}/'
                result = subprocess.run(command, input=input_data, capture_output=True, text=True, shell=True)
                
                if result.returncode != 0:
                    print(f"Warning: PumpKin returned error code {result.returncode} for cell {cell_index}")
                    print(f"current directory: {os.path.abspath(self.args.pumpkin_dir)} command: {command}")
                
                # Save output
                output_combined = result.stdout + "\n--- STDERR ---\n" + result.stderr
                output_filename = os.path.join(self.args.output_dir, f"{self.args.output_prefix}_{cell_index:03d}.txt")
                
                with open(output_filename, "w") as f:
                    f.write(output_combined)
                
                success_count += 1
            
            print(f"SUCCESS: PumpKin analysis completed for {success_count}/{self.args.num_cells} cells")
            return success_count > 0
            
        except Exception as e:
            print(f"Error in PumpKin analysis: {e}")
            traceback.print_exc()
            return False
    
    def _modify_pumpkin_input(self, input_file_path, cell_index):
        """Update PumpKin input file for specific cell."""
        try:
            with open(input_file_path, "r") as f:
                lines = f.readlines()
            
            new_lines = []
            for line in lines:
                if line.strip().startswith("t_init"):
                    new_lines.append(f"t_init = {cell_index}\n")
                elif line.strip().startswith("t_end"):
                    new_lines.append(f"t_end = {cell_index + 1}\n")
                else:
                    new_lines.append(line)
            
            with open(input_file_path, "w") as f:
                f.writelines(new_lines)
            
            return True
            
        except Exception as e:
            print(f"Error modifying PumpKin input: {e}")
            return False
    
    def generate_plots(self):
        """Generate comprehensive plots using the plotting tool."""
        print("=" * 60)
        print("STEP 3: Generating Plots and Visualizations")
        print("=" * 60)
        
        try:
            # Import plotting main
            plotting_script = os.path.join(os.path.dirname(__file__), 'Outputs', 'plotting_main.py')
            
            if not os.path.exists(plotting_script):
                print(f"Error: Plotting script not found at {plotting_script}")
                return False
            
            # Build command for plotting tool
            plot_cmd = [
                sys.executable, plotting_script,
                '--all',
                '--data-dir', self.args.output_dir,
                '--output-dir', os.path.join(self.args.output_dir, 'Plots'),
                '--num-cells', str(self.args.num_cells),
                '--species'] + self.args.species + [
                '--top-n', str(self.args.top_n),
                '--species-file', os.path.join(self.args.data_dir, self.args.species_file),
                '--densities-file', os.path.join(self.args.data_dir, self.args.densities_file),
                '--reactions-file', os.path.join(self.args.data_dir, self.args.reactions_file),
                '--rates-file', os.path.join(self.args.data_dir, self.args.rates_file),
                '--temp-file', os.path.join(self.args.data_dir, 'qt_conditions.txt')
            ]
            
            print("Running plotting tool...")
            result = subprocess.run(plot_cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print("SUCCESS: Plot generation completed successfully")
                print("Generated plots saved in:", os.path.join(self.args.output_dir, 'Plots'))
                return True
            else:
                print("Warning: Some plots may have failed to generate")
                print("STDOUT:", result.stdout)
                print("STDERR:", result.stderr)
                return False
                
        except Exception as e:
            print(f"Error generating plots: {e}")
            traceback.print_exc()
            return False
    
    def run_complete_pipeline(self):
        """Run the complete analysis pipeline."""
        print("Starting Complete AIOLOS-to-PumpKin Analysis Pipeline")
        print(f"Configuration: {self.args.num_cells} cells, Species: {self.args.species}")
        print(f"Data directory: {self.args.data_dir}")
        print(f"Output directory: {self.args.output_dir}")
        
        success_stages = 0
        total_stages = 3
        
        # Stage 1: Process AIOLOS data
        if self.args.process_data:
            if self.process_aiolos_data():
                success_stages += 1
            else:
                print("ERROR: AIOLOS data processing failed")
        else:
            print("SKIPPING: AIOLOS data processing")
            total_stages -= 1
        
        # Stage 2: Run PumpKin analysis
        if self.args.run_pumpkin:
            if self.run_pumpkin_analysis():
                success_stages += 1
            else:
                print("ERROR: PumpKin analysis failed")
        else:
            print("SKIPPING: PumpKin analysis")
            total_stages -= 1
        
        # Stage 3: Generate plots
        if self.args.generate_plots:
            if self.generate_plots():
                success_stages += 1
            else:
                print("ERROR: Plot generation failed")
        else:
            print("SKIPPING: Plot generation")
            total_stages -= 1
        
        # Summary
        print("\n" + "=" * 60)
        print("PIPELINE SUMMARY")
        print("=" * 60)
        print(f"Stages completed: {success_stages}/{total_stages}")
        
        if success_stages == total_stages:
            print("SUCCESS: Complete pipeline executed successfully!")
            print(f"Results available in: {self.args.output_dir}")
            return True
        else:
            print("WARNING: Pipeline completed with some failures")
            return False

def create_parser():
    """Create command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Complete AIOLOS-to-PumpKin Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Complete analysis with default settings
    python main.py --run-all

    # Custom configuration
    python main.py --run-all --num-cells 100 --species H H2 O OH H2O --output-dir ./results/

    # Run individual stages
    python main.py --process-data --run-pumpkin --num-cells 50
    python main.py --generate-plots --species H H2 O OH

    # Advanced configuration
    python main.py --run-all --pumpkin-dir "/custom/path/Examples/AIOLOS_New" --top-n 10

For more information, see README.md
        """)
    
    # Main operation modes
    parser.add_argument('--run-all', action='store_true',
                       help='Run complete pipeline (data processing + PumpKin + plots)')
    parser.add_argument('--process-data', action='store_true',
                       help='Process AIOLOS reaction files and simulation data')
    parser.add_argument('--run-pumpkin', action='store_true',
                       help='Run PumpKin analysis across cells')
    parser.add_argument('--generate-plots', action='store_true',
                       help='Generate comprehensive plots and visualizations')
    
    # Directory and file configuration
    parser.add_argument('--data-dir', type=str, default='./data/',
                       help='Directory for input/output data files (default: ./data/)')
    parser.add_argument('--output-dir', type=str, default='./results/',
                       help='Directory for analysis results (default: ./results/)')
    parser.add_argument('--pumpkin-dir', type=str, 
                       default='/mnt/d/OneDrive/Water Worlds/PumpKin/src',
                       help='Path to PumpKin directory')
    parser.add_argument('--pumpkin_in-dir', type=str, 
                       default='/Examples/AIOLOS_NEW',
                       help='Path to PumpKin example directory')
    parser.add_argument('--simulation-dir', type=str, default='../dynamic_cond0_data/',
                       help='Directory containing simulation data')
    parser.add_argument('--simulation-name', type=str, 
                       default='dynamic_ignoreleectrontest-cell200-newsol2-h2e+2-long-cond0',
                       help='Simulation name for data processing')
    
    # Analysis parameters
    parser.add_argument('--num-cells', type=int, default=233,
                       help='Number of spatial cells to process (default: 233)')
    parser.add_argument('--species', nargs='+', default=['H', 'H2', 'O', 'OH'],
                       help='Species to analyze (default: H H2 O OH)')
    parser.add_argument('--top-n', type=int, default=5,
                       help='Number of top pathways to show in plots (default: 5)')
    
    # File naming
    parser.add_argument('--output-prefix', type=str, default='pumpkin_output_cell',
                       help='Prefix for PumpKin output files (default: pumpkin_output_cell)')
    parser.add_argument('--species-file', type=str, default='qt_species_list.txt',
                       help='Species list filename (default: qt_species_list.txt)')
    parser.add_argument('--reactions-file', type=str, default='qt_reactions_list.txt',
                       help='Reactions list filename (default: qt_reactions_list.txt)')
    parser.add_argument('--densities-file', type=str, default='qt_densities.txt',
                       help='Densities filename (default: qt_densities.txt)')
    parser.add_argument('--rates-file', type=str, default='qt_rates.txt',
                       help='Rates filename (default: qt_rates.txt)')
    parser.add_argument('--reaction-file', type=str, default='steamfull_step3.reac',
                       help='Input reaction file (default: steamfull_step3.reac)')
    
    return parser

def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Set default behavior if no specific mode is selected
    if not any([args.run_all, args.process_data, args.run_pumpkin, args.generate_plots]):
        print("No operation mode selected. Use --help for options or --run-all for complete pipeline.")
        return 1
    
    # Set flags for --run-all
    if args.run_all:
        args.process_data = True
        args.run_pumpkin = True
        args.generate_plots = True
    
    try:
        # Create and run pipeline
        pipeline = AIOLOSPumpKinPipeline(args)
        success = pipeline.run_complete_pipeline()
        
        return 0 if success else 1
        
    except Exception as e:
        print(f"Fatal error: {e}")
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())