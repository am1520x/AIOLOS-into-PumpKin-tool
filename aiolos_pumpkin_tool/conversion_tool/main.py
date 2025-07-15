#!/usr/bin/env python3
"""
main.py - Primary execution script for processing Aiolos reaction files and simulation data
"""

import argparse
import os
from .processing_aiolos_reac_file import (
    parse_reaction_data,
    process_reaction_line,
    transform_species_set,
)
from .making_densities_file import process_timesteps, process_radial_profile
from .making_rates_file import make_rates


def process_reaction_files(input_file="steamfull_step3.reac"):
    """Process Aiolos reaction files and generate output files"""
    try:
        with open(input_file, "r") as f:
            reac_text = f.read()

        reaction_list, stoich_list, species_list = parse_reaction_data(reac_text)
        species_list = transform_species_set(species_list)

        # Generate output files
        with open("qt_species_list.txt", "w") as f:
            for i, sp in enumerate(species_list, start=1):
                f.write(f"{i} {sp}\n")

        with open("qt_reactions_list.txt", "w") as f:
            for idx, reaction in enumerate(reaction_list, start=1):
                modified = process_reaction_line(reaction)
                f.write(f"{idx} {modified}\n")

        with open("qt_matrix.txt", "w") as f:
            for sp in species_list:
                row = [
                    (
                        str(int(coeff))
                        if isinstance(coeff, (int, float)) and coeff.is_integer()
                        else str(coeff)
                    )
                    for coeff in (r.get(sp, 0) for r in stoich_list)
                ]
                f.write("\t".join(row) + "\n")

        print("Successfully generated reaction files:")
        print("- qt_species_list.txt\n- qt_reactions_list.txt\n- qt_matrix.txt")
        return True

    except FileNotFoundError:
        print(f"Error: Reaction file {input_file} not found.")
        return False
    except Exception as e:
        print(f"Error processing reaction files: {str(e)}")
        return False


def process_simulation_data(
    processing_type="timesteps", output_folder="Testing", timestep=None
):
    """Process simulation timesteps and generate density/rate files

    Args:
        processing_type (str): Either "timesteps" or "radial_profile"
        output_folder (str): Name of the output folder (default: "Testing")
        timestep (int): Single timestep for radial profile processing
    """
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
        ["e-", 5e-4, -1, r"$\rm e^-$"],
    ]

    # Create output directory
    output_dir = f"/mnt/d/OneDrive/Water Worlds/PumpKin/src/Examples/{output_folder}"
    os.makedirs(output_dir, exist_ok=True)

    # Define output files with full paths
    output_files = {
        "densities": os.path.join(output_dir, "qt_densities.txt"),
        "rates": os.path.join(output_dir, "qt_rates.txt"),
        "conditions": os.path.join(output_dir, "qt_conditions.txt"),
    }

    try:
        print(f"\nProcessing simulation data using {processing_type} method...")

        if processing_type == "timesteps":
            timesteps = range(10, 19, 1)
            avg_T, num_den = process_timesteps(
                directory="../dynamic_cond0_data/",
                sim="dynamic_ignoreleectrontest-cell200-newsol2-h2e+2-long-cond0",
                timesteps=timesteps,
                species=species,
                rplanet=1.37e8,
                index=10,
                output_file=output_files["densities"],
                output_file_2=output_files["conditions"],
            )
        elif processing_type == "radial_profile":
            if timestep is None:
                timestep = 15  # Default timestep
            timesteps = [timestep]
            avg_T, num_den = process_radial_profile(
                directory="../dynamic_cond0_data/",
                sim="dynamic_ignoreleectrontest-cell200-newsol2-h2e+2-long-cond0",
                timestep=timestep,
                species=species,
                rplanet=1.37e8,
                output_file=output_files["densities"],
                output_file_2=output_files["conditions"],
            )
        else:
            raise ValueError(f"Unknown processing type: {processing_type}")

        make_rates(
            output_file=output_files["rates"],
            avg_T_at_index=avg_T,
            number_densities=num_den,
            timesteps=timesteps,
        )
        print(f"\nSuccessfully generated files in {output_dir}:")
        print(f"- {os.path.basename(output_files['densities'])}")
        print(f"- {os.path.basename(output_files['rates'])}")
        print(f"- {os.path.basename(output_files['conditions'])}")
        return True

    except Exception as e:
        print(f"Error processing simulation data: {str(e)}")
        return False


def main():
    """Main function with command-line interface"""
    parser = argparse.ArgumentParser(
        description="Process Aiolos reaction files and simulation data"
    )
    parser.add_argument(
        "--process-data", action="store_true", help="Process simulation data only"
    )
    parser.add_argument(
        "--processing-type",
        choices=["timesteps", "radial_profile"],
        default="timesteps",
        help="Type of processing: timesteps or radial_profile (default: timesteps)",
    )
    parser.add_argument(
        "--output-folder",
        default="Testing",
        help="Output folder name (default: Testing)",
    )
    parser.add_argument(
        "--timestep", type=int, help="Timestep for radial profile processing"
    )

    args = parser.parse_args()

    if args.process_data:
        print("=== Processing Simulation Data Only ===")
        if not process_simulation_data(
            args.processing_type, args.output_folder, args.timestep
        ):
            print("Error: Simulation data processing failed")
            exit(1)
    else:
        print("=== Aiolos Data Processing Pipeline ===")

        # Step 1: Process reaction files
        if not process_reaction_files():
            print(
                "Warning: Reaction file processing failed, "
                "continuing with simulation data..."
            )

        # Step 2: Process simulation data
        if not process_simulation_data(
            args.processing_type, args.output_folder, args.timestep
        ):
            print("Error: Simulation data processing failed")
            exit(1)

    print("\nProcessing completed successfully!")


if __name__ == "__main__":
    main()
