"""
Main script to process chemical reaction data from AIOLOS and generate inputs \nfor PumpKin analysis.

This script:
- Extracts reactions from a log file
- Loads chemistry data from a `.dat` file
- Merges reaction rates with their corresponding reaction strings
- Processes and outputs:
    - Species list (`qt_species_list.txt`)
    - Formatted reaction list (`qt_reactions_list.txt`)
    - Stoichiometry matrix (`qt_matrix.txt`)
    - Radial temperature and density profile (`qt_densities.txt`, `qt_conditions.txt`)
    - Reaction rate matrix (`qt_rates.txt`)

Usage (via command line):
    python main.py --logfile LOGFILE --chemfile CHEMFILE --sim SIM_NAME \\
                   --directory DIR --timestep TIMESTEP

Dependencies:
- pandas
- argparse
- numpy
- Local modules: getting_reactions_from_log, processing_aiolos_reac_file,
                 making_densities_file, reading_rates_from_new_chem_files
"""

import argparse
import pandas as pd
from getting_reactions_from_log import extract_all_reactions
from processing_aiolos_reac_file import process_reaction_line
from making_densities_file import process_radial_profile
from reading_rates_from_new_chem_files import read_new_file, parse_reaction_column


def main():
    """
    Main entry point for processing chemistry files and generating required outputs.

    Steps:
        1. Extracts reactions from the log file.
        2. Loads reaction rates from the specified chemistry file.
        3. Merges rate and reaction data for further processing.
        4. Outputs species list, reactions, and stoichiometric matrix.
        5. Applies Aiolos-specific species name formatting.
        6. Processes temperature and density profiles for the given timestep.
        7. Outputs a matrix of reaction rates by cell.

    Command-line Args:
        --logfile (str): Path to the log file.
        --chemfile (str): Path to the chemistry .dat file.
        --sim (str): Base simulation file name (no extension).
        --directory (str): Directory containing simulation output files.
        --timestep (int): Timestep to process.

    Outputs:
        qt_species_list.txt    : List of species in transformed form
        qt_reactions_list.txt  : List of reactions in processed form
        qt_matrix.txt          : Stoichiometric matrix
        qt_densities.txt       : Species number densities (radial)
        qt_conditions.txt      : Temperature profile (radial)
        qt_rates.txt           : Rate matrix by cell
    """
    parser = argparse.ArgumentParser(
        description="Process chemistry and reaction data."
    )
    parser.add_argument("--logfile", required=True, help="Path to the log file")
    parser.add_argument(
        "--chemfile", required=True, help="Path to the chemistry .dat file"
    )
    parser.add_argument(
        "--sim", required=True, help="Simulation file base name (without extension)"
    )
    parser.add_argument(
        "--directory", required=True, help="Directory containing the simulation files"
    )
    parser.add_argument("--timestep", type=int, required=True, help="Timestep value")

    args = parser.parse_args()

    # Step 1: Extract reactions
    reactions_df = extract_all_reactions(args.logfile)

    # Step 2: Read the chemistry file
    rates_df = read_new_file(args.chemfile)
    print(rates_df.head())

    # Step 3: Melt the rates DataFrame
    rates_df_melted = pd.melt(
        rates_df,
        id_vars=["cell_number", "physical_radius"],
        var_name="reaction_number",
        value_name="rate",
    )
    rates_df_melted["reaction_number"] = rates_df_melted["reaction_number"].astype(int)
    reactions_df["number"] = reactions_df["number"].astype(int)

    # Step 4: Merge to get reaction equations
    merged_df = pd.merge(
        rates_df_melted,
        reactions_df,
        how="left",
        left_on="reaction_number",
        right_on="number",
    )
    merged_df = merged_df.drop(columns=["number"])
    print(merged_df)

    # Step 5: Process reactions
    merged_df["processed_reaction"] = merged_df["reaction"].apply(
        process_reaction_line
    )
    print(merged_df[["reaction", "processed_reaction"]].head())

    # Only process cell 0
    df = merged_df[merged_df["cell_number"] == 0.0]
    stoich_list, species_list = parse_reaction_column(df["reaction"])

    # Write species list
    with open("qt_species_list.txt", "w") as f:
        for i, sp in enumerate(species_list, start=1):
            f.write(f"{i} {sp}\n")

    # Write reactions
    with open("qt_reactions_list.txt", "w") as f:
        for idx, reaction in enumerate(df["processed_reaction"], start=1):
            f.write(f"{idx} {reaction}\n")

    # Write stoichiometry matrix
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

    print("Files generated: qt_species_list.txt, qt_reactions_list.txt, qt_matrix.txt")

    # Process for Aiolos
    replacements = [("E", "e-"), ("^+", "p"), ("^-", "-")]
    aiolos_species = []
    for species in species_list:
        modified_species = species
        for old, new in replacements:
            modified_species = modified_species.replace(old, new)
        aiolos_species.append(modified_species)

    print(aiolos_species)

    # Radial profile processing
    avg_T_list, radial_data = process_radial_profile(
        args.directory, args.sim, args.timestep, aiolos_species
    )

    # Pivot and save rate matrix
    pivoted = merged_df.pivot(
        index="cell_number", columns="reaction_number", values="rate"
    )
    pivoted.reset_index(inplace=True)
    pivoted.to_csv("qt_rates.txt", sep="\t", index=False)

    print("Generated qt_rates.txt")


if __name__ == "__main__":
    main()
