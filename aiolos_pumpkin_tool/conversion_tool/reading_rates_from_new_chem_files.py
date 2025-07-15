"""
Module for reading rate of reaction from chem output files and parsing 
reaction stoichiometry.

Functions:
- read_new_file: Loads tabular reaction rate data with reaction numbers as 
  column headers.
- process_side: Helper function to extract stoichiometric coefficients from 
  one side of a reaction.
- parse_reaction_column: Parses a series of reaction strings and returns 
  species-level stoichiometry.

"""

from .processing_aiolos_reac_file import transform_species_set
import pandas as pd
import numpy as np
import re


def read_new_file(file_path):
    """
    Reads a reaction rate output file with a header containing reaction numbers.

    Args:
        file_path (str): Path to the reaction rate file.

    Returns:
        pd.DataFrame: DataFrame with columns 
            ['cell_number', 'physical_radius', <reaction_ids>].

    Example:
        >>> df = read_new_file("data/reaction_rates.dat")
    """
    # Load the data using np.loadtxt and capture the header (first row) for 
    # reaction numbers
    data = np.loadtxt(file_path, skiprows=1)  # Skip the first row

    # Read the header line separately
    header = np.loadtxt(file_path, max_rows=1, dtype=int)

    # Create a DataFrame for the data
    columns = ["cell_number", "physical_radius"] + [f"{i}" for i in header]

    # Create the DataFrame
    df = pd.DataFrame(data, columns=columns)
    return df


def process_side(side, multiplier, reaction_stoich, species_set):
    """
    Processes one side (reactants or products) of a reaction string.

    Args:
        side (str): Reaction side string (e.g., "2 H + 1 O2").
        multiplier (int): -1 for reactants, +1 for products.
        reaction_stoich (dict): Dictionary to accumulate net stoichiometric coefficients.
        species_set (set): Set to accumulate all species encountered.

    Returns:
        None (modifies `reaction_stoich` and `species_set` in place)

    Example:
        >>> reaction_stoich = {}
        >>> species_set = set()
        >>> process_side("2 H + O2", -1, reaction_stoich, species_set)
    """
    term_pattern = re.compile(r"([\d\.]+)\s*([A-Za-z0-9\-\+\(\)^`]+)")
    terms = side.split("+")

    for term in terms:
        term = term.strip()
        match = term_pattern.match(term)
        if match:
            coeff_str, sp = match.groups()
            if sp.lower() == "ev":
                continue
            try:
                coeff = float(coeff_str)
            except ValueError:
                coeff = 0.0
            net_coeff = multiplier * coeff
            reaction_stoich[sp] = reaction_stoich.get(sp, 0) + net_coeff
            species_set.add(sp)


def parse_reaction_column(reaction_series):
    """
    Parses a column of reaction strings and returns stoichiometry per reaction.

    Args:
        reaction_series (pd.Series): Series of reaction strings like 
            '1 A + 2 B -> 1 C'.

    Returns:
        tuple:
            - stoich_list (list[dict]): List of dictionaries mapping species to 
              net coefficients.
            - transformed_species_list (list[str]): Sorted list of unique species 
              (transformed).

    Example:
        >>> series = pd.Series(["2 H + O2 -> 2 H2O"])
        >>> stoich, species = parse_reaction_column(series)
    """
    stoich_list = []
    species_set = set()

    for line in reaction_series:
        line = line.strip()
        if "->" in line:
            left, right = line.split("->")
        elif "=>" in line:
            left, right = line.split("=>")
        else:
            continue  # Skip if no recognizable arrow

        reaction_stoich = {}
        process_side(left.strip(), -1, reaction_stoich, species_set)
        process_side(right.strip(), +1, reaction_stoich, species_set)
        stoich_list.append(reaction_stoich)

    transformed_species_list = sorted(transform_species_set(species_set))
    return stoich_list, transformed_species_list
