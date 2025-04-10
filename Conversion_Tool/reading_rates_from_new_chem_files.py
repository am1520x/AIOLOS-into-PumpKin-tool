from getting_reactions_from_log import extract_all_reactions
from processing_aiolos_reac_file import process_reaction_line, transform_species_set
from making_densities_file import process_radial_profile
import pandas as pd
import numpy as np
import re

# Function to read the new .txt file
def read_new_file(file_path):
    # Load the data using np.loadtxt and capture the header (first row) for reaction numbers
    data = np.loadtxt(file_path, skiprows=1)  # Skip the first row (header) for now
    
    # Read the header line separately
    header = np.loadtxt(file_path, max_rows=1, dtype=int)

    # Create a DataFrame for the data
    # Assuming columns: CellNumber, PhysicalRadius, Rate1, Rate2, ...
    columns = ['cell_number', 'physical_radius'] + [f'{i}' for i in header] #+ header.tolist()
    
    # Create the DataFrame
    df = pd.DataFrame(data, columns=columns)
    #df['reaction_number'] = header
    return df

def process_side(side, multiplier, reaction_stoich, species_set):
    """Process one side of the reaction equation (reactants or products)."""
    term_pattern = re.compile(r'([\d\.]+)\s*([A-Za-z0-9\-\+\(\)^`]+)')
    terms = side.split('+')
    
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
    Process a column of reaction strings from a DataFrame and return stoichiometries and species.
    
    Args:
        reaction_series (pd.Series): Series of reaction strings like '1 A + 2 B -> 1 C'
    
    Returns:
        tuple:
            - stoich_list: List of dicts mapping species to net stoichiometric coefficients
            - transformed_species_list: Sorted list of transformed species
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