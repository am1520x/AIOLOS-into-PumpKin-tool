"""
This script provides functions for parsing and transforming chemical reaction \ndata from the AIOLOS simulation for feeding into the PumpKin analysis code.

The functions here:
- Parse reaction lines from text into structured representations.
- Convert species and reaction strings to standardized notation.
- Generate stoichiometry mappings and a complete species list.

Useful for preprocessing input for chemical solvers or simulations.

Functions:
- transform_species_set: Standardizes species notation.
- process_reaction_line: Converts a raw reaction line into a cleaned-up format.
- process_side: Parses either side of a reaction and extracts stoichiometry.
- parse_reaction_data: Parses full reaction data from annotated text into \n  structured lists.
"""

import re


def transform_species_set(species_set):
    """
    Apply replacement rules to all species in the set.

    Args:
        species_set (set): Original set of species names.

    Returns:
        set: New set with transformed species names.

    Example:
        >>> transform_species_set({"e-", "p", "H2", "CO", "eV"})
        {'E', '^+', 'H2', 'CO', ''}
    """
    replacements = [
        ("e-", "E"),
        ("p", "^+"),
        ("->", "=>"),
        ("M", "ANY_NEUTRAL"),
        (" ", ""),
        ("eV", ""),
        ("-", "^-"),
    ]

    transformed_set = set()
    for species in species_set:
        modified_species = species
        for old, new in replacements:
            modified_species = modified_species.replace(old, new)
        transformed_set.add(modified_species)

    return transformed_set


def process_reaction_line(line):
    """
    Process and standardize a single reaction line string.

    Removes stoichiometric numbers and applies standardized replacements to \n    make it like PumpKin inputs.

    Args:
        line (str): A raw reaction string (e.g., from a log or input file).

    Returns:
        str: A standardized reaction string.

    Example:
        >>> process_reaction_line("2 H2 + O2 -> 2 H2O")
        'H2+O2=>H2O'
        >>> process_reaction_line("2 Hp + O- -> 2 H2O")
        'H^++O^-=>H2O'
    """
    # Remove stoichiometric numbers
    modified = re.sub(r"\b[\d\.]+\s+", "", line)

    # Replace various notations
    replacements = [
        ("e-", "E"),
        ("p", "^+"),
        ("->", "=>"),
        ("M", "ANY_NEUTRAL"),
        ("+ eV", ""),
        ("-", "^-"),
    ]
    for old, new in replacements:
        modified = modified.replace(old, new)

    # Fix ion notation
    modified = re.sub(r"([A-Za-z0-9]+)([+-])", r"\1^\2", modified)
    modified = modified.replace(" ", "")
    return modified


def process_side(side, multiplier, reaction_stoich, species_set):
    """
    Process one side of a reaction to extract stoichiometric coefficients and species.

    Updates the reaction_stoich dictionary and species_set in place.

    Args:
        side (str): The reaction side string (reactants or products).
        multiplier (int): -1 for reactants, +1 for products.
        reaction_stoich (dict): Dictionary of stoichiometry (modified in place).
        species_set (set): Set of species encountered (modified in place).

    Returns:
        None

    Example:
        >>> stoich, species = {}, set()
        >>> process_side("2 H2 + O2", -1, stoich, species)
        >>> stoich
        {'H2': -2.0, 'O2': -1.0}
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


def parse_reaction_data(reac_text):
    """
    Parse reaction data from text and return reaction components.
    Handles reaction lines starting with `$` (thermo) or `%` (photo), splits equations,
    calculates net stoichiometry, and identifies unique species.

    Args:
        reac_text (str): The input text containing reaction data

    Returns:
        tuple: (reaction_list, stoich_list, species_list)
            - reaction_list: List of reaction equations
            - stoich_list: List of stoichiometry dictionaries
            - species_list: Sorted list of unique species
    Example:
        >>> txt = '''
        ... $ H2 + O -> OH + H | example
        ... % CO + photon -> C + O | photolysis
        ... '''
        >>> r, s, species = parse_reaction_data(txt)
        >>> r[0]
        'H2 + O -> OH + H'
        >>> s[0]['H2']
        -1.0
        >>> 'O' in species
        True
    """
    reaction_list = []
    stoich_list = []
    species_set = set()

    for line in reac_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        if line.startswith(("$", "%")):
            clean_line = line[1:].strip()
            parts = clean_line.split("|")
            if not parts:
                continue

            eq_str = parts[0].strip()
            reaction_list.append(eq_str)

            # Process stoichiometry
            reaction_stoich = {}
            if "->" in eq_str:
                left, right = eq_str.split("->")
            elif "=>" in eq_str:
                left, right = eq_str.split("=>")
            else:
                continue

            process_side(left.strip(), -1, reaction_stoich, species_set)
            process_side(right.strip(), +1, reaction_stoich, species_set)
            stoich_list.append(reaction_stoich)

    species_list = sorted(species_set)
    return reaction_list, stoich_list, species_list
