import re

# --- Helper Functions ---
def process_reaction_line(line):
    """Process a single reaction line and return modified components."""
    # Remove stoichiometric numbers
    modified = re.sub(r'\b[\d\.]+\s+', '', line)
    
    # Replace various notations
    replacements = [
        ('e-', 'E'),
        ('p', '^+'),
        ('->', '=>'),
        ('M', 'ANY_NEUTRAL'),
        (' ', '')
    ]
    for old, new in replacements:
        modified = modified.replace(old, new)
    
    # Fix ion notation
    modified = re.sub(r'([A-Za-z0-9]+)([+-])', r'\1^\2', modified)
    
    return modified

def process_side(side, multiplier, reaction_stoich, species_set):
    """Process one side of reaction equation (reactants or products)."""
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

def parse_reaction_data(reac_text):
    """
    Parse reaction data from text and return reaction components.
    
    Args:
        reac_text (str): The input text containing reaction data
        
    Returns:
        tuple: (reaction_list, stoich_list, species_list)
            - reaction_list: List of reaction equations
            - stoich_list: List of stoichiometry dictionaries
            - species_list: Sorted list of unique species
    """
    reaction_list = []
    stoich_list = []
    species_set = set()

    for line in reac_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
            
        if line.startswith(("$", "%")):
            reaction_type = line[0]
            clean_line = line[1:].strip()
            parts = clean_line.split("|")
            if not parts:
                continue

            eq_str = parts[0].strip()
            reaction_list.append(eq_str)  # Store original for processing
            
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