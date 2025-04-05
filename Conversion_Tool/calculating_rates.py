import numpy as np

def calc_rate_coeff(A, beta, gamma, T):
  try:     
    k = A * (T ** beta) * np.exp(-gamma / T)
    return k
  except (ValueError, KeyError) as e:
     print(f"Error calculating rate coefficient: {e}")
     return np.nan 

def parse_reactions(file_path):
    """
    Parses the reaction data from a .txt file.
    Returns a list of dictionaries, each containing:
      - reactants: List of tuples (species, stoichiometric_coeff)
      - alpha, beta, gamma: Modified Arrhenius parameters
    """
    reactions = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('%'):
                # Split reaction and parameters
                parts = line[1:].split('|')
                reaction_part = parts[0].strip()
                params = parts[1].strip().split()
                # Extract Arrhenius parameters
                alpha = float(params[0])
                beta = float(params[1])
                gamma = float(params[2])
                
                # Split reactants and products
                reactants_products = reaction_part.split('->')
                reactants_str = reactants_products[0].strip()
                
                # Parse reactants (species and coefficients)
                reactants = []
                for species_part in reactants_str.split('+'):
                    species_part = species_part.strip()
                    coeff, species = species_part.split(' ', 1)
                    reactants.append((species.strip(), int(coeff)))
                reactions.append({
                    'reactants': reactants,
                    'alpha': alpha,
                    'beta': beta,
                    'gamma': gamma
                })
    return reactions

def calculate_reaction_rates(reactions, number_densities, calc_rate_coeff, T=5000):
    """
    Calculates reaction rates for all reactions.
    
    Args:
        reactions: List of parsed reactions from parse_reactions()
        number_densities: Dict {species: concentration}
        calc_rate_coeff: Function to compute k from (alpha, beta, gamma)
    
    Returns:
        List of reaction rates (one per reaction)
    """
    rates = []
    for rxn in reactions:
        # Calculate rate coefficient
        k = calc_rate_coeff(rxn['alpha'], rxn['beta'], rxn['gamma'], T)
        
        # Calculate rate = k * product([species]^coeff)
        rate = k
        for species, coeff in rxn['reactants']:
            conc = number_densities.get(species, 0.0)  # Default to 0 if missing
            rate *= (conc ** coeff)
        
        rates.append(rate)
    
    return rates