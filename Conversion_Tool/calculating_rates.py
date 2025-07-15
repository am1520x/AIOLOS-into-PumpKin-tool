"""
Chemical reaction rate calculation module.

This module provides functions to calculate reaction rate coefficients using the 
modified Arrhenius equation and to compute overall reaction rates from 
reactant concentrations.

Functions:
    calc_rate_coeff: Calculate rate coefficient using modified Arrhenius equation
    parse_reactions: Parse reaction data from text files  
    calculate_reaction_rates: Calculate reaction rates from concentrations
"""

import numpy as np


def calc_rate_coeff(A, beta, gamma, T):
    """
    Calculate rate coefficient using the modified Arrhenius equation.
    
    The modified Arrhenius equation is: k = A * T^beta * exp(-gamma/T)
    where:
    - A is the pre-exponential factor (units depend on reaction order)
    - beta is the temperature exponent (dimensionless)  
    - gamma is the activation temperature (K)
    - T is the temperature (K)
    
    Args:
        A (float): Pre-exponential factor
        beta (float): Temperature exponent
        gamma (float): Activation temperature in Kelvin
        T (float): Temperature in Kelvin
        
    Returns:
        float: Rate coefficient, or np.nan if calculation fails
        
    Examples:
        >>> import numpy as np
        >>> # Simple reaction with A=1e-10, beta=0, gamma=0 at T=300K
        >>> k = calc_rate_coeff(1e-10, 0, 0, 300)
        >>> abs(k - 1e-10) < 1e-15
        True
        
        >>> # Temperature dependent reaction
        >>> k = calc_rate_coeff(1e-10, 0.5, 1000, 300)
        >>> isinstance(k, float) and k > 0
        True
        
        >>> # Invalid input should return NaN
        >>> k = calc_rate_coeff("invalid", 0, 0, 300)
        Error calculating rate coefficient: ...
        >>> np.isnan(k)
        True
    """
    try:     
        k = A * (T ** beta) * np.exp(-gamma / T)
        return k
    except (ValueError, TypeError, KeyError, ZeroDivisionError) as e:
        print(f"Error calculating rate coefficient: {e}")
        return np.nan 

def parse_reactions(file_path):
    """
    Parse reaction data from a text file containing reaction definitions.
    
    The file format expected has lines starting with '%' containing:
    % reaction_equation | alpha beta gamma
    
    For example:
    % 1 H + 1 O2 -> 1 OH + 1 O | 1.2e-10 0.5 1000
    
    Args:
        file_path (str): Path to the reaction file
        
    Returns:
        list: List of dictionaries, each containing:
            - reactants: List of tuples (species_name, stoichiometric_coeff)
            - alpha: Pre-exponential factor (float)
            - beta: Temperature exponent (float) 
            - gamma: Activation temperature (float)
            
    Raises:
        FileNotFoundError: If the specified file doesn't exist
        ValueError: If the file format is invalid
        
    Examples:
        >>> import tempfile, os
        >>> # Create a temporary reaction file
        >>> with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        ...     f.write("% 1 H + 1 O2 -> 1 OH + 1 O | 1.2e-10 0.5 1000\\n")
        ...     f.write("% 2 H2 + 1 O -> 1 H2O + 1 H | 3.4e-11 0.0 500\\n")
        ...     temp_file = f.name
        >>> reactions = parse_reactions(temp_file)
        >>> len(reactions)
        2
        >>> reactions[0]['reactants']
        [('H', 1), ('O2', 1)]
        >>> reactions[0]['alpha']
        1.2e-10
        >>> os.unlink(temp_file)  # Clean up
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
    Calculate reaction rates for all reactions using rate law kinetics.
    
    For each reaction, the rate is calculated as:
    rate = k * product([species_i]^coeff_i for all reactants)
    
    where k is the rate coefficient from the modified Arrhenius equation.
    
    Args:
        reactions (list): List of reaction dictionaries from parse_reactions()
        number_densities (dict): Species concentrations {species_name: concentration}
        calc_rate_coeff (function): Function to compute rate coefficient 
        T (float, optional): Temperature in Kelvin. Defaults to 5000.
        
    Returns:
        list: Reaction rates (float) for each reaction in molecules/cm³/s
        
    Examples:
        >>> import numpy as np
        >>> # Example reaction data
        >>> reactions = [
        ...     {'reactants': [('H', 1), ('O2', 1)], 'alpha': 1e-10, 'beta': 0, 'gamma': 0},
        ...     {'reactants': [('H2', 2)], 'alpha': 2e-11, 'beta': 0, 'gamma': 0}
        ... ]
        >>> densities = {'H': 1e12, 'O2': 1e11, 'H2': 1e10}
        >>> rates = calculate_reaction_rates(reactions, densities, calc_rate_coeff, T=300)
        >>> len(rates)
        2
        >>> rates[0] > 0  # Should be positive
        True
        >>> isinstance(rates[1], float)
        True
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