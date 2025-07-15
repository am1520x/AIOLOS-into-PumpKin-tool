"""
plot_styles.py - Consistent styling configuration for all plotting scripts

This module provides standardized colors, line styles, and formatting
for all plots to ensure visual consistency across the analysis pipeline.
"""

import matplotlib.pyplot as plt
import re

# Define consistent color palette for species
SPECIES_COLORS = {
    # Hydrogen species - Blues
    'H': '#1f77b4',      # blue
    'H2': '#aec7e8',     # light blue
    'H-': '#0066cc',     # dark blue
    'H2p': '#3399ff',    # medium blue
    'H3p': '#6bb6ff',    # light medium blue
    'H2O': '#87ceeb',    # sky blue
    'H2O2': '#4682b4',   # steel blue
    'H2Op': '#5dade2',   # medium sky blue
    'H3Op': '#7fb3d3',   # light steel blue
    
    # Oxygen species - Reds
    'O': '#d62728',      # red
    'O2': '#ff7f7f',     # light red
    'O-': '#cc0000',     # dark red
    'Op': '#ff3333',     # medium red
    'O2p': '#ff6666',    # light medium red
    'O2-': '#990000',    # dark dark red
    'OH': '#ff9999',     # light light red
    'OH-': '#cc3333',    # medium dark red
    'OHp': '#e74c3c',    # tomato red
    'O2H': '#c0392b',    # dark tomato
    'O2Hp': '#ec7063',   # light tomato
    
    # Carbon species - Greens
    'C': '#2ca02c',      # green
    'CO': '#98df8a',     # light green
    'CO2': '#1f5f1f',    # dark green
    'CH': '#52c41a',     # bright green
    'CH2': '#73d13d',    # medium green
    'CH3': '#95de64',    # light medium green
    'CH4': '#b7eb8f',    # very light green
    'HCO': '#389e0d',    # dark bright green
    'HCOp': '#52c41a',   # bright green
    
    # Nitrogen species - Purples
    'N': '#9467bd',      # purple
    'N2': '#c5b0d5',     # light purple
    'NH': '#7b68ee',     # medium slate blue
    'NH2': '#9370db',    # medium purple
    'NH3': '#ba55d3',    # medium orchid
    'NO': '#8a2be2',     # blue violet
    'NO2': '#9932cc',    # dark orchid
    'N2O': '#dda0dd',    # plum
    
    # Sulfur species - Yellows/Oranges
    'S': '#ff7f0e',      # orange
    'S2': '#ffbb78',     # light orange
    'SO': '#ff8c00',     # dark orange
    'SO2': '#ffa500',    # orange
    'H2S': '#ff6347',    # tomato orange
    'HS': '#ff4500',     # orange red
    'OCS': '#ff8c69',    # salmon
    
    # Metal species - Browns
    'Fe': '#8c564b',     # brown
    'Mg': '#c49c94',     # light brown
    'Ca': '#d2691e',     # chocolate
    'Na': '#daa520',     # golden rod
    'K': '#b8860b',      # dark golden rod
    'Al': '#cd853f',     # peru
    'Si': '#a0522d',     # sienna
    'Ti': '#8b4513',     # saddle brown
    
    # Other species - Grays and others
    'M': '#7f7f7f',      # gray (third body)
    'He': '#bcbd22',     # olive
    'Ar': '#17becf',     # cyan
    'Ne': '#e377c2',     # pink
    'e-': '#000000',     # black (electrons)
    'PHOTON': '#ffff00', # yellow
    'CRPHOT': '#ffd700', # gold
    'ANY_NEUTRAL': '#808080', # gray
    
    # Common ions - based on parent species but modified
    'Hp': '#0066cc',     # dark blue (H+)
    'Cp': '#1f5f1f',     # dark green (C+)
    'Np': '#6a5acd',     # slate blue (N+)
    'Sp': '#ff4500',     # orange red (S+)
    'Fep': '#654321',    # dark brown (Fe+)
    'Mgp': '#8b4513',    # saddle brown (Mg+)
    'Sip': '#696969',    # dim gray (Si+)
}

def get_species_style(species_name):
    """
    Get consistent color and line style for a species.
    
    Args:
        species_name (str): Name of the species
        
    Returns:
        dict: Dictionary with 'color' and 'linestyle' keys
    """
    # Determine if species is ionic
    is_positive_ion = species_name.endswith('p') or species_name.endswith('+')
    is_negative_ion = species_name.endswith('-') or species_name.endswith('m')
    is_neutral = not (is_positive_ion or is_negative_ion)
    
    # Get color
    color = SPECIES_COLORS.get(species_name, '#333333')  # default dark gray
    
    # If species not in palette, try to infer from base species
    if species_name not in SPECIES_COLORS:
        # Try removing charge indicators
        base_species = species_name.replace('p', '').replace('+', '').replace('-', '').replace('m', '')
        if base_species in SPECIES_COLORS:
            color = SPECIES_COLORS[base_species]
    
    # Set line style based on charge
    if is_positive_ion or is_negative_ion:
        linestyle = '--'  # dashed for ionic species
    else:
        linestyle = '-'   # solid for neutral species
    
    return {
        'color': color,
        'linestyle': linestyle
    }

def get_species_color(species_name):
    """Get only the color for a species."""
    return get_species_style(species_name)['color']

def get_species_linestyle(species_name):
    """Get only the line style for a species."""
    return get_species_style(species_name)['linestyle']

def apply_plot_style():
    """Apply consistent styling to matplotlib plots."""
    plt.style.use('default')
    
    # Set consistent plot parameters
    plt.rcParams.update({
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'figure.figsize': (10, 6),
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'legend.frameon': True,
        'legend.fancybox': True,
        'legend.shadow': True,
        'grid.alpha': 0.3,
        'axes.grid': True,
        'axes.axisbelow': True,
        'lines.linewidth': 2,
        'lines.markersize': 6,
    })

def format_species_name(species_name):
    """
    Format species name for display in plots.
    
    Args:
        species_name (str): Raw species name
        
    Returns:
        str: Formatted species name for display
    """
    # Replace common formatting
    formatted = species_name
    
    # Handle ionic species
    if formatted.endswith('p'):
        formatted = formatted[:-1] + '⁺'
    elif formatted.endswith('+'):
        formatted = formatted[:-1] + '⁺'
    elif formatted.endswith('-'):
        formatted = formatted[:-1] + '⁻'
    elif formatted.endswith('m'):
        formatted = formatted[:-1] + '⁻'
    
    # Handle subscripts for common molecules
    formatted = re.sub(r'(\d+)', r'₍\1₎', formatted)  # Convert numbers to subscripts
    
    # Common replacements
    replacements = {
        'H2': 'H₂',
        'O2': 'O₂',
        'CO2': 'CO₂',
        'H2O': 'H₂O',
        'H2O2': 'H₂O₂',
        'NH3': 'NH₃',
        'CH4': 'CH₄',
        'SO2': 'SO₂',
        'H2S': 'H₂S',
        'N2O': 'N₂O',
        'NO2': 'NO₂',
        'H3Op': 'H₃O⁺',
        'H2Op': 'H₂O⁺',
        'H2p': 'H₂⁺',
        'H3p': 'H₃⁺',
        'O2p': 'O₂⁺',
        'O2-': 'O₂⁻',
        'OH-': 'OH⁻',
        'OHp': 'OH⁺',
        'Hp': 'H⁺',
        'Op': 'O⁺',
        'e-': 'e⁻',
    }
    
    for old, new in replacements.items():
        formatted = formatted.replace(old, new)
    
    return formatted

def create_species_legend(species_list, ax=None):
    """
    Create a legend showing species colors and line styles.
    
    Args:
        species_list (list): List of species names
        ax (matplotlib.axes.Axes, optional): Axes to add legend to
        
    Returns:
        matplotlib.legend.Legend: The created legend
    """
    if ax is None:
        ax = plt.gca()
    
    # Create legend elements
    legend_elements = []
    for species in sorted(set(species_list)):
        style = get_species_style(species)
        formatted_name = format_species_name(species)
        
        legend_elements.append(
            plt.Line2D([0], [0], 
                      color=style['color'], 
                      linestyle=style['linestyle'],
                      linewidth=2,
                      label=formatted_name)
        )
    
    return ax.legend(handles=legend_elements, loc='best')

# Color palette for pathway types (for pathway-specific plots)
PATHWAY_COLORS = {
    'production': '#2ecc71',      # green
    'consumption': '#e74c3c',     # red
    'net': '#3498db',             # blue
    'deleted': '#f39c12',         # orange
    'important': '#9b59b6',       # purple
    'minor': '#95a5a6',           # gray
}

def get_pathway_color(pathway_type):
    """Get color for pathway type."""
    return PATHWAY_COLORS.get(pathway_type, '#333333')

# Apply the style when module is imported
apply_plot_style()