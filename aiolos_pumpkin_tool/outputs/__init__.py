"""
Outputs Module

Provides comprehensive plotting and visualization capabilities for chemical pathway analysis.

Key Functions:
- plotting_main: Main plotting orchestration
- species_specific_plots: Species-focused visualizations
- pathway_rates_plots: Pathway rate analysis plots
- rates_of_reactions_plots: Reaction rate visualizations
- deleted_pathways_plots: Analysis of removed pathways
"""

try:
    from .plotting_main import create_all_plots, main
    from .species_specific_plots import (
        plot_species_all_plots,
        load_all_cells_as_dict,
    )
    __all__ = [
        "create_all_plots",
        "main", 
        "plot_species_all_plots",
        "load_all_cells_as_dict",
    ]
except ImportError:
    # Handle missing dependencies gracefully
    __all__ = []