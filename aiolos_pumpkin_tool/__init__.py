"""
AIOLOS into PumpKin Tool

A comprehensive pipeline for converting AIOLOS simulation data into PumpKin analysis format
and generating comprehensive chemical pathway visualizations.

This package provides:
- AIOLOS data conversion utilities
- PumpKin integration and automation
- Comprehensive plotting and visualization tools
- Complete pipeline management

Main Components:
- conversion_tool: AIOLOS data processing and conversion
- outputs: Plotting and visualization modules  
- pumpkin: PumpKin integration utilities
- main: Complete pipeline orchestration
"""

__version__ = "1.0.0"
__author__ = "Astrochemistry Team"
__email__ = "astrochemistry@example.com"

# Import main pipeline class for easy access
from .main import AIOLOSPumpKinPipeline

# Import key conversion functions
from .conversion_tool.processing_aiolos_reac_file import (
    transform_species_set,
    process_reaction_line,
    parse_reaction_data,
)

# Import plotting utilities (handle missing dependencies gracefully)
try:
    from .outputs.plotting_main import create_all_plots
except ImportError:
    create_all_plots = None

__all__ = [
    "AIOLOSPumpKinPipeline",
    "transform_species_set", 
    "process_reaction_line",
    "parse_reaction_data",
]

# Add create_all_plots to __all__ only if it was imported successfully
if create_all_plots is not None:
    __all__.append("create_all_plots")