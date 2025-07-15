"""
PumpKin Integration Module

Provides integration utilities for running PumpKin analysis and managing PumpKin workflows.

Key Functions:
- running_pumpkin: Single cell PumpKin execution
- running_all_pumpkin: Multi-cell PumpKin batch processing
"""

try:
    from .running_pumpkin import run_pumpkin_single_cell
    from .running_all_pumpkin import run_pumpkin_all_cells
    __all__ = [
        "run_pumpkin_single_cell",
        "run_pumpkin_all_cells",
    ]
except ImportError:
    # Handle cases where PumpKin utilities aren't available
    __all__ = []