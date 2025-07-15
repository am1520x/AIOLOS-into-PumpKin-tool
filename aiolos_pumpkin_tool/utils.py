"""
Utility functions for the AIOLOS-into-PumpKin-tool package.

Includes path management, PumpKin location detection, and configuration helpers.
"""

import os
import sys
from pathlib import Path


def find_pumpkin_directory():
    """
    Locate the PumpKin directory relative to the package installation.
    
    Returns:
        str: Path to PumpKin directory if found, None otherwise
    """
    # Possible locations for PumpKin relative to package
    possible_locations = [
        # In the package itself (bundled)
        os.path.join(os.path.dirname(__file__), "pumpkin"),
        # Sibling directory to the package
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "PumpKin"),
        # Parent directory (for development)
        os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "PumpKin"),
        # Common system locations
        "/usr/local/bin/PumpKin",
        "/opt/PumpKin",
        "~/PumpKin",
    ]
    
    for location in possible_locations:
        expanded_path = os.path.expanduser(location)
        if os.path.exists(expanded_path):
            return expanded_path
    
    return None


def setup_pumpkin_environment(pumpkin_dir=None):
    """
    Set up environment for PumpKin execution.
    
    Args:
        pumpkin_dir (str, optional): Path to PumpKin directory
        
    Returns:
        dict: Environment configuration for PumpKin
    """
    if pumpkin_dir is None:
        pumpkin_dir = find_pumpkin_directory()
    
    if pumpkin_dir is None:
        print("Warning: PumpKin directory not found. PumpKin analysis will not be available.")
        return None
    
    env_config = {
        "pumpkin_dir": pumpkin_dir,
        "pumpkin_executable": os.path.join(pumpkin_dir, "pumpkin"),
        "available": True
    }
    
    return env_config


def get_package_data_path(filename):
    """
    Get path to a data file within the package.
    
    Args:
        filename (str): Name of the data file
        
    Returns:
        str: Full path to the data file
    """
    package_dir = os.path.dirname(__file__)
    data_dir = os.path.join(package_dir, "data")
    return os.path.join(data_dir, filename)


def validate_dependencies():
    """
    Check if all required dependencies are available.
    
    Returns:
        dict: Status of each dependency
    """
    deps = {}
    
    # Check Python packages
    try:
        import numpy
        deps["numpy"] = True
    except ImportError:
        deps["numpy"] = False
    
    try:
        import pandas
        deps["pandas"] = True
    except ImportError:
        deps["pandas"] = False
    
    try:
        import matplotlib
        deps["matplotlib"] = True
    except ImportError:
        deps["matplotlib"] = False
    
    # Check PumpKin availability
    deps["pumpkin"] = find_pumpkin_directory() is not None
    
    return deps