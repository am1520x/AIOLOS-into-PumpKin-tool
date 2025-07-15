# AIOLOS-into-PumpKin-tool Package Guide

## Package Structure

The project has been restructured as a proper Python package with the following organization:

```
AIOLOS-into-PumpKin-tool/
├── aiolos_pumpkin_tool/           # Main package directory
│   ├── __init__.py                # Package initialization
│   ├── main.py                    # Main pipeline class
│   ├── cli.py                     # Command-line interface
│   ├── utils.py                   # Utility functions
│   ├── conversion_tool/           # AIOLOS data conversion
│   ├── outputs/                   # Plotting and visualization
│   ├── pumpkin/                   # PumpKin integration
│   ├── tests/                     # Test suite
│   └── data/                      # Package data files
├── docs/                          # Sphinx documentation
├── pyproject.toml                 # Modern Python packaging config
├── setup.py                       # Backward compatibility
├── MANIFEST.in                    # Additional file inclusion
└── README.md                      # Project documentation
```

## Installation Options

### Option 1: Development Installation (Recommended for Development)
```bash
# Clone and install in editable mode
cd AIOLOS-into-PumpKin-tool
pip install -e .
```

### Option 2: Standard Installation
```bash
# Install from directory
pip install .
```

### Option 3: With Development Dependencies
```bash
# Install with additional dev tools
pip install -e ".[dev]"
```

## Usage After Installation

### Command Line Interface
Once installed, you can use the package from the command line:

```bash
# Run complete pipeline
aiolos-pumpkin --run-all

# Process only AIOLOS data
aiolos-pumpkin --process-data --num-cells 50

# Run only PumpKin analysis
aiolos-pumpkin --run-pumpkin

# Generate only plots
aiolos-pumpkin --generate-plots --species H H2 O OH

# Custom configuration
aiolos-pumpkin --run-all --data-dir ./my_data/ --output-dir ./my_results/
```

### Python API
You can also use the package programmatically:

```python
import aiolos_pumpkin_tool

# Use individual functions
species_set = aiolos_pumpkin_tool.transform_species_set({"H2", "O2"})
reaction = aiolos_pumpkin_tool.process_reaction_line("H2 + O2 -> H2O")

# Use the main pipeline class
from aiolos_pumpkin_tool import AIOLOSPumpKinPipeline

# Create pipeline with configuration
args = type('Args', (), {
    'run_all': True,
    'num_cells': 100,
    'data_dir': './data/',
    'output_dir': './results/'
})()

pipeline = AIOLOSPumpKinPipeline(args)
success = pipeline.run_complete_pipeline()
```

## Package Features

### Dependency Management
The package gracefully handles missing dependencies:
- **Core dependencies** (numpy, pandas): Required for basic functionality
- **Plotting dependencies** (matplotlib, seaborn): Optional, plotting features disabled if missing
- **PumpKin**: Auto-detected, warns if not available

### PumpKin Integration
The package handles PumpKin in several ways:
1. **Bundled**: If PumpKin is copied into the package
2. **Sibling directory**: If PumpKin is in the parent directory
3. **System installation**: Common installation paths
4. **Manual specification**: Via `--pumpkin-dir` argument

### Entry Points
Three command-line tools are available:
- `aiolos-pumpkin`: Main pipeline interface
- `aiolos-convert`: Conversion tools only
- `aiolos-plot`: Plotting tools only

## Distribution

### Creating a Distribution Package
```bash
# Build distribution packages
pip install build
python -m build

# This creates:
# dist/aiolos_pumpkin_tool-1.0.0.tar.gz        # Source distribution
# dist/aiolos_pumpkin_tool-1.0.0-py3-none-any.whl  # Wheel distribution
```

### Installing from Distribution
```bash
# Install from wheel
pip install dist/aiolos_pumpkin_tool-1.0.0-py3-none-any.whl

# Install from source
pip install dist/aiolos_pumpkin_tool-1.0.0.tar.gz
```

## Testing

```bash
# Run tests
python -m pytest aiolos_pumpkin_tool/tests/

# Run with coverage
pip install pytest-cov
python -m pytest aiolos_pumpkin_tool/tests/ --cov=aiolos_pumpkin_tool
```

## Documentation

```bash
# Build documentation
cd docs
make html

# View documentation
open _build/html/index.html
```

## Development Workflow

1. **Make changes** to the code
2. **Test changes**: `python -m pytest`
3. **Update version** in `pyproject.toml` if needed
4. **Build package**: `python -m build`
5. **Install and test**: `pip install dist/aiolos_pumpkin_tool-*.whl`

## Troubleshooting

### Import Warnings
Warnings about missing plot modules are normal if matplotlib/seaborn aren't installed:
```
Warning: Species plots module not available: No module named 'matplotlib'
```

### PumpKin Not Found
The package will warn if PumpKin can't be located:
```
Warning: PumpKin directory not found. PumpKin analysis will not be available.
```

### Command Not Found
If `aiolos-pumpkin` command isn't found, ensure the installation directory is in PATH:
```bash
export PATH=$PATH:~/.local/bin  # For user installations
```

## Benefits of Package Structure

1. **Proper dependency management**: Clear specification of required and optional dependencies
2. **Easy installation**: Standard `pip install` workflow
3. **Modular imports**: Import only what you need
4. **Command-line tools**: Installed globally after package installation
5. **Documentation**: Integrated Sphinx documentation
6. **Testing**: Standardized test structure
7. **Distribution**: Easy to share and deploy
8. **Version management**: Single source of truth for version information