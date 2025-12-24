# MCP GROMACS Scripts

Clean, self-contained scripts extracted from verified use cases for MCP tool wrapping.

## Design Principles

1. **Minimal Dependencies**: Only essential packages imported
2. **Self-Contained**: Functions inlined where possible
3. **Configurable**: Parameters in config files, not hardcoded
4. **MCP-Ready**: Each script has a main function ready for MCP wrapping

## Scripts

| Script | Description | Repo Dependent | Config | Example |
|--------|-------------|----------------|--------|---------|
| `md_simulation.py` | Run MD simulations with GROMACS | No | `configs/md_simulation_config.json` | `python md_simulation.py -i sample.tpr -o results/` |
| `tpr_analysis.py` | Analyze TPR files | No | `configs/tpr_analysis_config.json` | `python tpr_analysis.py -i sample.tpr -o analysis.txt` |
| `gromacs_command.py` | Execute GROMACS commands | No | `configs/gromacs_command_config.json` | `python gromacs_command.py --list-commands` |
| `gromacs_workflow.py` | Manage GROMACS workflows | No | `configs/gromacs_workflow_config.json` | `python gromacs_workflow.py --list-workflows` |

## Dependencies

All scripts use only standard library imports:
- `argparse` - Command line parsing
- `sys`, `os`, `pathlib` - System operations
- `logging` - Logging
- `subprocess` - GROMACS command execution
- `json` - Configuration loading
- `re` (tpr_analysis only) - Pattern matching

External dependency: GROMACS (accessed via `gmx` command)

## Usage

### Environment Setup

```bash
# Activate environment (prefer mamba over conda)
mamba activate ./env  # or: conda activate ./env
```

### Basic Usage

```bash
# Run MD simulation
python scripts/md_simulation.py --input examples/data/sample.tpr --output results/md_sim

# Analyze TPR file
python scripts/tpr_analysis.py --input examples/data/sample.tpr --output results/analysis.txt

# List available GROMACS commands
python scripts/gromacs_command.py --list-commands

# List available workflows
python scripts/gromacs_workflow.py --list-workflows

# Run demo workflow
python scripts/gromacs_workflow.py --workflow demo --verbose
```

### With Custom Configuration

```bash
# Use custom config
python scripts/md_simulation.py --input FILE --output DIR --config configs/custom.json
```

## Shared Library

Common functions are in `scripts/lib/`:
- `gromacs_utils.py`: GROMACS command utilities and parsing functions

## Configuration Files

Configuration files in `configs/` allow customization of:
- Simulation parameters
- Command timeouts
- Output formats
- File paths

## For MCP Wrapping (Step 6)

Each script exports a main function that can be wrapped:

```python
from scripts.md_simulation import run_md_simulation

# In MCP tool:
@mcp.tool()
def run_simulation(input_file: str, output_dir: str = None):
    return run_md_simulation(input_file, output_dir)
```

## Testing

All scripts have been tested with example data:

```bash
# Test with demo data
python scripts/md_simulation.py --input examples/data/sample.tpr --output results/test_md
python scripts/tpr_analysis.py --input examples/data/sample.tpr --output results/test_analysis.txt
python scripts/gromacs_workflow.py --workflow demo --verbose
```

## Script Architecture

Each script follows a consistent structure:

1. **Minimal Imports**: Only essential packages
2. **Configuration**: Default config with override capability
3. **Core Function**: Main logic (e.g., `run_md_simulation()`)
4. **CLI Interface**: `argparse` with help and examples
5. **Error Handling**: Comprehensive error checking
6. **Logging**: Structured logging throughout

This design makes scripts easy to:
- Test independently
- Wrap as MCP tools
- Configure for different use cases
- Debug and maintain