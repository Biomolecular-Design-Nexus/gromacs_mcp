# gromacs-2025.4 MCP

> MCP server for GROMACS 2025.4 providing molecular dynamics simulation, analysis, and workflow management tools

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Local Usage (Scripts)](#local-usage-scripts)
- [MCP Server Installation](#mcp-server-installation)
- [Using with Claude Code](#using-with-claude-code)
- [Using with Gemini CLI](#using-with-gemini-cli)
- [Available Tools](#available-tools)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)

## Overview

This MCP server provides comprehensive access to GROMACS 2025.4 molecular dynamics tools through both local scripts and MCP tools. It supports synchronous operations for quick analysis and asynchronous job submission for long-running simulations.

### Features
- **Fast Analysis Tools**: TPR file analysis, GROMACS commands, workflow management
- **Long-Running Simulations**: Background MD simulation with job tracking
- **Batch Processing**: Multiple file analysis in single operations
- **Job Management**: Full lifecycle management with status tracking, logs, and results
- **Comprehensive Workflows**: Predefined simulation and analysis workflows

### Directory Structure
```
./
├── README.md               # This file
├── env/                    # Conda environment
├── src/
│   ├── server.py           # MCP server
│   ├── jobs/
│   │   └── manager.py      # Job management system
│   └── utils.py            # Server utilities
├── scripts/
│   ├── md_simulation.py    # MD simulation runner
│   ├── tpr_analysis.py     # TPR file analyzer
│   ├── gromacs_command.py  # GROMACS command executor
│   ├── gromacs_workflow.py # Workflow manager
│   └── lib/                # Shared utilities
├── examples/
│   └── data/               # Demo data
│       ├── sample.tpr      # Sample simulation file (1000 Argon atoms)
│       ├── sample.gro      # Sample structure file
│       ├── sample.mdp      # Sample MD parameters
│       ├── sample.top      # Sample topology
│       ├── mdout.mdp       # Processed parameters
│       └── testdata.json   # Test data definitions
├── configs/                # Configuration files
├── jobs/                   # Runtime job storage
└── repo/                   # Original GROMACS repository
```

---

## Installation

### Prerequisites
- Conda or Mamba (mamba recommended for faster installation)
- Python 3.10+
- GROMACS 2025.4

### Create Environment
Please follow the environment setup procedure. An example workflow is shown below:

```bash
# Navigate to the MCP directory
cd /home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/gromacs_mcp

# Create conda environment (use mamba if available)
mamba create -p ./env python=3.10 -y
# or: conda create -p ./env python=3.10 -y

# Activate environment
mamba activate ./env
# or: conda activate ./env

# Install GROMACS
mamba install -c conda-forge gromacs=2025.4

# Install MCP dependencies
pip install fastmcp loguru --ignore-installed
```

---

## Local Usage (Scripts)

You can use the scripts directly without MCP for local processing.

### Available Scripts

| Script | Description | Example |
|--------|-------------|---------|
| `scripts/md_simulation.py` | Run MD simulations using GROMACS mdrun | See below |
| `scripts/tpr_analysis.py` | Analyze TPR files using gmx dump | See below |
| `scripts/gromacs_command.py` | Execute GROMACS command-line tools | See below |
| `scripts/gromacs_workflow.py` | Manage GROMACS simulation workflows | See below |

### Script Examples

#### MD Simulation

```bash
# Activate environment
mamba activate ./env

# Run simulation with 100 steps
python scripts/md_simulation.py \
  --input examples/data/sample.tpr \
  --output results/md_sim \
  --nsteps 100
```

**Parameters:**
- `--input, -i`: Path to TPR file (required)
- `--output, -o`: Output directory (default: current directory)
- `--nsteps, -n`: Override simulation steps (optional)
- `--config, -c`: Configuration file (optional)

#### TPR Analysis

```bash
python scripts/tpr_analysis.py \
  --input examples/data/sample.tpr \
  --output analysis.txt \
  --format text
```

**Parameters:**
- `--input, -i`: TPR file to analyze (required)
- `--output, -o`: Output file (optional, prints to console if not specified)
- `--format, -f`: Output format (text or json, default: text)

#### GROMACS Commands

```bash
# List available commands
python scripts/gromacs_command.py --list-commands

# Run specific command
python scripts/gromacs_command.py \
  --command dump \
  --input s examples/data/sample.tpr
```

#### Workflow Management

```bash
# List available workflows
python scripts/gromacs_workflow.py --list-workflows

# Run demo workflow
python scripts/gromacs_workflow.py \
  --workflow demo \
  --verbose
```

---

## MCP Server Installation

### Option 1: Using fastmcp (Recommended)

```bash
# Install MCP server for Claude Code
fastmcp install src/server.py --name gromacs-2025.4
```

### Option 2: Manual Installation for Claude Code

```bash
# Add MCP server to Claude Code
claude mcp add gromacs-2025.4 -- $(pwd)/env/bin/python $(pwd)/src/server.py

# Verify installation
claude mcp list
```

### Option 3: Configure in settings.json

Add to `~/.claude/settings.json`:

```json
{
  "mcpServers": {
    "gromacs-2025.4": {
      "command": "/home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/gromacs_mcp/env/bin/python",
      "args": ["/home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/gromacs_mcp/src/server.py"]
    }
  }
}
```

---

## Using with Claude Code

After installing the MCP server, you can use it directly in Claude Code.

### Quick Start

```bash
# Start Claude Code
claude
```

### Example Prompts

#### Tool Discovery
```
What tools are available from gromacs-2025.4?
```

#### Basic Analysis
```
Use analyze_tpr with input file @examples/data/sample.tpr
```

#### GROMACS Commands
```
Run the dump command on @examples/data/sample.tpr to see what information is in the file
```

#### Long-Running Simulations
```
Submit an MD simulation for @examples/data/sample.tpr with 10000 steps
Then check the job status
```

#### Batch Processing
```
Submit batch analysis for these TPR files:
- @examples/data/sample.tpr
- @examples/data/sample.tpr
```

### Using @ References

In Claude Code, use `@` to reference files and directories:

| Reference | Description |
|-----------|-------------|
| `@examples/data/sample.tpr` | Reference a specific TPR file |
| `@configs/md_simulation_config.json` | Reference a config file |
| `@results/` | Reference output directory |

---

## Using with Gemini CLI

### Configuration

Add to `~/.gemini/settings.json`:

```json
{
  "mcpServers": {
    "gromacs-2025.4": {
      "command": "/home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/gromacs_mcp/env/bin/python",
      "args": ["/home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/gromacs_mcp/src/server.py"]
    }
  }
}
```

### Example Prompts

```bash
# Start Gemini CLI
gemini

# Example prompts (same as Claude Code)
> What tools are available?
> Use analyze_tpr with file examples/data/sample.tpr
```

---

## Available Tools

### Quick Operations (Sync API)

These tools return results immediately (< 10 minutes):

| Tool | Description | Parameters |
|------|-------------|------------|
| `analyze_tpr` | Analyze TPR files and extract simulation parameters | `input_file`, `output_format`, `output_file` |
| `run_gromacs_command` | Execute GROMACS command-line tools | `command_name`, `input_files`, `output_files`, `parameters` |
| `run_gromacs_workflow` | Create and manage GROMACS workflows | `workflow_type`, `input_files`, `output_file`, `verbose` |

### Long-Running Tasks (Submit API)

These tools return a job_id for tracking (> 10 minutes):

| Tool | Description | Parameters |
|------|-------------|------------|
| `submit_md_simulation` | Submit MD simulation for background processing | `input_file`, `nsteps`, `output_dir`, `job_name` |
| `submit_batch_analysis` | Submit batch analysis for multiple files | `input_files`, `analysis_type`, `output_dir`, `job_name` |

### Job Management Tools

| Tool | Description |
|------|-------------|
| `get_job_status` | Check job progress and status |
| `get_job_result` | Get results when completed |
| `get_job_log` | View execution logs (supports tail) |
| `cancel_job` | Cancel running job |
| `list_jobs` | List all jobs with optional status filter |

### Information Tools

| Tool | Description |
|------|-------------|
| `get_server_info` | Get server version and capabilities |
| `list_available_commands` | List all supported commands and workflows |

---

## Examples

### Example 1: Quick TPR Analysis

**Goal:** Analyze a TPR file to understand simulation parameters

**Using Script:**
```bash
python scripts/tpr_analysis.py \
  --input examples/data/sample.tpr \
  --output results/analysis.txt
```

**Using MCP (in Claude Code):**
```
Use analyze_tpr to examine @examples/data/sample.tpr and save results to results/analysis.txt
```

**Expected Output:**
- System information: 1000 Argon atoms
- Simulation parameters: 500k steps, temperature 87K, Nose-Hoover thermostat
- Integrator: md-vv with 0.002 ps timestep

### Example 2: Run Short MD Simulation

**Goal:** Run a 100-step molecular dynamics simulation

**Using Script:**
```bash
python scripts/md_simulation.py \
  --input examples/data/sample.tpr \
  --output results/short_sim \
  --nsteps 100
```

**Using MCP (in Claude Code):**
```
Run a 100-step MD simulation on @examples/data/sample.tpr and save outputs to results/short_sim/
```

**Expected Output:**
- `confout.gro` - Final configuration
- `ener.edr` - Energy trajectory
- `md.log` - Simulation log with performance stats

### Example 3: Long Simulation (Asynchronous)

**Goal:** Submit a long simulation in the background

**Using MCP (in Claude Code):**
```
Submit an MD simulation for @examples/data/sample.tpr with 50000 steps, name it "long_test"
```

**Then monitor progress:**
```
Check the status of job long_test
Get the log output for the latest job
```

### Example 4: GROMACS Command Execution

**Goal:** Use GROMACS dump command to inspect file contents

**Using Script:**
```bash
python scripts/gromacs_command.py \
  --command dump \
  --input s examples/data/sample.tpr
```

**Using MCP (in Claude Code):**
```
Run the GROMACS dump command on @examples/data/sample.tpr to see the file structure
```

### Example 5: Workflow Demo

**Goal:** Demonstrate workflow management capabilities

**Using Script:**
```bash
python scripts/gromacs_workflow.py \
  --workflow demo \
  --verbose
```

**Using MCP (in Claude Code):**
```
Run the demo workflow with verbose output to see how workflow management works
```

---

## Demo Data

The `examples/data/` directory contains sample data for testing:

| File | Description | Use With | Atoms | Steps |
|------|-------------|----------|--------|-------|
| `sample.tpr` | Compiled simulation file | All tools | 1000 Ar | 500k |
| `sample.gro` | Structure file | Input preparation | 1000 Ar | - |
| `sample.mdp` | MD parameters | Input preparation | - | 500k |
| `sample.top` | Topology file | Input preparation | 1000 Ar | - |
| `mdout.mdp` | Processed parameters | Analysis | - | - |
| `testdata.json` | Test data definitions | Development | - | - |

### System Details
- **System**: Liquid Argon simulation
- **Atoms**: 1000 Argon atoms
- **Temperature**: 87 K
- **Integrator**: md-vv (velocity Verlet)
- **Thermostat**: Nose-Hoover
- **Barostat**: MTTK pressure coupling

---

## Configuration Files

The `configs/` directory contains configuration templates:

| Config | Description | Parameters |
|--------|-------------|------------|
| `md_simulation_config.json` | MD simulation settings | nsteps, timeout, output files |
| `tpr_analysis_config.json` | TPR analysis options | output format, parsing settings |
| `gromacs_command_config.json` | Command execution settings | timeouts, supported commands |
| `gromacs_workflow_config.json` | Workflow definitions | step dependencies, execution order |
| `default_config.json` | Global default settings | paths, logging, error handling |

### Config Example

```json
{
  "md_simulation": {
    "default_nsteps": 100,
    "timeout": 3600,
    "output_files": ["confout.gro", "ener.edr", "md.log"],
    "performance_report": true
  }
}
```

---

## Troubleshooting

### Environment Issues

**Problem:** Environment not found
```bash
# Recreate environment
mamba create -p ./env python=3.10 -y
mamba activate ./env
mamba install -c conda-forge gromacs=2025.4
pip install fastmcp loguru
```

**Problem:** GROMACS not found
```bash
# Test GROMACS installation
mamba activate ./env
gmx --version

# If not found, reinstall
mamba install -c conda-forge gromacs=2025.4
```

### MCP Issues

**Problem:** Server not found in Claude Code
```bash
# Check MCP registration
claude mcp list

# Re-add if needed
claude mcp remove gromacs-2025.4
claude mcp add gromacs-2025.4 -- $(pwd)/env/bin/python $(pwd)/src/server.py
```

**Problem:** Tools not working
```bash
# Test server directly
python src/server.py --help

# Check if tools are loadable
python -c "
import sys
sys.path.append('src')
from server import mcp
print('Available tools:', list(mcp.list_tools().keys()))
"
```

### Job Issues

**Problem:** Job stuck in pending
```bash
# Check job directory
ls -la jobs/

# Check server logs
cat gromacs_mcp.log
```

**Problem:** Job failed
```bash
# In Claude Code:
# Use get_job_log with job_id "<job_id>" and tail 50 to see error details

# Or check manually:
cat jobs/<job_id>/job.log
```

### Script Issues

**Problem:** Script execution fails
```bash
# Check environment activation
which python
which gmx

# Run with debug logging
python scripts/md_simulation.py --input examples/data/sample.tpr --output test_output --debug
```

### Performance Issues

**Problem:** Slow simulation performance
- Reduce number of steps for testing
- Check CPU usage during simulation
- Verify system has sufficient memory
- Use smaller test systems

---

## Development

### Running Tests

```bash
# Activate environment
mamba activate ./env

# Test individual scripts
python scripts/md_simulation.py --input examples/data/sample.tpr --output test_output --nsteps 10

# Test MCP server
python src/server.py
```

### Starting Dev Server

```bash
# Run MCP server in development mode
fastmcp dev src/server.py

# Test with MCP inspector
npx @anthropic/mcp-inspector src/server.py
```

### Adding New Tools

To add a new MCP tool:

1. Create script in `scripts/`
2. Add configuration in `configs/`
3. Add MCP tool wrapper in `src/server.py`
4. Update this README with examples

---

## License

Based on GROMACS 2025.4 (GPL-2.0 license)

## Credits

Based on [GROMACS 2025.4](https://github.com/gromacs/gromacs) molecular dynamics simulation package