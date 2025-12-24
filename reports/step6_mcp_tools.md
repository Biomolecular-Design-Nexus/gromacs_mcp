# Step 6: MCP Tools Documentation

## Server Information
- **Server Name**: gromacs-2025.4
- **Version**: 1.0.0
- **Created Date**: 2024-12-20
- **Server Path**: `src/server.py`
- **Package Manager**: mamba (preferred) / conda
- **Dependencies**: fastmcp, loguru, GROMACS 2025.4

## Architecture Overview

The GROMACS MCP server provides both synchronous and asynchronous APIs:

```
src/
├── server.py              # Main MCP server entry point
├── tools/
│   ├── __init__.py        # Tools package initialization
│   ├── sync_tools.py      # Future: Synchronous tools (if needed)
│   └── async_tools.py     # Future: Async tools (if needed)
├── jobs/
│   ├── __init__.py        # Jobs package initialization
│   └── manager.py         # Job queue management and persistence
└── utils.py               # Shared utilities and validation

jobs/                      # Runtime job storage
├── <job_id>/
│   ├── metadata.json      # Job status and metadata
│   ├── output.json        # Job results (when completed)
│   └── job.log           # Execution logs
```

## API Classification

Based on runtime analysis from Step 5:

| Script | API Type | Runtime | Reason |
|--------|----------|---------|--------|
| `tpr_analysis.py` | **Sync** | ~3.1s | Fast analysis operation |
| `gromacs_command.py` | **Sync** | 1-3s | Quick command execution |
| `gromacs_workflow.py` | **Sync** | ~1.5s | Demo workflows are fast |
| `md_simulation.py` | **Submit** | >15s | Variable length, can be hours |

## Job Management Tools

| Tool | Description | Returns |
|------|-------------|---------|
| `get_job_status` | Check job progress and status | Job metadata with status |
| `get_job_result` | Get completed job results | Job output data |
| `get_job_log` | View job execution logs | Log lines (tail support) |
| `cancel_job` | Cancel running job | Success/error message |
| `list_jobs` | List all jobs (with optional filter) | Array of job summaries |

### Job Status Lifecycle

```
PENDING → RUNNING → COMPLETED
                 → FAILED
                 → CANCELLED
```

### Job Management Examples

```bash
# Submit a long-running job
submit_md_simulation with input_file "examples/data/sample.tpr" nsteps 1000000
→ Returns: {"job_id": "abc12345", "status": "submitted"}

# Check job progress
get_job_status with job_id "abc12345"
→ Returns: {"status": "running", "started_at": "2024-12-20T10:30:00"}

# Get results when completed
get_job_result with job_id "abc12345"
→ Returns: {"status": "success", "result": {...}, "output_files": [...]}

# View execution logs
get_job_log with job_id "abc12345" tail 20
→ Returns: {"log_lines": [...], "total_lines": 150}
```

---

## Synchronous Tools (Fast Operations < 10 min)

### analyze_tpr
- **Description**: Analyze GROMACS TPR files and extract simulation parameters
- **Source Script**: `scripts/tpr_analysis.py`
- **Estimated Runtime**: ~3 seconds

**Parameters:**
| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| input_file | str | Yes | - | Path to TPR file to analyze |
| output_format | str | No | "text" | Output format ('text' or 'json') |
| output_file | str | No | None | Optional path to save results |

**Returns:**
- Parsed TPR information (atoms, steps, simulation parameters)
- Analysis results in specified format
- List of generated output files

**Example:**
```bash
analyze_tpr with input_file "examples/data/sample.tpr" output_format "json"
```

---

### run_gromacs_command
- **Description**: Execute GROMACS command-line tools with flexible parameter handling
- **Source Script**: `scripts/gromacs_command.py`
- **Estimated Runtime**: 1-3 seconds typical

**Parameters:**
| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| command_name | str | Yes | - | GROMACS command name (e.g., 'rms', 'energy') |
| input_files | dict | No | {} | Dict of input files with flags |
| output_files | dict | No | {} | Dict of output files with flags |
| parameters | list | No | [] | Additional parameters |

**Supported Commands:**
- `rms` - RMSD analysis
- `rmsf` - RMSF analysis
- `energy` - Energy analysis
- `sasa` - Surface area analysis
- `hbond` - Hydrogen bond analysis
- `trjconv` - Trajectory conversion
- `dump` - Dump file contents
- `check` - Check file integrity

**Example:**
```bash
run_gromacs_command with command_name "dump" input_files {"s": "examples/data/sample.tpr"}
```

---

### run_gromacs_workflow
- **Description**: Create and manage GROMACS simulation workflows
- **Source Script**: `scripts/gromacs_workflow.py`
- **Estimated Runtime**: ~1.5 seconds for demo workflows

**Parameters:**
| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| workflow_type | str | Yes | - | Type of workflow to run |
| input_files | dict | No | {} | Dict of input files for workflow |
| output_file | str | No | None | Path to save workflow results |
| verbose | bool | No | false | Enable verbose output |

**Workflow Types:**
- `simple_md` - grompp → mdrun workflow
- `analysis` - Trajectory analysis workflow
- `preprocessing` - System preparation workflow
- `demo` - Demonstration workflow patterns

**Example:**
```bash
run_gromacs_workflow with workflow_type "demo" verbose true
```

---

## Submit Tools (Long Operations > 10 min)

### submit_md_simulation
- **Description**: Submit molecular dynamics simulation for background processing
- **Source Script**: `scripts/md_simulation.py`
- **Estimated Runtime**: Minutes to hours (depends on simulation length)
- **Supports Batch**: No (single simulation per job)

**Parameters:**
| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| input_file | str | Yes | - | Path to TPR file for simulation |
| nsteps | int | No | None | Number of simulation steps (overrides TPR) |
| output_dir | str | No | None | Directory to save simulation outputs |
| job_name | str | No | auto | Custom job name for tracking |

**Generated Files:**
- `confout.gro` - Final configuration
- `ener.edr` - Energy trajectory
- `md.log` - Simulation log
- `traj.xtc` - Trajectory (if configured)

**Example:**
```bash
submit_md_simulation with input_file "examples/data/sample.tpr" nsteps 100000 job_name "test_simulation"
```

**Workflow:**
1. **Submit**: Returns job_id immediately
2. **Monitor**: Use `get_job_status(job_id)` to track progress
3. **Retrieve**: Use `get_job_result(job_id)` when status is "completed"
4. **Debug**: Use `get_job_log(job_id)` for troubleshooting

---

### submit_batch_analysis
- **Description**: Submit batch analysis for multiple input files
- **Source Script**: `scripts/tpr_analysis.py` (batch mode)
- **Estimated Runtime**: Depends on number of files and analysis type
- **Supports Batch**: ✅ Yes

**Parameters:**
| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| input_files | list | Yes | - | List of file paths to analyze |
| analysis_type | str | No | "tpr" | Type of analysis to perform |
| output_dir | str | No | None | Directory to save all results |
| job_name | str | No | auto | Custom job name |

**Analysis Types:**
- `tpr` - TPR file analysis

**Example:**
```bash
submit_batch_analysis with input_files ["file1.tpr", "file2.tpr", "file3.tpr"] analysis_type "tpr" output_dir "batch_results"
```

---

## Information and Help Tools

### get_server_info
- **Description**: Get information about the GROMACS MCP server
- **Returns**: Server version, capabilities, and example usage

### list_available_commands
- **Description**: List all available GROMACS commands and workflows
- **Returns**: Categorized lists of tools, commands, and workflows

---

## Error Handling

All tools return structured error responses:

```json
{
  "status": "error",
  "error": "Detailed error message explaining what went wrong"
}
```

### Common Error Types

1. **File Not Found**
   ```json
   {"status": "error", "error": "File not found: /path/to/file.tpr"}
   ```

2. **Invalid File Format**
   ```json
   {"status": "error", "error": "Expected .tpr file, got .gro"}
   ```

3. **GROMACS Command Failed**
   ```json
   {"status": "error", "error": "Process exited with code 1"}
   ```

4. **Job Not Found**
   ```json
   {"status": "error", "error": "Job abc12345 not found"}
   ```

---

## Workflow Examples

### Quick Analysis (Synchronous)
For fast operations that complete in seconds:

```bash
1. List available tools:
   list_available_commands

2. Analyze a TPR file:
   analyze_tpr with input_file "examples/data/sample.tpr"
   → Returns results immediately

3. Run GROMACS command:
   run_gromacs_command with command_name "dump" input_files {"s": "examples/data/sample.tpr"}
   → Returns command output immediately
```

### Long-Running Simulation (Asynchronous)
For operations that take more than 10 minutes:

```bash
1. Submit simulation:
   submit_md_simulation with input_file "examples/data/sample.tpr" nsteps 1000000
   → Returns: {"job_id": "abc12345", "status": "submitted"}

2. Check status periodically:
   get_job_status with job_id "abc12345"
   → Returns: {"status": "running", "started_at": "..."}

3. View logs if needed:
   get_job_log with job_id "abc12345" tail 30
   → Returns: {"log_lines": [...recent log entries...]}

4. Get results when completed:
   get_job_status with job_id "abc12345"
   → Returns: {"status": "completed", "completed_at": "..."}

   get_job_result with job_id "abc12345"
   → Returns: {"status": "success", "result": {...}, "output_files": [...]}
```

### Batch Processing
For processing multiple files efficiently:

```bash
1. Submit batch job:
   submit_batch_analysis with input_files ["file1.tpr", "file2.tpr", "file3.tpr"]
   → Returns: {"job_id": "batch123", "status": "submitted"}

2. Monitor batch progress:
   get_job_status with job_id "batch123"
   → Returns progress for entire batch

3. Get batch results:
   get_job_result with job_id "batch123"
   → Returns: {"status": "success", "result": {"file1.tpr": {...}, "file2.tpr": {...}, ...}}
```

---

## Performance Characteristics

### Synchronous Tools
- **Startup Time**: < 1 second
- **TPR Analysis**: ~3 seconds for typical files
- **GROMACS Commands**: 1-3 seconds typical
- **Memory Usage**: Minimal (standard library only)

### Asynchronous Tools
- **Job Submission**: < 1 second
- **MD Simulation**: 15.2s for 100 steps, scales with nsteps
- **Batch Analysis**: Linear scaling with number of files
- **Job Persistence**: Survives server restarts

### Resource Management
- **Concurrent Jobs**: Unlimited (managed by threading)
- **Job Storage**: Persistent on disk in `jobs/` directory
- **Log Rotation**: Manual cleanup required
- **Timeout Handling**: Configurable per operation

---

## Configuration and Customization

### Environment Setup
```bash
# Preferred: Use mamba for environment management
mamba activate ./env

# Fallback: Use conda
conda activate ./env

# Install dependencies
pip install fastmcp loguru
```

### Server Configuration
The server automatically detects:
- GROMACS installation via `gmx` command
- Environment path for mamba/conda execution
- Available example data in `examples/data/`

### Custom Configurations
Each script uses configuration files in `configs/`:
- `md_simulation_config.json` - MD simulation settings
- `tpr_analysis_config.json` - Analysis output options
- `gromacs_command_config.json` - Command execution settings
- `gromacs_workflow_config.json` - Workflow definitions

---

## Integration Examples

### With Claude Desktop
Add to Claude Desktop configuration:

```json
{
  "mcpServers": {
    "gromacs-2025.4": {
      "command": "mamba",
      "args": ["run", "-p", "./env", "python", "src/server.py"]
    }
  }
}
```

### With FastMCP CLI
```bash
# Install server
fastmcp install claude-code src/server.py

# Run in development mode
fastmcp dev src/server.py

# Test with MCP inspector
npx @anthropic/mcp-inspector src/server.py
```

---

## Success Criteria Assessment

- [x] MCP server created at `src/server.py`
- [x] Job manager implemented for async operations
- [x] Sync tools created for fast operations (<10 min)
- [x] Submit tools created for long-running operations (>10 min)
- [x] Batch processing support for applicable tools
- [x] Job management tools working (status, result, log, cancel, list)
- [x] All tools have clear descriptions for LLM use
- [x] Error handling returns structured responses
- [x] Server starts without errors: `fastmcp dev src/server.py`
- [x] Tool documentation complete in this file

---

## Tool Summary

**Total Tools**: 12
- **Job Management**: 5 tools (get_job_status, get_job_result, get_job_log, cancel_job, list_jobs)
- **Synchronous**: 3 tools (analyze_tpr, run_gromacs_command, run_gromacs_workflow)
- **Submit**: 2 tools (submit_md_simulation, submit_batch_analysis)
- **Information**: 2 tools (get_server_info, list_available_commands)

**API Coverage**:
- ✅ All Step 5 scripts wrapped as MCP tools
- ✅ Consistent parameter naming and validation
- ✅ Structured error responses for reliable LLM interaction
- ✅ Comprehensive job lifecycle management
- ✅ Batch processing capabilities
- ✅ Example data and usage patterns provided

This implementation successfully converts all GROMACS scripts from Step 5 into a fully functional MCP server with both synchronous and asynchronous capabilities, providing a robust foundation for LLM-driven molecular dynamics workflows.