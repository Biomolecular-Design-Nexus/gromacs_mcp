#!/usr/bin/env python3
"""MCP Server for GROMACS 2025.4

Provides both synchronous and asynchronous (submit) APIs for GROMACS tools.
"""

from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import sys

# Setup paths
SCRIPT_DIR = Path(__file__).parent.resolve()
MCP_ROOT = SCRIPT_DIR.parent
SCRIPTS_DIR = MCP_ROOT / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(SCRIPTS_DIR))

from jobs.manager import job_manager
from utils import validate_input_file, setup_output_directory, format_tool_result, configure_logging
from loguru import logger

# Configure logging
configure_logging()

# Create MCP server
mcp = FastMCP("gromacs-2025.4")

# ==============================================================================
# Job Management Tools (for async operations)
# ==============================================================================

@mcp.tool()
def get_job_status(job_id: str) -> dict:
    """
    Get the status of a submitted job.

    Args:
        job_id: The job ID returned from a submit_* function

    Returns:
        Dictionary with job status, timestamps, and any errors
    """
    return job_manager.get_job_status(job_id)

@mcp.tool()
def get_job_result(job_id: str) -> dict:
    """
    Get the results of a completed job.

    Args:
        job_id: The job ID of a completed job

    Returns:
        Dictionary with the job results or error if not completed
    """
    return job_manager.get_job_result(job_id)

@mcp.tool()
def get_job_log(job_id: str, tail: int = 50) -> dict:
    """
    Get log output from a running or completed job.

    Args:
        job_id: The job ID to get logs for
        tail: Number of lines from end (default: 50, use 0 for all)

    Returns:
        Dictionary with log lines and total line count
    """
    return job_manager.get_job_log(job_id, tail)

@mcp.tool()
def cancel_job(job_id: str) -> dict:
    """
    Cancel a running job.

    Args:
        job_id: The job ID to cancel

    Returns:
        Success or error message
    """
    return job_manager.cancel_job(job_id)

@mcp.tool()
def list_jobs(status: Optional[str] = None) -> dict:
    """
    List all submitted jobs.

    Args:
        status: Filter by status (pending, running, completed, failed, cancelled)

    Returns:
        List of jobs with their status
    """
    return job_manager.list_jobs(status)

# ==============================================================================
# Synchronous Tools (for fast operations < 10 min)
# ==============================================================================

@mcp.tool()
def analyze_tpr(
    input_file: str,
    output_format: str = "text",
    output_file: Optional[str] = None
) -> dict:
    """
    Analyze GROMACS TPR files and extract simulation parameters.

    Fast operation that completes in seconds. Use for examining TPR file contents,
    simulation parameters, and system information.

    Args:
        input_file: Path to TPR file to analyze
        output_format: Output format ('text' or 'json')
        output_file: Optional path to save analysis results

    Returns:
        Dictionary with parsed TPR information and analysis results
    """
    from tpr_analysis import run_tpr_analysis

    try:
        # Validate input
        validation = validate_input_file(input_file, ".tpr")
        if not validation["valid"]:
            return format_tool_result(False, error_message=validation["error"])

        # Setup output
        output_setup = setup_output_directory(output_file, "tpr_analysis")
        if not output_setup["success"]:
            return format_tool_result(False, error_message=output_setup["error"])

        # Run analysis
        result = run_tpr_analysis(
            input_file=input_file,
            output_file=output_setup.get("output_file"),
            format=output_format
        )

        return format_tool_result(
            success=result.get("result", False),
            data=result,
            output_files=result.get("output_files", [])
        )

    except FileNotFoundError as e:
        return format_tool_result(False, error_message=f"File not found: {e}")
    except ValueError as e:
        return format_tool_result(False, error_message=f"Invalid input: {e}")
    except Exception as e:
        logger.error(f"TPR analysis failed: {e}")
        return format_tool_result(False, error_message=str(e))

@mcp.tool()
def run_gromacs_command(
    command_name: str,
    input_files: Optional[dict] = None,
    output_files: Optional[dict] = None,
    parameters: Optional[List[str]] = None
) -> dict:
    """
    Execute GROMACS command-line tools with flexible parameter handling.

    Fast operation for running GROMACS analysis commands like rms, rmsf, energy.
    Supports all common GROMACS commands with automatic parameter handling.

    Args:
        command_name: GROMACS command name (e.g., 'rms', 'energy', 'dump')
        input_files: Dict of input files with flags (e.g., {'s': 'topology.tpr'})
        output_files: Dict of output files with flags (e.g., {'o': 'rmsd.xvg'})
        parameters: List of additional parameters

    Returns:
        Dictionary with command execution results and output files
    """
    from gromacs_command import run_gromacs_command

    try:
        # Run GROMACS command
        result = run_gromacs_command(
            command_name=command_name,
            input_files=input_files or {},
            output_files=output_files or {},
            parameters=parameters or []
        )

        return format_tool_result(
            success=result.get("result", False),
            data=result,
            output_files=result.get("output_files", [])
        )

    except Exception as e:
        logger.error(f"GROMACS command failed: {e}")
        return format_tool_result(False, error_message=str(e))

@mcp.tool()
def run_gromacs_workflow(
    workflow_type: str,
    input_files: Optional[dict] = None,
    output_file: Optional[str] = None,
    verbose: bool = False
) -> dict:
    """
    Create and manage GROMACS simulation workflows.

    Fast operation for running predefined workflows like 'simple_md', 'analysis',
    or 'demo'. Manages workflow dependencies and execution order automatically.

    Args:
        workflow_type: Type of workflow ('simple_md', 'analysis', 'demo')
        input_files: Dict of input files for the workflow
        output_file: Path to save workflow results
        verbose: Enable verbose output

    Returns:
        Dictionary with workflow execution results and generated files
    """
    from gromacs_workflow import run_gromacs_workflow

    try:
        # Run workflow
        result = run_gromacs_workflow(
            workflow_type=workflow_type,
            input_files=input_files,
            output_file=output_file,
            verbose=verbose
        )

        return format_tool_result(
            success=result.get("result", False),
            data=result,
            output_files=result.get("output_files", [])
        )

    except Exception as e:
        logger.error(f"GROMACS workflow failed: {e}")
        return format_tool_result(False, error_message=str(e))

# ==============================================================================
# Submit Tools (for long-running operations > 10 min)
# ==============================================================================

@mcp.tool()
def submit_md_simulation(
    input_file: str,
    nsteps: Optional[int] = None,
    output_dir: Optional[str] = None,
    job_name: Optional[str] = None
) -> dict:
    """
    Submit molecular dynamics simulation for background processing.

    This operation can take significant time (minutes to hours) depending on
    simulation length. Returns a job_id for tracking progress.

    Args:
        input_file: Path to TPR file for simulation
        nsteps: Number of simulation steps (overrides TPR setting)
        output_dir: Directory to save simulation outputs
        job_name: Optional name for tracking this simulation

    Returns:
        Dictionary with job_id. Use:
        - get_job_status(job_id) to check progress
        - get_job_result(job_id) to get results when completed
        - get_job_log(job_id) to see simulation logs
    """
    # Validate input first
    validation = validate_input_file(input_file, ".tpr")
    if not validation["valid"]:
        return format_tool_result(False, error_message=validation["error"])

    script_path = str(SCRIPTS_DIR / "md_simulation.py")

    return job_manager.submit_job(
        script_path=script_path,
        args={
            "input": input_file,
            "nsteps": nsteps,
            "output": output_dir
        },
        job_name=job_name or f"md_sim_{Path(input_file).stem}"
    )

@mcp.tool()
def submit_batch_analysis(
    input_files: List[str],
    analysis_type: str = "tpr",
    output_dir: Optional[str] = None,
    job_name: Optional[str] = None
) -> dict:
    """
    Submit batch analysis for multiple input files.

    Processes multiple TPR files or trajectories in a single background job.
    Suitable for analyzing many simulation files at once.

    Args:
        input_files: List of file paths to analyze
        analysis_type: Type of analysis ('tpr', 'trajectory')
        output_dir: Directory to save all analysis results
        job_name: Optional name for the batch job

    Returns:
        Dictionary with job_id for tracking the batch analysis
    """
    # Validate all input files first
    for input_file in input_files:
        validation = validate_input_file(input_file)
        if not validation["valid"]:
            return format_tool_result(False, error_message=f"Invalid file {input_file}: {validation['error']}")

    # Choose appropriate script based on analysis type
    if analysis_type == "tpr":
        script_path = str(SCRIPTS_DIR / "tpr_analysis.py")
    else:
        return format_tool_result(False, error_message=f"Unsupported analysis type: {analysis_type}")

    # For batch processing, we'll create a simple wrapper
    # In a full implementation, you'd want a dedicated batch script
    return job_manager.submit_job(
        script_path=script_path,
        args={
            "inputs": ",".join(input_files),  # Pass as comma-separated list
            "output": output_dir,
            "batch": True
        },
        job_name=job_name or f"batch_{analysis_type}_{len(input_files)}_files"
    )

# ==============================================================================
# Information and Help Tools
# ==============================================================================

@mcp.tool()
def list_available_commands() -> dict:
    """
    List all available GROMACS commands supported by the server.

    Returns:
        Dictionary with categorized lists of available commands and workflows
    """
    try:
        from gromacs_command import run_gromacs_command

        # Get available commands (this would be enhanced to read from config)
        sync_commands = [
            "analyze_tpr",
            "run_gromacs_command",
            "run_gromacs_workflow"
        ]

        submit_commands = [
            "submit_md_simulation",
            "submit_batch_analysis"
        ]

        job_commands = [
            "get_job_status",
            "get_job_result",
            "get_job_log",
            "cancel_job",
            "list_jobs"
        ]

        gromacs_commands = [
            "rms", "rmsf", "energy", "sasa", "hbond",
            "trjconv", "dump", "check"
        ]

        workflow_types = [
            "simple_md", "analysis", "preprocessing", "demo"
        ]

        return format_tool_result(True, {
            "sync_tools": sync_commands,
            "submit_tools": submit_commands,
            "job_management": job_commands,
            "gromacs_commands": gromacs_commands,
            "workflow_types": workflow_types,
            "example_data": str(MCP_ROOT / "examples" / "data")
        })

    except Exception as e:
        return format_tool_result(False, error_message=str(e))

@mcp.tool()
def get_server_info() -> dict:
    """
    Get information about the GROMACS MCP server.

    Returns:
        Dictionary with server version, capabilities, and example usage
    """
    return format_tool_result(True, {
        "server_name": "gromacs-2025.4",
        "version": "1.0.0",
        "description": "MCP server providing GROMACS molecular dynamics tools",
        "capabilities": {
            "sync_operations": "Fast analysis and command execution (<10 min)",
            "async_operations": "Long-running simulations with job tracking",
            "job_management": "Full job lifecycle management",
            "batch_processing": "Multiple file processing in single jobs"
        },
        "example_usage": {
            "sync": "analyze_tpr with input_file 'examples/data/sample.tpr'",
            "async": "submit_md_simulation with input_file 'examples/data/sample.tpr'"
        },
        "job_workflow": [
            "1. Submit: Use submit_* tools to start long jobs",
            "2. Monitor: Use get_job_status(job_id) to check progress",
            "3. Retrieve: Use get_job_result(job_id) when completed",
            "4. Debug: Use get_job_log(job_id) for troubleshooting"
        ]
    })

# ==============================================================================
# Entry Point
# ==============================================================================

if __name__ == "__main__":
    mcp.run()