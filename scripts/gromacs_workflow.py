#!/usr/bin/env python3
"""
Script: gromacs_workflow.py
Description: Create and manage GROMACS simulation workflows

Original Use Case: examples/use_case_5_workflow_management.py
Dependencies Removed: gmxapi (uses subprocess-based workflow patterns)

Usage:
    python scripts/gromacs_workflow.py --workflow <workflow_type> --input <files> --output <dir>

Example:
    python scripts/gromacs_workflow.py --workflow simple_md --structure input.gro --topology input.top --mdp input.mdp --output workflow_results
"""

# ==============================================================================
# Minimal Imports (only essential packages)
# ==============================================================================
import argparse
import sys
import os
import logging
import json
import time
from pathlib import Path
from typing import Union, Optional, Dict, Any, List

# Import shared utilities
try:
    from lib.gromacs_utils import (
        find_gromacs_command, run_command, ensure_output_directory,
        check_file_exists, copy_input_file
    )
except ImportError:
    # Fallback if lib not available
    import shutil
    import subprocess

    def find_gromacs_command():
        gmx_cmd = shutil.which("gmx")
        if gmx_cmd:
            return ["gmx"]
        script_dir = Path(__file__).parent.parent
        env_path = script_dir / "env"
        if env_path.exists():
            return ["mamba", "run", "-p", str(env_path), "gmx"]
        return ["gmx"]

    def run_command(cmd, working_dir=None, timeout=300):
        try:
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    cwd=working_dir, timeout=timeout)
            return result.returncode == 0, result.stdout, result.stderr
        except Exception as e:
            return False, "", str(e)

    def check_file_exists(file_path, description="File"):
        if not Path(file_path).exists():
            logging.getLogger(__name__).error(f"{description} not found: {file_path}")
            return False
        return True

    def copy_input_file(input_file, output_dir, new_name=None):
        src_path = Path(input_file)
        dst_dir = Path(output_dir)
        dst_dir.mkdir(parents=True, exist_ok=True)
        dst_path = dst_dir / (new_name or src_path.name)
        shutil.copy2(src_path, dst_path)
        return dst_path

# ==============================================================================
# Configuration
# ==============================================================================
DEFAULT_CONFIG = {
    "workflows": {
        "simple_md": {
            "description": "Simple MD workflow: grompp -> mdrun",
            "required_inputs": ["structure", "topology", "mdp"],
            "steps": ["grompp", "mdrun"],
            "default_nsteps": 1000
        },
        "analysis": {
            "description": "Trajectory analysis workflow",
            "required_inputs": ["structure", "trajectory"],
            "steps": ["rms", "rmsf", "energy"],
            "default_selections": ["Protein", "Protein", "Potential"]
        },
        "preprocessing": {
            "description": "System preprocessing workflow",
            "required_inputs": ["structure"],
            "steps": ["pdb2gmx", "editconf", "solvate"],
            "force_field": "oplsaa",
            "water_model": "spce"
        },
        "demo": {
            "description": "Demo workflow showing patterns",
            "required_inputs": [],
            "steps": ["show_patterns"],
            "patterns": ["linear", "parallel", "branched", "iterative"]
        }
    },
    "execution": {
        "timeout": 600,  # 10 minutes for workflows
        "step_timeout": 300,  # 5 minutes per step
        "continue_on_error": False,
        "log_level": "INFO"
    }
}

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# ==============================================================================
# Core Workflow Functions
# ==============================================================================
class WorkflowStep:
    """Represents a single workflow step."""

    def __init__(self, name: str, command: List[str], inputs: Dict[str, str] = None,
                 outputs: Dict[str, str] = None, depends_on: List[str] = None):
        self.name = name
        self.command = command
        self.inputs = inputs or {}
        self.outputs = outputs or {}
        self.depends_on = depends_on or []
        self.status = "pending"  # pending, running, completed, failed
        self.start_time = None
        self.end_time = None
        self.result = None

    def execute(self, working_dir: Path, timeout: int = 300) -> bool:
        """Execute this workflow step."""
        self.status = "running"
        self.start_time = time.time()

        logger.info(f"Executing step: {self.name}")
        logger.debug(f"Command: {' '.join(self.command)}")

        success, stdout, stderr = run_command(
            self.command,
            working_dir=working_dir,
            timeout=timeout
        )

        self.end_time = time.time()
        self.result = {
            "success": success,
            "stdout": stdout,
            "stderr": stderr,
            "duration": self.end_time - self.start_time
        }

        if success:
            self.status = "completed"
            logger.info(f"Step '{self.name}' completed successfully ({self.result['duration']:.1f}s)")
        else:
            self.status = "failed"
            logger.error(f"Step '{self.name}' failed")
            if stderr:
                logger.error(f"Error: {stderr}")

        return success

class WorkflowManager:
    """Manages workflow execution."""

    def __init__(self, output_dir: Path, config: Dict[str, Any] = None):
        self.output_dir = output_dir
        self.config = config or DEFAULT_CONFIG
        self.steps = []
        self.gmx_cmd = find_gromacs_command()

    def add_step(self, step: WorkflowStep):
        """Add a step to the workflow."""
        self.steps.append(step)

    def execute_workflow(self, continue_on_error: bool = False) -> Dict[str, Any]:
        """Execute all workflow steps."""
        logger.info(f"Starting workflow execution with {len(self.steps)} steps")

        results = {
            "total_steps": len(self.steps),
            "completed_steps": 0,
            "failed_steps": 0,
            "step_results": [],
            "overall_success": True
        }

        for step in self.steps:
            # Check dependencies
            if not self._check_dependencies(step):
                logger.error(f"Dependencies not met for step: {step.name}")
                step.status = "failed"
                results["failed_steps"] += 1
                results["overall_success"] = False
                if not continue_on_error:
                    break
                continue

            # Execute step
            step_timeout = self.config.get('execution', {}).get('step_timeout', 300)
            success = step.execute(self.output_dir, timeout=step_timeout)

            results["step_results"].append({
                "name": step.name,
                "status": step.status,
                "duration": step.result["duration"] if step.result else 0,
                "success": success
            })

            if success:
                results["completed_steps"] += 1
            else:
                results["failed_steps"] += 1
                results["overall_success"] = False
                if not continue_on_error:
                    break

        return results

    def _check_dependencies(self, step: WorkflowStep) -> bool:
        """Check if step dependencies are satisfied."""
        for dep_name in step.depends_on:
            dep_step = next((s for s in self.steps if s.name == dep_name), None)
            if not dep_step or dep_step.status != "completed":
                return False
        return True

def run_gromacs_workflow(
    workflow_type: str,
    input_files: Optional[Dict[str, Union[str, Path]]] = None,
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Execute a predefined GROMACS workflow.

    Args:
        workflow_type: Type of workflow to run
        input_files: Dict mapping input types to file paths
        output_file: Output directory for workflow results
        config: Configuration dict (uses DEFAULT_CONFIG if not provided)
        **kwargs: Override specific config parameters

    Returns:
        Dict containing workflow execution results

    Example:
        >>> result = run_gromacs_workflow(
        ...     "simple_md",
        ...     input_files={
        ...         "structure": "input.gro",
        ...         "topology": "input.top",
        ...         "mdp": "input.mdp"
        ...     },
        ...     output_file="workflow_output"
        ... )
    """
    # Setup
    config = {**DEFAULT_CONFIG, **(config or {}), **kwargs}

    # Extract parameters
    verbose = kwargs.get('verbose', False)
    continue_on_error = kwargs.get('continue_on_error',
                                   config['execution']['continue_on_error'])

    if verbose:
        logger.setLevel(logging.DEBUG)

    logger.info(f"Running GROMACS workflow: {workflow_type}")

    # Get workflow configuration
    workflow_config = config['workflows'].get(workflow_type)
    if not workflow_config:
        available = ", ".join(config['workflows'].keys())
        raise ValueError(f"Unknown workflow type: {workflow_type}. Available: {available}")

    # Setup output directory
    if output_file:
        output_dir = ensure_output_directory(output_file)
    else:
        output_dir = ensure_output_directory(f"{workflow_type}_output")

    logger.info(f"Workflow output directory: {output_dir}")

    try:
        # Create workflow manager
        workflow_manager = WorkflowManager(output_dir, config)

        # Build workflow based on type
        if workflow_type == "simple_md":
            _build_simple_md_workflow(workflow_manager, input_files, workflow_config, **kwargs)
        elif workflow_type == "analysis":
            _build_analysis_workflow(workflow_manager, input_files, workflow_config, **kwargs)
        elif workflow_type == "preprocessing":
            _build_preprocessing_workflow(workflow_manager, input_files, workflow_config, **kwargs)
        elif workflow_type == "demo":
            _build_demo_workflow(workflow_manager, workflow_config, **kwargs)
        else:
            raise ValueError(f"Workflow type '{workflow_type}' not implemented")

        # Execute workflow
        logger.info("Starting workflow execution...")
        results = workflow_manager.execute_workflow(continue_on_error=continue_on_error)

        # Add metadata
        results.update({
            "workflow_type": workflow_type,
            "output_dir": str(output_dir),
            "input_files": input_files or {},
            "config": workflow_config
        })

        if results["overall_success"]:
            logger.info("Workflow completed successfully")
        else:
            logger.warning(f"Workflow completed with {results['failed_steps']} failed steps")

        return results

    except Exception as e:
        logger.error(f"Error during workflow execution: {e}")
        return {
            "workflow_type": workflow_type,
            "overall_success": False,
            "error": str(e),
            "output_dir": str(output_dir) if 'output_dir' in locals() else None
        }

def _build_simple_md_workflow(manager: WorkflowManager, input_files: Dict[str, str],
                             config: Dict[str, Any], **kwargs):
    """Build simple MD workflow: grompp -> mdrun."""
    gmx_cmd = manager.gmx_cmd

    # Validate required inputs
    required = ["structure", "topology", "mdp"]
    for req in required:
        if not input_files or req not in input_files:
            raise ValueError(f"Required input '{req}' not provided")
        if not check_file_exists(input_files[req], f"{req} file"):
            raise FileNotFoundError(f"{req} file not found: {input_files[req]}")

    # Copy input files to workflow directory
    for input_type, filepath in input_files.items():
        copy_input_file(filepath, manager.output_dir)

    # Step 1: grompp
    grompp_cmd = gmx_cmd + [
        "grompp",
        "-f", Path(input_files["mdp"]).name,
        "-c", Path(input_files["structure"]).name,
        "-p", Path(input_files["topology"]).name,
        "-o", "topol.tpr"
    ]
    manager.add_step(WorkflowStep("grompp", grompp_cmd))

    # Step 2: mdrun
    nsteps = kwargs.get('nsteps', config.get('default_nsteps', 1000))
    mdrun_cmd = gmx_cmd + [
        "mdrun",
        "-s", "topol.tpr",
        "-nsteps", str(nsteps),
        "-c", "confout.gro",
        "-e", "ener.edr",
        "-g", "md.log"
    ]
    manager.add_step(WorkflowStep("mdrun", mdrun_cmd, depends_on=["grompp"]))

def _build_analysis_workflow(manager: WorkflowManager, input_files: Dict[str, str],
                           config: Dict[str, Any], **kwargs):
    """Build trajectory analysis workflow."""
    gmx_cmd = manager.gmx_cmd

    # Validate required inputs
    required = ["structure", "trajectory"]
    for req in required:
        if not input_files or req not in input_files:
            raise ValueError(f"Required input '{req}' not provided")
        if not check_file_exists(input_files[req], f"{req} file"):
            raise FileNotFoundError(f"{req} file not found: {input_files[req]}")

    # Copy input files
    struct_file = copy_input_file(input_files["structure"], manager.output_dir)
    traj_file = copy_input_file(input_files["trajectory"], manager.output_dir)

    # Add analysis steps
    analyses = ["rms", "rmsf"]
    for analysis in analyses:
        cmd = gmx_cmd + [
            analysis,
            "-s", struct_file.name,
            "-f", traj_file.name,
            "-o", f"{analysis}.xvg"
        ]
        # Add default selections via echo
        full_cmd = ["bash", "-c", f"echo '4 4' | {' '.join(cmd)}"]
        manager.add_step(WorkflowStep(analysis, full_cmd))

def _build_preprocessing_workflow(manager: WorkflowManager, input_files: Dict[str, str],
                                config: Dict[str, Any], **kwargs):
    """Build system preprocessing workflow."""
    logger.info("Demo: Preprocessing workflow would include:")
    logger.info("  1. pdb2gmx - Add hydrogens and generate topology")
    logger.info("  2. editconf - Define simulation box")
    logger.info("  3. solvate - Add water molecules")
    logger.info("  4. genion - Add ions for neutralization")

    # For demo, just create a simple step
    manager.add_step(WorkflowStep(
        "demo_preprocessing",
        ["echo", "Preprocessing workflow demo completed"]
    ))

def _build_demo_workflow(manager: WorkflowManager, config: Dict[str, Any], **kwargs):
    """Build demo workflow showing different patterns."""
    patterns = config.get('patterns', ["linear", "parallel", "branched"])

    logger.info("Demo workflow patterns:")
    for i, pattern in enumerate(patterns):
        logger.info(f"  {i+1}. {pattern.capitalize()} workflow pattern")

        # Add demo step
        manager.add_step(WorkflowStep(
            f"demo_{pattern}",
            ["echo", f"Demonstrating {pattern} workflow pattern"]
        ))

# ==============================================================================
# CLI Interface
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--workflow', '-w',
        choices=['simple_md', 'analysis', 'preprocessing', 'demo'],
        help='Workflow type to execute'
    )

    parser.add_argument(
        '--structure', '-s',
        help='Structure file (.gro, .pdb)'
    )

    parser.add_argument(
        '--topology', '-p',
        help='Topology file (.top)'
    )

    parser.add_argument(
        '--mdp', '-f',
        help='MD parameter file (.mdp)'
    )

    parser.add_argument(
        '--trajectory', '-t',
        help='Trajectory file (.xtc, .trr)'
    )

    parser.add_argument(
        '--output', '-o',
        help='Output directory (default: <workflow_type>_output)'
    )

    parser.add_argument(
        '--nsteps', '-n',
        type=int,
        help='Number of simulation steps for MD workflows'
    )

    parser.add_argument(
        '--config', '-c',
        help='Config file (JSON format)'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )

    parser.add_argument(
        '--continue-on-error',
        action='store_true',
        help='Continue workflow execution even if steps fail'
    )

    parser.add_argument(
        '--list-workflows',
        action='store_true',
        help='List available workflow types'
    )

    args = parser.parse_args()

    # Handle list workflows
    if args.list_workflows:
        print("Available workflow types:")
        for workflow_type, info in DEFAULT_CONFIG['workflows'].items():
            print(f"  {workflow_type:<15} - {info.get('description', 'No description')}")
            required = info.get('required_inputs', [])
            if required:
                print(f"    Required inputs: {', '.join(required)}")
        return 0

    # Validate required workflow argument
    if not args.workflow:
        parser.error("--workflow is required unless using --list-workflows")

    # Load config if provided
    config = None
    if args.config:
        with open(args.config) as f:
            config = json.load(f)

    # Prepare input files
    input_files = {}
    if args.structure:
        input_files['structure'] = args.structure
    if args.topology:
        input_files['topology'] = args.topology
    if args.mdp:
        input_files['mdp'] = args.mdp
    if args.trajectory:
        input_files['trajectory'] = args.trajectory

    # Prepare kwargs
    kwargs = {}
    if args.verbose:
        kwargs['verbose'] = args.verbose
    if args.nsteps:
        kwargs['nsteps'] = args.nsteps
    if args.continue_on_error:
        kwargs['continue_on_error'] = args.continue_on_error

    try:
        # Run workflow
        result = run_gromacs_workflow(
            workflow_type=args.workflow,
            input_files=input_files,
            output_file=args.output,
            config=config,
            **kwargs
        )

        if result.get('overall_success', False):
            print(f"✅ Success: Workflow '{args.workflow}' completed")
            print(f"   Output directory: {result.get('output_dir', 'Unknown')}")
            print(f"   Completed steps: {result.get('completed_steps', 0)}/{result.get('total_steps', 0)}")
            return 0
        else:
            print(f"❌ Failed: {result.get('error', 'Unknown error')}")
            if 'failed_steps' in result:
                print(f"   Failed steps: {result['failed_steps']}/{result.get('total_steps', 0)}")
            return 1

    except Exception as e:
        print(f"❌ Error: {e}")
        return 1

if __name__ == '__main__':
    sys.exit(main())