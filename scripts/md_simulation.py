#!/usr/bin/env python3
"""
Script: md_simulation.py
Description: Run basic molecular dynamics simulation using GROMACS

Original Use Case: examples/use_case_1_basic_md_simulation.py
Dependencies Removed: None (uses standard library + GROMACS command)

Usage:
    python scripts/md_simulation.py --input <tpr_file> --output <output_dir>

Example:
    python scripts/md_simulation.py --input examples/data/sample.tpr --output results/md_sim
"""

# ==============================================================================
# Minimal Imports (only essential packages)
# ==============================================================================
import argparse
import sys
import os
import logging
from pathlib import Path
from typing import Union, Optional, Dict, Any

# Import shared utilities
try:
    from lib.gromacs_utils import (
        find_gromacs_command, run_command, ensure_output_directory,
        check_file_exists, get_file_size, format_file_info, copy_input_file
    )
except ImportError:
    # Fallback if lib not available - inline essential functions
    import shutil
    import subprocess

    def find_gromacs_command():
        gmx_cmd = shutil.which("gmx")
        if gmx_cmd:
            return ["gmx"]
        # Try environment fallback
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

    def ensure_output_directory(output_dir):
        path = Path(output_dir)
        path.mkdir(parents=True, exist_ok=True)
        return path

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
    "simulation": {
        "default_nsteps": 100,  # Short demo simulation
        "timeout": 300,         # 5 minute timeout
        "tpr_name": "topol.tpr"
    },
    "output_files": {
        "trajectory": "traj.trr",
        "configuration": "confout.gro",
        "energy": "ener.edr",
        "log": "md.log"
    }
}

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# ==============================================================================
# Core Function
# ==============================================================================
def run_md_simulation(
    input_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Run a basic molecular dynamics simulation using GROMACS mdrun.

    Args:
        input_file: Path to TPR input file
        output_file: Path to save main output (directory will be created)
        config: Configuration dict (uses DEFAULT_CONFIG if not provided)
        **kwargs: Override specific config parameters
            - nsteps: Number of simulation steps
            - verbose: Enable verbose logging
            - timeout: Command timeout in seconds

    Returns:
        Dict containing:
            - result: Simulation success status
            - output_files: List of generated output files
            - metadata: Execution metadata

    Example:
        >>> result = run_md_simulation("input.tpr", "output/")
        >>> print(result['output_files'])
    """
    # Setup
    input_path = Path(input_file)
    config = {**DEFAULT_CONFIG, **(config or {}), **kwargs}

    # Extract parameters
    nsteps = kwargs.get('nsteps', config['simulation']['default_nsteps'])
    verbose = kwargs.get('verbose', False)
    timeout = kwargs.get('timeout', config['simulation']['timeout'])

    if verbose:
        logger.setLevel(logging.DEBUG)

    logger.info(f"Running MD simulation from TPR file: {input_path}")

    # Validate input
    if not check_file_exists(input_path, "TPR file"):
        raise FileNotFoundError(f"TPR file not found: {input_path}")

    # Setup output directory
    if output_file:
        output_dir = ensure_output_directory(output_file)
    else:
        output_dir = ensure_output_directory("md_output")

    logger.info(f"Output directory: {output_dir}")

    try:
        # Find GROMACS command
        gmx_cmd = find_gromacs_command()
        logger.debug(f"Using GROMACS command: {' '.join(gmx_cmd)}")

        # Copy TPR file to output directory
        tpr_output = copy_input_file(
            input_path,
            output_dir,
            config['simulation']['tpr_name']
        )

        # Prepare mdrun command
        mdrun_cmd = gmx_cmd + ["mdrun"]
        mdrun_cmd.extend(["-s", config['simulation']['tpr_name']])
        mdrun_cmd.extend(["-o", config['output_files']['trajectory']])
        mdrun_cmd.extend(["-c", config['output_files']['configuration']])
        mdrun_cmd.extend(["-e", config['output_files']['energy']])
        mdrun_cmd.extend(["-g", config['output_files']['log']])

        # Add number of steps
        mdrun_cmd.extend(["-nsteps", str(nsteps)])
        logger.info(f"Running simulation with {nsteps} steps")

        if verbose:
            mdrun_cmd.extend(["-v"])

        # Run simulation
        logger.info("Starting MD simulation...")
        logger.info(f"Command: {' '.join(mdrun_cmd)}")

        success, stdout, stderr = run_command(mdrun_cmd, working_dir=output_dir, timeout=timeout)

        if not success:
            logger.error("MD simulation failed!")
            if stdout:
                logger.error(f"STDOUT:\n{stdout}")
            if stderr:
                logger.error(f"STDERR:\n{stderr}")
            return {
                "result": False,
                "error": "Simulation failed",
                "stdout": stdout,
                "stderr": stderr,
                "metadata": {
                    "input_file": str(input_path),
                    "output_dir": str(output_dir),
                    "config": config
                }
            }

        # Check and collect output files
        output_files = []
        expected_files = list(config['output_files'].values())

        for filename in expected_files:
            filepath = output_dir / filename
            if filepath.exists():
                size = get_file_size(filepath)
                file_info = f"{filename} ({size:,} bytes)"
                output_files.append(file_info)
                logger.info(f"Generated: {file_info}")
            else:
                logger.warning(f"Expected output file not found: {filename}")

        logger.info("MD simulation completed successfully!")

        if verbose:
            if stdout:
                logger.debug(f"STDOUT:\n{stdout}")
            if stderr:
                logger.debug(f"STDERR:\n{stderr}")

        return {
            "result": True,
            "output_files": output_files,
            "output_dir": str(output_dir),
            "metadata": {
                "input_file": str(input_path),
                "nsteps": nsteps,
                "config": config,
                "total_output_files": len(output_files)
            }
        }

    except Exception as e:
        logger.error(f"Error during MD simulation: {e}")
        return {
            "result": False,
            "error": str(e),
            "metadata": {
                "input_file": str(input_path),
                "config": config
            }
        }

# ==============================================================================
# CLI Interface
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input TPR file path'
    )

    parser.add_argument(
        '--output', '-o',
        help='Output directory path (default: md_output)'
    )

    parser.add_argument(
        '--nsteps', '-n',
        type=int,
        help='Number of simulation steps (default: 100 for demo)'
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
        '--timeout',
        type=int,
        help='Command timeout in seconds (default: 300)'
    )

    args = parser.parse_args()

    # Load config if provided
    config = None
    if args.config:
        import json
        with open(args.config) as f:
            config = json.load(f)

    # Prepare kwargs
    kwargs = {}
    if args.nsteps:
        kwargs['nsteps'] = args.nsteps
    if args.verbose:
        kwargs['verbose'] = args.verbose
    if args.timeout:
        kwargs['timeout'] = args.timeout

    try:
        # Run simulation
        result = run_md_simulation(
            input_file=args.input,
            output_file=args.output,
            config=config,
            **kwargs
        )

        if result['result']:
            print(f"✅ Success: MD simulation completed")
            print(f"   Output directory: {result['output_dir']}")
            print(f"   Generated files: {len(result['output_files'])}")
            for file_info in result['output_files']:
                print(f"   - {file_info}")
            return 0
        else:
            print(f"❌ Failed: {result.get('error', 'Unknown error')}")
            return 1

    except Exception as e:
        print(f"❌ Error: {e}")
        return 1

if __name__ == '__main__':
    sys.exit(main())