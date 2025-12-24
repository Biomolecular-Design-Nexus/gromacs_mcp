#!/usr/bin/env python3
"""
Script: gromacs_command.py
Description: Execute GROMACS command-line operations

Original Use Case: examples/use_case_4_commandline_operations.py
Dependencies Removed: gmxapi (uses direct subprocess calls)

Usage:
    python scripts/gromacs_command.py --command <gromacs_tool> --input <files> --output <files>

Example:
    python scripts/gromacs_command.py --command rms --input trajectory.xtc --structure reference.gro --output rmsd.xvg
"""

# ==============================================================================
# Minimal Imports (only essential packages)
# ==============================================================================
import argparse
import sys
import os
import logging
import json
from pathlib import Path
from typing import Union, Optional, Dict, Any, List

# Import shared utilities
try:
    from lib.gromacs_utils import (
        find_gromacs_command, run_command, ensure_output_directory,
        check_file_exists
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

# ==============================================================================
# Configuration
# ==============================================================================
DEFAULT_CONFIG = {
    "commands": {
        # Common GROMACS analysis tools with their typical parameters
        "rms": {
            "description": "RMSD analysis",
            "required_inputs": ["structure", "trajectory"],
            "default_output": "rmsd.xvg",
            "parameters": ["-tu", "ns"]
        },
        "rmsf": {
            "description": "RMSF analysis",
            "required_inputs": ["structure", "trajectory"],
            "default_output": "rmsf.xvg",
            "parameters": ["-res"]
        },
        "energy": {
            "description": "Energy analysis",
            "required_inputs": ["energy"],
            "default_output": "energy.xvg",
            "parameters": []
        },
        "sasa": {
            "description": "Surface area analysis",
            "required_inputs": ["structure", "trajectory"],
            "default_output": "area.xvg",
            "parameters": []
        },
        "hbond": {
            "description": "Hydrogen bond analysis",
            "required_inputs": ["structure", "trajectory"],
            "default_output": "hbnum.xvg",
            "parameters": []
        },
        "trjconv": {
            "description": "Trajectory conversion",
            "required_inputs": ["structure", "trajectory"],
            "default_output": "trajout.xtc",
            "parameters": ["-pbc", "mol", "-center"]
        },
        "dump": {
            "description": "Dump file contents",
            "required_inputs": ["structure"],
            "default_output": "dump.txt",
            "parameters": []
        },
        "check": {
            "description": "Check file integrity",
            "required_inputs": ["structure"],
            "default_output": None,
            "parameters": []
        }
    },
    "execution": {
        "timeout": 300,
        "interactive_input": "echo '0' |",  # Default selection for interactive prompts
        "working_directory": None
    }
}

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# ==============================================================================
# Core Functions
# ==============================================================================
def run_gromacs_command(
    command_name: str,
    input_files: Optional[Dict[str, Union[str, Path]]] = None,
    output_files: Optional[Dict[str, Union[str, Path]]] = None,
    parameters: Optional[List[str]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Execute a GROMACS command with specified inputs, outputs, and parameters.

    Args:
        command_name: Name of GROMACS tool (e.g., 'rms', 'energy', 'trjconv')
        input_files: Dict mapping input types to file paths
        output_files: Dict mapping output types to file paths
        parameters: List of additional command parameters
        config: Configuration dict (uses DEFAULT_CONFIG if not provided)
        **kwargs: Override specific config parameters
            - verbose: Enable verbose logging
            - timeout: Command timeout in seconds
            - interactive_input: Input for interactive prompts

    Returns:
        Dict containing:
            - result: Command success status
            - command: Full command that was executed
            - output: Command output
            - output_files: List of generated files

    Example:
        >>> result = run_gromacs_command(
        ...     "rms",
        ...     input_files={"s": "ref.gro", "f": "traj.xtc"},
        ...     output_files={"o": "rmsd.xvg"}
        ... )
    """
    # Setup
    config = {**DEFAULT_CONFIG, **(config or {}), **kwargs}

    # Extract parameters
    verbose = kwargs.get('verbose', False)
    timeout = kwargs.get('timeout', config['execution']['timeout'])
    interactive_input = kwargs.get('interactive_input', config['execution']['interactive_input'])
    working_dir = kwargs.get('working_directory', config['execution']['working_directory'])

    if verbose:
        logger.setLevel(logging.DEBUG)

    logger.info(f"Running GROMACS command: gmx {command_name}")

    # Get command configuration
    cmd_config = config['commands'].get(command_name, {})

    try:
        # Find GROMACS command
        gmx_cmd = find_gromacs_command()

        # Build command
        full_cmd = gmx_cmd + [command_name]

        # Add input files
        if input_files:
            for flag, filepath in input_files.items():
                if filepath and check_file_exists(filepath, f"Input file (-{flag})"):
                    full_cmd.extend([f"-{flag}", str(filepath)])
                    logger.debug(f"Input: -{flag} {filepath}")

        # Add output files
        if output_files:
            for flag, filepath in output_files.items():
                if filepath:
                    # Ensure output directory exists
                    output_path = Path(filepath)
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    full_cmd.extend([f"-{flag}", str(filepath)])
                    logger.debug(f"Output: -{flag} {filepath}")

        # Add default parameters from config
        default_params = cmd_config.get('parameters', [])
        if default_params:
            full_cmd.extend(default_params)

        # Add custom parameters
        if parameters:
            full_cmd.extend(parameters)

        # Add verbose flag if requested
        if verbose:
            if command_name not in ['dump', 'check']:  # These don't support -v
                full_cmd.extend(["-v"])

        # Prepare command with interactive input if needed
        if interactive_input and command_name in ['rms', 'rmsf', 'energy', 'sasa', 'hbond']:
            # For commands that typically require interactive input
            cmd_str = f"{interactive_input} {' '.join(full_cmd)}"
            execute_cmd = ["bash", "-c", cmd_str]
        else:
            execute_cmd = full_cmd

        logger.info(f"Executing: {' '.join(execute_cmd)}")

        # Run command
        success, stdout, stderr = run_command(
            execute_cmd,
            working_dir=working_dir,
            timeout=timeout
        )

        # Collect output files
        output_file_info = []
        if output_files:
            for flag, filepath in output_files.items():
                if filepath and Path(filepath).exists():
                    size = Path(filepath).stat().st_size
                    output_file_info.append(f"{Path(filepath).name} ({size:,} bytes)")
                    logger.info(f"Generated: {Path(filepath).name} ({size:,} bytes)")

        if success:
            logger.info(f"GROMACS command '{command_name}' completed successfully")

            if verbose and stdout:
                logger.debug(f"STDOUT:\n{stdout}")
            if verbose and stderr:
                logger.debug(f"STDERR:\n{stderr}")

            return {
                "result": True,
                "command": ' '.join(full_cmd),
                "output": stdout,
                "output_files": output_file_info,
                "metadata": {
                    "command_name": command_name,
                    "input_files": input_files or {},
                    "output_files": output_files or {},
                    "parameters": parameters or []
                }
            }
        else:
            logger.error(f"GROMACS command '{command_name}' failed")
            if stdout:
                logger.error(f"STDOUT:\n{stdout}")
            if stderr:
                logger.error(f"STDERR:\n{stderr}")

            return {
                "result": False,
                "command": ' '.join(full_cmd),
                "error": stderr or "Command failed",
                "output": stdout,
                "metadata": {
                    "command_name": command_name,
                    "input_files": input_files or {},
                    "timeout": timeout
                }
            }

    except Exception as e:
        logger.error(f"Error executing GROMACS command '{command_name}': {e}")
        return {
            "result": False,
            "error": str(e),
            "command": command_name,
            "metadata": {
                "command_name": command_name,
                "config": config
            }
        }

def list_available_commands() -> Dict[str, str]:
    """
    List available GROMACS commands with descriptions.

    Returns:
    --------
    dict
        Command names mapped to descriptions
    """
    commands = {}
    for cmd, info in DEFAULT_CONFIG['commands'].items():
        commands[cmd] = info.get('description', f'GROMACS {cmd} command')
    return commands

def get_command_help(command_name: str) -> str:
    """
    Get help information for a specific command.

    Parameters:
    -----------
    command_name : str
        Name of GROMACS command

    Returns:
    --------
    str
        Help text for the command
    """
    cmd_config = DEFAULT_CONFIG['commands'].get(command_name, {})
    if not cmd_config:
        return f"Unknown command: {command_name}"

    help_text = f"GROMACS {command_name} - {cmd_config.get('description', 'No description')}\n\n"
    help_text += f"Required inputs: {', '.join(cmd_config.get('required_inputs', []))}\n"
    help_text += f"Default output: {cmd_config.get('default_output', 'None')}\n"
    help_text += f"Default parameters: {' '.join(cmd_config.get('parameters', []))}\n"

    return help_text

# ==============================================================================
# CLI Interface
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--command', '-cmd',
        help='GROMACS command to execute (e.g., rms, energy, trjconv)'
    )

    parser.add_argument(
        '--input', '-i',
        action='append',
        nargs=2,
        metavar=('FLAG', 'FILE'),
        help='Input file: -i FLAG FILE (can be used multiple times)'
    )

    parser.add_argument(
        '--output', '-o',
        action='append',
        nargs=2,
        metavar=('FLAG', 'FILE'),
        help='Output file: -o FLAG FILE (can be used multiple times)'
    )

    parser.add_argument(
        '--param', '-p',
        action='append',
        help='Additional parameter (can be used multiple times)'
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

    parser.add_argument(
        '--list-commands',
        action='store_true',
        help='List available GROMACS commands'
    )

    parser.add_argument(
        '--help-command',
        help='Show help for specific command'
    )

    args = parser.parse_args()

    # Handle special options
    if args.list_commands:
        print("Available GROMACS commands:")
        for cmd, desc in list_available_commands().items():
            print(f"  {cmd:<12} - {desc}")
        return 0

    if args.help_command:
        print(get_command_help(args.help_command))
        return 0

    # Validate required command for normal operations
    if not args.command:
        parser.error("--command is required unless using --list-commands or --help-command")

    # Load config if provided
    config = None
    if args.config:
        with open(args.config) as f:
            config = json.load(f)

    # Parse input files
    input_files = {}
    if args.input:
        for flag, filepath in args.input:
            input_files[flag] = filepath

    # Parse output files
    output_files = {}
    if args.output:
        for flag, filepath in args.output:
            output_files[flag] = filepath

    # Prepare kwargs
    kwargs = {}
    if args.verbose:
        kwargs['verbose'] = args.verbose
    if args.timeout:
        kwargs['timeout'] = args.timeout

    try:
        # Run command
        result = run_gromacs_command(
            command_name=args.command,
            input_files=input_files or None,
            output_files=output_files or None,
            parameters=args.param,
            config=config,
            **kwargs
        )

        if result['result']:
            print(f"✅ Success: GROMACS {args.command} completed")
            if result['output_files']:
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