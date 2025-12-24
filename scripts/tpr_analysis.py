#!/usr/bin/env python3
"""
Script: tpr_analysis.py
Description: Read and analyze GROMACS TPR files

Original Use Case: examples/use_case_2_read_tpr_file.py
Dependencies Removed: None (uses standard library + GROMACS command)

Usage:
    python scripts/tpr_analysis.py --input <tpr_file> --output <output_file>

Example:
    python scripts/tpr_analysis.py --input examples/data/sample.tpr --output results/analysis.txt
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
from typing import Union, Optional, Dict, Any

# Import shared utilities
try:
    from lib.gromacs_utils import (
        find_gromacs_command, run_command, ensure_output_directory,
        check_file_exists, get_file_size, parse_tpr_dump_output
    )
except ImportError:
    # Fallback if lib not available
    import shutil
    import subprocess
    import re

    def find_gromacs_command():
        gmx_cmd = shutil.which("gmx")
        if gmx_cmd:
            return ["gmx"]
        script_dir = Path(__file__).parent.parent
        env_path = script_dir / "env"
        if env_path.exists():
            return ["mamba", "run", "-p", str(env_path), "gmx"]
        return ["gmx"]

    def run_command(cmd, working_dir=None, timeout=60):
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

    def get_file_size(file_path):
        try:
            return Path(file_path).stat().st_size
        except:
            return 0

    def parse_tpr_dump_output(dump_output):
        info = {}
        patterns = {
            'nsteps': r'nsteps\s*=\s*(\d+)',
            'dt': r'delta-t\s*=\s*([\d\.e\-\+]+)',
            'natoms': r'#atoms\s*=\s*(\d+)',
            'title': r'title:\s*"([^"]*)"',
            'integrator': r'integrator\s*=\s*(\w+)',
            'nstlist': r'nstlist\s*=\s*(\d+)',
            'temperature': r'ref-t.*?(\d+\.?\d*)',
            'pressure': r'ref-p.*?(\d+\.?\d*)',
        }

        for key, pattern in patterns.items():
            match = re.search(pattern, dump_output, re.IGNORECASE)
            if match:
                try:
                    value = match.group(1)
                    if key in ['nsteps', 'natoms', 'nstlist']:
                        info[key] = int(value)
                    elif key in ['dt', 'temperature', 'pressure']:
                        info[key] = float(value)
                    else:
                        info[key] = value
                except:
                    info[key] = match.group(1)

        atom_lines = [line for line in dump_output.split('\n') if 'atom[' in line.lower()]
        info['total_atoms'] = len(atom_lines)
        return info

# ==============================================================================
# Configuration
# ==============================================================================
DEFAULT_CONFIG = {
    "analysis": {
        "timeout": 60,          # 1 minute timeout for dump command
        "include_raw_dump": True,  # Include raw gmx dump output
        "output_format": "text"    # or "json"
    },
    "output_sections": [
        "file_info",
        "simulation_parameters",
        "system_information",
        "raw_dump"
    ]
}

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# ==============================================================================
# Core Function
# ==============================================================================
def run_tpr_analysis(
    input_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Read and analyze a GROMACS TPR file using gmx dump.

    Args:
        input_file: Path to TPR input file
        output_file: Path to save analysis results (if None, return dict only)
        config: Configuration dict (uses DEFAULT_CONFIG if not provided)
        **kwargs: Override specific config parameters
            - verbose: Enable verbose logging
            - format: Output format ("text", "json")
            - include_raw: Include raw dump output

    Returns:
        Dict containing:
            - result: Analysis success status
            - tpr_info: Parsed TPR file information
            - analysis_text: Formatted analysis text
            - output_file: Path to output file (if saved)

    Example:
        >>> result = run_tpr_analysis("input.tpr", "analysis.txt")
        >>> print(result['tpr_info']['natoms'])
    """
    # Setup
    input_path = Path(input_file)
    config = {**DEFAULT_CONFIG, **(config or {}), **kwargs}

    # Extract parameters
    verbose = kwargs.get('verbose', False)
    output_format = kwargs.get('format', config['analysis']['output_format'])
    include_raw = kwargs.get('include_raw', config['analysis']['include_raw_dump'])
    timeout = kwargs.get('timeout', config['analysis']['timeout'])

    if verbose:
        logger.setLevel(logging.DEBUG)

    logger.info(f"Analyzing TPR file: {input_path}")

    # Validate input
    if not check_file_exists(input_path, "TPR file"):
        raise FileNotFoundError(f"TPR file not found: {input_path}")

    try:
        # Find GROMACS command
        gmx_cmd = find_gromacs_command()
        logger.debug(f"Using GROMACS command: {' '.join(gmx_cmd)}")

        # Run gmx dump command
        dump_cmd = gmx_cmd + ["dump", "-s", str(input_path)]

        logger.info("Extracting TPR file information using gmx dump...")
        success, stdout, stderr = run_command(dump_cmd, timeout=timeout)

        if not success:
            logger.error("Failed to run gmx dump")
            if stderr:
                logger.error(f"Error: {stderr}")
            return {
                "result": False,
                "error": "Failed to run gmx dump",
                "stderr": stderr
            }

        # Parse the output
        logger.info("Parsing TPR file information...")
        tpr_info = parse_tpr_dump_output(stdout)

        # Add file information
        file_stat = input_path.stat()
        tpr_info['file_size'] = file_stat.st_size
        tpr_info['file_path'] = str(input_path.absolute())

        # Generate analysis text
        analysis_text = format_analysis_text(tpr_info, include_raw_dump=include_raw, raw_dump=stdout)

        # Save output if requested
        output_path = None
        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            if output_format.lower() == "json":
                # Save as JSON
                output_data = {
                    "tpr_info": tpr_info,
                    "analysis_text": analysis_text
                }
                if include_raw:
                    output_data["raw_dump"] = stdout

                with open(output_path, 'w') as f:
                    json.dump(output_data, f, indent=2)
            else:
                # Save as text
                with open(output_path, 'w') as f:
                    f.write(analysis_text)

            logger.info(f"Analysis saved to: {output_path}")

        # Print analysis to console
        print(analysis_text)

        logger.info("TPR file analysis completed successfully")

        return {
            "result": True,
            "tpr_info": tpr_info,
            "analysis_text": analysis_text,
            "output_file": str(output_path) if output_path else None,
            "metadata": {
                "input_file": str(input_path),
                "output_format": output_format,
                "file_size": tpr_info['file_size']
            }
        }

    except Exception as e:
        logger.error(f"Error during TPR analysis: {e}")
        return {
            "result": False,
            "error": str(e),
            "metadata": {
                "input_file": str(input_path),
                "config": config
            }
        }

def format_analysis_text(tpr_info: Dict[str, Any], include_raw_dump: bool = False, raw_dump: str = "") -> str:
    """
    Format TPR analysis information as readable text.

    Parameters:
    -----------
    tpr_info : dict
        Parsed TPR information
    include_raw_dump : bool
        Whether to include raw dump output
    raw_dump : str
        Raw dump output

    Returns:
    --------
    str
        Formatted analysis text
    """
    lines = []
    lines.append("=== TPR File Analysis Results ===")
    lines.append(f"File: {tpr_info.get('file_path', 'Unknown')}")
    lines.append(f"Size: {tpr_info.get('file_size', 0):,} bytes")
    lines.append("")

    # Simulation Parameters
    lines.append("Simulation Parameters:")
    if 'title' in tpr_info:
        lines.append(f"  Title: {tpr_info['title']}")
    if 'integrator' in tpr_info:
        lines.append(f"  Integrator: {tpr_info['integrator']}")
    if 'nsteps' in tpr_info:
        lines.append(f"  Number of steps: {tpr_info['nsteps']:,}")
    if 'dt' in tpr_info:
        lines.append(f"  Time step: {tpr_info['dt']} ps")
        if 'nsteps' in tpr_info:
            total_time = tpr_info['nsteps'] * tpr_info['dt']
            lines.append(f"  Total simulation time: {total_time:.1f} ps")

    lines.append("")

    # System Information
    lines.append("System Information:")
    if 'natoms' in tpr_info:
        lines.append(f"  Number of atoms: {tpr_info['natoms']:,}")
    if 'total_atoms' in tpr_info and tpr_info['total_atoms'] != tpr_info.get('natoms', 0):
        lines.append(f"  Atoms found in dump: {tpr_info['total_atoms']:,}")

    if 'temperature' in tpr_info:
        lines.append(f"  Reference temperature: {tpr_info['temperature']} K")
    if 'pressure' in tpr_info:
        lines.append(f"  Reference pressure: {tpr_info['pressure']} bar")

    # Add raw dump if requested
    if include_raw_dump and raw_dump:
        lines.append("")
        lines.append("=== Raw gmx dump output ===")
        lines.append(raw_dump)

    return '\n'.join(lines)

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
        help='Output file path (optional)'
    )

    parser.add_argument(
        '--format', '-f',
        choices=['text', 'json'],
        default='text',
        help='Output format (default: text)'
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
        '--include-raw',
        action='store_true',
        help='Include raw gmx dump output'
    )

    parser.add_argument(
        '--info',
        action='store_true',
        help='Show TPR file format information'
    )

    args = parser.parse_args()

    if args.info:
        print("TPR File Format Information:")
        print("- TPR (Topology and Parameters Run) files contain all simulation settings")
        print("- Generated by gmx grompp from .mdp, .gro/.pdb, and .top files")
        print("- Binary format containing topology, parameters, and coordinates")
        print("- Required input for gmx mdrun")
        print("- Can be analyzed using gmx dump, gmx check, etc.")
        return 0

    # Load config if provided
    config = None
    if args.config:
        with open(args.config) as f:
            config = json.load(f)

    # Prepare kwargs
    kwargs = {}
    if args.verbose:
        kwargs['verbose'] = args.verbose
    if args.format:
        kwargs['format'] = args.format
    if args.include_raw:
        kwargs['include_raw'] = args.include_raw

    try:
        # Run analysis
        result = run_tpr_analysis(
            input_file=args.input,
            output_file=args.output,
            config=config,
            **kwargs
        )

        if result['result']:
            if result['output_file']:
                print(f"✅ Success: Analysis saved to {result['output_file']}")
            else:
                print("✅ Success: Analysis completed")
            return 0
        else:
            print(f"❌ Failed: {result.get('error', 'Unknown error')}")
            return 1

    except Exception as e:
        print(f"❌ Error: {e}")
        return 1

if __name__ == '__main__':
    sys.exit(main())