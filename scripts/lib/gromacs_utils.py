#!/usr/bin/env python3
"""
GROMACS Utilities Library

Shared functions for GROMACS operations extracted from use case scripts.
These utilities minimize dependencies and provide self-contained functionality.
"""

import os
import shutil
import subprocess
import logging
import re
from pathlib import Path
from typing import Tuple, Dict, List, Optional, Union

logger = logging.getLogger(__name__)

def find_gromacs_command() -> List[str]:
    """
    Find the GROMACS command, trying different approaches.

    Returns:
    --------
    list
        Command to use for GROMACS (either ['gmx'] or ['mamba', 'run', '-p', env_path, 'gmx'])
    """
    # Check if gmx is available in PATH
    gmx_cmd = shutil.which("gmx")
    if gmx_cmd:
        return ["gmx"]

    # Try mamba run approach with environment path
    try:
        # Calculate environment path relative to script location
        script_dir = Path(__file__).parent.parent.parent  # Go up from scripts/lib/
        env_path = script_dir / "env"

        if env_path.exists():
            return ["mamba", "run", "-p", str(env_path), "gmx"]
    except Exception:
        pass

    # Fallback to plain gmx (will fail if not in PATH)
    return ["gmx"]

def run_command(cmd: List[str], working_dir: Optional[Union[str, Path]] = None, timeout: int = 300) -> Tuple[bool, str, str]:
    """
    Run a system command and capture output.

    Parameters:
    -----------
    cmd : list
        Command and arguments to run
    working_dir : str or Path, optional
        Working directory for command execution
    timeout : int
        Timeout in seconds (default: 5 minutes)

    Returns:
    --------
    tuple
        (success: bool, stdout: str, stderr: str)
    """
    try:
        logger.debug(f"Running command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=working_dir,
            timeout=timeout
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after {timeout} seconds")
        return False, "", f"Command timed out after {timeout} seconds"
    except Exception as e:
        logger.error(f"Command failed with exception: {e}")
        return False, "", str(e)

def ensure_output_directory(output_dir: Union[str, Path]) -> Path:
    """
    Ensure output directory exists and return as Path object.

    Parameters:
    -----------
    output_dir : str or Path
        Output directory path

    Returns:
    --------
    Path
        Output directory as Path object
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    return output_path

def check_file_exists(file_path: Union[str, Path], description: str = "File") -> bool:
    """
    Check if file exists and log appropriate message.

    Parameters:
    -----------
    file_path : str or Path
        Path to check
    description : str
        Description of file type for logging

    Returns:
    --------
    bool
        True if file exists
    """
    path = Path(file_path)
    if not path.exists():
        logger.error(f"{description} not found: {file_path}")
        return False
    return True

def get_file_size(file_path: Union[str, Path]) -> int:
    """
    Get file size in bytes.

    Parameters:
    -----------
    file_path : str or Path
        Path to file

    Returns:
    --------
    int
        File size in bytes (0 if file doesn't exist)
    """
    try:
        return Path(file_path).stat().st_size
    except:
        return 0

def parse_tpr_dump_output(dump_output: str) -> Dict[str, Union[str, int, float]]:
    """
    Parse the output from gmx dump command to extract key information.
    Extracted and simplified from use_case_2_read_tpr_file.py

    Parameters:
    -----------
    dump_output : str
        Output from gmx dump command

    Returns:
    --------
    dict
        Parsed TPR information
    """
    info = {}

    # Extract basic simulation parameters
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
                # Try to convert to appropriate type
                value = match.group(1)
                if key in ['nsteps', 'natoms', 'nstlist']:
                    info[key] = int(value)
                elif key in ['dt', 'temperature', 'pressure']:
                    info[key] = float(value)
                else:
                    info[key] = value
            except:
                info[key] = match.group(1)

    # Count different atom types
    atom_lines = [line for line in dump_output.split('\n') if 'atom[' in line.lower()]
    info['total_atoms'] = len(atom_lines)

    return info

def format_file_info(file_path: Union[str, Path]) -> str:
    """
    Format file information for display.

    Parameters:
    -----------
    file_path : str or Path
        Path to file

    Returns:
    --------
    str
        Formatted file information
    """
    path = Path(file_path)
    if path.exists():
        size = get_file_size(path)
        return f"{path.name} ({size:,} bytes)"
    else:
        return f"{path.name} (not found)"

def copy_input_file(input_file: Union[str, Path], output_dir: Union[str, Path],
                   new_name: Optional[str] = None) -> Path:
    """
    Copy input file to output directory with optional renaming.

    Parameters:
    -----------
    input_file : str or Path
        Source file
    output_dir : str or Path
        Destination directory
    new_name : str, optional
        New filename (uses original if not provided)

    Returns:
    --------
    Path
        Path to copied file
    """
    src_path = Path(input_file)
    dst_dir = Path(output_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    if new_name:
        dst_path = dst_dir / new_name
    else:
        dst_path = dst_dir / src_path.name

    shutil.copy2(src_path, dst_path)
    logger.info(f"Copied {src_path} to {dst_path}")
    return dst_path