#!/usr/bin/env python3
"""Shared utilities for GROMACS MCP server."""

from pathlib import Path
from typing import Dict, Any, Optional
import logging

def validate_input_file(file_path: str, expected_extension: str = None) -> Dict[str, Any]:
    """Validate input file exists and has correct extension.

    Args:
        file_path: Path to the file to validate
        expected_extension: Expected file extension (e.g., '.tpr', '.gro')

    Returns:
        Dict with validation result and error message if any
    """
    try:
        path = Path(file_path)

        if not path.exists():
            return {
                "valid": False,
                "error": f"File not found: {file_path}"
            }

        if not path.is_file():
            return {
                "valid": False,
                "error": f"Path is not a file: {file_path}"
            }

        if expected_extension and not path.suffix.lower() == expected_extension.lower():
            return {
                "valid": False,
                "error": f"Expected {expected_extension} file, got {path.suffix}"
            }

        return {
            "valid": True,
            "file_size": path.stat().st_size,
            "absolute_path": str(path.resolve())
        }

    except Exception as e:
        return {
            "valid": False,
            "error": f"Validation error: {str(e)}"
        }

def setup_output_directory(output_path: Optional[str], default_name: str = "output") -> Dict[str, Any]:
    """Setup output directory for results.

    Args:
        output_path: User-specified output path (file or directory)
        default_name: Default directory name if none specified

    Returns:
        Dict with output directory path and any setup errors
    """
    try:
        if output_path:
            path = Path(output_path)
            if path.suffix:  # It's a file path
                output_dir = path.parent
                output_file = path
            else:  # It's a directory path
                output_dir = path
                output_file = None
        else:
            # Create default output directory
            output_dir = Path.cwd() / "results" / default_name
            output_file = None

        # Create directory if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)

        return {
            "success": True,
            "output_dir": str(output_dir.resolve()),
            "output_file": str(output_file.resolve()) if output_file else None
        }

    except Exception as e:
        return {
            "success": False,
            "error": f"Failed to setup output directory: {str(e)}"
        }

def format_tool_result(
    success: bool,
    data: Dict[str, Any] = None,
    error_message: str = None,
    output_files: list = None
) -> Dict[str, Any]:
    """Format consistent result for MCP tools.

    Args:
        success: Whether the operation succeeded
        data: Result data
        error_message: Error message if failed
        output_files: List of output files created

    Returns:
        Formatted result dictionary
    """
    result = {
        "status": "success" if success else "error"
    }

    if success:
        if data:
            result.update(data)
        if output_files:
            result["output_files"] = output_files
    else:
        result["error"] = error_message or "Unknown error"

    return result

def get_example_data_path() -> Path:
    """Get path to example data directory."""
    return Path(__file__).parent.parent / "examples" / "data"

def configure_logging(level: str = "INFO") -> None:
    """Configure logging for the MCP server."""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('gromacs_mcp.log')
        ]
    )