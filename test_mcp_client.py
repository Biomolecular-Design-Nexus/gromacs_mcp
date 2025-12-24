#!/usr/bin/env python3
"""Simple test client for GROMACS MCP server."""

import json
import subprocess
import sys
from pathlib import Path

def test_server_basic():
    """Test basic server functionality."""
    try:
        # Test server startup and tool discovery
        cmd = [
            "mamba", "run", "-p", str(Path.cwd() / "env"),
            "python", "-c",
            """
import sys
sys.path.insert(0, 'src')
from server import mcp

# Test that tools are registered
tools = [tool.name for tool in mcp._tools.values()]
print(f'Registered tools: {len(tools)}')
print('Tools:', tools)

# Test a simple utility function
from utils import validate_input_file
result = validate_input_file('examples/data/sample.tpr', '.tpr')
print(f'Validation result: {result}')
"""
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print("✅ Server basic test PASSED")
            print(result.stdout)
        else:
            print("❌ Server basic test FAILED")
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)

    except Exception as e:
        print(f"❌ Test failed with exception: {e}")

def test_fastmcp_dev():
    """Test using fastmcp dev mode."""
    try:
        print("Testing fastmcp dev mode (will timeout after 5 seconds)...")
        cmd = [
            "mamba", "run", "-p", str(Path.cwd() / "env"),
            "timeout", "5s",
            "fastmcp", "dev", "src/server.py"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Expect timeout (exit code 124) since we're testing startup
        if result.returncode == 124:
            print("✅ Server startup test PASSED (timed out as expected)")
        elif result.returncode == 0:
            print("✅ Server startup test PASSED")
        else:
            print("❌ Server startup test FAILED")
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)

    except Exception as e:
        print(f"❌ Test failed with exception: {e}")

if __name__ == "__main__":
    print("=== GROMACS MCP Server Tests ===")

    print("\n1. Testing basic server functionality...")
    test_server_basic()

    print("\n2. Testing server startup...")
    test_fastmcp_dev()

    print("\n=== Test Summary ===")
    print("If both tests passed, the server is working correctly!")
    print("Next steps:")
    print("- Run: fastmcp dev src/server.py")
    print("- Test with Claude Desktop or MCP inspector")
    print("- Use example: analyze_tpr with input_file 'examples/data/sample.tpr'")