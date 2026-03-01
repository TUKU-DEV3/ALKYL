import subprocess
import json
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"
SCRIPT_TIMEOUT = 30  # seconds

def run_script(script_name: str, args: list[str]) -> dict:
    """Run a chem_*.py script and parse its JSON stdout."""
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / script_name)] + args,
        capture_output=True,
        text=True,
        timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0, f"Script failed:\n{result.stderr}"
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise AssertionError(
            f"Script output is not valid JSON.\nstdout: {result.stdout!r}"
        ) from exc

# Common test molecules
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE_SMILES = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
