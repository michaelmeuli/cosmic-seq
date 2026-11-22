"""MANIFEST.txt of an output directory via the shared `kansasii-docs/bin/write_manifest.sh`."""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from . import REPO
from .kansasii_paths import MANIFEST_SCRIPT


def write_manifest(out: Path, script: str, *inputs: str) -> None:
    """Write `out/MANIFEST.txt`; skipped with a warning where bash or the script is missing (e.g. on Windows)."""
    bash = shutil.which("bash")
    if bash is None or not MANIFEST_SCRIPT.is_file():
        print(f"warning: no MANIFEST.txt (bash or {MANIFEST_SCRIPT} missing)", file=sys.stderr)
        return
    subprocess.run([bash, str(MANIFEST_SCRIPT), str(out), str(REPO), script, *inputs], check=True)
