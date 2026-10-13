"""Repository path helpers.

Override with environment variable MYSCRIPTS_ROOT if the repo is not at ~/myScripts.
"""

from __future__ import annotations

import os
from pathlib import Path


def _find_root() -> Path:
    env = os.environ.get("MYSCRIPTS_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    # src/myscripts/paths.py → repo root is parents[2]
    here = Path(__file__).resolve()
    candidate = here.parents[2]
    if (candidate / "src" / "md").is_dir():
        return candidate
    home = Path.home() / "myScripts"
    if (home / "src" / "md").is_dir():
        return home
    return candidate


ROOT = _find_root()
DATA = ROOT / "data"
REF = DATA / "reference" / "KinaseStructures"
MD = ROOT / "src" / "md"
BIOINFO = ROOT / "src" / "bioinformatics"
ML = ROOT / "src" / "machine_learning"
UTILS = ROOT / "utils"
