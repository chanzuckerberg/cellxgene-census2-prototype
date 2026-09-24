"""Make `tools/` importable, so the tests can build a catalog with the builder script."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
