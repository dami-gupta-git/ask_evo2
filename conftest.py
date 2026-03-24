import sys
from pathlib import Path

# Make the ask_evo2 root importable so that `scorer`, `scripts`, and `src`
# are all resolvable as top-level packages.
sys.path.insert(0, str(Path(__file__).resolve().parent))
