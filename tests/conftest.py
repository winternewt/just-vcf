import sys
from pathlib import Path

# scripts/ holds dev tools that import each other by module name (pgp_participants imports pgp_leech).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
