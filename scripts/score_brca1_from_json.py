"""
Read brca1_variants.json, score each variant via the Modal endpoint,
and write results (all original fields + scoring output) to brca1_scored.json.

Usage:
    python scripts/score_brca1_from_json.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.client import score_variant

DATA_FILE = Path(__file__).resolve().parents[1] / "data" / "brca1_variants.json"
OUT_FILE  = Path(__file__).resolve().parents[1] / "data" / "brca1_scored.json"


def main():
    with open(DATA_FILE) as f:
        variants = json.load(f)

    print(f"Loaded {len(variants)} variants from {DATA_FILE.name}\n")
    print(f"{'#':<4} {'Sig':<12} {'Name':<50} {'delta':>10}  interpretation")
    print("-" * 100)

    results = []
    for i, v in enumerate(variants, 1):
        name = v["name"]
        sig  = v["clinical_significance"]
        try:
            scored = score_variant(v["ref_seq"], v["alt_seq"])
            record = {**v, **scored}
            print(f"{i:<4} {sig:<12} {name[:49]:<50} {scored['delta']:>10.4f}  {scored['interpretation']}")
        except Exception as e:
            record = {**v, "ref_ll": None, "alt_ll": None, "delta": None, "interpretation": f"ERROR: {e}"}
            print(f"{i:<4} {sig:<12} {name[:49]:<50} {'SKIP':>10}  {e}")
        results.append(record)

    with open(OUT_FILE, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults written to {OUT_FILE.name} ({len(results)} records)")


if __name__ == "__main__":
    main()
