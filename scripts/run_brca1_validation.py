"""
Batch BRCA1 variant validation using Evo2 log-likelihood scoring.

Usage:
    python scripts/run_brca1_validation.py --clinvar path/to/variant_summary.txt

Downloads variant_summary.txt from:
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

Output:
    brca1_validation_results.csv  (40 rows: 20 Pathogenic + 20 Benign)
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Ensure repo root is on the path when run as a script
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scorer.client import score_variant
from scorer.clinvar import load_brca1_variants
from scorer.sequence import RefMismatchError, fetch_sequence_context

N_PER_CLASS = 20
OUT_FILE = Path(__file__).resolve().parents[1] / "brca1_validation_results.csv"


def score_variants_for_class(variants: pd.DataFrame, label: str) -> list[dict]:
    records = []
    skipped = 0

    for _, row in variants.iterrows():
        if len(records) >= N_PER_CLASS:
            break

        chrom = str(row["chromosome"])
        position = int(row["position"])
        ref_allele = row["ref_allele"]
        alt_allele = row["alt_allele"]
        name = row["name"]

        try:
            ref_seq, alt_seq = fetch_sequence_context(
                chrom, position, ref_allele, alt_allele
            )
        except RefMismatchError as e:
            print(f"  SKIP (ref mismatch): {name} — {e}")
            skipped += 1
            continue
        except Exception as e:
            print(f"  SKIP (sequence fetch error): {name} — {e}")
            skipped += 1
            continue

        try:
            result = score_variant(ref_seq, alt_seq)
        except Exception as e:
            print(f"  SKIP (scoring error): {name} — {e}")
            skipped += 1
            continue

        records.append(
            {
                "name": name,
                "clinical_significance": label,
                "chromosome": chrom,
                "position": position,
                "ref_allele": ref_allele,
                "alt_allele": alt_allele,
                "ref_ll": result["ref_ll"],
                "alt_ll": result["alt_ll"],
                "delta": result["delta"],
                "interpretation": result["interpretation"],
            }
        )
        print(
            f"  [{len(records)}/{N_PER_CLASS}] {name}  delta={result['delta']:.4f}"
            f"  {result['interpretation']}"
        )

    print(f"  Done: {len(records)} scored, {skipped} skipped\n")
    return records


def main():
    parser = argparse.ArgumentParser(description="Score BRCA1 ClinVar variants with Evo2")
    parser.add_argument(
        "--clinvar",
        required=True,
        help="Path to ClinVar variant_summary.txt (or .txt.gz)",
    )
    args = parser.parse_args()

    print("Loading ClinVar data…")
    pathogenic = load_brca1_variants(args.clinvar, clinical_significance="Pathogenic")
    benign = load_brca1_variants(args.clinvar, clinical_significance="Benign")
    print(f"  Pathogenic BRCA1 SNVs: {len(pathogenic)}")
    print(f"  Benign BRCA1 SNVs:     {len(benign)}\n")

    if len(pathogenic) < N_PER_CLASS:
        print(
            f"ERROR: only {len(pathogenic)} Pathogenic variants available, need {N_PER_CLASS}",
            file=sys.stderr,
        )
        sys.exit(1)
    if len(benign) < N_PER_CLASS:
        print(
            f"ERROR: only {len(benign)} Benign variants available, need {N_PER_CLASS}",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"=== Scoring {N_PER_CLASS} Pathogenic variants ===")
    path_records = score_variants_for_class(pathogenic, "Pathogenic")

    print(f"=== Scoring {N_PER_CLASS} Benign variants ===")
    benign_records = score_variants_for_class(benign, "Benign")

    if len(path_records) < N_PER_CLASS:
        print(
            f"ERROR: only scored {len(path_records)} Pathogenic variants "
            f"(need {N_PER_CLASS}). Increase ClinVar input or investigate skips.",
            file=sys.stderr,
        )
        sys.exit(1)
    if len(benign_records) < N_PER_CLASS:
        print(
            f"ERROR: only scored {len(benign_records)} Benign variants "
            f"(need {N_PER_CLASS}). Increase ClinVar input or investigate skips.",
            file=sys.stderr,
        )
        sys.exit(1)

    all_records = path_records[:N_PER_CLASS] + benign_records[:N_PER_CLASS]
    results_df = pd.DataFrame(all_records)
    results_df.to_csv(OUT_FILE, index=False)
    print(f"Results written to: {OUT_FILE}")
    print(f"Total rows: {len(results_df)}\n")

    mean_delta = results_df.groupby("clinical_significance")["delta"].mean()
    print("=== Mean delta by clinical significance ===")
    for sig, val in mean_delta.items():
        print(f"  {sig}: {val:.4f}")

    path_mean = mean_delta.get("Pathogenic")
    benign_mean = mean_delta.get("Benign")
    if path_mean is not None and benign_mean is not None:
        if path_mean < benign_mean:
            print("\nPASSED: Mean delta for Pathogenic < Mean delta for Benign")
        else:
            print(
                f"\nWARNING: Expected Pathogenic mean ({path_mean:.4f}) < "
                f"Benign mean ({benign_mean:.4f}), but this was not observed."
            )


if __name__ == "__main__":
    main()
