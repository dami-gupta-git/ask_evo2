# Usage: python map/draw_roc.py <scored_variants.json> [--output roc.png] [--title "My Title"]
# Example: python map/draw_roc.py data/brca1_scored.json --output map/brca1_roc.png --title "BRCA1 ROC Curve"

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

# draws roc
def main():
    parser = argparse.ArgumentParser(description="Draw ROC curve from scored variants JSON.")
    parser.add_argument("results", type=Path, help="Path to scored variants JSON file")
    parser.add_argument("--output", type=Path, help="Output PNG path (default: <results_stem>_roc.png)")
    parser.add_argument("--title", help="Plot title (default: derived from results filename)")
    args = parser.parse_args()

    with open(args.results) as f:
        records = json.load(f)

    records = [r for r in records if r.get("delta") is not None]

    scores = [-r["delta"] for r in records]
    labels = [1 if r["clinical_significance"] == "Pathogenic" else 0 for r in records]

    fpr, tpr, thresholds = roc_curve(labels, scores)
    roc_auc = auc(fpr, tpr)

    title = args.title or f"{args.results.stem} — ROC Curve"
    output = args.output or args.results.with_name(args.results.stem + "_roc.png")

    plt.figure(figsize=(6, 6))
    plt.plot(fpr, tpr, color='#1D9E75', lw=2, label=f'Evo2 (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--', label='random')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    plt.show()

    print(f"AUC: {roc_auc:.3f}")
    print(f"Saved: {output}")
    print(f"\nThreshold analysis:")
    for f, t, thresh in zip(fpr, tpr, thresholds):
        print(f"  threshold={thresh:.1f}  TPR={t:.2f}  FPR={f:.2f}")


if __name__ == "__main__":
    main()