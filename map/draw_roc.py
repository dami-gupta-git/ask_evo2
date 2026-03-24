import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

results_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parents[1] / "data" / "brca1_scored.json"

with open(results_file) as f:
    records = json.load(f)

records = [r for r in records if r.get("delta") is not None]

scores = [-r["delta"] for r in records]
labels = [1 if r["clinical_significance"] == "Pathogenic" else 0 for r in records]

fpr, tpr, thresholds = roc_curve(labels, scores)
roc_auc = auc(fpr, tpr)

plt.figure(figsize=(6, 6))
plt.plot(fpr, tpr, color='#1D9E75', lw=2, label=f'Evo2 (AUC = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--', label='random')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('BRCA1 Variant Scoring — ROC Curve')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig('brca1_roc.png', dpi=150)
plt.show()

print(f"AUC: {roc_auc:.3f}")
print(f"\nThreshold analysis:")
for f, t, thresh in zip(fpr, tpr, thresholds):
    print(f"  threshold={thresh:.1f}  TPR={t:.2f}  FPR={f:.2f}")