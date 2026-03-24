import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
from sklearn.metrics import roc_curve, auc

results_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parents[1] / "data" / "brca1_scored.json"

with open(results_file) as f:
    records = json.load(f)

records = [r for r in records if r.get("delta") is not None]

pathogenic = [r for r in records if r["clinical_significance"] == "Pathogenic"]
benign     = [r for r in records if r["clinical_significance"] == "Benign"]

scores = [-r["delta"] for r in records]
labels = [1 if r["clinical_significance"] == "Pathogenic" else 0 for r in records]

fpr, tpr, _ = roc_curve(labels, scores)
roc_auc = auc(fpr, tpr)

# ── colour palette ──────────────────────────────────────────────────────────
BG      = "#F5F0E8"
SALMON  = "#E07B6A"
TEAL    = "#6BBFB0"
TEXT    = "#2A2A2A"
SUBTEXT = "#555555"

# ── bin deltas ───────────────────────────────────────────────────────────────
bins = [(-np.inf, -120), (-120, -80), (-80, -40), (-40, -20),
        (-20, -12), (-12, -8), (-8, -4), (-4, 0), (0, np.inf)]
bin_labels = ["≤-120", "-120 to -80", "-80 to -40", "-40 to -20",
              "-20 to -12", "-12 to -8", "-8 to -4", "-4 to 0", "0"]

def bin_deltas(recs):
    counts = []
    for lo, hi in bins:
        counts.append(sum(1 for r in recs if lo < r["delta"] <= hi))
    return counts

path_counts   = bin_deltas(pathogenic)
benign_counts = bin_deltas(benign)

x = np.arange(len(bins))
width = 0.4

# ── layout ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 7), facecolor=BG)
fig.patch.set_facecolor(BG)

gs = gridspec.GridSpec(
    2, 3,
    figure=fig,
    height_ratios=[0.22, 0.78],
    hspace=0.08,
    wspace=0.35,
    left=0.05, right=0.97, top=0.96, bottom=0.12,
)

# ── stat cards ───────────────────────────────────────────────────────────────
card_data = [
    ("AUC", f"{roc_auc:.2f}"),
    ("pathogenic variants", str(len(pathogenic))),
    ("benign variants",     str(len(benign))),
]

for col, (label, value) in enumerate(card_data):
    ax = fig.add_subplot(gs[0, col])
    ax.set_facecolor(BG)
    for spine in ax.spines.values():
        spine.set_edgecolor("#DDDDCC")
        spine.set_linewidth(1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(0.06, 0.72, label, transform=ax.transAxes,
            fontsize=10, color=SUBTEXT, va="top")
    ax.text(0.06, 0.38, value, transform=ax.transAxes,
            fontsize=26, fontweight="bold", color=TEXT, va="top")

# ── ROC curve ────────────────────────────────────────────────────────────────
ax_roc = fig.add_subplot(gs[1, 0])
ax_roc.set_facecolor(BG)
ax_roc.plot(fpr, tpr, color=SALMON, lw=2)
ax_roc.plot([0, 1], [0, 1], color="#AAAAAA", lw=1, linestyle="--")
ax_roc.set_xlim(0, 1); ax_roc.set_ylim(0, 1.02)
ax_roc.set_xlabel("FPR", color=SUBTEXT, fontsize=10)
ax_roc.set_ylabel("TPR", color=SUBTEXT, fontsize=10)
ax_roc.set_title("ROC curve", loc="left", color=TEXT, fontsize=11, pad=8)
ax_roc.tick_params(colors=SUBTEXT, labelsize=9)
for spine in ax_roc.spines.values():
    spine.set_edgecolor("#CCCCBB")

# ── delta histogram ──────────────────────────────────────────────────────────
ax_hist = fig.add_subplot(gs[1, 1:])
ax_hist.set_facecolor(BG)
ax_hist.bar(x - width/2, path_counts,   width=width, color=SALMON, label="Pathogenic")
ax_hist.bar(x + width/2, benign_counts, width=width, color=TEAL,   label="Benign")
ax_hist.set_xticks(x)
ax_hist.set_xticklabels(bin_labels, rotation=35, ha="right", fontsize=8.5, color=SUBTEXT)
ax_hist.set_ylabel("count", color=SUBTEXT, fontsize=10)
ax_hist.set_title("delta distribution", loc="left", color=TEXT, fontsize=11, pad=8)
ax_hist.tick_params(colors=SUBTEXT, labelsize=9)
ax_hist.legend(fontsize=9, framealpha=0)
for spine in ax_hist.spines.values():
    spine.set_edgecolor("#CCCCBB")

fig.text(0.05, 0.01,
         "score = −delta  (more negative delta → higher score → more likely pathogenic)",
         fontsize=8.5, color=SUBTEXT)

out = Path(__file__).parent / "brca1_roc_dashboard.png"
plt.savefig(out, dpi=150, facecolor=BG)
plt.show()
print(f"Saved to {out}")
print(f"AUC: {roc_auc:.3f}")
