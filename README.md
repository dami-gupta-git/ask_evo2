---
title: AskEvo2
emoji: 🧬
colorFrom: green
colorTo: blue
sdk: gradio
sdk_version: "5.29.0"
app_file: app.py
pinned: false
license: mit
tags:
  - biology
  - genomics
  - variant-effect-prediction
  - evo2
  - bioinformatics
---

# AskEvo2

Zero-shot variant effect scoring using Evo2 7B. Given a reference and alternate DNA sequence, AskEvo2 computes a delta log-likelihood score to estimate whether a variant is likely deleterious or neutral.

> **For research use only. Not a clinical tool.**

---

## How it works

AskEvo2 runs both sequences through Evo2 7B using a teacher-forcing log-likelihood calculation:

1. Each sequence is tokenized and passed through Evo2 7B in a single forward pass.
2. Per-token log-probabilities are extracted via log-softmax on the logits.
3. The sequence log-likelihood is the sum of per-token log-probs (i.e., the model's probability of generating each nucleotide given all preceding context).
4. Delta is computed as: **delta = alt_LL − ref_LL**

### Interpreting delta

| Delta | Interpretation |
|-------|----------------|
| < −2.0 | Strongly deleterious |
| −2.0 to −1.0 | Likely deleterious |
| −1.0 to −0.5 | Possibly deleterious |
| −0.5 to 0.5 | Uncertain |
| > 0.5 | Likely neutral |

A negative delta means the model assigns lower likelihood to the alternate sequence — the mutation makes the sequence less consistent with patterns learned from natural genomes, which correlates with functional disruption.

---

## Prerequisites

- Python 3.10+
- A [Modal](https://modal.com) account (free tier is sufficient for testing)
- Modal CLI: `pip install modal && modal setup`

---

## Setup and deployment

**1. Install local dependencies**
```bash
pip install -r requirements.txt
```

**2. Authenticate Modal CLI**
```bash
modal setup
```

**3. Deploy the GPU backend**
```bash
modal deploy modal_app.py
```

Modal will print a URL like:
```
https://your-username--askevo2-scorer-score-variant.modal.run
```

**4. Set the endpoint URL in `app.py`**

Open [app.py](app.py) and replace the placeholder:
```python
MODAL_ENDPOINT_URL = "https://your-username--askevo2-scorer-score-variant.modal.run"
```

**5. Run the Gradio frontend**
```bash
python app.py
```

Open http://localhost:7860 in your browser.

---

## API usage

The Modal backend accepts POST requests directly:

```bash
curl -X POST https://<your-endpoint-url> \
  -H "Content-Type: application/json" \
  -d '{"ref_sequence": "ATGGATTTATCTGCTC", "alt_sequence": "ATGGATTTATCTGCTG"}'
```

Response:
```json
{
  "ref_ll": -12.3456,
  "alt_ll": -13.8901,
  "delta": -1.5445,
  "interpretation": "Likely deleterious"
}
```

---

## Sequence requirements

- Uppercase A, C, G, T only (lowercase is auto-normalized)
- Minimum 10 bp
- Maximum 10,000 bp
- No FASTA headers, spaces, or ambiguity codes (N, R, Y, etc.)

---

## Hosting on HuggingFace Spaces (free CPU)

1. Create a new Space on HuggingFace with the **Gradio** SDK and **CPU Basic** (free) hardware.
2. Upload `app.py` and `requirements.txt`.
3. Set `MODAL_ENDPOINT_URL` as a Space secret or hardcode it in `app.py`.
4. The Space will install dependencies and launch the Gradio app automatically.

---

## Cost notes

- Modal charges for A10G GPU time only while a request is being processed.
- `scaledown_window=300` keeps the container warm for 5 minutes after the last request, avoiding cold starts during active sessions.
- First request after idle: ~30–60 second cold start (model loading).
- Subsequent requests within the warm window: ~5–15 seconds each.

---

## Repository layout

```
ask_evo2/
├── src/
│   ├── modal_app.py          # Modal GPU backend: Evo2 model, log-likelihood, POST endpoint
│   └── app.py                # Gradio frontend: input validation, calls Modal endpoint
│
├── scorer/
│   ├── client.py             # Re-exports score_variant() from scripts/client.py
│   ├── clinvar.py            # Load ClinVar variant_summary.txt, filter BRCA1 SNVs on GRCh38
│   └── sequence.py           # Fetch hg38 flanking sequence from UCSC; validate ref base
│
├── scripts/
│   ├── client.py             # score_variant(ref_seq, alt_seq) — POSTs to Modal endpoint
│   ├── score_brca1_from_json.py   # Reads data/brca1_variants.json, scores all 40 variants,
│   │                              # writes data/brca1_scored.json
│   └── run_brca1_validation.py    # Alternative: loads from ClinVar TSV, fetches sequences
│                                  # on-the-fly, scores, writes brca1_validation_results.csv
│
├── data/
│   ├── brca1_variants.json   # 40 BRCA1 SNVs (20 Pathogenic, 20 Benign) from ClinVar + UCSC
│   │                         # Fields: name, clinical_significance, chromosome, position,
│   │                         #         ref_allele, alt_allele, ref_seq (1001bp), alt_seq (1001bp)
│   └── brca1_scored.json     # Same 40 variants with scoring results appended
│                             # Additional fields: ref_ll, alt_ll, delta, interpretation
│
├── map/
│   ├── draw_roc.py           # Simple ROC curve — reads from brca1_scored.json
│   │                         # Output: map/brca1_roc.png
│   ├── draw_roc_dashboard.py # Dashboard layout: stat cards + ROC + delta histogram
│   │                         # Output: map/brca1_roc_dashboard.png
│   ├── brca1_roc.png         # ROC curve plot
│   └── brca1_roc_dashboard.png  # Full dashboard plot (AUC 0.92)
│
├── tests/
│   ├── unit/
│   │   ├── test_scorer.py    # Tests for sequence fetch and ClinVar loading (no Modal calls)
│   │   ├── test_modal_app.py # Unit tests for compute_log_likelihood and validate_sequence
│   │   └── test_app.py       # Gradio frontend unit tests
│   └── integration/
│       ├── test_endpoint.py  # End-to-end test: hits live Modal endpoint
│       └── test_brca1_validation.py  # End-to-end: scores variants, checks Pathogenic < Benign mean delta
│
├── requirements.txt          # Local deps: gradio, requests, modal (NOT evo2/torch)
└── runtime.txt               # Python version for HuggingFace Spaces
```

### Key data flows

**Scoring a variant:**
```
brca1_variants.json
  → scripts/score_brca1_from_json.py
  → scripts/client.py → POST → src/modal_app.py (Modal GPU, Evo2 7B)
  → brca1_scored.json
```

**Building the variant dataset:**
```
ClinVar API (esearch/esummary) → 40 BRCA1 SNVs
UCSC hg38 (getData/sequence)   → 1001bp flanking sequences, ref base verified
  → data/brca1_variants.json
```

**Visualisation:**
```
data/brca1_scored.json → map/draw_roc_dashboard.py → map/brca1_roc_dashboard.png
```

---

## Validation results (BRCA1, 40 variants)

| Class | n | Mean delta |
|-------|---|-----------|
| Pathogenic | 20 | −38.1 |
| Benign | 20 | −3.5 |

AUC: **0.92**

Variants: 20 Pathogenic nonsense/splice SNVs, 20 Benign intronic/synonymous SNVs — all GRCh38, 1001bp context (500bp flanking each side).

---

## Limitations

- **Scores are relative, not absolute.** Log-likelihood is sequence-context-dependent; comparing scores across very different sequences is not meaningful.
- **No strand awareness.** The model scores the sequence as given. Reverse-complement effects are not considered.
- **No population frequency data.** Common variants can score as deleterious; rare variants can score as neutral.
- **Context window matters.** For long sequences, the local nucleotide context around a variant drives the score. Results may differ depending on how much flanking sequence is included.
- **Not calibrated.** Delta log-likelihood reflects sequence-level probability, not functional impact directly. Synonymous variants may still receive non-zero deltas. Interpretation thresholds are approximate and not clinically validated.
- **Research use only.** AskEvo2 is a research prototype. It has not been validated for clinical, diagnostic, or regulatory purposes.
