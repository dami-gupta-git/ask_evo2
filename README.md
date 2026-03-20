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

## Limitations

- **Scores are relative, not absolute.** Log-likelihood is sequence-context-dependent; comparing scores across very different sequences is not meaningful.
- **No strand awareness.** The model scores the sequence as given. Reverse-complement effects are not considered.
- **No population frequency data.** Common variants can score as deleterious; rare variants can score as neutral.
- **Context window matters.** For long sequences, the local nucleotide context around a variant drives the score. Results may differ depending on how much flanking sequence is included.
- **Not calibrated.** The −0.5 / +0.5 thresholds are heuristics. Do not use these scores for clinical decision-making.
- **Research use only.** AskEvo2 is a research prototype. It has not been validated for clinical, diagnostic, or regulatory purposes.
