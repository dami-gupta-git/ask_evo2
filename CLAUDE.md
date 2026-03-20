# AskEvo2 — CLAUDE.md

## Project Overview

**AskEvo2** is a bioinformatics tool for zero-shot variant effect scoring using Evo2 7B. It takes a reference and alternate DNA sequence, computes per-sequence log-likelihoods via a teacher-forcing forward pass, and returns a delta log-likelihood score indicating whether the variant is likely deleterious or neutral.

---

## Repository Layout

```
ask_evo2/
├── modal_app.py      # Modal GPU backend: Evo2 inference, log-likelihood, web endpoint
├── app.py            # Gradio frontend: UI, input validation, calls Modal endpoint
├── requirements.txt  # Local-only deps (gradio, requests, modal — NOT evo2/torch)
└── README.md         # Deployment guide, interpretation docs, limitations
```

---

## Infrastructure

- **GPU backend**: Modal (A10G GPU, `container_idle_timeout=300`)
- **Model**: `evo2_7b` via `from evo2 import Evo2`; loaded once per container via `@modal.enter()`
- **Weight caching**: Modal persistent Volume `evo2-weights` mounted at `/weights`; `HF_HOME=/weights` set on the image
- **Frontend**: Gradio, intended to be hosted on HuggingFace Spaces (free CPU tier)

---

## Core Logic

### Log-likelihood (teacher-forcing)
- Tokenize sequence → forward pass through Evo2 → log_softmax on logits → gather per-token log-probs (shift-by-1) → sum
- `compute_log_likelihood(model, sequence)` in `modal_app.py`

### Delta and interpretation
- `delta = alt_ll - ref_ll`
- `delta < -0.5` → "Likely deleterious"
- `delta > 0.5` → "Likely neutral"
- otherwise → "Uncertain"

---

## Coding Conventions

- Sequences are validated as uppercase ACGT only, 10–10,000 bp; lowercase is auto-normalized with `.upper()`.
- No fallback or default values for scores — if computation fails, raise with a clear error.
- `evo2` and `torch` are container-only dependencies; never add them to `requirements.txt`.
- The `torch.load` `weights_only=False` workaround must be applied before `Evo2(...)` is instantiated.
- Use `@modal.fastapi_endpoint` (not deprecated `@modal.web_endpoint`).
- No git in this repo.
