import re

import modal
import torch
import torch.nn.functional as F

app = modal.App("askevo2")

vol = modal.Volume.from_name("evo2-weights", create_if_missing=True)

image = (
    modal.Image.debian_slim()
    .pip_install("evo2", "torch", "fastapi[standard]")
    .env({"HF_HOME": "/weights"})
)

VALID_SEQ = re.compile(r"^[ACGT]+$")
MAX_SEQ_LEN = 10_000
MIN_SEQ_LEN = 10


def validate_sequence(seq: str, name: str) -> str:
    seq = seq.strip().upper()
    if not seq:
        raise ValueError(f"{name} is empty.")
    if not VALID_SEQ.match(seq):
        raise ValueError(
            f"{name} contains invalid characters. Only A, C, G, T are allowed."
        )
    if len(seq) < MIN_SEQ_LEN:
        raise ValueError(f"{name} must be at least {MIN_SEQ_LEN} bp.")
    if len(seq) > MAX_SEQ_LEN:
        raise ValueError(
            f"{name} exceeds {MAX_SEQ_LEN} bp limit (got {len(seq)} bp)."
        )
    return seq


def interpret(delta: float) -> str:
    if delta < -0.5:
        return "Likely deleterious"
    elif delta > 0.5:
        return "Likely neutral"
    return "Uncertain"


def compute_log_likelihood(model, sequence: str) -> float:
    input_ids = torch.tensor(
        model.tokenizer.tokenize(sequence),
        dtype=torch.int,
    ).unsqueeze(0).to("cuda:0")  # (1, N)

    with torch.no_grad():
        logits, _ = model(input_ids)  # (1, N, vocab)

    log_probs = F.log_softmax(logits, dim=-1)
    token_ids = input_ids[0, 1:].long()  # (N-1,)
    per_token_lp = log_probs[0, :-1, :].gather(
        dim=-1,
        index=token_ids.unsqueeze(-1),
    ).squeeze(-1)  # (N-1,)
    return per_token_lp.sum().item()


@app.cls(
    gpu="A10G",
    volumes={"/weights": vol},
    scaledown_window=300,
    image=image,
)
class Scorer:
    @modal.enter()
    def load_model(self):
        import torch as _torch
        _orig = _torch.load
        _torch.load = lambda *a, **kw: _orig(*a, **{**kw, "weights_only": False})

        from evo2 import Evo2
        self.model = Evo2("evo2_7b")
        if hasattr(self.model, "eval"):
            self.model.eval()

    @modal.fastapi_endpoint(method="POST")
    def score_variant(self, body: dict) -> dict:
        ref_seq = validate_sequence(body.get("ref_sequence", ""), "ref_sequence")
        alt_seq = validate_sequence(body.get("alt_sequence", ""), "alt_sequence")

        try:
            ref_ll = compute_log_likelihood(self.model, ref_seq)
            alt_ll = compute_log_likelihood(self.model, alt_seq)
        except RuntimeError as e:
            if "out of memory" in str(e).lower():
                raise ValueError("GPU OOM — try a shorter sequence.")
            raise

        delta = alt_ll - ref_ll
        return {
            "ref_ll": ref_ll,
            "alt_ll": alt_ll,
            "delta": delta,
            "interpretation": interpret(delta),
        }
