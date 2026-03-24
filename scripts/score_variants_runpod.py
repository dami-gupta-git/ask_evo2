"""
Score variants from a JSON file directly on RunPod (no Modal).

Loads Evo2 7B once, then scores all variants in a loop.
Outputs a scored JSON file with ref_ll, alt_ll, delta, interpretation.

Usage:
    python score_variants_runpod.py --input data/brca1_variants_200.json --output data/brca1_scored_200.json
    python score_variants_runpod.py --input data/mlh1_variants_200.json  --output data/mlh1_scored_200.json
"""

import argparse
import json
import logging
import math
import warnings
from pathlib import Path

import torch
import torch.nn.functional as F


def interpret(delta: float) -> str:
    if delta < -2.0:
        return "Strongly deleterious"
    elif delta < -1.0:
        return "Likely deleterious"
    elif delta < -0.5:
        return "Possibly deleterious"
    elif delta > 0.5:
        return "Likely neutral"
    return "Uncertain"


def compute_log_likelihood(model, sequence: str) -> float:
    input_ids = torch.tensor(
        model.tokenizer.tokenize(sequence),
        dtype=torch.int,
    ).unsqueeze(0).to("cuda:0")  # (1, N)

    with torch.no_grad():
        result = model(input_ids)
        logits = result
        while isinstance(logits, (tuple, list)):
            logits = logits[0]  # (1, N, vocab)

    log_probs = F.log_softmax(logits.float(), dim=-1)
    token_ids = input_ids[0, 1:].long()  # (N-1,)
    per_token_lp = log_probs[0, :-1, :].gather(
        dim=-1,
        index=token_ids.unsqueeze(-1),
    ).squeeze(-1)  # (N-1,)
    return per_token_lp.sum().item()


def load_model():
    logging.disable(logging.WARNING)
    warnings.filterwarnings("ignore")

    _orig = torch.load
    torch.load = lambda *a, **kw: _orig(*a, **{**kw, "weights_only": False})

    from evo2 import Evo2
    model = Evo2("evo2_7b")
    if hasattr(model, "eval"):
        model.eval()

    logging.disable(logging.NOTSET)
    return model


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  required=True, help="Input variants JSON file")
    parser.add_argument("--output", required=True, help="Output scored JSON file")
    args = parser.parse_args()

    in_path  = Path(args.input)
    out_path = Path(args.output)

    with open(in_path) as f:
        variants = json.load(f)

    print(f"Loaded {len(variants)} variants from {in_path.name}")
    print("Loading Evo2 7B...")
    model = load_model()
    print("Model ready.\n")

    print(f"{'#':<4} {'Sig':<12} {'Name':<50} {'ref_ll':>10} {'alt_ll':>10} {'delta':>10}  interpretation")
    print("-" * 110)

    results = []
    for i, v in enumerate(variants, 1):
        name = v["name"]
        sig  = v["clinical_significance"]
        try:
            ref_ll = compute_log_likelihood(model, v["ref_seq"])
            alt_ll = compute_log_likelihood(model, v["alt_seq"])
            delta  = alt_ll - ref_ll
            interp = interpret(delta)
            record = {
                "name": name,
                "clinical_significance": sig,
                "ref_ll": ref_ll,
                "alt_ll": alt_ll,
                "delta": delta,
                "interpretation": interp,
            }
            print(f"{i:<4} {sig:<12} {name[:49]:<50} {ref_ll:>10.3f} {alt_ll:>10.3f} {delta:>10.4f}  {interp}")
        except RuntimeError as e:
            if "out of memory" in str(e).lower():
                torch.cuda.empty_cache()
                msg = "OOM"
            else:
                msg = str(e)
            record = {
                "name": name,
                "clinical_significance": sig,
                "ref_ll": None,
                "alt_ll": None,
                "delta": None,
                "interpretation": f"ERROR: {msg}",
            }
            print(f"{i:<4} {sig:<12} {name[:49]:<50} {'SKIP':>10} {'':>10} {'':>10}  {msg}")
        except Exception as e:
            record = {
                "name": name,
                "clinical_significance": sig,
                "ref_ll": None,
                "alt_ll": None,
                "delta": None,
                "interpretation": f"ERROR: {e}",
            }
            print(f"{i:<4} {sig:<12} {name[:49]:<50} {'SKIP':>10} {'':>10} {'':>10}  {e}")

        results.append(record)

        # flush output so progress is visible over SSH
        import sys
        sys.stdout.flush()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    scored = [r for r in results if r["delta"] is not None]
    print(f"\nDone. {len(scored)}/{len(results)} scored successfully.")
    print(f"Results written to {out_path}")

    if scored:
        path_deltas   = [r["delta"] for r in scored if r["clinical_significance"] == "Pathogenic"]
        benign_deltas = [r["delta"] for r in scored if r["clinical_significance"] == "Benign"]
        if path_deltas:
            print(f"Pathogenic mean delta: {sum(path_deltas)/len(path_deltas):.3f}")
        if benign_deltas:
            print(f"Benign     mean delta: {sum(benign_deltas)/len(benign_deltas):.3f}")


if __name__ == "__main__":
    main()
