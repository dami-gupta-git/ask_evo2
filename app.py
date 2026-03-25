import csv
import io
import os
import re
import tempfile
import time
from datetime import datetime, timezone

import gradio as gr
import requests

MODAL_ENDPOINT_URL = "https://dami-gupta-git--askevo2-scorer-score-variant.modal.run"

VALID_SEQ = re.compile(r"^[ACGT]+$")
MAX_SEQ_LEN = 10_000
MIN_SEQ_LEN = 10
MAX_BATCH = 5

BRCA1_REF = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAA"
BRCA1_ALT = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAG"

if "<your-modal-app>" in MODAL_ENDPOINT_URL:
    raise RuntimeError(
        "Set MODAL_ENDPOINT_URL in app.py to your deployed Modal endpoint URL "
        "before running (run `modal deploy modal_app.py` to get the URL)."
    )


def validate_inputs(ref_seq: str, alt_seq: str):
    for name, seq in [("Reference", ref_seq), ("Alternate", alt_seq)]:
        seq = seq.strip().upper()
        if not seq:
            return None, None, f"Error: {name} sequence is empty."
        if not VALID_SEQ.match(seq):
            return None, None, (
                f"Error: {name} sequence contains invalid characters. "
                "Only A, C, G, T are allowed."
            )
        if len(seq) < MIN_SEQ_LEN:
            return None, None, (
                f"Error: {name} sequence must be at least {MIN_SEQ_LEN} bp."
            )
        if len(seq) > MAX_SEQ_LEN:
            return None, None, (
                f"Error: {name} sequence exceeds {MAX_SEQ_LEN} bp limit."
            )
    return ref_seq.strip().upper(), alt_seq.strip().upper(), None


def predict(ref_seq: str, alt_seq: str, request: gr.Request):
    t0 = time.time()
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    print(f"[{now}] ip={request.client.host} ref={ref_seq.strip()[:50]!r}")
    ref_clean, alt_clean, err = validate_inputs(ref_seq, alt_seq)
    if err:
        return gr.update(value=f'<p style="color:#b91c1c;margin:0">⚠ {err}</p>'), "", "", ""

    try:
        resp = requests.post(
            MODAL_ENDPOINT_URL,
            json={"ref_sequence": ref_clean, "alt_sequence": alt_clean},
            timeout=300,
        )
        resp.raise_for_status()
    except requests.exceptions.Timeout:
        return (
            gr.update(value='<p style="color:#b91c1c;margin:0">⚠ Request timed out. The model may be cold-starting — retry in 30s.</p>'),
            "", "", "",
        )
    except requests.exceptions.RequestException as e:
        return gr.update(value=f'<p style="color:#b91c1c;margin:0">⚠ {e}</p>'), "", "", ""

    data = resp.json()
    elapsed = time.time() - t0
    print(f"[done] elapsed={elapsed:.2f}s delta={data['delta']}")
    return (
        gr.update(value=""),
        f"{data['ref_ll']:.4f}",
        f"{data['alt_ll']:.4f}",
        f"{data['delta']:.4f}",
    )


def predict_batch(file, request: gr.Request):
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    if file is None:
        return gr.update(value='<p style="color:#b91c1c;margin:0">⚠ No file uploaded.</p>'), None, gr.update(value=None)

    try:
        content = open(file, "r").read()
        reader = csv.DictReader(io.StringIO(content))
    except Exception as e:
        return gr.update(value=f'<p style="color:#b91c1c;margin:0">⚠ Could not read file: {e}</p>'), None, gr.update(value=None)

    if "ref_sequence" not in (reader.fieldnames or []) or "alt_sequence" not in (reader.fieldnames or []):
        return gr.update(value='<p style="color:#b91c1c;margin:0">⚠ CSV must have columns: ref_sequence, alt_sequence</p>'), None, gr.update(value=None)

    rows = list(reader)
    if len(rows) == 0:
        return gr.update(value='<p style="color:#b91c1c;margin:0">⚠ CSV is empty.</p>'), None, gr.update(value=None)
    if len(rows) > MAX_BATCH:
        return gr.update(value=f'<p style="color:#b91c1c;margin:0">⚠ Too many rows. Maximum is {MAX_BATCH} sequence pairs.</p>'), None, gr.update(value=None)

    print(f"[{now}] ip={request.client.host} batch={len(rows)}")

    results = []
    for i, row in enumerate(rows, start=1):
        ref_raw = row.get("ref_sequence", "")
        alt_raw = row.get("alt_sequence", "")
        ref_clean, alt_clean, err = validate_inputs(ref_raw, alt_raw)
        if err:
            return gr.update(value=f'<p style="color:#b91c1c;margin:0">⚠ Row {i}: {err}</p>'), None, gr.update(value=None)

        try:
            resp = requests.post(
                MODAL_ENDPOINT_URL,
                json={"ref_sequence": ref_clean, "alt_sequence": alt_clean},
                timeout=300,
            )
            resp.raise_for_status()
        except requests.exceptions.Timeout:
            return gr.update(value=f'<p style="color:#b91c1c;margin:0">⚠ Row {i}: Request timed out. The model may be cold-starting — retry in 30s.</p>'), None, gr.update(value=None)
        except requests.exceptions.RequestException as e:
            return gr.update(value=f'<p style="color:#b91c1c;margin:0">⚠ Row {i}: {e}</p>'), None, gr.update(value=None)

        data = resp.json()
        results.append({
            "ref_sequence": ref_clean,
            "alt_sequence": alt_clean,
            "ref_ll": f"{data['ref_ll']:.4f}",
            "alt_ll": f"{data['alt_ll']:.4f}",
            "delta": f"{data['delta']:.4f}",
        })
        print(f"[batch row {i}] delta={data['delta']:.4f}")

    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, newline="")
    writer = csv.DictWriter(tmp, fieldnames=["ref_sequence", "alt_sequence", "ref_ll", "alt_ll", "delta"])
    writer.writeheader()
    writer.writerows(results)
    tmp.close()

    table = [
        [
            i,
            r["ref_sequence"][:20] + ("…" if len(r["ref_sequence"]) > 20 else ""),
            r["alt_sequence"][:20] + ("…" if len(r["alt_sequence"]) > 20 else ""),
            r["ref_ll"],
            r["alt_ll"],
            r["delta"],
        ]
        for i, r in enumerate(results, start=1)
    ]
    return gr.update(value=""), table, gr.update(value=tmp.name)


css = """
#score-btn { background: #60a5fa !important; border-color: #60a5fa !important; }
#score-btn:hover { background: #3b82f6 !important; border-color: #3b82f6 !important; }
#batch-btn { background: #60a5fa !important; border-color: #60a5fa !important; }
#batch-btn:hover { background: #3b82f6 !important; border-color: #3b82f6 !important; }
#load-example-btn { background: #dcffeb !important; border-color: #dcffeb !important; color: #14532d !important; }
#load-example-btn:hover { background: #bce8cc !important; border-color: #bce8cc !important; }
#download-example-btn { background: #dcffeb !important; border-color: #dcffeb !important; color: #14532d !important; width: 50% !important; }
#download-example-btn:hover { background: #bce8cc !important; border-color: #bce8cc !important; }
#download-results-btn { background: #60a5fa !important; border-color: #60a5fa !important; }
#download-results-btn:hover { background: #3b82f6 !important; border-color: #3b82f6 !important; }
.compact-error { display: contents !important; }
#error-out { display: none; color: #b91c1c; margin: 0; padding: 0; }
#error-out:not(:empty) { display: block !important; }
#batch-error-out { display: none; color: #b91c1c; margin: 0; padding: 0; }
#batch-error-out:not(:empty) { display: block !important; }
"""

with gr.Blocks(title="AskEvo2", css=css) as demo:
    gr.Markdown("# AskEvo2 — Zero-shot variant effect scoring powered by Evo2 7B")
    gr.Markdown("Arc Institute Evo2 (https://arcinstitute.org/tools/evo)")
    gr.Markdown("> **For research use only. Not a clinical tool.**")

    with gr.Tabs():
        with gr.Tab("Single Variant"):
            gr.Markdown("Enter reference and alternate DNA sequences to see how Evo2 scores the difference.")

            with gr.Row(elem_classes="compact-error"):
                error_out = gr.HTML(value="", elem_id="error-out")

            with gr.Row():
                ref_input = gr.Textbox(
                    label="Reference Sequence",
                    placeholder="Enter reference DNA sequence (A, C, G, T only)",
                    lines=4,
                )
                alt_input = gr.Textbox(
                    label="Alternate Sequence",
                    placeholder="Enter alternate DNA sequence (A, C, G, T only)",
                    lines=4,
                )

            with gr.Row():
                run_btn = gr.Button("Score Variant", variant="primary", elem_id="score-btn")
                clear_btn = gr.Button("Clear")

            with gr.Row():
                ref_ll_out = gr.Textbox(label="Reference Log-Likelihood")
                alt_ll_out = gr.Textbox(label="Alternate Log-Likelihood")
                delta_out = gr.Textbox(label="Delta (Alt − Ref)")

            with gr.Accordion("► How to read the scores", open=False):
                gr.Markdown(
                    "Log-likelihood measures how probable a sequence is under Evo2's learned model of genomic DNA.\n\n"
                    "**Δ (Alt − Ref):** Difference in log-likelihood between the alternate and reference sequences. "
                    "More negative values are associated with pathogenic variants.\n\n"
                    "- More negative Δ → Greater likelihood of pathogenicity \n"
                    "- Near zero Δ → Model sees little difference between ref and alt\n\n"
                    "*For research use only — not a clinical prediction.*"
                )

            gr.Markdown("**Example:** BRCA1 synonymous SNV (C→G at position 60)")
            load_btn = gr.Button("Load Example", elem_id="load-example-btn")

            run_btn.click(
                fn=predict,
                inputs=[ref_input, alt_input],
                outputs=[error_out, ref_ll_out, alt_ll_out, delta_out],
            )
            load_btn.click(
                fn=lambda: (BRCA1_REF, BRCA1_ALT),
                outputs=[ref_input, alt_input],
            )
            clear_btn.click(
                fn=lambda: ("", "", "", "", "", ""),
                outputs=[error_out, ref_input, alt_input, ref_ll_out, alt_ll_out, delta_out],
            )

        with gr.Tab("Batch Upload"):
            gr.Markdown("Upload a CSV file with up to 5 sequence pairs. Required columns: `ref_sequence`, `alt_sequence`.")

            with gr.Row(elem_classes="compact-error"):
                batch_error_out = gr.HTML(value="", elem_id="batch-error-out")

            with gr.Row():
                gr.Markdown("**Upload your CSV file** — need a template?")
                gr.DownloadButton(label="Example File", value="example_batch.csv", elem_id="download-example-btn")

            batch_file = gr.File(label="Upload CSV", file_types=[".csv"])

            with gr.Row():
                batch_btn = gr.Button("Score Batch", variant="primary", elem_id="batch-btn")

            batch_results = gr.Dataframe(
                headers=["#", "Ref (preview)", "Alt (preview)", "Ref LL", "Alt LL", "Delta"],
                datatype=["number", "str", "str", "str", "str", "str"],
                label="Results",
                interactive=False,
            )

            with gr.Row():
                batch_clear_btn = gr.Button("Clear")
                download_results_btn = gr.DownloadButton(label="Download Results", value=None, elem_id="download-results-btn")

            batch_btn.click(
                fn=predict_batch,
                inputs=[batch_file],
                outputs=[batch_error_out, batch_results, download_results_btn],
            )
            batch_clear_btn.click(
                fn=lambda: ("", None, None, gr.update(value=None)),
                outputs=[batch_error_out, batch_file, batch_results, download_results_btn],
            )

if __name__ == "__main__":
    demo.launch()
