import re

import gradio as gr
import requests

MODAL_ENDPOINT_URL = "https://dami-gupta-git--askevo2-scorer-score-variant.modal.run"

VALID_SEQ = re.compile(r"^[ACGT]+$")
MAX_SEQ_LEN = 10_000
MIN_SEQ_LEN = 10

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


def predict(ref_seq: str, alt_seq: str):
    ref_clean, alt_clean, err = validate_inputs(ref_seq, alt_seq)
    if err:
        return gr.update(value=err, visible=True), "", "", "", ""



    try:
        resp = requests.post(
            MODAL_ENDPOINT_URL,
            json={"ref_sequence": ref_clean, "alt_sequence": alt_clean},
            timeout=300,
        )
        resp.raise_for_status()
    except requests.exceptions.Timeout:
        return (
            gr.update(value="Error: Request timed out. The model may be cold-starting — retry in 30s.", visible=True),
            "", "", "", "",
        )
    except requests.exceptions.RequestException as e:
        return gr.update(value=f"Error: {e}", visible=True), "", "", "", ""

    data = resp.json()
    return (
        gr.update(value="", visible=False),
        f"{data['ref_ll']:.4f}",
        f"{data['alt_ll']:.4f}",
        f"{data['delta']:.4f}",
        data["interpretation"],
    )


with gr.Blocks(title="AskEvo2") as demo:
    gr.Markdown("# AskEvo2 — Zero-shot variant effect scoring powered by Evo2 7B")
    gr.Markdown("> **For research use only. Not a clinical tool.**")

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

    gr.HTML("<style>#score-btn { background: #1d4ed8 !important; border-color: #1d4ed8 !important; } #score-btn:hover { background: #1e40af !important; border-color: #1e40af !important; } #error-box { color: #7f1d1d; background: #fee2e2; border: 1px solid #f87171; border-radius: 6px; padding: 16px 20px; font-weight: 600; font-size: 1rem; min-height: 56px; }</style>")
    with gr.Row():
        error_box = gr.Markdown(
            value="",
            visible=False,
            elem_id="error-box",
        )

    with gr.Row():
        ref_ll_out = gr.Textbox(label="Reference Log-Likelihood")
        alt_ll_out = gr.Textbox(label="Alternate Log-Likelihood")
        delta_out = gr.Textbox(label="Delta (Alt − Ref)")
        interp_out = gr.Textbox(label="Interpretation")

    gr.Examples(
        examples=[[BRCA1_REF, BRCA1_ALT]],
        inputs=[ref_input, alt_input],
        label="Example: BRCA1 synonymous SNV (C→G at position 60)",
    )

    run_btn.click(
        fn=predict,
        inputs=[ref_input, alt_input],
        outputs=[error_box, ref_ll_out, alt_ll_out, delta_out, interp_out],
    )
    clear_btn.click(
        fn=lambda: (gr.update(value="", visible=False), "", "", "", "", ""),
        outputs=[error_box, ref_input, alt_input, ref_ll_out, alt_ll_out, delta_out, interp_out],
    )

if __name__ == "__main__":
    demo.launch()
