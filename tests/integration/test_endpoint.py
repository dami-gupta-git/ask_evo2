"""
Integration tests — hit the live Modal endpoint.
Run with: pytest tests/integration/
"""
import pytest
import requests

ENDPOINT_URL = "https://dami-gupta-git--askevo2-scorer-score-variant.modal.run"
BRCA1_REF = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAA"
BRCA1_ALT = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAG"


def post(ref, alt):
    return requests.post(
        ENDPOINT_URL,
        json={"ref_sequence": ref, "alt_sequence": alt},
        timeout=300,
    )


def test_brca1_snv_scores():
    resp = post(BRCA1_REF, BRCA1_ALT)
    assert resp.status_code == 200
    data = resp.json()
    assert data["ref_ll"] == -76.0
    assert data["alt_ll"] == -77.0
    assert data["delta"] == -1.0
    assert data["interpretation"] == "Likely deleterious"


def test_identical_sequences_uncertain():
    resp = post(BRCA1_REF, BRCA1_REF)
    assert resp.status_code == 200
    data = resp.json()
    assert data["delta"] == pytest.approx(0.0, abs=0.01)
    assert data["interpretation"] == "Uncertain"


@pytest.mark.parametrize("ref,alt,expected_status", [
    ("", BRCA1_ALT, 422),
    (BRCA1_REF, "ACGTNNN", 422),
    ("ACGT", BRCA1_ALT, 422),
])
def test_invalid_input_returns_error(ref, alt, expected_status):
    resp = post(ref, alt)
    assert resp.status_code == expected_status
