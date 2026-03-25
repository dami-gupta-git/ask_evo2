import csv
import os
import tempfile
import pytest
from unittest.mock import patch, MagicMock

from app import validate_inputs, predict, predict_batch


VALID_REF = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAA"
VALID_ALT = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAG"


# --- validate_inputs ---

@pytest.mark.parametrize("ref,alt,exp_ref,exp_alt", [
    (VALID_REF, VALID_ALT, VALID_REF, VALID_ALT),
    (VALID_REF.lower(), VALID_ALT.lower(), VALID_REF, VALID_ALT),  # normalised
    (f"  {VALID_REF}  ", VALID_ALT, VALID_REF, VALID_ALT),         # stripped
])
def test_validate_inputs_valid(ref, alt, exp_ref, exp_alt):
    r, a, err = validate_inputs(ref, alt)
    assert r == exp_ref
    assert a == exp_alt
    assert err is None


@pytest.mark.parametrize("ref,alt,expected_fragment", [
    ("", VALID_ALT, "empty"),
    (VALID_REF, "", "empty"),
    ("ACGTNNN", VALID_ALT, "invalid characters"),
    (VALID_REF, "ACGTNNN", "invalid characters"),
    ("ACGT", VALID_ALT, "at least 10"),
    ("A" * 10_001, VALID_ALT, "10000 bp"),
])
def test_validate_inputs_invalid(ref, alt, expected_fragment):
    r, a, err = validate_inputs(ref, alt)
    assert r is None
    assert a is None
    assert expected_fragment in err


# --- predict ---

def test_predict_returns_validation_error_on_bad_input():
    ref_ll, alt_ll, delta, interp = predict("ACGTNNN", VALID_ALT)
    assert "invalid characters" in ref_ll
    assert alt_ll == ""
    assert delta == ""
    assert interp == ""


def test_predict_success():
    mock_response = MagicMock()
    mock_response.json.return_value = {
        "ref_ll": -76.0,
        "alt_ll": -77.0,
        "delta": -1.0,
        "interpretation": "Likely deleterious",
    }
    mock_response.raise_for_status = MagicMock()

    with patch("app.requests.post", return_value=mock_response):
        ref_ll, alt_ll, delta, interp = predict(VALID_REF, VALID_ALT)

    assert ref_ll == "-76.0000"
    assert alt_ll == "-77.0000"
    assert delta == "-1.0000"
    assert interp == "Likely deleterious"


def test_predict_timeout():
    import requests as req
    with patch("app.requests.post", side_effect=req.exceptions.Timeout):
        ref_ll, alt_ll, delta, interp = predict(VALID_REF, VALID_ALT)

    assert "timed out" in ref_ll
    assert alt_ll == ""
    assert delta == ""
    assert interp == ""


def test_predict_request_error():
    import requests as req
    with patch("app.requests.post", side_effect=req.exceptions.ConnectionError("refused")):
        ref_ll, alt_ll, delta, interp = predict(VALID_REF, VALID_ALT)

    assert "Error" in ref_ll
    assert alt_ll == ""
    assert delta == ""
    assert interp == ""


# --- predict_batch ---

def test_predict_batch_success():
    rows = [
        {"ref_sequence": VALID_REF, "alt_sequence": VALID_ALT},
        {"ref_sequence": VALID_REF, "alt_sequence": VALID_REF},
    ]
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, newline="")
    writer = csv.DictWriter(tmp, fieldnames=["ref_sequence", "alt_sequence"])
    writer.writeheader()
    writer.writerows(rows)
    tmp.close()

    mock_response = MagicMock()
    mock_response.raise_for_status = MagicMock()
    mock_response.json.side_effect = [
        {"ref_ll": -76.0, "alt_ll": -77.0, "delta": -1.0},
        {"ref_ll": -76.0, "alt_ll": -76.0, "delta": 0.0},
    ]

    try:
        with patch("app.requests.post", return_value=mock_response):
            error_html, table, download_btn = predict_batch(tmp.name)

        assert error_html["value"] == ""
        assert len(table) == 2
        assert table[0][0] == 1
        assert table[0][5] == "-1.0000"
        assert table[1][5] == "0.0000"
        assert download_btn["value"] is not None
    finally:
        os.unlink(tmp.name)
