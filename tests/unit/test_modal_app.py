import pytest

from modal_app import validate_sequence, interpret


# --- validate_sequence ---

@pytest.mark.parametrize("seq,name,expected", [
    ("ACGTACGTAC", "ref", "ACGTACGTAC"),
    ("acgtacgtac", "ref", "ACGTACGTAC"),       # lowercase normalised
    ("  ACGTACGT  ", "ref", "ACGTACGT"),        # whitespace stripped
    ("ACGTACGTACGT", "alt", "ACGTACGTACGT"),
])
def test_validate_sequence_valid(seq, name, expected):
    assert validate_sequence(seq, name) == expected


@pytest.mark.parametrize("seq,name,expected_fragment", [
    ("", "ref", "empty"),
    ("ACGTNNN", "ref", "invalid characters"),
    ("ACGT ACG", "alt", "invalid characters"),
    ("ACGU", "ref", "invalid characters"),
    ("ACGT", "ref", "at least 10"),                          # too short
    ("A" * 10_001, "ref", "exceeds 10000 bp"),
])
def test_validate_sequence_invalid(seq, name, expected_fragment):
    with pytest.raises(ValueError, match=expected_fragment):
        validate_sequence(seq, name)


# --- interpret ---

@pytest.mark.parametrize("delta,expected", [
    (-1.0, "Likely deleterious"),
    (-0.51, "Likely deleterious"),
    (-0.5, "Uncertain"),
    (0.0, "Uncertain"),
    (0.5, "Uncertain"),
    (0.51, "Likely neutral"),
    (1.0, "Likely neutral"),
])
def test_interpret(delta, expected):
    assert interpret(delta) == expected
