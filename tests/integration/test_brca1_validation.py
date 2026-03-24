"""
End-to-end validation test.

Checks that brca1_validation_results.csv exists with 40 rows and that
mean delta for Pathogenic variants is lower than for Benign variants.

Run after scripts/run_brca1_validation.py has completed:
    pytest tests/integration/test_brca1_validation.py
"""

from pathlib import Path

import pandas as pd
import pytest

RESULTS_CSV = Path(__file__).resolve().parents[2] / "brca1_validation_results.csv"


@pytest.fixture(scope="module")
def results():
    if not RESULTS_CSV.exists():
        pytest.skip(
            f"brca1_validation_results.csv not found at {RESULTS_CSV}. "
            "Run scripts/run_brca1_validation.py first."
        )
    return pd.read_csv(RESULTS_CSV)


def test_end_to_end_brca1_validation(results):
    # 40 rows total
    assert len(results) == 40, f"Expected 40 rows, got {len(results)}"

    # 20 Pathogenic, 20 Benign
    counts = results["clinical_significance"].value_counts()
    assert counts.get("Pathogenic", 0) == 20
    assert counts.get("Benign", 0) == 20

    # Required columns present
    required = {"name", "clinical_significance", "chromosome", "position",
                "ref_allele", "alt_allele", "ref_ll", "alt_ll", "delta", "interpretation"}
    assert required.issubset(set(results.columns))

    # No null scores
    assert results["delta"].notna().all()
    assert results["ref_ll"].notna().all()
    assert results["alt_ll"].notna().all()

    # Mean delta for Pathogenic < Mean delta for Benign
    mean_path = results[results["clinical_significance"] == "Pathogenic"]["delta"].mean()
    mean_benign = results[results["clinical_significance"] == "Benign"]["delta"].mean()
    assert mean_path < mean_benign, (
        f"Expected mean Pathogenic delta ({mean_path:.4f}) < "
        f"mean Benign delta ({mean_benign:.4f})"
    )
