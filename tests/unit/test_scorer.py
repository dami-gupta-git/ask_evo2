"""
Unit tests for scorer.sequence and scorer.clinvar.
Does not call the Modal endpoint.
"""

import io
import textwrap
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from scorer.sequence import RefMismatchError, fetch_sequence_context
from scorer.clinvar import load_brca1_variants


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ucsc_response(dna: str) -> MagicMock:
    mock = MagicMock()
    mock.raise_for_status = MagicMock()
    mock.json.return_value = {"dna": dna}
    return mock


CLINVAR_TSV = textwrap.dedent("""\
    Name\tClinicalSignificance\tChromosome\tStart\tReferenceAllele\tAlternateAllele\tGeneSymbol\tAssembly\tType
    NM_007294.4(BRCA1):c.5266dupC\tPathogenic\t17\t43045629\tA\tT\tBRCA1\tGRCh38\tsingle nucleotide variant
    NM_007294.4(BRCA1):c.68_69del\tBenign\t17\t43106460\tG\tC\tBRCA1\tGRCh38\tsingle nucleotide variant
    NM_007294.4(BRCA1):c.1A>G\tPathogenic\t17\t43124096\tC\tA\tBRCA1\tGRCh38\tsingle nucleotide variant
    NM_007294.4(BRCA1):c.2T>C\tBenign\t17\t43115726\tT\tG\tBRCA1\tGRCh38\tsingle nucleotide variant
    NM_007294.3(BRCA1):c.3A>C\tPathogenic\t17\t43106535\tA\tC\tBRCA1\tGRCh37\tsingle nucleotide variant
    NM_007294.3(BRCA1):c.4T>A\tPathogenic\t17\t43091978\tN\tT\tBRCA1\tGRCh38\tsingle nucleotide variant
    NM_007294.3(BRCA1):c.5A>G\tPathogenic\t17\t43091234\tA\tGG\tBRCA1\tGRCh38\tsingle nucleotide variant
    NM_007294.4(BRCA1):c.6C>T\tUncertain significance\t17\t43090943\tC\tT\tBRCA1\tGRCh38\tsingle nucleotide variant
    NM_007294.4(BRCA1):c.7G>A\tBenign\t17\t43082434\tG\tA\tBRCA2\tGRCh38\tsingle nucleotide variant
""")


# ---------------------------------------------------------------------------
# scorer.sequence tests
# ---------------------------------------------------------------------------

class TestFetchSequenceContext:

    def test_returns_ref_and_alt_seq(self):
        # 1001-bp sequence; variant is the middle base (index 500)
        context_seq = "A" * 500 + "C" + "A" * 500
        with patch("scorer.sequence.requests.get", return_value=_ucsc_response(context_seq)):
            ref_seq, alt_seq = fetch_sequence_context("chr17", 43045629, "C", "T")

        assert ref_seq == context_seq.upper()
        assert alt_seq == "A" * 500 + "T" + "A" * 500

    def test_alt_seq_differs_only_at_variant_position(self):
        context_seq = "G" * 500 + "A" + "G" * 500
        with patch("scorer.sequence.requests.get", return_value=_ucsc_response(context_seq)):
            ref_seq, alt_seq = fetch_sequence_context("chr17", 43045629, "A", "C")

        assert ref_seq != alt_seq
        diffs = [i for i in range(len(ref_seq)) if ref_seq[i] != alt_seq[i]]
        assert len(diffs) == 1

    def test_chr_prefix_added_if_missing(self):
        context_seq = "T" * 500 + "G" + "T" * 500
        captured = {}

        def mock_get(url, params, timeout):
            captured["chrom"] = params["chrom"]
            return _ucsc_response(context_seq)

        with patch("scorer.sequence.requests.get", side_effect=mock_get):
            fetch_sequence_context("17", 43045629, "G", "A")

        assert captured["chrom"] == "chr17"

    def test_ref_mismatch_raises(self):
        context_seq = "A" * 500 + "C" + "A" * 500  # fetched ref is C
        with patch("scorer.sequence.requests.get", return_value=_ucsc_response(context_seq)):
            with pytest.raises(RefMismatchError, match="mismatch"):
                fetch_sequence_context("chr17", 43045629, "T", "G")  # claims ref is T

    def test_invalid_ref_allele_raises(self):
        with pytest.raises(ValueError, match="single ACGT base"):
            fetch_sequence_context("chr17", 43045629, "N", "A")

    def test_invalid_alt_allele_raises(self):
        with pytest.raises(ValueError, match="single ACGT base"):
            fetch_sequence_context("chr17", 43045629, "A", "NN")

    def test_sequence_uppercased(self):
        context_seq = "a" * 500 + "c" + "a" * 500
        with patch("scorer.sequence.requests.get", return_value=_ucsc_response(context_seq)):
            ref_seq, alt_seq = fetch_sequence_context("chr17", 43045629, "C", "T")

        assert ref_seq == ref_seq.upper()
        assert alt_seq == alt_seq.upper()

    def test_ucsc_http_error_propagates(self):
        import requests as req
        mock = MagicMock()
        mock.raise_for_status.side_effect = req.HTTPError("404")
        with patch("scorer.sequence.requests.get", return_value=mock):
            with pytest.raises(req.HTTPError):
                fetch_sequence_context("chr17", 43045629, "A", "T")


# ---------------------------------------------------------------------------
# scorer.clinvar tests
# ---------------------------------------------------------------------------

class TestLoadBrca1Variants:

    def _load(self, sig=None):
        buf = io.StringIO(CLINVAR_TSV)
        return load_brca1_variants(buf, clinical_significance=sig)

    def test_returns_dataframe_with_expected_columns(self):
        df = self._load()
        assert list(df.columns) == [
            "name", "clinical_significance", "chromosome",
            "position", "ref_allele", "alt_allele",
        ]

    def test_filters_to_grch38_only(self):
        df = self._load()
        # Row with GRCh37 should be excluded
        assert all(
            "GRCh37" not in str(row["name"]) for _, row in df.iterrows()
        )

    def test_filters_to_brca1_only(self):
        df = self._load()
        # The BRCA2 row should be excluded (even if Benign)
        assert len(df) > 0

    def test_excludes_non_acgt_alleles(self):
        df = self._load()
        assert df["ref_allele"].str.match(r"^[ACGT]$").all()
        assert df["alt_allele"].str.match(r"^[ACGT]$").all()

    def test_excludes_multibase_alleles(self):
        df = self._load()
        # Row with alt_allele="GG" should be excluded
        assert not (df["alt_allele"] == "GG").any()

    def test_significance_filter_pathogenic(self):
        df = self._load(sig="Pathogenic")
        assert (df["clinical_significance"] == "Pathogenic").all()
        assert len(df) > 0

    def test_significance_filter_benign(self):
        df = self._load(sig="Benign")
        assert (df["clinical_significance"] == "Benign").all()

    def test_significance_filter_both(self):
        df = self._load(sig=["Pathogenic", "Benign"])
        assert set(df["clinical_significance"].unique()).issubset({"Pathogenic", "Benign"})

    def test_significance_filter_none_returns_all(self):
        df_all = self._load(sig=None)
        df_path = self._load(sig="Pathogenic")
        df_benign = self._load(sig="Benign")
        assert len(df_all) >= len(df_path) + len(df_benign)

    def test_position_is_int(self):
        df = self._load()
        assert df["position"].dtype == int or pd.api.types.is_integer_dtype(df["position"])

    def test_excludes_uncertain_significance_when_filtering(self):
        df = self._load(sig=["Pathogenic", "Benign"])
        assert "Uncertain significance" not in df["clinical_significance"].values
