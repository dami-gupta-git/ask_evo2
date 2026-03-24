"""
Load and filter ClinVar variant_summary.txt for BRCA1 SNVs on GRCh38.

Download variant_summary.txt from:
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
"""

from __future__ import annotations

import pandas as pd

BRCA1_GENE = "BRCA1"
ASSEMBLY = "GRCh38"
VALID_SIGNIFICANCE = {"Pathogenic", "Benign"}

# Columns we care about in variant_summary.txt
_USECOLS = [
    "Name",
    "ClinicalSignificance",
    "Chromosome",
    "Start",
    "ReferenceAllele",
    "AlternateAllele",
    "GeneSymbol",
    "Assembly",
    "Type",
]


def load_brca1_variants(
    path: str,
    clinical_significance: str | list[str] | None = None,
) -> pd.DataFrame:
    """
    Load ClinVar variant_summary.txt and return BRCA1 SNVs on GRCh38.

    Parameters
    ----------
    path : str
        Path to variant_summary.txt (tab-separated, may be gzipped).
    clinical_significance : str, list[str], or None
        Filter to one or more significance labels (e.g. "Pathogenic", "Benign",
        or ["Pathogenic", "Benign"]). None returns all.

    Returns
    -------
    pd.DataFrame with columns:
        name, clinical_significance, chromosome, position, ref_allele, alt_allele

    Only exact single-value significance labels are included (no "Pathogenic/Likely
    pathogenic" composites unless explicitly requested).
    """
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=_USECOLS,
        low_memory=False,
    )

    # Filter: GRCh38, BRCA1, SNV
    df = df[
        (df["Assembly"] == ASSEMBLY)
        & (df["GeneSymbol"] == BRCA1_GENE)
        & (df["Type"] == "single nucleotide variant")
    ].copy()

    # Normalise alleles — drop rows with missing or multi-base alleles
    df["ReferenceAllele"] = df["ReferenceAllele"].astype(str).str.strip().str.upper()
    df["AlternateAllele"] = df["AlternateAllele"].astype(str).str.strip().str.upper()

    single_base = df["ReferenceAllele"].str.match(r"^[ACGT]$") & df[
        "AlternateAllele"
    ].str.match(r"^[ACGT]$")
    df = df[single_base].copy()

    # Drop rows with missing position
    df = df.dropna(subset=["Start"])
    df["Start"] = df["Start"].astype(int)

    # Significance filter
    if clinical_significance is not None:
        if isinstance(clinical_significance, str):
            clinical_significance = [clinical_significance]
        df = df[df["ClinicalSignificance"].isin(clinical_significance)].copy()

    df = df.rename(
        columns={
            "Name": "name",
            "ClinicalSignificance": "clinical_significance",
            "Chromosome": "chromosome",
            "Start": "position",
            "ReferenceAllele": "ref_allele",
            "AlternateAllele": "alt_allele",
        }
    )

    return df[
        ["name", "clinical_significance", "chromosome", "position", "ref_allele", "alt_allele"]
    ].reset_index(drop=True)
