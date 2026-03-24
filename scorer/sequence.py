"""
Fetch sequence context from UCSC hg38 for a given variant and build
ref_seq / alt_seq ready to pass to score_variant().
"""

import re

import requests

UCSC_API = "https://api.genome.ucsc.edu/getData/sequence"
CONTEXT = 500  # bp of flanking sequence on each side
VALID_BASE = re.compile(r"^[ACGT]$")


class RefMismatchError(ValueError):
    pass


def fetch_sequence_context(
    chrom: str,
    position: int,
    ref_allele: str,
    alt_allele: str,
    context: int = CONTEXT,
) -> tuple[str, str]:
    """
    Fetch hg38 sequence context around a SNV from UCSC.

    Parameters
    ----------
    chrom : str
        Chromosome, e.g. "chr17". Prefixed with "chr" if absent.
    position : int
        1-based genomic position of the variant.
    ref_allele : str
        Expected reference base (single nucleotide, uppercase).
    alt_allele : str
        Alternate base (single nucleotide, uppercase).
    context : int
        Flanking bp on each side of the variant (default 500).

    Returns
    -------
    ref_seq, alt_seq : (str, str)
        Uppercase sequences ready for score_variant().

    Raises
    ------
    RefMismatchError
        If the fetched reference base does not match ref_allele.
    ValueError
        If alleles are not single uppercase ACGT bases.
    requests.HTTPError
        On UCSC API failure.
    """
    ref_allele = ref_allele.strip().upper()
    alt_allele = alt_allele.strip().upper()

    if not VALID_BASE.match(ref_allele):
        raise ValueError(f"ref_allele must be a single ACGT base, got {ref_allele!r}")
    if not VALID_BASE.match(alt_allele):
        raise ValueError(f"alt_allele must be a single ACGT base, got {alt_allele!r}")

    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"

    # UCSC API uses 0-based half-open coordinates
    start = position - 1 - context
    end = position + context  # includes the variant base

    if start < 0:
        start = 0

    resp = requests.get(
        UCSC_API,
        params={"genome": "hg38", "chrom": chrom, "start": start, "end": end},
        timeout=30,
    )
    resp.raise_for_status()

    dna = resp.json().get("dna", "")
    if not dna:
        raise ValueError(f"UCSC returned empty sequence for {chrom}:{start}-{end}")

    seq = dna.upper()

    # The variant base is at index (position - 1 - start) within the fetched region
    var_idx = (position - 1) - start
    if var_idx < 0 or var_idx >= len(seq):
        raise ValueError(
            f"Variant index {var_idx} out of range for fetched sequence of length {len(seq)}"
        )

    fetched_ref = seq[var_idx]
    if fetched_ref != ref_allele:
        raise RefMismatchError(
            f"Reference mismatch at {chrom}:{position}: "
            f"ClinVar says {ref_allele!r}, UCSC hg38 has {fetched_ref!r}"
        )

    ref_seq = seq
    alt_seq = seq[:var_idx] + alt_allele + seq[var_idx + 1 :]

    return ref_seq, alt_seq
