#!/usr/bin/env python3
"""
Fetch 100 Pathogenic and 100 Benign BRCA1 SNVs from ClinVar + UCSC flanking sequences.
Filters strictly for exact germline_classification.description.
"""

import json
import time
import requests
from urllib.parse import quote

OUTPUT_PATH = "/Users/dgupta/code/portfolio/EvoMine/ask_evo2/data/brca1_variants_200.json"
CHROM = "chr17"
FLANK = 500  # 500bp each side → 1001bp total

HEADERS = {"User-Agent": "Mozilla/5.0 (research use)"}

def clinvar_esearch(term, retmax=500):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "clinvar",
        "term": term,
        "retmax": retmax,
        "retmode": "json",
    }
    r = requests.get(url, params=params, timeout=30, headers=HEADERS)
    r.raise_for_status()
    data = r.json()
    ids = data["esearchresult"]["idlist"]
    count = data["esearchresult"]["count"]
    return ids, int(count)

def clinvar_esummary_batch(ids):
    """Fetch esummary for up to 100 IDs."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "clinvar",
        "id": ",".join(str(i) for i in ids),
        "retmode": "json",
    }
    r = requests.get(url, params=params, timeout=60, headers=HEADERS)
    r.raise_for_status()
    return r.json()

def ucsc_sequence(chrom, start, end):
    """Fetch sequence from UCSC. start/end are 0-based half-open."""
    url = "https://api.genome.ucsc.edu/getData/sequence"
    params = {
        "genome": "hg38",
        "chrom": chrom,
        "start": start,
        "end": end,
    }
    r = requests.get(url, params=params, timeout=30, headers=HEADERS)
    r.raise_for_status()
    data = r.json()
    return data.get("dna", "").upper()

def parse_spdi(spdi_str):
    """
    Parse SPDI string like 'NC_000017.11:43063902:G:A'
    Returns (seq_acc, pos_0based, ref, alt) or None.
    """
    if not spdi_str:
        return None
    parts = spdi_str.split(":")
    if len(parts) != 4:
        return None
    seq_acc = parts[0]
    try:
        pos_0based = int(parts[1])
    except ValueError:
        return None
    ref_allele = parts[2].upper()
    alt_allele = parts[3].upper()
    return seq_acc, pos_0based, ref_allele, alt_allele

def parse_variant(uid, summary):
    """
    Extract variant info from esummary record.
    Returns dict or None if not valid SNV.
    """
    # Get variation_set - canonical_spdi lives here
    variation_set = summary.get("variation_set", [])
    if not variation_set:
        return None

    vs = variation_set[0]

    # Check variant type
    vtype = vs.get("variant_type", "").lower()
    if "single nucleotide" not in vtype:
        return None

    # Get SPDI
    spdi_str = vs.get("canonical_spdi", "")
    parsed = parse_spdi(spdi_str)
    if parsed is None:
        return None

    seq_acc, pos_0based, ref_allele, alt_allele = parsed

    # Must be chr17
    if not seq_acc.startswith("NC_000017"):
        return None

    # Single base SNV only
    if len(ref_allele) != 1 or len(alt_allele) != 1:
        return None
    if ref_allele not in "ACGT" or alt_allele not in "ACGT":
        return None
    if ref_allele == alt_allele:
        return None

    # 1-based position
    position_1based = pos_0based + 1

    # Name
    name = summary.get("title", f"VCV{uid}")

    return {
        "uid": uid,
        "name": name,
        "chromosome": CHROM,
        "position": position_1based,
        "pos_0based": pos_0based,
        "ref_allele": ref_allele,
        "alt_allele": alt_allele,
    }

def fetch_and_verify_sequence(variant):
    """
    Fetch 1001bp flanking from UCSC, verify ref base at index 500.
    Returns (ref_seq, alt_seq) or (None, None).
    start = pos_0based - 500 (0-based)
    end = pos_0based + 501 (exclusive)
    → sequence length = 1001, ref base at index 500
    """
    pos_0based = variant["pos_0based"]
    start = pos_0based - FLANK       # pos_0based - 500
    end = pos_0based + FLANK + 1     # pos_0based + 501

    if start < 0:
        print(f"  Skipping uid {variant['uid']}: start < 0")
        return None, None

    try:
        seq = ucsc_sequence(CHROM, start, end)
    except Exception as e:
        print(f"  UCSC error for {variant['uid']}: {e}")
        return None, None

    if len(seq) != 1001:
        print(f"  Unexpected seq length {len(seq)} for uid {variant['uid']}")
        return None, None

    ref_at_500 = seq[500]
    if ref_at_500 != variant["ref_allele"]:
        print(f"  Ref mismatch uid {variant['uid']}: expected {variant['ref_allele']} got {ref_at_500} at pos {variant['position']}")
        return None, None

    alt_seq = seq[:500] + variant["alt_allele"] + seq[501:]
    return seq, alt_seq

def get_all_candidate_ids(search_terms, retmax=500):
    """Collect unique IDs from multiple search terms."""
    all_ids = []
    seen = set()
    for term in search_terms:
        ids, count = clinvar_esearch(term, retmax=retmax)
        print(f"  Search '{term[:60]}...' → count={count}, returned {len(ids)} IDs")
        for id_ in ids:
            if id_ not in seen:
                seen.add(id_)
                all_ids.append(id_)
        time.sleep(0.4)
    return all_ids

def fetch_variants_for_clinsig(search_terms, exact_clinsig, needed=100):
    """
    Collect enough candidates, filter for exact classification, fetch sequences.
    """
    print(f"\n{'='*60}")
    print(f"Collecting {needed} {exact_clinsig} BRCA1 SNVs")
    print(f"{'='*60}")

    all_ids = get_all_candidate_ids(search_terms, retmax=500)
    print(f"Total unique candidates: {len(all_ids)}")

    valid_variants = []
    processed = 0
    batch_size = 100
    dedup_positions = set()

    for batch_start in range(0, len(all_ids), batch_size):
        if len(valid_variants) >= needed:
            break

        batch = all_ids[batch_start:batch_start + batch_size]
        print(f"\nBatch {batch_start//batch_size + 1} (IDs {batch_start}–{batch_start+len(batch)-1}), "
              f"valid so far: {len(valid_variants)}/{needed}")

        try:
            data = clinvar_esummary_batch(batch)
        except Exception as e:
            print(f"  esummary error: {e}")
            time.sleep(2)
            continue

        time.sleep(0.4)

        result = data.get("result", {})
        uids = result.get("uids", batch)

        for uid in uids:
            if len(valid_variants) >= needed:
                break
            if uid == "uids":
                continue

            summary = result.get(str(uid), result.get(uid, {}))
            if not summary:
                continue

            processed += 1

            # Strict classification check
            germ = summary.get("germline_classification", {})
            clinsig = germ.get("description", "")
            if clinsig != exact_clinsig:
                continue

            # Parse variant
            variant = parse_variant(uid, summary)
            if variant is None:
                continue

            # Dedup by position
            pos_key = (variant["pos_0based"], variant["ref_allele"], variant["alt_allele"])
            if pos_key in dedup_positions:
                print(f"  Skip dup position: uid={uid} pos={variant['position']}")
                continue
            dedup_positions.add(pos_key)

            print(f"  [{len(valid_variants)+1}/{needed}] {exact_clinsig}: uid={uid}, "
                  f"chr17:{variant['position']} {variant['ref_allele']}>{variant['alt_allele']}")

            # Fetch and verify UCSC sequence
            ref_seq, alt_seq = fetch_and_verify_sequence(variant)
            time.sleep(0.25)

            if ref_seq is None:
                continue

            record = {
                "name": variant["name"],
                "clinical_significance": exact_clinsig,
                "chromosome": CHROM,
                "position": variant["position"],
                "ref_allele": variant["ref_allele"],
                "alt_allele": variant["alt_allele"],
                "ref_seq": ref_seq,
                "alt_seq": alt_seq,
            }
            valid_variants.append(record)

    print(f"\nProcessed {processed} records, found {len(valid_variants)} valid {exact_clinsig}")
    return valid_variants

def main():
    # Search terms for Pathogenic - multiple strategies to get enough candidates
    path_terms = [
        'BRCA1[gene] AND "pathogenic"[clinsig] AND "reviewed by expert panel"[review_status]',
        'BRCA1[gene] AND "pathogenic"[clinsig] AND "practice guideline"[review_status]',
        'BRCA1[gene] AND pathogenic[clinsig] AND "criteria provided, single submitter"[review_status]',
        'BRCA1[gene] AND pathogenic[clinsig]',
    ]

    benign_terms = [
        'BRCA1[gene] AND "benign"[clinsig] AND "reviewed by expert panel"[review_status]',
        'BRCA1[gene] AND "benign"[clinsig] AND "practice guideline"[review_status]',
        'BRCA1[gene] AND benign[clinsig] AND "criteria provided, single submitter"[review_status]',
        'BRCA1[gene] AND benign[clinsig]',
    ]

    pathogenic_variants = fetch_variants_for_clinsig(path_terms, "Pathogenic", needed=100)
    benign_variants = fetch_variants_for_clinsig(benign_terms, "Benign", needed=100)

    all_variants = pathogenic_variants + benign_variants

    with open(OUTPUT_PATH, "w") as f:
        json.dump(all_variants, f, indent=2)

    print(f"\n{'='*60}")
    print(f"FINAL SUMMARY")
    print(f"{'='*60}")
    print(f"Pathogenic variants: {len(pathogenic_variants)}")
    print(f"Benign variants:     {len(benign_variants)}")
    print(f"Total:               {len(all_variants)}")
    print(f"Output: {OUTPUT_PATH}")

    if len(pathogenic_variants) < 100:
        print(f"WARNING: Only {len(pathogenic_variants)} Pathogenic variants (needed 100)")
    if len(benign_variants) < 100:
        print(f"WARNING: Only {len(benign_variants)} Benign variants (needed 100)")

if __name__ == "__main__":
    main()
