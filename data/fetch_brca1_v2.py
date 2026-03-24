#!/usr/bin/env python3
"""
Fetch 100 Pathogenic and 100 Benign BRCA1 SNVs from ClinVar + UCSC flanking sequences.
"""

import json
import time
import requests

OUTPUT_PATH = "/Users/dgupta/code/portfolio/EvoMine/ask_evo2/data/brca1_variants_200.json"
CHROM = "chr17"
FLANK = 500  # 500bp each side → 1001bp total

def clinvar_esearch(term, retmax=500):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "clinvar",
        "term": term,
        "retmax": retmax,
        "retmode": "json",
    }
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    return data["esearchresult"]["idlist"]

def clinvar_esummary_batch(ids):
    """Fetch esummary for up to 100 IDs."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "clinvar",
        "id": ",".join(ids),
        "retmode": "json",
    }
    r = requests.get(url, params=params, timeout=60)
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
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    return data.get("dna", "").upper()

def parse_variant(uid, summary):
    """
    Extract variant info from esummary record.
    Returns dict or None if not valid.
    """
    # Check exact germline classification
    germ = summary.get("germline_classification", {})
    clinsig = germ.get("description", "")

    # Check variation_set for the name
    variation_set = summary.get("variation_set", [])
    name = summary.get("title", f"VCV{uid}")

    # Get canonical SPDI
    spdi_list = summary.get("canonical_spdi", [])
    if not spdi_list:
        return None

    # canonical_spdi may be a list or a single dict
    if isinstance(spdi_list, list):
        spdi = spdi_list[0] if spdi_list else None
    else:
        spdi = spdi_list

    if not spdi:
        return None

    # Extract SPDI fields
    # Format: NC_000017.11:position:deleted:inserted
    if isinstance(spdi, str):
        # Parse string format: "NC_000017.11:43045757:C:A"
        parts = spdi.split(":")
        if len(parts) != 4:
            return None
        seq_acc = parts[0]
        try:
            pos_0based = int(parts[1])
        except ValueError:
            return None
        ref_allele = parts[2].upper()
        alt_allele = parts[3].upper()
    elif isinstance(spdi, dict):
        seq_acc = spdi.get("sequence", "")
        try:
            pos_0based = int(spdi.get("position", -1))
        except (ValueError, TypeError):
            return None
        ref_allele = spdi.get("deleted_sequence", "").upper()
        alt_allele = spdi.get("inserted_sequence", "").upper()
    else:
        return None

    # Must be on chr17
    if not seq_acc.startswith("NC_000017"):
        return None

    # Single base SNV only
    if len(ref_allele) != 1 or len(alt_allele) != 1:
        return None
    if ref_allele not in "ACGT" or alt_allele not in "ACGT":
        return None
    if ref_allele == alt_allele:
        return None

    # SPDI position is 0-based; 1-based position = pos_0based + 1
    position_1based = pos_0based + 1

    return {
        "uid": uid,
        "name": name,
        "clinical_significance": clinsig,
        "chromosome": CHROM,
        "position": position_1based,
        "pos_0based": pos_0based,
        "ref_allele": ref_allele,
        "alt_allele": alt_allele,
    }

def fetch_and_verify_sequence(variant):
    """
    Fetch 1001bp flanking from UCSC, verify ref base at index 500.
    Returns (ref_seq, alt_seq) or (None, None) if mismatch.
    """
    pos_0based = variant["pos_0based"]
    start = pos_0based - FLANK  # 0-based start for UCSC (position - 501 in 0-based = pos_0based - 500)
    end = pos_0based + FLANK + 1  # exclusive end → covers pos_0based + 500

    # UCSC: start=position-501 end=position+500 (from task description)
    # Task says: start=position-501, end=position+500 where position is 1-based
    # So: start = pos_1based - 501 = pos_0based - 500
    # end = pos_1based + 500 = pos_0based + 501
    # That's 1001 bases total: [pos_0based-500 .. pos_0based+500] inclusive
    # In 0-based half-open: start=pos_0based-500, end=pos_0based+501
    start = pos_0based - FLANK  # = pos_0based - 500
    end = pos_0based + FLANK + 1  # = pos_0based + 501

    if start < 0:
        return None, None

    try:
        seq = ucsc_sequence(CHROM, start, end)
    except Exception as e:
        print(f"  UCSC error for {variant['uid']}: {e}")
        return None, None

    if len(seq) != 1001:
        print(f"  Unexpected sequence length {len(seq)} for uid {variant['uid']}")
        return None, None

    # Verify ref base at index 500
    ref_at_500 = seq[500]
    if ref_at_500 != variant["ref_allele"]:
        print(f"  Ref mismatch uid {variant['uid']}: expected {variant['ref_allele']} got {ref_at_500}")
        return None, None

    # Build alt_seq
    alt_seq = seq[:500] + variant["alt_allele"] + seq[501:]

    return seq, alt_seq

def fetch_candidates(search_term, exact_clinsig, needed=100):
    """
    Search ClinVar with broad term, fetch esummary in batches,
    filter for exact clinsig, and return valid variants with sequences.
    """
    print(f"\n=== Searching for {exact_clinsig} variants ===")

    # Get candidate IDs - fetch more to have enough after filtering
    all_ids = clinvar_esearch(search_term, retmax=500)
    print(f"Found {len(all_ids)} candidate IDs from esearch")
    time.sleep(0.34)

    valid_variants = []
    batch_size = 100

    for batch_start in range(0, len(all_ids), batch_size):
        if len(valid_variants) >= needed:
            break

        batch = all_ids[batch_start:batch_start + batch_size]
        print(f"Processing batch {batch_start//batch_size + 1} (IDs {batch_start}-{batch_start+len(batch)-1})")

        try:
            data = clinvar_esummary_batch(batch)
        except Exception as e:
            print(f"  esummary error: {e}")
            time.sleep(1)
            continue

        time.sleep(0.34)

        result = data.get("result", {})
        uids = result.get("uids", batch)

        for uid in uids:
            if uid == "uids":
                continue
            summary = result.get(uid, {})

            # Check exact classification
            germ = summary.get("germline_classification", {})
            clinsig = germ.get("description", "")

            if clinsig != exact_clinsig:
                continue

            # Parse variant
            variant = parse_variant(uid, summary)
            if variant is None:
                continue

            print(f"  Found {exact_clinsig}: uid={uid}, chr17:{variant['position']} {variant['ref_allele']}>{variant['alt_allele']}")

            # Fetch UCSC sequence and verify
            ref_seq, alt_seq = fetch_and_verify_sequence(variant)
            time.sleep(0.2)

            if ref_seq is None:
                continue

            # Build final record
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
            print(f"  => Added record #{len(valid_variants)}/{needed}")

            if len(valid_variants) >= needed:
                break

    return valid_variants

def main():
    # Search terms - broad enough to get many candidates
    path_term = "BRCA1[gene]+AND+Pathogenic[clinsig]+AND+%22single+nucleotide+variant%22[Type]+AND+GRCh38[Assembly]"
    benign_term = "BRCA1[gene]+AND+Benign[clinsig]+AND+%22single+nucleotide+variant%22[Type]+AND+GRCh38[Assembly]"

    # Also try broader terms to get more candidates
    path_term2 = 'BRCA1[Gene Name] AND "Pathogenic"[Clinical significance] AND "single nucleotide variant"[Variant type]'
    benign_term2 = 'BRCA1[Gene Name] AND "Benign"[Clinical significance] AND "single nucleotide variant"[Variant type]'

    pathogenic_variants = fetch_candidates(path_term2, "Pathogenic", needed=100)
    print(f"\nPathogenic: collected {len(pathogenic_variants)} valid variants")

    if len(pathogenic_variants) < 100:
        print("Not enough pathogenic variants from first search, trying additional sources...")

    benign_variants = fetch_candidates(benign_term2, "Benign", needed=100)
    print(f"\nBenign: collected {len(benign_variants)} valid variants")

    all_variants = pathogenic_variants + benign_variants

    # Save output
    with open(OUTPUT_PATH, "w") as f:
        json.dump(all_variants, f, indent=2)

    print(f"\n=== SUMMARY ===")
    print(f"Pathogenic variants: {len(pathogenic_variants)}")
    print(f"Benign variants: {len(benign_variants)}")
    print(f"Total: {len(all_variants)}")
    print(f"Output written to: {OUTPUT_PATH}")

if __name__ == "__main__":
    main()
