#!/usr/bin/env python3
"""
Fetch 100 Pathogenic and 100 Benign BRCA1 SNVs from ClinVar + UCSC flanking sequences.
Uses pagination to get more candidates.
"""

import json
import time
import requests

OUTPUT_PATH = "/Users/dgupta/code/portfolio/EvoMine/ask_evo2/data/brca1_variants_200.json"
CHROM = "chr17"
FLANK = 500

HEADERS = {"User-Agent": "BRCAResearch/1.0"}

def clinvar_esearch_paged(term, retmax=500, retstart=0):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "clinvar",
        "term": term,
        "retmax": retmax,
        "retstart": retstart,
        "retmode": "json",
    }
    r = requests.get(url, params=params, timeout=30, headers=HEADERS)
    r.raise_for_status()
    data = r.json()
    ids = data["esearchresult"]["idlist"]
    count = int(data["esearchresult"]["count"])
    return ids, count

def get_ids_for_term(term, max_ids=2000):
    """Get up to max_ids IDs using pagination."""
    all_ids = []
    seen = set()
    retstart = 0
    retmax = 500

    while len(all_ids) < max_ids:
        ids, total_count = clinvar_esearch_paged(term, retmax=retmax, retstart=retstart)
        if not ids:
            break
        for id_ in ids:
            if id_ not in seen:
                seen.add(id_)
                all_ids.append(id_)
        retstart += len(ids)
        if retstart >= total_count or retstart >= max_ids:
            break
        time.sleep(0.4)

    return all_ids, total_count

def clinvar_esummary_batch(ids):
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
    url = "https://api.genome.ucsc.edu/getData/sequence"
    params = {"genome": "hg38", "chrom": chrom, "start": start, "end": end}
    r = requests.get(url, params=params, timeout=30, headers=HEADERS)
    r.raise_for_status()
    return r.json().get("dna", "").upper()

def parse_spdi(spdi_str):
    if not spdi_str:
        return None
    parts = spdi_str.split(":")
    if len(parts) != 4:
        return None
    try:
        pos_0based = int(parts[1])
    except ValueError:
        return None
    return parts[0], pos_0based, parts[2].upper(), parts[3].upper()

def parse_variant(uid, summary):
    variation_set = summary.get("variation_set", [])
    if not variation_set:
        return None
    vs = variation_set[0]

    vtype = vs.get("variant_type", "").lower()
    if "single nucleotide" not in vtype:
        return None

    spdi_str = vs.get("canonical_spdi", "")
    parsed = parse_spdi(spdi_str)
    if parsed is None:
        return None

    seq_acc, pos_0based, ref_allele, alt_allele = parsed

    if not seq_acc.startswith("NC_000017"):
        return None
    if len(ref_allele) != 1 or len(alt_allele) != 1:
        return None
    if ref_allele not in "ACGT" or alt_allele not in "ACGT":
        return None
    if ref_allele == alt_allele:
        return None

    return {
        "uid": uid,
        "name": summary.get("title", f"VCV{uid}"),
        "chromosome": CHROM,
        "position": pos_0based + 1,
        "pos_0based": pos_0based,
        "ref_allele": ref_allele,
        "alt_allele": alt_allele,
    }

def fetch_and_verify_sequence(variant):
    pos_0based = variant["pos_0based"]
    start = pos_0based - FLANK
    end = pos_0based + FLANK + 1

    if start < 0:
        return None, None

    try:
        seq = ucsc_sequence(CHROM, start, end)
    except Exception as e:
        print(f"  UCSC error uid {variant['uid']}: {e}")
        return None, None

    if len(seq) != 1001:
        print(f"  Wrong seq length {len(seq)} uid {variant['uid']}")
        return None, None

    if seq[500] != variant["ref_allele"]:
        print(f"  Ref mismatch uid {variant['uid']}: expected {variant['ref_allele']} got {seq[500]} @ pos {variant['position']}")
        return None, None

    alt_seq = seq[:500] + variant["alt_allele"] + seq[501:]
    return seq, alt_seq

def collect_variants(search_terms, exact_clinsig, needed=100):
    print(f"\n{'='*60}")
    print(f"Collecting {needed} {exact_clinsig} BRCA1 SNVs")
    print(f"{'='*60}")

    # Collect unique IDs from all search terms with pagination
    all_ids = []
    seen_ids = set()
    for term in search_terms:
        ids, total = get_ids_for_term(term, max_ids=3000)
        print(f"  Term: '{term[:70]}' → total={total}, got {len(ids)} IDs")
        for id_ in ids:
            if id_ not in seen_ids:
                seen_ids.add(id_)
                all_ids.append(id_)
        time.sleep(0.5)

    print(f"Total unique candidates: {len(all_ids)}")

    valid = []
    dedup_pos = set()
    batch_size = 100

    for batch_start in range(0, len(all_ids), batch_size):
        if len(valid) >= needed:
            break

        batch = all_ids[batch_start:batch_start + batch_size]
        print(f"\nBatch {batch_start//batch_size + 1}: IDs {batch_start}–{batch_start+len(batch)-1}, valid={len(valid)}/{needed}")

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
            if len(valid) >= needed:
                break
            if uid == "uids":
                continue

            summary = result.get(str(uid), result.get(uid, {}))
            if not summary:
                continue

            # Strict classification filter
            clinsig = summary.get("germline_classification", {}).get("description", "")
            if clinsig != exact_clinsig:
                continue

            variant = parse_variant(uid, summary)
            if variant is None:
                continue

            pos_key = (variant["pos_0based"], variant["ref_allele"], variant["alt_allele"])
            if pos_key in dedup_pos:
                continue
            dedup_pos.add(pos_key)

            print(f"  [{len(valid)+1}/{needed}] uid={uid} chr17:{variant['position']} {variant['ref_allele']}>{variant['alt_allele']}")

            ref_seq, alt_seq = fetch_and_verify_sequence(variant)
            time.sleep(0.25)

            if ref_seq is None:
                continue

            valid.append({
                "name": variant["name"],
                "clinical_significance": exact_clinsig,
                "chromosome": CHROM,
                "position": variant["position"],
                "ref_allele": variant["ref_allele"],
                "alt_allele": variant["alt_allele"],
                "ref_seq": ref_seq,
                "alt_seq": alt_seq,
            })

    print(f"Found {len(valid)} valid {exact_clinsig} variants")
    return valid

def main():
    path_terms = [
        'BRCA1[gene] AND "pathogenic"[clinsig] AND "reviewed by expert panel"[review_status]',
        'BRCA1[gene] AND "pathogenic"[clinsig] AND "criteria provided, single submitter"[review_status]',
        'BRCA1[gene] AND "pathogenic"[clinsig] AND "criteria provided, multiple submitters, no conflicts"[review_status]',
        'BRCA1[gene] AND "pathogenic"[clinsig]',
    ]

    benign_terms = [
        'BRCA1[gene] AND "benign"[clinsig] AND "reviewed by expert panel"[review_status]',
        'BRCA1[gene] AND "benign"[clinsig] AND "criteria provided, single submitter"[review_status]',
        'BRCA1[gene] AND "benign"[clinsig] AND "criteria provided, multiple submitters, no conflicts"[review_status]',
        'BRCA1[gene] AND "benign"[clinsig]',
    ]

    pathogenic_variants = collect_variants(path_terms, "Pathogenic", needed=100)
    benign_variants = collect_variants(benign_terms, "Benign", needed=100)

    all_variants = pathogenic_variants + benign_variants

    with open(OUTPUT_PATH, "w") as f:
        json.dump(all_variants, f, indent=2)

    print(f"\n{'='*60}")
    print("FINAL SUMMARY")
    print(f"{'='*60}")
    print(f"Pathogenic: {len(pathogenic_variants)}")
    print(f"Benign:     {len(benign_variants)}")
    print(f"Total:      {len(all_variants)}")
    print(f"Output:     {OUTPUT_PATH}")
    if len(pathogenic_variants) < 100:
        print(f"WARNING: Only {len(pathogenic_variants)} Pathogenic (needed 100)")
    if len(benign_variants) < 100:
        print(f"WARNING: Only {len(benign_variants)} Benign (needed 100)")

if __name__ == "__main__":
    main()
