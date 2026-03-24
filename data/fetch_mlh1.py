#!/usr/bin/env python3
"""
Fetch 100 Pathogenic + 100 Benign MLH1 SNVs from ClinVar + UCSC hg38 flanking sequences.
"""

import json
import time
import urllib.request
import urllib.parse
import sys

OUTPUT_PATH = "/Users/dgupta/code/portfolio/EvoMine/ask_evo2/data/mlh1_variants_200.json"
FLANK = 500  # bp each side → total 1001bp

def fetch_url(url, retries=3, delay=1.0):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read().decode("utf-8")
        except Exception as e:
            print(f"  [retry {attempt+1}/{retries}] {e}", file=sys.stderr)
            time.sleep(delay * (attempt + 1))
    raise RuntimeError(f"Failed to fetch: {url}")

def esearch_ids(term, retmax=500):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = urllib.parse.urlencode({"db": "clinvar", "term": term, "retmax": retmax, "retmode": "json"})
    url = f"{base}?{params}"
    print(f"esearch: {url}")
    data = json.loads(fetch_url(url))
    ids = data["esearchresult"]["idlist"]
    total = data["esearchresult"]["count"]
    print(f"  total={total}, returned={len(ids)}")
    return ids

def esummary_batch(ids):
    """Fetch esummary for a list of IDs (up to 100 at a time)."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    id_str = ",".join(ids)
    params = urllib.parse.urlencode({"db": "clinvar", "id": id_str, "retmode": "json"})
    url = f"{base}?{params}"
    data = json.loads(fetch_url(url))
    return data.get("result", {})

def parse_variant(uid, result):
    """Parse a single esummary result. Returns dict or None if not valid."""
    rec = result.get(str(uid))
    if rec is None or not isinstance(rec, dict):
        return None

    # Classification
    gc = rec.get("germline_classification", {})
    classification = gc.get("description", "")
    if classification not in ("Pathogenic", "Benign"):
        return None

    # Variant type - must be single nucleotide variant
    variant_type = rec.get("obj_type", "")
    if "single nucleotide variant" not in variant_type.lower():
        return None

    # Get canonical SPDI
    spdi = rec.get("canonical_spdi", "")
    if not spdi:
        return None

    # Parse SPDI: accession:position:ref:alt  (0-based position)
    parts = spdi.split(":")
    if len(parts) != 4:
        return None
    _, pos_str, ref_allele, alt_allele = parts

    # Single base only
    if len(ref_allele) != 1 or ref_allele not in "ACGT":
        return None
    if len(alt_allele) != 1 or alt_allele not in "ACGT":
        return None

    # Position: SPDI is 0-based, convert to 1-based
    try:
        pos_0based = int(pos_str)
    except ValueError:
        return None
    position = pos_0based + 1  # 1-based

    # Verify chromosome from variation_loc with GRCh38
    chromosome = None
    vset = rec.get("variation_set", [])
    for vs in vset:
        vlocs = vs.get("variation_loc", [])
        for vloc in vlocs:
            if vloc.get("assembly_name", "") == "GRCh38":
                chromosome = vloc.get("chr", "")
                break
        if chromosome:
            break

    if not chromosome:
        return None

    # Must be on chr3 for MLH1
    if chromosome != "3":
        return None

    name = rec.get("title", f"variant_{uid}")

    return {
        "uid": str(uid),
        "name": name,
        "clinical_significance": classification,
        "chromosome": chromosome,
        "position": position,
        "ref_allele": ref_allele,
        "alt_allele": alt_allele,
    }

def fetch_ucsc_sequence(chrom, position_1based, flank=500):
    """
    Fetch flanking sequence from UCSC.
    position_1based: 1-based genomic position
    Returns 1001bp string (flank + ref + flank), or None on error.
    """
    # 0-based half-open: start = position - 1 - flank, end = position + flank
    start = position_1based - 1 - flank   # 0-based start (inclusive)
    end = position_1based + flank          # 0-based end (exclusive) → 1001bp total

    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38&chrom=chr{chrom}&start={start}&end={end}"
    try:
        raw = fetch_url(url)
        data = json.loads(raw)
        seq = data.get("dna", "")
        return seq.upper() if seq else None
    except Exception as e:
        print(f"  UCSC fetch error for chr{chrom}:{position_1based}: {e}", file=sys.stderr)
        return None

def collect_variants(search_term, target_count, label):
    """Collect target_count valid variants matching label (Pathogenic/Benign)."""
    print(f"\n=== Collecting {target_count} {label} variants ===")

    # Fetch enough candidate IDs
    all_ids = esearch_ids(search_term, retmax=500)
    print(f"  Candidate pool: {len(all_ids)} IDs")

    valid = []
    batch_size = 100
    i = 0

    while len(valid) < target_count and i < len(all_ids):
        batch = all_ids[i:i+batch_size]
        i += batch_size
        print(f"  Processing IDs {i-len(batch)+1}–{i} (have {len(valid)}/{target_count} valid so far)...")

        time.sleep(0.34)  # NCBI rate limit: ~3 req/sec
        result = esummary_batch(batch)

        for uid in batch:
            parsed = parse_variant(uid, result)
            if parsed is None:
                continue

            # Fetch UCSC sequence
            print(f"    Fetching UCSC for {parsed['name']} chr{parsed['chromosome']}:{parsed['position']}")
            time.sleep(0.2)
            seq = fetch_ucsc_sequence(parsed["chromosome"], parsed["position"])
            if seq is None or len(seq) != 1001:
                print(f"    SKIP: bad sequence length {len(seq) if seq else 'None'}")
                continue

            # Verify ref base at index 500
            ref_at_500 = seq[500]
            if ref_at_500 != parsed["ref_allele"]:
                print(f"    SKIP: ref mismatch at index 500: got {ref_at_500}, expected {parsed['ref_allele']}")
                continue

            # Build alt_seq
            alt_seq = seq[:500] + parsed["alt_allele"] + seq[501:]

            parsed["ref_seq"] = seq
            parsed["alt_seq"] = alt_seq
            del parsed["uid"]

            valid.append(parsed)
            print(f"    VALID [{len(valid)}/{target_count}]: {parsed['name']}")

            if len(valid) >= target_count:
                break

    print(f"  Collected {len(valid)} {label} variants")
    return valid

def main():
    pathogenic_term = "MLH1[gene] AND Pathogenic[ClinSig]"
    benign_term = "MLH1[gene] AND Benign[ClinSig]"

    pathogenic = collect_variants(pathogenic_term, 100, "Pathogenic")
    benign = collect_variants(benign_term, 100, "Benign")

    if len(pathogenic) < 100:
        print(f"WARNING: Only collected {len(pathogenic)} Pathogenic variants (target: 100)")
    if len(benign) < 100:
        print(f"WARNING: Only collected {len(benign)} Benign variants (target: 100)")

    combined = pathogenic + benign
    print(f"\nTotal variants: {len(combined)} ({len(pathogenic)} Pathogenic + {len(benign)} Benign)")

    with open(OUTPUT_PATH, "w") as f:
        json.dump(combined, f, indent=2)

    print(f"\nWrote {len(combined)} records to {OUTPUT_PATH}")
    print("\n=== SUMMARY ===")
    print(f"  Pathogenic: {len(pathogenic)}")
    print(f"  Benign:     {len(benign)}")
    print(f"  Total:      {len(combined)}")

    # Show first few of each
    for label in ("Pathogenic", "Benign"):
        subset = [v for v in combined if v["clinical_significance"] == label][:3]
        print(f"\n  First 3 {label}:")
        for v in subset:
            print(f"    {v['name']}  chr{v['chromosome']}:{v['position']}  {v['ref_allele']}>{v['alt_allele']}")

if __name__ == "__main__":
    main()
