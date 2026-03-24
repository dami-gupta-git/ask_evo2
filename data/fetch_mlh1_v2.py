#!/usr/bin/env python3
"""
Fetch 100 Pathogenic + 100 Benign MLH1 SNVs from ClinVar + UCSC hg38 flanking sequences.

Key fixes from debug:
- canonical_spdi is in variation_set[0]["canonical_spdi"]
- variant_type is in variation_set[0]["variant_type"]
- Need to scan many more IDs since ~29% are SNVs and many won't be exact Pathogenic/Benign
"""

import json
import time
import urllib.request
import urllib.parse
import sys

OUTPUT_PATH = "/Users/dgupta/code/portfolio/EvoMine/ask_evo2/data/mlh1_variants_200.json"
FLANK = 500  # bp each side → total 1001bp

def fetch_url(url, retries=4, delay=1.5):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read().decode("utf-8")
        except Exception as e:
            wait = delay * (attempt + 1)
            print(f"  [retry {attempt+1}/{retries}] {e} (waiting {wait:.1f}s)", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"Failed to fetch: {url}")

def esearch_ids(term, retmax=500, retstart=0):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = urllib.parse.urlencode({
        "db": "clinvar", "term": term,
        "retmax": retmax, "retstart": retstart, "retmode": "json"
    })
    url = f"{base}?{params}"
    print(f"esearch (start={retstart}): {term[:60]}...")
    data = json.loads(fetch_url(url))
    ids = data["esearchresult"]["idlist"]
    total = int(data["esearchresult"]["count"])
    print(f"  total available={total}, returned={len(ids)}")
    return ids, total

def esummary_batch(ids):
    """Fetch esummary for a list of IDs (up to 100 at a time)."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    id_str = ",".join(str(i) for i in ids)
    params = urllib.parse.urlencode({"db": "clinvar", "id": id_str, "retmode": "json"})
    url = f"{base}?{params}"
    data = json.loads(fetch_url(url))
    return data.get("result", {})

def parse_variant(uid, result, required_classification):
    """
    Parse a single esummary result.
    Returns dict or None if not valid SNV with exact classification.
    """
    rec = result.get(str(uid))
    if rec is None or not isinstance(rec, dict):
        return None

    # Classification must be exactly "Pathogenic" or "Benign"
    gc = rec.get("germline_classification", {})
    classification = gc.get("description", "")
    if classification != required_classification:
        return None

    # Get variation_set[0] for type and SPDI
    vset = rec.get("variation_set", [])
    if not vset:
        return None
    vs0 = vset[0]

    # Variant type must be single nucleotide variant
    variant_type = vs0.get("variant_type", "")
    if "single nucleotide variant" not in variant_type.lower():
        return None

    # Get canonical SPDI from variation_set[0]
    spdi = vs0.get("canonical_spdi", "")
    if not spdi:
        return None

    # Parse SPDI: accession:position:ref:alt (0-based position)
    parts = spdi.split(":")
    if len(parts) != 4:
        return None
    _, pos_str, ref_allele, alt_allele = parts

    # Single base only
    if len(ref_allele) != 1 or ref_allele not in "ACGT":
        return None
    if len(alt_allele) != 1 or alt_allele not in "ACGT":
        return None
    if ref_allele == alt_allele:
        return None

    # SPDI position is 0-based → convert to 1-based
    try:
        pos_0based = int(pos_str)
    except ValueError:
        return None
    position = pos_0based + 1  # 1-based

    # Verify chromosome from variation_loc with GRCh38
    chromosome = None
    vlocs = vs0.get("variation_loc", [])
    for vloc in vlocs:
        if vloc.get("assembly_name", "") == "GRCh38":
            chromosome = vloc.get("chr", "")
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
    Returns 1001bp string, or None on error.
    UCSC API uses 0-based half-open coordinates.
    """
    start = position_1based - 1 - flank   # 0-based inclusive
    end = position_1based + flank          # 0-based exclusive → 1001bp total

    if start < 0:
        return None

    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38&chrom=chr{chrom}&start={start}&end={end}"
    try:
        raw = fetch_url(url)
        data = json.loads(raw)
        seq = data.get("dna", "")
        return seq.upper() if seq else None
    except Exception as e:
        print(f"  UCSC error for chr{chrom}:{position_1based}: {e}", file=sys.stderr)
        return None

def collect_variants(search_term, target_count, required_classification):
    """Collect target_count valid variants."""
    print(f"\n=== Collecting {target_count} {required_classification} variants ===")

    valid = []
    seen_positions = set()  # avoid duplicate positions
    offset = 0
    batch_size = 100  # esummary batch
    search_batch = 500  # how many IDs to fetch per esearch call

    while len(valid) < target_count:
        # Fetch IDs
        all_ids, total_available = esearch_ids(search_term, retmax=search_batch, retstart=offset)
        if not all_ids:
            print(f"  No more IDs available (fetched {offset} of {total_available})")
            break

        offset += len(all_ids)

        # Process in batches of 100
        i = 0
        while len(valid) < target_count and i < len(all_ids):
            batch = all_ids[i:i+batch_size]
            i += batch_size
            abs_pos = offset - len(all_ids) + i
            print(f"  Esummary batch (pos ~{abs_pos}/{total_available}), have {len(valid)}/{target_count}...")

            time.sleep(0.4)  # NCBI rate limit
            result = esummary_batch(batch)

            for uid in batch:
                parsed = parse_variant(uid, result, required_classification)
                if parsed is None:
                    continue

                # Deduplicate by position
                pos_key = (parsed["chromosome"], parsed["position"])
                if pos_key in seen_positions:
                    continue
                seen_positions.add(pos_key)

                # Fetch UCSC sequence
                print(f"    Fetching UCSC for chr{parsed['chromosome']}:{parsed['position']} {parsed['ref_allele']}>{parsed['alt_allele']}")
                time.sleep(0.2)
                seq = fetch_ucsc_sequence(parsed["chromosome"], parsed["position"])

                if seq is None or len(seq) != 1001:
                    print(f"    SKIP: bad sequence (got len={len(seq) if seq else 'None'})")
                    continue

                # Verify ref base at index 500
                ref_at_500 = seq[500]
                if ref_at_500 != parsed["ref_allele"]:
                    print(f"    SKIP: ref mismatch at index 500: got {ref_at_500}, expected {parsed['ref_allele']}")
                    continue

                # Build alt_seq (substitute at index 500)
                alt_seq = seq[:500] + parsed["alt_allele"] + seq[501:]

                parsed["ref_seq"] = seq
                parsed["alt_seq"] = alt_seq
                del parsed["uid"]

                valid.append(parsed)
                print(f"    VALID [{len(valid)}/{target_count}]: {parsed['name']}")

                if len(valid) >= target_count:
                    break

        if offset >= total_available:
            print(f"  Exhausted all {total_available} available IDs")
            break

        time.sleep(0.5)

    print(f"  Collected {len(valid)} {required_classification} variants")
    return valid

def main():
    pathogenic_term = "MLH1[gene] AND Pathogenic[ClinSig]"
    benign_term = "MLH1[gene] AND Benign[ClinSig]"

    pathogenic = collect_variants(pathogenic_term, 100, "Pathogenic")
    benign = collect_variants(benign_term, 100, "Benign")

    if len(pathogenic) < 100:
        print(f"\nWARNING: Only collected {len(pathogenic)} Pathogenic variants (target: 100)")
    if len(benign) < 100:
        print(f"\nWARNING: Only collected {len(benign)} Benign variants (target: 100)")

    combined = pathogenic + benign
    print(f"\nTotal variants: {len(combined)} ({len(pathogenic)} Pathogenic + {len(benign)} Benign)")

    with open(OUTPUT_PATH, "w") as f:
        json.dump(combined, f, indent=2)

    print(f"\nWrote {len(combined)} records to {OUTPUT_PATH}")
    print("\n=== SUMMARY ===")
    print(f"  Pathogenic: {len(pathogenic)}")
    print(f"  Benign:     {len(benign)}")
    print(f"  Total:      {len(combined)}")

    for label in ("Pathogenic", "Benign"):
        subset = [v for v in combined if v["clinical_significance"] == label][:3]
        print(f"\n  First 3 {label}:")
        for v in subset:
            print(f"    {v['name']}")
            print(f"      chr{v['chromosome']}:{v['position']}  {v['ref_allele']}>{v['alt_allele']}")
            print(f"      ref_seq[498:503] = {v['ref_seq'][498:503]}")

if __name__ == "__main__":
    main()
