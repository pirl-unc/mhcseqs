#!/usr/bin/env python3
"""Validate signal peptide inference against UniProt curated annotations.

Fetches UniProt entries for well-characterized MHC proteins and compares
their curated signal peptide lengths against our Cys-pair-based inference.
"""

import csv
import json
import sys
import urllib.request
from pathlib import Path

# Representative UniProt accessions for major HLA genes
# These are the canonical/reviewed Swiss-Prot entries
UNIPROT_MHC_ENTRIES = {
    # Class I
    "P01892": "HLA-A*02:01",  # HLA-A
    "P30461": "HLA-A*03:01",  # HLA-A
    "P30685": "HLA-B*07:02",  # HLA-B
    "P10321": "HLA-B*08:01",  # HLA-B
    "P30504": "HLA-C*04:01",  # HLA-C
    "P13747": "HLA-E*01:01",  # HLA-E
    "P30511": "HLA-F*01:01",  # HLA-F
    "P17693": "HLA-G*01:01",  # HLA-G
    # Class II alpha
    "P01903": "HLA-DRA*01:01",  # HLA-DRA
    "P01920": "HLA-DQA1*01:01",  # HLA-DQA1
    "P20036": "HLA-DPA1*01:03",  # HLA-DPA1
    "P28067": "HLA-DMA*01:01",  # HLA-DMA
    "P06340": "HLA-DOA*01:01",  # HLA-DOA
    # Class II beta
    "P01911": "HLA-DRB1*01:01",  # HLA-DRB1
    "P04229": "HLA-DRB1*04:01",  # HLA-DRB1
    "P01918": "HLA-DQB1*02:01",  # HLA-DQB1
    "P04440": "HLA-DPB1*04:01",  # HLA-DPB1
    "P28068": "HLA-DMB*01:01",  # HLA-DMB
    "P06341": "HLA-DOB*01:01",  # HLA-DOB
}


def fetch_uniprot_sp(accession: str) -> dict:
    """Fetch signal peptide and sequence info from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read())

    result = {"accession": accession, "sp_start": None, "sp_end": None, "seq_len": None}

    features = data.get("features", [])
    for f in features:
        if f.get("type") == "Signal":
            loc = f.get("location", {})
            result["sp_start"] = loc.get("start", {}).get("value")
            result["sp_end"] = loc.get("end", {}).get("value")
            break

    seq_info = data.get("sequence", {})
    result["seq_len"] = seq_info.get("length")
    result["sequence"] = seq_info.get("value", "")

    genes = data.get("genes", [])
    result["gene"] = genes[0].get("geneName", {}).get("value", "") if genes else ""

    return result


def main():
    raw_csv = Path("mhc-seqs-raw.csv")
    if not raw_csv.exists():
        print("Run 'python build.py' first to generate mhc-seqs-raw.csv", file=sys.stderr)
        sys.exit(1)

    # Load our inferred signal peptide data
    our_data = {}
    with open(raw_csv, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            our_data[row["allele_normalized"]] = row

    print("Signal Peptide Validation: Cys-pair inference vs UniProt curated annotations")
    print("=" * 90)
    cols = f"{'Allele':<25} {'UniProt':>8} {'Ours':>6} {'Delta':>7}"
    cols += f" {'UP SP len':>15} {'Our SP len':>11} {'Match':>6}"
    print(cols)
    print("-" * 90)

    matches = 0
    mismatches = 0
    skipped = 0

    seen = set()
    for accession, allele_hint in sorted(UNIPROT_MHC_ENTRIES.items(), key=lambda x: x[1]):
        if accession in seen:
            continue
        seen.add(accession)

        try:
            up = fetch_uniprot_sp(accession)
        except Exception as e:
            print(f"  {allele_hint:<25} FETCH ERROR: {e}")
            skipped += 1
            continue

        up_sp_len = up["sp_end"] if up["sp_end"] else 0
        up_seq_len = up["seq_len"] or 0

        # Find matching entry in our raw CSV
        # Try exact allele match first, then prefix match
        our_row = our_data.get(allele_hint)
        if not our_row:
            # Try finding by sequence match
            for key, row in our_data.items():
                if up["sequence"] and row["sequence"] == up["sequence"]:
                    our_row = row
                    break
        if not our_row:
            # Try by allele prefix (two-field)
            for key, row in our_data.items():
                if key.startswith(allele_hint.replace("*", "*")):
                    our_row = row
                    break

        if not our_row:
            print(f"  {allele_hint:<25} {up_sp_len:>8} {'N/A':>6} {'':>7} {f'1-{up_sp_len}':>15} {'not found':>11}")
            skipped += 1
            continue

        our_sp_len = int(our_row.get("signal_peptide_len", "0") or "0")
        our_seq_len = int(our_row.get("seq_len", "0") or "0")
        delta = our_sp_len - up_sp_len

        if delta == 0:
            matches += 1
            match_str = "OK"
        else:
            mismatches += 1
            match_str = "DIFF"

        up_range = f"1-{up_sp_len}"
        row_str = f"  {allele_hint:<25} {up_seq_len:>8} {our_seq_len:>6} {delta:>+7}"
        row_str += f" {up_range:>15} {our_sp_len:>11} {match_str:>6}"
        print(row_str)

    print("-" * 90)
    print(f"Summary: {matches} match, {mismatches} differ, {skipped} skipped")
    if mismatches > 0:
        print("Note: small differences (±1-3 aa) may reflect different allele variants or")
        print("      the Cys-pair heuristic inferring a slightly different boundary.")


if __name__ == "__main__":
    main()
