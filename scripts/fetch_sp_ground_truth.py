#!/usr/bin/env python3
"""Fetch signal peptide ground truth from UniProt.

Downloads UniProt Signal feature annotations for MHC proteins across all
vertebrate clades. These are primarily SignalP computational predictions
from unreviewed TrEMBL entries, with a few experimentally validated
annotations from reviewed Swiss-Prot entries.

Output: data/sp_ground_truth.csv

Usage:
    python scripts/fetch_sp_ground_truth.py
"""

from __future__ import annotations

import csv
import json
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "data" / "sp_ground_truth.csv"

FIELDS = [
    "accession",
    "organism",
    "taxon_id",
    "sp_length",
    "reviewed",
    "sequence",
]

CLADES = [
    ("Mammalia", 40674),
    ("Aves", 8782),
    ("Actinopterygii", 7898),
    ("Reptilia", 8504),
    ("Amphibia", 8292),
    ("Chondrichthyes", 7777),
]


def fetch_clade(name: str, taxid: int, max_results: int = 500) -> list[dict]:
    """Fetch MHC entries with Signal features for a clade."""
    query = (
        f"(taxonomy_id:{taxid}) AND "
        f"(keyword:KW-0491 OR protein_name:histocompatibility OR "
        f'protein_name:"class I" OR protein_name:"class II") AND '
        f"(ft_signal:*) AND (length:[100 TO 500])"
    )
    encoded = urllib.parse.urlencode(
        {
            "query": query,
            "format": "json",
            "size": str(max_results),
            "fields": "accession,organism_name,organism_id,ft_signal,sequence,reviewed",
        }
    )
    url = f"https://rest.uniprot.org/uniprotkb/search?{encoded}"
    req = urllib.request.Request(url)
    req.add_header("User-Agent", "mhcseqs-sp-ground-truth/1.0")

    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read())

    rows = []
    for entry in data.get("results", []):
        features = entry.get("features", [])
        sp = [f for f in features if f.get("type") == "Signal"]
        if not sp:
            continue
        sp_end = sp[0].get("location", {}).get("end", {}).get("value")
        if not sp_end:
            continue

        seq = entry.get("sequence", {}).get("value", "")
        if not seq:
            continue

        is_reviewed = entry.get("entryType", "").startswith("UniProtKB reviewed")
        organism = entry.get("organism", {}).get("scientificName", "")
        taxon = entry.get("organism", {}).get("taxonId", "")

        rows.append(
            {
                "accession": entry.get("primaryAccession", ""),
                "organism": organism,
                "taxon_id": str(taxon),
                "sp_length": str(sp_end),
                "reviewed": "Y" if is_reviewed else "N",
                "sequence": seq,
            }
        )

    return rows


def main():
    all_rows = []
    for name, taxid in CLADES:
        rows = fetch_clade(name, taxid)
        print(f"{name} (taxid {taxid}): {len(rows)} entries with SP")
        all_rows.extend(rows)

    # Deduplicate by accession
    seen = set()
    unique = []
    for row in all_rows:
        if row["accession"] not in seen:
            seen.add(row["accession"])
            unique.append(row)

    with open(OUTPUT, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for row in unique:
            writer.writerow(row)

    print(f"\nTotal: {len(unique)} entries written to {OUTPUT}")


if __name__ == "__main__":
    main()
