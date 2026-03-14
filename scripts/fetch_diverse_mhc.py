#!/usr/bin/env python3
"""Fetch diverse MHC sequences from UniProt.

Downloads MHC protein sequences from UniProt for taxonomic groups
underrepresented in IMGT/HLA and IPD-MHC (reptiles, amphibians, birds,
sharks, bony fish, marsupials, monotremes, bats).

Applies only the minimum filtering needed to exclude obvious non-MHC
proteins. All further curation is done by curate_diverse_mhc.py.

Output: data/diverse_mhc_raw.csv (intermediate; not shipped with package)

Usage:
    python scripts/fetch_diverse_mhc.py [--output PATH]

Re-run when upstream UniProt data is updated to refresh the dataset.
"""
from __future__ import annotations

import argparse
import csv
import io
import re
import time
import urllib.request
from collections import Counter
from pathlib import Path

# Taxonomic groups to query
QUERY_GROUPS = {
    # (label, taxonomy_id, extra_query_filters)
    "reptile_lepidosauria": ("8504", ""),  # lizards, snakes, tuatara
    "reptile_crocodylia": ("8493", ""),  # crocodilians
    "reptile_testudines": ("8459", ""),  # turtles
    "amphibian": ("8292", ""),  # frogs, salamanders, caecilians
    "bird_non_chicken": ("8782", "+NOT+taxonomy_id:9031"),
    "chicken": ("9031", ""),
    "shark_ray": ("7777", ""),  # Chondrichthyes
    "bony_fish": ("7898", "+NOT+taxonomy_id:8030+NOT+taxonomy_id:8022"),
    "marsupial": ("9263", ""),
    "monotreme": ("9255", ""),
    "bat": ("9397", ""),  # Chiroptera
}

# UniProt keyword-based MHC query
MHC_QUERY = (
    "(keyword:KW-0490+OR+keyword:KW-0491"
    "+OR+protein_name:%22histocompatibility%22"
    "+OR+protein_name:%22MHC+class%22)"
)

# Patterns to EXCLUDE — clearly not MHC peptide-binding chains
EXCLUDE_PATTERNS = re.compile(
    r"("
    r"minor histocompatibility|histocompatibility \(minor\)|histocompatibility 13|HM13\b|"
    r"Ig-like domain-containing|invariant chain|gamma chain|"
    r"CIITA|CD74|tapasin|TAP[12][ -]|"
    r"class IV|B-G antigen|butyrophilin|BTN\b|BTNL\b|"
    r"MHC class I polypeptide-related|MICA\b|MICB\b|"
    r"class I-related gene protein-like|"
    r"class I-like antigen recognition-like|"
    r"I-related gene protein\b|"
    r"transactivator|"
    r"regulatory factor RFX|RFX[1-5]\b|"
    r"neonatal Fc receptor|FcRn\b|"
    r"\bMR1\b|CD1[a-e]?\b|"
    r"proteasome|"
    r"ULBP|NKG2|raet1"
    r")",
    re.IGNORECASE,
)


def fetch_group(label: str, tax_id: str, extra: str) -> list[dict]:
    """Fetch one taxonomic group from UniProt REST API."""
    url = (
        f"https://rest.uniprot.org/uniprotkb/stream?"
        f"query={MHC_QUERY}+AND+taxonomy_id:{tax_id}{extra}"
        f"&format=tsv"
        f"&fields=accession,id,protein_name,gene_names,organism_name,organism_id,length,sequence,fragment"
    )
    print(f"  Fetching {label} (tax:{tax_id})...", end="", flush=True)
    try:
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "mhcseqs/0.5 (https://github.com/openvax/mhcseqs)"},
        )
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = resp.read().decode("utf-8")
        reader = csv.DictReader(io.StringIO(data), delimiter="\t")
        rows = list(reader)
        print(f" {len(rows)} entries")
        return rows
    except Exception as e:
        print(f" ERROR: {e}")
        return []


def is_mhc(protein_name: str, gene_names: str) -> bool:
    """Return True if entry looks like a real MHC peptide-binding chain or B2M."""
    combined = f"{protein_name} {gene_names}"
    if EXCLUDE_PATTERNS.search(combined):
        return False
    # Must mention histocompatibility, MHC, or a known MHC gene
    if re.search(r"(histocompat|MHC|class\s*I|beta.?2.?microglobulin|B2M\b)", combined, re.IGNORECASE):
        return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data" / "diverse_mhc_raw.csv",
        help="Output CSV path (default: data/diverse_mhc_raw.csv)",
    )
    args = parser.parse_args()

    print("Fetching diverse MHC sequences from UniProt...")
    all_rows = []
    for label, (tax_id, extra) in QUERY_GROUPS.items():
        rows = fetch_group(label, tax_id, extra)
        for r in rows:
            r["_group"] = label
        all_rows.extend(rows)
        time.sleep(0.5)  # rate limit

    print(f"\nTotal fetched: {len(all_rows)}")

    # Minimal filtering: exclude obvious non-MHC
    kept = []
    excluded = 0
    for r in all_rows:
        protein_name = r.get("Protein names", "")
        gene_names = r.get("Gene Names", "")
        if not is_mhc(protein_name, gene_names):
            excluded += 1
            continue
        kept.append({
            "uniprot_accession": r["Entry"],
            "protein_name": protein_name,
            "gene_names": gene_names,
            "organism": r.get("Organism", ""),
            "organism_id": r.get("Organism (ID)", ""),
            "length": r.get("Length", "0"),
            "is_fragment": str(r.get("Fragment", "").strip().lower() == "fragment"),
            "source_group": r["_group"],
            "sequence": r.get("Sequence", ""),
        })

    print(f"  Excluded (non-MHC): {excluded}")
    print(f"  Kept: {len(kept)}")

    by_group = Counter(r["source_group"] for r in kept)
    print("\nBy source group:")
    for g, n in by_group.most_common():
        print(f"  {g}: {n}")

    # Write raw CSV
    fields = [
        "uniprot_accession", "protein_name", "gene_names",
        "organism", "organism_id", "length", "is_fragment",
        "source_group", "sequence",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(kept)

    print(f"\nWrote {len(kept)} entries to {args.output}")
    print("Next step: python scripts/curate_diverse_mhc.py")


if __name__ == "__main__":
    main()
