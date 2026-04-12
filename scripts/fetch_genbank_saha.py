#!/usr/bin/env python3
"""Fetch Sarcophilus harrisii (Tasmanian devil) SahaI MHC class I allele
sequences from NCBI GenBank.

These alleles (SahaI*1 ... SahaI*98 plus pa/pb variants) were deposited by
Siddle, Cheng, and colleagues as part of multiple Tasmanian devil MHC
diversity studies (2010-2015). They include all classical class I alleles
discussed in Caldwell et al. 2018 (eLife, DFT2) -- e.g. SahaI*27, *32, *35,
*74, *88, *90 -- which are not available through IMGT/HLA or IPD-MHC and
not yet assigned a locus (Saha-UA/UB/UC) in UniProt.

Output: mhcseqs/genbank_saha_sequences.csv (shipped with the package;
re-run only when NCBI adds new SahaI entries).

Usage:
    python scripts/fetch_genbank_saha.py
"""

from __future__ import annotations

import argparse
import csv
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SEARCH_TERM = "Sarcophilus harrisii[Organism] AND SahaI[Title]"
BATCH_SIZE = 50
USER_AGENT = "mhcseqs/0.5 (https://github.com/openvax/mhcseqs)"


def _urlopen(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return resp.read().decode("utf-8")


def esearch_ids() -> list[str]:
    url = f"{EUTILS}/esearch.fcgi?db=nuccore&term={urllib.parse.quote(SEARCH_TERM)}&retmax=500"
    xml = _urlopen(url)
    return re.findall(r"<Id>(\d+)</Id>", xml)


def efetch_gb(uids: list[str]) -> str:
    url = f"{EUTILS}/efetch.fcgi?db=nuccore&id={','.join(uids)}&rettype=gb&retmode=text"
    return _urlopen(url)


def parse_gb_records(text: str) -> list[dict]:
    records = []
    for rec in text.split("\n//\n"):
        rec = rec.strip()
        if not rec or not rec.startswith("LOCUS"):
            continue
        acc = re.search(r"VERSION\s+(\S+)", rec)
        gene = re.search(r'/gene="([^"]+)"', rec)
        allele = re.search(r'/allele="([^"]+)"', rec)
        note = re.search(r'/note="([^"]*)"', rec)
        pseudo = "/pseudogene=" in rec or "/pseudo" in rec
        # /translation can span multiple lines with embedded whitespace
        trans_m = re.search(r'/translation="([^"]+)"', rec, re.DOTALL)
        translation = re.sub(r"\s+", "", trans_m.group(1)) if trans_m else ""
        if not (acc and allele):
            continue
        records.append(
            {
                "accession": acc.group(1),
                "allele_name": allele.group(1),
                "gene": gene.group(1) if gene else "SahaI",
                "sequence": translation,
                "length": len(translation),
                "is_pseudogene": pseudo,
                "note": note.group(1).strip().replace("\n", " ") if note else "",
            }
        )
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "mhcseqs" / "genbank_saha_sequences.csv",
    )
    args = parser.parse_args()

    print(f"Searching NCBI: {SEARCH_TERM}")
    uids = esearch_ids()
    print(f"  Found {len(uids)} UIDs")

    all_records: list[dict] = []
    for i in range(0, len(uids), BATCH_SIZE):
        batch = uids[i : i + BATCH_SIZE]
        print(f"  Fetching batch {i // BATCH_SIZE + 1} ({len(batch)} records)...")
        gb_text = efetch_gb(batch)
        all_records.extend(parse_gb_records(gb_text))
        time.sleep(0.4)

    with_translation = [r for r in all_records if r["sequence"]]
    without_translation = [r for r in all_records if not r["sequence"]]
    print(f"Parsed {len(all_records)} records ({len(with_translation)} with translation, {len(without_translation)} without)")
    if without_translation:
        print("  No-translation (skipped):")
        for r in without_translation:
            print(f"    {r['accession']} {r['allele_name']}")

    with_translation.sort(key=lambda r: r["accession"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        fields = ["accession", "allele_name", "gene", "length", "is_pseudogene", "note", "sequence"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in with_translation:
            w.writerow(
                {
                    "accession": r["accession"],
                    "allele_name": r["allele_name"],
                    "gene": r["gene"],
                    "length": r["length"],
                    "is_pseudogene": str(r["is_pseudogene"]),
                    "note": r["note"],
                    "sequence": r["sequence"],
                }
            )
    print(f"Wrote {len(with_translation)} records to {args.output}")


if __name__ == "__main__":
    main()
