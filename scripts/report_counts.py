#!/usr/bin/env python3
"""Report build counts and mhcgnomes parse rate.

Reads the built CSV from ~/.cache/mhcseqs/ and prints a summary.
Intended to run after `mhcseqs build` or in CI after merge to main.

Usage:
    python scripts/report_counts.py
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

BUILT_CSV = Path.home() / ".cache" / "mhcseqs" / "mhc-full-seqs.csv"
DIVERSE_CSV = Path(__file__).resolve().parent.parent / "mhcseqs" / "diverse_mhc_sequences.csv"
B2M_CSV = Path(__file__).resolve().parent.parent / "mhcseqs" / "b2m_sequences.csv"


def main():
    try:
        import mhcgnomes

        mhcgnomes_version = mhcgnomes.__version__
    except ImportError:
        mhcgnomes_version = "not installed"

    from mhcseqs.species import extract_latin_binomial
    from mhcseqs.version import __version__

    print(f"mhcseqs {__version__} | mhcgnomes {mhcgnomes_version}")

    if not BUILT_CSV.exists():
        print(f"\nNo built CSV at {BUILT_CSV} — run `mhcseqs build` first.")
        sys.exit(1)

    # Read full-seqs CSV
    total = 0
    with_sp = 0
    groove = Counter()
    species_cat = Counter()
    class_chain = Counter()

    with open(BUILT_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            total += 1
            groove[row.get("groove_status", "")] += 1
            species_cat[row.get("species_category", "")] += 1
            cc = f"{row.get('mhc_class', '')} {row.get('chain', '')}".strip()
            class_chain[cc] += 1
            ms = int(row.get("mature_start", "0") or "0")
            if ms >= 15:
                with_sp += 1

    # B2M count
    b2m_count = 0
    if B2M_CSV.exists():
        with open(B2M_CSV, "r", encoding="utf-8") as f:
            b2m_count = sum(1 for _ in csv.DictReader(f))

    # Diverse parse rate (deduplicate by gene+organism, not gene alone)
    diverse_total = 0
    seen_pairs = set()
    parsed_ok = 0
    if DIVERSE_CSV.exists() and mhcgnomes_version != "not installed":
        import inspect

        parse_params = inspect.signature(mhcgnomes.parse).parameters
        sp_kwarg = "species" if "species" in parse_params else "default_species"

        with open(DIVERSE_CSV, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                diverse_total += 1
                gene = row.get("gene", "")
                organism = row.get("organism", "")
                if not gene or (gene, organism) in seen_pairs:
                    continue
                seen_pairs.add((gene, organism))

                latin = extract_latin_binomial(organism)

                # Ortholog-transferred genes: parse with source species
                # Format: {actual}-ortho:{source_prefix}:{ortholog_gene}
                if "-ortho:" in gene:
                    rest = gene.split("-ortho:", 1)[1]
                    if ":" in rest:
                        source_prefix, bare = rest.split(":", 1)
                        try:
                            sp_obj = mhcgnomes.Species.get(source_prefix)
                            parse_species = sp_obj.name if sp_obj else latin
                        except Exception:
                            parse_species = latin
                    else:
                        bare = rest
                        parse_species = latin
                else:
                    bare = gene.split("-", 1)[1] if "-" in gene else gene
                    prefix = gene.split("-")[0] if "-" in gene else ""
                    if bare.lower().startswith(prefix.lower()) and prefix:
                        bare = bare[len(prefix) :].lstrip("-_")
                    parse_species = latin

                if parse_species:
                    try:
                        r = mhcgnomes.parse(bare, **{sp_kwarg: parse_species})
                        tp = type(r).__name__
                        if tp in ("Gene", "Allele", "AlleleWithoutGene"):
                            if "-ortho:" in gene:
                                parsed_ok += 1
                            else:
                                sp = getattr(getattr(r, "species", None), "name", "")
                                if sp and sp.lower() in organism.lower():
                                    parsed_ok += 1
                    except Exception:
                        pass

    # Print report
    print(f"\n{'=' * 50}")
    print("BUILD COUNTS")
    print(f"{'=' * 50}")
    print(f"Total alleles:          {total:,}")
    print(f"With signal peptide:    {with_sp:,}")
    print(f"B2M references:         {b2m_count}")
    print(f"Diverse MHC entries:    {diverse_total:,}")
    if seen_pairs:
        print(f"Diverse parse rate:     {parsed_ok}/{len(seen_pairs)} ({parsed_ok / len(seen_pairs) * 100:.1f}%)")

    cats = [k for k, _ in species_cat.most_common() if k]
    print(f"\nSpecies categories ({len(cats)}):")
    for cat in cats:
        print(f"  {cat:20s} {species_cat[cat]:6,}")

    print("\nGroove status:")
    for status, n in groove.most_common():
        if not status:
            status = "(empty)"
        print(f"  {status:25s} {n:6,} ({n / total * 100:.1f}%)")

    print("\nClass/chain:")
    for cc, n in class_chain.most_common():
        if not cc:
            cc = "(empty)"
        print(f"  {cc:20s} {n:6,}")


if __name__ == "__main__":
    main()
