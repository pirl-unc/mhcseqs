#!/usr/bin/env python3
"""Migrate the packaged diverse corpus to full-binomial gene aliases.

The historical raw input is not complete enough to regenerate every packaged
row, so this script performs the narrow, lossless name migration in place. It
also removes MHC-region genes that mhcgnomes/local curation identify as not
encoding MHC groove proteins.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from mhcseqs.alleles import is_non_mhc_gene, parse_allele_name
from mhcseqs.species import extract_latin_binomial, get_established_mhc_prefixes

try:
    from scripts.curate_diverse_mhc import (
        GENOMIC_LOCUS_RE,
        canonicalize_gene_species_alias,
        derive_species_tag,
    )
except ModuleNotFoundError:  # Running this file directly puts scripts/ on sys.path.
    from curate_diverse_mhc import GENOMIC_LOCUS_RE, canonicalize_gene_species_alias, derive_species_tag

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CSV = ROOT / "mhcseqs" / "diverse_mhc_sequences.csv"


def migrate_rows(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], Counter]:
    """Canonicalize names and return migrated rows plus counters."""
    result: list[dict[str, str]] = []
    stats: Counter = Counter()
    for row in rows:
        gene = row.get("gene", "")
        organism = row.get("organism", "")
        species = extract_latin_binomial(organism)
        if gene and GENOMIC_LOCUS_RE.match(gene):
            migrated = dict(row)
            migrated["gene"] = f"~loc:{derive_species_tag(organism)}|{gene}"
            migrated["gene_status"] = "loc"
            migrated.setdefault("raw_gene_label", gene)
            result.append(migrated)
            stats["encoded_genomic_locus"] += 1
            continue
        if gene and is_non_mhc_gene(gene, species=species):
            stats["removed_non_mhc_gene"] += 1
            continue

        canonical_gene = canonicalize_gene_species_alias(gene, organism)
        if canonical_gene != gene:
            stats["canonicalized_gene"] += 1
        migrated = dict(row)
        body = canonical_gene.split("-", 1)[-1]
        source_aliases = {prefix.casefold() for prefix in get_established_mhc_prefixes(species)}
        is_source_identifier = body.casefold() in source_aliases
        if (
            canonical_gene
            and not canonical_gene.startswith("~")
            and is_source_identifier
            and parse_allele_name(canonical_gene, species=species) is None
        ):
            # Retain paper/database identifiers as provenance rather than
            # presenting a species abbreviation as though it were a locus.
            accession = row.get("uniprot_accession", "")
            migrated["gene"] = f"~ref:{derive_species_tag(organism)}|{accession}:{body}"
            migrated["gene_status"] = "paper_specific"
            migrated["raw_gene_label"] = body
            stats["encoded_paper_specific"] += 1
        else:
            migrated["gene"] = canonical_gene
            migrated.setdefault("raw_gene_label", "")
        result.append(migrated)
    stats["kept"] = len(result)
    return result, stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--output", type=Path, default=DEFAULT_CSV)
    args = parser.parse_args()

    with open(args.input, encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    if "raw_gene_label" not in fields:
        fields.insert(fields.index("gene_status") + 1, "raw_gene_label")

    migrated, stats = migrate_rows(rows)
    output = args.output
    temporary = output.with_suffix(output.suffix + ".tmp")
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(temporary, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(migrated)
    temporary.replace(output)

    print(f"Migrated {args.input} -> {output}")
    for key, value in sorted(stats.items()):
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
