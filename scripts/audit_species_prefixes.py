#!/usr/bin/env python3
"""Audit whether species prefixes were generated or source-attested.

This is primarily a migration tool for the pre-full-binomial diverse MHC
corpus. It compares the curated gene prefix with the raw UniProt gene/protein
labels and the explicitly sourced short-prefix allowlist in ``mhcseqs``.

Usage:
    python scripts/audit_species_prefixes.py --output data/species_prefix_audit.csv
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from mhcseqs.species import (
    ESTABLISHED_MHC_PREFIX_SOURCES,
    extract_latin_binomial,
    full_species_name_alias,
    get_established_mhc_prefix,
)

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CURATED = ROOT / "mhcseqs" / "diverse_mhc_sequences.csv"
DEFAULT_RAW = ROOT / "data" / "diverse_mhc_raw.csv"

AUDIT_FIELDS = [
    "species",
    "legacy_prefix",
    "status",
    "row_count",
    "source_attested_row_count",
    "evidence",
    "example_accessions",
    "example_source_labels",
]


def legacy_2plus2_prefix(species: str) -> str:
    """Reproduce the retired mhcseqs 2+2 generator for audit purposes."""
    parts = species.split()
    if len(parts) >= 2:
        return parts[0][:2].capitalize() + parts[1][:2].lower()
    return parts[0][:4].capitalize() if parts else ""


def prefix_is_present(prefix: str, source_text: str) -> bool:
    """Return whether a source label explicitly contains ``prefix``."""
    if len(prefix) < 2:
        return False
    pattern = rf"(?<![A-Za-z]){re.escape(prefix)}(?=[-_A-Z0-9]|\b)"
    return re.search(pattern, source_text) is not None


def classify_prefix(species: str, prefix: str, source_accessions: list[str]) -> tuple[str, str]:
    """Classify one species/prefix pair and return status plus evidence."""
    if prefix.casefold() == full_species_name_alias(species).casefold():
        return "canonical_full_binomial", "self-describing full scientific-name alias"

    established = get_established_mhc_prefix(species)
    if established and prefix.casefold() == established.casefold():
        source = ESTABLISHED_MHC_PREFIX_SOURCES[established]
        return "external_registry_or_literature", source

    if source_accessions:
        urls = [f"https://rest.uniprot.org/uniprotkb/{accession}" for accession in source_accessions[:3]]
        return "external_uniprot_label", ";".join(urls)

    if prefix.casefold() == legacy_2plus2_prefix(species).casefold():
        return (
            "generated_by_mhcseqs_2plus2",
            "retired 2+2 rule; absent from available raw UniProt labels",
        )

    return (
        "unverified_historical",
        "not in the sourced allowlist or available raw UniProt labels",
    )


def audit_prefixes(curated_path: Path, raw_path: Path) -> list[dict[str, str | int]]:
    """Return one audit row per species/prefix pair in a curated CSV."""
    raw_by_accession: dict[str, dict[str, str]] = {}
    if raw_path.exists():
        with open(raw_path, encoding="utf-8") as handle:
            raw_by_accession = {row["uniprot_accession"]: row for row in csv.DictReader(handle)}

    groups: dict[tuple[str, str], dict[str, object]] = defaultdict(lambda: {"rows": 0, "source_rows": 0, "accessions": [], "labels": []})
    with open(curated_path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            gene = row.get("gene", "")
            if not gene or gene.startswith("~") or "-" not in gene:
                continue
            prefix = gene.split("-", 1)[0]
            species = extract_latin_binomial(row.get("organism", ""))
            group = groups[(species, prefix)]
            group["rows"] = int(group["rows"]) + 1

            accession = row.get("uniprot_accession", "")
            raw = raw_by_accession.get(accession)
            if raw is None:
                continue
            label = " ".join((raw.get("gene_names", ""), raw.get("protein_name", ""))).strip()
            if not prefix_is_present(prefix, label):
                continue
            group["source_rows"] = int(group["source_rows"]) + 1
            accessions = group["accessions"]
            labels = group["labels"]
            if isinstance(accessions, list) and accession not in accessions and len(accessions) < 3:
                accessions.append(accession)
            if isinstance(labels, list) and label not in labels and len(labels) < 3:
                labels.append(label)

    result: list[dict[str, str | int]] = []
    for (species, prefix), group in groups.items():
        accessions = group["accessions"]
        labels = group["labels"]
        assert isinstance(accessions, list)
        assert isinstance(labels, list)
        status, evidence = classify_prefix(species, prefix, accessions)
        result.append(
            {
                "species": species,
                "legacy_prefix": prefix,
                "status": status,
                "row_count": int(group["rows"]),
                "source_attested_row_count": int(group["source_rows"]),
                "evidence": evidence,
                "example_accessions": ";".join(accessions),
                "example_source_labels": ";".join(labels),
            }
        )
    return sorted(result, key=lambda row: (str(row["status"]), str(row["species"]), str(row["legacy_prefix"])))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curated", type=Path, default=DEFAULT_CURATED)
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows = audit_prefixes(args.curated, args.raw)
    counts = Counter(str(row["status"]) for row in rows)
    print(f"Audited {len(rows)} species/prefix pairs from {args.curated}")
    for status, count in sorted(counts.items()):
        print(f"  {status}: {count}")

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=AUDIT_FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
