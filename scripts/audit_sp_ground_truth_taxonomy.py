#!/usr/bin/env python3
"""Audit SP taxa against the release-pinned taxonomy artifact."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.enrich_sp_ground_truth import CONTROL_FIELDS, ENRICHED_FIELDS, GT_ENRICHED_CSV, NEGATIVE_CONTROL_CSV
from scripts.fetch_sp_ground_truth import FIELDS as RAW_FIELDS
from scripts.fetch_sp_ground_truth import OUTPUT as RAW_CSV
from scripts.sp_corpus_artifacts import taxonomy_by_id, validate_artifact_bundle
from scripts.sp_ground_truth_taxonomy import (
    SOURCE_CLADES,
    TAXONOMY_AUDIT_CSV,
    classification_from_taxonomy_row,
)

AUDIT_FIELDS = [
    "taxon_id",
    "organism",
    "taxonomy_name",
    "source_clade",
    "source_clade_taxon_id",
    "species_category",
    "category_basis_taxon_id",
    "taxonomy_release",
    "source_url",
    "audited_on",
]


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _verify_source_clade_roots(taxonomy: dict[str, dict[str, str]]) -> None:
    """Check every source-clade label against the pinned taxonomy snapshot."""
    for name, taxon_id in SOURCE_CLADES:
        scientific_name = taxonomy[str(taxon_id)]["Scientific name"]
        if scientific_name != name:
            raise ValueError(f"Source clade {name!r} does not match UniProt name {scientific_name!r} for taxon {taxon_id}")
    print(f"Verified {len(SOURCE_CLADES)} source-clade root names against the taxonomy snapshot")


def _audit_taxon(item: tuple[str, str], taxonomy: dict[str, dict[str, str]]) -> dict[str, str]:
    taxon_id_text, organism = item
    snapshot = taxonomy.get(taxon_id_text)
    if snapshot is None:
        raise ValueError(f"Taxonomy artifact has no record for taxon {taxon_id_text}")
    scientific_name = snapshot["Scientific name"]
    if scientific_name != organism:
        raise ValueError(f"Taxon {taxon_id_text} is stored as {organism!r} but the artifact reports {scientific_name!r}")
    source_clade, source_root, category, category_root = classification_from_taxonomy_row(snapshot)
    return {
        "taxon_id": taxon_id_text,
        "organism": organism,
        "taxonomy_name": scientific_name,
        "source_clade": source_clade,
        "source_clade_taxon_id": str(source_root),
        "species_category": category,
        "category_basis_taxon_id": str(category_root),
        "taxonomy_release": snapshot["Taxonomy release"],
        "source_url": snapshot["Source URL"],
        "audited_on": snapshot["Retrieved on"],
    }


def _write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _apply_audit(rows: list[dict[str, str]], audit: dict[str, dict[str, str]]) -> None:
    """Reconcile taxonomy fields against the audit, leaving labels untouched.

    A stored source clade that disagrees with the lineage-derived one means the
    fetch query and the returned lineage have diverged (a reclassified taxon, or
    a root alias change like Reptilia -> Lepidosauria). Overwriting it silently
    would erase exactly the signal this audit exists to surface, so it raises.
    """
    for row in rows:
        decision = audit.get(row["taxon_id"])
        if decision is None:
            raise ValueError(f"Missing taxonomy audit for taxon {row['taxon_id']}")
        stored_clade = row.get("source_clade", "")
        if stored_clade and stored_clade != decision["source_clade"]:
            raise ValueError(
                f"Source clade drift for taxon {row['taxon_id']}: stored={stored_clade!r}, lineage={decision['source_clade']!r}. "
                f"Re-fetch the corpus if the clade genuinely moved."
            )
        row["source_clade"] = decision["source_clade"]
        if "species_category" in row:
            row["species_category"] = decision["species_category"]


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=RAW_CSV)
    parser.add_argument("--enriched", type=Path, default=GT_ENRICHED_CSV)
    parser.add_argument("--controls", type=Path, default=NEGATIVE_CONTROL_CSV)
    parser.add_argument("--output", type=Path, default=TAXONOMY_AUDIT_CSV)
    parser.add_argument("--data-dir", type=Path, help="Artifact cache root.")
    parser.add_argument(
        "--reuse-audit",
        action="store_true",
        help="Apply the existing audit instead of rebuilding it from the artifact bundle.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)

    raw_rows = _read_csv(args.raw)
    organisms_by_taxon: dict[str, str] = {}
    for row in raw_rows:
        taxon_id = row["taxon_id"].strip()
        organism = row["organism"].strip()
        previous = organisms_by_taxon.setdefault(taxon_id, organism)
        if previous != organism:
            raise ValueError(f"Taxon {taxon_id} has multiple organism names: {previous!r}, {organism!r}")

    items = sorted(organisms_by_taxon.items(), key=lambda item: int(item[0]))
    if args.reuse_audit:
        audit_rows = _read_csv(args.output)
    else:
        validate_artifact_bundle(args.data_dir)
        taxonomy = taxonomy_by_id(args.data_dir)
        _verify_source_clade_roots(taxonomy)
        audit_rows = [_audit_taxon(item, taxonomy) for item in items]
    releases = {row["taxonomy_release"] for row in audit_rows}
    if len(releases) != 1 or "" in releases:
        raise ValueError(f"Expected one non-empty UniProt release, found {releases!r}")

    audit_by_taxon = {row["taxon_id"]: row for row in audit_rows}
    if set(audit_by_taxon) != set(organisms_by_taxon):
        raise ValueError("Taxonomy audit and raw corpus contain different taxon IDs")
    audited_clades = {row["source_clade"] for row in audit_rows}
    known_clades = {name for name, _taxon_id in SOURCE_CLADES}
    if not audited_clades <= known_clades:
        raise ValueError(f"Taxonomy audit uses unknown source clades: {sorted(audited_clades - known_clades)!r}")

    enriched_rows = _read_csv(args.enriched)
    control_rows = _read_csv(args.controls)
    _apply_audit(raw_rows, audit_by_taxon)
    _apply_audit(enriched_rows, audit_by_taxon)
    _apply_audit(control_rows, audit_by_taxon)

    if not args.reuse_audit:
        _write_csv(args.output, audit_rows, AUDIT_FIELDS)
    _write_csv(args.raw, raw_rows, RAW_FIELDS)
    _write_csv(args.enriched, enriched_rows, ENRICHED_FIELDS)
    _write_csv(args.controls, control_rows, CONTROL_FIELDS)

    print(
        f"Audited {len(audit_rows)} taxa across {len(raw_rows)} raw and {len(enriched_rows)} enriched rows "
        f"against the pinned UniProt {releases.pop()} snapshot"
    )
    for clade, _taxon_id in SOURCE_CLADES:
        count = sum(row["source_clade"] == clade for row in raw_rows)
        print(f"  {clade}: {count}")
    print(f"{'Used' if args.reuse_audit else 'Wrote'} {args.output}")
    print(f"Updated {args.raw}")
    print(f"Updated {args.enriched}")
    print(f"Updated {args.controls}")


if __name__ == "__main__":
    main()
