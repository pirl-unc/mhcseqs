#!/usr/bin/env python3
"""Audit stored SP taxa against UniProt taxonomy and persist source clades."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.enrich_sp_ground_truth import CONTROL_FIELDS, ENRICHED_FIELDS, GT_ENRICHED_CSV, NEGATIVE_CONTROL_CSV
from scripts.fetch_sp_ground_truth import FIELDS as RAW_FIELDS
from scripts.fetch_sp_ground_truth import OUTPUT as RAW_CSV
from scripts.sp_ground_truth_taxonomy import (
    SOURCE_CLADES,
    TAXONOMY_AUDIT_CSV,
    category_from_lineage,
    source_clade_from_lineage,
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


def _fetch_taxonomy(taxon_id: int, attempts: int = 4) -> tuple[str, dict, str]:
    """Return the request URL, decoded payload, and UniProt release for a taxon.

    Retries with exponential backoff: an audit is 253 requests against a
    rate-limited public endpoint, and a single transient 429/503 would
    otherwise discard the whole run.
    """
    url = f"https://rest.uniprot.org/taxonomy/{taxon_id}"
    request = urllib.request.Request(url, headers={"User-Agent": "mhcseqs-sp-taxonomy-audit/1.0"})
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = json.loads(response.read())
                release = response.headers.get("X-UniProt-Release", "")
            return url, payload, release
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            if attempt == attempts - 1:
                raise RuntimeError(f"UniProt taxonomy lookup failed for {taxon_id} after {attempts} attempts") from exc
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def _verify_source_clade_roots() -> None:
    """Check every hardcoded clade name against the UniProt scientific name.

    The leaf taxa are resolved from UniProt, so the six lineage roots that
    define the benchmark must be held to the same standard instead of being
    trusted as handwritten labels.
    """
    for name, taxon_id in SOURCE_CLADES:
        _url, payload, _release = _fetch_taxonomy(taxon_id)
        scientific_name = str(payload.get("scientificName", ""))
        if scientific_name != name:
            raise ValueError(f"Source clade {name!r} does not match UniProt name {scientific_name!r} for taxon {taxon_id}")
    print(f"Verified {len(SOURCE_CLADES)} source-clade root names against UniProt")


def _fetch_taxon(item: tuple[str, str]) -> dict[str, str]:
    taxon_id_text, organism = item
    taxon_id = int(taxon_id_text)
    url, payload, release = _fetch_taxonomy(taxon_id)
    scientific_name = str(payload.get("scientificName", ""))
    if scientific_name != organism:
        raise ValueError(f"Taxon {taxon_id} is stored as {organism!r} but UniProt now reports {scientific_name!r}")
    lineage_ids = {int(node["taxonId"]) for node in payload.get("lineage", [])}
    source_clade, source_root = source_clade_from_lineage(taxon_id, lineage_ids)
    category, category_root = category_from_lineage(taxon_id, lineage_ids)
    return {
        "taxon_id": taxon_id_text,
        "organism": organism,
        "taxonomy_name": scientific_name,
        "source_clade": source_clade,
        "source_clade_taxon_id": str(source_root),
        "species_category": category,
        "category_basis_taxon_id": str(category_root),
        "taxonomy_release": release,
        "source_url": url,
        "audited_on": date.today().isoformat(),
    }


def _write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=RAW_CSV)
    parser.add_argument("--enriched", type=Path, default=GT_ENRICHED_CSV)
    parser.add_argument("--controls", type=Path, default=NEGATIVE_CONTROL_CSV)
    parser.add_argument("--output", type=Path, default=TAXONOMY_AUDIT_CSV)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--reuse-audit",
        action="store_true",
        help="Apply the existing audit without querying UniProt again.",
    )
    args = parser.parse_args()

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
        _verify_source_clade_roots()
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            audit_rows = list(pool.map(_fetch_taxon, items))
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

    print(f"Audited {len(audit_rows)} taxa across {len(raw_rows)} raw and {len(enriched_rows)} enriched rows against UniProt {releases.pop()}")
    for clade, _taxon_id in SOURCE_CLADES:
        count = sum(row["source_clade"] == clade for row in raw_rows)
        print(f"  {clade}: {count}")
    print(f"{'Used' if args.reuse_audit else 'Wrote'} {args.output}")
    print(f"Updated {args.raw}")
    print(f"Updated {args.enriched}")
    print(f"Updated {args.controls}")


if __name__ == "__main__":
    main()
