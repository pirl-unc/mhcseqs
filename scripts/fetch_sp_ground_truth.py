#!/usr/bin/env python3
"""Derive signal-peptide ground truth from saved UniProt artifacts.

Network access is explicit: ``--refresh`` downloads one broad, release-pinned
Vertebrata candidate artifact and complete UniSave records for historical
benchmark accessions absent from that current stream. Every run then derives
the ground-truth CSV from local artifacts only.

Output: data/sp_ground_truth.csv

Usage:
    python scripts/fetch_sp_ground_truth.py --refresh
    python scripts/fetch_sp_ground_truth.py
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections.abc import Sequence
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.sp_corpus_artifacts import (
    download_corpus_artifact,
    download_deleted_artifact,
    load_corpus_artifact,
    load_deleted_artifact,
    source_clade_from_artifact_row,
    write_manifest,
)

OUTPUT = ROOT / "data" / "sp_ground_truth.csv"

FIELDS = [
    "accession",
    "organism",
    "taxon_id",
    "source_clade",
    "sp_length",
    "reviewed",
    "sequence",
]

# Preserve the benchmark's historical local eligibility window while moving it
# out of the network query. Changing this is now an offline, reviewable choice.
MIN_LENGTH = 100
MAX_LENGTH = 500
_EXACT_SIGNAL_RE = re.compile(r"(?:^|;\s*)SIGNAL\s+1\.\.(\d+)(?:;|$)")


def _read_ground_truth(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _signal_end(feature: str) -> int | None:
    """Return an exact N-terminal SIGNAL endpoint, rejecting fuzzy bounds."""
    match = _EXACT_SIGNAL_RE.search(feature)
    return int(match.group(1)) if match else None


def artifact_row_to_ground_truth(row: dict[str, str]) -> dict[str, str] | None:
    """Convert one eligible artifact row to the established raw GT schema."""
    sequence = row.get("Sequence", "").strip()
    if not sequence or not MIN_LENGTH <= len(sequence) <= MAX_LENGTH:
        return None
    sp_length = _signal_end(row.get("Signal peptide", ""))
    if sp_length is None or not 0 < sp_length < len(sequence):
        return None
    accession = row.get("Entry", "").strip()
    organism = row.get("Organism", "").strip()
    taxon_id = row.get("Organism (ID)", "").strip()
    if not accession or not organism or not taxon_id:
        raise ValueError(f"Artifact row lacks benchmark identity fields: {accession or '<missing accession>'}")
    return {
        "accession": accession,
        "organism": organism,
        "taxon_id": taxon_id,
        "source_clade": source_clade_from_artifact_row(row),
        "sp_length": str(sp_length),
        "reviewed": "Y" if row.get("Reviewed", "").strip().lower() == "reviewed" else "N",
        "sequence": sequence,
    }


def build_ground_truth_rows(artifact_rows: Sequence[dict[str, str]]) -> list[dict[str, str]]:
    """Derive and validate benchmark rows from current plus archived records."""
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for artifact_row in artifact_rows:
        row = artifact_row_to_ground_truth(artifact_row)
        if row is None:
            continue
        accession = row["accession"]
        if accession in seen:
            raise ValueError(f"Duplicate artifact accession in ground truth: {accession}")
        seen.add(accession)
        rows.append(row)
    return sorted(rows, key=lambda row: row["accession"])


def refresh_artifacts(
    preserved_rows: Sequence[dict[str, str]],
    data_dir: Path | None = None,
) -> None:
    """Download current candidates and full archives for disappeared GT rows."""
    corpus = download_corpus_artifact(data_dir)
    current_accessions = {row["Entry"] for row in load_corpus_artifact(data_dir)}
    missing_rows = [row for row in preserved_rows if row["accession"] not in current_accessions]
    source_clades = {row["accession"]: row.get("source_clade", "") for row in missing_rows}
    deleted = download_deleted_artifact(
        [row["accession"] for row in missing_rows],
        data_dir,
        source_clades=source_clades,
    )
    manifest = write_manifest([corpus, deleted], data_dir)
    print(f"Downloaded {corpus.rows:,} current candidates from UniProt {corpus.release or '<unknown>'}")
    print(f"Archived {deleted.rows:,} disappeared benchmark entries from UniSave")
    print(f"Wrote artifact manifest to {manifest}")


def write_ground_truth(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write the established raw benchmark schema."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Download current UniProt and missing-entry UniSave artifacts before deriving the CSV.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        help="Cache root (default: MHCSEQS_DATA or ~/.cache/mhcseqs).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT,
        help=f"Derived CSV path (default: {OUTPUT}).",
    )
    parser.add_argument(
        "--preserve-from",
        type=Path,
        default=OUTPUT,
        help="Existing benchmark whose disappeared accessions must be archived.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = _parse_args(argv)
    preserved_rows = _read_ground_truth(args.preserve_from)
    if args.refresh:
        refresh_artifacts(preserved_rows, args.data_dir)

    artifact_rows = load_corpus_artifact(args.data_dir) + load_deleted_artifact(args.data_dir)
    rows = build_ground_truth_rows(artifact_rows)
    write_ground_truth(args.output, rows)
    print(f"Wrote {len(rows):,} ground-truth entries to {args.output}")


if __name__ == "__main__":
    main()
