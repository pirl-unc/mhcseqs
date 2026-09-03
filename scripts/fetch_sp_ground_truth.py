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
import tempfile
from collections import Counter
from collections.abc import Sequence
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.sp_corpus_artifacts import (
    artifact_dir,
    download_corpus_artifact,
    download_deleted_artifact,
    download_taxonomy_artifact,
    load_corpus_artifact,
    load_deleted_artifact,
    publish_artifact_bundle,
    required_taxonomy_ids,
    taxonomy_by_id,
    validate_artifact_bundle,
    write_manifest,
)
from scripts.sp_ground_truth_eligibility import (
    GT_LABEL_CURATION_CSV,
    load_label_curation,
    resolve_mhc_label,
)
from scripts.sp_ground_truth_taxonomy import classification_from_taxonomy_row

OUTPUT = ROOT / "data" / "sp_ground_truth.csv"
UNLABELLED_OUTPUT = ROOT / "data" / "sp_unlabelled_sequences.csv"

FIELDS = [
    "accession",
    "organism",
    "taxon_id",
    "source_clade",
    "sp_length",
    "reviewed",
    "sequence",
]
UNLABELLED_FIELDS = [
    "accession",
    "organism",
    "taxon_id",
    "source_clade",
    "reviewed",
    "mhc_class",
    "chain",
    "protein_name",
    "gene_names",
    "metadata_source",
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


def _identity_fields(
    row: dict[str, str],
    taxonomy: dict[str, dict[str, str]],
) -> dict[str, str]:
    accession = row.get("Entry", "").strip()
    taxon_id = row.get("Organism (ID)", "").strip()
    decision = taxonomy.get(taxon_id)
    if not accession or not taxon_id or decision is None:
        raise ValueError(f"Artifact row lacks benchmark identity or taxonomy: {accession or '<missing accession>'}")
    organism = decision.get("Scientific name", "").strip()
    if not organism:
        raise ValueError(f"Taxonomy artifact has no scientific name for taxon {taxon_id}")
    return {
        "accession": accession,
        "organism": organism,
        "taxon_id": taxon_id,
        "source_clade": classification_from_taxonomy_row(decision)[0],
        "reviewed": "Y" if row.get("Reviewed", "").strip().lower() == "reviewed" else "N",
        "sequence": row.get("Sequence", "").strip(),
    }


def artifact_row_to_ground_truth(
    row: dict[str, str],
    taxonomy: dict[str, dict[str, str]],
    curation: dict[str, dict[str, str]] | None = None,
) -> dict[str, str] | None:
    """Convert one eligible artifact row to the established raw GT schema."""
    identity = _identity_fields(row, taxonomy)
    sequence = identity["sequence"]
    if not sequence or not MIN_LENGTH <= len(sequence) <= MAX_LENGTH:
        return None
    if row.get("Fragment", "").strip():
        return None
    if not resolve_mhc_label(row, curation or {}).eligible:
        return None
    sp_length = _signal_end(row.get("Signal peptide", ""))
    if sp_length is None or not 0 < sp_length < len(sequence):
        return None
    return {
        **{field: identity[field] for field in FIELDS if field not in {"sp_length"}},
        "sp_length": str(sp_length),
    }


def derive_corpus_rows(
    artifact_rows: Sequence[dict[str, str]],
    taxonomy: dict[str, dict[str, str]],
    curation: dict[str, dict[str, str]] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]], Counter[str]]:
    """Derive labelled and unlabelled MHC chains with disjoint audit counts."""
    curation = curation or {}
    labelled: list[dict[str, str]] = []
    unlabelled: list[dict[str, str]] = []
    stats: Counter[str] = Counter()
    seen: set[str] = set()
    for artifact_row in artifact_rows:
        stats["candidates"] += 1
        identity = _identity_fields(artifact_row, taxonomy)
        accession = identity["accession"]
        if accession in seen:
            raise ValueError(f"Duplicate artifact accession in corpus: {accession}")
        seen.add(accession)

        sequence = identity["sequence"]
        if artifact_row.get("Fragment", "").strip():
            stats["excluded_fragment"] += 1
            continue
        if not sequence or not MIN_LENGTH <= len(sequence) <= MAX_LENGTH:
            stats["excluded_length"] += 1
            continue

        label = resolve_mhc_label(artifact_row, curation)
        if not label.eligible:
            key = "excluded_non_mhc" if label.disposition == "exclude_non_mhc" else "excluded_unresolved"
            stats[key] += 1
            continue

        sp_length = _signal_end(artifact_row.get("Signal peptide", ""))
        if sp_length is not None and 0 < sp_length < len(sequence):
            labelled.append(
                {
                    **{field: identity[field] for field in FIELDS if field != "sp_length"},
                    "sp_length": str(sp_length),
                }
            )
            stats["labelled"] += 1
            continue
        if artifact_row.get("Signal peptide", "").strip():
            stats["excluded_invalid_signal"] += 1
            continue
        if not sequence.startswith("M"):
            stats["excluded_incomplete_n_terminus"] += 1
            continue
        unlabelled.append(
            {
                **identity,
                "mhc_class": label.mhc_class,
                "chain": label.chain,
                "protein_name": artifact_row.get("Protein names", "").strip(),
                "gene_names": artifact_row.get("Gene Names", "").strip(),
                "metadata_source": "unisave_artifact" if artifact_row.get("Archive source URL", "").strip() else "uniprot_corpus_artifact",
            }
        )
        stats["unlabelled"] += 1

    return (
        sorted(labelled, key=lambda row: row["accession"]),
        sorted(unlabelled, key=lambda row: row["accession"]),
        stats,
    )


def build_ground_truth_rows(
    artifact_rows: Sequence[dict[str, str]],
    taxonomy: dict[str, dict[str, str]],
    curation: dict[str, dict[str, str]] | None = None,
) -> list[dict[str, str]]:
    """Derive and validate benchmark rows from current plus archived records."""
    rows, _unlabelled, _stats = derive_corpus_rows(artifact_rows, taxonomy, curation)
    return rows


def refresh_artifacts(
    preserved_rows: Sequence[dict[str, str]],
    data_dir: Path | None = None,
) -> None:
    """Stage, validate, and publish one release-consistent artifact bundle."""
    destination = artifact_dir(data_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".sp-refresh-", dir=destination.parent) as tmp:
        staged_data_dir = Path(tmp)
        corpus = download_corpus_artifact(staged_data_dir)
        corpus_rows = load_corpus_artifact(staged_data_dir)
        current_accessions = {row["Entry"] for row in corpus_rows}
        missing_rows = [row for row in preserved_rows if row["accession"] not in current_accessions]
        source_clades = {row["accession"]: row.get("source_clade", "") for row in missing_rows}
        deleted = download_deleted_artifact(
            [row["accession"] for row in missing_rows],
            staged_data_dir,
            source_clades=source_clades,
            expected_release=corpus.release,
        )
        deleted_rows = load_deleted_artifact(staged_data_dir)
        taxonomy = download_taxonomy_artifact(
            required_taxonomy_ids(corpus_rows + deleted_rows),
            staged_data_dir,
            expected_release=corpus.release,
        )
        write_manifest([corpus, deleted, taxonomy], staged_data_dir)
        validate_artifact_bundle(staged_data_dir)
        manifest = publish_artifact_bundle(staged_data_dir, data_dir)
    print(f"Downloaded {corpus.rows:,} current candidates from UniProt {corpus.release or '<unknown>'}")
    print(f"Archived {deleted.rows:,} disappeared benchmark entries from UniSave")
    print(f"Snapshotted {taxonomy.rows:,} taxonomy records")
    print(f"Wrote artifact manifest to {manifest}")


def write_ground_truth(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write the established raw benchmark schema."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_unlabelled(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write full-length MHC chains without a usable signal annotation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=UNLABELLED_FIELDS, lineterminator="\n")
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
        "--unlabelled-output",
        type=Path,
        default=UNLABELLED_OUTPUT,
        help=f"Unlabelled inference-set path (default: {UNLABELLED_OUTPUT}).",
    )
    parser.add_argument(
        "--label-curation",
        type=Path,
        default=GT_LABEL_CURATION_CSV,
        help="Accession-level source-backed class/chain decisions.",
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

    validate_artifact_bundle(args.data_dir)
    artifact_rows = load_corpus_artifact(args.data_dir) + load_deleted_artifact(args.data_dir)
    taxonomy = taxonomy_by_id(args.data_dir)
    curation = load_label_curation(args.label_curation)
    rows, unlabelled, stats = derive_corpus_rows(artifact_rows, taxonomy, curation)
    write_ground_truth(args.output, rows)
    write_unlabelled(args.unlabelled_output, unlabelled)
    print(f"Wrote {len(rows):,} ground-truth entries to {args.output}")
    print(f"Wrote {len(unlabelled):,} unlabelled MHC entries to {args.unlabelled_output}")
    print("Corpus decisions:")
    for key in (
        "candidates",
        "labelled",
        "unlabelled",
        "excluded_fragment",
        "excluded_length",
        "excluded_non_mhc",
        "excluded_unresolved",
        "excluded_invalid_signal",
        "excluded_incomplete_n_terminus",
    ):
        print(f"  {key}: {stats[key]:,}")


if __name__ == "__main__":
    main()
