#!/usr/bin/env python3
"""Build the versioned, record-complete UniProt MHC protein dataset.

The source bundle is deliberately broader than any benchmark view. Every
current or historically preserved candidate is emitted exactly once. Source
annotations remain nullable source facts, while normalized MHC labels and
mhcseqs parser outputs live in separate column families.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import sys
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mhcseqs.domain_grammar import SP_BOUNDARY_MODEL_PATH, SP_SEQUENCE_CUE_MODEL_PATH
from mhcseqs.domain_parsing import AlleleRecord, decompose_domains
from mhcseqs.version import __version__
from scripts.enrich_sp_ground_truth import _classify_from_names
from scripts.sp_corpus_artifacts import (
    load_corpus_artifact,
    load_deleted_artifact,
    read_manifest,
    taxonomy_by_id,
    validate_artifact_bundle,
)
from scripts.sp_ground_truth_eligibility import GT_LABEL_CURATION_CSV, load_label_curation, resolve_mhc_label
from scripts.sp_ground_truth_taxonomy import classification_from_taxonomy_row

DATASET_NAME = "mhc-proteins"
DATASET_REVISION = 1
DATASET_SCHEMA_VERSION = 1
DATASET_FILENAME = "mhc-proteins-{version}.csv.gz"
MANIFEST_FILENAME = "mhc-proteins-{version}.manifest.json"

PROTEIN_FIELDS = (
    "dataset_version",
    "accession",
    "source_kind",
    "source_release",
    "source_url",
    "source_entry_version",
    "source_modified_date",
    "source_archive_release",
    "source_archive_release_date",
    "source_retrieved_on",
    "reviewed",
    "entry_name",
    "organism",
    "taxon_id",
    "taxonomic_lineage_ids",
    "taxonomic_lineage",
    "source_clade",
    "species_category",
    "protein_name",
    "gene_names",
    "protein_existence",
    "keywords",
    "source_fragment",
    "source_length",
    "sequence",
    "source_signal_feature",
    "source_signal_status",
    "source_signal_start",
    "source_signal_end",
    "source_chain_feature",
    "inferred_mhc_class",
    "inferred_chain",
    "inferred_protein_type",
    "inferred_gene",
    "inferred_gene_status",
    "inferred_label_status",
    "inferred_disposition",
    "parse_status",
    "parse_ok",
    "parsed_mhc_class",
    "parsed_chain",
    "parsed_protein_type",
    "parsed_signal_status",
    "parsed_signal_end",
    "parsed_signal_sequence",
    "mature_sequence",
    "groove1",
    "groove2",
    "groove_seq",
    "groove1_len",
    "groove2_len",
    "ig_domain",
    "ig_domain_len",
    "tail",
    "tail_len",
    "domain_architecture",
    "domain_spans",
    "parse_score",
    "sp_subscore",
    "groove_subscore",
    "ig_subscore",
    "tail_subscore",
    "anchor_type",
    "anchor_cys1",
    "anchor_cys2",
    "secondary_cys1",
    "secondary_cys2",
    "candidate_type",
    "nterm_state",
    "support_state",
    "tail_state",
    "parse_flags",
    "parse_error",
)

_SOURCE_SIGNAL_RE = re.compile(r"(?:^|;\s*)SIGNAL\s+([<>?]?\d+)\.\.([<>?]?\d+)(?:;|$)")


def dataset_version(source_release: str, revision: int = DATASET_REVISION) -> str:
    """Return the independent data version for one source release."""
    release = str(source_release).strip()
    if not release:
        raise ValueError("Cannot version an MHC protein dataset without a source release")
    return f"uniprot-{release}-r{revision}"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def _protein_type(mhc_class: str, chain: str) -> str:
    mhc_class = str(mhc_class or "").strip().upper()
    chain = str(chain or "").strip().lower()
    if mhc_class not in {"I", "II"}:
        return ""
    return f"mhc_{mhc_class.lower()}_{chain or 'unknown'}"


def _source_signal_fields(feature: str) -> dict[str, str]:
    """Normalize a nullable UniProt SIGNAL feature without inventing a label."""
    raw = str(feature or "").strip()
    if not raw:
        return {
            "source_signal_feature": "",
            "source_signal_status": "not_annotated",
            "source_signal_start": "",
            "source_signal_end": "",
        }
    match = _SOURCE_SIGNAL_RE.search(raw)
    if match is None:
        return {
            "source_signal_feature": raw,
            "source_signal_status": "unparsed",
            "source_signal_start": "",
            "source_signal_end": "",
        }
    start_text, end_text = match.groups()
    if not start_text.isdigit() or not end_text.isdigit():
        return {
            "source_signal_feature": raw,
            "source_signal_status": "fuzzy",
            "source_signal_start": "",
            "source_signal_end": "",
        }
    start = int(start_text)
    end = int(end_text)
    if start < 1 or end < start:
        return {
            "source_signal_feature": raw,
            "source_signal_status": "invalid",
            "source_signal_start": "",
            "source_signal_end": "",
        }
    return {
        "source_signal_feature": raw,
        "source_signal_status": "exact",
        # All normalized coordinate columns use zero-based half-open offsets.
        "source_signal_start": str(start - 1),
        "source_signal_end": str(end),
    }


def _format_score(value: float) -> str:
    return format(float(value), ".8g")


def _optional_index(value: int | None) -> str:
    return "" if value is None else str(value)


def _parser_fields(result: AlleleRecord | None, error: str = "") -> dict[str, str]:
    if result is None:
        return {
            "parse_status": "parser_error",
            "parse_ok": "false",
            "parsed_mhc_class": "",
            "parsed_chain": "",
            "parsed_protein_type": "",
            "parsed_signal_status": "unresolved",
            "parsed_signal_end": "",
            "parsed_signal_sequence": "",
            "mature_sequence": "",
            "groove1": "",
            "groove2": "",
            "groove_seq": "",
            "groove1_len": "",
            "groove2_len": "",
            "ig_domain": "",
            "ig_domain_len": "",
            "tail": "",
            "tail_len": "",
            "domain_architecture": "",
            "domain_spans": "",
            "parse_score": "",
            "sp_subscore": "",
            "groove_subscore": "",
            "ig_subscore": "",
            "tail_subscore": "",
            "anchor_type": "",
            "anchor_cys1": "",
            "anchor_cys2": "",
            "secondary_cys1": "",
            "secondary_cys2": "",
            "candidate_type": "",
            "nterm_state": "",
            "support_state": "",
            "tail_state": "",
            "parse_flags": "",
            "parse_error": error,
        }

    ok = result.ok
    if result.mature_start > 0:
        signal_status = "present"
        signal_end = str(result.mature_start)
    elif ok:
        signal_status = "not_detected"
        signal_end = "0"
    else:
        signal_status = "unresolved"
        signal_end = ""
    signal_sequence = result.sequence[: result.mature_start] if result.mature_start > 0 else ""
    mature_sequence = result.mature_sequence if signal_status != "unresolved" else ""

    return {
        "parse_status": result.status,
        "parse_ok": str(ok).lower(),
        "parsed_mhc_class": result.mhc_class,
        "parsed_chain": result.chain,
        "parsed_protein_type": _protein_type(result.mhc_class, result.chain),
        "parsed_signal_status": signal_status,
        "parsed_signal_end": signal_end,
        "parsed_signal_sequence": signal_sequence,
        "mature_sequence": mature_sequence,
        "groove1": result.groove1,
        "groove2": result.groove2,
        "groove_seq": result.groove_seq,
        "groove1_len": str(result.groove1_len),
        "groove2_len": str(result.groove2_len),
        "ig_domain": result.ig_domain,
        "ig_domain_len": str(result.ig_domain_len),
        "tail": result.tail,
        "tail_len": str(result.tail_len),
        "domain_architecture": result.domain_architecture,
        "domain_spans": result.domain_spans,
        "parse_score": _format_score(result.parse_score),
        "sp_subscore": _format_score(result.sp_subscore),
        "groove_subscore": _format_score(result.groove_subscore),
        "ig_subscore": _format_score(result.ig_subscore),
        "tail_subscore": _format_score(result.tail_subscore),
        "anchor_type": result.anchor_type,
        "anchor_cys1": _optional_index(result.anchor_cys1),
        "anchor_cys2": _optional_index(result.anchor_cys2),
        "secondary_cys1": _optional_index(result.secondary_cys1),
        "secondary_cys2": _optional_index(result.secondary_cys2),
        "candidate_type": result.candidate_type,
        "nterm_state": result.nterm_state,
        "support_state": result.support_state,
        "tail_state": result.tail_state,
        "parse_flags": ",".join(result.flags),
        "parse_error": error,
    }


def _parse_prepared_record(prepared: dict[str, str]) -> dict[str, str]:
    """Run mhcseqs parsing for one source-normalized record."""
    row = {key: value for key, value in prepared.items() if not key.startswith("_parse_")}
    sequence = row["sequence"]
    try:
        result = decompose_domains(
            sequence,
            mhc_class=prepared["_parse_mhc_class"],
            chain=prepared["_parse_chain"] or None,
            allele=row["accession"],
            gene=row["inferred_gene"],
            species=row["organism"] or None,
            use_early_shortcuts=prepared["_parse_use_early_shortcuts"] == "true",
        )
        row.update(_parser_fields(result))
    except Exception as exc:  # Preserve the record and make parser failures auditable.
        row.update(_parser_fields(None, f"{type(exc).__name__}: {exc}"))
    return row


def prepare_record(
    artifact: dict[str, str],
    *,
    source_release: str,
    version: str,
    taxonomy: dict[str, dict[str, str]],
    curation: dict[str, dict[str, str]],
    use_early_shortcuts: bool = False,
) -> dict[str, str]:
    """Preserve source fields and add normalized metadata for parser input."""
    accession = artifact.get("Entry", "").strip()
    taxon_id = artifact.get("Organism (ID)", "").strip()
    taxonomy_row = taxonomy.get(taxon_id)
    if not accession or taxonomy_row is None:
        raise ValueError(f"Artifact record lacks accession or taxonomy: {accession or '<missing>'}")

    source_clade, _source_root, species_category, _category_root = classification_from_taxonomy_row(taxonomy_row)
    organism = taxonomy_row["Scientific name"].strip()
    protein_name = artifact.get("Protein names", "").strip()
    gene_names = artifact.get("Gene Names", "").strip()
    label = resolve_mhc_label(artifact, curation)
    classified = _classify_from_names(organism=organism, protein_name=protein_name, gene_names=gene_names)
    inferred_gene = classified["gene"]
    archived = bool(artifact.get("Archive source URL", "").strip())
    source_url = artifact.get("Archive source URL", "").strip() if archived else f"https://rest.uniprot.org/uniprotkb/{accession}"
    row = {
        "dataset_version": version,
        "accession": accession,
        "source_kind": "uniprot_unisave" if archived else "uniprotkb",
        "source_release": source_release,
        "source_url": source_url,
        "source_entry_version": artifact.get("Entry version", "").strip(),
        "source_modified_date": artifact.get("Date of last modification", "").strip(),
        "source_archive_release": artifact.get("Archive release", "").strip(),
        "source_archive_release_date": artifact.get("Archive release date", "").strip(),
        "source_retrieved_on": taxonomy_row.get("Retrieved on", "").strip(),
        "reviewed": "true" if artifact.get("Reviewed", "").strip().lower() == "reviewed" else "false",
        "entry_name": artifact.get("Entry Name", "").strip(),
        "organism": organism,
        "taxon_id": taxon_id,
        "taxonomic_lineage_ids": taxonomy_row.get("Taxonomic lineage (Ids)", "").strip(),
        "taxonomic_lineage": taxonomy_row.get("Taxonomic lineage", "").strip(),
        "source_clade": source_clade,
        "species_category": species_category,
        "protein_name": protein_name,
        "gene_names": gene_names,
        "protein_existence": artifact.get("Protein existence", "").strip(),
        "keywords": artifact.get("Keywords", "").strip(),
        "source_fragment": str(bool(artifact.get("Fragment", "").strip())).lower(),
        "source_length": artifact.get("Length", "").strip(),
        "sequence": artifact.get("Sequence", "").strip().upper(),
        **_source_signal_fields(artifact.get("Signal peptide", "")),
        "source_chain_feature": artifact.get("Chain", "").strip(),
        "inferred_mhc_class": label.mhc_class,
        "inferred_chain": label.chain,
        "inferred_protein_type": _protein_type(label.mhc_class, label.chain),
        "inferred_gene": inferred_gene,
        "inferred_gene_status": classified["gene_status"],
        "inferred_label_status": label.label_status,
        "inferred_disposition": label.disposition,
        "_parse_mhc_class": label.mhc_class if label.mhc_class in {"I", "II"} else "",
        "_parse_chain": label.chain if label.chain in {"alpha", "beta"} else "",
        "_parse_use_early_shortcuts": str(use_early_shortcuts).lower(),
    }
    if len(row["sequence"]) != int(row["source_length"] or 0):
        raise ValueError(f"Artifact sequence length mismatch for {accession}")
    return row


def build_records(
    artifacts: Sequence[dict[str, str]],
    *,
    source_release: str,
    version: str,
    taxonomy: dict[str, dict[str, str]],
    curation: dict[str, dict[str, str]],
    workers: int = 1,
    use_early_shortcuts: bool = False,
) -> list[dict[str, str]]:
    """Build one source-complete normalized and parsed record per artifact."""
    prepared = [
        prepare_record(
            artifact,
            source_release=source_release,
            version=version,
            taxonomy=taxonomy,
            curation=curation,
            use_early_shortcuts=use_early_shortcuts,
        )
        for artifact in artifacts
    ]
    if workers <= 1:
        rows = [_parse_prepared_record(row) for row in prepared]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_parse_prepared_record, prepared, chunksize=64))
    rows.sort(key=lambda row: row["accession"])
    accessions = [row["accession"] for row in rows]
    if len(set(accessions)) != len(accessions):
        raise ValueError("MHC protein dataset contains duplicate accessions")
    return rows


def _write_gzip_csv(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write deterministic gzip CSV bytes, publishing the completed file atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}-", suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                    writer = csv.DictWriter(text, fieldnames=PROTEIN_FIELDS, lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(rows)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}-", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def write_dataset_manifest(
    path: Path,
    *,
    output: Path,
    rows: Sequence[dict[str, str]],
    source_manifest: dict[str, object],
    curation_path: Path,
    version: str,
    use_early_shortcuts: bool,
) -> dict[str, object]:
    """Write complete source, schema, model, and content provenance."""
    payload: dict[str, object] = {
        "dataset": DATASET_NAME,
        "dataset_version": version,
        "schema_version": DATASET_SCHEMA_VERSION,
        "scope": (
            "Every current or historically preserved record from the release-pinned broad UniProt Vertebrata MHC candidate query. "
            "inferred_disposition identifies query contaminants and unresolved candidates; parse_ok is structural parser status, "
            "not an identity label."
        ),
        "annotation_families": {
            "source": "Verbatim or losslessly normalized UniProt/UniSave/taxonomy facts; missing annotations remain empty.",
            "inferred": "MHC identity and label-policy fields inferred from source metadata and accession curation.",
            "parsed": "mhcseqs sequence-parser outputs; partial parts may be populated when parse_ok is false.",
        },
        "source": {
            "database": "UniProtKB and UniSave",
            "release": source_manifest["release"],
            "bundle_schema_version": source_manifest["schema_version"],
            "artifacts": source_manifest["artifacts"],
            "provenance": source_manifest["provenance"],
        },
        "license": {
            "identifier": "CC-BY-4.0",
            "attribution": "UniProt Consortium",
            "url": "https://www.uniprot.org/help/license/",
        },
        "generator": {
            "name": "mhcseqs",
            "version": __version__,
            "curation_file": curation_path.name,
            "curation_sha256": _sha256(curation_path),
            "sp_boundary_model_sha256": _sha256(SP_BOUNDARY_MODEL_PATH),
            "sp_sequence_cue_model_sha256": _sha256(SP_SEQUENCE_CUE_MODEL_PATH),
            "use_early_shortcuts": use_early_shortcuts,
        },
        "coordinates": "zero-based, half-open for normalized numeric columns; domain_spans is one-based, inclusive",
        "records": {
            "filename": output.name,
            "rows": len(rows),
            "bytes": output.stat().st_size,
            "sha256": _sha256(output),
            "columns": list(PROTEIN_FIELDS),
            "source_kinds": dict(sorted(Counter(row["source_kind"] for row in rows).items())),
            "source_signal_status": dict(sorted(Counter(row["source_signal_status"] for row in rows).items())),
            "inferred_disposition": dict(sorted(Counter(row["inferred_disposition"] for row in rows).items())),
            "inferred_protein_type": dict(sorted(Counter(row["inferred_protein_type"] or "unresolved" for row in rows).items())),
            "parse_status": dict(sorted(Counter(row["parse_status"] for row in rows).items())),
            "parsed_protein_type": dict(sorted(Counter(row["parsed_protein_type"] or "unresolved" for row in rows).items())),
        },
    }
    _write_json(path, payload)
    return payload


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, help="Cache root containing the validated UniProt source bundle.")
    parser.add_argument("--output", type=Path, help="Output .csv.gz path (default: versioned directory under the cache root).")
    parser.add_argument("--manifest-output", type=Path, help="Output manifest path (default: beside the dataset).")
    parser.add_argument("--label-curation", type=Path, default=GT_LABEL_CURATION_CSV)
    parser.add_argument("--revision", type=int, default=DATASET_REVISION)
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument(
        "--use-early-shortcuts",
        action="store_true",
        help="Opt into production exact/lexical SP shortcuts; disabled by default to keep parsed annotations independent of source boundaries.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = _parse_args(argv)
    if not args.label_curation.is_file():
        raise FileNotFoundError(f"Missing pinned MHC label curation: {args.label_curation}")
    source_manifest = validate_artifact_bundle(args.data_dir)
    version = dataset_version(str(source_manifest["release"]), args.revision)
    cache_root = args.data_dir or Path(os.environ.get("MHCSEQS_DATA", Path.home() / ".cache" / "mhcseqs"))
    output = args.output or cache_root / "datasets" / DATASET_NAME / version / DATASET_FILENAME.format(version=version)
    manifest_output = args.manifest_output or output.with_name(MANIFEST_FILENAME.format(version=version))

    artifacts = load_corpus_artifact(args.data_dir) + load_deleted_artifact(args.data_dir)
    taxonomy = taxonomy_by_id(args.data_dir)
    curation = load_label_curation(args.label_curation)
    use_early_shortcuts = args.use_early_shortcuts
    print(f"Parsing {len(artifacts):,} source records with {args.workers} worker(s)...", flush=True)
    rows = build_records(
        artifacts,
        source_release=str(source_manifest["release"]),
        version=version,
        taxonomy=taxonomy,
        curation=curation,
        workers=args.workers,
        use_early_shortcuts=use_early_shortcuts,
    )
    _write_gzip_csv(output, rows)
    manifest = write_dataset_manifest(
        manifest_output,
        output=output,
        rows=rows,
        source_manifest=read_manifest(args.data_dir),
        curation_path=args.label_curation,
        version=version,
        use_early_shortcuts=use_early_shortcuts,
    )

    print(f"Wrote {len(rows):,} MHC protein records to {output}")
    print(f"Wrote data manifest to {manifest_output}")
    print(f"Data version: {version}")
    print(f"SHA-256: {manifest['records']['sha256']}")


if __name__ == "__main__":
    main()
