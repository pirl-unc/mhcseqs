#!/usr/bin/env python3
"""Enrich the SP ground-truth corpus with MHC metadata and control sequences.

The raw SP corpus only carries UniProt accessions, organisms, sequences, and
annotated signal-peptide lengths. This script upgrades it into a benchmark
that can be dispatched by source-backed MHC class and chain while preserving
the release-pinned metadata provenance.

Outputs:
  - data/sp_ground_truth_enriched.csv
  - data/sp_negative_controls.csv

The enriched CSV adds:
  - source_clade / lineage-backed species_category
  - mhc_class / chain / gene
  - protein_name / gene_names
  - is_fragment / source_group
  - metadata_source / gene_status / raw_gene_label / label_status

The default control CSV contains mature-only sequences derived from the GT
rows by removing the annotated leader.  Optional fragment-like controls are
available behind ``--include-fragment-controls``, but they are slower because
they require a full domain decomposition of each parent sequence.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
if str(SCRIPTS_DIR.parent) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR.parent))

from curate_diverse_mhc import (
    classify_mhc,
    derive_species_alias,
    resolve_gene_annotation,
)
from evaluate_sp_ground_truth import GT_RAW_CSV

from scripts.sp_corpus_artifacts import (
    load_corpus_artifact,
    load_deleted_artifact,
    taxonomy_by_id,
    validate_artifact_bundle,
)
from scripts.sp_ground_truth_eligibility import (
    GT_LABEL_CURATION_CSV,
    load_label_curation,
    resolve_mhc_label,
)
from scripts.sp_ground_truth_taxonomy import classification_from_taxonomy_row

ROOT = Path(__file__).resolve().parent.parent
GT_ENRICHED_CSV = ROOT / "data" / "sp_ground_truth_enriched.csv"
NEGATIVE_CONTROL_CSV = ROOT / "data" / "sp_negative_controls.csv"
ENRICHED_FIELDS = [
    "accession",
    "organism",
    "taxon_id",
    "source_clade",
    "species_category",
    "sp_length",
    "reviewed",
    "mhc_class",
    "chain",
    "gene",
    "protein_name",
    "gene_names",
    "is_fragment",
    "source_group",
    "metadata_source",
    "label_status",
    "gene_status",
    "raw_gene_label",
    "sequence",
]

CONTROL_FIELDS = [
    "control_id",
    "source_accession",
    "control_type",
    "organism",
    "taxon_id",
    "source_clade",
    "species_category",
    "reviewed",
    "mhc_class",
    "chain",
    "gene",
    "expected_sp_length",
    "parent_sp_length",
    "metadata_source",
    "sequence",
]


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _load_label_curation(path: Path = GT_LABEL_CURATION_CSV) -> dict[str, dict[str, str]]:
    """Load accession-level benchmark label decisions and provenance."""
    return load_label_curation(path)


def _apply_label_curation(
    row: dict[str, str],
    curation: dict[str, dict[str, str]],
) -> dict[str, str]:
    """Apply a source-backed class/chain decision to one enriched row."""
    decision = curation.get(row["accession"])
    if decision is None:
        return row

    disposition = decision.get("disposition", "")
    if disposition == "include":
        mhc_class = decision.get("mhc_class", "")
        chain = decision.get("chain", "")
        if mhc_class not in {"I", "II"} or chain not in {"alpha", "beta"}:
            raise ValueError(f"Invalid curated MHC label for {row['accession']}: {(mhc_class, chain)!r}")
        row["mhc_class"] = mhc_class
        row["chain"] = chain
    elif disposition == "exclude_non_mhc":
        row["mhc_class"] = ""
        row["chain"] = ""
    elif disposition == "retain_unresolved":
        row["mhc_class"] = decision.get("mhc_class", "")
        row["chain"] = decision.get("chain", "")
    else:
        raise ValueError(f"Invalid curation disposition for {row['accession']}: {disposition!r}")

    row["label_status"] = decision.get("label_status", "")
    return row


def _load_artifact_metadata(data_dir: Path | None = None) -> dict[str, dict[str, str]]:
    """Load current and archived UniProt metadata from one validated bundle."""
    validate_artifact_bundle(data_dir)
    result: dict[str, dict[str, str]] = {}
    for rows, source in (
        (load_corpus_artifact(data_dir), "uniprot_corpus_artifact"),
        (load_deleted_artifact(data_dir), "unisave_artifact"),
    ):
        for row in rows:
            accession = row["Entry"]
            if accession in result:
                raise ValueError(f"Duplicate artifact metadata for {accession}")
            result[accession] = {**row, "_metadata_source": source}
    return result


def _gene_from_protein_name(protein_name: str, organism: str, mhc_class: str, chain: str) -> tuple[str, str]:
    """Recover a gene-like token from UniProt protein names when gene_names are missing.

    This is intentionally conservative: it only translates explicit tokens in
    UniProt names such as ``F10 alpha chain`` or ``DR beta 1 chain`` into the
    gene namespace used by the rest of the package.
    """
    prefix = derive_species_alias(organism)
    text = protein_name or ""

    m = re.search(r"\(MHC class [^)]* antigen ([A-Z0-9-]+)\)", text, re.IGNORECASE)
    if m:
        raw = m.group(1).upper()
        bare = raw.replace("-", "")
        if bare.startswith("HLA"):
            bare = bare[3:]
        if bare.startswith(("DRA", "DRB", "DQA", "DQB", "DPA", "DPB", "DMA", "DMB", "DOA", "DOB", "DXA", "DXB")):
            return (f"{prefix}-{bare}" if prefix else bare, "protein_name_parenthetical")
        return (f"{prefix}-{bare}" if prefix else bare, "protein_name_parenthetical")

    m = re.search(
        r"antigen,\s*([A-Za-z0-9-]+)\s*(alpha|beta)\s*(\d+)?\s*chain(?:-like)?",
        text,
        re.IGNORECASE,
    )
    if not m:
        return "", "protein_name_unresolved"

    token = m.group(1).upper()
    arm = m.group(2).lower()
    suffix = m.group(3) or ""

    if token in {"B-L", "BL", "BLA", "BLB"}:
        bare = "BLA" if arm == "alpha" else "BLB"
    elif token in {"DR", "DQ", "DP", "DM", "DO", "DX", "DY", "DN"}:
        bare = f"{token}{'A' if arm == 'alpha' else 'B'}{suffix}"
    elif token == "HLA":
        return "", "protein_name_unresolved"
    else:
        bare = f"{token}{suffix}"

    if prefix:
        return f"{prefix}-{bare}", "protein_name_token"
    return bare, "protein_name_token"


def _classify_from_names(
    *,
    organism: str,
    protein_name: str,
    gene_names: str,
) -> dict[str, str]:
    cls = classify_mhc(protein_name, gene_names) or ("", "")
    mhc_class, chain = cls
    prefix = derive_species_alias(organism)
    gene, raw_gene_label, gene_status = resolve_gene_annotation(gene_names, protein_name, prefix)

    if not gene and mhc_class in {"I", "II"}:
        gene, fallback_status = _gene_from_protein_name(protein_name, organism, mhc_class, chain)
        if gene:
            raw_gene_label = raw_gene_label or gene
            gene_status = fallback_status

    label_status = "gold" if mhc_class else "unresolved"
    if mhc_class and chain in {"alpha", "beta", "B2M", "unknown"}:
        label_status = "gold"

    return {
        "mhc_class": mhc_class,
        "chain": chain,
        "gene": gene,
        "raw_gene_label": raw_gene_label,
        "gene_status": gene_status,
        "label_status": label_status,
    }


def _merge_row(
    row: dict[str, str],
    *,
    artifacts: dict[str, dict[str, str]],
    taxonomy: dict[str, dict[str, str]],
    label_curation: dict[str, dict[str, str]],
) -> dict[str, str]:
    accession = row["accession"]
    organism = row["organism"]
    taxon_id = row.get("taxon_id", "")
    source_clade = row.get("source_clade", "")
    taxonomy_row = taxonomy.get(taxon_id)
    if taxonomy_row is None:
        raise ValueError(f"Missing artifact taxonomy for {accession} taxon {taxon_id}")
    artifact_clade, _source_root, species_category, _category_root = classification_from_taxonomy_row(taxonomy_row)
    if source_clade != artifact_clade:
        raise ValueError(f"Source clade mismatch for {accession}: raw={source_clade!r}, artifact={artifact_clade!r}")

    artifact = artifacts.get(accession)
    if artifact is None:
        raise ValueError(f"Missing artifact metadata for ground-truth accession {accession}")
    protein_name = artifact.get("Protein names", "").strip()
    gene_names = artifact.get("Gene Names", "").strip()
    label = resolve_mhc_label(artifact, label_curation)
    if not label.eligible:
        raise ValueError(f"Ineligible accession {accession} reached enriched ground truth")
    classified = _classify_from_names(
        organism=organism,
        protein_name=protein_name,
        gene_names=gene_names,
    )

    merged = {
        "accession": accession,
        "organism": organism,
        "taxon_id": taxon_id,
        "source_clade": source_clade,
        "species_category": species_category,
        "sp_length": row["sp_length"],
        "reviewed": row["reviewed"],
        "sequence": row["sequence"],
        "mhc_class": label.mhc_class,
        "chain": label.chain,
        "gene": classified["gene"],
        "protein_name": protein_name,
        "gene_names": gene_names,
        "is_fragment": artifact.get("Fragment", "").strip(),
        "source_group": "",
        "metadata_source": artifact["_metadata_source"],
        "label_status": label.label_status,
        "gene_status": classified["gene_status"],
        "raw_gene_label": classified["raw_gene_label"],
    }

    return _apply_label_curation(merged, label_curation)


def _write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _build_controls(rows: list[dict[str, str]], *, include_fragment_controls: bool = False) -> list[dict[str, str]]:
    controls: list[dict[str, str]] = []
    for row in rows:
        mhc_class = row.get("mhc_class", "")
        chain = row.get("chain", "")
        label_status = row.get("label_status", "")
        if label_status not in {"gold", "curated"} or mhc_class not in {"I", "II"} or chain not in {"alpha", "beta"}:
            continue

        sequence = row["sequence"]
        sp_length = int(row["sp_length"])
        if sp_length <= 0 or sp_length >= len(sequence):
            continue

        base = {
            "source_accession": row["accession"],
            "organism": row["organism"],
            "taxon_id": row["taxon_id"],
            "source_clade": row["source_clade"],
            "species_category": row["species_category"],
            "reviewed": row["reviewed"],
            "mhc_class": mhc_class,
            "chain": chain,
            "gene": row.get("gene", ""),
            "expected_sp_length": "0",
            "parent_sp_length": row["sp_length"],
            "metadata_source": row.get("metadata_source", ""),
        }

        mature_only = sequence[sp_length:]
        controls.append(
            {
                **base,
                "control_id": f"{row['accession']}:mature_only",
                "control_type": "mature_only",
                "sequence": mature_only,
            }
        )

        if not include_fragment_controls:
            continue

        from mhcseqs.domain_parsing import decompose_domains

        try:
            parsed = decompose_domains(
                sequence,
                mhc_class=mhc_class,
                chain=chain or None,
                gene=row.get("gene", ""),
                species=row.get("organism", ""),
            )
        except Exception:
            continue
        if not parsed.ok:
            continue
        if mhc_class == "I" and parsed.groove_seq:
            controls.append(
                {
                    **base,
                    "control_id": f"{row['accession']}:class_i_exon23_like",
                    "control_type": "class_i_exon23_like",
                    "sequence": parsed.groove_seq,
                }
            )
        elif mhc_class == "II" and chain == "alpha" and parsed.groove1:
            controls.append(
                {
                    **base,
                    "control_id": f"{row['accession']}:class_ii_exon2_like",
                    "control_type": "class_ii_exon2_like",
                    "sequence": parsed.groove1,
                }
            )
        elif mhc_class == "II" and chain == "beta" and parsed.groove2:
            controls.append(
                {
                    **base,
                    "control_id": f"{row['accession']}:class_ii_exon2_like",
                    "control_type": "class_ii_exon2_like",
                    "sequence": parsed.groove2,
                }
            )

    return controls


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--include-fragment-controls",
        action="store_true",
        help="Also generate synthetic exon-like groove fragments (slower).",
    )
    parser.add_argument("--data-dir", type=Path, help="Artifact cache root.")
    parser.add_argument("--raw", type=Path, default=GT_RAW_CSV)
    parser.add_argument("--output", type=Path, default=GT_ENRICHED_CSV)
    parser.add_argument("--controls", type=Path, default=NEGATIVE_CONTROL_CSV)
    parser.add_argument("--label-curation", type=Path, default=GT_LABEL_CURATION_CSV)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)

    gt_rows = _read_csv(args.raw)
    label_curation = _load_label_curation(args.label_curation)
    artifacts = _load_artifact_metadata(args.data_dir)
    taxonomy = taxonomy_by_id(args.data_dir)

    enriched = [
        _merge_row(
            row,
            artifacts=artifacts,
            taxonomy=taxonomy,
            label_curation=label_curation,
        )
        for row in gt_rows
    ]
    controls = _build_controls(enriched, include_fragment_controls=args.include_fragment_controls)

    _write_csv(args.output, enriched, ENRICHED_FIELDS)
    _write_csv(args.controls, controls, CONTROL_FIELDS)

    by_source = Counter(row["metadata_source"] for row in enriched)
    by_label = Counter((row["mhc_class"], row["chain"]) for row in enriched if row["mhc_class"])
    by_gene_status = Counter(row["gene_status"] for row in enriched)
    by_control = Counter(row["control_type"] for row in controls)

    print(f"Loaded {len(gt_rows)} raw GT rows from {args.raw.name}")
    print(f"Loaded release-pinned artifact metadata for {len(artifacts)} accessions")
    print(f"Wrote enriched GT to {args.output}")
    print(f"Wrote negative controls to {args.controls}")
    print("\nMetadata sources:")
    for source, count in sorted(by_source.items()):
        print(f"  {source}: {count}")
    print("\nGold class / chain labels:")
    for key, count in sorted(by_label.items()):
        print(f"  {key[0]} {key[1]}: {count}")
    print("\nGene status:")
    for status, count in sorted(by_gene_status.items()):
        print(f"  {status}: {count}")
    print("\nControls:")
    for kind, count in sorted(by_control.items()):
        print(f"  {kind}: {count}")


if __name__ == "__main__":
    main()
