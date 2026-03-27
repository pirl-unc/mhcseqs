#!/usr/bin/env python3
"""Enrich the SP ground-truth corpus with MHC metadata and control sequences.

The raw SP corpus only carries UniProt accessions, organisms, sequences, and
annotated signal-peptide lengths.  This script upgrades it into a benchmark
that can be dispatched by gold MHC class / chain whenever UniProt provides
enough metadata, while still preserving provenance for ambiguous rows.

Outputs:
  - data/sp_ground_truth_enriched.csv
  - data/sp_negative_controls.csv

The enriched CSV adds:
  - species_category
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
import io
import re
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path
from urllib.error import HTTPError

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from curate_diverse_mhc import (
    classify_mhc,
    derive_prefix,
    resolve_gene_annotation,
)
from evaluate_sp_ground_truth import GT_RAW_CSV, _species_category

ROOT = Path(__file__).resolve().parent.parent
GT_ENRICHED_CSV = ROOT / "data" / "sp_ground_truth_enriched.csv"
NEGATIVE_CONTROL_CSV = ROOT / "data" / "sp_negative_controls.csv"
DIVERSE_CURATED_CSV = ROOT / "mhcseqs" / "diverse_mhc_sequences.csv"
MOUSE_H2_CSV = ROOT / "mhcseqs" / "mouse_h2_sequences.csv"
DIVERSE_RAW_CSV = ROOT / "data" / "diverse_mhc_raw.csv"

ENRICHED_FIELDS = [
    "accession",
    "organism",
    "taxon_id",
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


def _load_local_metadata() -> tuple[dict[str, dict], dict[str, dict], dict[str, dict]]:
    diverse_rows = _read_csv(DIVERSE_CURATED_CSV)
    mouse_rows = _read_csv(MOUSE_H2_CSV)
    raw_rows = _read_csv(DIVERSE_RAW_CSV)

    curated = {
        row["uniprot_accession"]: row
        for row in diverse_rows
    }
    mouse = {
        row["uniprot_accession"]: row
        for row in mouse_rows
    }
    raw = {
        row["uniprot_accession"]: row
        for row in raw_rows
    }
    return curated, mouse, raw


def _chunked(values: list[str], size: int) -> list[list[str]]:
    return [values[i : i + size] for i in range(0, len(values), size)]


def _fetch_uniprot_chunk(chunk: list[str]) -> dict[str, dict[str, str]]:
    query = "(" + " OR ".join(f"accession:{acc}" for acc in chunk) + ")"
    params = urllib.parse.urlencode(
        {
            "query": query,
            "format": "tsv",
            "fields": "accession,protein_name,gene_names,organism_name,organism_id,reviewed,fragment",
            "size": str(len(chunk)),
        }
    )
    url = f"https://rest.uniprot.org/uniprotkb/search?{params}"
    req = urllib.request.Request(url, headers={"User-Agent": "mhcseqs-sp-gt/1.0"})
    with urllib.request.urlopen(req, timeout=60) as response:
        payload = response.read().decode("utf-8")
    reader = csv.DictReader(io.StringIO(payload), delimiter="\t")
    rows: dict[str, dict[str, str]] = {}
    for row in reader:
        accession = row.get("Entry", "").strip()
        if accession:
            rows[accession] = {
                "protein_name": row.get("Protein names", "").strip(),
                "gene_names": row.get("Gene Names", "").strip(),
                "organism": row.get("Organism", "").strip(),
                "organism_id": row.get("Organism (ID)", "").strip(),
                "fragment": row.get("Fragment", "").strip(),
                "reviewed": row.get("Reviewed", "").strip(),
            }
    return rows


def _fetch_uniprot_metadata(accessions: list[str], chunk_size: int = 64) -> dict[str, dict[str, str]]:
    """Fetch protein/gene metadata for accessions via UniProt search."""
    results: dict[str, dict[str, str]] = {}
    for chunk in _chunked(accessions, chunk_size):
        try:
            results.update(_fetch_uniprot_chunk(chunk))
            continue
        except HTTPError as exc:
            if exc.code != 400 or len(chunk) <= 8:
                raise

        midpoint = len(chunk) // 2
        results.update(_fetch_uniprot_metadata(chunk[:midpoint], chunk_size=max(midpoint, 8)))
        results.update(_fetch_uniprot_metadata(chunk[midpoint:], chunk_size=max(len(chunk) - midpoint, 8)))
    return results


def _gene_from_protein_name(protein_name: str, organism: str, mhc_class: str, chain: str) -> tuple[str, str]:
    """Recover a gene-like token from UniProt protein names when gene_names are missing.

    This is intentionally conservative: it only translates explicit tokens in
    UniProt names such as ``F10 alpha chain`` or ``DR beta 1 chain`` into the
    gene namespace used by the rest of the package.
    """
    prefix = derive_prefix(organism)
    text = protein_name or ""

    m = re.search(r"\(MHC class [^)]* antigen ([A-Z0-9-]+)\)", text, re.IGNORECASE)
    if m:
        raw = m.group(1).upper()
        bare = raw.replace("-", "")
        if bare.startswith("HLA"):
            if prefix.upper() == "HOSA":
                return (f"HLA-{bare[3:]}", "protein_name_parenthetical")
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
    prefix = derive_prefix(organism)
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
    curated: dict[str, dict],
    mouse: dict[str, dict],
    raw: dict[str, dict],
    fetched: dict[str, dict[str, str]],
) -> dict[str, str]:
    accession = row["accession"]
    organism = row["organism"]
    taxon_id = row.get("taxon_id", "")
    species_category = _species_category(organism, taxon_id)

    local_curated = curated.get(accession)
    local_mouse = mouse.get(accession)
    local_raw = raw.get(accession)
    remote = fetched.get(accession, {})

    protein_name = ""
    gene_names = ""
    if local_raw:
        protein_name = local_raw.get("protein_name", "")
        gene_names = local_raw.get("gene_names", "")
    if remote:
        protein_name = remote.get("protein_name", "") or protein_name
        gene_names = remote.get("gene_names", "") or gene_names

    merged = {
        "accession": accession,
        "organism": organism,
        "taxon_id": taxon_id,
        "species_category": species_category,
        "sp_length": row["sp_length"],
        "reviewed": row["reviewed"],
        "sequence": row["sequence"],
        "mhc_class": "",
        "chain": "",
        "gene": "",
        "protein_name": protein_name,
        "gene_names": gene_names,
        "is_fragment": remote.get("fragment", ""),
        "source_group": "",
        "metadata_source": "",
        "label_status": "unresolved",
        "gene_status": "missing",
        "raw_gene_label": "",
    }

    if local_curated:
        merged.update(
            {
                "mhc_class": local_curated.get("mhc_class", ""),
                "chain": local_curated.get("chain", ""),
                "gene": local_curated.get("gene", ""),
                "is_fragment": local_curated.get("is_fragment", ""),
                "source_group": local_curated.get("source_group", ""),
                "metadata_source": "local_curated_accession",
                "label_status": "gold",
                "gene_status": "ok" if local_curated.get("gene") else "missing",
                "raw_gene_label": local_curated.get("gene", ""),
            }
        )
        return merged

    if local_mouse:
        merged.update(
            {
                "mhc_class": local_mouse.get("mhc_class", ""),
                "chain": local_mouse.get("chain", ""),
                "gene": local_mouse.get("gene", ""),
                "is_fragment": local_mouse.get("is_fragment", ""),
                "source_group": "murine_curated",
                "metadata_source": "local_mouse_accession",
                "label_status": "gold",
                "gene_status": "ok" if local_mouse.get("gene") else "missing",
                "raw_gene_label": local_mouse.get("gene", ""),
            }
        )
        return merged

    if local_raw:
        merged["is_fragment"] = local_raw.get("is_fragment", "")
        merged["source_group"] = local_raw.get("source_group", "")

    classified = _classify_from_names(
        organism=organism,
        protein_name=protein_name,
        gene_names=gene_names,
    )
    merged.update(classified)
    if remote:
        merged["metadata_source"] = "uniprot_accession_tsv"
    elif local_raw:
        merged["metadata_source"] = "local_raw_accession"
    else:
        merged["metadata_source"] = "none"

    return merged


def _write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _build_controls(rows: list[dict[str, str]], *, include_fragment_controls: bool = False) -> list[dict[str, str]]:
    controls: list[dict[str, str]] = []
    for row in rows:
        mhc_class = row.get("mhc_class", "")
        chain = row.get("chain", "")
        if mhc_class not in {"I", "II"}:
            continue

        sequence = row["sequence"]
        sp_length = int(row["sp_length"])
        if sp_length <= 0 or sp_length >= len(sequence):
            continue

        base = {
            "source_accession": row["accession"],
            "organism": row["organism"],
            "taxon_id": row["taxon_id"],
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--include-fragment-controls",
        action="store_true",
        help="Also generate synthetic exon-like groove fragments (slower).",
    )
    args = parser.parse_args()

    gt_rows = _read_csv(GT_RAW_CSV)
    accessions = [row["accession"] for row in gt_rows]
    curated, mouse, raw = _load_local_metadata()
    fetched = _fetch_uniprot_metadata(accessions)

    enriched = [
        _merge_row(row, curated=curated, mouse=mouse, raw=raw, fetched=fetched)
        for row in gt_rows
    ]
    controls = _build_controls(enriched, include_fragment_controls=args.include_fragment_controls)

    _write_csv(GT_ENRICHED_CSV, enriched, ENRICHED_FIELDS)
    _write_csv(NEGATIVE_CONTROL_CSV, controls, CONTROL_FIELDS)

    by_source = Counter(row["metadata_source"] for row in enriched)
    by_label = Counter((row["mhc_class"], row["chain"]) for row in enriched if row["mhc_class"])
    by_gene_status = Counter(row["gene_status"] for row in enriched)
    by_control = Counter(row["control_type"] for row in controls)

    print(f"Loaded {len(gt_rows)} raw GT rows from {GT_RAW_CSV.name}")
    print(f"Fetched UniProt metadata for {len(fetched)} accessions")
    print(f"Wrote enriched GT to {GT_ENRICHED_CSV}")
    print(f"Wrote negative controls to {NEGATIVE_CONTROL_CSV}")
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
