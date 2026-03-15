"""Main pipeline: parse FASTA → raw CSV → full-seqs CSV.

This module contains the core logic for building both output files.
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from .alleles import (
    allele_suffix_flags,
    infer_gene,
    infer_mhc_class,
    infer_species_identity,
    normalize_allele_name,
    normalize_mhc_class,
    parse_allele_name,
)
from .groove import (
    NON_GROOVE_GENES,
    AlleleRecord,
    is_class_ii_alpha_gene,
    parse_class_i,
    parse_class_ii_alpha,
    parse_class_ii_beta,
)
from .species import get_latin_name, normalize_mhc_species

# Minimum protein length to include (allows groove-bearing fragments)
MIN_MHC_SEQUENCE_LEN = 70

# Curated reference CSVs (shipped with this repo)
_B2M_CSV = Path(__file__).resolve().parent / "b2m_sequences.csv"
_MOUSE_H2_CSV = Path(__file__).resolve().parent / "mouse_h2_sequences.csv"
_DIVERSE_MHC_CSV = Path(__file__).resolve().parent / "diverse_mhc_sequences.csv"

_NUCLEOTIDE_LIKE_CHARS = set("ACGTUNWSMKRYBDHV")

# Statuses considered functional for groove extraction
FUNCTIONAL_GROOVE_STATUSES = {
    "ok",
    "alpha3_fallback",
    "beta1_only_fallback",
    "fragment_fallback",
}

# ---------------------------------------------------------------------------
# CSV field definitions
# ---------------------------------------------------------------------------

RAW_FIELDS = [
    "allele_raw",
    "allele_normalized",
    "two_field_allele",
    "gene",
    "mhc_class",
    "chain",
    "species",
    "species_category",
    "species_prefix",
    "source",
    "source_id",
    "seq_len",
    "sequence",
    "has_signal_peptide",
    "signal_peptide_len",
    "signal_peptide_seq",
    "is_null",
    "is_questionable",
    "is_pseudogene",
]

FULL_FIELDS = [
    "two_field_allele",
    "representative_allele",
    "protein_seq_selection",
    "gene",
    "mhc_class",
    "chain",
    "species",
    "species_category",
    "species_prefix",
    "source",
    "source_id",
    "seq_len",
    "sequence",
    "mature_start",
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
    "groove_status",
    "groove_flags",
    "anchor_type",
    "is_null",
    "is_questionable",
    "is_pseudogene",
    "is_functional",
]


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------


def _iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) from a FASTA file."""
    header = None
    seq_parts: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def _looks_like_nucleotide(sequence: str) -> bool:
    seq = re.sub(r"\s+", "", str(sequence or "").strip().upper())
    if not seq:
        return False
    chars = {ch for ch in seq if ch.isalpha()}
    nucleotide_chars = chars & set("ACGTU")
    return bool(chars) and chars <= _NUCLEOTIDE_LIKE_CHARS and len(nucleotide_chars) >= 3


def _extract_source_id(header: str, source_label: str) -> str:
    """Extract database accession from a FASTA header.

    IMGT/HLA headers start with ``HLA:HLA00001``.
    IPD-MHC headers start with ``IPD-MHC:NHP00001``.
    Returns empty string if no accession is found.
    """
    if source_label == "imgt":
        # Format: "HLA:HLA00001 A*01:01:01:01 365 bp"
        m = re.match(r"HLA:(HLA\d+)", header)
        return m.group(1) if m else ""
    if source_label == "ipd_mhc":
        # Format: "IPD-MHC:NHP00001 Aona-DQA1*27:01 73 bp"
        m = re.match(r"IPD-MHC:(\w+)", header)
        return m.group(1) if m else ""
    return ""


def _candidate_tokens(header: str) -> List[str]:
    """Score and rank tokens from a FASTA header for allele parsing."""
    cleaned = header.replace("|", " ").replace(";", " ")
    tokens = [t.strip() for t in cleaned.split() if t.strip()]
    if not tokens:
        return []
    scored = []
    for tok in tokens:
        score = 0
        if "*" in tok:
            score += 3
        if tok.upper().startswith("HLA-"):
            score += 2
        if tok.upper().startswith("H-2"):
            score += 2
        if tok.upper().startswith("MAMU"):
            score += 2
        if ":" in tok:
            score += 1
        scored.append((score, tok))
    scored.sort(key=lambda x: x[0], reverse=True)
    return [tok for _, tok in scored]


def _resolve_header_allele(header: str):
    """Parse allele metadata from a FASTA header using mhcgnomes."""
    last_error = None
    for token in _candidate_tokens(header):
        try:
            parsed = parse_allele_name(token)
            if parsed is None:
                continue
            normalized = parsed.to_string()
            gene = getattr(parsed.gene, "name", None) if getattr(parsed, "gene", None) else None
            species = getattr(parsed.species, "name", None) if getattr(parsed, "species", None) else None
            mhc_class = normalize_mhc_class(getattr(parsed, "mhc_class", None))
            return normalized, gene, mhc_class, species, token
        except Exception as exc:
            last_error = exc
            continue
    raise RuntimeError(f"Failed to parse allele from FASTA header: '{header}'.") from last_error


def _infer_species_prefix(allele: Optional[str]) -> str:
    """Extract the species MHC prefix from an allele name (e.g. HLA, SLA, Mamu)."""
    if not allele:
        return ""
    try:
        parsed = parse_allele_name(allele)
    except Exception:
        parsed = None
    if parsed is not None:
        species = getattr(parsed, "species", None)
        prefix = getattr(species, "mhc_prefix", None)
        if prefix:
            return str(prefix)
    # Fallback: extract prefix before the dash
    token = str(allele).strip()
    if "-" in token:
        return token.split("-", 1)[0]
    return ""


def _load_b2m_references() -> List[dict]:
    """Load curated B2M reference sequences from data/b2m_sequences.csv."""
    from .species import get_canonical_prefix

    if not _B2M_CSV.exists():
        return []
    rows = []
    with open(_B2M_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            species_key = row["species_key"]
            seq = row["sequence"]
            accession = row.get("uniprot_accession", "")
            allele_name = f"B2M_{species_key}"
            latin = get_latin_name(species_key)
            category = normalize_mhc_species(species_key) or ""
            prefix = get_canonical_prefix(species_key)
            rows.append(
                {
                    "allele_raw": accession or allele_name,
                    "allele_normalized": allele_name,
                    "two_field_allele": allele_name,
                    "gene": "B2M",
                    "mhc_class": "I",
                    "chain": "B2M",
                    "species": latin,
                    "species_category": category,
                    "species_prefix": prefix,
                    "source": "uniprot_reference",
                    "source_id": accession,
                    "seq_len": str(len(seq)),
                    "sequence": seq,
                    "has_signal_peptide": "False",
                    "signal_peptide_len": "",
                    "signal_peptide_seq": "",
                    "is_null": "False",
                    "is_questionable": "False",
                    "is_pseudogene": "False",
                }
            )
    return rows


def _load_mouse_h2_references() -> List[dict]:
    """Load curated mouse H-2 sequences from mouse_h2_sequences.csv.

    These are SwissProt-reviewed reference sequences covering the major
    H-2 loci (K, D, L class I; Aa, Ab1, Ea, Eb1 class II) across
    common inbred mouse haplotypes (b, d, k, q, s, etc.).
    """
    if not _MOUSE_H2_CSV.exists():
        return []
    rows = []
    with open(_MOUSE_H2_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            allele_name = row["allele_name"]
            gene = row["gene"]
            mhc_class = row["mhc_class"]
            chain = row["chain"]
            accession = row.get("uniprot_accession", "")
            seq = row["sequence"]
            is_fragment = row.get("is_fragment", "False") == "True"

            # Groove parse to infer signal peptide / mature start
            groove = _try_groove_parse(seq, mhc_class=mhc_class, gene=gene, allele=allele_name)
            mature_start = groove.mature_start if groove and groove.ok and groove.mature_start > 0 else 0
            has_sp = mature_start >= 15 and seq[:1].upper() == "M"
            sp_seq = seq[:mature_start] if has_sp else ""

            rows.append(
                {
                    "allele_raw": accession or allele_name,
                    "allele_normalized": allele_name,
                    "two_field_allele": allele_name,
                    "gene": gene,
                    "mhc_class": mhc_class,
                    "chain": chain,
                    "species": "Mus musculus",
                    "species_category": "murine",
                    "species_prefix": "H2",
                    "source": "uniprot_curated",
                    "source_id": accession,
                    "seq_len": str(len(seq)),
                    "sequence": seq,
                    "has_signal_peptide": str(has_sp),
                    "signal_peptide_len": str(mature_start),
                    "signal_peptide_seq": sp_seq,
                    "is_null": "False",
                    "is_questionable": str(is_fragment),
                    "is_pseudogene": "False",
                }
            )
    return rows


def _infer_chain(gene: str, mhc_class: str) -> str:
    """Infer chain type from gene name and MHC class."""
    mc = normalize_mhc_class(mhc_class)
    if mc == "I":
        return "alpha"
    if mc == "II":
        if is_class_ii_alpha_gene(gene):
            return "alpha"
        return "beta"
    return ""


# Source group → species_category mapping (for diverse_mhc_sequences.csv)
_DIVERSE_GROUP_TO_CATEGORY = {
    "reptile_lepidosauria": "other_vertebrate",
    "reptile_crocodylia": "other_vertebrate",
    "reptile_testudines": "other_vertebrate",
    "amphibian": "other_vertebrate",
    "bird_non_chicken": "bird",
    "chicken": "bird",
    "shark_ray": "fish",
    "bony_fish": "fish",
    "marsupial": "other_mammal",
    "monotreme": "other_mammal",
    "bat": "other_mammal",
}


def _load_diverse_mhc_references() -> List[dict]:
    """Load curated diverse MHC sequences from diverse_mhc_sequences.csv.

    These are UniProt sequences from taxonomic groups underrepresented in
    IMGT/HLA and IPD-MHC: reptiles, amphibians, birds, fish, sharks,
    marsupials, monotremes, and bats.
    """
    if not _DIVERSE_MHC_CSV.exists():
        return []
    rows = []
    with open(_DIVERSE_MHC_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            accession = row.get("uniprot_accession", "")
            gene = row.get("gene", "")
            mhc_class = row.get("mhc_class", "")
            chain = row.get("chain", "")
            organism = row.get("organism", "")
            source_group = row.get("source_group", "")
            seq = row.get("sequence", "")
            is_fragment = row.get("is_fragment", "False") == "True"
            length = len(seq)

            if not seq or length < MIN_MHC_SEQUENCE_LEN:
                continue

            # Use accession as the allele name (no standard nomenclature)
            allele_name = accession

            # Determine species_category from source_group (reliable)
            species_category = _DIVERSE_GROUP_TO_CATEGORY.get(source_group, "")
            if not species_category:
                species_category = normalize_mhc_species(organism) or ""

            # Fix chain for class II unknowns: try class_ii_alpha_gene check
            if mhc_class == "II" and chain in ("unknown", ""):
                if gene and is_class_ii_alpha_gene(gene):
                    chain = "alpha"
                elif gene:
                    chain = "beta"

            # Groove parse for signal peptide detection
            groove = _try_groove_parse(seq, mhc_class=mhc_class, gene=gene, allele=allele_name)
            mature_start = groove.mature_start if groove and groove.ok and groove.mature_start > 0 else 0
            has_sp = mature_start >= 15 and seq[:1].upper() == "M"
            sp_seq = seq[:mature_start] if has_sp else ""

            rows.append(
                {
                    "allele_raw": accession,
                    "allele_normalized": allele_name,
                    "two_field_allele": allele_name,
                    "gene": gene,
                    "mhc_class": mhc_class,
                    "chain": chain,
                    "species": organism,
                    "species_category": species_category,
                    "species_prefix": "",
                    "source": "uniprot_diverse",
                    "source_id": accession,
                    "seq_len": str(length),
                    "sequence": seq,
                    "has_signal_peptide": str(has_sp),
                    "signal_peptide_len": str(mature_start),
                    "signal_peptide_seq": sp_seq,
                    "is_null": "False",
                    "is_questionable": str(is_fragment),
                    "is_pseudogene": "False",
                }
            )
    return rows


# ---------------------------------------------------------------------------
# Step 1: Raw index
# ---------------------------------------------------------------------------


def build_raw_index(
    fasta_paths: List[Tuple[Path, str]],
    out_csv: Path,
) -> Dict[str, int]:
    """Parse all FASTA files into the raw index CSV.

    Returns stats dict.
    """
    records: Dict[str, dict] = {}
    stats = Counter(
        total=0,
        parsed=0,
        skipped=0,
        skipped_nucleotide=0,
        skipped_short=0,
        duplicates=0,
        replaced=0,
    )

    for path, source_label in fasta_paths:
        for header, seq in _iter_fasta(path):
            stats["total"] += 1
            if _looks_like_nucleotide(seq):
                stats["skipped_nucleotide"] += 1
                continue
            # Allow even very short/empty for the raw CSV — but skip nucleotide
            source_id = _extract_source_id(header, source_label)
            try:
                normalized, gene, mhc_class, species_raw, allele_token = _resolve_header_allele(header)
            except RuntimeError:
                stats["skipped"] += 1
                continue

            if not gene:
                gene = infer_gene(normalized)
            if not mhc_class:
                mhc_class = infer_mhc_class(normalized) or ""
            mhc_class = normalize_mhc_class(mhc_class, default=infer_mhc_class(normalized)) or ""
            if not species_raw:
                species_raw = infer_species_identity(normalized) or ""

            chain = _infer_chain(gene or "", mhc_class)
            species_fine = species_raw
            species_category = normalize_mhc_species(species_fine) or ""
            species_prefix = _infer_species_prefix(normalized)

            # Groove parse to infer mature protein start.
            # signal_peptide_len stores the mature_start offset from the
            # groove parser (useful for alignment even when no true SP is
            # present).  has_signal_peptide is only True when the sequence
            # actually begins with Met and the offset is plausible.
            groove = _try_groove_parse(seq, mhc_class=mhc_class, gene=gene or "", allele=normalized)
            mature_start = groove.mature_start if groove and groove.ok and groove.mature_start > 0 else 0
            has_sp = mature_start >= 15 and seq[:1].upper() == "M"
            sp_seq = seq[:mature_start] if has_sp else ""

            suffix = allele_suffix_flags(normalized)

            try:
                two_field = normalize_allele_name(normalized)
            except Exception:
                two_field = normalized

            row = {
                "allele_raw": allele_token,
                "allele_normalized": normalized,
                "two_field_allele": two_field,
                "gene": gene or "",
                "mhc_class": mhc_class,
                "chain": chain,
                "species": species_fine,
                "species_category": species_category,
                "species_prefix": species_prefix,
                "source": source_label,
                "source_id": source_id,
                "seq_len": str(len(seq)),
                "sequence": seq,
                "has_signal_peptide": str(has_sp),
                "signal_peptide_len": str(mature_start),
                "signal_peptide_seq": sp_seq,
                "is_null": str(suffix["is_null"]),
                "is_questionable": str(suffix["is_questionable"]),
                "is_pseudogene": str(suffix["is_pseudogene"]),
            }

            # Dedup: prefer IMGT, prefer protein, prefer longer
            if normalized in records:
                existing = records[normalized]
                if _looks_like_nucleotide(existing["sequence"]):
                    records[normalized] = row
                    stats["replaced"] += 1
                elif existing["source"] != "imgt" and source_label == "imgt":
                    records[normalized] = row
                    stats["replaced"] += 1
                else:
                    stats["duplicates"] += 1
                continue

            records[normalized] = row
            stats["parsed"] += 1

    # Inject curated B2M reference sequences
    b2m_rows = _load_b2m_references()
    for b2m in b2m_rows:
        key = b2m["allele_normalized"]
        if key not in records:
            records[key] = b2m
            stats["parsed"] += 1
            stats["b2m_references"] = stats.get("b2m_references", 0) + 1

    # Inject curated mouse H-2 reference sequences
    h2_rows = _load_mouse_h2_references()
    for h2 in h2_rows:
        key = h2["allele_normalized"]
        if key not in records:
            records[key] = h2
            stats["parsed"] += 1
            stats["h2_references"] = stats.get("h2_references", 0) + 1

    # Inject diverse MHC sequences (reptiles, amphibians, birds, fish, etc.)
    diverse_rows = _load_diverse_mhc_references()
    for dr in diverse_rows:
        key = dr["allele_normalized"]
        if key not in records:
            records[key] = dr
            stats["parsed"] += 1
            stats["diverse_references"] = stats.get("diverse_references", 0) + 1

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=RAW_FIELDS)
        writer.writeheader()
        for key in sorted(records.keys()):
            writer.writerow(records[key])

    return dict(stats)


def _try_groove_parse(
    seq: str,
    *,
    mhc_class: str,
    gene: str,
    allele: str,
) -> Optional[AlleleRecord]:
    try:
        gene_upper = (gene or "").strip().upper()
        # B2M and non-groove genes (MICA/MICB/HFE) don't have peptide-binding grooves
        if gene_upper in ("B2M", "BETA-2-MICROGLOBULIN") or gene_upper in NON_GROOVE_GENES:
            return None
        nc = normalize_mhc_class(mhc_class)
        if nc == "I":
            return parse_class_i(seq, allele=allele, gene=gene)
        if nc == "II":
            if is_class_ii_alpha_gene(gene):
                return parse_class_ii_alpha(seq, allele=allele, gene=gene)
            return parse_class_ii_beta(seq, allele=allele, gene=gene)
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Step 2: Full sequences (two-field representative selection)
# ---------------------------------------------------------------------------


def _unique_by_sequence(members: List[dict]) -> List[dict]:
    """Deduplicate by sequence content, keeping the longest normalized name."""
    chosen: Dict[str, dict] = {}
    for m in members:
        seq = m["sequence"]
        cur = chosen.get(seq)
        if cur is None:
            chosen[seq] = m
        elif len(m["allele_normalized"]) > len(cur["allele_normalized"]):
            chosen[seq] = m
    return list(chosen.values())


def _merge_sequences_with_exact_overlap(
    left: str,
    right: str,
    *,
    min_overlap: int,
) -> Optional[str]:
    """Try to merge two sequences that share an exact overlap region."""
    left_seq = left.strip().upper()
    right_seq = right.strip().upper()
    if not left_seq or not right_seq:
        return None
    if right_seq in left_seq:
        return left_seq
    if left_seq in right_seq:
        return right_seq

    candidates: List[Tuple[int, int, str]] = []
    min_offset = -(len(right_seq) - min_overlap)
    max_offset = len(left_seq) - min_overlap
    for offset in range(min_offset, max_offset + 1):
        left_start = max(0, offset)
        right_start = max(0, -offset)
        overlap_len = min(
            len(left_seq) - left_start,
            len(right_seq) - right_start,
        )
        if overlap_len < min_overlap:
            continue
        if left_seq[left_start : left_start + overlap_len] != right_seq[right_start : right_start + overlap_len]:
            continue
        merged_start = min(0, offset)
        merged_end = max(len(left_seq), offset + len(right_seq))
        merged_chars: List[str] = []
        conflict = False
        for pos in range(merged_start, merged_end):
            left_char = left_seq[pos] if 0 <= pos < len(left_seq) else ""
            right_pos = pos - offset
            right_char = right_seq[right_pos] if 0 <= right_pos < len(right_seq) else ""
            if left_char and right_char and left_char != right_char:
                conflict = True
                break
            merged_chars.append(left_char or right_char)
        if conflict:
            continue
        candidates.append((overlap_len, len(merged_chars), "".join(merged_chars)))

    if not candidates:
        return None
    candidates.sort(key=lambda item: (-item[0], item[1], item[2]))
    best_overlap, best_len, best_seq = candidates[0]
    tied = [s for o, ln, s in candidates if o == best_overlap and ln == best_len]
    if len(set(tied)) != 1:
        return None
    return best_seq


def _try_assemble_overlap(
    unique_records: List[dict],
    group_key: str,
    *,
    min_overlap: int = 100,
) -> Optional[Tuple[str, dict]]:
    """Try to assemble a full-length sequence from overlapping fragments.

    Returns (assembled_sequence, anchor_record) or None.
    """
    sorted_recs = sorted(unique_records, key=lambda r: (-len(r["sequence"]), r["allele_normalized"]))
    if len(sorted_recs) <= 1:
        return None

    anchor = sorted_recs[0]
    consensus = anchor["sequence"].strip().upper()
    if not consensus:
        return None

    for rec in sorted_recs[1:]:
        seq = rec["sequence"].strip().upper()
        if not seq:
            return None
        merged = _merge_sequences_with_exact_overlap(consensus, seq, min_overlap=min_overlap)
        if merged is None:
            return None
        consensus = merged

    # Verify the assembled sequence produces a valid groove
    groove = _try_groove_parse(
        consensus,
        mhc_class=anchor.get("mhc_class", ""),
        gene=anchor.get("gene", ""),
        allele=group_key,
    )
    if groove is None or not groove.ok:
        return None

    return consensus, anchor


def _groove_signature(row: dict) -> Optional[Tuple[str, str]]:
    """Compute groove signature for a raw record (parses on the fly)."""
    groove = _try_groove_parse(
        row["sequence"],
        mhc_class=row.get("mhc_class", ""),
        gene=row.get("gene", ""),
        allele=row.get("allele_normalized", ""),
    )
    if groove is None or not groove.ok:
        return None
    return (groove.groove1.upper(), groove.groove2.upper())


def _sp_len_from_row(row: dict) -> int:
    """Get signal peptide length (= mature start offset) from a raw CSV row."""
    try:
        return int(row.get("signal_peptide_len", "0") or "0")
    except (ValueError, TypeError):
        return 0


def _pairwise_mismatches(a: str, b: str) -> List[Tuple[int, str, str]]:
    """Find mismatches between two sequences aligned at position 0."""
    mm = []
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]:
            mm.append((i, a[i], b[i]))
    return mm


def _classify_group(
    group_key: str,
    members: List[dict],
) -> Tuple[Optional[dict], Optional[str], str, Optional[dict]]:
    """Classify a two-field allele group and select a representative.

    Returns (representative_row_or_None, assembled_seq_or_None, policy, details).
    Policy is one of: unique, nested_longest, nested_longest_mature,
    assembled_overlap, groove_equivalent, conflict_same_length, conflict_non_nested.
    """
    unique = _unique_by_sequence(members)

    # 1. All map to one protein → trivial
    if len(unique) <= 1:
        rep = unique[0] if unique else members[0]
        return rep, None, "unique", None

    # 2. Check nested containment: does the longest contain all others?
    sorted_by_len = sorted(unique, key=lambda r: (-len(r["sequence"]), r["allele_normalized"]))
    longest = sorted_by_len[0]
    longest_seq = longest["sequence"]
    all_nested = all(r["sequence"] in longest_seq for r in sorted_by_len[1:])
    if all_nested:
        return longest, None, "nested_longest", None

    # 2b. Mature-protein alignment: sequences may differ only in signal peptide.
    #     Strip signal peptides and re-check nesting on the mature portions.
    mature_data = []
    has_any_sp = False
    for r in sorted_by_len:
        ms = _sp_len_from_row(r)
        if ms > 0:
            has_any_sp = True
        mature = r["sequence"][ms:].upper()
        mature_data.append((r, mature, ms))

    if has_any_sp:
        mature_sorted = sorted(
            mature_data,
            key=lambda x: (-len(x[1]), x[0]["allele_normalized"]),
        )
        longest_mature_seq = mature_sorted[0][1]

        if longest_mature_seq:
            all_nested_mature = all(m_seq in longest_mature_seq for _, m_seq, _ in mature_sorted[1:] if m_seq)
            if all_nested_mature:
                details = {
                    "type": "mature_aligned",
                    "alleles": [(r["allele_normalized"], len(r["sequence"]), ms) for r, _, ms in mature_data],
                }
                return longest, None, "nested_longest_mature", details

    # 3. Try overlap assembly
    assembled = _try_assemble_overlap(unique, group_key)
    if assembled is not None:
        assembled_seq, anchor = assembled
        return anchor, assembled_seq, "assembled_overlap", None

    # 4. Groove equivalence: all unique sequences produce the same groove?
    signatures: Dict[Tuple[str, str], List[dict]] = {}
    all_have_sig = True
    for rec in unique:
        sig = _groove_signature(rec)
        if sig is None:
            all_have_sig = False
            break
        signatures.setdefault(sig, []).append(rec)

    if all_have_sig and len(signatures) == 1:
        # All produce the same groove — pick the longest
        exemplar = sorted_by_len[0]
        return exemplar, None, "groove_equivalent", None

    # 5. Conflict — cannot resolve.  Collect details for the merge report.
    conflict_alleles = []
    for r in sorted_by_len:
        try:
            two_field_check = normalize_allele_name(r["allele_normalized"])
        except Exception:
            two_field_check = r.get("two_field_allele", "")
        conflict_alleles.append(
            {
                "allele": r["allele_normalized"],
                "length": len(r["sequence"]),
                "mature_start": _sp_len_from_row(r),
                "two_field_check": two_field_check,
            }
        )
    # Find the first non-nesting pair (mature-aligned) for the report
    mature_mismatches = []
    mm_pair = ("", "")
    mm_offset = 0
    if len(mature_data) >= 2:
        mature_sorted_for_mm = sorted(
            mature_data,
            key=lambda x: (-len(x[1]), x[0]["allele_normalized"]),
        )
        m1_rec, m1, _ = mature_sorted_for_mm[0]
        # Find the first member that doesn't nest in the longest
        m2_rec = None
        m2 = ""
        for rec, mat, _ in mature_sorted_for_mm[1:]:
            if mat and mat not in m1:
                m2_rec = rec
                m2 = mat
                break
        if not m2:
            # All nest — shouldn't happen for conflicts, but just compare top two
            m2_rec = mature_sorted_for_mm[1][0]
            m2 = mature_sorted_for_mm[1][1]
        if m2:
            mm_pair = (m1_rec["allele_normalized"], m2_rec["allele_normalized"])
            # Try substring alignment — find best offset
            best_offset = 0
            best_count = sum(1 for i in range(min(len(m1), len(m2))) if m1[i] != m2[i])
            for off in range(1, max(1, len(m1) - len(m2) + 1)):
                count = sum(1 for i in range(len(m2)) if m1[off + i] != m2[i])
                if count < best_count:
                    best_count = count
                    best_offset = off
            mm_offset = best_offset
            for i in range(min(len(m2), len(m1) - best_offset)):
                if m1[best_offset + i] != m2[i]:
                    mature_mismatches.append((best_offset + i, m1[best_offset + i], m2[i]))

    details = {
        "type": "conflict",
        "alleles": conflict_alleles,
        "n_unique": len(unique),
        "unique_lengths": sorted([len(r["sequence"]) for r in unique], reverse=True),
        "mature_mismatches": mature_mismatches[:30],
        "mature_mismatch_pair": mm_pair,
        "mature_mismatch_offset": mm_offset,
        "two_field_consistent": all(a["two_field_check"] == group_key for a in conflict_alleles),
    }

    lengths = sorted({len(r["sequence"]) for r in unique})
    if len(lengths) == 1:
        return None, None, "conflict_same_length", details
    return None, None, "conflict_non_nested", details


def _emit_full_row(
    group_key: str,
    representative: dict,
    policy: str,
    seq: str,
    groove: Optional[AlleleRecord],
) -> dict:
    """Build a row for mhc-full-seqs.csv (includes groove decomposition)."""
    mature_start = groove.mature_start if groove and groove.ok else 0
    mature_seq = seq[mature_start:] if mature_start > 0 else seq

    gene = representative.get("gene", "")
    gene_upper = gene.strip().upper()
    is_non_groove = gene_upper in NON_GROOVE_GENES or gene_upper in ("B2M", "BETA-2-MICROGLOBULIN")

    if groove:
        groove_status = groove.status
    elif is_non_groove:
        groove_status = "not_applicable"
    else:
        groove_status = ""

    suffix = allele_suffix_flags(representative.get("allele_normalized", ""))
    is_functional = (
        groove_status in FUNCTIONAL_GROOVE_STATUSES and not suffix["is_null"] and not suffix["is_pseudogene"]
    )
    return {
        "two_field_allele": group_key,
        "representative_allele": representative.get("allele_normalized", ""),
        "protein_seq_selection": policy,
        "gene": gene,
        "mhc_class": representative.get("mhc_class", ""),
        "chain": representative.get("chain", ""),
        "species": representative.get("species", ""),
        "species_category": representative.get("species_category", ""),
        "species_prefix": representative.get("species_prefix", ""),
        "source": representative.get("source", ""),
        "source_id": representative.get("source_id", ""),
        "seq_len": str(len(seq)),
        "sequence": seq,
        "mature_start": str(mature_start),
        "mature_sequence": mature_seq,
        "groove1": groove.groove1 if groove else "",
        "groove2": groove.groove2 if groove else "",
        "groove_seq": groove.groove_seq if groove else "",
        "groove1_len": str(groove.groove1_len) if groove else "0",
        "groove2_len": str(groove.groove2_len) if groove else "0",
        "ig_domain": groove.ig_domain if groove else "",
        "ig_domain_len": str(groove.ig_domain_len) if groove else "0",
        "tail": groove.tail if groove else "",
        "tail_len": str(groove.tail_len) if groove else "0",
        "groove_status": groove_status,
        "groove_flags": ",".join(groove.flags) if groove and groove.flags else "",
        "anchor_type": groove.anchor_type if groove else "",
        "is_null": str(suffix["is_null"]),
        "is_questionable": str(suffix["is_questionable"]),
        "is_pseudogene": str(suffix["is_pseudogene"]),
        "is_functional": str(is_functional),
    }


def _write_merge_report(
    report_path: Path,
    stats: Counter,
    mature_aligned_entries: List[Tuple[str, dict]],
    conflict_entries: List[Tuple[str, dict]],
) -> None:
    """Write a human-readable merge report alongside the full-seqs CSV."""
    from datetime import date

    lines: List[str] = []
    lines.append("mhcseqs Merge Report")
    lines.append("=" * 60)
    lines.append(f"Generated: {date.today().isoformat()}")
    lines.append("")

    lines.append("Summary")
    lines.append("-" * 40)
    lines.append(f"Total two-field groups processed:     {stats['total_groups']}")
    lines.append(f"  unique:                             {stats['policy_unique']}")
    lines.append(f"  nested_longest:                     {stats['policy_nested_longest']}")
    lines.append(f"  nested_longest [mature-aligned]:    {stats['policy_nested_longest_mature']}")
    lines.append(f"  assembled_overlap:                  {stats['policy_assembled_overlap']}")
    lines.append(f"  groove_equivalent:                  {stats['policy_groove_equivalent']}")
    n_conflict = stats["policy_conflict_same_length"] + stats["policy_conflict_non_nested"]
    lines.append(f"  conflict (unresolved):              {n_conflict}")
    lines.append(f"  longest_no_groove:                  {stats['policy_longest_no_groove']}")
    lines.append(f"  excluded_all_null:                  {stats['excluded_all_null']}")
    lines.append(f"  excluded_all_short (<{MIN_MHC_SEQUENCE_LEN} aa):       {stats['excluded_all_short']}")
    lines.append(f"  with_representative:                {stats['with_representative']}")
    lines.append("")

    # Mature-aligned groups
    if mature_aligned_entries:
        lines.append(f"Mature-aligned groups ({len(mature_aligned_entries)} signal peptide offset conflicts resolved)")
        lines.append("-" * 60)
        for group_key, details in mature_aligned_entries:
            lines.append(f"  {group_key}:")
            for allele, raw_len, sp_len in details.get("alleles", []):
                mature_len = raw_len - sp_len
                lines.append(f"    {allele} — {raw_len} aa (signal: {sp_len}, mature: {mature_len})")
        lines.append("")

    # Conflict groups
    if conflict_entries:
        lines.append(f"Conflict groups ({len(conflict_entries)} remaining)")
        lines.append("-" * 60)
        for group_key, details in conflict_entries:
            n_unique = details.get("n_unique", "?")
            lengths = details.get("unique_lengths", [])
            consistent = details.get("two_field_consistent", True)
            lines.append(f"  {group_key}:")
            lines.append(f"    {n_unique} unique sequences (lengths: {', '.join(str(x) for x in lengths)})")
            if not consistent:
                lines.append("    *** TWO-FIELD NORMALIZATION MISMATCH — alleles may not belong to this group ***")
            for a in details.get("alleles", []):
                flag = ""
                if a["two_field_check"] != group_key:
                    flag = f"  *** normalizes to {a['two_field_check']} ***"
                lines.append(f"    {a['allele']} — {a['length']} aa (mature_start={a['mature_start']}){flag}")
            mismatches = details.get("mature_mismatches", [])
            mm_offset = details.get("mature_mismatch_offset", 0)
            mm_pair = details.get("mature_mismatch_pair", ("", ""))
            if mismatches:
                offset_note = f", offset {mm_offset}" if mm_offset else ""
                pair_note = f" ({mm_pair[0]} vs {mm_pair[1]})" if mm_pair[0] else ""
                lines.append(f"    Mature-aligned mismatches{pair_note}{offset_note}: {len(mismatches)} positions")
                for pos, c1, c2 in mismatches[:10]:
                    lines.append(f"      mature pos {pos}: {c1} vs {c2}")
                if len(mismatches) > 10:
                    lines.append(f"      ... and {len(mismatches) - 10} more")
            lines.append("")

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def build_full_seqs(
    raw_csv: Path,
    out_csv: Path,
    report_path: Optional[Path] = None,
    min_seq_len: int = MIN_MHC_SEQUENCE_LEN,
) -> Dict[str, int]:
    """Select the best two-field representative for each allele group.

    Selection cascade:
    1. Filter out null alleles (suffix N)
    2. Deduplicate by sequence content
    3. If all map to one protein → unique
    4. If the longest contains all shorter seqs → nested_longest
    4b. If nested after stripping signal peptides → nested_longest_mature
    5. Try assembling overlapping fragments → assembled_overlap
    6. If all produce the same binding groove → groove_equivalent
    7. Otherwise → conflict (recorded but still uses longest as fallback)

    Writes a merge report to report_path if provided.
    Returns stats dict.
    """
    raw_records: List[dict] = []
    with open(raw_csv, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            raw_records.append(row)

    # Group by two-field allele
    groups: Dict[str, List[dict]] = {}
    for row in raw_records:
        seq = row.get("sequence", "")
        if len(seq) < min_seq_len:
            continue
        if _looks_like_nucleotide(seq):
            continue
        key = row.get("two_field_allele", "")
        if not key:
            continue
        groups.setdefault(key, []).append(row)

    full_rows: List[dict] = []
    stats = Counter(
        total_groups=0,
        with_representative=0,
        no_representative=0,
        excluded_all_null=0,
        excluded_all_short=0,
        policy_unique=0,
        policy_nested_longest=0,
        policy_nested_longest_mature=0,
        policy_assembled_overlap=0,
        policy_groove_equivalent=0,
        policy_conflict_same_length=0,
        policy_conflict_non_nested=0,
        policy_longest_no_groove=0,
    )

    # Report data
    mature_aligned_entries: List[Tuple[str, dict]] = []
    conflict_entries: List[Tuple[str, dict]] = []

    for group_key, all_members in sorted(groups.items()):
        stats["total_groups"] += 1

        # Filter out null alleles — they should not be representatives
        members = [m for m in all_members if m.get("is_null") != "True"]
        if not members:
            stats["excluded_all_null"] += 1
            continue

        rep, assembled_seq, policy, details = _classify_group(group_key, members)

        # Collect report data
        if details:
            if details.get("type") == "mature_aligned":
                mature_aligned_entries.append((group_key, details))
            elif details.get("type") == "conflict":
                conflict_entries.append((group_key, details))

        if rep is not None:
            # Use the assembled sequence if available, otherwise the rep's own
            seq = assembled_seq if assembled_seq else rep["sequence"]
            groove = _try_groove_parse(
                seq,
                mhc_class=rep.get("mhc_class", ""),
                gene=rep.get("gene", ""),
                allele=rep.get("allele_normalized", ""),
            )
            stats["with_representative"] += 1
            stats[f"policy_{policy}"] += 1
            full_rows.append(_emit_full_row(group_key, rep, policy, seq, groove))
        else:
            # Conflict — use longest as fallback anyway
            sorted_members = sorted(members, key=lambda r: -len(r["sequence"]))
            fallback = sorted_members[0]
            seq = fallback["sequence"]
            groove = _try_groove_parse(
                seq,
                mhc_class=fallback.get("mhc_class", ""),
                gene=fallback.get("gene", ""),
                allele=fallback.get("allele_normalized", ""),
            )
            if groove and groove.ok:
                stats["with_representative"] += 1
                stats[f"policy_{policy}"] += 1
                full_rows.append(_emit_full_row(group_key, fallback, policy, seq, groove))
            else:
                # True failure: no valid groove even on longest
                stats["with_representative"] += 1
                stats["policy_longest_no_groove"] += 1
                full_rows.append(_emit_full_row(group_key, fallback, "longest_no_groove", seq, groove))

    # Count groups excluded because all entries were too short
    all_keys_in_raw = set()
    for row in raw_records:
        key = row.get("two_field_allele", "")
        if key:
            all_keys_in_raw.add(key)
    short_excluded = all_keys_in_raw - set(groups.keys())
    stats["excluded_all_short"] = len(short_excluded)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FULL_FIELDS)
        writer.writeheader()
        writer.writerows(full_rows)

    # Write merge report
    if report_path is None:
        report_path = out_csv.parent / "mhc-merge-report.txt"
    _write_merge_report(report_path, stats, mature_aligned_entries, conflict_entries)

    return dict(stats)
