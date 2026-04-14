#!/usr/bin/env python3
"""Curate diverse MHC sequences downloaded from UniProt.

Reads the raw download from fetch_diverse_mhc.py and produces a clean CSV
suitable for shipping with the mhcseqs package.

Curation steps:
  1. Classify MHC class and chain from protein/gene names
  2. Extract or construct a species prefix from the organism name
  3. Normalize gene names (split concatenated prefixes, fix bare genes)
  4. Filter out entries where species, gene, or class can't be determined
  5. Filter out genomic loci that aren't real MHC gene names
  6. Apply minimum length threshold (default 80 aa)

Output: mhcseqs/diverse_mhc_sequences.csv (shipped with the package)

Usage:
    python scripts/curate_diverse_mhc.py [--input PATH] [--output PATH] [--min-length N]
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path

# ---------------------------------------------------------------------------
# MHC class / chain classification
# ---------------------------------------------------------------------------

CLASS_I_PATTERN = re.compile(
    r"(class\s*I[^IV]|class\s*Ia\b|class I alpha|class\s*1\b|"
    r"\bBF[12]?\b|\bYF[12]?\b|MHC-Y|"
    r"MHC class I heavy chain|MHC I\b|MhcI\b)",
    re.IGNORECASE,
)
CLASS_II_ALPHA_PATTERN = re.compile(
    r"(class\s*II\S*\s*alpha|class\s*II\s+\bA\b|"
    r"\bD[A-Z]A\d*\b|\bBLA\d*\b|"
    r"MHC-II.*alpha)",
    re.IGNORECASE,
)
CLASS_II_BETA_PATTERN = re.compile(
    r"(class\s*II\S*\s*beta|class\s*II\s+\bB\b|class\s*IIB\b|"
    r"\bD[A-Z]B\d*\b|\bBLB\d*\b|"
    r"MHC-II.*beta)",
    re.IGNORECASE,
)
B2M_PATTERN = re.compile(r"beta.?2.?microglobulin|B2M\b", re.IGNORECASE)

# Genomic locus patterns — not real MHC gene names
GENOMIC_LOCUS_RE = re.compile(
    r"^("
    r"si:|zmp:|zgc:|wu:|"  # zebrafish genomic loci
    r"LOC\d|"  # NCBI locus tags
    r"[A-Z]\d{3,}_\d|"  # e.g., D623_10000991, L345_18537
    r"[A-Z]{3,}\d{2}_\d{5,}|"  # e.g., ROHU_025598, AMELA_G00181650
    r"[A-Z]+_G\d{5,}|"  # e.g., AMEX_G25088, HHUSO_G36862
    r"[A-Z]+_LOCUS\d+|"  # e.g., SPARVUS_LOCUS3216176
    r"BAC[-_]\d|"  # BAC clone names (not gene names)
    r"FSCOSCO|KUDE01|IRJ41|NDU88|DR999|EPR50|PFLUV|"  # genome project IDs
    r"[A-Z0-9]{4,}_[A-Z]?\d{5,}|"  # generic genome scaffold loci
    r"[A-Z]{3,}_[A-Z0-9]{2,}\d{5,}"  # e.g., PODLI_1B036551
    r")",
    re.IGNORECASE,
)


def classify_mhc(protein_name: str, gene_names: str) -> tuple[str, str] | None:
    """Classify MHC class and chain.

    Returns (mhc_class, chain) or None if unclassifiable.
    """
    combined = f"{protein_name} {gene_names}"

    if B2M_PATTERN.search(combined):
        return ("I", "B2M")
    if CLASS_II_ALPHA_PATTERN.search(combined):
        return ("II", "alpha")
    if CLASS_II_BETA_PATTERN.search(combined):
        return ("II", "beta")
    if CLASS_I_PATTERN.search(combined):
        return ("I", "alpha")

    # Generic "MHC class II" or "mhc2" prefix or "MHCII" — infer chain from gene name
    if re.search(r"(class\s*II|\bmhc2|MHCII)", combined, re.IGNORECASE):
        chain = _infer_class_ii_chain(gene_names, protein_name)
        return ("II", chain)

    # Generic histocompatibility / MHC (with or without "major")
    if re.search(r"(histocompat|MHC\b)", combined, re.IGNORECASE):
        # Check for class II indicators first
        if re.search(r"\bII\b|class\s*2\b|DP\s*beta|DQ\s*beta|DR\s*beta", combined, re.IGNORECASE):
            chain = _infer_class_ii_chain(gene_names, protein_name)
            return ("II", chain)
        if re.search(r"alpha.*chain|alpha 1.*alpha 2", combined, re.IGNORECASE):
            return ("I", "alpha")
        if re.search(r"class\s*[1Ia]\b", combined, re.IGNORECASE):
            return ("I", "alpha")
        return ("unknown", "unknown")

    return None


def _infer_class_ii_chain(gene_names: str, protein_name: str) -> str:
    """Infer alpha vs beta for class II from gene/protein names."""
    combined_text = f"{protein_name} {gene_names}"
    if re.search(r"\balpha\d*\b", combined_text, re.IGNORECASE):
        return "alpha"
    if re.search(r"\bbeta\d*\b", combined_text, re.IGNORECASE):
        return "beta"
    # "class IIB" or "class IIA" or "class II, B" or "class II, A" in protein name
    m = re.search(r"class\s*II\s*,?\s*([AB])\b", protein_name, re.IGNORECASE)
    if m:
        return "alpha" if m.group(1).upper() == "A" else "beta"

    combined = f"{gene_names} {protein_name}"
    # Known beta gene patterns (with optional digits/suffix, allow B-LB hyphenated form)
    if re.search(r"(^|\b|mhc2)(BLB|B-LB|DAB|DRB|DQB|DPB|DMB|DOB|DXB|DBB|DEB)[\dI]*\b", combined, re.IGNORECASE):
        return "beta"
    # Known alpha gene patterns (with optional digits)
    if re.search(r"(^|\b|mhc2)(DAA|DRA|DQA|DPA|DMA|DOA|DXA|DNA)\d*\b", combined, re.IGNORECASE):
        return "alpha"

    if gene_names:
        for tok in gene_names.replace(";", " ").split():
            tok_upper = tok.upper()
            # Gene names containing D?A or D?B anywhere (e.g., Crmi-DA-Ex3 → alpha)
            m = re.search(r"[-]D([AB])[-]", tok_upper)
            if m:
                return "alpha" if m.group(1) == "A" else "beta"
            # Terminal A/B after hyphen: Cosp-B1, Sppu-A2
            m = re.search(r"-([AB])\d*$", tok_upper)
            if m:
                return "alpha" if m.group(1) == "A" else "beta"
            # Standard pattern: D[letter][AB] possibly with digits/suffix
            m = re.search(r"D[A-Z]?([AB])", tok_upper)
            if m:
                return "alpha" if m.group(1) == "A" else "beta"
            m = re.match(r"^[A-Z]{2,6}([AB])\d*$", tok_upper)
            if m and len(tok_upper) >= 3:
                return "alpha" if m.group(1) == "A" else "beta"
    return "unknown"


# ---------------------------------------------------------------------------
# Species prefix derivation
# ---------------------------------------------------------------------------


def derive_prefix(organism: str) -> str:
    """Derive 4-letter MHC prefix from Latin name.

    Convention: first 2 letters of genus + first 2 of species epithet.
    E.g., Phasianus colchicus → Phco, Danio rerio → Dare.
    Returns capitalized prefix or empty string.
    """
    # Strip parenthetical common names
    latin = organism.split("(")[0].strip()
    parts = latin.split()
    if len(parts) >= 2:
        return (parts[0][:2] + parts[1][:2]).capitalize()
    if len(parts) == 1:
        return parts[0][:4].capitalize()
    return ""


def derive_species_tag(organism: str) -> str:
    """Derive a CamelCase species tag for non-standard gene identifiers.

    E.g., Mus musculus → MusMusculus, Rattus norvegicus → RattusNorvegicus.
    Used in ~keyword:SpeciesTag|id format.
    """
    latin = organism.split("(")[0].strip()
    parts = latin.split()
    if len(parts) >= 2:
        return parts[0].capitalize() + parts[1].capitalize()
    if len(parts) == 1:
        return parts[0].capitalize()
    return ""


# ---------------------------------------------------------------------------
# Gene name normalization
# ---------------------------------------------------------------------------

# Model-organism MHC prefixes that genome annotation pipelines transfer
# onto unrelated species via homology (e.g. HLA-DRA on a rodent).
_TRANSFERRED_PREFIXES = {"HLA", "PATR", "POPY", "MAMU", "MAFA", "GOGO", "SLA"}

_BARE_GENE_RE = (
    r"^(BF|BLB|BLA|YF|"
    r"DAB|DAA|DRB|DRA|DQA|DQB|DPA|DPB|DMA|DMB|DOA|DOB|DXB|DXA|DYA|DYB|DNA|"
    r"DBA|DBB|DCA|DCB|DDA|DDB|DEA|DEB|"
    r"UBA|UAA|UCA|UDA|UEA|UFA|UGA|UHA|ULA|"
    r"ZAA|ZBA|ZCA|ZDA|SAA)"
)


def _is_opaque_numbering(tok: str, organism_prefix: str) -> bool:
    """Detect species-specific numbering like Crpo94 (prefix + digits only)."""
    if not organism_prefix:
        return False
    # Extract gene part after any prefix separator
    gene_part = tok.split("-", 1)[-1] if "-" in tok else tok
    gene_part = gene_part.split("_", 1)[-1] if "_" in gene_part else gene_part
    gp_upper = gene_part.upper()
    prefix_upper = organism_prefix.upper()
    if gp_upper.startswith(prefix_upper):
        rest = gp_upper[len(prefix_upper) :]
        if rest and rest.isdigit():
            return True
    return False


def _is_canonical_gene(gene_part: str) -> bool:
    """Check if a gene part (after prefix stripping) uses standard MHC nomenclature."""
    return bool(
        re.match(_BARE_GENE_RE, gene_part, re.IGNORECASE)
        or re.match(
            r"^(UA|UB|UC|UD|UE|UF|UG|UH|UK|UL|UM|"
            r"A|B|C|E|F|G|"  # single-letter class I genes (valid with species prefix)
            r"I|II|Ia|Ib|B2[Mm]|DRB|beta|alpha|BF|BLB|MHC|Mhc)\d*$",
            gene_part,
        )
        or re.match(r"^[A-Z][a-z]{1,5}[AB]\d*$", gene_part)
    )


def normalize_gene(gene_names: str, protein_name: str, organism_prefix: str) -> tuple[str, bool]:
    """Extract and normalize the MHC gene name.

    Handles:
    - Standard Prefix-Gene format (Acsc-UA, orni-dba) → Acsc-UA, Orni-DBA
    - Transferred HLA prefixes (hla-dqa1 on non-human) → Xetr-DQA1
    - Concatenated prefix+gene (XimuDAB, SahaI) → Ximu-DAB, Saha-I
    - Concatenated HLA+gene (HLADRB1) → Spto-DRB1
    - Literature prefixes (PochUA with different organism) → Poch-UA
    - Underscore separator (Hyam_DAB1) → Hyam-DAB1
    - Bare gene with no prefix (DAB1*06) → Prefix-DAB1*06 using organism_prefix
    - Opaque species numbering (Crpo94) → empty string (filtered)
    - Genomic loci (si:xxx, LOC123) → empty string (filtered)

    Returns (normalized_gene, is_canonical) where is_canonical is True when
    the gene name uses standard MHC nomenclature (UA, DAB, DRB, etc.) vs
    paper-specific identifiers (dila_a1, MumuTL, Secit, etc.).
    """
    if not gene_names and not protein_name:
        return "", False

    tokens = gene_names.replace(";", " ").split() if gene_names else []

    # Pass 1: find a token in Prefix-Gene format (hyphen or underscore separator)
    for tok in tokens:
        if _is_opaque_numbering(tok, organism_prefix):
            continue

        m = re.match(r"^([A-Za-z]{2,5})([-_])([A-Za-z].*)$", tok)
        if m:
            tok_prefix, _sep, gene_part = m.group(1), m.group(2), m.group(3)
            # Check if the gene part is opaque numbering
            if _is_opaque_numbering(gene_part, organism_prefix):
                continue
            # Transferred model-organism prefixes on unrelated species →
            # re-prefix with the actual organism prefix.  Legitimate uses
            # (e.g. Gogo on Gobio gobio) pass through because capitalize()
            # matches the organism_prefix.
            if tok_prefix.upper() in _TRANSFERRED_PREFIXES and organism_prefix and tok_prefix.capitalize() != organism_prefix:
                return f"{organism_prefix}-{gene_part.upper()}", _is_canonical_gene(gene_part)
            canonical = _is_canonical_gene(gene_part)
            return f"{tok_prefix.capitalize()}-{gene_part.upper()}", canonical

    # Pass 2: find a known MHC gene token and normalize it
    for tok in tokens:
        if GENOMIC_LOCUS_RE.match(tok):
            continue
        if _is_opaque_numbering(tok, organism_prefix):
            continue

        tok_stripped = tok.strip()
        if not tok_stripped:
            continue

        # Concatenated HLA + gene (HLADRB1 → DRB1)
        if tok_stripped.upper().startswith("HLA") and len(tok_stripped) > 3:
            gene_part = tok_stripped[3:]
            if re.match(_BARE_GENE_RE, gene_part, re.IGNORECASE):
                if organism_prefix and organism_prefix.upper() != "HOSA":
                    return f"{organism_prefix}-{gene_part.upper()}", True
                return f"HLA-{gene_part.upper()}", True

        # Concatenated prefix+gene: 4-letter prefix (matching organism) + gene
        # E.g., XimuDAB → Ximu-DAB, MaeuDBB → Maeu-DBB, SahaI → Saha-I
        if organism_prefix and len(tok_stripped) > 4:
            candidate_prefix = tok_stripped[:4].capitalize()
            if candidate_prefix == organism_prefix and tok_stripped[4:5].isupper():
                gene_part = tok_stripped[4:]
                return f"{candidate_prefix}-{gene_part}", _is_canonical_gene(gene_part)

        # Literature prefix: 3-4 letter capitalized prefix + uppercase gene,
        # where prefix differs from organism (PochUA on Zhch → Poch-UA)
        m = re.match(r"^([A-Z][a-z]{2,3})([A-Z][A-Za-z0-9*]+)$", tok_stripped)
        if m:
            lit_prefix, gene_part = m.group(1), m.group(2)
            if lit_prefix != (organism_prefix or ""):
                return f"{lit_prefix}-{gene_part}", _is_canonical_gene(gene_part)

        # Bare MHC gene name (no prefix): DAB1*06, UBA*01, BF1, DRB1, etc.
        if re.match(_BARE_GENE_RE, tok_stripped, re.IGNORECASE):
            if organism_prefix:
                return f"{organism_prefix}-{tok_stripped.upper()}", True
            return tok_stripped.upper(), True

        # Single-letter or short gene names that look MHC-ish: UA, I, B2m, Mhc, etc.
        if re.match(r"^(UA|UB|UC|UD|UE|UF|UG|UM|I|II|Ia|Ib|Mhc|B2[Mm]|DRB|beta|alpha)\d*$", tok_stripped):
            if organism_prefix:
                return f"{organism_prefix}-{tok_stripped}", True
            return tok_stripped, True

        # Gene names ending in standard class II A/B pattern
        if re.match(r"^[A-Z][a-z]{1,5}[AB]\d*$", tok_stripped):
            return tok_stripped, True

    # Pass 3: fall back to first non-locus, non-opaque token
    # These are paper-specific identifiers (dila_a1, MumuTL, Secit, etc.)
    for tok in tokens:
        if GENOMIC_LOCUS_RE.match(tok):
            continue
        if _is_opaque_numbering(tok, organism_prefix):
            continue
        if re.match(r"^[A-Za-z][\w.*:-]*$", tok) and len(tok) <= 30:
            if organism_prefix:
                return f"{organism_prefix}-{tok}", False
            return tok, False

    return "", False


def resolve_gene_annotation(gene_names: str, protein_name: str, organism_prefix: str) -> tuple[str, str, str]:
    """Resolve gene annotation with provenance tracking.

    Returns (normalized_gene, raw_label, gene_status) where:
    - normalized_gene: clean gene name or ""
    - raw_label: original gene identifier (for provenance)
    - gene_status: one of:
        "ok"                   — standard MHC nomenclature (e.g. Sppu-UA)
        "ortholog_transferred" — gene: ~ortho:Species|Source:gene
        "paper_specific"       — gene: ~ref:Species|Accession:id
        "loc"                  — gene: ~loc:Species|LOC_id
        "opaque_unassigned"    — opaque species-specific numbering (Crpo94)
        "inferred"             — B2M inferred from protein name
        "missing"              — no gene name found
    """
    if not gene_names:
        if B2M_PATTERN.search(protein_name or ""):
            return ("B2M", "", "inferred")
        return ("", "", "missing")

    tokens = gene_names.replace(";", " ").split()

    # Check if all non-locus tokens are opaque numbering
    opaque_label = ""
    all_opaque_or_locus = True
    for tok in tokens:
        if GENOMIC_LOCUS_RE.match(tok):
            continue
        if _is_opaque_numbering(tok, organism_prefix):
            # Preserve the raw identifier (gene part after organism prefix)
            gene_part = tok.split("-", 1)[-1] if "-" in tok else tok
            opaque_label = opaque_label or gene_part
            continue
        all_opaque_or_locus = False

    if all_opaque_or_locus and (opaque_label or tokens):
        return ("", opaque_label, "opaque_unassigned")

    gene, is_canonical = normalize_gene(gene_names, protein_name, organism_prefix)
    if gene:
        status = "ok" if is_canonical else "paper_specific"
        return (gene, tokens[0] if tokens else "", status)

    return ("", tokens[0] if tokens else "", "missing")


# Species-specific MHC nomenclature systems and their source species.
# When these patterns appear on other species, it's because NCBI's automated
# annotation pipeline named the gene after its closest reference ortholog.
# H-2 and RT1 nomenclature is legitimate for Mus and Rattus (any species).
# Match by genus prefix (first 2 chars) to avoid enumerating every species.
_MURIDAE_GENUS_PREFIXES = {"Mu", "Ra"}  # Mus, Rattus

_ORTHOLOG_NOMENCLATURE = [
    # (gene regex, mhcgnomes source prefix — must resolve via Species.get())
    (re.compile(r"(^|-)H2[-.]", re.IGNORECASE), "H2"),
    (re.compile(r"(^|-)H-2", re.IGNORECASE), "H2"),
    (re.compile(r"(^|-)RT1[-.]", re.IGNORECASE), "Rano"),
]


def _detect_ortholog_transfer(gene: str, organism_prefix: str) -> str | None:
    """Detect gene names from another species' nomenclature system.

    Returns the mhcgnomes source prefix (e.g. "H2" for H-2 genes)
    if the gene uses a foreign nomenclature, or None if legitimate.
    H-2/RT1 nomenclature is legitimate for any Mus or Rattus species.
    """
    if organism_prefix[:2] in _MURIDAE_GENUS_PREFIXES:
        return None  # Mus/Rattus — H2 and RT1 are their own system
    for pattern, canonical in _ORTHOLOG_NOMENCLATURE:
        if pattern.search(gene):
            return canonical
    return None


# ---------------------------------------------------------------------------
# Structural validation (for gene-less rescue)
# ---------------------------------------------------------------------------


def _structural_rescue(seq: str, mhc_class: str, chain: str) -> bool:
    """Check if a sequence is structurally valid MHC (domain parse succeeds).

    Requires at least 2 Cys residues (conserved disulfide bond) to avoid
    rescuing poly-repeat or non-MHC sequences that pass the lenient
    fragment-fallback path.
    """
    if seq.upper().count("C") < 2:
        return False
    try:
        from mhcseqs.domain_parsing import (
            decompose_class_i,
            decompose_class_ii_alpha,
            decompose_class_ii_beta,
        )

        if mhc_class == "I":
            return decompose_class_i(seq).ok
        if mhc_class == "II":
            if chain == "alpha":
                return decompose_class_ii_alpha(seq).ok
            return decompose_class_ii_beta(seq).ok
    except Exception:
        pass
    return False


# ---------------------------------------------------------------------------
# Single-row curation (testable entry point)
# ---------------------------------------------------------------------------


def curate_row(row: dict, min_length: int = 80) -> tuple[dict | None, Counter]:
    """Curate a single row from the raw diverse MHC download.

    Returns (curated_dict_or_None, stats) where stats tracks what happened.
    """
    stats: Counter = Counter()

    protein_name = row.get("protein_name", "")
    gene_names = row.get("gene_names", "")
    organism = row.get("organism", "")
    seq = row.get("sequence", "")
    length = len(seq)

    if length < min_length:
        stats["short"] += 1
        return None, stats

    cls = classify_mhc(protein_name, gene_names)
    if cls is None:
        stats["unclassifiable"] += 1
        return None, stats
    mhc_class, chain = cls

    combined = f"{protein_name} {gene_names}"
    has_i = CLASS_I_PATTERN.search(combined)
    has_ii = CLASS_II_ALPHA_PATTERN.search(combined) or CLASS_II_BETA_PATTERN.search(combined)
    if has_i and has_ii:
        stats["ambiguous_class"] += 1
        return None, stats

    prefix = derive_prefix(organism)
    if not prefix:
        stats["no_prefix"] += 1
        return None, stats

    gene, raw_gene_label, gene_status = resolve_gene_annotation(gene_names, protein_name, prefix)

    species_tag = derive_species_tag(organism)

    # Detect ortholog-transferred names: gene named after a reference
    # species ortholog by automated annotation (e.g. H2-K1 on a fish).
    # Format: ~ortho:Species|SourcePrefix:ortholog_gene
    if gene and gene_status in ("ok", "paper_specific"):
        source_prefix = _detect_ortholog_transfer(gene, prefix)
        if source_prefix:
            ortholog_name = gene.split("-", 1)[1] if "-" in gene else gene
            gene = f"~ortho:{species_tag}|{source_prefix}:{ortholog_name}"
            gene_status = "ortholog_transferred"

    if not gene:
        # Try structural rescue: keep if sequence parses as valid MHC
        if _structural_rescue(seq, mhc_class, chain):
            stats["rescued_no_gene"] += 1
            if gene_status == "opaque_unassigned":
                stats["opaque_gene_label"] += 1
        else:
            stats["no_gene"] += 1
            return None, stats
    else:
        # Genomic loci that slipped through normalization → encode as ~loc:
        gene_part = gene.split("-", 1)[-1] if "-" in gene else gene
        if GENOMIC_LOCUS_RE.match(gene_part):
            gene = f"~loc:{species_tag}|{gene_part}"
            gene_status = "loc"
            stats["genomic_locus"] += 1

        # Filter compound gene names with both class I and class II indicators
        gene_upper = gene.upper()
        has_ci_gene = bool(re.search(r"\b(UA|UB|UC|UD|UE|UF|UG|BF)\d", gene_upper))
        has_cii_gene = bool(re.search(r"(DAB|DAA|DRB|DRA|DQB|DQA|DPB|DPA|BLB)\d", gene_upper))
        if has_ci_gene and has_cii_gene:
            stats["ambiguous_gene"] += 1
            return None, stats

        # Encode paper-specific genes with ~ref: format
        if gene_status == "paper_specific":
            # Extract the identifier (gene part after species prefix)
            ref_id = gene.split("-", 1)[1] if "-" in gene else gene
            accession = row.get("uniprot_accession", "")
            gene = f"~ref:{species_tag}|{accession}:{ref_id}"

        stats["kept"] += 1

    curated = {
        "uniprot_accession": row.get("uniprot_accession", ""),
        "gene": gene,
        "raw_gene_label": raw_gene_label,
        "gene_status": gene_status,
        "organism": organism,
        "organism_id": row.get("organism_id", ""),
        "length": str(length),
        "mhc_class": mhc_class,
        "chain": chain,
        "is_fragment": row.get("is_fragment", "False"),
        "source_group": row.get("source_group", ""),
        "sequence": seq,
    }
    return curated, stats


# ---------------------------------------------------------------------------
# Main curation pipeline
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data" / "diverse_mhc_raw.csv",
        help="Input raw CSV from fetch_diverse_mhc.py",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "mhcseqs" / "diverse_mhc_sequences.csv",
        help="Output curated CSV (default: mhcseqs/diverse_mhc_sequences.csv)",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=80,
        help="Minimum sequence length (default: 80)",
    )
    args = parser.parse_args()

    if not args.input.exists():
        print(f"Input not found: {args.input}")
        print("Run scripts/fetch_diverse_mhc.py first.")
        return

    with open(args.input, "r", encoding="utf-8") as f:
        raw_rows = list(csv.DictReader(f))
    print(f"Loaded {len(raw_rows)} raw entries from {args.input}")

    kept = []
    stats = Counter()

    for r in raw_rows:
        curated, row_stats = curate_row(r, min_length=args.min_length)
        stats.update(row_stats)
        if curated is not None:
            kept.append(curated)

    # Report
    print(f"\nCuration results ({len(raw_rows)} → {len(kept)}):")
    for k, v in stats.most_common():
        print(f"  {k}: {v}")

    by_class = Counter((r["mhc_class"], r["chain"]) for r in kept)
    print("\nBy class/chain:")
    for (mc, ch), n in by_class.most_common():
        print(f"  class={mc} chain={ch}: {n}")

    by_group = Counter(r["source_group"] for r in kept)
    print("\nBy source group:")
    for g, n in by_group.most_common():
        print(f"  {g}: {n}")

    prefixes = Counter(r["gene"].split("-")[0] if "-" in r["gene"] else "??" for r in kept)
    print(f"\nUnique species prefixes: {len(prefixes)}")

    fields = [
        "uniprot_accession",
        "gene",
        "gene_status",
        "organism",
        "organism_id",
        "length",
        "mhc_class",
        "chain",
        "is_fragment",
        "source_group",
        "sequence",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(kept)

    print(f"\nWrote {len(kept)} entries to {args.output}")


if __name__ == "__main__":
    main()
