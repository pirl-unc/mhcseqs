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
    r"\bBF[12]?\b|MHC-Y|H-2 class I|"
    r"MHC class I heavy chain|MHC I\b|MhcI\b)",
    re.IGNORECASE,
)
CLASS_II_ALPHA_PATTERN = re.compile(
    r"(class\s*II\S*\s*alpha|class\s*II\s+\bA\b|"
    r"\bD[A-Z]A\d*\b|"
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


# ---------------------------------------------------------------------------
# Gene name normalization
# ---------------------------------------------------------------------------


def normalize_gene(gene_names: str, protein_name: str, organism_prefix: str) -> str:
    """Extract and normalize the MHC gene name.

    Handles:
    - Standard Prefix-Gene format (Acsc-UA) → returns as-is
    - Concatenated prefix+gene (XimuDAB, SahaI, PochUA) → Ximu-DAB, Saha-I, Poch-UA
    - Bare gene with no prefix (DAB1*06) → Prefix-DAB1*06 using organism_prefix
    - Genomic loci (si:xxx, LOC123) → returns empty string (filtered)

    Returns the normalized gene name, or empty string if unusable.
    """
    if not gene_names and not protein_name:
        return ""

    tokens = gene_names.replace(";", " ").split() if gene_names else []

    # Pass 1: find a token in standard Prefix-Gene format
    for tok in tokens:
        if re.match(r"^[A-Za-z]{2,5}-[A-Z]", tok):
            # Skip HLA- gene names on non-human organisms (UniProt uses human gene
            # names as homolog identifiers for non-human species)
            if tok.startswith("HLA-") and organism_prefix and organism_prefix != "Hosa":
                # Re-prefix with organism prefix: HLA-DRA → Sppu-DRA
                gene_part = tok[4:]
                return f"{organism_prefix}-{gene_part}"
            return tok

    # Pass 2: find a known MHC gene token and normalize it
    for tok in tokens:
        # Skip genomic loci
        if GENOMIC_LOCUS_RE.match(tok):
            continue

        tok_stripped = tok.strip()
        if not tok_stripped:
            continue

        # Concatenated prefix+gene: 4-letter prefix (matching organism) + gene
        # E.g., XimuDAB → Ximu-DAB, MaeuDBB → Maeu-DBB, SahaI → Saha-I
        if organism_prefix and len(tok_stripped) > 4:
            candidate_prefix = tok_stripped[:4].capitalize()
            if candidate_prefix == organism_prefix and tok_stripped[4:5].isupper():
                gene_part = tok_stripped[4:]
                return f"{candidate_prefix}-{gene_part}"

        # Bare MHC gene name (no prefix): DAB1*06, UBA*01, BF1, DRB1, etc.
        bare_gene_re = r"^(BF|BLB|DAB|DAA|DRB|DRA|DQA|DQB|DMB|DXB|DXA|UBA|UAA|UCA|UDA|UEA|UFA|UGA|DBA|DBB)"
        if re.match(bare_gene_re, tok_stripped, re.IGNORECASE):
            if organism_prefix:
                return f"{organism_prefix}-{tok_stripped}"
            return tok_stripped

        # Single-letter or short gene names that look MHC-ish: UA, I, B2m, Mhc, etc.
        if re.match(r"^(UA|UB|UC|UD|UE|UF|UG|UM|I|II|Ia|Ib|Mhc|B2[Mm]|DRB|beta|alpha)\d*$", tok_stripped):
            if organism_prefix:
                return f"{organism_prefix}-{tok_stripped}"
            return tok_stripped

        # Gene names ending in standard class II A/B pattern
        if re.match(r"^[A-Z][a-z]{1,5}[AB]\d*$", tok_stripped):
            # Could be concatenated prefix, but prefix doesn't match organism
            # Still usable as a gene name
            return tok_stripped

    # Pass 3: fall back to first non-locus token if nothing else matched
    for tok in tokens:
        if not GENOMIC_LOCUS_RE.match(tok):
            # If it looks like a gene name at all (starts with letter, no weird chars)
            if re.match(r"^[A-Za-z][\w.*:-]*$", tok) and len(tok) <= 30:
                if organism_prefix:
                    return f"{organism_prefix}-{tok}"
                return tok

    return ""


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

    # Load raw data
    with open(args.input, "r", encoding="utf-8") as f:
        raw_rows = list(csv.DictReader(f))
    print(f"Loaded {len(raw_rows)} raw entries from {args.input}")

    # Curate
    kept = []
    stats = Counter()

    for r in raw_rows:
        protein_name = r.get("protein_name", "")
        gene_names = r.get("gene_names", "")
        organism = r.get("organism", "")
        seq = r.get("sequence", "")
        length = len(seq)

        # Length filter
        if length < args.min_length:
            stats["short"] += 1
            continue

        # Classify MHC class/chain
        cls = classify_mhc(protein_name, gene_names)
        if cls is None:
            stats["unclassifiable"] += 1
            continue
        mhc_class, chain = cls

        # Filter ambiguous entries that match both class I and class II patterns
        combined = f"{protein_name} {gene_names}"
        has_i = CLASS_I_PATTERN.search(combined)
        has_ii = CLASS_II_ALPHA_PATTERN.search(combined) or CLASS_II_BETA_PATTERN.search(combined)
        if has_i and has_ii:
            stats["ambiguous_class"] += 1
            continue

        # Derive species prefix
        prefix = derive_prefix(organism)
        if not prefix:
            stats["no_prefix"] += 1
            continue

        # Normalize gene name
        gene = normalize_gene(gene_names, protein_name, prefix)
        if not gene:
            stats["no_gene"] += 1
            continue

        # Filter remaining genomic loci that slipped through
        if GENOMIC_LOCUS_RE.match(gene.split("-", 1)[-1] if "-" in gene else gene):
            stats["genomic_locus"] += 1
            continue

        # Filter compound gene names containing both class I and class II indicators
        gene_upper = gene.upper()
        has_ci_gene = bool(re.search(r"\b(UA|UB|UC|UD|UE|UF|UG|BF)\d", gene_upper))
        has_cii_gene = bool(re.search(r"(DAB|DAA|DRB|DRA|DQB|DQA|DPB|DPA|BLB)\d", gene_upper))
        if has_ci_gene and has_cii_gene:
            stats["ambiguous_gene"] += 1
            continue

        kept.append({
            "uniprot_accession": r["uniprot_accession"],
            "gene": gene,
            "organism": organism,
            "organism_id": r.get("organism_id", ""),
            "length": str(length),
            "mhc_class": mhc_class,
            "chain": chain,
            "is_fragment": r.get("is_fragment", "False"),
            "source_group": r.get("source_group", ""),
            "sequence": seq,
        })
        stats["kept"] += 1

    # Report
    print(f"\nCuration results ({len(raw_rows)} → {stats['kept']}):")
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

    # Write output
    fields = [
        "uniprot_accession", "gene", "organism", "organism_id",
        "length", "mhc_class", "chain", "is_fragment",
        "source_group", "sequence",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(kept)

    print(f"\nWrote {stats['kept']} entries to {args.output}")


if __name__ == "__main__":
    main()
