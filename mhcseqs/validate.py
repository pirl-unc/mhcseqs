"""Post-build sanity checks for mhcseqs output CSVs.

Validates signal peptides, B2M sequences, groove properties, amino acid
composition, and sequence length distributions. Runs automatically at the
end of the build pipeline and produces warnings for anything unexpected.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------
# Constants: expected properties by class / chain
# ---------------------------------------------------------------------------

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYX")

# Expected signal peptide lengths by MHC class (wide ranges to cover species variation)
SP_RANGE_CLASS_I = (15, 35)  # typical 20-24, some species 15-17
SP_RANGE_CLASS_II_ALPHA = (15, 42)  # typical 24-30, DMA sometimes 39
SP_RANGE_CLASS_II_BETA = (15, 35)  # typical 25-31, some species 15-17

# Expected mature protein lengths (wide ranges to accommodate fragments)
MATURE_RANGE_CLASS_I = (70, 365)  # full ~340 aa, fragments down to groove length
MATURE_RANGE_CLASS_II_ALPHA = (70, 265)  # full ~230 aa, fragments shorter
MATURE_RANGE_CLASS_II_BETA = (70, 270)  # full ~237 aa, fragments shorter

# Expected groove lengths
GROOVE_RANGE_CLASS_I_H1 = (80, 100)  # alpha1 half, typical ~91 aa
GROOVE_RANGE_CLASS_I_H2 = (80, 105)  # alpha2 half, typical ~93 aa
GROOVE_RANGE_CLASS_II_ALPHA_H1 = (70, 95)  # alpha1 groove, typical ~84 aa
GROOVE_RANGE_CLASS_II_BETA_H2 = (80, 105)  # beta1 groove, typical ~94 aa

# B2M expected properties (including signal peptide)
B2M_LENGTH_RANGE = (95, 130)  # mature ~99 aa (human), with signal ~119 aa; varies by species
B2M_EXPECTED_CYS_COUNT = (2, 4)  # B2M has a disulfide bond (2 Cys minimum)

# Hydrophobic residues (for signal peptide check)
HYDROPHOBIC = set("AILMFVW")


# ---------------------------------------------------------------------------
# Individual check functions
# ---------------------------------------------------------------------------


def _check_valid_aa(seq: str, label: str) -> List[str]:
    """Check for invalid amino acid characters."""
    invalid = set(seq.upper()) - VALID_AA
    if invalid:
        return [f"{label}: invalid characters in sequence: {sorted(invalid)}"]
    return []


def _check_signal_peptide(row: dict) -> List[str]:
    """Validate signal peptide properties for a raw CSV row."""
    warnings = []
    sp_len = int(row.get("signal_peptide_len", "0") or "0")
    sp_seq = row.get("signal_peptide_seq", "")
    allele = row.get("allele_normalized", "?")
    mhc_class = row.get("mhc_class", "")
    chain = row.get("chain", "")
    has_sp = row.get("has_signal_peptide", "") == "True"

    if not has_sp or sp_len == 0:
        return []

    # SP should start with M (methionine)
    if sp_seq and not sp_seq.upper().startswith("M"):
        warnings.append(f"SP({allele}): signal peptide does not start with M: '{sp_seq[:5]}...'")

    # SP length should be in expected range
    if mhc_class == "I" and chain == "alpha":
        lo, hi = SP_RANGE_CLASS_I
    elif mhc_class == "II" and chain == "alpha":
        lo, hi = SP_RANGE_CLASS_II_ALPHA
    elif mhc_class == "II" and chain == "beta":
        lo, hi = SP_RANGE_CLASS_II_BETA
    else:
        return warnings  # Skip length check for unknown class/chain

    if not (lo <= sp_len <= hi):
        warnings.append(f"SP({allele}): unusual SP length {sp_len} (expected {lo}-{hi})")

    # SP should have hydrophobic core
    if sp_seq:
        hydro_frac = sum(1 for c in sp_seq.upper() if c in HYDROPHOBIC) / len(sp_seq)
        if hydro_frac < 0.30:
            warnings.append(f"SP({allele}): low hydrophobicity ({hydro_frac:.0%}) in signal peptide '{sp_seq[:20]}...'")

    return warnings


def _check_b2m(row: dict) -> List[str]:
    """Validate B2M sequence properties."""
    warnings = []
    seq = row.get("sequence", "")
    allele = row.get("allele_normalized", "?")

    if not seq:
        warnings.append(f"B2M({allele}): empty sequence")
        return warnings

    # Length check
    lo, hi = B2M_LENGTH_RANGE
    if not (lo <= len(seq) <= hi):
        warnings.append(f"B2M({allele}): unusual length {len(seq)} (expected {lo}-{hi})")

    # Cysteine count (B2M has one Ig-fold disulfide bond)
    cys_count = seq.upper().count("C")
    clo, chi = B2M_EXPECTED_CYS_COUNT
    if not (clo <= cys_count <= chi):
        warnings.append(f"B2M({allele}): unexpected Cys count {cys_count} (expected {clo}-{chi})")

    # B2M should NOT start with M if it's a curated mature sequence
    # (our B2M references include signal peptide, so M is expected)

    return warnings


def _check_groove_row(row: dict) -> List[str]:
    """Validate groove extraction results."""
    warnings = []
    allele = row.get("two_field_allele", "?")
    status = row.get("groove_status", "")
    mhc_class = row.get("mhc_class", "")
    chain = row.get("chain", "")
    h1 = row.get("groove1", "")
    h2 = row.get("groove2", "")

    if status in ("not_applicable", "too_short", ""):
        return []

    if status == "ok":
        if mhc_class == "I" and chain == "alpha":
            lo1, hi1 = GROOVE_RANGE_CLASS_I_H1
            lo2, hi2 = GROOVE_RANGE_CLASS_I_H2
            if h1 and not (lo1 <= len(h1) <= hi1):
                warnings.append(f"GROOVE({allele}): class I alpha1 half length {len(h1)} (expected {lo1}-{hi1})")
            if h2 and not (lo2 <= len(h2) <= hi2):
                warnings.append(f"GROOVE({allele}): class I alpha2 half length {len(h2)} (expected {lo2}-{hi2})")

        elif mhc_class == "II" and chain == "alpha":
            lo, hi = GROOVE_RANGE_CLASS_II_ALPHA_H1
            if h1 and not (lo <= len(h1) <= hi):
                warnings.append(f"GROOVE({allele}): class II alpha groove length {len(h1)} (expected {lo}-{hi})")

        elif mhc_class == "II" and chain == "beta":
            lo, hi = GROOVE_RANGE_CLASS_II_BETA_H2
            if h2 and not (lo <= len(h2) <= hi):
                warnings.append(f"GROOVE({allele}): class II beta groove length {len(h2)} (expected {lo}-{hi})")

    # Check for invalid characters in groove sequences
    for label, groove_seq in [("h1", h1), ("h2", h2)]:
        if groove_seq:
            warnings.extend(_check_valid_aa(groove_seq, f"GROOVE({allele}).{label}"))

    return warnings


def _check_mature_sequence(row: dict) -> List[str]:
    """Validate mature protein properties for a full-seqs row."""
    warnings = []
    allele = row.get("two_field_allele", "?")
    mature = row.get("mature_sequence", "")
    mhc_class = row.get("mhc_class", "")
    chain = row.get("chain", "")

    if not mature:
        return []

    # Check sequence content
    warnings.extend(_check_valid_aa(mature, f"MATURE({allele})"))

    # Check excessive X content
    x_count = mature.upper().count("X")
    if x_count > 0 and x_count / len(mature) > 0.05:
        warnings.append(f"MATURE({allele}): {x_count} unknown residues (X) = {x_count / len(mature):.1%}")

    # Length checks
    if mhc_class == "I" and chain == "alpha":
        lo, hi = MATURE_RANGE_CLASS_I
        if not (lo <= len(mature) <= hi):
            warnings.append(f"MATURE({allele}): class I mature length {len(mature)} (expected {lo}-{hi})")
    elif mhc_class == "II" and chain == "alpha":
        lo, hi = MATURE_RANGE_CLASS_II_ALPHA
        if not (lo <= len(mature) <= hi):
            warnings.append(f"MATURE({allele}): class II alpha mature length {len(mature)} (expected {lo}-{hi})")
    elif mhc_class == "II" and chain == "beta":
        lo, hi = MATURE_RANGE_CLASS_II_BETA
        if not (lo <= len(mature) <= hi):
            warnings.append(f"MATURE({allele}): class II beta mature length {len(mature)} (expected {lo}-{hi})")

    return warnings


# ---------------------------------------------------------------------------
# Aggregate statistics
# ---------------------------------------------------------------------------


def _aa_composition_check(sequences: List[str], label: str) -> List[str]:
    """Check amino acid composition of a collection of sequences."""
    warnings = []
    total = Counter()
    for seq in sequences:
        total.update(seq.upper())

    n = sum(total.values())
    if n == 0:
        return []

    # Expected amino acid frequencies (approximate, from typical proteins)
    # Flag if any amino acid is >2x or <0.5x its expected frequency
    expected_freq = {
        "A": 0.074,
        "R": 0.042,
        "N": 0.044,
        "D": 0.059,
        "C": 0.033,
        "E": 0.058,
        "Q": 0.037,
        "G": 0.074,
        "H": 0.029,
        "I": 0.038,
        "L": 0.076,
        "K": 0.072,
        "M": 0.018,
        "F": 0.040,
        "P": 0.050,
        "S": 0.081,
        "T": 0.062,
        "W": 0.013,
        "Y": 0.033,
        "V": 0.068,
    }

    for aa, expected in expected_freq.items():
        observed = total.get(aa, 0) / n
        ratio = observed / expected if expected > 0 else 0
        if ratio > 3.0:
            warnings.append(f"AA({label}): {aa} enriched {ratio:.1f}x ({observed:.1%} vs expected {expected:.1%})")
        elif ratio < 0.2 and observed > 0:
            warnings.append(f"AA({label}): {aa} depleted {ratio:.1f}x ({observed:.1%} vs expected {expected:.1%})")

    return warnings


# ---------------------------------------------------------------------------
# Main validation entry point
# ---------------------------------------------------------------------------


def validate_build(
    raw_csv: Path,
    full_csv: Path,
    groove_csv: Path,
) -> Tuple[List[str], Dict[str, int]]:
    """Run all sanity checks on the three output CSVs.

    Returns (warnings, stats) where warnings is a list of warning strings
    and stats is a dict of check counts.
    """
    warnings: List[str] = []
    stats = Counter()

    # --- Raw CSV checks ---
    raw_rows = []
    sp_sequences = []
    b2m_rows = []
    with open(raw_csv, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            raw_rows.append(row)

    stats["raw_total"] = len(raw_rows)

    for row in raw_rows:
        # Sequence content
        seq = row.get("sequence", "")
        warnings.extend(_check_valid_aa(seq, f"RAW({row.get('allele_normalized', '?')})"))

        # Signal peptide checks
        sp_warnings = _check_signal_peptide(row)
        warnings.extend(sp_warnings)
        if row.get("has_signal_peptide") == "True":
            stats["raw_with_sp"] += 1
            sp_seq = row.get("signal_peptide_seq", "")
            if sp_seq:
                sp_sequences.append(sp_seq)

        # B2M checks
        if row.get("gene", "").upper() == "B2M":
            b2m_rows.append(row)
            warnings.extend(_check_b2m(row))

    stats["raw_b2m"] = len(b2m_rows)

    # SP aggregate: should most start with M?
    if sp_sequences:
        m_start = sum(1 for sp in sp_sequences if sp.upper().startswith("M"))
        m_pct = m_start / len(sp_sequences) * 100
        stats["sp_start_with_M"] = m_start
        stats["sp_total"] = len(sp_sequences)
        if m_pct < 95:
            warnings.append(
                f"SP_AGGREGATE: only {m_pct:.1f}% of signal peptides start with M ({m_start}/{len(sp_sequences)})"
            )

    # SP composition: signal peptides are expected to be enriched in
    # hydrophobic residues and depleted in charged/polar ones.  Skip the
    # standard composition check (which uses general protein frequencies)
    # and instead verify the hydrophobic enrichment.
    if sp_sequences:
        all_sp = "".join(sp_sequences).upper()
        hydro_frac = sum(1 for c in all_sp if c in HYDROPHOBIC) / len(all_sp) if all_sp else 0
        stats["sp_hydrophobic_fraction"] = round(hydro_frac * 100, 1)
        if hydro_frac < 0.35:
            warnings.append(
                f"SP_AGGREGATE: low hydrophobic fraction ({hydro_frac:.1%}); expected >35% for signal peptides"
            )

    # SP length distribution
    sp_lens = [len(sp) for sp in sp_sequences]
    if sp_lens:
        stats["sp_len_min"] = min(sp_lens)
        stats["sp_len_max"] = max(sp_lens)
        stats["sp_len_mean"] = round(sum(sp_lens) / len(sp_lens), 1)

    # --- Full seqs CSV checks ---
    full_rows = []
    mature_by_class: Dict[str, List[str]] = {"I_alpha": [], "II_alpha": [], "II_beta": []}
    with open(full_csv, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            full_rows.append(row)

    stats["full_total"] = len(full_rows)

    for row in full_rows:
        mature_warnings = _check_mature_sequence(row)
        warnings.extend(mature_warnings)

        mature = row.get("mature_sequence", "")
        mc = row.get("mhc_class", "")
        ch = row.get("chain", "")
        if mature:
            if mc == "I" and ch == "alpha":
                mature_by_class["I_alpha"].append(mature)
            elif mc == "II" and ch == "alpha":
                mature_by_class["II_alpha"].append(mature)
            elif mc == "II" and ch == "beta":
                mature_by_class["II_beta"].append(mature)

    # Amino acid composition of mature sequences by class
    for class_label, seqs in mature_by_class.items():
        if seqs:
            warnings.extend(_aa_composition_check(seqs, f"mature_{class_label}"))
            stats[f"mature_{class_label}_count"] = len(seqs)
            stats[f"mature_{class_label}_mean_len"] = round(sum(len(s) for s in seqs) / len(seqs), 1)

    # --- Groove CSV checks ---
    groove_rows = []
    with open(groove_csv, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            groove_rows.append(row)

    stats["groove_total"] = len(groove_rows)
    groove_ok = 0

    for row in groove_rows:
        groove_warnings = _check_groove_row(row)
        warnings.extend(groove_warnings)
        if row.get("groove_status") == "ok":
            groove_ok += 1

    stats["groove_ok"] = groove_ok

    # Groove length distributions
    for mc, ch, half_key, label in [
        ("I", "alpha", "groove1", "class_I_h1"),
        ("I", "alpha", "groove2", "class_I_h2"),
        ("II", "alpha", "groove1", "class_II_alpha"),
        ("II", "beta", "groove2", "class_II_beta"),
    ]:
        lengths = []
        for row in groove_rows:
            if row.get("mhc_class") == mc and row.get("chain") == ch:
                g = row.get(half_key, "")
                if g:
                    lengths.append(len(g))
        if lengths:
            stats[f"groove_{label}_count"] = len(lengths)
            stats[f"groove_{label}_mean"] = round(sum(lengths) / len(lengths), 1)
            stats[f"groove_{label}_min"] = min(lengths)
            stats[f"groove_{label}_max"] = max(lengths)

    return warnings, dict(stats)


def format_validation_report(
    warnings: List[str],
    stats: Dict[str, int],
) -> str:
    """Format validation results as a human-readable report string."""
    lines = []
    lines.append("Validation Summary")
    lines.append("=" * 60)
    lines.append("")

    # Counts
    lines.append(f"Raw entries:         {stats.get('raw_total', 0)}")
    lines.append(f"  with signal pep:   {stats.get('raw_with_sp', 0)}")
    lines.append(f"  B2M entries:       {stats.get('raw_b2m', 0)}")
    lines.append(f"Full-seq entries:    {stats.get('full_total', 0)}")
    lines.append(f"Groove entries:      {stats.get('groove_total', 0)} ({stats.get('groove_ok', 0)} ok)")
    lines.append("")

    # SP stats
    sp_total = stats.get("sp_total", 0)
    if sp_total:
        sp_min = stats.get("sp_len_min")
        sp_max = stats.get("sp_len_max")
        sp_mean = stats.get("sp_len_mean")
        lines.append(f"Signal peptide lengths: min={sp_min}, max={sp_max}, mean={sp_mean}")
        m_start = stats.get("sp_start_with_M", 0)
        lines.append(f"  Start with M: {m_start}/{sp_total} ({m_start / sp_total * 100:.1f}%)")
        hydro = stats.get("sp_hydrophobic_fraction", 0)
        if hydro:
            lines.append(f"  Hydrophobic fraction: {hydro}%")
        lines.append("")

    # Mature sequence stats
    for label in ["I_alpha", "II_alpha", "II_beta"]:
        count = stats.get(f"mature_{label}_count", 0)
        if count:
            lines.append(f"Mature {label}: n={count}, mean_len={stats.get(f'mature_{label}_mean_len')}")

    lines.append("")

    # Groove length stats
    for label in ["class_I_h1", "class_I_h2", "class_II_alpha", "class_II_beta"]:
        count = stats.get(f"groove_{label}_count", 0)
        if count:
            lines.append(
                f"Groove {label}: n={count}, "
                f"mean={stats.get(f'groove_{label}_mean')}, "
                f"range=[{stats.get(f'groove_{label}_min')}, {stats.get(f'groove_{label}_max')}]"
            )

    lines.append("")

    # Warnings
    if warnings:
        # Aggregate by type
        by_type: Dict[str, List[str]] = {}
        for w in warnings:
            prefix = w.split("(")[0] if "(" in w else w.split(":")[0]
            by_type.setdefault(prefix, []).append(w)

        lines.append(f"Warnings: {len(warnings)} total")
        lines.append("-" * 40)
        for prefix, ws in sorted(by_type.items()):
            if len(ws) <= 5:
                for w in ws:
                    lines.append(f"  {w}")
            else:
                for w in ws[:3]:
                    lines.append(f"  {w}")
                lines.append(f"  ... and {len(ws) - 3} more {prefix} warnings")
    else:
        lines.append("No warnings!")

    return "\n".join(lines)
