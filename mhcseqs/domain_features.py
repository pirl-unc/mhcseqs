"""Structural domain feature extraction, Cys-pair classification, boundaries, and mutations."""

from __future__ import annotations

import re
from dataclasses import replace
from typing import Optional, Sequence

from .domain_grammar import (
    AA_PROPERTY as _AA_PROPERTY,
)
from .domain_grammar import (
    BOUNDARY_CLASS1_M5 as _BOUNDARY_CLASS1_M5,
)
from .domain_grammar import (
    BOUNDARY_CLASS1_P1 as _BOUNDARY_CLASS1_P1,
)
from .domain_grammar import (
    BOUNDARY_CLASS1_P3 as _BOUNDARY_CLASS1_P3,
)
from .domain_grammar import (
    BOUNDARY_CLASS2A_0 as _BOUNDARY_CLASS2A_0,
)
from .domain_grammar import (
    BOUNDARY_CLASS2A_P2 as _BOUNDARY_CLASS2A_P2,
)
from .domain_grammar import (
    BOUNDARY_CLASS2A_P3 as _BOUNDARY_CLASS2A_P3,
)
from .domain_grammar import (
    BOUNDARY_CLASS2B_M1 as _BOUNDARY_CLASS2B_M1,
)
from .domain_grammar import (
    BOUNDARY_CLASS2B_P3 as _BOUNDARY_CLASS2B_P3,
)
from .domain_grammar import (
    BOUNDARY_CLASS2B_P5 as _BOUNDARY_CLASS2B_P5,
)
from .domain_grammar import (
    CYS1_BEFORE_CLIKE as _CYS1_BEFORE_CLIKE,
)
from .domain_grammar import (
    CYS1_BEFORE_GDOMAIN as _CYS1_BEFORE_GDOMAIN,
)
from .domain_grammar import (
    CYS2_AFTER_CLIKE as _CYS2_AFTER_CLIKE,
)
from .domain_grammar import (
    CYS2_AFTER_GDOMAIN as _CYS2_AFTER_GDOMAIN,
)
from .domain_grammar import (
    GROOVE_I_CYS1_M1 as _GROOVE_I_CYS1_M1,
)
from .domain_grammar import (
    GROOVE_I_CYS1_M3 as _GROOVE_I_CYS1_M3,
)
from .domain_grammar import (
    GROOVE_I_CYS1_M5 as _GROOVE_I_CYS1_M5,
)
from .domain_grammar import (
    GROOVE_I_CYS2_P3 as _GROOVE_I_CYS2_P3,
)
from .domain_grammar import (
    GROOVE_I_CYS2_P4 as _GROOVE_I_CYS2_P4,
)
from .domain_grammar import (
    GROOVE_II_CYS1_M1 as _GROOVE_II_CYS1_M1,
)
from .domain_grammar import (
    GROOVE_II_CYS1_M3 as _GROOVE_II_CYS1_M3,
)
from .domain_grammar import (
    GROOVE_II_CYS1_M5 as _GROOVE_II_CYS1_M5,
)
from .domain_grammar import (
    GROOVE_II_CYS2_P3 as _GROOVE_II_CYS2_P3,
)
from .domain_grammar import (
    GROOVE_II_CYS2_P4 as _GROOVE_II_CYS2_P4,
)
from .domain_grammar import (
    IG_CYS2_M2 as _IG_CYS2_M2,
)
from .domain_grammar import (
    IG_CYS2_P2 as _IG_CYS2_P2,
)
from .domain_grammar import (
    IG_I_CYS1_M1 as _IG_I_CYS1_M1,
)
from .domain_grammar import (
    IG_I_CYS1_M2 as _IG_I_CYS1_M2,
)
from .domain_grammar import (
    IG_I_CYS1_P1 as _IG_I_CYS1_P1,
)
from .domain_grammar import (
    IG_I_CYS1_P2 as _IG_I_CYS1_P2,
)
from .domain_grammar import (
    IG_II_CYS1_M1 as _IG_II_CYS1_M1,
)
from .domain_grammar import (
    IG_II_CYS2_P4 as _IG_II_CYS2_P4,
)
from .domain_grammar import (
    IG_SEP_MAX,
    IG_SEP_MIN,
)
from .domain_grammar import (
    JUNCTION_POS0 as _JUNCTION_POS0,
)
from .domain_grammar import (
    JUNCTION_POS1 as _JUNCTION_POS1,
)
from .domain_grammar import (
    JUNCTION_POS2 as _JUNCTION_POS2,
)
from .domain_grammar import (
    JUNCTION_POS3 as _JUNCTION_POS3,
)
from .domain_grammar import (
    JUNCTION_POS5 as _JUNCTION_POS5,
)
from .domain_grammar import (
    SP_CHARGED as _SP_CHARGED,
)
from .domain_grammar import (
    SP_HYDROPHOBIC as _SP_HYDROPHOBIC,
)
from .domain_grammar import (
    TM_CHARGED as _TM_CHARGED,
)
from .domain_grammar import (
    TM_HYDROPHOBIC as _TM_HYDROPHOBIC,
)
from .domain_types import (
    AlleleRecord,
    CysPairAnnotation,
    DomainFoldAnnotation,
    SequenceFeatures,
    StructuralDomain,
    _clean_seq,
    _gaussian_score,
    _grammar_for_record,
)

# ---------------------------------------------------------------------------
# Cysteine pair detection
# ---------------------------------------------------------------------------


def find_cys_pairs(
    seq: str,
    min_sep: int = IG_SEP_MIN,
    max_sep: int = IG_SEP_MAX,
) -> list[tuple[int, int, int]]:
    """Find all Cys-Cys pairs with plausible Ig-fold separation.

    Returns list of (cys1_pos, cys2_pos, separation) tuples.
    """
    cleaned = _clean_seq(seq)
    cys_positions = [idx for idx, aa in enumerate(cleaned) if aa == "C"]
    pairs: list[tuple[int, int, int]] = []
    for i, c1 in enumerate(cys_positions):
        for c2 in cys_positions[i + 1 :]:
            sep = c2 - c1
            if sep < min_sep:
                continue
            if sep > max_sep:
                break
            pairs.append((c1, c2, sep))
    return pairs


# ---- Groove-domain Cys flanking scoring tables ----
def _score_groove_flanking(seq: str, c1: int, c2: int) -> tuple[float, list[str]]:
    """Score how well flanking residues match a groove-domain Cys pair.

    Scores against BOTH class I α2 and class II β1 groove patterns,
    taking the max at each position.  This serves double duty:
    1. Classifying the pair as groove vs Ig (combined with _score_ig_flanking)
    2. Inferring class I vs class II when the MHC class is unknown

    Class I α2 groove: G@-1, Q@-5, M@-3, W@+3, L@+4
    Class II β1 groove: E@-1, K@-3, N@+3, Y@+4
    """
    score = 0.0
    evidence: list[str] = []
    n = len(seq)

    # Separation: groove pairs cluster at ~63 (both class I and II)
    sep = c2 - c1
    sep_s = _gaussian_score(sep, 63, 4) * 2.0
    score += sep_s
    evidence.append(f"sep={sep}({sep_s:+.1f})")

    # Cys1 flanking — take max of class I and class II scores at each position
    for offset, tables_I, tables_II, label in [
        (-1, _GROOVE_I_CYS1_M1, _GROOVE_II_CYS1_M1, "c1-1"),
        (-5, _GROOVE_I_CYS1_M5, _GROOVE_II_CYS1_M5, "c1-5"),
        (-3, _GROOVE_I_CYS1_M3, _GROOVE_II_CYS1_M3, "c1-3"),
    ]:
        pos = c1 + offset
        if 0 <= pos < n:
            aa = seq[pos]
            s = max(tables_I.get(aa, 0.0), tables_II.get(aa, 0.0))
            if s:
                score += s
                evidence.append(f"{aa}@{label}({s:+.1f})")

    # Cys2 flanking — same dual-class scoring
    for offset, tables_I, tables_II, label in [
        (3, _GROOVE_I_CYS2_P3, _GROOVE_II_CYS2_P3, "c2+3"),
        (4, _GROOVE_I_CYS2_P4, _GROOVE_II_CYS2_P4, "c2+4"),
    ]:
        pos = c2 + offset
        if 0 <= pos < n:
            aa = seq[pos]
            s = max(tables_I.get(aa, 0.0), tables_II.get(aa, 0.0))
            if s:
                score += s
                evidence.append(f"{aa}@{label}({s:+.1f})")

    # Groove domains almost never have W at Cys2+15 (14% groove vs 51% Ig)
    if c2 + 15 < n and seq[c2 + 15] != "W":
        score += 1.0
        evidence.append("no_W@c2+15(+1.0)")

    return score, evidence


def _score_ig_flanking(seq: str, c1: int, c2: int) -> tuple[float, list[str]]:
    """Score how well flanking residues match an Ig-domain Cys pair.

    Scores universal Ig signals (shared across all Ig domain types) plus
    class-I-specific α3 signals.  Class II Ig domains (α2, β2) match on
    the universal positions; α3 additionally matches on its distinctive
    LR-C-WA motif.

    Universal Ig signals:
      Separation ~56, Cys2-2=Y (64-97%), Cys2+2=V (74-93%)
    Class I α3 bonus: LR-C-WA at Cys1 (L@-2 92%, R@-1 78%, W@+1 78%, A@+2 96%)
    Class II β2 bonus: Cys2+4=H (73%), Cys2+15=W (72%)
    W@Cys2+15: 51% overall Ig, 14% groove — useful but not definitive
    """
    score = 0.0
    evidence: list[str] = []
    n = len(seq)

    # Separation: Ig pairs cluster at ~56
    sep = c2 - c1
    sep_s = _gaussian_score(sep, 56, 3) * 2.0
    score += sep_s
    evidence.append(f"sep={sep}({sep_s:+.1f})")

    # W at Cys2+15 — moderate Ig signal (51% Ig vs 14% groove)
    if c2 + 15 < n:
        if seq[c2 + 15] == "W":
            score += 2.0
            evidence.append("W@c2+15(+2.0)")
        elif seq[c2 + 15] in "LF":
            score += 0.3
            evidence.append(f"{seq[c2 + 15]}@c2+15(+0.3)")

    # Universal Ig Cys2 flanking (shared across all Ig types)
    if c2 >= 2:
        s = _IG_CYS2_M2.get(seq[c2 - 2], 0.0)
        if s:
            score += s
            evidence.append(f"{seq[c2 - 2]}@c2-2({s:+.1f})")
    if c2 + 2 < n:
        s = _IG_CYS2_P2.get(seq[c2 + 2], 0.0)
        if s:
            score += s
            evidence.append(f"{seq[c2 + 2]}@c2+2({s:+.1f})")

    # Class I α3 Cys1 flanking: LR-C-WA (take max of class I and class II)
    if c1 >= 2:
        for tables_I, tables_II, offset, label in [
            (_IG_I_CYS1_M2, {}, -2, "c1-2"),
            (_IG_I_CYS1_M1, _IG_II_CYS1_M1, -1, "c1-1"),
        ]:
            pos = c1 + offset
            if 0 <= pos < n:
                aa = seq[pos]
                s = max(tables_I.get(aa, 0.0), tables_II.get(aa, 0.0) if tables_II else 0.0)
                if s:
                    score += s
                    evidence.append(f"{aa}@{label}({s:+.1f})")
    if c1 + 2 < n:
        for table, offset, label in [
            (_IG_I_CYS1_P1, 1, "c1+1"),
            (_IG_I_CYS1_P2, 2, "c1+2"),
        ]:
            s = table.get(seq[c1 + offset], 0.0)
            if s:
                score += s
                evidence.append(f"{seq[c1 + offset]}@{label}({s:+.1f})")

    # Class II β2 Cys2+4 bonus
    if c2 + 4 < n:
        s = _IG_II_CYS2_P4.get(seq[c2 + 4], 0.0)
        if s:
            score += s
            evidence.append(f"{seq[c2 + 4]}@c2+4({s:+.1f})")

    return score, evidence


def classify_cys_pair(seq: str, c1: int, c2: int) -> CysPairAnnotation:
    """Classify a Cys-Cys pair as groove-domain or Ig-domain using flanking motifs.

    Returns a CysPairAnnotation with groove_score, ig_score, and domain_type.
    Does not use raw position ranges — classification is purely from local
    sequence context around the Cys pair.
    """
    seq = seq.upper()
    groove_score, groove_ev = _score_groove_flanking(seq, c1, c2)
    ig_score, ig_ev = _score_ig_flanking(seq, c1, c2)

    if groove_score > ig_score + 1.5:
        domain_type = "groove"
    elif ig_score > groove_score + 1.5:
        domain_type = "ig"
    else:
        domain_type = "ambiguous"

    evidence = tuple(groove_ev + ig_ev)
    return CysPairAnnotation(
        c1=c1,
        c2=c2,
        separation=c2 - c1,
        groove_score=groove_score,
        ig_score=ig_score,
        domain_type=domain_type,
        evidence=evidence,
    )


def classify_domain_fold(seq: str, c1: int, c2: int) -> DomainFoldAnnotation:
    """Classify a Cys pair as G-domain or C-like by the conserved Trp41 marker.

    The IMGT-conserved Trp41 in C-like (Ig/C1-set) domains sits at
    exactly Cys1 + 14 in sequence space.  Validated across 6145 Cys pairs
    from 2406 ground-truth sequences:

        W at c1+14:  91.1% sensitivity for Ig, 98.8% specificity vs groove
        W at c1+13:  secondary — captures insertions/deletions but noisier

    G-domains (groove/peptide-binding) lack Trp at this position.

    This is a species-independent, gene-independent structural classifier
    that uses the Ig fold hydrophobic core Trp as the primary signal and
    Cys pair separation as secondary confirmation. See [M2] in the module
    docstring for the underlying IMGT domain grammar.
    """
    seq = seq.upper()
    sep = c2 - c1
    n = len(seq)
    evidence: list[str] = [f"sep={sep}"]

    # --- Primary signal: conserved Trp at the canonical position -----------
    # Position c1+14 is the dominant Trp41 site across vertebrate MHC Ig
    # domains (1792/1968 Ig pairs = 91.1%, vs 39/3216 groove = 1.2%).
    # Positions c1+13 and c1+15 are secondary (insertions/deletions).
    trp_pos = -1
    trp_tier = 0  # 0=none, 1=canonical, 2=near-canonical

    # Tier 1: canonical c1+14 Trp (91.1% of Ig pairs, 1.2% of groove pairs)
    pos14 = c1 + 14
    if pos14 < c2 and pos14 < n and seq[pos14] == "W":
        trp_pos = pos14
        trp_tier = 1
        evidence.append("W@c1+14")
    else:
        # Tier 2: c1+13 Trp (captures insertions)
        pos13 = c1 + 13
        if pos13 < c2 and pos13 < n and seq[pos13] == "W":
            trp_pos = pos13
            trp_tier = 2
            evidence.append("W@c1+13")
        # Tier 3: Ile at c1+14 — divergent Ig domains in fish/reptiles
        # substitute Ile for the conserved Trp (4.5% of Ig pairs).
        # W+I at c1+14: 95.6% sensitivity, 95.5% specificity.
        elif pos14 < c2 and pos14 < n and seq[pos14] == "I" and sep <= 60:
            trp_pos = pos14
            trp_tier = 3
            evidence.append("I@c1+14(trp_substitute)")

    # --- Classification ----------------------------------------------------
    if trp_tier == 1:
        # Canonical Trp41 — very strong C-like signal
        fold_type = "c_like"
        confidence = 0.95
    elif trp_tier == 2 and sep <= 62:
        # Near-canonical Trp with Ig-like separation
        fold_type = "c_like"
        confidence = 0.75
    elif trp_tier in (2, 3) and sep <= 62:
        # Trp shift or Ile substitute with Ig-like separation
        fold_type = "c_like"
        confidence = 0.65
    elif trp_tier == 2:
        # Near-canonical Trp but groove-like separation — ambiguous
        fold_type = "ambiguous"
        confidence = 0.50
        evidence.append("trp_near_but_wide_sep")
    elif sep >= 58:
        # No Trp at core positions + groove-like separation
        fold_type = "g_domain"
        confidence = 0.90
    elif sep >= 48:
        # No Trp + short separation — could be degraded Ig or unusual groove
        fold_type = "ambiguous"
        confidence = 0.40
        evidence.append("short_sep_no_trp")
    else:
        fold_type = "ambiguous"
        confidence = 0.30

    return DomainFoldAnnotation(
        c1=c1,
        c2=c2,
        separation=sep,
        fold_type=fold_type,
        trp_position=trp_pos,
        confidence=confidence,
        evidence=tuple(evidence),
    )


def score_cys_flanking_properties(seq: str, c1: int, c2: int) -> tuple[float, float]:
    """Score a Cys pair's flanking property patterns for domain-type classification.

    Returns (g_domain_score, c_like_score) based on the coarsened amino acid
    property patterns in the 3 residues before/after each Cys.  These patterns
    are species-independent and complement the Trp41 classifier.
    """
    n = len(seq)
    g_score = 0.0
    c_score = 0.0

    # 3 aa before Cys1
    if c1 >= 3:
        pat = tuple(_AA_PROPERTY.get(seq[c1 - 3 + i], "") for i in range(3))
        g_score += _CYS1_BEFORE_GDOMAIN.get(pat, 0.0)
        c_score += _CYS1_BEFORE_CLIKE.get(pat, 0.0)

    # 3 aa after Cys2
    if c2 + 3 < n:
        pat = tuple(_AA_PROPERTY.get(seq[c2 + 1 + i], "") for i in range(3))
        g_score += _CYS2_AFTER_GDOMAIN.get(pat, 0.0)
        c_score += _CYS2_AFTER_CLIKE.get(pat, 0.0)

    return g_score, c_score


def _score_pair_topology_support(
    seq: str,
    c1: int,
    c2: int,
    expected_fold: str,
) -> float:
    """Score how well a Cys pair matches the expected fold topology.

    expected_fold is ``"g_domain"`` or ``"c_like"``.  This combines the
    Trp41 fold classifier with coarse flanking-property patterns.  The Trp41
    signal is treated as the dominant C-like feature; the absence of Trp is
    only a mild G-domain bonus because separation is already encoded elsewhere.
    """
    fold = classify_domain_fold(seq, c1, c2)
    g_prop, c_prop = score_cys_flanking_properties(seq, c1, c2)

    if expected_fold == "g_domain":
        score = g_prop * 0.8 - c_prop * 0.25
        if fold.fold_type == "g_domain":
            score += 0.6 * fold.confidence
        elif fold.fold_type == "c_like":
            penalty = 2.2 if fold.trp_position == c1 + 14 else 1.3
            score -= penalty * fold.confidence
        return score

    score = c_prop * 0.8 - g_prop * 0.25
    if fold.fold_type == "c_like":
        bonus = 2.2 if fold.trp_position == c1 + 14 else 1.3
        score += bonus * fold.confidence
    elif fold.fold_type == "g_domain":
        score -= 0.8 * fold.confidence
    return score


def _cached_topology_support(
    seq: str,
    c1: int,
    c2: int,
    expected: str,
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Read a precomputed pair-topology score when available."""
    if features is not None:
        cached = features.topology_scores.get((c1, c2, expected))
        if cached is not None:
            return cached
    return _score_pair_topology_support(seq, c1, c2, expected)


def _cached_between_cys(
    seq: str,
    c1: int,
    c2: int,
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Read a precomputed between-Cys score when available."""
    if features is not None:
        cached = features.between_cys_scores.get((c1, c2))
        if cached is not None:
            return cached
    return _score_between_cys(seq, c1, c2)


def score_junction(seq: str, pos: int) -> float:
    """Score a candidate α1/α2 junction position (class I groove1/groove2 boundary).

    pos is the 0-indexed position of the first α2 residue.  The junction
    motif is G..H..Q at positions 0, +2, +5 (strongest signals), with
    supporting residues at +1 and +3.
    """
    if pos < 0 or pos + 5 >= len(seq):
        return -999.0
    score = _JUNCTION_POS0.get(seq[pos], -0.5)
    score += _JUNCTION_POS2.get(seq[pos + 2], -1.0)
    score += _JUNCTION_POS5.get(seq[pos + 5], -0.5)
    score += _JUNCTION_POS1.get(seq[pos + 1], 0.0)
    score += _JUNCTION_POS3.get(seq[pos + 3], 0.0)
    return score


def find_junction(seq: str, groove_c1: int) -> tuple[int, float]:
    """Find the α1/α2 junction position by scanning near the groove Cys1.

    The junction is structurally ~10 residues before the groove-domain
    Cys1.  Scans ±5 around that expected position and picks the
    highest-scoring candidate.

    Returns (junction_pos, junction_score).
    """
    expected = groove_c1 - 10
    best_pos = expected
    best_score = -999.0
    for pos in range(max(0, expected - 5), min(len(seq), expected + 6)):
        s = score_junction(seq, pos)
        if s > best_score:
            best_score = s
            best_pos = pos
    return best_pos, best_score


def score_groove_ig_boundary(seq: str, pos: int, mhc_class: str, chain: str = "") -> float:
    """Score a candidate groove→Ig domain boundary.

    pos is the first residue of the Ig domain.  Class-specific because
    the structural linker entering the Ig support domain differs:
      Class I:   hydrophobic...Pro linker
      Class II α: Asn...Pro-Pro linker
      Class II β: Arg...Pro linker
    """
    n = len(seq)
    if pos < 5 or pos + 5 >= n:
        return 0.0
    score = 0.0
    if mhc_class == "I":
        score += _BOUNDARY_CLASS1_M5.get(seq[pos - 5], 0.0)
        score += _BOUNDARY_CLASS1_P1.get(seq[pos + 1], -0.5)
        score += _BOUNDARY_CLASS1_P3.get(seq[pos + 3], 0.0)
    elif chain == "alpha":
        score += _BOUNDARY_CLASS2A_0.get(seq[pos], 0.0)
        score += _BOUNDARY_CLASS2A_P2.get(seq[pos + 2], 0.0)
        score += _BOUNDARY_CLASS2A_P3.get(seq[pos + 3], 0.0)
    else:  # beta
        score += _BOUNDARY_CLASS2B_M1.get(seq[pos - 1], 0.0)
        score += _BOUNDARY_CLASS2B_P3.get(seq[pos + 3], -0.5)
        score += _BOUNDARY_CLASS2B_P5.get(seq[pos + 5], 0.0)
    return score


# ---------------------------------------------------------------------------
# Between-cys properties
# ---------------------------------------------------------------------------


def _score_between_cys(seq: str, c1: int, c2: int) -> float:
    """Score properties of the region between the two Cys residues.

    MHC Ig-fold domains have a characteristic mix of hydrophobic (core packing)
    and charged/polar (surface) residues.  Extremely hydrophobic or extremely
    polar regions are unlikely to be real Ig-fold domains.
    """
    if c2 <= c1 + 5 or c2 >= len(seq):
        return 0.0
    region = seq[c1 + 1 : c2]
    n = len(region)
    if n == 0:
        return 0.0

    hydro = sum(1 for c in region if c in _SP_HYDROPHOBIC) / n
    charged = sum(1 for c in region if c in _SP_CHARGED) / n

    score = 0.0
    # Expect 25-45% hydrophobic (Ig fold core)
    if 0.25 <= hydro <= 0.45:
        score += 1.0
    elif hydro > 0.60 or hydro < 0.10:
        score -= 2.0

    # Expect 10-30% charged (surface residues)
    if 0.10 <= charged <= 0.30:
        score += 0.5
    elif charged > 0.45:
        score -= 1.0

    return score


# ---------------------------------------------------------------------------
# Tail properties
# ---------------------------------------------------------------------------


def _score_tail(
    seq: str,
    tail_start: int,
    *,
    tm_span: tuple[int, int] = (0, 0),
) -> float:
    """Score tail region — MHC tails have a transmembrane helix (hydrophobic).

    A stretch of ≥15 hydrophobic residues in the first 40 aa of tail is
    expected for membrane-anchored MHC molecules.  Soluble/truncated
    forms lack this.

    When *tm_span* is provided (from SequenceFeatures), skips the
    expensive windowed TM scan.
    """
    if tail_start >= len(seq) or tail_start <= 0:
        return 0.0
    tm_start, tm_end = tm_span
    if tm_end <= tm_start:
        tm_start, tm_end = _find_tm_span(seq[tail_start:], raw_start=tail_start)
    tm_len = tm_end - tm_start
    if tm_len >= 20:
        return 1.5
    if tm_len >= 15:
        return 0.5
    return 0.0


def _find_tm_span(tail: str, *, raw_start: int = 0) -> tuple[int, int]:
    """Find a plausible single TM helix near the start of the tail region.

    Uses a windowed hydrophobicity scan rather than a pure contiguous-run
    heuristic so it still recognizes helices with a small number of polar
    interruptions, which are common in real MHC membrane anchors.
    """
    cleaned = _clean_seq(tail)
    search = cleaned[:45]
    if len(search) < 15:
        return 0, 0

    best_score = -999.0
    best_span = (0, 0)
    for win_len in range(17, 25):
        if len(search) < win_len:
            continue
        for start in range(0, len(search) - win_len + 1):
            window = search[start : start + win_len]
            hydro = sum(1 for aa in window if aa in _TM_HYDROPHOBIC)
            charged = sum(1 for aa in window if aa in _TM_CHARGED)
            frac = hydro / win_len
            score = frac * 4.0 - charged * 0.8 + win_len * 0.03
            if score > best_score:
                best_score = score
                best_span = (start, start + win_len)

    start, end = best_span
    if end <= start:
        return 0, 0
    window = search[start:end]
    hydro = sum(1 for aa in window if aa in _TM_HYDROPHOBIC)
    charged = sum(1 for aa in window if aa in _TM_CHARGED)
    if hydro / len(window) < 0.68 or charged > 3:
        return 0, 0
    return raw_start + start, raw_start + end


def _cys_evidence(prefix: str, c1: Optional[int], c2: Optional[int]) -> tuple[str, ...]:
    if c1 is None or c2 is None:
        return ()
    return (f"{prefix}_disulfide=raw_C{c1 + 1}-C{c2 + 1}",)


def _append_domain(
    domains: list[StructuralDomain],
    *,
    total_len: int,
    kind: str,
    family: str,
    role: str,
    start: int,
    end: int,
    evidence: Sequence[str] = (),
) -> None:
    start = max(0, min(start, total_len))
    end = max(start, min(end, total_len))
    if end <= start:
        return
    domains.append(
        StructuralDomain(
            kind=kind,
            family=family,
            role=role,
            start=start,
            end=end,
            evidence=tuple(str(item) for item in evidence if item),
        )
    )


def infer_structural_domains(record: AlleleRecord) -> tuple[StructuralDomain, ...]:
    """Materialize a source-backed domain grammar from a parsed record.

    The grammar is intentionally modular:
      signal peptide -> G-domain groove module(s) -> C-like support domain -> TM -> cytoplasmic tail

    This matches the conserved architecture seen across vertebrate MHC proteins
    more closely than short linear motifs do; see [M1]-[M5] in the module
    docstring.
    """
    total_len = record.seq_len or (max(record.mature_start, 0) + len(record.mature_sequence))
    if total_len <= 0:
        return ()

    domains: list[StructuralDomain] = []
    grammar = _grammar_for_record(record)
    mature_start = max(0, min(record.mature_start, total_len))
    cursor = mature_start

    if mature_start > 0:
        _append_domain(
            domains,
            total_len=total_len,
            kind="signal_peptide",
            family="signal_peptide",
            role="signal_peptide",
            start=0,
            end=mature_start,
            evidence=("signal_peptide_boundary",),
        )

    if record.mhc_class == "I":
        _append_domain(
            domains,
            total_len=total_len,
            kind="g_domain",
            family="G-domain",
            role=(grammar.groove_segments[0][1] if grammar is not None else "g_alpha1"),
            start=cursor,
            end=cursor + record.groove1_len,
            evidence=("peptide_binding_module", "junction_scored"),
        )
        cursor += record.groove1_len
        groove2_evidence = ["peptide_binding_module"]
        if record.anchor_type == "alpha2_cys":
            groove2_evidence.append("g_domain_anchor")
            groove2_evidence.extend(_cys_evidence("g_domain", record.anchor_cys1, record.anchor_cys2))
        elif record.anchor_type == "alpha3_cys":
            groove2_evidence.append("inferred_from_downstream_c_like_anchor")
        _append_domain(
            domains,
            total_len=total_len,
            kind="g_domain",
            family="G-domain",
            role=(grammar.groove_segments[1][1] if grammar is not None and len(grammar.groove_segments) > 1 else "g_alpha2"),
            start=cursor,
            end=cursor + record.groove2_len,
            evidence=groove2_evidence,
        )
        cursor += record.groove2_len
        ig_evidence = ["c_like_support_domain"]
        ig_evidence.extend(_cys_evidence("c_like", record.secondary_cys1, record.secondary_cys2))
        if record.anchor_type == "alpha3_cys":
            ig_evidence.extend(_cys_evidence("c_like", record.anchor_cys1, record.anchor_cys2))
        _append_domain(
            domains,
            total_len=total_len,
            kind="c_like_domain",
            family="C-like",
            role=(grammar.support_role if grammar is not None else "c1_alpha3"),
            start=cursor,
            end=cursor + record.ig_domain_len,
            evidence=ig_evidence,
        )
        cursor += record.ig_domain_len
    elif record.chain == "alpha":
        _append_domain(
            domains,
            total_len=total_len,
            kind="g_domain",
            family="G-domain",
            role=(grammar.groove_segments[0][1] if grammar is not None else "g_alpha1"),
            start=cursor,
            end=cursor + record.groove1_len,
            evidence=("peptide_binding_module", "bounded_by_c_like_anchor"),
        )
        cursor += record.groove1_len
        ig_evidence = ["c_like_support_domain", "c_like_anchor"]
        ig_evidence.extend(_cys_evidence("c_like", record.anchor_cys1, record.anchor_cys2))
        _append_domain(
            domains,
            total_len=total_len,
            kind="c_like_domain",
            family="C-like",
            role=(grammar.support_role if grammar is not None else "c1_alpha2"),
            start=cursor,
            end=cursor + record.ig_domain_len,
            evidence=ig_evidence,
        )
        cursor += record.ig_domain_len
    elif record.chain == "beta":
        groove_evidence = ["peptide_binding_module"]
        if record.anchor_type == "beta1_cys":
            groove_evidence.append("g_domain_anchor")
            groove_evidence.extend(_cys_evidence("g_domain", record.anchor_cys1, record.anchor_cys2))
        elif record.secondary_cys1 is not None and record.secondary_cys2 is not None:
            groove_evidence.append("g_domain_anchor")
            groove_evidence.extend(_cys_evidence("g_domain", record.secondary_cys1, record.secondary_cys2))
        else:
            groove_evidence.append("bounded_by_c_like_anchor")
        _append_domain(
            domains,
            total_len=total_len,
            kind="g_domain",
            family="G-domain",
            role=(grammar.groove_segments[0][1] if grammar is not None else "g_beta1"),
            start=cursor,
            end=cursor + record.groove2_len,
            evidence=groove_evidence,
        )
        cursor += record.groove2_len
        ig_evidence = ["c_like_support_domain"]
        if record.anchor_type == "beta2_cys":
            ig_evidence.append("c_like_anchor")
            ig_evidence.extend(_cys_evidence("c_like", record.anchor_cys1, record.anchor_cys2))
        _append_domain(
            domains,
            total_len=total_len,
            kind="c_like_domain",
            family="C-like",
            role=(grammar.support_role if grammar is not None else "c1_beta2"),
            start=cursor,
            end=cursor + record.ig_domain_len,
            evidence=ig_evidence,
        )
        cursor += record.ig_domain_len

    tail_end = min(total_len, cursor + record.tail_len)
    if tail_end > cursor:
        tm_start, tm_end = _find_tm_span(record.tail, raw_start=cursor)
        if tm_start > cursor:
            _append_domain(
                domains,
                total_len=total_len,
                kind="tail_region",
                family="tail_linker",
                role="tail_linker",
                start=cursor,
                end=tm_start,
                evidence=("between_c_like_and_tm",),
            )
        if tm_end > tm_start:
            _append_domain(
                domains,
                total_len=total_len,
                kind="transmembrane",
                family="TM-helix",
                role="transmembrane",
                start=tm_start,
                end=tm_end,
                evidence=("hydrophobic_tail_helix",),
            )
            if tail_end > tm_end:
                _append_domain(
                    domains,
                    total_len=total_len,
                    kind="cytoplasmic_tail",
                    family="cytoplasmic",
                    role="cytoplasmic_tail",
                    start=tm_end,
                    end=tail_end,
                    evidence=("post_tm_tail",),
                )
        else:
            _append_domain(
                domains,
                total_len=total_len,
                kind="tail_region",
                family="tail_support",
                role="tail_region",
                start=cursor,
                end=tail_end,
                evidence=("no_tm_detected",),
            )

    return tuple(domains)


# ---------------------------------------------------------------------------
# Mutation parsing and application
# ---------------------------------------------------------------------------

_MUTATION_RE = re.compile(r"^([A-Y])(\d+)([A-Y])$", re.IGNORECASE)


def _parse_mutation(mut: object) -> tuple[int, str, str]:
    """Parse a mutation spec into (mature_pos, original_aa, mutant_aa).

    Accepts:
      - str: "K66A" (original + position + mutant, IEDB/mhcgnomes format)
      - tuple/list: (66, "K", "A") or (66, "A") — position, [original], mutant
      - mhcgnomes Mutation object: has .pos, .aa_original, .aa_mutant

    Positions are 1-indexed mature protein numbering.
    """
    if isinstance(mut, str):
        m = _MUTATION_RE.match(mut.strip())
        if not m:
            raise ValueError(f"Cannot parse mutation string: {mut!r} (expected format like 'K66A')")
        return (int(m.group(2)), m.group(1).upper(), m.group(3).upper())

    # mhcgnomes Mutation object or similar duck-typed object
    if hasattr(mut, "pos") and hasattr(mut, "aa_mutant"):
        return (int(mut.pos), str(getattr(mut, "aa_original", "") or "").upper(), str(mut.aa_mutant).upper())  # type: ignore[union-attr]

    if isinstance(mut, (tuple, list)):
        if len(mut) == 3:
            pos, orig, mutant = mut
            return (int(pos), str(orig).upper(), str(mutant).upper())
        if len(mut) == 2:
            pos, mutant = mut
            return (int(pos), "", str(mutant).upper())
        raise ValueError(f"Mutation tuple must have 2-3 elements, got {len(mut)}: {mut!r}")

    raise TypeError(f"Unsupported mutation type: {type(mut).__name__}")


def _normalize_mutations(mutations: str | Sequence[object] | None) -> tuple[object, ...]:
    """Normalize a mutation specification input into a tuple of individual specs.

    Accepts ``None``, a single mutation string (e.g. ``"K66A"``), or any
    sequence of mutation specs. Rejects raw ``bytes``/``bytearray`` with a
    clear ``TypeError`` and avoids iterating a single string character by
    character.
    """
    if mutations is None:
        return ()
    if isinstance(mutations, (bytes, bytearray)):
        raise TypeError("mutations must be a string or sequence of mutation specs, not bytes")
    if isinstance(mutations, str):
        stripped = mutations.strip()
        return (stripped,) if stripped else ()
    return tuple(mutations)


def _apply_requested_mutations(result: AlleleRecord, mutations: str | Sequence[object]) -> AlleleRecord:
    """Apply mutations after a successful parse, including fragment and salvage results."""
    normalized = _normalize_mutations(mutations)
    if normalized and result.ok:
        return apply_mutations(result, normalized)
    return result


def _infer_class_i_mature_offset(result: AlleleRecord) -> int:
    """Infer leading mature residues missing from a partial Class I record.

    Many IMGT/HLA and IPD-MHC partial exon-2 entries begin at mature residue 2
    (e.g. 89-aa ``SHSMRYF...`` missing the leading ``G`` at mature position 1),
    with ``mature_start == 0``. When the record has no signal peptide and
    starts cleanly inside G-ALPHA1 (or G-ALPHA2 for α2-only fragments), return
    the 0-based mature offset (e.g. ``1`` when starting at mature position 2)
    so 1-indexed mature mutations index the right residue.
    """
    if result.mhc_class != "I" or result.mature_start > 0:
        return 0

    from .imgt import align_domain_to_imgt

    if result.groove1:
        aln = align_domain_to_imgt(result.groove1, domain="G-ALPHA1")
        if aln.n_term_missing > 0 and not aln.n_term_overhang and aln.identity >= 0.25:
            return aln.n_term_missing
        return 0

    if result.groove2:
        aln = align_domain_to_imgt(result.groove2, domain="G-ALPHA2")
        start_pos = aln.canonical_start_mature_pos
        if start_pos is not None and not aln.n_term_overhang and aln.identity >= 0.25:
            return start_pos - 1
    return 0


def apply_mutations(result: AlleleRecord, mutations: str | Sequence[object] | None) -> AlleleRecord:
    """Apply mutations to an AlleleRecord, returning a new result with mutated sequences.

    Mutations are specified in mature protein numbering (1-indexed).
    Accepts a single mutation string (``"K66A"``) or any sequence of mutation
    specs: strings (``"K66A"``), tuples, or mhcgnomes ``Mutation`` objects.

    For partial Class I records missing leading N-terminal mature residues
    (such as 89-aa partial exon-2 alleles starting at mature residue 2,
    ``SHSMRYF...``), mature coordinates are automatically offset using the
    IMGT G-DOMAIN alignment so ``"K66A"`` still targets canonical mature
    position 66.

    The original groove extraction must have succeeded (result.ok must be True).
    Mutations are applied to the mature sequence reconstructed from domain parts,
    then re-sliced into the same domain boundaries.

    >>> from mhcseqs.domain_parsing import decompose_class_i, apply_mutations
    >>> wt = decompose_class_i("G" * 400, allele="test")  # doctest: +SKIP
    """
    if not result.ok:
        raise ValueError(f"Cannot apply mutations to failed result (status={result.status!r})")

    normalized = _normalize_mutations(mutations)
    parsed = [_parse_mutation(m) for m in normalized]
    if not parsed:
        return result

    # Reconstruct mature sequence from domain parts
    if result.mhc_class == "I":
        mature = result.groove1 + result.groove2 + result.ig_domain + result.tail
    elif result.chain == "alpha":
        mature = result.groove1 + result.ig_domain + result.tail
    else:  # beta
        mature = result.groove2 + result.ig_domain + result.tail

    mature_offset = _infer_class_i_mature_offset(result)

    # Apply mutations
    mature_list = list(mature)
    mutation_strs: list[str] = []
    for pos, orig, mutant in parsed:
        idx = pos - 1 - mature_offset  # convert 1-indexed mature pos to 0-indexed slice index
        if idx < 0 or idx >= len(mature_list):
            if mature_offset > 0:
                first_pos = mature_offset + 1
                last_pos = mature_offset + len(mature_list)
                raise ValueError(
                    f"Mutation position {pos} out of range for mature sequence covering "
                    f"mature positions {first_pos}-{last_pos} (length {len(mature_list)})"
                )
            raise ValueError(f"Mutation position {pos} out of range for mature sequence (length {len(mature_list)})")
        actual = mature_list[idx]
        if orig and actual != orig:
            raise ValueError(f"Mutation {orig}{pos}{mutant}: expected {orig} at mature position {pos}, found {actual}")
        mutation_strs.append(f"{actual}{pos}{mutant}")
        mature_list[idx] = mutant

    mutated = "".join(mature_list)

    # Re-slice into domain parts using original lengths
    g1_len = result.groove1_len
    g2_len = result.groove2_len
    ig_len = result.ig_domain_len

    if result.mhc_class == "I":
        g1 = mutated[:g1_len]
        g2 = mutated[g1_len : g1_len + g2_len]
        ig = mutated[g1_len + g2_len : g1_len + g2_len + ig_len]
        t = mutated[g1_len + g2_len + ig_len :]
    elif result.chain == "alpha":
        g1 = mutated[:g1_len]
        g2 = ""
        ig = mutated[g1_len : g1_len + ig_len]
        t = mutated[g1_len + ig_len :]
    else:  # beta
        g1 = ""
        g2 = mutated[:g2_len]
        ig = mutated[g2_len : g2_len + ig_len]
        t = mutated[g2_len + ig_len :]

    updated_sequence = result.sequence
    if result.sequence and result.sequence[result.mature_start :] == mature:
        updated_sequence = result.sequence[: result.mature_start] + mutated

    return replace(
        result,
        sequence=updated_sequence,
        groove_seq=g1 + g2,
        groove1=g1,
        groove2=g2,
        ig_domain=ig,
        tail=t,
        tail_len=len(t),
        mutations=tuple(mutation_strs),
    )
