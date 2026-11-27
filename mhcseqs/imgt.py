"""IMGT G-DOMAIN unique numbering for MHC groove domains.

Maps between mature protein positions (1-indexed) and IMGT G-DOMAIN positions
for MHC class I and class II molecules.

The IMGT numbering is per-domain, not per-chain.  Each groove domain
(G-ALPHA1, G-ALPHA2, G-ALPHA, G-BETA) has its own numbering from 1 to ~92.

For Class I (single alpha chain):
    G-ALPHA1 [D1]: alpha-1 domain (mature positions 1-90)
    G-ALPHA2 [D2]: alpha-2 domain (mature positions 91-182)

For Class II (separate alpha and beta chains):
    G-ALPHA [D1]: alpha-1 domain on the alpha chain
    G-BETA  [D1]: beta-1 domain on the beta chain

Reference: Lefranc et al. (2005) Dev Comp Immunol 29:917-938.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

# =============================================================================
# IMGT structural elements
# =============================================================================

#: Structural element boundaries (inclusive) in the IMGT G-DOMAIN numbering.
#: The groove floor comprises four antiparallel beta strands and three turns.
#: The helix (positions 50-92) forms the wall of the groove.
STRUCTURAL_ELEMENTS = {
    "A-STRAND": (1, 14),
    "AB-TURN": (15, 17),
    "B-STRAND": (18, 28),
    "BC-TURN": (29, 30),
    "C-STRAND": (31, 38),
    "CD-TURN": (39, 41),
    "D-STRAND": (42, 49),
    "HELIX": (50, 92),
}

#: Helix insertion positions defined in the IMGT numbering.
#: 54A is only occupied in some G-ALPHA domains (DQA1, DOA, DMA).
#: 61A and 72A are occupied in G-ALPHA2 and G-BETA domains.
#: 61B is occupied in G-BETA domains.
#: 92A is occupied only in HLA-DMA and H2-DMA G-ALPHA domains.
HELIX_INSERTIONS = ("54A", "61A", "61B", "72A", "92A")

#: Conserved disulfide-bridge Cys positions (IMGT numbering).
#: CYS 11 (A-strand) and CYS 74 (helix) form the Ig-fold disulfide bridge
#: in G-ALPHA2, G-BETA, and their G-LIKE-DOMAIN counterparts.
CONSERVED_CYS_POSITIONS = (11, 74)


def structural_element(imgt_pos: str) -> str:
    """Return the structural element for an IMGT G-DOMAIN position.

    >>> structural_element("11")
    'A-STRAND'
    >>> structural_element("61A")
    'HELIX'
    >>> structural_element("1.2")
    'A-STRAND'
    """
    # Handle dot-notation insertions (1.1, 1.2, 49.1, 51.1, 61A.1, etc.)
    if "." in imgt_pos:
        base_str = imgt_pos.split(".", 1)[0]
        if base_str in HELIX_INSERTIONS:
            return "HELIX"
        try:
            base = int(base_str.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
        except ValueError:
            return "unknown"
        if base <= 1:
            return "A-STRAND"
        if base >= 50:
            return "HELIX"
        for element, (start, end) in STRUCTURAL_ELEMENTS.items():
            if element == "HELIX":
                continue
            if start <= base <= end:
                return element
        return "unknown"

    # Handle letter-suffix insertions (19A, 54A-54C, 61A-61C, 72A, 40A-40C, etc.)
    num_part = imgt_pos.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    try:
        num = int(num_part)
    except ValueError:
        return "unknown"

    if num >= 50 or imgt_pos in HELIX_INSERTIONS:
        return "HELIX"

    for element, (start, end) in STRUCTURAL_ELEMENTS.items():
        if element == "HELIX":
            continue
        if start <= num <= end:
            return element

    return "unknown"


# =============================================================================
# Class I G-ALPHA1 [D1] mapping
# =============================================================================
# Trivial: IMGT position = mature position for positions 1-90.
# Some alleles may use insertion position 7A, but for HLA-A the mapping
# is 1:1.  Domain covers the alpha-1 groove half.

G_ALPHA1_SIZE = 90  # residues in the G-ALPHA1 domain


# =============================================================================
# Class I G-ALPHA2 [D2] mapping (HLA-A*02:01 reference)
# =============================================================================
#
# From Lefranc et al. 2005, p.924:
#   "The IMGT unique numbering for positions 1-39 and 73-92 of the
#    G-ALPHA2 domains can be obtained very easily by subtracting 90
#    from the mature protein numbering (91-129 and 163-182).  Between
#    these positions, the two gaps (at positions 40 and 41) and the
#    two insertions (at positions 61A and 72A) are necessary."
#
# Mapping verified against Table 3 of the paper.

#: IMGT G-ALPHA2 positions that are gaps (no corresponding mature residue)
#: in HLA class I molecules.
G_ALPHA2_GAP_POSITIONS = frozenset({"40", "41"})


def _build_g_alpha2_maps() -> tuple[dict[str, int], dict[int, str]]:
    """Build bidirectional IMGT<->mature maps for Class I G-ALPHA2."""
    fwd: dict[str, int] = {}

    # IMGT 1-39 -> mature 91-129 (offset +90)
    for i in range(1, 40):
        fwd[str(i)] = i + 90

    # IMGT 40, 41 -> gaps (CD-TURN unoccupied)

    # IMGT 42-61 -> mature 130-149 (offset +88, because 2 gap positions were skipped)
    for i in range(42, 62):
        fwd[str(i)] = i + 88

    # IMGT 61A -> mature 150 (helix insertion)
    fwd["61A"] = 150

    # IMGT 62-72 -> mature 151-161 (offset +89, shifted by 61A insertion)
    for i in range(62, 73):
        fwd[str(i)] = i + 89

    # IMGT 72A -> mature 162 (helix insertion)
    fwd["72A"] = 162

    # IMGT 73-92 -> mature 163-182 (offset +90, back to original)
    for i in range(73, 93):
        fwd[str(i)] = i + 90

    rev = {v: k for k, v in fwd.items()}
    return fwd, rev


G_ALPHA2_IMGT_TO_MATURE, G_ALPHA2_MATURE_TO_IMGT = _build_g_alpha2_maps()

#: Total residues mapped in G-ALPHA2 (for HLA-A*02:01): 92
G_ALPHA2_SIZE = len(G_ALPHA2_IMGT_TO_MATURE)

#: Ordered IMGT positions for G-ALPHA1 (1 through 90).
G_ALPHA1_POSITIONS: tuple[str, ...] = tuple(str(i) for i in range(1, G_ALPHA1_SIZE + 1))

#: Ordered IMGT positions for G-ALPHA2 (occupied only, excludes gaps).
G_ALPHA2_POSITIONS: tuple[str, ...] = tuple(k for k, v in sorted(G_ALPHA2_IMGT_TO_MATURE.items(), key=lambda x: x[1]))

#: Canonical HLA-A*02:01 reference sequences (UniProt P01892) for IMGT G-DOMAIN
#: alignment: G-ALPHA1 (mature 1-90, 90 aa) and G-ALPHA2 (mature 91-182, 92 aa).
HLA_A0201_G_ALPHA1 = "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEA"
HLA_A0201_G_ALPHA2 = "GSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRT"

# =============================================================================
# Class II G-ALPHA [D1] (HLA-DRA*01:01) and G-BETA [D1] (HLA-DRB1*01:01) maps
# =============================================================================
#
# Verified via 3D Kabsch superposition of PDB 1DLH (HLA-DRA*01:01 / HLA-DRB1*01:01),
# 1JK8 (HLA-DQA1*03:01), and 3LQZ (HLA-DPB1*02:01) against 1HHK (HLA-A*02:01),
# following Lefranc et al. (2005) Dev Comp Immunol 29:917-938.

#: IMGT G-ALPHA positions that are structural gaps in canonical HLA-DRA*01:01.
G_ALPHA_GAP_POSITIONS: frozenset[str] = frozenset({"16", "17", "18", "19", "40", "41", "57", "58", "59", "60"})

#: Ordered occupied IMGT positions for G-ALPHA (HLA-DRA*01:01 reference, 82 aa).
G_ALPHA_POSITIONS: tuple[str, ...] = (
    ("1.2", "1.1")
    + tuple(str(i) for i in range(1, 16))
    + tuple(str(i) for i in range(20, 40))
    + tuple(str(i) for i in range(42, 57))
    + tuple(str(i) for i in range(61, 91))
)

#: Canonical HLA-DRA*01:01 reference sequence (UniProt P01903 mature 1-82, 82 aa).
HLA_DRA0101_G_ALPHA = "IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYTPI"
G_ALPHA_SIZE = len(G_ALPHA_POSITIONS)  # 82
G_ALPHA_MATURE_POSITIONS: tuple[int, ...] = tuple(range(1, G_ALPHA_SIZE + 1))

#: IMGT G-BETA positions that are structural gaps in canonical HLA-DRB1*01:01.
G_BETA_GAP_POSITIONS: frozenset[str] = frozenset({"40", "41"})

#: Ordered occupied IMGT positions for G-BETA (HLA-DRB1*01:01 reference, 89 aa, mature 5-93).
G_BETA_POSITIONS: tuple[str, ...] = (
    tuple(str(i) for i in range(1, 40))
    + tuple(str(i) for i in range(42, 62))
    + ("61A", "61B")
    + tuple(str(i) for i in range(62, 73))
    + ("72A",)
    + tuple(str(i) for i in range(73, 89))
)

#: Canonical HLA-DRB1*01:01 reference sequence (UniProt P04229 mature 5-93, 89 aa).
HLA_DRB1_0101_G_BETA = "PRFLWQLKFECHFFNGTERVRLLERCIYNQEESVRFDSDVGEYRAVTELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVGESFTVQR"
G_BETA_SIZE = len(G_BETA_POSITIONS)  # 89
G_BETA_MATURE_POSITIONS: tuple[int, ...] = tuple(range(5, 5 + G_BETA_SIZE))  # mature 5..93


# =============================================================================
# Public API
# =============================================================================


def mature_to_imgt_class_i(mature_pos: int) -> tuple[str, str]:
    """Convert a 1-indexed mature protein position to an IMGT G-DOMAIN position.

    Returns (domain, imgt_pos) for Class I alpha chain.

    >>> mature_to_imgt_class_i(66)
    ('G-ALPHA1', '66')
    >>> mature_to_imgt_class_i(101)
    ('G-ALPHA2', '11')
    >>> mature_to_imgt_class_i(150)
    ('G-ALPHA2', '61A')
    >>> mature_to_imgt_class_i(200)
    ('C-LIKE', '18')
    """
    if mature_pos < 1:
        raise ValueError(f"mature_pos must be >= 1, got {mature_pos}")

    if mature_pos <= G_ALPHA1_SIZE:
        return ("G-ALPHA1", str(mature_pos))

    if mature_pos in G_ALPHA2_MATURE_TO_IMGT:
        return ("G-ALPHA2", G_ALPHA2_MATURE_TO_IMGT[mature_pos])

    # Beyond G-ALPHA2: C-LIKE domain (alpha-3 Ig-fold)
    if mature_pos > 182:
        return ("C-LIKE", str(mature_pos - 182))

    raise ValueError(f"mature_pos {mature_pos} has no G-ALPHA2 IMGT mapping")


def imgt_to_mature_class_i(domain: str, imgt_pos: str) -> Optional[int]:
    """Convert an IMGT G-DOMAIN position to a 1-indexed mature protein position.

    Returns the mature position, or None if the IMGT position is a gap.

    >>> imgt_to_mature_class_i("G-ALPHA1", "66")
    66
    >>> imgt_to_mature_class_i("G-ALPHA2", "11")
    101
    >>> imgt_to_mature_class_i("G-ALPHA2", "40")  # gap
    >>> imgt_to_mature_class_i("G-ALPHA2", "61A")
    150
    """
    if domain == "G-ALPHA1":
        try:
            pos = int(imgt_pos)
        except ValueError:
            return None
        return pos if 1 <= pos <= G_ALPHA1_SIZE else None

    if domain == "G-ALPHA2":
        if imgt_pos in G_ALPHA2_GAP_POSITIONS:
            return None
        return G_ALPHA2_IMGT_TO_MATURE.get(imgt_pos)

    if domain == "C-LIKE":
        try:
            return int(imgt_pos) + 182
        except ValueError:
            return None

    return None


def mature_to_imgt(mature_pos: int, mhc_class: str = "I", chain: str = "") -> tuple[str, str]:
    """Convert a 1-indexed mature protein position to an IMGT G-DOMAIN position.

    For Class I, both groove domains are on the same alpha chain:
        mature 1-90   -> G-ALPHA1 (trivial 1:1 mapping)
        mature 91-182 -> G-ALPHA2 (with gaps at IMGT 40,41 and insertions at 61A,72A)
        mature 183+   -> C-LIKE (alpha-3 Ig-fold, not groove)

    For Class II, each chain has one groove domain.  The ``chain`` parameter
    ("alpha" or "beta") selects which mapping to use.  Within each chain,
    the IMGT mapping is approximate (exact alignment varies by allele).

    Returns:
        (domain, imgt_pos) tuple.
    """
    if mhc_class == "I":
        return mature_to_imgt_class_i(mature_pos)

    if mhc_class == "II":
        if chain == "alpha":
            # G-ALPHA [D1]: approximate 1:1 for the groove domain
            # Exact mapping requires structural alignment (varies by allele)
            return ("G-ALPHA", str(mature_pos))
        if chain == "beta":
            return ("G-BETA", str(mature_pos))
        raise ValueError("Class II requires chain='alpha' or chain='beta'")

    raise ValueError(f"Unknown mhc_class: {mhc_class!r}")


def imgt_to_mature(domain: str, imgt_pos: str) -> Optional[int]:
    """Convert an IMGT G-DOMAIN position to a 1-indexed mature protein position.

    For Class I domains (G-ALPHA1, G-ALPHA2, C-LIKE), returns the exact
    mature position.  For Class II domains (G-ALPHA, G-BETA), returns an
    approximate position (exact mapping varies by allele).

    Returns None if the IMGT position is a gap.
    """
    if domain in ("G-ALPHA1", "G-ALPHA2", "C-LIKE"):
        return imgt_to_mature_class_i(domain, imgt_pos)

    # Class II: approximate 1:1 mapping within the domain
    if domain in ("G-ALPHA", "G-BETA"):
        try:
            return int(imgt_pos)
        except ValueError:
            return None

    return None


# =============================================================================
# Reference-anchored Class I IMGT G-DOMAIN alignment and indel annotation
# =============================================================================

_BLOSUM62_ALPHABET = "ARNDCQEGHILKMFPSTWYVBZX*"
_BLOSUM62_ROWS: tuple[tuple[int, ...], ...] = (
    (4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4),
    (-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4),
    (-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4),
    (-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4),
    (0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4),
    (-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4),
    (-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4),
    (0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4),
    (-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4),
    (-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4),
    (-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4),
    (-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4),
    (-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4),
    (-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4),
    (-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4),
    (1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4),
    (0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4),
    (-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4),
    (-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4),
    (0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4),
    (-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4),
    (-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4),
    (0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4),
    (-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1),
)
_BLOSUM62: dict[str, dict[str, int]] = {
    aa: {bb: _BLOSUM62_ROWS[i][j] for j, bb in enumerate(_BLOSUM62_ALPHABET)} for i, aa in enumerate(_BLOSUM62_ALPHABET)
}
_BLOSUM62_UNKNOWN = _BLOSUM62["X"]

G_ALPHA2_REF_C1_IDX = 10  # mature 101 = IMGT G-ALPHA2 11
G_ALPHA2_REF_C2_IDX = 73  # mature 164 = IMGT G-ALPHA2 74
G_ALPHA2_CYS_ANCHOR_BONUS = 6.0

G_BETA_REF_C1_IDX = 10  # mature 15 = IMGT G-BETA 11
G_BETA_REF_C2_IDX = 74  # mature 79 = IMGT G-BETA 74
G_BETA_CYS_ANCHOR_BONUS = 6.0

#: Shipped multi-clade Class I and Class II IMGT G-DOMAIN structural PSSM + loop gap profiles.
IMGT_CLASS_I_PROFILE_PATH = Path(__file__).resolve().parent / "imgt_class_i_profile.json"
IMGT_CLASS_II_PROFILE_PATH = Path(__file__).resolve().parent / "imgt_class_ii_profile.json"

#: 22-token alphabet for fixed-index IMGT frame strings:
#: 20 canonical amino acids + '-' (structural gap) + 'X' (missing/unknown).
FIXED_FRAME_ALPHABET = "ARNDCQEGHILKMFPSTWYV-X"
CANONICAL_AA20 = frozenset("ARNDCQEGHILKMFPSTWYV")


@dataclass(frozen=True)
class InsertionHotspot:
    """Named structural loop or helix-hinge insertion hotspot in an IMGT fixed frame.

    Each hotspot sits strictly between ``anchor_imgt_pos`` (left flank) and
    ``right_flank_imgt_pos`` (right flank) in its target slot table, and ONLY
    absorbs insertions with ``after_imgt_pos == anchor_imgt_pos`` — guaranteeing
    monotonic left-to-right sequence ordering with zero residue transposition.
    """

    name: str
    structural_element: str
    anchor_imgt_pos: str
    right_flank_imgt_pos: str
    slots: tuple[str, ...]
    fill_mode: str = "left"
    lineage: str = ""


#: Fixed-index IMGT G-ALPHA1 slot schema (96 slots):
#: 90 canonical match columns + AB-turn slots (19A-19C) + N-helix hinge slots (54A-54C).
IMGT_G_ALPHA1_SLOTS: tuple[str, ...] = (
    tuple(str(i) for i in range(1, 20))
    + ("19A", "19B", "19C")
    + tuple(str(i) for i in range(20, 55))
    + ("54A", "54B", "54C")
    + tuple(str(i) for i in range(55, 91))
)

#: Fixed-index IMGT G-ALPHA2 slot schema (102 slots):
#: 94 IMGT backbone positions (including CD-turn 40, 41 and helix 61A, 72A)
#: + AB-turn slots (13A-13C) + CD-turn slots (40A-40C) + helix-kink slots (61B-61C).
IMGT_G_ALPHA2_SLOTS: tuple[str, ...] = (
    tuple(str(i) for i in range(1, 14))
    + ("13A", "13B", "13C")
    + tuple(str(i) for i in range(14, 40))
    + ("40", "40A", "40B", "40C", "41")
    + tuple(str(i) for i in range(42, 62))
    + ("61A", "61B", "61C")
    + tuple(str(i) for i in range(62, 73))
    + ("72A",)
    + tuple(str(i) for i in range(73, 93))
)

#: Combined Class I fixed-index IMGT groove frame (198 slots: 96 G-ALPHA1 + 102 G-ALPHA2).
IMGT_CLASS_I_GROOVE_SLOTS: tuple[tuple[str, str], ...] = tuple(("G-ALPHA1", s) for s in IMGT_G_ALPHA1_SLOTS) + tuple(
    ("G-ALPHA2", s) for s in IMGT_G_ALPHA2_SLOTS
)

#: Fixed-index IMGT G-ALPHA slot schema (87 slots, zero dead Class I-only columns):
#: 82 canonical match columns + A-strand 7A (DQA1) + CD-turn 40-41 + N-helix hinge 54A-54B.
IMGT_G_ALPHA_SLOTS: tuple[str, ...] = (
    ("1.2", "1.1")
    + tuple(str(i) for i in range(1, 8))
    + ("7A",)
    + tuple(str(i) for i in range(8, 16))
    + tuple(str(i) for i in range(20, 40))
    + ("40", "41")
    + tuple(str(i) for i in range(42, 55))
    + ("54A", "54B")
    + ("55", "56")
    + tuple(str(i) for i in range(61, 91))
)

#: Fixed-index IMGT G-BETA slot schema (98 slots, zero dead columns):
#: 89 canonical match columns + AB-turn 13A-13C + CD-turn 40/40A-40C/41 + helix-kink 61C.
IMGT_G_BETA_SLOTS: tuple[str, ...] = (
    tuple(str(i) for i in range(1, 14))
    + ("13A", "13B", "13C")
    + tuple(str(i) for i in range(14, 40))
    + ("40", "40A", "40B", "40C", "41")
    + tuple(str(i) for i in range(42, 62))
    + ("61A", "61B", "61C")
    + tuple(str(i) for i in range(62, 73))
    + ("72A",)
    + tuple(str(i) for i in range(73, 89))
)

#: Combined Class II fixed-index IMGT groove frame (185 slots: 87 G-ALPHA + 98 G-BETA).
IMGT_CLASS_II_GROOVE_SLOTS: tuple[tuple[str, str], ...] = tuple(("G-ALPHA", s) for s in IMGT_G_ALPHA_SLOTS) + tuple(
    ("G-BETA", s) for s in IMGT_G_BETA_SLOTS
)

#: Shared 3D-homologous Pan-Class Half 1 slot schema (100 slots, shared by G-ALPHA1 and G-ALPHA).
IMGT_PAN_CLASS_HALF1_SLOTS: tuple[str, ...] = (
    ("1.2", "1.1")
    + tuple(str(i) for i in range(1, 8))
    + ("7A",)
    + tuple(str(i) for i in range(8, 20))
    + ("19A", "19B", "19C")
    + tuple(str(i) for i in range(20, 55))
    + ("54A", "54B", "54C")
    + tuple(str(i) for i in range(55, 91))
    + ("92A",)
)

#: Shared 3D-homologous Pan-Class Half 2 slot schema (102 slots, shared by G-ALPHA2 and G-BETA).
IMGT_PAN_CLASS_HALF2_SLOTS: tuple[str, ...] = IMGT_G_ALPHA2_SLOTS

#: Unified Pan-Class fixed-index IMGT groove frame (202 slots: 100 HALF1 + 102 HALF2).
IMGT_PAN_CLASS_GROOVE_SLOTS: tuple[tuple[str, str], ...] = tuple(("HALF1", s) for s in IMGT_PAN_CLASS_HALF1_SLOTS) + tuple(
    ("HALF2", s) for s in IMGT_PAN_CLASS_HALF2_SLOTS
)

#: 34 Class I NetMHCpan-4.1 peptide-contacting positions expressed in 3D-homologous
#: pan-class IMGT coordinates (18 HALF1 + 16 HALF2).
IMGT_PAN_CLASS_NETMHCPAN_SLOTS: tuple[tuple[str, str], ...] = (
    ("HALF1", "7"),
    ("HALF1", "9"),
    ("HALF1", "24"),
    ("HALF1", "45"),
    ("HALF1", "59"),
    ("HALF1", "62"),
    ("HALF1", "63"),
    ("HALF1", "66"),
    ("HALF1", "67"),
    ("HALF1", "69"),
    ("HALF1", "70"),
    ("HALF1", "73"),
    ("HALF1", "74"),
    ("HALF1", "76"),
    ("HALF1", "77"),
    ("HALF1", "80"),
    ("HALF1", "81"),
    ("HALF1", "84"),
    ("HALF2", "5"),
    ("HALF2", "7"),
    ("HALF2", "9"),
    ("HALF2", "24"),
    ("HALF2", "26"),
    ("HALF2", "28"),
    ("HALF2", "55"),
    ("HALF2", "59"),
    ("HALF2", "61A"),
    ("HALF2", "63"),
    ("HALF2", "67"),
    ("HALF2", "69"),
    ("HALF2", "70"),
    ("HALF2", "73"),
    ("HALF2", "77"),
    ("HALF2", "81"),
)

#: 34 Class II NetMHCIIpan-3.0/4.0 peptide-contacting positions (Karosiene et al. 2013:
#: 15 DRA alpha mature [9,11,22,24,31,52,53,58,59,61,65,66,68,72,73] + 19 DRB1 beta mature
#: [9,11,13,26,28,30,47,57,67,70,71,74,77,78,81,85,86,89,90]) expressed in 3D-homologous
#: pan-class IMGT coordinates (15 HALF1 + 19 HALF2). Exact 22/34 overlap with Class I
#: NetMHCpan-34 (11 HALF1 + 11 HALF2), plus 12 Class II-specific P1/P4/P7/P9 slots.
IMGT_PAN_CLASS_NETMHCIIPAN_SLOTS: tuple[tuple[str, str], ...] = (
    ("HALF1", "7"),
    ("HALF1", "9"),
    ("HALF1", "24"),
    ("HALF1", "26"),
    ("HALF1", "33"),
    ("HALF1", "56"),
    ("HALF1", "61"),
    ("HALF1", "66"),
    ("HALF1", "67"),
    ("HALF1", "69"),
    ("HALF1", "73"),
    ("HALF1", "74"),
    ("HALF1", "76"),
    ("HALF1", "80"),
    ("HALF1", "81"),
    ("HALF2", "5"),
    ("HALF2", "7"),
    ("HALF2", "9"),
    ("HALF2", "22"),
    ("HALF2", "24"),
    ("HALF2", "26"),
    ("HALF2", "45"),
    ("HALF2", "55"),
    ("HALF2", "63"),
    ("HALF2", "66"),
    ("HALF2", "67"),
    ("HALF2", "70"),
    ("HALF2", "72A"),
    ("HALF2", "73"),
    ("HALF2", "76"),
    ("HALF2", "80"),
    ("HALF2", "81"),
    ("HALF2", "84"),
    ("HALF2", "85"),
)

#: 46-slot exact 3D-aligned union of Class I NetMHCpan-34 and Class II NetMHCIIpan-34
#: pocket contacts (22 shared + 12 Class I-only + 12 Class II-only), ordered along
#: ``IMGT_PAN_CLASS_GROOVE_SLOTS`` (22 HALF1 + 24 HALF2).
IMGT_PAN_CLASS_POCKET_SLOTS: tuple[tuple[str, str], ...] = tuple(
    slot for slot in IMGT_PAN_CLASS_GROOVE_SLOTS if slot in (set(IMGT_PAN_CLASS_NETMHCPAN_SLOTS) | set(IMGT_PAN_CLASS_NETMHCIIPAN_SLOTS))
)

#: Class-specific IMGT G-DOMAIN insertion hotspot schemas.
#: Each :class:`InsertionHotspot` sits strictly between ``anchor_imgt_pos`` and
#: ``right_flank_imgt_pos`` and ONLY absorbs insertions with
#: ``after_imgt_pos == anchor_imgt_pos`` so fixed frames never transpose residues
#: across canonical columns.
G_ALPHA1_INSERTION_HOTSPOTS: tuple[InsertionHotspot, ...] = (
    InsertionHotspot(
        name="ab_turn",
        structural_element="AB-TURN",
        anchor_imgt_pos="19",
        right_flank_imgt_pos="20",
        slots=("19A", "19B", "19C"),
        fill_mode="left",
        lineage="Class I G-ALPHA1 AB-turn loop insertions across teleosts, sharks, birds, and non-classical mammals (Lefranc et al. 2005)",
    ),
    InsertionHotspot(
        name="n_helix_hinge",
        structural_element="HELIX",
        anchor_imgt_pos="54",
        right_flank_imgt_pos="55",
        slots=("54A", "54B", "54C"),
        fill_mode="left",
        lineage=(
            "3_10-to-alpha1 helix hinge insertions after IMGT 54 in marsupial Class I "
            "(Trvu-UB*01:01 VGE), bony fish, and reptiles (Lefranc et al. 2005)"
        ),
    ),
)

G_ALPHA2_INSERTION_HOTSPOTS: tuple[InsertionHotspot, ...] = (
    InsertionHotspot(
        name="ab_turn",
        structural_element="AB-TURN",
        anchor_imgt_pos="13",
        right_flank_imgt_pos="14",
        slots=("13A", "13B", "13C"),
        fill_mode="left",
        lineage="Class I G-ALPHA2 AB-turn loop insertions downstream of conserved Cys11 across ray-finned fish, sharks, amphibians, and reptiles",
    ),
    InsertionHotspot(
        name="cd_turn",
        structural_element="CD-TURN",
        anchor_imgt_pos="39",
        right_flank_imgt_pos="42",
        slots=("40", "40A", "40B", "40C", "41"),
        fill_mode="cd_turn",
        lineage=(
            "Class I G-ALPHA2 CD-turn insertions filling canonical IMGT 40/41 before central 40A-40C in non-mammalian vertebrates "
            "(gapped in human HLA, Lefranc et al. 2005)"
        ),
    ),
    InsertionHotspot(
        name="helix_kink",
        structural_element="HELIX",
        anchor_imgt_pos="61A",
        right_flank_imgt_pos="62",
        slots=("61B", "61C"),
        fill_mode="left",
        lineage=(
            "Class I G-ALPHA2 alpha2-helix kink insertions downstream of IMGT 61A in non-mammalian vertebrates (homologous to Class II G-BETA 61B)"
        ),
    ),
)

G_ALPHA_INSERTION_HOTSPOTS: tuple[InsertionHotspot, ...] = (
    InsertionHotspot(
        name="a_strand_bulge",
        structural_element="A-STRAND",
        anchor_imgt_pos="7",
        right_flank_imgt_pos="8",
        slots=("7A",),
        fill_mode="left",
        lineage="Class II G-ALPHA A-strand 1-aa bulge in HLA-DQA1, H2-Aa, RT1-Ba, and teleost DAA (PDB 1JK8/1S9V, Lefranc et al. 2005)",
    ),
    InsertionHotspot(
        name="cd_turn",
        structural_element="CD-TURN",
        anchor_imgt_pos="39",
        right_flank_imgt_pos="42",
        slots=("40", "41"),
        fill_mode="left",
        lineage=(
            "Class II G-ALPHA CD-turn insertions occupying IMGT 40-41 in teleost/amphibian DAA (structurally gapped in human DRA/DQA/DPA, PDB 1DLH)"
        ),
    ),
    InsertionHotspot(
        name="n_helix_hinge",
        structural_element="HELIX",
        anchor_imgt_pos="54",
        right_flank_imgt_pos="55",
        slots=("54A", "54B"),
        fill_mode="left",
        lineage="Class II G-ALPHA N-helix hinge insertions in HLA-DQA/DOA/DMA and teleost DAA (Lefranc et al. 2005)",
    ),
)

G_BETA_INSERTION_HOTSPOTS: tuple[InsertionHotspot, ...] = (
    InsertionHotspot(
        name="ab_turn",
        structural_element="AB-TURN",
        anchor_imgt_pos="13",
        right_flank_imgt_pos="14",
        slots=("13A", "13B", "13C"),
        fill_mode="left",
        lineage="Class II G-BETA AB-turn loop insertions downstream of conserved Cys11 across ray-finned fish, sharks, amphibians, and reptiles",
    ),
    InsertionHotspot(
        name="cd_turn",
        structural_element="CD-TURN",
        anchor_imgt_pos="39",
        right_flank_imgt_pos="42",
        slots=("40", "40A", "40B", "40C", "41"),
        fill_mode="cd_turn",
        lineage=(
            "Class II G-BETA CD-turn insertions filling canonical IMGT 40/41 before central 40A-40C in non-mammalian vertebrates "
            "(gapped in human HLA, Lefranc et al. 2005)"
        ),
    ),
    InsertionHotspot(
        name="helix_kink",
        structural_element="HELIX",
        anchor_imgt_pos="61B",
        right_flank_imgt_pos="62",
        slots=("61C",),
        fill_mode="left",
        lineage=("Class II G-BETA beta1-helix kink insertions downstream of canonical IMGT 61B in non-mammalian DAB/DBB and rare HLA-DRB variants"),
    ),
)


def _load_profile_file(path: Path) -> dict[str, Any]:
    """Load one shipped multi-clade IMGT structural PSSM + loop gap profile JSON."""
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    alphabet = raw.get("alphabet", "ARNDCQEGHILKMFPSTWYV")
    domains: dict[str, dict[str, Any]] = {}
    for dom_name, dom_data in raw.get("domains", {}).items():
        pssm_dicts = [{aa: float(val) for aa, val in zip(alphabet, row)} for row in dom_data["pssm"]]
        domains[dom_name] = {
            "length": int(dom_data["length"]),
            "consensus": str(dom_data["consensus"]),
            "gap_open_del": [float(x) for x in dom_data["gap_open_del"]],
            "gap_open_ins": [float(x) for x in dom_data["gap_open_ins"]],
            "pssm": pssm_dicts,
        }
    return {
        "version": str(raw.get("version", "1")),
        "alphabet": alphabet,
        "gap_ext": float(raw.get("gap_ext", -1.0)),
        "end_gap": float(raw.get("end_gap", -0.5)),
        "training_source": raw.get("training_source", {}),
        "clade_representatives": raw.get("clade_representatives", {}),
        "representative_ids": raw.get("representative_ids", {}),
        "domains": domains,
    }


@lru_cache(maxsize=1)
def load_class_i_profile() -> dict[str, Any]:
    """Load the shipped multi-clade Class I IMGT structural PSSM + loop gap profile."""
    return _load_profile_file(IMGT_CLASS_I_PROFILE_PATH)


@lru_cache(maxsize=1)
def load_class_ii_profile() -> dict[str, Any]:
    """Load the shipped multi-clade Class II IMGT structural PSSM + loop gap profile."""
    return _load_profile_file(IMGT_CLASS_II_PROFILE_PATH)


def _fill_hotspot_slots(seq: str, slots: tuple[str, ...], mode: str) -> tuple[dict[str, str], bool]:
    """Fill a fixed hotspot slot tuple with *seq*, returning (slot_map, overflowed)."""
    cleaned = "".join(ch if ch in CANONICAL_AA20 else "X" for ch in seq)
    k = len(slots)
    n = len(cleaned)
    out = {s: "-" for s in slots}
    if n == 0:
        return out, False
    overflow = n > k
    used = cleaned[:k]
    m = len(used)

    if mode == "left" or m == k:
        for idx, ch in enumerate(used):
            out[slots[idx]] = ch
        return out, overflow

    if mode == "cd_turn" and k == 5:
        # ("40", "40A", "40B", "40C", "41"): 1->40, 2->40+41, 3->40+40B+41, 4->40+40A+40C+41
        patterns = {
            1: (0,),
            2: (0, 4),
            3: (0, 2, 4),
            4: (0, 1, 3, 4),
            5: (0, 1, 2, 3, 4),
        }
        for pos_idx, ch in zip(patterns[m], used):
            out[slots[pos_idx]] = ch
        return out, overflow

    # Default "center_out" (e.g. 3 slots: 1->mid, 2->flanks, 3->all)
    if m == 1:
        out[slots[k // 2]] = used[0]
    elif m == 2:
        out[slots[0]] = used[0]
        out[slots[-1]] = used[1]
    else:
        for idx, ch in enumerate(used):
            out[slots[idx]] = ch
    return out, overflow


@dataclass(frozen=True)
class AlignedResidue:
    """One actual residue in a variable-length groove with its IMGT coordinates.

    For canonical match residues, ``mature_pos`` is the 1-indexed canonical
    mature position and ``is_insertion`` is ``False``.  For internal insertions,
    ``mature_pos`` is the preceding anchor mature position and ``is_insertion``
    is ``True``.  For N- or C-terminal overhangs outside the canonical groove
    boundaries, ``mature_pos`` is ``None``, ``is_insertion`` is ``True``, and
    ``structural_element`` is ``"N-OVERHANG"`` or ``"C-OVERHANG"``.
    """

    aa: str
    domain: str
    query_index: int
    mature_pos: Optional[int]
    imgt_pos: str
    is_insertion: bool = False
    insertion_rank: int = 0
    structural_element: str = "unknown"


@dataclass(frozen=True)
class GrooveInsertion:
    """Structured annotation for an internal insertion relative to the IMGT reference."""

    domain: str
    after_mature_pos: int
    after_imgt_pos: str
    sequence: str
    query_start: int
    imgt_codes: tuple[str, ...]
    structural_element: str

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class GrooveDeletion:
    """Structured annotation for an internal deletion (gapped canonical position)."""

    domain: str
    mature_pos: int
    imgt_pos: str
    structural_element: str


@dataclass(frozen=True)
class DomainAlignment:
    """Alignment of one MHC groove domain (G-ALPHA1, G-ALPHA2, G-ALPHA, or G-BETA) to IMGT coordinates."""

    domain: str
    query_seq: str
    ref_seq: str
    aligned_ref: str
    aligned_query: str
    mature_to_aa: dict[int, str]
    imgt_to_aa: dict[str, str]
    mature_to_query_index: dict[int, int]
    query_index_to_mature: dict[int, int]
    query_index_to_imgt: dict[int, str]
    insertions: tuple[GrooveInsertion, ...] = ()
    deletions: tuple[GrooveDeletion, ...] = ()
    n_term_missing: int = 0
    c_term_missing: int = 0
    n_term_overhang: str = ""
    c_term_overhang: str = ""
    score: float = 0.0
    identity: float = 0.0
    residues: tuple[AlignedResidue, ...] = ()
    a3m_seq: str = ""
    fixed_seq: str = ""
    fixed_slots: tuple[str, ...] = ()
    fixed_slot_to_aa: dict[str, str] = field(default_factory=dict)
    overflow_insertions: tuple[GrooveInsertion, ...] = ()
    pan_class_fixed_seq: str = ""
    pan_class_fixed_slots: tuple[str, ...] = ()
    pan_class_slot_to_aa: dict[str, str] = field(default_factory=dict)
    pan_class_overflow_insertions: tuple[GrooveInsertion, ...] = ()
    pan_class_insertion_counts: tuple[int, ...] = ()

    @property
    def has_insertions(self) -> bool:
        return bool(self.insertions)

    @property
    def has_deletions(self) -> bool:
        return bool(self.deletions)

    @property
    def has_indels(self) -> bool:
        return bool(self.insertions or self.deletions)

    @property
    def has_overflow_insertions(self) -> bool:
        return bool(self.overflow_insertions)

    @property
    def has_pan_class_overflow_insertions(self) -> bool:
        return bool(self.pan_class_overflow_insertions)

    @property
    def is_truncated(self) -> bool:
        return self.n_term_missing > 0 or self.c_term_missing > 0

    @property
    def canonical_start_mature_pos(self) -> Optional[int]:
        """First canonical 1-indexed mature position aligned to a query residue."""
        if not self.mature_to_query_index:
            return None
        return min(self.mature_to_query_index)


@dataclass(frozen=True)
class GrooveAlignment:
    """Combined Class I (G-ALPHA1/G-ALPHA2) or Class II (G-ALPHA/G-BETA) IMGT alignment.

    Provides all sequence granularities from a single exact profile alignment:
      - Variable-length lossless: ``residues`` (:class:`AlignedResidue` tokens) and ``a3m_seq``
      - Class-specific fixed-index IMGT frame (zero dead columns): ``fixed_seq``
        (198 slots for Class I, 185 slots for Class II, 22-token alphabet with ``'-'`` vs ``'X'``)
      - Pan-class shared 3D-homologous IMGT frame: ``pan_class_fixed_seq``
        (202 slots over ``IMGT_PAN_CLASS_GROOVE_SLOTS``, shared by Class I and Class II)
      - Compact pocket / hotspot extraction: :meth:`extract_positions`,
        :meth:`extract_imgt_slots`, and :meth:`extract_pan_class_slots`
    """

    alpha1: Optional[DomainAlignment] = None
    alpha2: Optional[DomainAlignment] = None
    mhc_class: Optional[str] = None

    @property
    def resolved_mhc_class(self) -> str:
        """Return ``"I"`` or ``"II"`` inferred from explicit ``mhc_class`` or aligned domains."""
        if self.mhc_class in ("I", "II"):
            return self.mhc_class
        if (self.alpha1 is not None and self.alpha1.domain == "G-ALPHA") or (self.alpha2 is not None and self.alpha2.domain == "G-BETA"):
            return "II"
        return "I"

    @property
    def half1(self) -> Optional[DomainAlignment]:
        """First groove half (``G-ALPHA1`` for Class I, ``G-ALPHA`` for Class II)."""
        return self.alpha1

    @property
    def half2(self) -> Optional[DomainAlignment]:
        """Second groove half (``G-ALPHA2`` for Class I, ``G-BETA`` for Class II)."""
        return self.alpha2

    @property
    def g_alpha(self) -> Optional[DomainAlignment]:
        """Class II ``G-ALPHA`` half (or ``None`` for Class I)."""
        return self.alpha1 if self.resolved_mhc_class == "II" else None

    @property
    def g_beta(self) -> Optional[DomainAlignment]:
        """Class II ``G-BETA`` half (or ``None`` for Class I)."""
        return self.alpha2 if self.resolved_mhc_class == "II" else None

    @property
    def insertions(self) -> tuple[GrooveInsertion, ...]:
        parts: list[GrooveInsertion] = []
        if self.alpha1 is not None:
            parts.extend(self.alpha1.insertions)
        if self.alpha2 is not None:
            parts.extend(self.alpha2.insertions)
        return tuple(parts)

    @property
    def deletions(self) -> tuple[GrooveDeletion, ...]:
        parts: list[GrooveDeletion] = []
        if self.alpha1 is not None:
            parts.extend(self.alpha1.deletions)
        if self.alpha2 is not None:
            parts.extend(self.alpha2.deletions)
        return tuple(parts)

    @property
    def overflow_insertions(self) -> tuple[GrooveInsertion, ...]:
        parts: list[GrooveInsertion] = []
        if self.alpha1 is not None:
            parts.extend(self.alpha1.overflow_insertions)
        if self.alpha2 is not None:
            parts.extend(self.alpha2.overflow_insertions)
        return tuple(parts)

    @property
    def pan_class_overflow_insertions(self) -> tuple[GrooveInsertion, ...]:
        """Internal insertions that did not fit into the 202-slot pan-class IMGT frame."""
        parts: list[GrooveInsertion] = []
        if self.alpha1 is not None:
            parts.extend(self.alpha1.pan_class_overflow_insertions)
        if self.alpha2 is not None:
            parts.extend(self.alpha2.pan_class_overflow_insertions)
        return tuple(parts)

    @property
    def has_pan_class_overflow_insertions(self) -> bool:
        return bool(self.pan_class_overflow_insertions)

    @property
    def pan_class_insertion_counts(self) -> tuple[int, ...]:
        """202-length tuple giving unrepresented overflow insertion length after each pan-class slot."""
        c1 = (
            self.alpha1.pan_class_insertion_counts
            if (self.alpha1 is not None and self.alpha1.pan_class_insertion_counts)
            else (0,) * len(IMGT_PAN_CLASS_HALF1_SLOTS)
        )
        c2 = (
            self.alpha2.pan_class_insertion_counts
            if (self.alpha2 is not None and self.alpha2.pan_class_insertion_counts)
            else (0,) * len(IMGT_PAN_CLASS_HALF2_SLOTS)
        )
        return c1 + c2

    @property
    def has_insertions(self) -> bool:
        return bool(self.insertions)

    @property
    def has_deletions(self) -> bool:
        return bool(self.deletions)

    @property
    def has_indels(self) -> bool:
        return bool(self.insertions or self.deletions)

    @property
    def has_overflow_insertions(self) -> bool:
        return bool(self.overflow_insertions)

    @property
    def is_truncated(self) -> bool:
        return (self.alpha1 is not None and self.alpha1.is_truncated) or (self.alpha2 is not None and self.alpha2.is_truncated)

    @property
    def residues(self) -> tuple[AlignedResidue, ...]:
        """Lossless variable-length per-residue tokens across both groove halves."""
        parts: list[AlignedResidue] = []
        if self.alpha1 is not None:
            parts.extend(self.alpha1.residues)
        if self.alpha2 is not None:
            parts.extend(self.alpha2.residues)
        return tuple(parts)

    @property
    def a3m_seq(self) -> str:
        """Lossless variable-length IMGT A3M string.

        Uppercase ``A-Z`` = canonical match column (including observed ``'X'``),
        lowercase ``a-z`` = internal insertion or N/C-terminal overhang,
        ``'-'`` = unoccupied match column. Each domain's ``deletions``,
        ``n_term_missing``, and ``c_term_missing`` distinguish internal
        deletions from terminal truncations. An absent domain is all gaps.
        """
        is_c2 = self.resolved_mhc_class == "II"
        len1 = G_ALPHA_SIZE if is_c2 else G_ALPHA1_SIZE
        len2 = G_BETA_SIZE if is_c2 else G_ALPHA2_SIZE
        s1 = self.alpha1.a3m_seq if self.alpha1 is not None else ("-" * len1)
        s2 = self.alpha2.a3m_seq if self.alpha2 is not None else ("-" * len2)
        return s1 + s2

    @property
    def fixed_slots(self) -> tuple[tuple[str, str], ...]:
        """Ordered ``(domain, imgt_slot)`` coordinates corresponding to :attr:`fixed_seq`."""
        return IMGT_CLASS_II_GROOVE_SLOTS if self.resolved_mhc_class == "II" else IMGT_CLASS_I_GROOVE_SLOTS

    @property
    def fixed_seq(self) -> str:
        """Class-specific fixed-index IMGT string (198 chars for Class I, 185 for Class II)."""
        is_c2 = self.resolved_mhc_class == "II"
        n1 = len(IMGT_G_ALPHA_SLOTS) if is_c2 else len(IMGT_G_ALPHA1_SLOTS)
        n2 = len(IMGT_G_BETA_SLOTS) if is_c2 else len(IMGT_G_ALPHA2_SLOTS)
        s1 = self.alpha1.fixed_seq if self.alpha1 is not None else ("X" * n1)
        s2 = self.alpha2.fixed_seq if self.alpha2 is not None else ("X" * n2)
        return s1 + s2

    @property
    def pan_class_fixed_slots(self) -> tuple[tuple[str, str], ...]:
        """Ordered ``(half, imgt_slot)`` coordinates for :attr:`pan_class_fixed_seq` (202 slots)."""
        return IMGT_PAN_CLASS_GROOVE_SLOTS

    @property
    def pan_class_fixed_seq(self) -> str:
        """Shared 3D-homologous pan-class IMGT string (202 chars, 22-token alphabet with '-' vs 'X')."""
        s1 = self.alpha1.pan_class_fixed_seq if self.alpha1 is not None else ("X" * len(IMGT_PAN_CLASS_HALF1_SLOTS))
        s2 = self.alpha2.pan_class_fixed_seq if self.alpha2 is not None else ("X" * len(IMGT_PAN_CLASS_HALF2_SLOTS))
        return s1 + s2

    @property
    def mature_to_aa(self) -> dict[int, str]:
        merged: dict[int, str] = {}
        if self.alpha1 is not None:
            merged.update(self.alpha1.mature_to_aa)
        if self.alpha2 is not None:
            merged.update(self.alpha2.mature_to_aa)
        return merged

    @property
    def imgt_to_aa(self) -> dict[tuple[str, str], str]:
        merged: dict[tuple[str, str], str] = {}
        if self.alpha1 is not None:
            for pos, aa in self.alpha1.imgt_to_aa.items():
                merged[(self.alpha1.domain, pos)] = aa
        if self.alpha2 is not None:
            for pos, aa in self.alpha2.imgt_to_aa.items():
                merged[(self.alpha2.domain, pos)] = aa
        return merged

    def residue_at(self, mature_pos: int) -> Optional[str]:
        """Return the aligned amino acid at canonical 1-indexed *mature_pos*.

        For Class I, positions 1–90 query ``G-ALPHA1`` and >90 query ``G-ALPHA2``.
        For Class II (where alpha and beta are separate chains with independent
        mature numbering), use :meth:`residue_at_imgt` or query ``.g_alpha`` /
        ``.g_beta`` directly when both chains are present.
        """
        if self.resolved_mhc_class == "II":
            if self.alpha1 is not None and self.alpha2 is None:
                return self.alpha1.mature_to_aa.get(mature_pos)
            if self.alpha2 is not None and self.alpha1 is None:
                return self.alpha2.mature_to_aa.get(mature_pos)
            return self.alpha1.mature_to_aa.get(mature_pos) if self.alpha1 is not None else None
        if mature_pos <= G_ALPHA1_SIZE:
            return self.alpha1.mature_to_aa.get(mature_pos) if self.alpha1 is not None else None
        return self.alpha2.mature_to_aa.get(mature_pos) if self.alpha2 is not None else None

    def residue_at_imgt(self, domain: str, imgt_pos: str) -> Optional[str]:
        """Return the aligned amino acid at IMGT ``(domain, imgt_pos)``.

        Supports Class I (``G-ALPHA1``, ``G-ALPHA2``), Class II (``G-ALPHA``,
        ``G-BETA``), and pan-class half names (``HALF1``, ``HALF2``).
        Returns ``"-"`` for structural gaps / unoccupied insertion slots, or
        ``None`` when the position is missing / truncated.
        """
        token = str(domain or "").strip().upper().replace("_", "-")
        if token in ("HALF1", "H1"):
            dom = self.alpha1
            if dom is None:
                return None
            ch = dom.pan_class_slot_to_aa.get(imgt_pos)
            if ch is not None:
                return None if ch == "X" else ch
            return dom.imgt_to_aa.get(imgt_pos)
        if token in ("HALF2", "H2"):
            dom = self.alpha2
            if dom is None:
                return None
            ch = dom.pan_class_slot_to_aa.get(imgt_pos)
            if ch is not None:
                return None if ch == "X" else ch
            return dom.imgt_to_aa.get(imgt_pos)

        norm = _normalize_groove_domain(domain)
        dom = self.alpha1 if norm in ("G-ALPHA1", "G-ALPHA") else self.alpha2
        if dom is None or dom.domain != norm:
            return None
        if imgt_pos in dom.imgt_to_aa:
            return dom.imgt_to_aa[imgt_pos]
        if imgt_pos in dom.fixed_slot_to_aa:
            ch = dom.fixed_slot_to_aa[imgt_pos]
            return None if ch == "X" else ch
        return None

    def extract_positions(
        self,
        mature_positions: Sequence[int],
        *,
        missing: str = "X",
        deleted: str = "-",
    ) -> str:
        """Extract a fixed-position sequence (e.g. pseudosequence) at canonical mature positions.

        Internal deletions render as *deleted* (default ``"-"``) and
        truncated/unavailable positions render as *missing* (default ``"X"``).
        """
        chars: list[str] = []
        for pos in mature_positions:
            aa = self.residue_at(pos)
            if aa is None:
                chars.append(missing)
            elif aa == "-":
                chars.append(deleted)
            else:
                chars.append(aa)
        return "".join(chars)

    def extract_imgt_slots(
        self,
        slots: Sequence[tuple[str, str]],
        *,
        missing: str = "X",
        deleted: str = "-",
    ) -> str:
        """Extract a fixed-position sequence at arbitrary ``(domain_or_half, imgt_slot)`` coordinates.

        Accepts Class I (``G-ALPHA1``/``G-ALPHA2``), Class II (``G-ALPHA``/``G-BETA``),
        and pan-class (``HALF1``/``HALF2``) slot coordinates, keeping structural
        gaps (``deleted='-'``) distinct from missing/truncated positions
        (``missing='X'``).
        """
        chars: list[str] = []
        for domain, slot in slots:
            token = str(domain or "").strip().upper().replace("_", "-")
            if token in ("HALF1", "H1"):
                dom = self.alpha1
                slot_lookup = dom.pan_class_slot_to_aa if dom is not None else {}
            elif token in ("HALF2", "H2"):
                dom = self.alpha2
                slot_lookup = dom.pan_class_slot_to_aa if dom is not None else {}
            else:
                norm = _normalize_groove_domain(domain)
                dom = self.alpha1 if norm in ("G-ALPHA1", "G-ALPHA") else self.alpha2
                if dom is not None and dom.domain != norm:
                    dom = None
                slot_lookup = dom.fixed_slot_to_aa if dom is not None else {}

            if dom is None:
                chars.append(missing)
                continue
            if slot in slot_lookup:
                ch = slot_lookup[slot]
                if ch == "-":
                    chars.append(deleted)
                elif ch == "X":
                    chars.append(missing)
                else:
                    chars.append(ch)
            elif slot in dom.imgt_to_aa:
                ch = dom.imgt_to_aa[slot]
                chars.append(deleted if ch == "-" else ch)
            else:
                chars.append(missing)
        return "".join(chars)

    def extract_pan_class_slots(
        self,
        slots: Sequence[tuple[str, str]],
        *,
        missing: str = "X",
        deleted: str = "-",
    ) -> str:
        """Extract a fixed-position sequence from the shared 202-slot pan-class IMGT frame."""
        return self.extract_imgt_slots(slots, missing=missing, deleted=deleted)

    def extract_pan_class_pocket_seq(
        self,
        *,
        missing: str = "X",
        deleted: str = "-",
    ) -> str:
        """Extract the 46-slot 3D-aligned union of NetMHCpan-34 and NetMHCIIpan-34 pocket contacts."""
        return self.extract_imgt_slots(IMGT_PAN_CLASS_POCKET_SLOTS, missing=missing, deleted=deleted)

    def to_pan_class_features(self) -> dict[str, Any]:
        """Export a predictor-ready feature dictionary on the 202-slot pan-class IMGT grid.

        Includes ``seq`` (202 chars, 22-token alphabet with ``'-'`` vs ``'X'``),
        ``pocket_seq`` (46 chars covering 100% of NetMHCpan-34 and NetMHCIIpan-34
        contacts), ``insertion_counts`` (202-tuple of per-slot unrepresented overflow
        insertion lengths), ``overflow_insertions``, and terminal overhangs.
        """
        return {
            "mhc_class": self.resolved_mhc_class,
            "slots": IMGT_PAN_CLASS_GROOVE_SLOTS,
            "seq": self.pan_class_fixed_seq,
            "pocket_slots": IMGT_PAN_CLASS_POCKET_SLOTS,
            "pocket_seq": self.extract_pan_class_pocket_seq(),
            "insertion_counts": self.pan_class_insertion_counts,
            "overflow_insertions": self.pan_class_overflow_insertions,
            "n_term_overhang_half1": self.alpha1.n_term_overhang if self.alpha1 is not None else "",
            "c_term_overhang_half1": self.alpha1.c_term_overhang if self.alpha1 is not None else "",
            "n_term_overhang_half2": self.alpha2.n_term_overhang if self.alpha2 is not None else "",
            "c_term_overhang_half2": self.alpha2.c_term_overhang if self.alpha2 is not None else "",
        }

    def summary(self) -> str:
        """Return a concise diagnostic summary of indels and truncations."""
        parts: list[str] = []
        labels = ("alpha", "beta") if self.resolved_mhc_class == "II" else ("alpha1", "alpha2")
        for label, dom in ((labels[0], self.alpha1), (labels[1], self.alpha2)):
            if dom is None:
                continue
            notes: list[str] = []
            if dom.n_term_missing:
                notes.append(f"N-trunc {dom.n_term_missing}aa")
            if dom.c_term_missing:
                notes.append(f"C-trunc {dom.c_term_missing}aa")
            for ins in dom.insertions:
                notes.append(f"+{ins.length}aa ins@{ins.after_mature_pos}({ins.sequence})")
            if dom.deletions:
                del_pos = ",".join(str(d.mature_pos) for d in dom.deletions)
                notes.append(f"-{len(dom.deletions)}aa del@{del_pos}")
            if notes:
                parts.append(f"{label}: " + ", ".join(notes))
        return "; ".join(parts) if parts else "canonical"


def _normalize_groove_domain(domain: str) -> str:
    token = str(domain or "").strip().upper().replace("_", "-")
    if token in ("G-ALPHA1", "ALPHA1", "GROOVE1", "A1", "I-ALPHA1", "CLASS-I-ALPHA1"):
        return "G-ALPHA1"
    if token in ("G-ALPHA2", "ALPHA2", "GROOVE2", "A2", "I-ALPHA2", "CLASS-I-ALPHA2"):
        return "G-ALPHA2"
    if token in ("G-ALPHA", "ALPHA", "II-ALPHA", "CLASS-II-ALPHA", "C2A"):
        return "G-ALPHA"
    if token in ("G-BETA", "BETA", "BETA1", "II-BETA", "CLASS-II-BETA", "C2B"):
        return "G-BETA"
    raise ValueError(f"Unsupported IMGT groove domain: {domain!r} (expected 'G-ALPHA1', 'G-ALPHA2', 'G-ALPHA', or 'G-BETA')")


def _normalize_class_i_domain(domain: str) -> str:
    norm = _normalize_groove_domain(domain)
    if norm not in ("G-ALPHA1", "G-ALPHA2"):
        raise ValueError(f"Unsupported Class I groove domain: {domain!r} (expected 'G-ALPHA1' or 'G-ALPHA2')")
    return norm


def _disulfide_cys_anchor_bonus(
    qry: str,
    ref_c1_idx: int,
    ref_c2_idx: int,
    bonus: float,
) -> Optional[dict[tuple[int, int], float]]:
    """Return a gentle alignment bonus for the conserved Cys11–Cys74 disulfide pair."""
    cys = [idx for idx, aa in enumerate(qry) if aa == "C"]
    if len(cys) < 2:
        return None
    target_sep = ref_c2_idx - ref_c1_idx
    best_pair: Optional[tuple[int, int]] = None
    best_dist = 10**9
    for i, c1 in enumerate(cys):
        for c2 in cys[i + 1 :]:
            sep = c2 - c1
            if 50 <= sep <= 76:
                dist = abs(sep - target_sep) * 2 + abs(c1 - ref_c1_idx)
                if dist < best_dist:
                    best_dist = dist
                    best_pair = (c1, c2)
    if best_pair is None:
        return None
    return {
        (ref_c1_idx, best_pair[0]): bonus,
        (ref_c2_idx, best_pair[1]): bonus,
    }


def g_alpha2_cys_anchor_bonus(qry: str) -> Optional[dict[tuple[int, int], float]]:
    """Return a gentle alignment bonus for the conserved G-ALPHA2 Cys11–Cys74 pair."""
    return _disulfide_cys_anchor_bonus(qry, G_ALPHA2_REF_C1_IDX, G_ALPHA2_REF_C2_IDX, G_ALPHA2_CYS_ANCHOR_BONUS)


def g_beta_cys_anchor_bonus(qry: str) -> Optional[dict[tuple[int, int], float]]:
    """Return a gentle alignment bonus for the conserved G-BETA Cys11–Cys74 pair."""
    return _disulfide_cys_anchor_bonus(qry, G_BETA_REF_C1_IDX, G_BETA_REF_C2_IDX, G_BETA_CYS_ANCHOR_BONUS)


def _align_groove_domain(ref: str, qry: str, domain: str) -> tuple[float, str, str]:
    """Run exact semi-global alignment for one Class I or Class II groove domain."""
    if domain == "G-ALPHA2":
        anchor_bonus = g_alpha2_cys_anchor_bonus(qry)
    elif domain == "G-BETA":
        anchor_bonus = g_beta_cys_anchor_bonus(qry)
    else:
        anchor_bonus = None

    profile = load_class_i_profile() if domain in ("G-ALPHA1", "G-ALPHA2") else load_class_ii_profile()
    dom_profile = profile.get("domains", {}).get(domain) if profile else None
    if dom_profile is not None:
        return _profile_gotoh(
            ref,
            dom_profile["pssm"],
            dom_profile["gap_open_del"],
            dom_profile["gap_open_ins"],
            qry,
            gap_ext=profile["gap_ext"],
            end_gap=profile["end_gap"],
            anchor_bonus=anchor_bonus,
        )
    return _semiglobal_gotoh(ref, qry, anchor_bonus=anchor_bonus)


# Cache only immutable score and traceback strings. Public coordinate maps are fresh per call.
_cached_groove_alignment = lru_cache(maxsize=4096)(_align_groove_domain)
_cached_class_i_alignment = _cached_groove_alignment


def _profile_gotoh(
    ref: str,
    pssm: list[dict[str, float]],
    go_del: list[float],
    go_ins: list[float],
    qry: str,
    *,
    gap_ext: float = -1.0,
    end_gap: float = -0.5,
    anchor_bonus: Optional[dict[tuple[int, int], float]] = None,
) -> tuple[float, str, str]:
    """Semi-global affine-gap (Gotoh) alignment against a multi-clade structural profile.

    Uses per-column log-odds match scores (*pssm*) and position-specific
    gap-open penalties (*go_del*, *go_ins*) lowered at the 6 conserved IMGT
    G-DOMAIN loop hotspots, with light end-gap penalties and no direct
    ``Ix <-> Iy`` zigzag transitions.
    """
    m = len(pssm)
    n = len(qry)
    neg = -1e9

    mat = [[neg] * (n + 1) for _ in range(m + 1)]
    ix = [[neg] * (n + 1) for _ in range(m + 1)]
    iy = [[neg] * (n + 1) for _ in range(m + 1)]
    tr_m = [[0] * (n + 1) for _ in range(m + 1)]
    tr_x = [[0] * (n + 1) for _ in range(m + 1)]
    tr_y = [[0] * (n + 1) for _ in range(m + 1)]

    mat[0][0] = 0.0
    for i in range(1, m + 1):
        ix[i][0] = end_gap * i
        tr_x[i][0] = 1 if i > 1 else 0
    for j in range(1, n + 1):
        iy[0][j] = end_gap * j
        tr_y[0][j] = 2 if j > 1 else 0

    for i in range(1, m + 1):
        col_scores = pssm[i - 1]
        del_open_i = go_del[i - 1]
        ins_open_i = go_ins[i - 1]
        m_i, m_p = mat[i], mat[i - 1]
        ix_i, ix_p = ix[i], ix[i - 1]
        iy_i, iy_p = iy[i], iy[i - 1]
        trm_i, trx_i, try_i = tr_m[i], tr_x[i], tr_y[i]
        is_last_ref = i == m

        for j in range(1, n + 1):
            qa = qry[j - 1]
            sub = col_scores.get(qa, -1.0)
            if anchor_bonus is not None:
                bonus = anchor_bonus.get((i - 1, j - 1))
                if bonus is not None:
                    sub += bonus

            v0 = m_p[j - 1]
            v1 = ix_p[j - 1]
            v2 = iy_p[j - 1]
            if v0 >= v1 and v0 >= v2:
                m_i[j] = v0 + sub
                trm_i[j] = 0
            elif v1 >= v2:
                m_i[j] = v1 + sub
                trm_i[j] = 1
            else:
                m_i[j] = v2 + sub
                trm_i[j] = 2

            go_x = end_gap if j == n else del_open_i
            ge_x = end_gap if j == n else gap_ext
            open_x = m_p[j] + go_x
            ext_x = ix_p[j] + ge_x
            if ext_x >= open_x:
                ix_i[j] = ext_x
                trx_i[j] = 1
            else:
                ix_i[j] = open_x
                trx_i[j] = 0

            go_y = end_gap if is_last_ref else ins_open_i
            ge_y = end_gap if is_last_ref else gap_ext
            open_y = m_i[j - 1] + go_y
            ext_y = iy_i[j - 1] + ge_y
            if ext_y >= open_y:
                iy_i[j] = ext_y
                try_i[j] = 2
            else:
                iy_i[j] = open_y
                try_i[j] = 0

    score, state = max(((mat[m][n], 0), (ix[m][n], 1), (iy[m][n], 2)))
    ar: list[str] = []
    aq: list[str] = []
    i, j = m, n
    while i > 0 or j > 0:
        if state == 0:
            prev = tr_m[i][j]
            ar.append(ref[i - 1])
            aq.append(qry[j - 1])
            i -= 1
            j -= 1
        elif state == 1:
            prev = tr_x[i][j]
            ar.append(ref[i - 1])
            aq.append("-")
            i -= 1
        else:
            prev = tr_y[i][j]
            ar.append("-")
            aq.append(qry[j - 1])
            j -= 1
        state = prev
    ar.reverse()
    aq.reverse()
    return score, "".join(ar), "".join(aq)


def _semiglobal_gotoh(
    ref: str,
    qry: str,
    *,
    gap_open: float = -13.0,
    gap_ext: float = -1.0,
    end_gap: float = -0.5,
    anchor_bonus: Optional[dict[tuple[int, int], float]] = None,
) -> tuple[float, str, str]:
    """Semi-global affine-gap (Gotoh) alignment with BLOSUM62 scoring.

    Uses light end-gap penalties on both leading and trailing overhangs/truncations
    and disallows direct gap-to-gap (``Ix <-> Iy``) transitions to prevent
    spurious zigzag indels.
    """
    m = len(ref)
    n = len(qry)
    neg = -1e9

    mat = [[neg] * (n + 1) for _ in range(m + 1)]
    ix = [[neg] * (n + 1) for _ in range(m + 1)]
    iy = [[neg] * (n + 1) for _ in range(m + 1)]
    tr_m = [[0] * (n + 1) for _ in range(m + 1)]
    tr_x = [[0] * (n + 1) for _ in range(m + 1)]
    tr_y = [[0] * (n + 1) for _ in range(m + 1)]

    mat[0][0] = 0.0
    for i in range(1, m + 1):
        ix[i][0] = end_gap * i
        tr_x[i][0] = 1 if i > 1 else 0
    for j in range(1, n + 1):
        iy[0][j] = end_gap * j
        tr_y[0][j] = 2 if j > 1 else 0

    for i in range(1, m + 1):
        ra = ref[i - 1]
        b62_row = _BLOSUM62.get(ra, _BLOSUM62_UNKNOWN)
        m_i, m_p = mat[i], mat[i - 1]
        ix_i, ix_p = ix[i], ix[i - 1]
        iy_i, iy_p = iy[i], iy[i - 1]
        trm_i, trx_i, try_i = tr_m[i], tr_x[i], tr_y[i]
        is_last_ref = i == m

        for j in range(1, n + 1):
            qa = qry[j - 1]
            sub = float(b62_row.get(qa, -1))
            if anchor_bonus is not None:
                bonus = anchor_bonus.get((i - 1, j - 1))
                if bonus is not None:
                    sub += bonus

            v0 = m_p[j - 1]
            v1 = ix_p[j - 1]
            v2 = iy_p[j - 1]
            if v0 >= v1 and v0 >= v2:
                m_i[j] = v0 + sub
                trm_i[j] = 0
            elif v1 >= v2:
                m_i[j] = v1 + sub
                trm_i[j] = 1
            else:
                m_i[j] = v2 + sub
                trm_i[j] = 2

            go_x = end_gap if j == n else gap_open
            ge_x = end_gap if j == n else gap_ext
            open_x = m_p[j] + go_x
            ext_x = ix_p[j] + ge_x
            if ext_x >= open_x:
                ix_i[j] = ext_x
                trx_i[j] = 1
            else:
                ix_i[j] = open_x
                trx_i[j] = 0

            go_y = end_gap if is_last_ref else gap_open
            ge_y = end_gap if is_last_ref else gap_ext
            open_y = m_i[j - 1] + go_y
            ext_y = iy_i[j - 1] + ge_y
            if ext_y >= open_y:
                iy_i[j] = ext_y
                try_i[j] = 2
            else:
                iy_i[j] = open_y
                try_i[j] = 0

    score, state = max(((mat[m][n], 0), (ix[m][n], 1), (iy[m][n], 2)))
    ar: list[str] = []
    aq: list[str] = []
    i, j = m, n
    while i > 0 or j > 0:
        if state == 0:
            prev = tr_m[i][j]
            ar.append(ref[i - 1])
            aq.append(qry[j - 1])
            i -= 1
            j -= 1
        elif state == 1:
            prev = tr_x[i][j]
            ar.append(ref[i - 1])
            aq.append("-")
            i -= 1
        else:
            prev = tr_y[i][j]
            ar.append("-")
            aq.append(qry[j - 1])
            j -= 1
        state = prev
    ar.reverse()
    aq.reverse()
    return score, "".join(ar), "".join(aq)


def _build_fixed_frame(
    domain: str,
    imgt_to_aa: dict[str, str],
    insertions: Sequence[GrooveInsertion],
) -> tuple[str, tuple[str, ...], dict[str, str], tuple[GrooveInsertion, ...]]:
    """Project an aligned domain onto its class-specific fixed-index IMGT frame with '-' vs 'X'."""
    if domain == "G-ALPHA1":
        slots = IMGT_G_ALPHA1_SLOTS
        hotspots = G_ALPHA1_INSERTION_HOTSPOTS
    elif domain == "G-ALPHA2":
        slots = IMGT_G_ALPHA2_SLOTS
        hotspots = G_ALPHA2_INSERTION_HOTSPOTS
    elif domain == "G-ALPHA":
        slots = IMGT_G_ALPHA_SLOTS
        hotspots = G_ALPHA_INSERTION_HOTSPOTS
    elif domain == "G-BETA":
        slots = IMGT_G_BETA_SLOTS
        hotspots = G_BETA_INSERTION_HOTSPOTS
    else:
        raise ValueError(f"Unsupported IMGT groove domain: {domain!r}")

    slot_map: dict[str, str] = {}
    for slot in slots:
        if slot in imgt_to_aa:
            aa = imgt_to_aa[slot]
            if aa == "-":
                slot_map[slot] = "-"
            elif aa in CANONICAL_AA20:
                slot_map[slot] = aa
            else:
                slot_map[slot] = "X"
        else:
            slot_map[slot] = "X"

    # Bucket internal insertions into hotspot slot tuples; track any overflow or
    # out-of-hotspot insertions so no structural insertion is silently lost.
    hotspot_seqs: list[list[str]] = [[] for _ in hotspots]
    hotspot_sources: list[list[GrooveInsertion]] = [[] for _ in hotspots]
    overflow: list[GrooveInsertion] = []

    for ins in insertions:
        matched_idx: Optional[int] = None
        for idx, hotspot in enumerate(hotspots):
            if ins.after_imgt_pos == hotspot.anchor_imgt_pos:
                matched_idx = idx
                break
        if matched_idx is None:
            overflow.append(ins)
        else:
            hotspot_seqs[matched_idx].append(ins.sequence)
            hotspot_sources[matched_idx].append(ins)

    for idx, hotspot in enumerate(hotspots):
        combined = "".join(hotspot_seqs[idx])
        flanks_observed = hotspot.anchor_imgt_pos in imgt_to_aa and hotspot.right_flank_imgt_pos in imgt_to_aa
        if not combined:
            fill_ch = "-" if flanks_observed else "X"
            for s in hotspot.slots:
                slot_map[s] = fill_ch
        else:
            filled, did_overflow = _fill_hotspot_slots(combined, hotspot.slots, hotspot.fill_mode)
            slot_map.update(filled)
            if did_overflow:
                overflow.extend(hotspot_sources[idx])

    fixed_seq = "".join(slot_map[s] for s in slots)
    return fixed_seq, slots, slot_map, tuple(overflow)


#: Flanking canonical backbone positions used to distinguish structural gaps ('-')
#: from missing/truncated ends ('X') for pan-class slots absent from a domain's native frame.
PAN_CLASS_EXTRA_SLOT_FLANKS: dict[str, dict[str, tuple[str, ...]]] = {
    "G-ALPHA1": {
        "1.2": ("1",),
        "1.1": ("1",),
        "7A": ("7", "8"),
        "92A": ("88",),
    },
    "G-ALPHA": {
        "16": ("15", "20"),
        "17": ("15", "20"),
        "18": ("15", "20"),
        "19": ("15", "20"),
        "19A": ("15", "20"),
        "19B": ("15", "20"),
        "19C": ("15", "20"),
        "54C": ("54", "55"),
        "57": ("56", "61"),
        "58": ("56", "61"),
        "59": ("56", "61"),
        "60": ("56", "61"),
        "92A": ("88",),
    },
    "G-ALPHA2": {},
    "G-BETA": {
        "89": ("88",),
        "90": ("88",),
        "91": ("88",),
        "92": ("88",),
    },
}


#: Pan-class Half 1 (G-ALPHA1 / G-ALPHA) insertion hotspot schemas for the 202-slot
#: shared pan-class IMGT frame. Every hotspot is strictly anchored at its exact
#: monotonic left-flanking IMGT column and uses uniform "left" fill across Class I and Class II.
PAN_CLASS_HALF1_INSERTION_HOTSPOTS: dict[str, tuple[InsertionHotspot, ...]] = {
    "G-ALPHA1": (
        InsertionHotspot(
            name="a_strand_bulge",
            structural_element="A-STRAND",
            anchor_imgt_pos="7",
            right_flank_imgt_pos="8",
            slots=("7A",),
            fill_mode="left",
            lineage="Non-mammalian Class I A-strand 1-aa bulge at IMGT 7A (shared pan-class slot with Class II G-ALPHA DQA1/H2-Aa)",
        ),
        *G_ALPHA1_INSERTION_HOTSPOTS,
    ),
    "G-ALPHA": (
        G_ALPHA_INSERTION_HOTSPOTS[0],  # 7A (a_strand_bulge)
        InsertionHotspot(
            name="ab_turn",
            structural_element="AB-TURN",
            anchor_imgt_pos="15",
            right_flank_imgt_pos="20",
            slots=("19A", "19B", "19C"),
            fill_mode="left",
            lineage="Class II G-ALPHA AB-turn insertions between IMGT 15 and 20 across the 16-19 structural deletion relative to Class I",
        ),
        G_ALPHA_INSERTION_HOTSPOTS[1],  # 40-41 (cd_turn)
        InsertionHotspot(
            name="n_helix_hinge",
            structural_element="HELIX",
            anchor_imgt_pos="54",
            right_flank_imgt_pos="55",
            slots=("54A", "54B", "54C"),
            fill_mode="left",
            lineage=(
                "Class II G-ALPHA N-helix hinge insertions (up to 3 aa in teleost DAA, plus HLA-DQA/DOA/DMA) sharing pan-class 54A-54C with Class I"
            ),
        ),
    ),
}

#: Pan-class Half 2 (G-ALPHA2 / G-BETA) insertion hotspot schemas for the 202-slot frame.
PAN_CLASS_HALF2_INSERTION_HOTSPOTS: dict[str, tuple[InsertionHotspot, ...]] = {
    "G-ALPHA2": G_ALPHA2_INSERTION_HOTSPOTS,
    "G-BETA": G_BETA_INSERTION_HOTSPOTS,
}


def validate_imgt_frame_schemas() -> None:
    """Validate structural integrity of all native and pan-class IMGT frame & hotspot schemas.

    Raises :class:`ValueError` if any hotspot slot is missing or non-contiguous,
    sits out of monotonic order relative to its anchor/flank columns, omits a
    native hotspot from its pan-class table, uses an invalid fill mode, or lacks
    documented biological lineage.
    """
    canonical_positions_by_domain: dict[str, frozenset[str]] = {
        "G-ALPHA1": frozenset(G_ALPHA1_POSITIONS),
        "G-ALPHA2": frozenset(G_ALPHA2_POSITIONS),
        "G-ALPHA": frozenset(G_ALPHA_POSITIONS),
        "G-BETA": frozenset(G_BETA_POSITIONS),
    }
    native_tables: dict[str, tuple[tuple[str, ...], tuple[InsertionHotspot, ...]]] = {
        "G-ALPHA1": (IMGT_G_ALPHA1_SLOTS, G_ALPHA1_INSERTION_HOTSPOTS),
        "G-ALPHA2": (IMGT_G_ALPHA2_SLOTS, G_ALPHA2_INSERTION_HOTSPOTS),
        "G-ALPHA": (IMGT_G_ALPHA_SLOTS, G_ALPHA_INSERTION_HOTSPOTS),
        "G-BETA": (IMGT_G_BETA_SLOTS, G_BETA_INSERTION_HOTSPOTS),
    }
    pan_tables: dict[str, tuple[tuple[str, ...], tuple[InsertionHotspot, ...]]] = {
        "G-ALPHA1": (IMGT_PAN_CLASS_HALF1_SLOTS, PAN_CLASS_HALF1_INSERTION_HOTSPOTS["G-ALPHA1"]),
        "G-ALPHA": (IMGT_PAN_CLASS_HALF1_SLOTS, PAN_CLASS_HALF1_INSERTION_HOTSPOTS["G-ALPHA"]),
        "G-ALPHA2": (IMGT_PAN_CLASS_HALF2_SLOTS, PAN_CLASS_HALF2_INSERTION_HOTSPOTS["G-ALPHA2"]),
        "G-BETA": (IMGT_PAN_CLASS_HALF2_SLOTS, PAN_CLASS_HALF2_INSERTION_HOTSPOTS["G-BETA"]),
    }

    def _check_table(label: str, domain: str, slots: tuple[str, ...], hotspots: tuple[InsertionHotspot, ...]) -> None:
        if len(slots) != len(set(slots)):
            raise ValueError(f"{label} ({domain}): duplicate slot names in slot table")
        pos_index = {slot: idx for idx, slot in enumerate(slots)}
        canonical = canonical_positions_by_domain[domain]
        seen_anchors: set[str] = set()
        seen_hotspot_slots: set[str] = set()
        for hs in hotspots:
            if not hs.name or not hs.lineage:
                raise ValueError(f"{label} ({domain}): hotspot {hs!r} must have non-empty name and lineage")
            if hs.fill_mode not in ("left", "cd_turn"):
                raise ValueError(f"{label} ({domain}): invalid fill_mode {hs.fill_mode!r} on {hs.name!r}")
            if hs.fill_mode == "cd_turn" and len(hs.slots) != 5:
                raise ValueError(f"{label} ({domain}): cd_turn fill requires 5 slots on {hs.name!r}")
            if hs.anchor_imgt_pos in seen_anchors:
                raise ValueError(f"{label} ({domain}): duplicate anchor {hs.anchor_imgt_pos!r}")
            seen_anchors.add(hs.anchor_imgt_pos)
            if hs.anchor_imgt_pos not in canonical or hs.right_flank_imgt_pos not in canonical:
                raise ValueError(f"{label} ({domain}): hotspot {hs.name!r} flanks must be canonical {domain} columns")
            if hs.anchor_imgt_pos not in pos_index or hs.right_flank_imgt_pos not in pos_index:
                raise ValueError(f"{label} ({domain}): hotspot {hs.name!r} flanks missing from slot table")
            anchor_idx = pos_index[hs.anchor_imgt_pos]
            right_idx = pos_index[hs.right_flank_imgt_pos]
            slot_indices: list[int] = []
            for s in hs.slots:
                if s not in pos_index:
                    raise ValueError(f"{label} ({domain}): hotspot slot {s!r} missing from slot table")
                if s in seen_hotspot_slots:
                    raise ValueError(f"{label} ({domain}): hotspot slot {s!r} assigned twice")
                seen_hotspot_slots.add(s)
                slot_indices.append(pos_index[s])
            if slot_indices != list(range(slot_indices[0], slot_indices[0] + len(slot_indices))):
                raise ValueError(f"{label} ({domain}): hotspot {hs.name!r} slots {hs.slots} are not contiguous")
            if not (anchor_idx < slot_indices[0] and slot_indices[-1] < right_idx):
                raise ValueError(
                    f"{label} ({domain}): hotspot {hs.name!r} slots {hs.slots} "
                    f"not strictly between {hs.anchor_imgt_pos} and {hs.right_flank_imgt_pos}"
                )
            # No canonical match column may sit between anchor_imgt_pos and slots[0] or between slots[-1] and right_flank_imgt_pos
            for between_idx in list(range(anchor_idx + 1, slot_indices[0])) + list(range(slot_indices[-1] + 1, right_idx)):
                if slots[between_idx] in canonical:
                    raise ValueError(
                        f"{label} ({domain}): canonical column {slots[between_idx]!r} trapped inside hotspot {hs.name!r} "
                        f"({hs.anchor_imgt_pos}..{hs.right_flank_imgt_pos}), which would scramble residue order"
                    )

    for dom, (slots, hotspots) in native_tables.items():
        _check_table("native", dom, slots, hotspots)
    for dom, (slots, hotspots) in pan_tables.items():
        _check_table("pan_class", dom, slots, hotspots)
        native_anchors = {hs.anchor_imgt_pos for hs in native_tables[dom][1]}
        pan_anchors = {hs.anchor_imgt_pos for hs in hotspots}
        missing_in_pan = native_anchors - pan_anchors
        if missing_in_pan:
            raise ValueError(f"pan_class ({dom}): omits native hotspot anchors {sorted(missing_in_pan)} (would double-count overflow)")
        pan_slot_set = set(slots)
        for extra_slot, flanks in PAN_CLASS_EXTRA_SLOT_FLANKS[dom].items():
            if extra_slot not in pan_slot_set:
                raise ValueError(f"PAN_CLASS_EXTRA_SLOT_FLANKS[{dom!r}] references unknown slot {extra_slot!r}")
            for flank in flanks:
                if flank not in canonical_positions_by_domain[dom]:
                    raise ValueError(f"PAN_CLASS_EXTRA_SLOT_FLANKS[{dom!r}][{extra_slot!r}] references non-canonical flank {flank!r}")


validate_imgt_frame_schemas()


def _build_pan_class_frame(
    domain: str,
    fixed_slot_to_aa: dict[str, str],
    imgt_to_aa: dict[str, str],
    insertions: Sequence[GrooveInsertion] = (),
) -> tuple[str, tuple[str, ...], dict[str, str], tuple[GrooveInsertion, ...], tuple[int, ...]]:
    """Project a Class I or Class II domain alignment onto the shared 3D pan-class IMGT frame.

    Buckets internal insertions directly into pan-class hotspot slots (including
    ``7A``, ``19A-19C``, ``40-41``, and ``54A-54C`` on ``HALF1``) and returns both
    any remaining ``pan_class_overflow_insertions`` and a per-slot
    ``pan_class_insertion_counts`` tuple giving the number of unrepresented
    inserted residues immediately following each pan-class slot.
    """
    is_half1 = domain in ("G-ALPHA1", "G-ALPHA")
    slots = IMGT_PAN_CLASS_HALF1_SLOTS if is_half1 else IMGT_PAN_CLASS_HALF2_SLOTS
    flanks_map = PAN_CLASS_EXTRA_SLOT_FLANKS[domain]
    pan_map: dict[str, str] = {}
    for slot in slots:
        if slot in fixed_slot_to_aa:
            pan_map[slot] = fixed_slot_to_aa[slot]
        else:
            flanks = flanks_map.get(slot, ())
            if flanks and all(f in imgt_to_aa for f in flanks):
                pan_map[slot] = "-"
            else:
                pan_map[slot] = "X"

    hotspots = PAN_CLASS_HALF1_INSERTION_HOTSPOTS[domain] if is_half1 else PAN_CLASS_HALF2_INSERTION_HOTSPOTS[domain]
    hotspot_seqs: list[list[str]] = [[] for _ in hotspots]
    hotspot_sources: list[list[GrooveInsertion]] = [[] for _ in hotspots]
    overflow: list[GrooveInsertion] = []
    extra_counts: dict[str, int] = {}

    for ins in insertions:
        matched_idx: Optional[int] = None
        for idx, hotspot in enumerate(hotspots):
            if ins.after_imgt_pos == hotspot.anchor_imgt_pos:
                matched_idx = idx
                break
        if matched_idx is None:
            overflow.append(ins)
            if ins.after_imgt_pos in pan_map:
                extra_counts[ins.after_imgt_pos] = extra_counts.get(ins.after_imgt_pos, 0) + ins.length
        else:
            hotspot_seqs[matched_idx].append(ins.sequence)
            hotspot_sources[matched_idx].append(ins)

    for idx, hotspot in enumerate(hotspots):
        combined = "".join(hotspot_seqs[idx])
        flanks_observed = hotspot.anchor_imgt_pos in imgt_to_aa and hotspot.right_flank_imgt_pos in imgt_to_aa
        if not combined:
            fill_ch = "-" if flanks_observed else "X"
            for s in hotspot.slots:
                pan_map[s] = fill_ch
        else:
            filled, did_overflow = _fill_hotspot_slots(combined, hotspot.slots, hotspot.fill_mode)
            pan_map.update(filled)
            if did_overflow:
                overflow.extend(hotspot_sources[idx])
                extra_len = len(combined) - len(hotspot.slots)
                anchor_slot = hotspot.slots[-1]
                extra_counts[anchor_slot] = extra_counts.get(anchor_slot, 0) + extra_len

    pan_seq = "".join(pan_map[s] for s in slots)
    insertion_counts = tuple(extra_counts.get(s, 0) for s in slots)
    return pan_seq, slots, pan_map, tuple(overflow), insertion_counts


def _domain_reference_spec(norm_domain: str) -> tuple[str, tuple[str, ...], tuple[int, ...]]:
    """Return ``(ref_seq, imgt_positions, mature_positions)`` for a normalized groove domain."""
    if norm_domain == "G-ALPHA1":
        return HLA_A0201_G_ALPHA1, G_ALPHA1_POSITIONS, tuple(range(1, G_ALPHA1_SIZE + 1))
    if norm_domain == "G-ALPHA2":
        return HLA_A0201_G_ALPHA2, G_ALPHA2_POSITIONS, tuple(G_ALPHA2_IMGT_TO_MATURE[pos] for pos in G_ALPHA2_POSITIONS)
    if norm_domain == "G-ALPHA":
        return HLA_DRA0101_G_ALPHA, G_ALPHA_POSITIONS, G_ALPHA_MATURE_POSITIONS
    if norm_domain == "G-BETA":
        return HLA_DRB1_0101_G_BETA, G_BETA_POSITIONS, G_BETA_MATURE_POSITIONS
    raise ValueError(f"Unsupported IMGT groove domain: {norm_domain!r}")


def align_domain_to_imgt(
    seq: str,
    domain: str = "G-ALPHA1",
    *,
    use_fast_path: bool = True,
) -> DomainAlignment:
    """Align a Class I or Class II groove half to canonical IMGT G-DOMAIN coordinates.

    Supports all four vertebrate MHC groove domains via separate per-domain
    multi-clade structural PSSM + loop gap-penalty profiles:
      - Class I: ``"G-ALPHA1"`` (α1, 90 cols) and ``"G-ALPHA2"`` (α2, 92 cols)
      - Class II: ``"G-ALPHA"`` (α1, 82 cols) and ``"G-BETA"`` (β1, 89 cols)

    Emits all sequence granularities from a single exact alignment:
      - Per-position coordinate maps (``mature_to_aa``, ``imgt_to_aa``) and
        structured ``insertions`` / ``deletions`` (distinct from terminal
        truncations and overhangs).
      - Lossless variable-length views: ``residues`` (:class:`AlignedResidue`
        tokens, including N- and C-terminal overhangs) and ``a3m_seq``
        (uppercase ``A-Z`` match including observed ``'X'``, lowercase ``a-z``
        internal insertion or terminal overhang, ``'-'`` for both internal
        deletions and missing match columns).
      - Class-specific fixed-index IMGT frame: ``fixed_seq`` (96 chars for
        G-ALPHA1, 102 for G-ALPHA2, 87 for G-ALPHA, 98 for G-BETA; 22-token
        alphabet with ``'-'`` vs ``'X'``) plus ``overflow_insertions``.
      - Shared 3D-homologous pan-class frame: ``pan_class_fixed_seq`` (100
        chars for Half 1, 102 chars for Half 2).

    Parameters
    ----------
    seq : str
        Amino acid sequence of the groove half.
    domain : str, optional
        ``"G-ALPHA1"``, ``"G-ALPHA2"``, ``"G-ALPHA"``, or ``"G-BETA"``
        (aliases ``"alpha1"``, ``"alpha2"``, ``"class_ii_alpha"``,
        ``"class_ii_beta"``, ``"beta1"``, etc. accepted).
    use_fast_path : bool, optional
        When True (default), reuse cached exact score/traceback results for
        repeated domain sequences. Set to False to recompute the exact same
        alignment without cache lookup.

    Returns
    -------
    DomainAlignment
        Per-position mappings, structured indels, lossless variable-length
        ``residues`` / ``a3m_seq``, class-specific ``fixed_seq``, and shared
        ``pan_class_fixed_seq``.
    """
    norm_domain = _normalize_groove_domain(domain)
    ref_seq, imgt_positions, mature_positions = _domain_reference_spec(norm_domain)

    qry = "".join(ch for ch in str(seq or "").strip().upper() if ch.isalpha())
    if not qry:
        fixed_seq, fixed_slots, slot_map, overflow = _build_fixed_frame(norm_domain, {}, ())
        pan_seq, pan_slots, pan_map, pan_overflow, pan_ins_counts = _build_pan_class_frame(norm_domain, slot_map, {}, ())
        return DomainAlignment(
            domain=norm_domain,
            query_seq="",
            ref_seq=ref_seq,
            aligned_ref=ref_seq,
            aligned_query="-" * len(ref_seq),
            mature_to_aa={},
            imgt_to_aa={},
            mature_to_query_index={},
            query_index_to_mature={},
            query_index_to_imgt={},
            n_term_missing=len(ref_seq),
            residues=(),
            a3m_seq="-" * len(ref_seq),
            fixed_seq=fixed_seq,
            fixed_slots=fixed_slots,
            fixed_slot_to_aa=slot_map,
            overflow_insertions=overflow,
            pan_class_fixed_seq=pan_seq,
            pan_class_fixed_slots=pan_slots,
            pan_class_slot_to_aa=pan_map,
            pan_class_overflow_insertions=pan_overflow,
            pan_class_insertion_counts=pan_ins_counts,
        )

    align = _cached_groove_alignment if use_fast_path else _align_groove_domain
    score, aligned_ref, aligned_qry = align(ref_seq, qry, norm_domain)

    # Identify the span of columns containing aligned (non-gap in both) pairs
    # so leading/trailing truncations and overhangs stay distinct from internal indels.
    paired_cols = [idx for idx, (r, q) in enumerate(zip(aligned_ref, aligned_qry)) if r != "-" and q != "-"]
    if not paired_cols:
        fixed_seq, fixed_slots, slot_map, overflow = _build_fixed_frame(norm_domain, {}, ())
        pan_seq, pan_slots, pan_map, pan_overflow, pan_ins_counts = _build_pan_class_frame(norm_domain, slot_map, {}, ())
        overhang_residues = tuple(
            AlignedResidue(
                aa=ch,
                domain=norm_domain,
                query_index=idx,
                mature_pos=None,
                imgt_pos=f"0.{idx + 1}",
                is_insertion=True,
                insertion_rank=idx + 1,
                structural_element="N-OVERHANG",
            )
            for idx, ch in enumerate(qry)
        )
        return DomainAlignment(
            domain=norm_domain,
            query_seq=qry,
            ref_seq=ref_seq,
            aligned_ref=aligned_ref,
            aligned_query=aligned_qry,
            mature_to_aa={},
            imgt_to_aa={},
            mature_to_query_index={},
            query_index_to_mature={},
            query_index_to_imgt={idx: f"0.{idx + 1}" for idx in range(len(qry))},
            n_term_missing=len(ref_seq),
            n_term_overhang=qry,
            score=score,
            identity=0.0,
            residues=overhang_residues,
            a3m_seq=qry.lower() + ("-" * len(ref_seq)),
            fixed_seq=fixed_seq,
            fixed_slots=fixed_slots,
            fixed_slot_to_aa=slot_map,
            overflow_insertions=overflow,
            pan_class_fixed_seq=pan_seq,
            pan_class_fixed_slots=pan_slots,
            pan_class_slot_to_aa=pan_map,
            pan_class_overflow_insertions=pan_overflow,
            pan_class_insertion_counts=pan_ins_counts,
        )

    first_paired = paired_cols[0]
    last_paired = paired_cols[-1]

    mature_to_aa: dict[int, str] = {}
    imgt_to_aa: dict[str, str] = {}
    mature_to_query_index: dict[int, int] = {}
    query_index_to_mature: dict[int, int] = {}
    query_index_to_imgt: dict[int, str] = {}
    insertions: list[GrooveInsertion] = []
    deletions: list[GrooveDeletion] = []
    residues: list[AlignedResidue] = []
    a3m_chars: list[str] = []

    n_term_missing = 0
    c_term_missing = 0
    n_term_overhang_chars: list[str] = []
    c_term_overhang_chars: list[str] = []
    n_overhang_rank = 0
    c_overhang_rank = 0
    last_imgt_pos = imgt_positions[-1]

    ref_idx = 0
    qry_idx = 0
    identical_pairs = 0
    aligned_pairs = 0

    col = 0
    total_cols = len(aligned_ref)
    while col < total_cols:
        r = aligned_ref[col]
        q = aligned_qry[col]

        if col < first_paired:
            if r != "-":
                n_term_missing += 1
                a3m_chars.append("-")
                ref_idx += 1
            if q != "-":
                n_overhang_rank += 1
                code = f"0.{n_overhang_rank}"
                query_index_to_imgt[qry_idx] = code
                residues.append(
                    AlignedResidue(
                        aa=q,
                        domain=norm_domain,
                        query_index=qry_idx,
                        mature_pos=None,
                        imgt_pos=code,
                        is_insertion=True,
                        insertion_rank=n_overhang_rank,
                        structural_element="N-OVERHANG",
                    )
                )
                a3m_chars.append(q.lower())
                n_term_overhang_chars.append(q)
                qry_idx += 1
            col += 1
            continue

        if col > last_paired:
            if r != "-":
                c_term_missing += 1
                a3m_chars.append("-")
                ref_idx += 1
            if q != "-":
                c_overhang_rank += 1
                code = f"{last_imgt_pos}.{c_overhang_rank}"
                query_index_to_imgt[qry_idx] = code
                residues.append(
                    AlignedResidue(
                        aa=q,
                        domain=norm_domain,
                        query_index=qry_idx,
                        mature_pos=None,
                        imgt_pos=code,
                        is_insertion=True,
                        insertion_rank=c_overhang_rank,
                        structural_element="C-OVERHANG",
                    )
                )
                a3m_chars.append(q.lower())
                c_term_overhang_chars.append(q)
                qry_idx += 1
            col += 1
            continue

        if r != "-" and q != "-":
            mpos = mature_positions[ref_idx]
            ipos = imgt_positions[ref_idx]
            elem = structural_element(ipos)
            mature_to_aa[mpos] = q
            imgt_to_aa[ipos] = q
            mature_to_query_index[mpos] = qry_idx
            query_index_to_mature[qry_idx] = mpos
            query_index_to_imgt[qry_idx] = ipos
            residues.append(
                AlignedResidue(
                    aa=q,
                    domain=norm_domain,
                    query_index=qry_idx,
                    mature_pos=mpos,
                    imgt_pos=ipos,
                    is_insertion=False,
                    insertion_rank=0,
                    structural_element=elem,
                )
            )
            a3m_chars.append(q)
            aligned_pairs += 1
            if r == q:
                identical_pairs += 1
            ref_idx += 1
            qry_idx += 1
            col += 1
        elif r != "-" and q == "-":
            mpos = mature_positions[ref_idx]
            ipos = imgt_positions[ref_idx]
            mature_to_aa[mpos] = "-"
            imgt_to_aa[ipos] = "-"
            a3m_chars.append("-")
            deletions.append(
                GrooveDeletion(
                    domain=norm_domain,
                    mature_pos=mpos,
                    imgt_pos=ipos,
                    structural_element=structural_element(ipos),
                )
            )
            ref_idx += 1
            col += 1
        else:
            # Internal insertion run in query (r == "-" and q != "-")
            after_ref_idx = ref_idx - 1
            after_mpos = mature_positions[after_ref_idx]
            after_ipos = imgt_positions[after_ref_idx]
            elem = structural_element(after_ipos)
            ins_start_qry = qry_idx
            ins_chars: list[str] = []
            ins_codes: list[str] = []
            k = 1
            while col <= last_paired and aligned_ref[col] == "-" and aligned_qry[col] != "-":
                ch = aligned_qry[col]
                code = f"{after_ipos}.{k}"
                ins_chars.append(ch)
                ins_codes.append(code)
                query_index_to_imgt[qry_idx] = code
                residues.append(
                    AlignedResidue(
                        aa=ch,
                        domain=norm_domain,
                        query_index=qry_idx,
                        mature_pos=after_mpos,
                        imgt_pos=code,
                        is_insertion=True,
                        insertion_rank=k,
                        structural_element=elem,
                    )
                )
                a3m_chars.append(ch.lower())
                qry_idx += 1
                k += 1
                col += 1
            insertions.append(
                GrooveInsertion(
                    domain=norm_domain,
                    after_mature_pos=after_mpos,
                    after_imgt_pos=after_ipos,
                    sequence="".join(ins_chars),
                    query_start=ins_start_qry,
                    imgt_codes=tuple(ins_codes),
                    structural_element=elem,
                )
            )

    fixed_seq, fixed_slots, slot_map, overflow = _build_fixed_frame(norm_domain, imgt_to_aa, insertions)
    pan_seq, pan_slots, pan_map, pan_overflow, pan_ins_counts = _build_pan_class_frame(norm_domain, slot_map, imgt_to_aa, insertions)
    identity = identical_pairs / aligned_pairs if aligned_pairs else 0.0
    return DomainAlignment(
        domain=norm_domain,
        query_seq=qry,
        ref_seq=ref_seq,
        aligned_ref=aligned_ref,
        aligned_query=aligned_qry,
        mature_to_aa=mature_to_aa,
        imgt_to_aa=imgt_to_aa,
        mature_to_query_index=mature_to_query_index,
        query_index_to_mature=query_index_to_mature,
        query_index_to_imgt=query_index_to_imgt,
        insertions=tuple(insertions),
        deletions=tuple(deletions),
        n_term_missing=n_term_missing,
        c_term_missing=c_term_missing,
        n_term_overhang="".join(n_term_overhang_chars),
        c_term_overhang="".join(c_term_overhang_chars),
        score=score,
        identity=identity,
        residues=tuple(residues),
        a3m_seq="".join(a3m_chars),
        fixed_seq=fixed_seq,
        fixed_slots=fixed_slots,
        fixed_slot_to_aa=slot_map,
        overflow_insertions=overflow,
        pan_class_fixed_seq=pan_seq,
        pan_class_fixed_slots=pan_slots,
        pan_class_slot_to_aa=pan_map,
        pan_class_overflow_insertions=pan_overflow,
        pan_class_insertion_counts=pan_ins_counts,
    )


def map_groove_positions(
    groove1: str = "",
    groove2: str = "",
    *,
    mhc_class: str = "I",
    use_fast_path: bool = True,
) -> GrooveAlignment:
    """Align Class I ``(groove1, groove2)`` or Class II ``(alpha_groove, beta_groove)`` onto IMGT coordinates.

    Returns a :class:`GrooveAlignment` with per-position mature and IMGT
    residue mappings, structured indels, lossless ``residues`` / ``a3m_seq``,
    class-specific ``fixed_seq`` (198 for Class I, 185 for Class II), and
    unified 202-slot ``pan_class_fixed_seq``.
    """
    cls = str(mhc_class or "I").strip().upper()
    if cls == "II":
        d1, d2 = "G-ALPHA", "G-BETA"
    elif cls == "I":
        d1, d2 = "G-ALPHA1", "G-ALPHA2"
    else:
        raise ValueError(f"Unsupported mhc_class: {mhc_class!r} (expected 'I' or 'II')")
    a1 = align_domain_to_imgt(groove1, domain=d1, use_fast_path=use_fast_path) if groove1 else None
    a2 = align_domain_to_imgt(groove2, domain=d2, use_fast_path=use_fast_path) if groove2 else None
    return GrooveAlignment(alpha1=a1, alpha2=a2, mhc_class=cls)


def map_class_ii_groove_positions(
    alpha_groove: str = "",
    beta_groove: str = "",
    *,
    use_fast_path: bool = True,
) -> GrooveAlignment:
    """Align a Class II α/β groove pair (or single chain) onto IMGT G-ALPHA / G-BETA coordinates."""
    return map_groove_positions(alpha_groove, beta_groove, mhc_class="II", use_fast_path=use_fast_path)


def align_groove_to_imgt(
    record_or_groove1: Any,
    groove2: Optional[str] = None,
    *,
    mhc_class: Optional[str] = None,
    use_fast_path: bool = True,
) -> GrooveAlignment:
    """Align a Class I or Class II allele record or ``(groove1, groove2)`` pair to IMGT G-DOMAIN coordinates.

    Accepts:
      - an :class:`~mhcseqs.domain_parsing.AlleleRecord` (or any object with
        ``.groove1`` / ``.groove2`` and optional ``.mhc_class`` attributes),
      - a mapping / CSV row dict with ``"groove1"`` and ``"groove2"`` keys
        (auto-detecting ``mhc_class="II"`` when present), or
      - explicit ``(groove1, groove2)`` strings (pass ``mhc_class="II"`` for
        Class II α/β heterodimer pairs).
    """
    if isinstance(record_or_groove1, Mapping):
        g1 = str(record_or_groove1.get("groove1") or "")
        g2 = str(record_or_groove1.get("groove2") or "")
        cls = mhc_class or str(record_or_groove1.get("mhc_class") or "I")
        return map_groove_positions(g1, g2, mhc_class=cls, use_fast_path=use_fast_path)

    if hasattr(record_or_groove1, "groove1") or hasattr(record_or_groove1, "groove2"):
        g1 = str(getattr(record_or_groove1, "groove1", "") or "")
        g2 = str(getattr(record_or_groove1, "groove2", "") or "")
        cls = mhc_class or str(getattr(record_or_groove1, "mhc_class", "") or "I")
        return map_groove_positions(g1, g2, mhc_class=cls, use_fast_path=use_fast_path)

    if isinstance(record_or_groove1, str):
        return map_groove_positions(record_or_groove1, groove2 or "", mhc_class=mhc_class or "I", use_fast_path=use_fast_path)

    raise TypeError(
        "align_groove_to_imgt expects an AlleleRecord, mapping with groove1/groove2 keys, or (groove1, groove2) strings; "
        f"got {type(record_or_groove1).__name__}"
    )


# Backwards-compatible aliases for earlier compact spellings (prefer the spaced G_ALPHA / G_BETA / PAN_CLASS names above).
GALPHA1_POSITIONS = G_ALPHA1_POSITIONS
GALPHA2_POSITIONS = G_ALPHA2_POSITIONS
GALPHA2_GAP_POSITIONS = G_ALPHA2_GAP_POSITIONS
GALPHA_POSITIONS = G_ALPHA_POSITIONS
GALPHA_GAP_POSITIONS = G_ALPHA_GAP_POSITIONS
GBETA_POSITIONS = G_BETA_POSITIONS
GBETA_GAP_POSITIONS = G_BETA_GAP_POSITIONS
HLA_A0201_GALPHA1 = HLA_A0201_G_ALPHA1
HLA_A0201_GALPHA2 = HLA_A0201_G_ALPHA2
HLA_DRA0101_GALPHA = HLA_DRA0101_G_ALPHA
HLA_DRB1_0101_GBETA = HLA_DRB1_0101_G_BETA
IMGT_GALPHA1_SLOTS = IMGT_G_ALPHA1_SLOTS
IMGT_GALPHA2_SLOTS = IMGT_G_ALPHA2_SLOTS
IMGT_GALPHA_SLOTS = IMGT_G_ALPHA_SLOTS
IMGT_GBETA_SLOTS = IMGT_G_BETA_SLOTS
IMGT_PAN_HALF1_SLOTS = IMGT_PAN_CLASS_HALF1_SLOTS
IMGT_PAN_HALF2_SLOTS = IMGT_PAN_CLASS_HALF2_SLOTS
