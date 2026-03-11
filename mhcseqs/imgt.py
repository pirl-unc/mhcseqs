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

from typing import Optional

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
    # Handle N-terminal insertions (1.1, 1.2, etc.)
    if "." in imgt_pos:
        base = int(imgt_pos.split(".")[0])
        if base <= 1:
            return "A-STRAND"
        if base == 49:
            return "D-STRAND"
        return "unknown"

    # Handle letter-suffix insertions (61A, 72A, etc.)
    num_part = imgt_pos.rstrip("AB")
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

_GALPHA1_SIZE = 90  # residues in the G-ALPHA1 domain


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
GALPHA2_GAP_POSITIONS = frozenset({"40", "41"})


def _build_galpha2_maps() -> tuple[dict[str, int], dict[int, str]]:
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


_GALPHA2_FWD, _GALPHA2_REV = _build_galpha2_maps()

#: Total residues mapped in G-ALPHA2 (for HLA-A*02:01): 92
_GALPHA2_SIZE = len(_GALPHA2_FWD)

#: Ordered IMGT positions for G-ALPHA2 (occupied only, excludes gaps).
GALPHA2_POSITIONS: tuple[str, ...] = tuple(k for k, v in sorted(_GALPHA2_FWD.items(), key=lambda x: x[1]))


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

    if mature_pos <= _GALPHA1_SIZE:
        return ("G-ALPHA1", str(mature_pos))

    if mature_pos in _GALPHA2_REV:
        return ("G-ALPHA2", _GALPHA2_REV[mature_pos])

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
        return pos if 1 <= pos <= _GALPHA1_SIZE else None

    if domain == "G-ALPHA2":
        if imgt_pos in GALPHA2_GAP_POSITIONS:
            return None
        return _GALPHA2_FWD.get(imgt_pos)

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
