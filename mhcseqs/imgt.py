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

from dataclasses import dataclass
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
            base = int(base_str.rstrip("AB"))
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

#: Ordered IMGT positions for G-ALPHA1 (1 through 90).
GALPHA1_POSITIONS: tuple[str, ...] = tuple(str(i) for i in range(1, _GALPHA1_SIZE + 1))

#: Ordered IMGT positions for G-ALPHA2 (occupied only, excludes gaps).
GALPHA2_POSITIONS: tuple[str, ...] = tuple(k for k, v in sorted(_GALPHA2_FWD.items(), key=lambda x: x[1]))

#: Canonical HLA-A*02:01 reference sequences (UniProt P01892) for IMGT G-DOMAIN
#: alignment: G-ALPHA1 (mature 1-90, 90 aa) and G-ALPHA2 (mature 91-182, 92 aa).
HLA_A0201_GALPHA1 = "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEA"
HLA_A0201_GALPHA2 = "GSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRT"


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

_GALPHA2_REF_C1_IDX = 10  # mature 101 = IMGT G-ALPHA2 11
_GALPHA2_REF_C2_IDX = 73  # mature 164 = IMGT G-ALPHA2 74
_GALPHA2_CYS_ANCHOR_BONUS = 6.0


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
    """Alignment of one Class I groove domain (G-ALPHA1 or G-ALPHA2) to IMGT coordinates."""

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
    """Combined Class I G-ALPHA1 / G-ALPHA2 alignment onto IMGT G-DOMAIN coordinates."""

    alpha1: Optional[DomainAlignment] = None
    alpha2: Optional[DomainAlignment] = None

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
    def has_insertions(self) -> bool:
        return bool(self.insertions)

    @property
    def has_deletions(self) -> bool:
        return bool(self.deletions)

    @property
    def has_indels(self) -> bool:
        return bool(self.insertions or self.deletions)

    @property
    def is_truncated(self) -> bool:
        return (self.alpha1 is not None and self.alpha1.is_truncated) or (self.alpha2 is not None and self.alpha2.is_truncated)

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
                merged[("G-ALPHA1", pos)] = aa
        if self.alpha2 is not None:
            for pos, aa in self.alpha2.imgt_to_aa.items():
                merged[("G-ALPHA2", pos)] = aa
        return merged

    def residue_at(self, mature_pos: int) -> Optional[str]:
        """Return the aligned amino acid at canonical 1-indexed *mature_pos*.

        Returns ``"-"`` for an internal deletion, or ``None`` if the position
        is outside the aligned domain (e.g. N/C-terminal truncation or missing domain).
        """
        if mature_pos <= _GALPHA1_SIZE:
            return self.alpha1.mature_to_aa.get(mature_pos) if self.alpha1 is not None else None
        return self.alpha2.mature_to_aa.get(mature_pos) if self.alpha2 is not None else None

    def residue_at_imgt(self, domain: str, imgt_pos: str) -> Optional[str]:
        """Return the aligned amino acid at canonical IMGT ``(domain, imgt_pos)``."""
        norm = _normalize_class_i_domain(domain)
        if norm == "G-ALPHA1":
            return self.alpha1.imgt_to_aa.get(imgt_pos) if self.alpha1 is not None else None
        return self.alpha2.imgt_to_aa.get(imgt_pos) if self.alpha2 is not None else None

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

    def summary(self) -> str:
        """Return a concise diagnostic summary of indels and truncations."""
        parts: list[str] = []
        for label, dom in (("alpha1", self.alpha1), ("alpha2", self.alpha2)):
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


def _normalize_class_i_domain(domain: str) -> str:
    token = str(domain or "").strip().upper().replace("_", "-")
    if token in ("G-ALPHA1", "ALPHA1", "GROOVE1", "A1"):
        return "G-ALPHA1"
    if token in ("G-ALPHA2", "ALPHA2", "GROOVE2", "A2"):
        return "G-ALPHA2"
    raise ValueError(f"Unsupported Class I groove domain: {domain!r} (expected 'G-ALPHA1' or 'G-ALPHA2')")


def _galpha2_cys_anchor_bonus(qry: str) -> Optional[dict[tuple[int, int], float]]:
    """Return a gentle alignment bonus for the conserved G-ALPHA2 Cys11–Cys74 pair."""
    cys = [idx for idx, aa in enumerate(qry) if aa == "C"]
    if len(cys) < 2:
        return None
    best_pair: Optional[tuple[int, int]] = None
    best_dist = 10**9
    for i, c1 in enumerate(cys):
        for c2 in cys[i + 1 :]:
            sep = c2 - c1
            if 50 <= sep <= 76:
                dist = abs(sep - 63) * 2 + abs(c1 - _GALPHA2_REF_C1_IDX)
                if dist < best_dist:
                    best_dist = dist
                    best_pair = (c1, c2)
    if best_pair is None:
        return None
    return {
        (_GALPHA2_REF_C1_IDX, best_pair[0]): _GALPHA2_CYS_ANCHOR_BONUS,
        (_GALPHA2_REF_C2_IDX, best_pair[1]): _GALPHA2_CYS_ANCHOR_BONUS,
    }


def _try_canonical_ungapped_alignment(
    ref: str,
    qry: str,
    domain: str,
    *,
    end_gap: float = -0.5,
) -> Optional[tuple[float, str, str]]:
    """Fast-path for canonical-length, high-similarity ungapped grooves."""
    m = len(ref)
    n = len(qry)
    if domain == "G-ALPHA1":
        if n != m:
            return None
        overhang = ""
        core = qry
    else:
        if n not in (m, m + 1):
            return None
        if qry[_GALPHA2_REF_C1_IDX] != "C" or qry[_GALPHA2_REF_C2_IDX] != "C":
            return None
        core = qry[:m]
        overhang = qry[m:]

    # Verify every overlapping 15-aa window scores well ungapped so hidden
    # compensating indels (+k/-k) always fall through to full Gotoh DP.
    scores = [_BLOSUM62.get(ra, _BLOSUM62_UNKNOWN).get(qa, -1) for ra, qa in zip(ref, core)]
    matches = [1 if ra == qa else 0 for ra, qa in zip(ref, core)]
    win = 15
    for start in range(0, m - win + 1, 5):
        if sum(scores[start : start + win]) < 15 or sum(matches[start : start + win]) < 4:
            return None
    if sum(scores[m - win :]) < 15 or sum(matches[m - win :]) < 4:
        return None

    total = float(sum(scores)) + (end_gap * len(overhang))
    if overhang:
        return total, ref + ("-" * len(overhang)), core + overhang
    return total, ref, core


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


def align_domain_to_imgt(
    seq: str,
    domain: str = "G-ALPHA1",
    *,
    use_fast_path: bool = True,
) -> DomainAlignment:
    """Align a Class I groove half (G-ALPHA1 or G-ALPHA2) to canonical IMGT coordinates.

    Handles N/C-terminal truncations (such as 89-aa partial exon-2 G-ALPHA1
    sequences starting at mature residue 2), N/C-terminal overhangs (such as
    the 93rd residue at mature 183 of mhcseqs ``groove2``), and internal
    insertions or deletions across divergent vertebrate species.

    Parameters
    ----------
    seq : str
        Amino acid sequence of the groove half (``groove1`` or ``groove2``).
    domain : str, optional
        ``"G-ALPHA1"`` (default) or ``"G-ALPHA2"`` (aliases ``"alpha1"`` /
        ``"alpha2"`` / ``"groove1"`` / ``"groove2"`` accepted).
    use_fast_path : bool, optional
        When True (default), bypass full DP for exact-length, high-similarity
        ungapped canonical grooves.

    Returns
    -------
    DomainAlignment
        Per-position mature and IMGT coordinate mappings plus structured
        internal ``insertions`` and ``deletions`` (distinct from terminal
        truncations and overhangs).
    """
    norm_domain = _normalize_class_i_domain(domain)
    if norm_domain == "G-ALPHA1":
        ref_seq = HLA_A0201_GALPHA1
        imgt_positions = GALPHA1_POSITIONS
        mature_positions = tuple(range(1, _GALPHA1_SIZE + 1))
    else:
        ref_seq = HLA_A0201_GALPHA2
        imgt_positions = GALPHA2_POSITIONS
        mature_positions = tuple(_GALPHA2_FWD[pos] for pos in GALPHA2_POSITIONS)

    qry = "".join(ch for ch in str(seq or "").strip().upper() if ch.isalpha())
    if not qry:
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
        )

    fast = _try_canonical_ungapped_alignment(ref_seq, qry, norm_domain) if use_fast_path else None
    if fast is not None:
        score, aligned_ref, aligned_qry = fast
    else:
        anchor_bonus = _galpha2_cys_anchor_bonus(qry) if norm_domain == "G-ALPHA2" else None
        score, aligned_ref, aligned_qry = _semiglobal_gotoh(ref_seq, qry, anchor_bonus=anchor_bonus)

    # Identify the span of columns containing aligned (non-gap in both) pairs
    # so leading/trailing truncations and overhangs stay distinct from internal indels.
    paired_cols = [idx for idx, (r, q) in enumerate(zip(aligned_ref, aligned_qry)) if r != "-" and q != "-"]
    if not paired_cols:
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
            query_index_to_imgt={},
            n_term_missing=len(ref_seq),
            n_term_overhang=qry,
            score=score,
            identity=0.0,
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

    n_term_missing = 0
    c_term_missing = 0
    n_term_overhang_chars: list[str] = []
    c_term_overhang_chars: list[str] = []

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
                ref_idx += 1
            if q != "-":
                n_term_overhang_chars.append(q)
                qry_idx += 1
            col += 1
            continue

        if col > last_paired:
            if r != "-":
                c_term_missing += 1
                ref_idx += 1
            if q != "-":
                c_term_overhang_chars.append(q)
                qry_idx += 1
            col += 1
            continue

        if r != "-" and q != "-":
            mpos = mature_positions[ref_idx]
            ipos = imgt_positions[ref_idx]
            mature_to_aa[mpos] = q
            imgt_to_aa[ipos] = q
            mature_to_query_index[mpos] = qry_idx
            query_index_to_mature[qry_idx] = mpos
            query_index_to_imgt[qry_idx] = ipos
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
            ins_start_qry = qry_idx
            ins_chars: list[str] = []
            ins_codes: list[str] = []
            k = 1
            while col <= last_paired and aligned_ref[col] == "-" and aligned_qry[col] != "-":
                code = f"{after_ipos}.{k}"
                ins_chars.append(aligned_qry[col])
                ins_codes.append(code)
                query_index_to_imgt[qry_idx] = code
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
                    structural_element=structural_element(after_ipos),
                )
            )

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
    )


def map_groove_positions(
    groove1: str = "",
    groove2: str = "",
    *,
    use_fast_path: bool = True,
) -> GrooveAlignment:
    """Align Class I ``groove1`` (α1) and ``groove2`` (α2) onto IMGT G-DOMAIN coordinates.

    Returns a :class:`GrooveAlignment` with per-position mature and IMGT
    residue mappings plus structured internal ``insertions`` and ``deletions``.
    """
    a1 = align_domain_to_imgt(groove1, domain="G-ALPHA1", use_fast_path=use_fast_path) if groove1 else None
    a2 = align_domain_to_imgt(groove2, domain="G-ALPHA2", use_fast_path=use_fast_path) if groove2 else None
    return GrooveAlignment(alpha1=a1, alpha2=a2)


def align_groove_to_imgt(
    record_or_groove1: Any,
    groove2: Optional[str] = None,
    *,
    use_fast_path: bool = True,
) -> GrooveAlignment:
    """Align a Class I allele record or ``(groove1, groove2)`` pair to IMGT G-DOMAIN coordinates.

    Accepts:
      - an :class:`~mhcseqs.domain_parsing.AlleleRecord` (or any object with
        ``.groove1`` / ``.groove2`` attributes),
      - a mapping / CSV row dict with ``"groove1"`` and ``"groove2"`` keys, or
      - explicit ``(groove1, groove2)`` strings.
    """
    if isinstance(record_or_groove1, Mapping):
        g1 = str(record_or_groove1.get("groove1") or "")
        g2 = str(record_or_groove1.get("groove2") or "")
        return map_groove_positions(g1, g2, use_fast_path=use_fast_path)

    if hasattr(record_or_groove1, "groove1") or hasattr(record_or_groove1, "groove2"):
        g1 = str(getattr(record_or_groove1, "groove1", "") or "")
        g2 = str(getattr(record_or_groove1, "groove2", "") or "")
        return map_groove_positions(g1, g2, use_fast_path=use_fast_path)

    if isinstance(record_or_groove1, str):
        return map_groove_positions(record_or_groove1, groove2 or "", use_fast_path=use_fast_path)

    raise TypeError(
        "align_groove_to_imgt expects an AlleleRecord, mapping with groove1/groove2 keys, or (groove1, groove2) strings; "
        f"got {type(record_or_groove1).__name__}"
    )
