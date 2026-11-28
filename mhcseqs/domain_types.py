"""Core dataclasses, grammar accessors, and sequence utilities for MHC domain parsing."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional, Sequence

if TYPE_CHECKING:
    from .domain_scoring import ParseScaffold

from .domain_grammar import (
    GRAMMAR_SPECS as _GRAMMAR_SPECS,
)
from .domain_grammar import (
    KD_SCALE as _KD_SCALE,
)
from .domain_grammar import (
    SP_CHARGED as _SP_CHARGED,
)
from .domain_grammar import (
    SP_HYDROPHOBIC as _SP_HYDROPHOBIC,
)
from .domain_grammar import (
    GrammarSpec,
)

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RawAllele:
    """Unparsed allele record: identity, species, and sequence data.

    This is what you have before groove extraction — the raw sequence
    plus metadata from the FASTA header or CSV.
    """

    allele: str = ""
    full_allele: str = ""
    gene: str = ""
    mhc_class: str = ""
    chain: str = ""  # "alpha" or "beta"
    species: str = ""
    species_category: str = ""
    species_prefix: str = ""
    source: str = ""
    source_id: str = ""  # database accession (HLA00001, NHP00001, P01901, …)
    sequence: str = ""  # full protein (may include signal peptide)
    seq_len: int = 0
    mature_start: int = 0

    @property
    def mature_sequence(self) -> str:
        """The mature protein sequence (signal peptide removed)."""
        if self.sequence:
            return self.sequence[self.mature_start :]
        return ""


@dataclass(frozen=True)
class AlleleRecord(RawAllele):
    """Fully parsed allele record with groove domains, Ig domain, and tail.

    The peptide-binding groove has two structural halves:
      groove1  (N-terminal):  α1 in class I, α1 in class II alpha
      groove2  (C-terminal):  α2 in class I, β1 in class II beta

    Class I alpha provides both halves; class II alpha and beta each provide one.
    """

    groove_seq: str = ""  # concatenation of groove halves
    groove1: str = ""  # α1 (class I & II alpha), empty (class II beta)
    groove2: str = ""  # α2 (class I), β1 (class II beta), empty (class II alpha)
    groove1_len: int = 0
    groove2_len: int = 0
    ig_domain: str = ""  # Ig support domain: α3 (class I), α2 (class II α), β2 (class II β)
    ig_domain_len: int = 0
    tail: str = ""  # everything after Ig domain (TM + cytoplasmic)
    tail_len: int = 0
    parse_score: float = 0.0
    sp_subscore: float = 0.0
    groove_subscore: float = 0.0
    ig_subscore: float = 0.0
    tail_subscore: float = 0.0
    status: str = "ok"
    anchor_type: str = ""
    anchor_cys1: Optional[int] = None
    anchor_cys2: Optional[int] = None
    secondary_cys1: Optional[int] = None
    secondary_cys2: Optional[int] = None
    candidate_type: str = ""
    nterm_state: str = ""
    support_state: str = ""
    tail_state: str = ""
    lexical_score: float = 0.0
    compatibility_score: float = 0.0
    candidate_score: float = 0.0
    mutations: tuple[str, ...] = field(default_factory=tuple)
    flags: tuple[str, ...] = field(default_factory=tuple)

    @property
    def ok(self) -> bool:
        return self.status in {
            "ok",
            "inferred_from_alpha3",
            "alpha1_only",
            "alpha2_only",
            "beta1_only_fallback",
            "fragment_fallback",
            "non_classical",
        }

    @property
    def mature_sequence(self) -> str:
        """The mature protein sequence, materialized from sequence or domain parts."""
        if self.sequence:
            return self.sequence[self.mature_start :]
        # Reconstruct from parsed domains
        if self.mhc_class == "I":
            return self.groove1 + self.groove2 + self.ig_domain + self.tail
        if self.chain == "alpha":
            return self.groove1 + self.ig_domain + self.tail
        if self.chain == "beta":
            return self.groove2 + self.ig_domain + self.tail
        return ""

    @property
    def domains(self) -> tuple["StructuralDomain", ...]:
        """Ordered structural domains implied by the parsed domain grammar."""
        from .domain_features import infer_structural_domains

        return infer_structural_domains(self)

    @property
    def domain_architecture(self) -> str:
        """Compact architecture string using source-backed domain roles."""
        return ">".join(domain.token for domain in self.domains)

    @property
    def domain_spans(self) -> str:
        """Human-readable raw-sequence spans for the inferred domain grammar."""
        return ";".join(domain.span_text for domain in self.domains)

    @property
    def parse_candidate(self) -> "ParseCandidate":
        """Structured candidate view of this parsed record."""
        from .domain_scoring import _build_parse_candidate

        return _build_parse_candidate(self)


# ---------------------------------------------------------------------------
# Structural template scoring — Cys pair classification & motif detection
# ---------------------------------------------------------------------------
#
# These functions classify Cys-Cys pairs as "groove domain" or "Ig domain"
# based on conserved flanking residues, and detect domain boundaries using
# structural motifs.  They replace per-gene/per-species position constants
# with universal biochemical signals.


@dataclass(frozen=True)
class CysPairAnnotation:
    """A Cys-Cys pair with structural classification scores."""

    c1: int
    c2: int
    separation: int
    groove_score: float  # higher = more likely groove domain (α2, β1)
    ig_score: float  # higher = more likely Ig domain (α3, β2, class II α2)
    domain_type: str  # "groove" | "ig" | "ambiguous"
    evidence: tuple[str, ...] = ()  # human-readable scoring details


@dataclass(frozen=True)
class DomainFoldAnnotation:
    """Cys pair classified by fold topology: G-domain vs C-like (Ig/C1-set).

    Uses the conserved Trp41 marker (IMGT position 41 in C-like domains) as
    the primary discriminator.  In sequence space, this Trp falls at
    approximately C1+14 (range C1+10 to C1+22) and is part of the Ig fold
    hydrophobic core.  G-domains (groove/peptide-binding) lack any Trp in
    this window.
    """

    c1: int
    c2: int
    separation: int
    fold_type: str  # "c_like" | "g_domain" | "ambiguous"
    trp_position: int  # position of the conserved Trp, or -1 if absent
    confidence: float  # 0.0-1.0
    evidence: tuple[str, ...] = ()


@dataclass
class ParseTrace:
    """Full debug trace of parsing decisions.

    Attached to an AlleleRecord when trace=True is passed to a parser.
    Contains all intermediate scoring decisions for debugging ground
    truth disagreements.
    """

    sequence_len: int = 0
    mhc_class: str = ""
    cys_pairs: tuple[CysPairAnnotation, ...] = ()
    junction_pos: int = 0
    junction_score: float = 0.0
    sp_candidates: tuple[tuple[int, float], ...] = ()  # (pos, score) top candidates
    selected_sp: int = 0
    boundary_scores: dict[str, float] = field(default_factory=dict)
    log: list[str] = field(default_factory=list)

    def summary(self) -> str:
        """Human-readable summary of parsing decisions."""
        lines = list(self.log)
        return "\n".join(lines)


@dataclass(frozen=True)
class ParseSubscores:
    """Explicit domain-level evidence for a candidate parse."""

    sp: float
    groove: float
    ig: float
    tail: float

    @property
    def total(self) -> float:
        return self.sp + self.groove + self.ig + self.tail


@dataclass(frozen=True)
class StructuralDomain:
    """Typed structural domain span in raw-sequence coordinates."""

    kind: str
    family: str
    role: str
    start: int
    end: int
    evidence: tuple[str, ...] = ()

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)

    @property
    def token(self) -> str:
        return self.role or self.kind

    @property
    def span_text(self) -> str:
        """1-based inclusive coordinates for compact human-readable output."""
        return f"{self.token}:{self.start + 1}-{self.end}"


@dataclass(frozen=True)
class BoundaryParseObjective:
    """Atomic SP / groove boundary evidence for a candidate parse.

    This is intentionally finer-grained than the public ParseSubscores
    buckets.  The model objective is SP- and groove-centric; Ig/tail terms
    are retained only as support signals when they help disambiguate the
    groove end.
    """

    sp_cleavage: float = 0.0
    sp_length: float = 0.0
    mature_nterm: float = 0.0
    groove_anchor: float = 0.0
    groove_internal: float = 0.0
    groove_junction: float = 0.0
    groove_end: float = 0.0
    groove1_length: float = 0.0
    groove2_length: float = 0.0
    groove_total_length: float = 0.0
    ig_anchor_support: float = 0.0
    ig_internal_support: float = 0.0
    ig_length_support: float = 0.0
    tail_support: float = 0.0

    @property
    def sp(self) -> float:
        return self.sp_cleavage + self.sp_length + self.mature_nterm

    @property
    def groove(self) -> float:
        return (
            self.groove_anchor
            + self.groove_internal
            + self.groove_junction
            + self.groove_end
            + self.groove1_length
            + self.groove2_length
            + self.groove_total_length
        )

    @property
    def ig(self) -> float:
        return self.ig_anchor_support + self.ig_internal_support + self.ig_length_support

    @property
    def tail(self) -> float:
        return self.tail_support

    @property
    def total(self) -> float:
        return self.sp + self.groove + self.ig + self.tail


@dataclass(frozen=True)
class ParseCandidate:
    """Whole-parse candidate with explicit architecture and evidence states."""

    candidate_type: str
    nterm_state: str
    support_state: str
    tail_state: str
    lexical_score: float
    compatibility_score: float
    total_score: float
    evidence: tuple[str, ...] = ()


def _sigmoid(x: float, center: float = 0.0, scale: float = 1.0) -> float:
    """Smooth mapping from additive score to [0, 1] factor.

    At x=center, returns 0.5.  Scores well above center approach 1.0,
    scores well below approach 0.0.  *scale* controls steepness.
    """
    import math

    try:
        return 1.0 / (1.0 + math.exp(-(x - center) / max(scale, 0.01)))
    except OverflowError:
        return 0.0 if x < center else 1.0


@dataclass
class SequenceFeatures:
    """Precomputed structural features for a single sequence.

    Computed once by ``analyze_sequence()``, then passed to parsing and
    scoring functions.  Eliminates redundant Cys-pair scanning, fold
    scoring, h-region detection, TM detection, and short-window SP
    rescoring that otherwise happen independently in each parser and in
    refinement.
    """

    seq: str
    seq_len: int

    # Cys pair annotations (computed once)
    cys_pairs: tuple[tuple[int, int, int], ...]
    pair_annotations: tuple["CysPairAnnotation", ...]
    topology_scores: dict[tuple[int, int, str], float]

    # SP features and window statistics (computed once from N-terminus)
    h_region: tuple[int, int]
    sp_estimate: int
    sp_estimate_confidence: float
    sp_estimate_candidates: tuple[tuple[int, float], ...]
    sp_shortcut_estimate: int
    sp_shortcut_confidence: float
    sp_shortcut_state: str
    sp_shortcut_kind: str
    sp_nterm_state_score: float
    sp_word_candidates: tuple[tuple[int, float], ...]
    hydrophobic_prefix: tuple[int, ...]
    charged_prefix: tuple[int, ...]
    kd_prefix: tuple[float, ...]

    # TM features (computed once from C-terminus)
    tm_span: tuple[int, int]

    # Between-cys composition (computed once per pair)
    between_cys_scores: dict[tuple[int, int], float]


def _is_truncated_class_i_sp(
    seq: str,
    mature_start: int,
    groove_boundary: int,
    *,
    cleavage_score: float,
    mature_nterm_score: float,
) -> bool:
    """Return whether ``seq[:mature_start]`` is a truncated Class I SP tail.

    Partial Class I deposits sometimes start mid-signal-peptide (without the
    initiator Met, leaving 5–16 aa of hydrophobic/C-region leader upstream of
    the canonical mature N-terminus). Recognizing these prevents the full-length
    SP length/Met penalties from leaving leader residues glued onto ``groove1``.
    """
    if seq[:1] == "M":
        return False
    if not (5 <= mature_start <= 22):
        return False
    if groove_boundary < 95:
        return False
    alpha1_len = groove_boundary - mature_start
    if not (86 <= alpha1_len <= 94):
        return False
    if cleavage_score < 5.0 or mature_nterm_score < 4.0:
        return False
    hydro_prefix = seq[: max(3, mature_start - 4)]
    hydro_frac = sum(1 for aa in hydro_prefix if aa in _SP_HYDROPHOBIC) / len(hydro_prefix)
    return hydro_frac >= 0.50


def _prefix_count(prefix: Sequence[int], start: int, end: int) -> int:
    """Return the count in [start, end) from a 1-based cumulative prefix array."""
    lo = max(0, start)
    hi = max(lo, min(len(prefix) - 1, end))
    return prefix[hi] - prefix[lo]


def _prefix_sum(prefix: Sequence[float], start: int, end: int) -> float:
    """Return the sum in [start, end) from a cumulative prefix array."""
    lo = max(0, start)
    hi = max(lo, min(len(prefix) - 1, end))
    return prefix[hi] - prefix[lo]


def _hydrophobic_fraction(
    seq: str,
    start: int,
    end: int,
    *,
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Return hydrophobic fraction for seq[start:end]."""
    lo = max(0, start)
    hi = max(lo, min(len(seq), end))
    width = hi - lo
    if width <= 0:
        return 0.0
    if features is not None:
        hydro = _prefix_count(features.hydrophobic_prefix, lo, hi)
    else:
        hydro = sum(1 for aa in seq[lo:hi] if aa in _SP_HYDROPHOBIC)
    return hydro / width


def _charged_fraction(
    seq: str,
    start: int,
    end: int,
    *,
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Return charged-residue fraction for seq[start:end]."""
    lo = max(0, start)
    hi = max(lo, min(len(seq), end))
    width = hi - lo
    if width <= 0:
        return 0.0
    if features is not None:
        charged = _prefix_count(features.charged_prefix, lo, hi)
    else:
        charged = sum(1 for aa in seq[lo:hi] if aa in _SP_CHARGED)
    return charged / width


def _kd_mean(
    seq: str,
    start: int,
    end: int,
    *,
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Return mean KD hydrophobicity over seq[start:end]."""
    lo = max(0, start)
    hi = max(lo, min(len(seq), end))
    width = hi - lo
    if width <= 0:
        return 0.0
    if features is not None:
        total = _prefix_sum(features.kd_prefix, lo, hi)
    else:
        total = sum(_KD_SCALE.get(seq[i], 0.0) for i in range(lo, hi))
    return total / width


def _gaussian_score(value: int, mean: int, sigma: int) -> float:
    """Gaussian-like score: 1.0 at mean, decaying with distance."""
    import math

    return math.exp(-0.5 * ((value - mean) / max(sigma, 1)) ** 2)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _clean_seq(sequence: Optional[str]) -> str:
    return "".join(ch for ch in str(sequence or "").strip().upper() if not ch.isspace())


# ---------------------------------------------------------------------------
# Domain-length soft priors
# ---------------------------------------------------------------------------
#
# Each entry is (typical, soft_range, hard_range).
# Soft prior: triangular bonus peaking at typical (~4 pts), tapering to 0 at
# soft range edges.  A strong motif signal (~6-8 pts) can override this.
# Hard range: steep penalty outside.  Catches truly implausible decompositions.
#
# Empirically derived from UniProt ground truth (p2–p98 of correctly-parsed).
def grammar_spec(mhc_class: str, chain: str) -> GrammarSpec:
    """Return the parameterized grammar for one class/chain pair."""
    try:
        return _GRAMMAR_SPECS[(mhc_class, chain)]
    except KeyError as exc:
        raise ValueError(f"Unsupported grammar key: {(mhc_class, chain)!r}") from exc


def _grammar_for_record(record: AlleleRecord) -> Optional[GrammarSpec]:
    """Return the grammar spec for a parsed record when one exists."""
    try:
        return grammar_spec(record.mhc_class, record.chain)
    except ValueError:
        return None


def _soft_prior_score(value: int, typical: int, soft: tuple[int, int], hard: tuple[int, int]) -> float:
    """Score a domain length against a soft prior.

    Within soft range: flat plateau at +2.0.  All plausible lengths score
    equally — a ±1 residue SP shift should NOT change which groove parse wins.
    This is critical for cross-species parsing where groove domains range
    from 75 to 110 aa across vertebrate lineages.

    Between soft and hard: linear ramp down.
    Beyond hard: steep penalty.
    """
    slo, shi = soft
    hlo, hhi = hard
    if slo <= value <= shi:
        if value <= typical:
            width = max(typical - slo, 1)
        else:
            width = max(shi - typical, 1)
        return 4.0 * (1.0 - abs(value - typical) / width)
    if hlo <= value <= hhi:
        if value < slo:
            return -(slo - value) * 1.5
        return -(value - shi) * 1.5
    dist = min(abs(value - hlo), abs(value - hhi))
    return -dist * 3.0 - 10.0


def _flags_to_tuple(flags: Sequence[str]) -> tuple[str, ...]:
    return tuple(str(flag) for flag in flags if flag)


def _slice_or_empty(seq: str, start: int, end: int) -> str:
    lo = max(0, int(start))
    hi = max(lo, min(len(seq), int(end)))
    return seq[lo:hi]


@dataclass(frozen=True)
class PrimaryParseSelection:
    """One scored anchor configuration under one chain grammar."""

    mature_start: int
    anchor: CysPairAnnotation
    secondary_pair: Optional[tuple[int, int]]
    scaffold: ParseScaffold
    selection_score: float


@dataclass(frozen=True)
class ExtractedDomains:
    """Concrete domain slices produced from a grammar-specific parse."""

    groove1: str
    groove2: str
    groove_seq: str
    ig_domain: str
    tail: str
