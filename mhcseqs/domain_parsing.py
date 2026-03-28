"""Holistic domain decomposition for MHC protein sequences.

The parser materializes the conserved vertebrate domain grammar rather than
matching one fixed mature-position constant:

  Class I alpha:   SP? -> G-alpha1 -> G-alpha2 -> C-like alpha3 -> TM? -> tail?
  Class II alpha:  SP? -> G-alpha1 -> C-like alpha2 -> TM? -> tail?
  Class II beta:   SP? -> G-beta1  -> C-like beta2 -> TM? -> tail?

Signal peptide (SP), groove, support-domain, and transmembrane evidence are
scored jointly. The live parser enumerates candidate cleavage sites, Cys pairs,
and partial architectures, then picks the highest-scoring whole parse. Active
fallback paths are still structural: class-I α3-only salvage and class-II β1-only
salvage both use relative offsets, boundary motifs, and length priors rather
than a single absolute mature-position constant.

What is scored
==============

1. SP boundary evidence
   - von Heijne -3/-1 cleavage compatibility
   - negative filters for impossible -3/-1 property combinations
   - h-region / c-region grammar and hydrophobicity drop
   - weak mature N-terminus cues used only in MHC context

2. Domain-fold evidence
   - Cys-Cys separation in the Ig/C-like range
   - Trp41-like fold topology around Cys1+14 for C-like domains
   - coarse amino-acid property patterns flanking the anchor cysteines

3. Groove-boundary evidence
   - class-I α1/α2 junction motifs
   - class-specific groove -> C-like boundary motifs
   - per-segment and total groove-length soft priors

4. Support evidence
   - C-like domain length and topology
   - transmembrane helix support to keep the groove end plausible

Fragments and SP-stripped entries
=================================

The parser supports partial sequences directly.

- SP-stripped proteins can still parse because the mature start may be zero.
- Class-I single-exon fragments can return `alpha1_only` or `alpha2_only`.
- Class-II exon-2-like fragments can return `fragment_fallback`.
- Full-length or near-full-length sequences with no credible groove now fail
  as `missing_groove` with explanatory flags instead of inventing boundaries.

References
==========

MHC domain grammar and topology:
  [M1] https://pmc.ncbi.nlm.nih.gov/articles/PMC3805034/
       Conserved MHC class I / class II domain organization across jawed vertebrates.
  [M2] https://pmc.ncbi.nlm.nih.gov/articles/PMC3913909/
       IMGT domain numbering and the G-domain versus C-like domain landmarks
       (for example the G-domain Cys11-Cys74 disulfide and C-like Cys23/Trp41/Cys104 grammar).
  [M3] https://pmc.ncbi.nlm.nih.gov/articles/PMC2386828/
       Salmonid class-II alpha/beta cysteine topology and lineage-specific extra cysteines.
  [M4] https://pmc.ncbi.nlm.nih.gov/articles/PMC4219347/
       Teleost class-II evolution, including cases where short local motifs drift but
       the modular domain grammar remains recognizable.
  [M5] https://pmc.ncbi.nlm.nih.gov/articles/PMC2434379/
       Classical class-I domain layout and conserved heavy-chain landmarks.

Signal peptide logic:
  [S1] https://pubmed.ncbi.nlm.nih.gov/6423828/
       von Heijne cleavage specificity and the canonical -3/-1 rule.
  [S2] https://pmc.ncbi.nlm.nih.gov/articles/PMC2638155/
       Flanking residues around the cleavage site, including weaker but real -2 / +1 effects.
  [S3] https://pubmed.ncbi.nlm.nih.gov/1544500/
       Strong disfavoring of Pro at +1 for signal peptidase cleavage.
  [S4] https://www.sciencedirect.com/science/article/pii/S1097276521006006
       Structural determinants of human signal peptidase recognition, including
       h-region and c-region context beyond the local -3/-1 motif.
"""

from __future__ import annotations

import json
import re
from bisect import bisect_right
from dataclasses import dataclass, field, replace
from functools import lru_cache
from typing import Optional, Sequence

from .alleles import infer_gene, normalize_mhc_class
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
    CLASS_I_ALPHA2_END_AFTER_CYS2,
    CLASS_II_ALPHA_GENE_PREFIXES,
    CLASS_II_BETA_GENE_PREFIXES,
    GENE_CLASS_I_PATTERNS,
    GENE_CLASS_II_ALPHA_PATTERNS,
    GENE_CLASS_II_BETA_PATTERNS,
    IG_DOMAIN_END_AFTER_CYS2,
    IG_SEP_MAX,
    IG_SEP_MIN,
    MAX_PLAUSIBLE_SP,
    MIN_FUNCTIONAL_GROOVE_HALF_LEN,
    MIN_GROOVE_SOURCE_LEN,
    NON_CLASSICAL_CLASS_I_GENE_PATTERNS,
    NON_GROOVE_GENES,
    NON_MHC_GENE_NAMES,
    PRIMARY_PARSE_CANDIDATE_KEEP,
    PRIMARY_PARSE_EXPANSION_MARGIN,
    PRIMARY_PARSE_LOW_CONFIDENCE,
    SP_ESTIMATE_CANDIDATE_KEEP,
    GrammarSpec,
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
    GRAMMAR_SPECS as _GRAMMAR_SPECS,
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
    GROOVE_PRIOR as _GROOVE_PRIOR,
)
from .domain_grammar import (
    GROOVE_SEGMENT_PRIOR as _GROOVE_SEGMENT_PRIOR,
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
    IG_PRIOR as _IG_PRIOR,
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
    KD_SCALE as _KD_SCALE,
)
from .domain_grammar import (
    MAMMAL_CATEGORIES as _MAMMAL_CATEGORIES,
)
from .domain_grammar import (
    MATURE_START_CLASS_I as _MATURE_START_CLASS_I,
)
from .domain_grammar import (
    MATURE_START_CLASS_II as _MATURE_START_CLASS_II,
)
from .domain_grammar import (
    MHC_BOUNDARY_MATURE_TRIADS as _MHC_BOUNDARY_MATURE_TRIADS,
)
from .domain_grammar import (
    MHC_MATURE_PREFIX_TRIADS as _MHC_MATURE_PREFIX_TRIADS,
)
from .domain_grammar import (
    MHC_SP_PRESENT_PREFIX_TRIADS as _MHC_SP_PRESENT_PREFIX_TRIADS,
)
from .domain_grammar import (
    PARSE_CANDIDATE_COMPAT_WEIGHTS as _PARSE_CANDIDATE_COMPAT_WEIGHTS,
)
from .domain_grammar import (
    PRIMARY_CANDIDATE_CONTRADICTION_WEIGHTS as _PRIMARY_CANDIDATE_CONTRADICTION_WEIGHTS,
)
from .domain_grammar import (
    REFINEMENT_DISTANCE_PENALTY as _REFINEMENT_DISTANCE_PENALTY,
)
from .domain_grammar import (
    REFINEMENT_MOTIFS as _REFINEMENT_MOTIFS,
)
from .domain_grammar import (
    REFINEMENT_WINDOW_BY_GROUP as _REFINEMENT_WINDOW_BY_GROUP,
)
from .domain_grammar import (
    RELATIVE_OFFSET_PRIOR as _RELATIVE_OFFSET_PRIOR,
)
from .domain_grammar import (
    SP_BOUNDARY_AA_CLASS as _SP_BOUNDARY_AA_CLASS,
)
from .domain_grammar import (
    SP_BOUNDARY_CLASS_WEIGHT as _SP_BOUNDARY_CLASS_WEIGHT,
)
from .domain_grammar import (
    SP_BOUNDARY_MODEL_PATH as _SP_BOUNDARY_MODEL_PATH,
)
from .domain_grammar import (
    SP_BOUNDARY_MODEL_WEIGHT_BY_GROUP as _SP_BOUNDARY_MODEL_WEIGHT_BY_GROUP,
)
from .domain_grammar import (
    SP_CHARGED as _SP_CHARGED,
)
from .domain_grammar import (
    SP_CLEAVAGE_STRONG as _SP_CLEAVAGE_STRONG,
)
from .domain_grammar import (
    SP_EARLIER_TIE_MARGIN as _SP_EARLIER_TIE_MARGIN,
)
from .domain_grammar import (
    SP_EXCLUDED_M3_M1 as _SP_EXCLUDED_M3_M1,
)
from .domain_grammar import (
    SP_H_REGION_WEIGHT as _SP_H_REGION_WEIGHT,
)
from .domain_grammar import (
    SP_HYDROPHOBIC as _SP_HYDROPHOBIC,
)
from .domain_grammar import (
    SP_LATE_DRIFT_SCALE as _SP_LATE_DRIFT_SCALE,
)
from .domain_grammar import (
    SP_PRIOR as _SP_PRIOR,
)
from .domain_grammar import (
    SP_SEQUENCE_CUE_MODEL_PATH as _SP_SEQUENCE_CUE_MODEL_PATH,
)
from .domain_grammar import (
    SP_SMALL_ALIPHATIC as _SP_SMALL_ALIPHATIC,
)
from .domain_grammar import (
    TM_CHARGED as _TM_CHARGED,
)
from .domain_grammar import (
    TM_HYDROPHOBIC as _TM_HYDROPHOBIC,
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


@dataclass(frozen=True)
class ParseScaffold:
    """Reusable parse-scoring scaffold for a fixed anchor configuration.

    The expensive structural evidence for a candidate anchor pair does not
    change while mature_start is enumerated.  This scaffold stores those
    pair-fixed terms once and only recomputes the start-dependent SP and
    length terms in the hot loop.
    """

    seq: str
    grammar: GrammarSpec
    groove_boundary: int
    groove_end: int
    fixed_groove: float
    fixed_ig: float
    fixed_tail: float
    features: Optional[SequenceFeatures] = None
    h_region: tuple[int, int] = (0, 0)
    class_i_alpha2_len: int = 0
    beta1_c1: int = 0
    beta1_start_weight: float = 0.0
    anchor_g_topology: float = 0.0
    anchor_c_topology: float = 0.0
    secondary_g_topology: float = 0.0
    secondary_c_topology: float = 0.0

    def score_components(self, mature_start: int) -> ParseSubscores:
        """Score a candidate mature start against this fixed scaffold."""
        sp = 0.0
        if mature_start > 0:
            sp = (
                _score_sp_cleavage(
                    self.seq,
                    mature_start,
                    h_region=self.h_region,
                    features=self.features,
                )
                + _score_sp_length(mature_start)
                + _score_mature_nterm(self.seq, mature_start, self.grammar.mhc_class)
                + _score_nterm_lexical_state(self.seq, mature_start)
            )

            # --- Joint SP evidence (continuous, combinatorial) ---
            # A real SP requires multiple features to be present together:
            #   - Initiator Met at position 0 (93.6% of true SPs, 0.6% of controls)
            #   - Strong h-region (hydro frac ≥0.65: 98.8% of SPs, 13.8% of controls)
            #   - Cleavage after the h-region (structural ordering)
            #
            # When features are absent, score degrades continuously.  A non-Met
            # start without a strong h-region is very unlikely to have an SP;
            # a Met start with weak h-region could be a partial or divergent SP
            # (6.4% of GT) and should require strong cleavage-site evidence.

            starts_met = self.seq[:1] == "M"
            h_start, h_end = self.h_region
            h_frac = 0.0
            if h_end > h_start:
                h_frac = _hydrophobic_fraction(
                    self.seq,
                    h_start,
                    h_end,
                    features=self.features,
                )

            cleavage_after_h = h_end > h_start and mature_start > h_end

            if not starts_met and h_frac < 0.55:
                sp -= 5.0  # no Met + no real h-region = almost certainly no SP
            elif not starts_met and h_frac < 0.65:
                sp -= 3.0  # no Met + weak h-region = very suspicious
            elif not starts_met:
                sp -= 1.0  # no Met but good h-region = suspicious, not impossible
            elif h_frac < 0.55:
                sp -= 2.0  # Met but no real h-region = suspicious
            elif h_frac < 0.65:
                sp -= 0.5  # Met but borderline h-region

            if not cleavage_after_h and h_end > h_start:
                sp -= 1.5  # cleavage before h-region = structurally wrong

        groove = self.fixed_groove
        if self.grammar.length_mode == "class_i_split":
            alpha1_len = self.groove_boundary - mature_start
            if alpha1_len > 0 and self.class_i_alpha2_len > 0:
                g1, g2, total = _score_class_i_groove_lengths(alpha1_len, self.class_i_alpha2_len)
                groove += g1 + g2 + total
            else:
                groove += -80.0
        else:
            groove_len = self.groove_boundary - mature_start
            groove += _soft_prior_score(groove_len, *_GROOVE_SEGMENT_PRIOR[self.grammar.groove_prior_key]) if groove_len > 0 else -30.0
            if self.beta1_c1 > 0:
                groove += (
                    _soft_prior_score(
                        self.beta1_c1 - mature_start,
                        *_RELATIVE_OFFSET_PRIOR[("II", "beta1_start_to_c1")],
                    )
                    * self.beta1_start_weight
                )

        return ParseSubscores(sp=sp, groove=groove, ig=self.fixed_ig, tail=self.fixed_tail)

    def score(self, mature_start: int) -> float:
        """Return the total score for *mature_start*.

        Uses factored scoring where three structural claims (SP grammar,
        domain architecture, completeness) each produce a [0,1] factor.
        The factors multiply: contradictory evidence in ANY factor drives
        the product toward zero, while missing evidence is a softer
        reduction.  This prevents a strong groove-length score from
        compensating for contradictory SP evidence or missing support.
        """
        components = self.score_components(mature_start)

        # --- Factor 1: SP grammar ---
        # Convert the additive sp score to a [0, 1] factor.
        # center=0: a sp score of 0 maps to 0.5 (neutral).
        # Positive sp = good evidence → factor > 0.5.
        # Negative sp = contradictory → factor < 0.5 (gates down).
        if mature_start > 0:
            sp_factor = _sigmoid(components.sp, center=0.0, scale=5.0)
        else:
            sp_factor = 1.0

        # --- Factor 2: Domain architecture ---
        # center=-5: only strongly negative groove+ig scores gate.
        # This is lenient — most valid parses have groove+ig > 0.
        arch_factor = _sigmoid(components.groove + components.ig, center=-5.0, scale=6.0)

        # --- Factor 3: Completeness ---
        # Full-length: missing tail is a mild reduction, not a gate.
        if self.features is not None and self.features.seq_len > 280:
            tail_factor = 0.7 + 0.3 * _sigmoid(components.tail, center=0.0, scale=0.5)
        else:
            tail_factor = 1.0

        # Combine: the additive total sets the base magnitude,
        # the multiplicative factors gate it down when evidence contradicts.
        # This preserves the additive ranking for consistent parses while
        # killing parses with contradictory structure.
        additive_total = components.total
        gate = sp_factor * arch_factor * tail_factor
        return additive_total * (0.3 + 0.7 * gate)


def _analyze_sequence_uncached(cleaned: str) -> SequenceFeatures:
    """Compute all structural features for one cleaned sequence."""
    pairs = find_cys_pairs(cleaned)
    annotations = tuple(classify_cys_pair(cleaned, c1, c2) for c1, c2, _ in pairs)

    topo: dict[tuple[int, int, str], float] = {}
    between: dict[tuple[int, int], float] = {}
    for c1, c2, _ in pairs:
        for expected in ("g_domain", "c_like"):
            topo[(c1, c2, expected)] = _score_pair_topology_support(cleaned, c1, c2, expected)
        between[(c1, c2)] = _score_between_cys(cleaned, c1, c2)

    h_start, h_end = detect_h_region(cleaned)
    sp_candidates = infer_signal_peptide_candidates(cleaned, h_region=(h_start, h_end))
    sp_est, sp_conf = infer_signal_peptide(cleaned, h_region=(h_start, h_end))
    shortcut_est, shortcut_conf, shortcut_state, shortcut_kind, shortcut_candidates = fast_sp_triage(
        cleaned,
        h_region=(h_start, h_end),
    )
    nterm_state_score = _score_sp_nterm_state(cleaned)

    tail_start = max(0, len(cleaned) - 80)
    tm_start, tm_end = _find_tm_span(cleaned[tail_start:], raw_start=tail_start)

    hydrophobic_prefix = [0]
    charged_prefix = [0]
    kd_prefix = [0.0]
    for aa in cleaned:
        hydrophobic_prefix.append(hydrophobic_prefix[-1] + int(aa in _SP_HYDROPHOBIC))
        charged_prefix.append(charged_prefix[-1] + int(aa in _SP_CHARGED))
        kd_prefix.append(kd_prefix[-1] + _KD_SCALE.get(aa, 0.0))

    return SequenceFeatures(
        seq=cleaned,
        seq_len=len(cleaned),
        cys_pairs=tuple(pairs),
        pair_annotations=annotations,
        topology_scores=topo,
        h_region=(h_start, h_end),
        sp_estimate=sp_est,
        sp_estimate_confidence=sp_conf,
        sp_estimate_candidates=sp_candidates,
        sp_shortcut_estimate=shortcut_est,
        sp_shortcut_confidence=shortcut_conf,
        sp_shortcut_state=shortcut_state,
        sp_shortcut_kind=shortcut_kind,
        sp_nterm_state_score=nterm_state_score,
        sp_word_candidates=shortcut_candidates,
        hydrophobic_prefix=tuple(hydrophobic_prefix),
        charged_prefix=tuple(charged_prefix),
        kd_prefix=tuple(kd_prefix),
        tm_span=(tm_start, tm_end),
        between_cys_scores=between,
    )


@lru_cache(maxsize=4096)
def _analyze_sequence_cached(cleaned: str) -> SequenceFeatures:
    """Cache exact-sequence feature extraction for repeated proteins."""
    return _analyze_sequence_uncached(cleaned)


def analyze_sequence(seq: str) -> SequenceFeatures:
    """Compute all structural features for *seq* in a single pass.

    The returned ``SequenceFeatures`` should be passed to
    ``decompose_class_i``, ``decompose_class_ii_alpha``,
    ``decompose_class_ii_beta``, and ``refine_signal_peptide`` via
    their ``features=`` keyword argument to avoid redundant work.
    Exact cleaned-sequence repeats reuse a cached feature bundle.
    """
    cleaned = _clean_seq(seq)
    return _analyze_sequence_cached(cleaned)


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


def _selected_sp_estimate(
    features: SequenceFeatures,
    *,
    use_early_shortcuts: bool,
) -> int:
    """Select which SP estimate to feed into structural parse enumeration.

    Exact-match shortcuts are allowed to override the grammar-based estimate
    when enabled because they are effectively dictionary lookups on highly
    repeated leaders or mature starts.  Softer lexical shortcuts still remain
    advisory only.
    """
    if use_early_shortcuts:
        if features.sp_shortcut_kind == "exact_sp_prefix30" and features.sp_shortcut_state == "sp_present" and features.sp_shortcut_estimate > 0:
            return features.sp_shortcut_estimate
        if features.sp_shortcut_kind == "exact_mature10" and features.sp_shortcut_state == "leaderless":
            return 0
    return features.sp_estimate


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


def _gaussian_score(value: int, mean: int, sigma: int) -> float:
    """Gaussian-like score: 1.0 at mean, decaying with distance."""
    import math

    return math.exp(-0.5 * ((value - mean) / max(sigma, 1)) ** 2)


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


# ---------------------------------------------------------------------------
# Universal boundary scoring — species-independent structural signals
# ---------------------------------------------------------------------------
#
# These functions use coarsened amino acid properties (_AA_PROPERTY) to
# detect domain boundaries from biophysical constraints alone, without
# clade-specific motif tables.


def sp_boundary_excluded(seq: str, pos: int) -> bool:
    """Check if a candidate SP cleavage site has an impossible -3/-1 property combo.

    Returns True if the property classes at positions -3 and -1 form a
    combination that is absent from all 2406 ground-truth MHC SP cleavages.
    This is a species-independent biophysical exclusion filter layered on top
    of the classic von Heijne cleavage rule; see [S1] and [S2] in the module
    docstring.

    Args:
        seq: full protein sequence (uppercase)
        pos: candidate cleavage position (index of first mature residue)
    """
    if pos < 3 or pos >= len(seq):
        return False
    p3 = _AA_PROPERTY.get(seq[pos - 3], "")
    p1 = _AA_PROPERTY.get(seq[pos - 1], "")
    if not p3 or not p1:
        return False
    return (p3, p1) in _SP_EXCLUDED_M3_M1


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


def detect_h_region(seq: str, *, limit: int = 50) -> tuple[int, int]:
    """Find the signal peptide hydrophobic core (h-region) in the N-terminal region.

    Uses a windowed hydrophobicity scan (similar to TM detection) to find the
    best hydrophobic window in the first *limit* residues.  The h-region
    tolerates occasional polar interruptions (e.g., Ser/Thr in the middle of
    an otherwise hydrophobic stretch), which is common in real MHC signal
    peptides.

    Returns (h_start, h_end) as 0-indexed positions, or (0, 0) if no
    plausible h-region is found. This follows the h-region -> c-region view of
    SPase substrates rather than using only a local -3/-1 motif; see [S1] and
    [S4] in the module docstring.
    """
    region = seq[: min(limit, len(seq))].upper()
    if len(region) < 8:
        return 0, 0

    best_score = -999.0
    best_span = (0, 0)

    for win_len in range(8, 16):
        if len(region) < win_len:
            continue
        for start in range(0, len(region) - win_len + 1):
            window = region[start : start + win_len]
            hydro = sum(1 for aa in window if aa in _SP_HYDROPHOBIC)
            charged = sum(1 for aa in window if aa in _SP_CHARGED)
            frac = hydro / win_len
            # Score: hydrophobic fraction drives the signal, charged residues
            # penalize.  Prefer windows that start after position ~2 (skip
            # initial Met + n-region).
            score = frac * 4.0 - charged * 1.5
            if start < 2:
                score -= 0.5
            if score > best_score:
                best_score = score
                best_span = (start, start + win_len)

    start, end = best_span
    if end <= start:
        return 0, 0

    window = region[start:end]
    hydro = sum(1 for aa in window if aa in _SP_HYDROPHOBIC)
    if hydro / len(window) < 0.50:
        return 0, 0

    # Trim trailing weakly-hydrophobic residues from the h-region end using
    # Kyte-Doolittle hydrophobicity.  The h-region core should end where the
    # hydrophobicity drops — GRAVY-like scoring is more nuanced than a binary
    # hydrophobic set.  Only trim if it leaves at least 7 residues.
    trimmed_end = end
    while trimmed_end > start + 7:
        kd = _KD_SCALE.get(region[trimmed_end - 1], 0.0)
        if kd >= 1.8:  # strongly hydrophobic (Ala=1.8 and above)
            break
        trimmed_end -= 1
    if trimmed_end > start + 6:
        end = trimmed_end

    return start, end


def estimate_sp_from_h_region(seq: str) -> int:
    """Estimate SP cleavage position from the hydrophobic core alone.

    Uses the biophysical rule that SP cleavage occurs ~5 residues after
    the h-region ends (in the polar c-region), at a small/neutral residue.
    Returns 0 if no plausible h-region is found.

    This provides a species-independent first-pass SP estimate that can
    anchor more detailed scoring. The intent is the same as the structural
    SP grammar described in [S1]-[S4] in the module docstring.
    """
    h_start, h_end = detect_h_region(seq)
    if h_end == 0:
        return 0

    seq_upper = seq.upper()
    # Search for a valid cleavage site in [h_end+1, h_end+12]
    # Prefer positions where -1 is A/G/S/C (von Heijne) and not excluded
    best_pos = 0
    best_score = -999.0
    for offset in range(1, 13):
        pos = h_end + offset
        if pos >= len(seq_upper) or pos < 3:
            continue
        if sp_boundary_excluded(seq_upper, pos):
            continue
        m1 = seq_upper[pos - 1]
        score = 0.0
        if m1 in "AG":
            score += 3.0
        elif m1 in "SC":
            score += 1.5
        elif m1 in "T":
            score += 0.5
        else:
            score -= 1.0
        # Prefer ~5 residues after h-end
        score -= 0.3 * abs(offset - 5)
        if score > best_score:
            best_score = score
            best_pos = pos

    return best_pos


def _score_h_region_alignment(
    seq: str,
    pos: int,
    *,
    h_region: tuple[int, int] | None = None,
) -> float:
    """Score a candidate SP cleavage site against the structural grammar:

        [n-region] → [h-region: 8-14 aa] → [c-region: 3-10 aa] → [cleavage]

    Signals scored (all species-independent):
      1. c-region length: soft prior peaked at 5-7, range [3, 10]
      2. c-region polarity: should be less hydrophobic than h-region
      3. Helix-breaking residue at c-region start (Pro/Gly = 26% in MHC SPs)
      4. Cleavage inside h-region: strong penalty
      5. Pro at -2: rare (0.6%) penalty

    Derived from 2373 ground truth MHC SPs: h-region median 8 aa,
    c-region median 6 aa (peak at 5, p5-p95 = [3, 10]). The scoring logic is
    the MHC-specific implementation of the general SP grammar summarized in
    [S1]-[S4] in the module docstring.
    """
    if pos <= 0 or pos >= len(seq):
        return 0.0

    if h_region is None:
        h_start, h_end = detect_h_region(seq)
    else:
        h_start, h_end = h_region
    if h_end <= h_start:
        return 0.0

    score = 0.0
    c_len = pos - h_end  # c-region length: h-region end to cleavage site

    # --- 1. c-region length prior ---
    # Peak at 5 (30.9%), 81% in [4, 9].  When multiple candidates remain
    # plausible, prefer the earlier coherent cut rather than drifting
    # farther downstream into a longer c-region.
    if c_len <= 0:
        # Cleavage inside or before h-region — biophysically implausible
        score -= 3.0
    elif 4 <= c_len <= 8:
        # Sweet spot: bonus peaked at 5-6, with a mild earlier preference.
        score += 1.6 - 0.2 * abs(c_len - 5)
    elif c_len == 9:
        score += 0.4
    elif 2 <= c_len <= 3:
        score += 0.3  # short but plausible (3.5% of cases)
    elif 10 <= c_len <= 12:
        score -= 0.4 + 0.25 * (c_len - 9)  # long c-regions drift late
    else:
        score -= 1.5  # very long or negative — suspicious

    # --- 2. c-region polarity ---
    # Mean hydrophobic fraction 0.38 in ground truth c-regions.
    # The c-region should be LESS hydrophobic than the h-region.
    if c_len > 0:
        c_region = seq[h_end:pos]
        c_hydro = sum(1 for aa in c_region if aa in _SP_HYDROPHOBIC) / len(c_region)
        if c_hydro <= 0.30:
            score += 0.5  # nicely polar
        elif c_hydro >= 0.65:
            score -= 1.0  # too hydrophobic for a c-region

    # --- 3. Helix-breaking at c-region start ---
    # Pro+Gly at the first c-region residue = 25.9% of MHC SPs.
    # Hydroxyl (S/T) = 29.9%.  Together these cover 55%.
    if h_end < len(seq):
        c_first = seq[h_end]
        if c_first in "PG":
            score += 0.6  # classic helix breaker
        elif c_first in "ST":
            score += 0.3  # polar, common at c-region start

    # --- 4. Pro at -2 ---
    # Only 0.6% of MHC SPs have Pro at -2 (barely above zero).
    if pos >= 2 and seq[pos - 2] == "P":
        score -= 0.8

    return score


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


_SP_SEQUENCE_CUE_MODEL_CACHE: Optional[dict[str, object]] = None


def _load_sp_sequence_cue_model() -> dict[str, object]:
    """Load the learned MHC-specific sequence-cue model if present."""
    global _SP_SEQUENCE_CUE_MODEL_CACHE
    if _SP_SEQUENCE_CUE_MODEL_CACHE is not None:
        return _SP_SEQUENCE_CUE_MODEL_CACHE
    try:
        with _SP_SEQUENCE_CUE_MODEL_PATH.open("r", encoding="utf-8") as f:
            _SP_SEQUENCE_CUE_MODEL_CACHE = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        _SP_SEQUENCE_CUE_MODEL_CACHE = {}
    return _SP_SEQUENCE_CUE_MODEL_CACHE


def _cue_config_value(model: dict[str, object], key: str, default: float) -> float:
    """Return one numeric cue-model config value."""
    config = model.get("config", {})
    try:
        return float(config.get(key, default))
    except (TypeError, ValueError):
        return default


def _score_sp_nterm_state(
    seq: str,
    model: Optional[dict[str, object]] = None,
) -> float:
    """Score the raw N-terminus as SP-like vs mature-like from learned stats."""
    if not seq:
        return 0.0
    model = model or _load_sp_sequence_cue_model()
    if not model:
        return 0.0
    nterm = model.get("nterm", {})
    score = 0.0
    start_m = nterm.get("start_m_log_odds", {})
    key = "present" if seq[:1] == "M" else "absent"
    try:
        score += float(start_m.get(key, 0.0))
    except (TypeError, ValueError):
        pass
    prefix3_scores = nterm.get("prefix3_log_odds", {})
    prefix3 = seq[:3] if len(seq) >= 3 else seq
    try:
        score += float(prefix3_scores.get(prefix3, 0.0))
    except (TypeError, ValueError):
        pass
    return score


def _exact_sp_shortcut(
    seq: str,
    model: dict[str, object],
) -> tuple[int, float, str, str]:
    """Return an exact-match early shortcut from learned repeated prefixes."""
    exact = model.get("exact_shortcuts", {})
    if not exact or not seq:
        return 0, 0.0, "", ""

    try:
        mature10_len = int(exact.get("mature10_len", 10))
    except (TypeError, ValueError):
        mature10_len = 10
    if len(seq) >= mature10_len:
        prefix10 = seq[:mature10_len]
        mature10 = exact.get("mature10_prefixes", {})
        if prefix10 in mature10:
            return 0, 0.995, "leaderless", "exact_mature10"

    if seq[:1] != "M":
        return 0, 0.0, "", ""

    try:
        prefix30_len = int(exact.get("sp_prefix30_len", 30))
    except (TypeError, ValueError):
        prefix30_len = 30
    if len(seq) < prefix30_len:
        return 0, 0.0, "", ""

    prefix30 = seq[:prefix30_len]
    payload = exact.get("sp_prefix30", {}).get(prefix30)
    if not isinstance(payload, dict):
        return 0, 0.0, "", ""
    try:
        split = int(payload.get("split", 0))
    except (TypeError, ValueError):
        split = 0
    if 3 <= split < len(seq):
        return split, 0.995, "sp_present", "exact_sp_prefix30"
    return 0, 0.0, "", ""


def _score_sp_boundary_words(
    seq: str,
    mature_start: int,
    model: Optional[dict[str, object]] = None,
) -> float:
    """Score a candidate cleavage boundary from exact and backed-off words."""
    if mature_start < 3 or mature_start + 3 > len(seq):
        return 0.0
    model = model or _load_sp_sequence_cue_model()
    if not model:
        return 0.0
    boundary = model.get("boundary", {})
    weights = model.get("weights", {})
    left3 = seq[mature_start - 3 : mature_start]
    right3 = seq[mature_start : mature_start + 3]
    exact6 = left3 + right3
    score = 0.0
    try:
        score += float(boundary.get("boundary6_log_odds", {}).get(exact6, 0.0)) * float(weights.get("boundary6", 1.0))
    except (TypeError, ValueError):
        pass
    try:
        score += float(boundary.get("sp_end3_log_odds", {}).get(left3, 0.0)) * float(weights.get("sp_end3", 0.4))
    except (TypeError, ValueError):
        pass
    try:
        score += float(boundary.get("mature3_log_odds", {}).get(right3, 0.0)) * float(weights.get("mature3", 0.4))
    except (TypeError, ValueError):
        pass
    return score


def _rank_sp_word_candidates(
    seq: str,
    *,
    model: Optional[dict[str, object]] = None,
    max_candidates: int = 5,
) -> tuple[tuple[int, float], ...]:
    """Return top lexical boundary candidates in the plausible SP window."""
    model = model or _load_sp_sequence_cue_model()
    if not model or len(seq) < 12:
        return ()
    config = model.get("config", {})
    try:
        search_min = int(config.get("shortcut_search_min", 8))
        search_max = int(config.get("shortcut_search_max", 45))
    except (TypeError, ValueError):
        search_min, search_max = 8, 45
    candidates: list[tuple[int, float]] = []
    for pos in range(max(3, search_min), min(search_max, len(seq) - 3) + 1):
        score = _score_sp_boundary_words(seq, pos, model)
        if score > 0.0:
            candidates.append((pos, score))
    candidates.sort(key=lambda item: item[1], reverse=True)
    return tuple(candidates[:max_candidates])


def infer_signal_peptide_candidates(
    seq: str,
    *,
    h_region: tuple[int, int] | None = None,
    max_candidates: int = SP_ESTIMATE_CANDIDATE_KEEP,
) -> tuple[tuple[int, float], ...]:
    """Return a ranked sequence-only SP cleavage zone.

    The current parser only needs a short, high-quality beam of plausible
    cleavage positions.  This keeps SP evidence flexible on divergent
    sequences without turning every parse into a wide search.
    """
    seq_upper = seq.upper() if seq else ""
    if len(seq_upper) < 15:
        return ()

    h_start, h_end = h_region if h_region is not None else detect_h_region(seq_upper)
    if h_end <= h_start or h_start > 25:
        return ()

    scored: list[tuple[int, float]] = []
    for pos in range(max(3, h_end + 3), min(len(seq_upper), h_end + 14)):
        if sp_boundary_excluded(seq_upper, pos):
            continue

        score = 0.0
        c_len = pos - h_end
        if 4 <= c_len <= 8:
            score += 2.0 - 0.3 * abs(c_len - 5)
        elif c_len == 9:
            score += 0.4
        elif c_len == 3:
            score += 0.8
        elif c_len <= 2:
            score -= 3.0
        else:
            score -= 0.8 * (c_len - 8)

        m1 = seq_upper[pos - 1]
        if m1 in _SP_CLEAVAGE_STRONG:
            score += 3.0
        elif m1 == "S":
            score += 2.0
        elif m1 == "C":
            score += 1.5
        elif m1 == "T":
            score += 1.0
        elif m1 in _SP_CHARGED or m1 == "P":
            score -= 3.0
        else:
            score -= 0.5

        if pos >= 3 and seq_upper[pos - 3] in _SP_SMALL_ALIPHATIC:
            score += 1.0
        elif pos >= 3 and seq_upper[pos - 3] in _SP_CHARGED:
            score -= 1.0

        if pos < len(seq_upper) and seq_upper[pos] == "P":
            score -= 2.0
        if pos >= 2 and seq_upper[pos - 2] == "P":
            score -= 0.8

        if h_end < len(seq_upper):
            c_first = seq_upper[h_end]
            if c_first in "PG":
                score += 0.4
            elif c_first in "ST":
                score += 0.2

        if c_len > 0:
            c_region = seq_upper[h_end:pos]
            c_hydro = sum(1 for aa in c_region if aa in _SP_HYDROPHOBIC)
            if c_hydro / len(c_region) <= 0.30:
                score += 0.3
            elif c_hydro / len(c_region) >= 0.65:
                score -= 0.8

        if pos >= 3 and pos + 3 <= len(seq_upper):
            pre_kd = sum(_KD_SCALE.get(seq_upper[i], 0.0) for i in range(pos - 3, pos)) / 3
            post_kd = sum(_KD_SCALE.get(seq_upper[i], 0.0) for i in range(pos, min(pos + 3, len(seq_upper)))) / 3
            kd_drop = pre_kd - post_kd
            if kd_drop > 0.5:
                score += min(kd_drop * 0.5, 1.5)

        if score > 0.0:
            scored.append((pos, score))

    scored.sort(key=lambda item: item[1], reverse=True)
    return tuple(scored[:max_candidates])


def fast_sp_triage(
    seq: str,
    *,
    h_region: tuple[int, int] | None = None,
) -> tuple[int, float, str, str, tuple[tuple[int, float], ...]]:
    """Return an optional early lexical SP shortcut, separate from general parsing.

    This path intentionally uses a stricter bar than the general parser:
    it only emits a shortcut when the N-terminus looks decisively SP-like
    or decisively mature-only. Otherwise callers should ignore the shortcut
    and let the full parser reason over the structural evidence.
    """
    model = _load_sp_sequence_cue_model()
    if not model or not seq:
        return 0, 0.0, "", "", ()

    exact_est, exact_conf, exact_state, exact_kind = _exact_sp_shortcut(seq, model)
    if exact_state:
        return exact_est, exact_conf, exact_state, exact_kind, ()

    nterm_score = _score_sp_nterm_state(seq, model)
    h_start, h_end = h_region if h_region is not None else detect_h_region(seq)
    h_frac = 0.0
    if h_end > h_start:
        h_frac = _hydrophobic_fraction(seq, h_start, h_end)

    leaderless_threshold = _cue_config_value(model, "leaderless_shortcut_threshold", -3.0)
    leaderless_max_hfrac = _cue_config_value(model, "leaderless_shortcut_max_hfrac", 0.45)
    if nterm_score <= leaderless_threshold and h_frac <= leaderless_max_hfrac:
        confidence = min(1.0, max(0.75, 0.55 + 0.08 * abs(nterm_score)))
        return 0, confidence, "leaderless", "lexical_leaderless", ()

    if seq[:1] != "M":
        return 0, 0.0, "", "", ()

    strong_nterm = _cue_config_value(model, "sp_shortcut_nterm_threshold", 1.0)
    min_hfrac = _cue_config_value(model, "sp_shortcut_hfrac_min", 0.65)
    min_score = _cue_config_value(model, "sp_shortcut_score_threshold", 5.0)
    min_margin = _cue_config_value(model, "sp_shortcut_margin_threshold", 1.0)
    max_c_region = int(_cue_config_value(model, "sp_shortcut_max_c_region", 12.0))
    min_cleavage_score = _cue_config_value(model, "sp_shortcut_min_cleavage_score", 2.5)

    if nterm_score < strong_nterm or h_frac < min_hfrac:
        return 0, 0.0, "", "", ()

    candidates = _rank_sp_word_candidates(seq, model=model)
    if not candidates:
        return 0, 0.0, "", "", ()

    best_pos, best_score = candidates[0]
    second_score = candidates[1][1] if len(candidates) > 1 else -999.0
    if best_score < min_score or (best_score - second_score) < min_margin:
        return 0, 0.0, "", "", candidates
    if h_end > h_start:
        c_len = best_pos - h_end
        if c_len < 1 or c_len > max_c_region:
            return 0, 0.0, "", "", candidates
    cleavage_score = _score_sp_cleavage(seq, best_pos, h_region=(h_start, h_end))
    if cleavage_score < min_cleavage_score:
        return 0, 0.0, "", "", candidates
    confidence = min(1.0, 0.55 + 0.05 * max(0.0, best_score - second_score) + 0.10 * max(0.0, h_frac - min_hfrac))
    return best_pos, confidence, "sp_present", "lexical_boundary", candidates


def infer_signal_peptide(
    seq: str,
    *,
    h_region: tuple[int, int] | None = None,
) -> tuple[int, float]:
    """Infer signal peptide cleavage position from sequence alone.

    This is the unified, class-independent SP detector.  It uses only
    biophysical signals — no Cys pair positions, no MHC class knowledge,
    no species-specific motifs.

    The inference follows the structural grammar:
        [n-region] → [h-region: hydrophobic core] → [c-region: polar] → [cleavage]

    Returns (cleavage_pos, confidence) where cleavage_pos is the index of
    the first mature residue, or (0, 0.0) if no plausible SP is found.

    Called once per sequence by all three class-specific parsers.  The result
    constrains the mature_start search range: candidate parses whose
    mature_start is far from this estimate get penalized.
    """
    candidates = infer_signal_peptide_candidates(seq, h_region=h_region, max_candidates=1)
    if not candidates:
        return 0, 0.0
    best_pos, best_score = candidates[0]

    # Confidence from score magnitude and h-region quality
    confidence = min(1.0, max(0.0, best_score / 8.0))
    return best_pos, confidence


def _score_junction(seq: str, pos: int) -> float:
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


def _find_junction(seq: str, groove_c1: int) -> tuple[int, float]:
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
        s = _score_junction(seq, pos)
        if s > best_score:
            best_score = s
            best_pos = pos
    return best_pos, best_score


def _score_groove_ig_boundary(seq: str, pos: int, mhc_class: str, chain: str = "") -> float:
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


def _score_sp_candidate(
    seq: str,
    pos: int,
    junction_pos: int,
    *,
    h_region: tuple[int, int] = (0, 0),
) -> float:
    """Score a candidate signal peptide cleavage site with junction confirmation.

    Combines von Heijne cleavage rules, hydrophobicity transition, mature
    start properties (H at +3), and structural consistency with the known
    α1/α2 junction position.

    Parameters
    ----------
    seq : str
        Full protein sequence (uppercase).
    pos : int
        Candidate mature protein start position.
    junction_pos : int
        Known or estimated α1/α2 junction position (from _find_junction).
        Used to verify that the implied α1 domain length is plausible.
    """
    if pos < 3 or pos >= len(seq) - 3:
        return -999.0

    # Fast exclusion: impossible -3/-1 property combos
    if sp_boundary_excluded(seq, pos):
        return -10.0

    score = 0.0
    n = len(seq)

    # Signal 1: von Heijne -1 cleavage residue
    m1 = seq[pos - 1]
    if m1 in _SP_CLEAVAGE_STRONG:  # A, G
        score += 3.0
    elif m1 == "S":
        score += 2.0
    elif m1 == "C":
        score += 1.5
    elif m1 == "T":
        score += 1.0
    elif m1 in _SP_CHARGED or m1 == "P":
        score -= 3.0
    else:
        score -= 0.5

    # Signal 2: von Heijne -3 residue
    if pos >= 3:
        if seq[pos - 3] in _SP_SMALL_ALIPHATIC:
            score += 1.0
        elif seq[pos - 3] in _SP_CHARGED:
            score -= 1.0

    # Signal 3: +1 not proline
    if seq[pos] == "P":
        score -= 2.0

    # Signal 4: upstream hydrophobic density (h-region)
    if pos >= 12:
        upstream = seq[pos - 12 : pos - 3]
        hydro_frac = sum(1 for c in upstream if c in _SP_HYDROPHOBIC) / max(len(upstream), 1)
        score += (hydro_frac - 0.4) * 2.0

    # Signal 5: c-region polarity
    if pos >= 5:
        c_hydro = sum(1 for c in seq[pos - 3 : pos] if c in _SP_HYDROPHOBIC) / 3
        if c_hydro > 0.66:
            score -= 1.5

    # Signal 6: H at +3 — conserved structural His in MHC α1 S1 strand
    if pos + 2 < n and seq[pos + 2] == "H":
        score += 2.0

    # Signal 7: first mature residue properties
    if pos < n:
        p1 = seq[pos]
        if p1 in _SP_CHARGED:
            score += 1.0
        elif p1 in "GAS":
            score += 0.5

    # Signal 8: KD hydrophobicity transition
    if pos >= 4 and pos + 3 <= n:
        pre_kd = sum(_KD_SCALE.get(seq[i], 0.0) for i in range(pos - 3, pos)) / 3
        post_kd = sum(_KD_SCALE.get(seq[i], 0.0) for i in range(pos, pos + 3)) / 3
        kd_drop = pre_kd - post_kd
        if kd_drop > 0.5:
            score += min(kd_drop * 0.5, 1.5)

    # Signal 9: agreement with the signal-peptide hydrophobic core
    score += _score_h_region_alignment(seq, pos, h_region=h_region) * _SP_H_REGION_WEIGHT

    # Signal 10: α1 domain length plausibility (from junction position)
    if junction_pos > 0:
        alpha1_len = junction_pos - pos
        if 80 <= alpha1_len <= 95:
            score += 2.0  # standard range
        elif 75 <= alpha1_len <= 125:
            score += 0.5  # extended range (birds, fish, DMA)
        else:
            score -= 3.0  # very unlikely

    return score


def _find_mature_start(
    seq: str,
    groove_c1: int,
    *,
    species_category: str = "",
    mhc_class: str = "I",
) -> int:
    """Find the mature protein start using structural scoring.

    Uses a scored search that combines SP cleavage signals with α1/α2
    junction confirmation and groove-length plausibility.

    For mammals: uses the traditional ±2 scan from the Cys-based prediction
    (safe, preserves existing accuracy).
    For non-mammals: wide scored search using junction discovery.

    Parameters
    ----------
    seq : str
        Full protein sequence (uppercase).
    groove_c1 : int
        Raw position of Cys1 of the groove-domain Cys pair.
    species_category : str
        Species category (used to decide mammal vs non-mammal path).
    mhc_class : str
        "I", "II" — affects junction motif search.
    """
    # Traditional estimate as fallback
    if mhc_class == "I":
        traditional = max(0, groove_c1 - 100)
    else:
        traditional = max(0, groove_c1 - 106)

    if not seq or groove_c1 <= 0:
        return traditional

    # For mammals: use the traditional estimate + ±2 refinement (safe)
    if species_category in _MAMMAL_CATEGORIES or not species_category:
        return traditional

    # For non-mammals: junction-based scored search
    # Step 1: Find the junction (α1/α2 boundary) near groove_c1 - 10
    if mhc_class == "I":
        junction_pos, junction_score = _find_junction(seq, groove_c1)
    else:
        # Class II doesn't have a groove1/groove2 junction on one chain
        # Use the Cys-based estimate for the groove-end position
        junction_pos = max(0, groove_c1 - 23)  # groove ends 23 before Ig Cys1

    # Step 2: Search for SP boundary
    # α1 domain is 75-125 aa across all vertebrates
    search_min = max(5, junction_pos - 125)
    search_max = min(junction_pos - 75, MAX_PLAUSIBLE_SP)

    if search_max <= search_min:
        # Can't find a valid range — probably SP-stripped
        return max(0, traditional)

    _h_region = detect_h_region(seq)

    best_pos = traditional
    best_score = _score_sp_candidate(seq, traditional, junction_pos, h_region=_h_region)

    for pos in range(search_min, search_max + 1):
        s = _score_sp_candidate(seq, pos, junction_pos, h_region=_h_region)
        if s > best_score:
            best_score = s
            best_pos = pos

    return best_pos


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _clean_seq(sequence: Optional[str]) -> str:
    return "".join(ch for ch in str(sequence or "").strip().upper() if not ch.isspace())


# ---------------------------------------------------------------------------
# Candidate-parse enumeration
# ---------------------------------------------------------------------------
#
# Enumerate candidate mature_start positions for each plausible Cys pair and
# score the full [SP]-[groove]-[Ig]-[tail] decomposition. Pick the highest-
# scoring whole parse rather than anchoring on one absolute position.
#
# Domain size expectations (aa):
#   Class I:   groove1 (α1) 80-100, groove2 (α2) ~93, Ig (α3) ~76
#   Class II α: groove1 (α1) 75-95, Ig (α2) ~76
#   Class II β: groove2 (β1) 85-100, Ig (β2) ~76
#   SP: 14-45 (0 = stripped)


def _score_sp_cleavage(
    seq: str,
    pos: int,
    *,
    h_region: tuple[int, int] | None = None,
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Score cleavage site quality at position *pos* (first mature residue).

    Universal biochemical signals: von Heijne -1/-3 rules, upstream
    hydrophobic density, c-region polarity, hydrophobicity transition,
    and conserved MHC mature-start properties (H at +3, charged +1).
    """
    if pos < 3 or pos >= len(seq) - 3:
        return -999.0

    # Fast exclusion: impossible -3/-1 property combos
    if sp_boundary_excluded(seq, pos):
        return -10.0

    score = 0.0
    n = len(seq)

    # -1 cleavage residue (von Heijne)
    # Among the strong residues (A/G), G is slightly preferred: SPase I
    # encounters the first valid site and cleaves there.  In 951 cases
    # where both GT and GT+1 have valid -1 residues, GT picks G at -1
    # 32% of the time but GT+1 has G at -1 50% of the time — the enzyme
    # prefers the earlier (G) cleavage.
    m1 = seq[pos - 1]
    if m1 == "G":
        score += 3.3  # slight G preference for "cleave earlier" tiebreaker
    elif m1 == "A":
        score += 3.0
    elif m1 == "S":
        score += 2.0
    elif m1 == "C":
        score += 1.5
    elif m1 == "T":
        score += 1.0
    elif m1 in _SP_CHARGED or m1 == "P":
        score -= 3.0
    else:
        score -= 0.5

    # -3 residue (von Heijne)
    if seq[pos - 3] in _SP_SMALL_ALIPHATIC:
        score += 1.0
    elif seq[pos - 3] in _SP_CHARGED:
        score -= 1.0

    # +1 not proline
    if seq[pos] == "P":
        score -= 2.0

    # Upstream hydrophobic density (h-region)
    if pos >= 12:
        hydro_frac = _hydrophobic_fraction(seq, pos - 12, pos - 3, features=features)
        score += (hydro_frac - 0.4) * 2.0

    # c-region polarity
    if pos >= 5:
        c_hydro = _hydrophobic_fraction(seq, pos - 3, pos, features=features)
        if c_hydro > 0.66:
            score -= 1.5

    # H at +3 — conserved structural His in MHC α1 S1 strand
    if pos + 2 < n and seq[pos + 2] == "H":
        score += 2.0

    # First mature residue properties
    if pos < n:
        p1 = seq[pos]
        if p1 in "AGCVIL":
            score += 0.6
        elif p1 in "ST":
            score += 0.2
        elif p1 == "D":
            score -= 0.8
        elif p1 == "E":
            score -= 0.2

    # KD hydrophobicity transition
    if pos >= 4 and pos + 3 <= n:
        pre_kd = _kd_mean(seq, pos - 3, pos, features=features)
        post_kd = _kd_mean(seq, pos, pos + 3, features=features)
        kd_drop = pre_kd - post_kd
        if kd_drop > 0.5:
            score += min(kd_drop * 0.5, 1.5)

    score += _score_h_region_alignment(seq, pos, h_region=h_region) * _SP_H_REGION_WEIGHT

    return score


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
def _grammar_spec(mhc_class: str, chain: str) -> GrammarSpec:
    """Return the parameterized grammar for one class/chain pair."""
    try:
        return _GRAMMAR_SPECS[(mhc_class, chain)]
    except KeyError as exc:
        raise ValueError(f"Unsupported grammar key: {(mhc_class, chain)!r}") from exc


def _grammar_for_record(record: AlleleRecord) -> Optional[GrammarSpec]:
    """Return the grammar spec for a parsed record when one exists."""
    try:
        return _grammar_spec(record.mhc_class, record.chain)
    except ValueError:
        return None


def _score_nterm_lexical_state(
    seq: str,
    mature_start: int,
) -> float:
    """Return an MHC-specific lexical cue score from learned sequence stats.

    This is the general lexical path, not the early shortcut.  It contributes
    moderate evidence to all parses and remains active even when early lexical
    shortcuts are disabled for benchmarking.
    """
    if not seq:
        return 0.0
    model = _load_sp_sequence_cue_model()
    if model:
        nterm_score = _score_sp_nterm_state(seq, model)
        if mature_start > 0:
            return nterm_score + _score_sp_boundary_words(seq, mature_start, model)
        return -nterm_score

    # Fallback if the learned cue model is unavailable.
    prefix3 = seq[:3] if len(seq) >= 3 else seq
    score = 0.0
    if mature_start > 0:
        if prefix3 in _MHC_SP_PRESENT_PREFIX_TRIADS:
            score += 2.0
        if prefix3 in _MHC_MATURE_PREFIX_TRIADS:
            score -= 3.0
        if mature_start + 3 <= len(seq):
            mature3 = seq[mature_start : mature_start + 3]
            if mature3 in _MHC_BOUNDARY_MATURE_TRIADS:
                score += 1.25
    else:
        if prefix3 in _MHC_MATURE_PREFIX_TRIADS:
            score += 2.5
        if prefix3 in _MHC_SP_PRESENT_PREFIX_TRIADS:
            score -= 3.0
    return score


def _full_candidate_type(grammar: GrammarSpec) -> str:
    """Return the canonical full-length candidate type for one grammar."""
    if grammar.mhc_class == "I":
        return "class_I_full"
    return f"class_II_{grammar.chain}_full"


def _missing_support_candidate_type(grammar: GrammarSpec) -> str:
    """Return the partial candidate type for a groove parse missing support."""
    if grammar.mhc_class == "I":
        return "class_I_missing_support"
    return f"class_II_{grammar.chain}_missing_support"


def _candidate_type_for_record(record: AlleleRecord) -> str:
    """Map a parsed record to a structured candidate type."""
    if record.status == "alpha1_only":
        return "class_I_alpha1_only"
    if record.status == "alpha2_only":
        return "class_I_alpha2_only"
    if record.status == "inferred_from_alpha3":
        return "class_I_alpha3_salvage"
    if record.status == "beta1_only_fallback":
        return "class_II_beta_beta1_only"
    if record.status == "fragment_fallback":
        return f"class_{record.mhc_class}_{record.chain}_fragment"
    if record.status == "missing_groove":
        return "missing_groove"
    if record.status == "non_groove":
        return "reject_non_groove_like"
    support_state = _support_state_for_record(record)
    if record.ok and support_state == "support_missing":
        if record.mhc_class == "I":
            return "class_I_missing_support"
        if record.mhc_class == "II" and record.chain == "alpha":
            return "class_II_alpha_missing_support"
        if record.mhc_class == "II" and record.chain == "beta":
            return "class_II_beta_missing_support"
    if record.mhc_class == "I" and record.ok:
        return "class_I_full"
    if record.mhc_class == "II" and record.chain == "alpha" and record.ok:
        return "class_II_alpha_full"
    if record.mhc_class == "II" and record.chain == "beta" and record.ok:
        return "class_II_beta_full"
    return record.status or "unknown"


def _nterm_state_for_record(record: AlleleRecord) -> str:
    """Classify the N-terminal state independently from architecture type."""
    if record.mature_start <= 0:
        candidate_type = record.candidate_type or _candidate_type_for_record(record)
        if record.status in {"alpha1_only", "alpha2_only", "fragment_fallback"} or "fragment" in record.status or candidate_type.endswith("_only"):
            return "sp_stripped_fragment"
        return "leaderless_mature_only"
    return "sp_present"


def _support_state_for_record(record: AlleleRecord) -> str:
    """Classify how much support-domain evidence the parse carries."""
    if record.ig_domain_len > 0:
        if record.status == "inferred_from_alpha3":
            return "support_inferred"
        return "support_present"
    if record.status in {"alpha1_only", "alpha2_only", "fragment_fallback"}:
        return "support_truncated"
    if "missing_support_anchor" in record.flags or record.status == "beta1_only_fallback":
        return "support_missing"
    if record.ok and (record.groove1_len > 0 or record.groove2_len > 0):
        return "support_missing"
    return "support_unknown"


def _tail_state_for_record(record: AlleleRecord) -> str:
    """Classify the C-terminal state using the parsed tail/TM evidence."""
    if record.tail_len <= 0:
        return "tail_absent"
    tm_start, tm_end = _find_tm_span(record.tail, raw_start=0)
    if tm_end > tm_start:
        return "tm_present"
    if record.status in {"alpha1_only", "alpha2_only", "fragment_fallback"}:
        return "cterm_truncated"
    return "tail_no_tm"


def _candidate_compatibility_score(
    *,
    candidate_type: str,
    nterm_state: str,
    support_state: str,
    tail_state: str,
    seq_len: int,
    mature_start: int,
    flags: Sequence[str] = (),
    seq: str = "",
    h_region: tuple[int, int] = (0, 0),
    sp_estimate: int = 0,
    features: Optional[SequenceFeatures] = None,
) -> tuple[float, tuple[str, ...]]:
    """Score higher-level compatibility between partial evidence sources.

    This layer distinguishes missing evidence from contradictory evidence.
    It should stay small and structural: the local residue-level evidence
    already lives in the parse score, while this function scores whether the
    chosen architecture, N-terminus state, support-domain state, and tail
    state agree with each other.
    """
    score = 0.0
    evidence: list[str] = []

    starts_met = seq[:1] == "M"
    h_start, h_end = h_region
    h_frac = 0.0
    if seq and h_end > h_start:
        h_frac = _hydrophobic_fraction(seq, h_start, h_end, features=features)

    if nterm_state == "leaderless_mature_only":
        if starts_met and h_frac >= 0.65 and sp_estimate > 0:
            score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["leaderless_conflict"]
            evidence.append("leaderless_contradicted_by_sp")
        elif (not starts_met and h_frac < 0.55) or sp_estimate <= 0:
            score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["leaderless_supported"]
            evidence.append("leaderless_supported")
    elif nterm_state == "sp_present":
        if starts_met and h_frac >= 0.65:
            score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["sp_supported"]
            evidence.append("sp_supported")
        elif (not starts_met and h_frac < 0.55) or mature_start <= 0:
            score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["sp_weak"]
            evidence.append("sp_weak")
        if sp_estimate > 0:
            delta = abs(mature_start - sp_estimate)
            if delta <= 2:
                score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["sp_estimate_close"]
                evidence.append("sp_estimate_consistent")
            elif h_frac >= 0.70 and delta >= 8:
                score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["sp_estimate_far"]
                evidence.append("sp_estimate_conflict")

    if candidate_type.endswith("_full") and support_state == "support_missing" and seq_len > 220:
        score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["full_missing_support"]
        evidence.append("full_parse_missing_support")
    if (
        support_state == "support_missing"
        and seq_len > 220
        and not candidate_type.endswith("_fragment")
        and not candidate_type.endswith("_missing_support")
        and candidate_type != "missing_groove"
        and candidate_type not in {"class_I_alpha1_only", "class_I_alpha2_only"}
    ):
        score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["support_missing_in_long_sequence"]
        evidence.append("support_missing_long_sequence")
    if candidate_type.endswith("_full") and support_state in {"support_present", "support_inferred"}:
        score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["full_support_consistent"]
        evidence.append("support_consistent")
    if candidate_type.endswith("_full") and tail_state == "tm_present":
        score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["full_tm_consistent"]
        evidence.append("tail_tm_consistent")
    if "fragment" in candidate_type:
        if seq_len > 220:
            score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["fragment_in_long_sequence"]
            evidence.append("fragment_in_long_sequence")
        else:
            score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["fragment_in_short_sequence"]
            evidence.append("fragment_length_consistent")
    if "non_groove_gene" in flags and candidate_type.endswith("_full"):
        score += _PARSE_CANDIDATE_COMPAT_WEIGHTS["non_classical_gene_context"]
        evidence.append("non_classical_gene_context")

    return score, tuple(evidence)


def _compatibility_score_for_record(record: AlleleRecord) -> tuple[float, tuple[str, ...]]:
    """Score higher-level compatibility between partial evidence sources."""
    h_region = (0, 0)
    sp_estimate = 0
    if record.sequence:
        h_region = detect_h_region(record.sequence)
        sp_estimate, _confidence, shortcut_state, _shortcut_kind, _shortcut_candidates = fast_sp_triage(
            record.sequence,
            h_region=h_region,
        )
        if shortcut_state != "sp_present" or sp_estimate == 0:
            sp_estimate, _confidence = infer_signal_peptide(record.sequence, h_region=h_region)
    return _candidate_compatibility_score(
        candidate_type=_candidate_type_for_record(record),
        nterm_state=_nterm_state_for_record(record),
        support_state=_support_state_for_record(record),
        tail_state=_tail_state_for_record(record),
        seq_len=record.seq_len,
        mature_start=record.mature_start,
        flags=record.flags,
        seq=record.sequence,
        h_region=h_region,
        sp_estimate=sp_estimate,
    )


def _build_parse_candidate(
    record: AlleleRecord,
    *,
    seq: str = "",
    h_region: tuple[int, int] | None = None,
    sp_estimate: Optional[int] = None,
    features: Optional[SequenceFeatures] = None,
) -> ParseCandidate:
    """Materialize a structured candidate view from one parsed record."""
    candidate_type = record.candidate_type or _candidate_type_for_record(record)
    nterm_state = record.nterm_state or _nterm_state_for_record(record)
    support_state = record.support_state or _support_state_for_record(record)
    tail_state = record.tail_state or _tail_state_for_record(record)
    lexical_score = record.lexical_score
    score_seq = seq or record.sequence
    if score_seq:
        lexical_score = _score_nterm_lexical_state(score_seq, record.mature_start)
    compatibility_score = record.compatibility_score
    evidence: tuple[str, ...] = ()
    if compatibility_score == 0.0:
        if score_seq:
            if h_region is None:
                h_region = features.h_region if features is not None else detect_h_region(score_seq)
            if sp_estimate is None:
                if features is not None:
                    sp_estimate = features.sp_estimate
                else:
                    shortcut_state = ""
                    shortcut_est = 0
                    if score_seq:
                        shortcut_est, _confidence, shortcut_state, _shortcut_kind, _shortcut_candidates = fast_sp_triage(
                            score_seq,
                            h_region=h_region,
                        )
                    if shortcut_state != "sp_present" or shortcut_est == 0:
                        shortcut_est, _confidence = infer_signal_peptide(score_seq, h_region=h_region)
                    sp_estimate = shortcut_est
            compatibility_score, evidence = _candidate_compatibility_score(
                candidate_type=candidate_type,
                nterm_state=nterm_state,
                support_state=support_state,
                tail_state=tail_state,
                seq_len=record.seq_len,
                mature_start=record.mature_start,
                flags=record.flags,
                seq=score_seq,
                h_region=h_region,
                sp_estimate=(sp_estimate or 0),
                features=features,
            )
        else:
            compatibility_score, evidence = _compatibility_score_for_record(record)
    total_score = record.candidate_score or (record.parse_score + compatibility_score)
    return ParseCandidate(
        candidate_type=candidate_type,
        nterm_state=nterm_state,
        support_state=support_state,
        tail_state=tail_state,
        lexical_score=lexical_score,
        compatibility_score=compatibility_score,
        total_score=total_score,
        evidence=evidence,
    )


def _attach_parse_candidate(
    record: AlleleRecord,
    *,
    seq: str = "",
    h_region: tuple[int, int] | None = None,
    sp_estimate: Optional[int] = None,
    features: Optional[SequenceFeatures] = None,
) -> AlleleRecord:
    """Populate structured candidate fields on an AlleleRecord."""
    candidate = _build_parse_candidate(
        record,
        seq=seq,
        h_region=h_region,
        sp_estimate=sp_estimate,
        features=features,
    )
    return replace(
        record,
        candidate_type=candidate.candidate_type,
        nterm_state=candidate.nterm_state,
        support_state=candidate.support_state,
        tail_state=candidate.tail_state,
        lexical_score=candidate.lexical_score,
        compatibility_score=candidate.compatibility_score,
        candidate_score=candidate.total_score,
    )


def _record_candidate_key(record: AlleleRecord) -> tuple[float, ...]:
    """Return a stable ranking key for whole-parse record candidates.

    This is intentionally type-aware but not status-hardcoded: full parses,
    support-anchored salvages, fragments, and abstentions all compete on the
    same total score, with only a small structural prior to break ties.
    """
    candidate = record.parse_candidate
    domain_count = sum(int(bool(value)) for value in (record.groove1, record.groove2, record.ig_domain, record.tail))
    candidate_type = candidate.candidate_type
    type_bonus = 0.0
    if candidate_type.endswith("_full"):
        type_bonus = 3.0
    elif candidate_type.endswith("_missing_support"):
        type_bonus = 1.5
    elif candidate_type in {"class_I_alpha3_salvage", "class_II_beta_beta1_only"}:
        type_bonus = 1.0
    elif candidate_type.endswith("_only") or candidate_type.endswith("_fragment"):
        type_bonus = 2.0
    elif candidate_type == "missing_groove":
        type_bonus = 0.0
    elif candidate_type == "reject_non_groove_like":
        type_bonus = -0.5
    structural_score = candidate.total_score + type_bonus + (0.1 * domain_count)
    return (
        float(structural_score),
        float(type_bonus),
        float(candidate.total_score),
        float(domain_count),
        float(record.parse_score),
    )


def _select_best_record_candidate(records: Sequence[AlleleRecord]) -> AlleleRecord:
    """Pick the best whole-parse record candidate from competing alternatives."""
    if not records:
        raise ValueError("No record candidates to select from")
    return max(records, key=_record_candidate_key)


def _append_record_candidate(records: list[AlleleRecord], record: Optional[AlleleRecord]) -> None:
    """Append a candidate record when it represents a usable parse alternative."""
    if record is None:
        return
    if record.status == "invalid_boundaries":
        return
    records.append(record)


def _select_best_viable_record(records: Sequence[AlleleRecord]) -> Optional[AlleleRecord]:
    """Pick the best viable record when at least one parser produced a parse."""
    viable = [record for record in records if record.ok and record.mature_start >= 0]
    if not viable:
        return None
    return _select_best_record_candidate(viable)


def _soft_prior_score(value: int, typical: int, soft: tuple[int, int], hard: tuple[int, int]) -> float:
    """Score a domain length against a soft prior.

    Within soft range: 0 to +4, peaking at typical.
    Between soft and hard: linear penalty.
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


def _score_sp_length(mature_start: int) -> float:
    """Score signal-peptide length with extra pressure against long leaders."""
    score = _soft_prior_score(mature_start, *_SP_PRIOR)
    if mature_start > 31:
        score -= (mature_start - 31) * 0.15
    if mature_start > 35:
        score -= (mature_start - 35) * 0.35
    return score


def _score_sp_estimate_consistency(
    pos: int,
    sp_estimate: int,
    confidence: float,
) -> float:
    """Score how compatible a boundary is with the sequence-only SP estimate.

    This is intentionally smooth and confidence-driven: strong sequence-only SP
    evidence should keep later structural scoring from drifting too far,
    especially downstream, while low-confidence estimates should stay weak.
    """
    if sp_estimate <= 0 or confidence <= 0.0:
        return 0.0
    delta = pos - sp_estimate
    span = max(4.0, 10.0 - 5.0 * min(confidence, 1.0))
    distance = abs(delta)
    if delta > 0:
        distance *= _SP_LATE_DRIFT_SCALE
    return max(-4.0, (1.0 + 2.5 * confidence) - (distance / span) * (2.0 + confidence))


def _score_sp_estimate_zone_consistency(
    pos: int,
    candidates: Sequence[tuple[int, float]],
    confidence: float,
) -> float:
    """Score compatibility with a short ranked zone of SP cleavage proposals.

    The sequence-only SP detector is most useful as a small beam of plausible
    cuts rather than a single hard point estimate.  This function rewards
    agreement with any strong candidate and turns into a contradiction only
    when *all* plausible cuts are far away.
    """
    if not candidates or confidence <= 0.0:
        return 0.0

    best_raw = max(float(score) for _pos, score in candidates)
    if best_raw <= 0.0:
        return 0.0

    span = max(3.0, 9.0 - 4.0 * min(confidence, 1.0))
    best_score = -999.0
    for cand_pos, cand_score in candidates:
        rel = max(0.25, min(1.0, float(cand_score) / best_raw))
        delta = pos - cand_pos
        distance = abs(delta)
        if delta > 0:
            distance *= _SP_LATE_DRIFT_SCALE
        candidate_score = rel * ((1.0 + 2.5 * confidence) - (distance / span) * (2.0 + confidence))
        if candidate_score > best_score:
            best_score = candidate_score
    return max(-4.0, best_score)


def _sp_zone_is_ambiguous(features: Optional[SequenceFeatures]) -> bool:
    """Whether the sequence-only SP zone has multiple close competitors."""
    if features is None:
        return False
    candidates = features.sp_estimate_candidates
    if len(candidates) < 2:
        return False
    (best_pos, best_score), (next_pos, next_score) = candidates[:2]
    if features.sp_estimate_confidence < PRIMARY_PARSE_LOW_CONFIDENCE:
        return False
    score_ratio = (float(next_score) / float(best_score)) if best_score > 0 else 0.0
    if abs(best_pos - next_pos) <= 4 and score_ratio >= 0.75 and (best_score - next_score) <= 1.5:
        return True
    return False


# ---------------------------------------------------------------------------
# N-terminal mature protein motifs (first ~5 residues after cleavage)
# ---------------------------------------------------------------------------
# Class I α1 N-terminus:
#   Mammals: G/S at +1, P/S at +2, H at +3 (GSH, GPH, SPH patterns)
#   Birds:   E at +1, L/P at +2, H at +3 (ELH, EPH patterns)
#   Fish:    K/D/E at +1, various at +2, H at +3
#   Universal: H at +3 is 33-52% across all vertebrates.
#
# Class II α1 N-terminus:
#   Mammals: variable first residue, often K/E/D/I
#   Less conserved than class I.
#
# Class II β1 N-terminus:
#   Mammals: often R/D/E at +1, F at +2
#   Less conserved than class I.


def _refinement_group(species_category: str) -> str:
    if species_category in _MAMMAL_CATEGORIES or not species_category:
        return "mammal"
    if species_category == "bird":
        return "bird"
    if species_category == "fish":
        return "fish"
    return "other_vertebrate"


_SP_BOUNDARY_MODEL_CACHE: Optional[dict[str, object]] = None


def _load_sp_boundary_model() -> dict[str, object]:
    """Load the learned SP motif-composition model if present."""
    global _SP_BOUNDARY_MODEL_CACHE
    if _SP_BOUNDARY_MODEL_CACHE is not None:
        return _SP_BOUNDARY_MODEL_CACHE
    try:
        with _SP_BOUNDARY_MODEL_PATH.open("r", encoding="utf-8") as f:
            _SP_BOUNDARY_MODEL_CACHE = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        _SP_BOUNDARY_MODEL_CACHE = {}
    return _SP_BOUNDARY_MODEL_CACHE


def _score_sp_boundary_composition(
    seq: str,
    mature_start: int,
    species_category: str = "",
) -> float:
    """Score a cleavage site using learned motif-composition log-odds.

    The model is trained from positives vs nearby decoys and uses both
    exact residues and coarse amino-acid classes around the boundary.
    """
    model = _load_sp_boundary_model()
    if not model:
        return 0.0
    groups = model.get("groups", {})
    group = _refinement_group(species_category)
    payload = groups.get(group) or groups.get("global")
    if not payload:
        return 0.0

    residue_scores = payload.get("residue_log_odds", {})
    class_scores = payload.get("class_log_odds", {})
    score = 0.0
    for offset_text, residue_table in residue_scores.items():
        try:
            offset = int(offset_text)
        except ValueError:
            continue
        pos = mature_start + offset
        if not (0 <= pos < len(seq)):
            continue
        aa = seq[pos]
        score += float(residue_table.get(aa, 0.0))
        aa_class = _SP_BOUNDARY_AA_CLASS.get(aa)
        if aa_class:
            score += float(class_scores.get(offset_text, {}).get(aa_class, 0.0)) * _SP_BOUNDARY_CLASS_WEIGHT
    return score * _SP_BOUNDARY_MODEL_WEIGHT_BY_GROUP.get(group, 0.0)


def _score_class_i_groove_lengths(alpha1_len: int, alpha2_len: int) -> tuple[float, float, float]:
    """Return class-I groove1 / groove2 / total length evidence."""
    g1 = _soft_prior_score(alpha1_len, *_GROOVE_SEGMENT_PRIOR[("I", "groove1")])
    g2 = _soft_prior_score(alpha2_len, *_GROOVE_SEGMENT_PRIOR[("I", "groove2")])
    total = _soft_prior_score(alpha1_len + alpha2_len, *_GROOVE_SEGMENT_PRIOR[("I", "total")]) * 0.5
    return g1, g2, total


def _score_mature_nterm(seq: str, mature_start: int, mhc_class: str) -> float:
    """Score the first ~5 residues of the mature protein."""
    if mature_start <= 0:
        return 0.0
    n = len(seq)
    tables = _MATURE_START_CLASS_I if mhc_class == "I" else _MATURE_START_CLASS_II
    score = 0.0
    for offset, residue_scores in tables.items():
        pos = mature_start + offset
        if pos < n:
            score += residue_scores.get(seq[pos], 0.0)
    return score


def _score_refinement_motif(
    seq: str,
    mature_start: int,
    species_category: str,
    mhc_class: str = "",
) -> float:
    """Score mature-start motif composition around a candidate boundary."""
    if mature_start < 0 or mature_start + 3 > len(seq):
        return 0.0
    group = _refinement_group(species_category)
    motif = seq[mature_start : mature_start + 3]
    score = _score_sp_boundary_composition(seq, mature_start, species_category)
    score += _REFINEMENT_MOTIFS.get(group, {}).get(motif, 0.0)
    if mhc_class in {"I", "II"}:
        score += _score_mature_nterm(seq, mature_start, mhc_class) * 0.5
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
# Comprehensive parse scorer
# ---------------------------------------------------------------------------


def _build_parse_scaffold(
    seq: str,
    anchor_c1: int,
    anchor_c2: int,
    mhc_class: str,
    chain: str,
    secondary_c1: int = 0,
    secondary_c2: int = 0,
    anchor_groove_score: float = 0.0,
    anchor_ig_score: float = 0.0,
    *,
    features: Optional[SequenceFeatures] = None,
) -> ParseScaffold:
    """Build pair-fixed scoring terms for a candidate anchor configuration."""
    grammar = _grammar_spec(mhc_class, chain)
    weights = grammar.weights
    n = len(seq)
    anchor_g_topology = _cached_topology_support(seq, anchor_c1, anchor_c2, "g_domain", features)
    anchor_c_topology = _cached_topology_support(seq, anchor_c1, anchor_c2, "c_like", features)
    anchor_between = _cached_between_cys(seq, anchor_c1, anchor_c2, features)

    secondary_g_topology = 0.0
    secondary_c_topology = 0.0
    secondary_between = 0.0
    if secondary_c1 > 0 and secondary_c2 > 0:
        secondary_g_topology = _cached_topology_support(seq, secondary_c1, secondary_c2, "g_domain", features)
        secondary_c_topology = _cached_topology_support(seq, secondary_c1, secondary_c2, "c_like", features)
        secondary_between = _cached_between_cys(seq, secondary_c1, secondary_c2, features)

    groove_boundary = anchor_c1 - grammar.groove_boundary_from_c1
    groove_end = groove_boundary if grammar.groove_end_from_c2 == 0 else anchor_c2 + grammar.groove_end_from_c2

    ig_end = 0
    if not grammar.uses_primary_anchor_for_ig and secondary_c1 > 0 and secondary_c2 > 0:
        ig_end = min(n, secondary_c2 + IG_DOMAIN_END_AFTER_CYS2)
    elif grammar.uses_primary_anchor_for_ig:
        ig_end = min(n, anchor_c2 + IG_DOMAIN_END_AFTER_CYS2)

    fixed_groove = 0.0
    fixed_ig = 0.0
    class_i_alpha2_len = 0
    beta1_c1 = 0
    beta1_start_weight = 0.0

    fixed_groove = (
        anchor_groove_score * weights.groove_anchor_score
        + anchor_g_topology * weights.groove_anchor_g_topology
        + anchor_between * weights.groove_anchor_between
    )
    if weights.groove_junction and groove_boundary > 0 and groove_boundary + 5 < n:
        fixed_groove += _score_junction(seq, groove_boundary) * weights.groove_junction
    if weights.groove_boundary and 5 < groove_end < n - 5:
        fixed_groove += _score_groove_ig_boundary(seq, groove_end, grammar.mhc_class, grammar.chain) * weights.groove_boundary

    if grammar.length_mode == "class_i_split":
        class_i_alpha2_len = groove_end - groove_boundary
        if secondary_c1 > 0 and secondary_c2 > 0:
            ig_start = groove_end
            ig_len = ig_end - ig_start
            fixed_ig = (
                anchor_ig_score * weights.ig_anchor_score
                + secondary_c_topology * weights.ig_secondary_c_topology
                + secondary_between * weights.ig_secondary_between
                + (_soft_prior_score(ig_len, *_IG_PRIOR) * weights.ig_length if ig_len > 0 else -8.0)
            )
    else:
        if secondary_c1 > 0 and secondary_c2 > 0:
            fixed_groove += secondary_g_topology * weights.groove_secondary_g_topology + secondary_between * weights.groove_secondary_between
            fixed_groove += (
                _soft_prior_score(
                    groove_boundary - secondary_c2,
                    *_RELATIVE_OFFSET_PRIOR[("II", "beta1_c2_to_end")],
                )
                * weights.groove_secondary_end_offset
            )
            if weights.beta1_start_weight:
                beta1_c1 = secondary_c1
                beta1_start_weight = weights.beta1_start_weight
        if ig_end > 0:
            ig_start = groove_boundary
            ig_len = ig_end - ig_start
            fixed_ig = (
                anchor_ig_score * weights.ig_anchor_score
                + anchor_c_topology * weights.ig_anchor_c_topology
                + anchor_between * weights.ig_anchor_between
                + (_soft_prior_score(ig_len, *_IG_PRIOR) * weights.ig_length if ig_len > 0 else -8.0)
            )

    tail_start = ig_end if ig_end > 0 else groove_end
    tm_span = features.tm_span if features is not None else (0, 0)
    fixed_tail = _score_tail(seq, tail_start, tm_span=tm_span) * weights.tail

    return ParseScaffold(
        seq=seq,
        grammar=grammar,
        groove_boundary=groove_boundary,
        groove_end=groove_end,
        fixed_groove=fixed_groove,
        fixed_ig=fixed_ig,
        fixed_tail=fixed_tail,
        features=features,
        h_region=(features.h_region if features is not None else detect_h_region(seq)),
        class_i_alpha2_len=class_i_alpha2_len,
        beta1_c1=beta1_c1,
        beta1_start_weight=beta1_start_weight,
        anchor_g_topology=anchor_g_topology,
        anchor_c_topology=anchor_c_topology,
        secondary_g_topology=secondary_g_topology,
        secondary_c_topology=secondary_c_topology,
    )


def _score_parse_components(
    seq: str,
    mature_start: int,
    anchor_c1: int,
    anchor_c2: int,
    mhc_class: str,
    chain: str,
    secondary_c1: int = 0,
    secondary_c2: int = 0,
    anchor_groove_score: float = 0.0,
    anchor_ig_score: float = 0.0,
    *,
    features: Optional[SequenceFeatures] = None,
    scaffold: Optional[ParseScaffold] = None,
) -> ParseSubscores:
    """Return explicit SP / groove / support evidence for a candidate parse."""
    if scaffold is None:
        scaffold = _build_parse_scaffold(
            seq,
            anchor_c1,
            anchor_c2,
            mhc_class,
            chain,
            secondary_c1,
            secondary_c2,
            anchor_groove_score,
            anchor_ig_score,
            features=features,
        )
    return scaffold.score_components(mature_start)


def _score_parse(
    seq: str,
    mature_start: int,
    anchor_c1: int,
    anchor_c2: int,
    mhc_class: str,
    chain: str,
    secondary_c1: int = 0,
    secondary_c2: int = 0,
    anchor_groove_score: float = 0.0,
    anchor_ig_score: float = 0.0,
    *,
    features: Optional[SequenceFeatures] = None,
    scaffold: Optional[ParseScaffold] = None,
) -> float:
    """Score a complete candidate parse.

    Evaluates the full domain decomposition implied by (mature_start, Cys pair):

    Class I:  SP? - groove1(α1) - groove2(α2) - Ig?(α3) - tail?
                    mature_start   c1-10        c2+20     ic2+20

    Class II: SP? - groove(α1/β1) - Ig?(α2/β2) - tail?
                    mature_start    c1-23        c2+20

    Signals:
    1. SP cleavage quality (von Heijne, hydrophobicity, KD transition)
    2. Mature protein N-terminal motifs (first 3-5 aa, class-specific)
    3. Anchor Cys pair classification (groove_score for class I, ig_score for II)
    4. Flanking motifs around anchor Cys pair (already in groove/ig scores)
    5. Between-cys compositional properties
    6. Junction motif at α1/α2 boundary (class I: G..H..Q)
    7. Groove/Ig boundary motif (class-specific linker residues)
    8. Secondary Cys pair quality (Ig for class I)
    9. Tail transmembrane helix detection
    10. Domain-length soft priors (groove, Ig, SP)
    """
    return _score_parse_components(
        seq,
        mature_start,
        anchor_c1,
        anchor_c2,
        mhc_class,
        chain,
        secondary_c1,
        secondary_c2,
        anchor_groove_score,
        anchor_ig_score,
        features=features,
        scaffold=scaffold,
    ).total


def _score_class_i_alpha3_salvage(
    seq: str,
    mature_start: int,
    alpha2_start: int,
    alpha3_start: int,
    alpha3_c1: int,
    alpha3_c2: int,
    *,
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Score a class-I parse salvaged from an α3 C-like anchor only.

    This is used when the α2 groove disulfide is missing or too divergent.
    The parse is still chosen holistically from relative structural evidence:
    SP boundary, α1/α2 lengths, α1/α2 junction, α2→α3 boundary, α3 topology,
    and α3 internal offset to its disulfide.
    """
    if not (0 <= mature_start < alpha2_start < alpha3_start < len(seq)):
        return -999.0
    if not (alpha3_start < alpha3_c1 < alpha3_c2 < len(seq)):
        return -999.0

    alpha1_len = alpha2_start - mature_start
    alpha2_len = alpha3_start - alpha2_start
    if alpha1_len <= 0 or alpha2_len <= 0:
        return -999.0

    ig_end = min(len(seq), alpha3_c2 + IG_DOMAIN_END_AFTER_CYS2)
    ig_len = ig_end - alpha3_start
    alpha3_offset = alpha3_c1 - alpha3_start

    score = 0.0
    if mature_start > 0:
        score += _score_sp_cleavage(
            seq,
            mature_start,
            h_region=(features.h_region if features is not None else None),
            features=features,
        )
        score += _score_sp_length(mature_start)
        score += _score_mature_nterm(seq, mature_start, "I")
        score += _score_nterm_lexical_state(seq, mature_start)

    g1_score, g2_score, total_score = _score_class_i_groove_lengths(alpha1_len, alpha2_len)
    score += g1_score + g2_score + total_score
    score += _soft_prior_score(alpha3_offset, *_RELATIVE_OFFSET_PRIOR[("I", "alpha3_start_to_c1")]) * 0.7
    score += _score_junction(seq, alpha2_start) * 0.8 if alpha2_start + 5 < len(seq) else -8.0
    score += _score_groove_ig_boundary(seq, alpha3_start, "I", "alpha") if 5 < alpha3_start < len(seq) - 5 else -2.0
    score += _cached_topology_support(seq, alpha3_c1, alpha3_c2, "c_like", features) * 0.9
    score += _cached_between_cys(seq, alpha3_c1, alpha3_c2, features) * 0.25
    score += _soft_prior_score(ig_len, *_IG_PRIOR) * 0.35 if ig_len > 0 else -8.0
    score += _score_tail(seq, ig_end, tm_span=(features.tm_span if features is not None else (0, 0))) * 0.25
    return score


def _enumerate_class_i_from_alpha3(
    seq: str,
    alpha3_c1: int,
    alpha3_c2: int,
    *,
    features: Optional[SequenceFeatures] = None,
) -> tuple[int, int, int, float] | None:
    """Enumerate class-I parses using only the downstream α3 C-like anchor."""
    best: tuple[int, int, int, float] | None = None
    best_score = -999.0

    _typical, _soft, alpha3_offset_hard = _RELATIVE_OFFSET_PRIOR[("I", "alpha3_start_to_c1")]
    off_lo, off_hi = alpha3_offset_hard
    g2_hlo, g2_hhi = _GROOVE_SEGMENT_PRIOR[("I", "groove2")][2]
    g1_hlo, g1_hhi = _GROOVE_SEGMENT_PRIOR[("I", "groove1")][2]

    alpha3_start_lo = max(1, alpha3_c1 - off_hi)
    alpha3_start_hi = min(alpha3_c1 - 1, alpha3_c1 - off_lo)
    for alpha3_start in range(alpha3_start_lo, alpha3_start_hi + 1):
        alpha2_start_lo = max(1, alpha3_start - g2_hhi)
        alpha2_start_hi = min(alpha3_start - 1, alpha3_start - g2_hlo)
        for alpha2_start in range(alpha2_start_lo, alpha2_start_hi + 1):
            search_lo = max(0, alpha2_start - g1_hhi)
            search_hi = min(MAX_PLAUSIBLE_SP, alpha2_start - g1_hlo)
            if search_hi < search_lo:
                continue
            for mature_start in range(search_lo, search_hi + 1):
                score = _score_class_i_alpha3_salvage(
                    seq,
                    mature_start,
                    alpha2_start,
                    alpha3_start,
                    alpha3_c1,
                    alpha3_c2,
                    features=features,
                )
                if score > best_score:
                    best_score = score
                    best = (mature_start, alpha2_start, alpha3_start, score)

    if best is None:
        return None
    return best


def _score_class_ii_beta1_salvage(
    seq: str,
    mature_start: int,
    groove_end: int,
    beta1_c1: int,
    beta1_c2: int,
    beta1_groove_score: float,
    *,
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Score a class-II beta parse salvaged from the β1 groove disulfide only."""
    if not (0 <= mature_start < beta1_c1 < beta1_c2 < groove_end <= len(seq)):
        return -999.0

    groove_len = groove_end - mature_start
    start_offset = beta1_c1 - mature_start
    end_offset = groove_end - beta1_c2

    score = 0.0
    if mature_start > 0:
        score += _score_sp_cleavage(
            seq,
            mature_start,
            h_region=(features.h_region if features is not None else None),
            features=features,
        )
        score += _score_sp_length(mature_start)
        score += _score_mature_nterm(seq, mature_start, "II")
        score += _score_nterm_lexical_state(seq, mature_start)

    score += _soft_prior_score(groove_len, *_GROOVE_SEGMENT_PRIOR[("II", "beta")])
    score += _soft_prior_score(start_offset, *_RELATIVE_OFFSET_PRIOR[("II", "beta1_start_to_c1")]) * 0.7
    score += _soft_prior_score(end_offset, *_RELATIVE_OFFSET_PRIOR[("II", "beta1_c2_to_end")]) * 0.7
    score += beta1_groove_score * 0.35
    score += _cached_topology_support(seq, beta1_c1, beta1_c2, "g_domain", features) * 0.9
    score += _cached_between_cys(seq, beta1_c1, beta1_c2, features) * 0.45
    score += _score_groove_ig_boundary(seq, groove_end, "II", "beta") if 5 < groove_end < len(seq) - 5 else -1.0
    score += _score_tail(seq, groove_end, tm_span=(features.tm_span if features is not None else (0, 0))) * 0.15
    return score


def _enumerate_class_ii_beta_from_beta1(
    seq: str,
    beta1_c1: int,
    beta1_c2: int,
    beta1_groove_score: float,
    *,
    features: Optional[SequenceFeatures] = None,
) -> tuple[int, int, float] | None:
    """Enumerate class-II beta parses using only the β1 groove disulfide."""
    best: tuple[int, int, float] | None = None
    best_score = -999.0

    start_hlo, start_hhi = _RELATIVE_OFFSET_PRIOR[("II", "beta1_start_to_c1")][2]
    end_hlo, end_hhi = _RELATIVE_OFFSET_PRIOR[("II", "beta1_c2_to_end")][2]
    groove_hlo, groove_hhi = _GROOVE_SEGMENT_PRIOR[("II", "beta")][2]

    start_lo = max(0, beta1_c1 - start_hhi)
    start_hi = min(MAX_PLAUSIBLE_SP, beta1_c1 - start_hlo)
    end_lo = max(beta1_c2 + 1, beta1_c2 + end_hlo)
    end_hi = min(len(seq), beta1_c2 + end_hhi)

    if start_hi < start_lo or end_hi < end_lo:
        return None

    for mature_start in range(start_lo, start_hi + 1):
        for groove_end in range(end_lo, end_hi + 1):
            groove_len = groove_end - mature_start
            if not (groove_hlo <= groove_len <= groove_hhi):
                continue
            score = _score_class_ii_beta1_salvage(
                seq,
                mature_start,
                groove_end,
                beta1_c1,
                beta1_c2,
                beta1_groove_score,
                features=features,
            )
            if score > best_score:
                best_score = score
                best = (mature_start, groove_end, score)

    if best is None:
        return None
    return best


# ---------------------------------------------------------------------------
# Candidate enumeration
# ---------------------------------------------------------------------------


def _enumerate_mature_starts(
    scaffold: ParseScaffold,
    *,
    sp_estimate: int = 0,
) -> int:
    """Find the best mature_start by scoring all positions in [0, MAX_PLAUSIBLE_SP].

    When *sp_estimate* is provided (from ``infer_signal_peptide``), it is
    used to expand the search range and add a proximity bonus.  This
    prevents the Cys-pair-derived search range from excluding the true SP
    when the wrong Cys pair is selected.

    Returns the best-scoring position, or -1 if no candidate produces a
    plausible groove length.
    """
    groove_end = scaffold.groove_boundary
    _typical, _soft, hard = _GROOVE_PRIOR.get(scaffold.grammar.groove_prior_key, (90, (75, 100), (65, 115)))
    hlo, hhi = hard

    search_lo = max(0, groove_end - hhi)
    search_hi = min(MAX_PLAUSIBLE_SP, groove_end - max(hlo, 10))

    sp_confidence = 0.0
    sp_candidates: tuple[tuple[int, float], ...] = ()
    if scaffold.features is not None:
        if scaffold.features.sp_estimate == sp_estimate:
            sp_confidence = scaffold.features.sp_estimate_confidence
        sp_candidates = scaffold.features.sp_estimate_candidates
    if not sp_candidates and sp_estimate > 0:
        sp_candidates = ((sp_estimate, 1.0),)

    # Use the sequence-only SP estimate as a proposal anchor, not a hard
    # restriction. Even high-confidence estimates can be wrong on divergent
    # sequences, so we widen the search around them rather than clipping the
    # structurally valid range.
    for cand_pos, cand_score in sp_candidates:
        cand_conf = min(1.0, max(sp_confidence, float(cand_score) / 8.0))
        anchor_radius = max(4, int(round(12 - 8 * cand_conf)))
        search_lo = min(search_lo, max(0, cand_pos - anchor_radius))
        search_hi = max(search_hi, min(MAX_PLAUSIBLE_SP, cand_pos + anchor_radius))

    # Even when we expand toward the h-region estimate, the mature start must
    # stay upstream of the inferred groove boundary. Otherwise the candidate
    # implies a zero- or negative-length groove half and should not survive
    # enumeration.
    search_hi = min(search_hi, groove_end - 1)

    if search_hi < search_lo:
        return -1

    best_pos = -1
    best_score = -999.0

    for pos in range(search_lo, search_hi + 1):
        s = scaffold.score(pos)
        if sp_candidates:
            s += _score_sp_estimate_zone_consistency(pos, sp_candidates, sp_confidence)
        elif sp_estimate > 0:
            s += _score_sp_estimate_consistency(pos, sp_estimate, sp_confidence)
        if s > best_score + _SP_EARLIER_TIE_MARGIN:
            best_score = s
            best_pos = pos
        elif best_pos < 0 or (pos < best_pos and s >= best_score - _SP_EARLIER_TIE_MARGIN):
            best_score = s
            best_pos = pos

    # Also check pos=0 (SP-stripped / leaderless).
    # The pos=0 candidate gets its groove/ig/tail scores from the scaffold
    # normally, plus a leaderless bonus based on the same SP evidence used
    # above (Met, h-region, charged N-terminus).  This ensures pos=0 and
    # pos>0 compete on the same evidence model.
    s = scaffold.score(0)

    seq = scaffold.seq
    h_start, h_end = scaffold.h_region

    # --- Leaderless evidence (mirrors the joint SP penalty above) ---
    starts_met = seq[:1] == "M"
    h_frac = 0.0
    if h_end > h_start:
        h_frac = _hydrophobic_fraction(seq, h_start, h_end, features=scaffold.features)

    if not starts_met and h_frac < 0.55:
        s += 6.0  # no Met + no h-region = almost certainly leaderless
    elif not starts_met and h_frac < 0.65:
        s += 4.0  # no Met + weak h-region
    elif not starts_met:
        s += 2.0  # no Met but decent h-region (could be partial)
    elif sp_estimate == 0:
        s += 2.0  # has Met but infer_signal_peptide found nothing
    elif starts_met and h_frac >= 0.65 and sp_estimate > 0:
        # Anti-leaderless: strong SP evidence (Met + h-region + sp_estimate)
        # contradicts the leaderless hypothesis.  Penalize pos=0 so the
        # SP candidates from the main loop can win.
        s -= 4.0

    # N-terminus looks mature (charged/polar first residues)
    if len(seq) > 4:
        n_charged = sum(1 for aa in seq[:5] if aa in _SP_CHARGED)
        if n_charged >= 2:
            s += 1.0

    s += _score_nterm_lexical_state(seq, 0)

    if s > best_score + _SP_EARLIER_TIE_MARGIN:
        best_pos = 0
    elif best_pos < 0 or (0 < best_pos and s >= best_score - _SP_EARLIER_TIE_MARGIN):
        best_pos = 0

    return best_pos


# Signal peptide cleavage site residues.
#
# The -1 position (last residue of the SP) is A/G/S/C across all jawed
# vertebrates: sharks 96%, mammals 98%, birds 100%, fish 88%, reptiles 84%,
# amphibians 95% (including C at 25%).
#
def _score_refinement_candidate(
    seq: str,
    pos: int,
    initial_start: int,
    species_category: str = "",
    mhc_class: str = "",
    *,
    h_region: tuple[int, int] = (0, 0),
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Score a candidate SP cleavage site around an anchor prediction."""
    if pos < 3 or pos >= len(seq) - 3:
        return -999.0

    # Fast exclusion: impossible -3/-1 property combos
    if sp_boundary_excluded(seq, pos):
        return -10.0

    score = 0.0
    group = _refinement_group(species_category)
    sp_estimate = features.sp_estimate if features is not None else 0
    sp_confidence = features.sp_estimate_confidence if features is not None else 0.0
    sp_candidates = features.sp_estimate_candidates if features is not None else ()

    # Signal 1: -1 cleavage residue
    m1 = seq[pos - 1]
    if m1 in _SP_CLEAVAGE_STRONG:
        score += 3.0
    elif m1 == "S":
        score += 2.0
    elif m1 == "C":
        score += 1.5
    elif m1 == "T":
        score += 1.0
    elif m1 in _SP_CHARGED or m1 == "P":
        score -= 3.0
    else:
        score -= 0.5

    # Signal 2: -3 residue and canonical AxA cleavage motif.
    if seq[pos - 3] in _SP_SMALL_ALIPHATIC:
        score += 1.0
    elif seq[pos - 3] in _SP_CHARGED:
        score -= 1.0
    if m1 in "AGS" and seq[pos - 3] in _SP_SMALL_ALIPHATIC:
        score += 1.0

    # Signal 3: +1 not proline
    if seq[pos] == "P":
        score -= 2.0

    # Signal 4: upstream hydrophobic density (h-region)
    if pos >= 12:
        hydro_frac = _hydrophobic_fraction(seq, pos - 12, pos - 3, features=features)
        h_weight = 2.0 if group == "mammal" else 1.0
        score += (hydro_frac - 0.4) * h_weight
        if hydro_frac < 0.3:
            score -= 2.0

    # Signal 5: c-region polarity
    if pos >= 5:
        c_hydro = _hydrophobic_fraction(seq, pos - 3, pos, features=features)
        if c_hydro > 0.66:
            score -= 1.5

    # Signal 6: distance from the anchor prediction.
    delta_from_initial = pos - initial_start
    distance_from_initial = abs(delta_from_initial)
    if delta_from_initial > 0:
        distance_from_initial *= _SP_LATE_DRIFT_SCALE
    score -= distance_from_initial * _REFINEMENT_DISTANCE_PENALTY[group]
    if sp_candidates:
        score += _score_sp_estimate_zone_consistency(pos, sp_candidates, sp_confidence)
    else:
        score += _score_sp_estimate_consistency(pos, sp_estimate, sp_confidence)

    # Signal 7: mature protein start properties and category-specific motifs.
    if pos + 2 < len(seq):
        if seq[pos + 2] == "H":
            score += 2.0
        p1 = seq[pos]
        if p1 in "AGCVIL":
            score += 0.8
        elif p1 in "ST":
            score += 0.3
        elif p1 == "D":
            score -= 0.8
        elif p1 == "E":
            if _refinement_group(species_category) == "bird":
                score += 0.3
            else:
                score -= 0.2
    score += _score_refinement_motif(seq, pos, species_category, mhc_class)

    # Signal 8: charged mature-start density.
    if pos + 5 <= len(seq):
        charged_frac = _charged_fraction(seq, pos, pos + 5, features=features)
        if charged_frac > 0.3:
            score += 1.0
        elif charged_frac < 0.1:
            score -= 1.0

    # Signal 9: hydrophobicity transition at cleavage boundary.
    if pos >= 4 and pos + 3 <= len(seq):
        pre_kd = _kd_mean(seq, pos - 3, pos, features=features)
        post_kd = _kd_mean(seq, pos, pos + 3, features=features)
        kd_drop = pre_kd - post_kd
        if kd_drop > 0.5:
            score += min(kd_drop * 0.5, 1.5)

    # Signal 10: alignment with the signal-peptide hydrophobic core.
    score += _score_h_region_alignment(seq, pos, h_region=h_region) * _SP_H_REGION_WEIGHT

    # Signal 10b: MHC-specific statistical boundary words.
    score += 0.20 * _score_sp_boundary_words(seq, pos)

    # Signal 11: don't let the cleavage drift into the SP n-region.
    if seq[:1] == "M" and pos < 12:
        score -= 2.0

    return score


def refine_signal_peptide(
    sequence: str,
    mature_start: int,
    species_category: str = "",
    mhc_class: str = "",
    *,
    features: Optional["SequenceFeatures"] = None,
    groove_anchor: Optional[tuple[int, int]] = None,
) -> int:
    """Refine signal peptide cleavage site using sequence features.

    Parameters
    ----------
    sequence : str
        Full protein sequence (including signal peptide).
    mature_start : int
        Predicted mature protein start from Cys-pair heuristic.
    species_category : str
        One of the mhcseqs species categories.  Used for motif and
        distance-prior calibration.
    mhc_class : str
        Optional MHC class hint ("I" or "II") so class-specific mature
        N-terminal signals can participate in the score.
    features : SequenceFeatures, optional
        Pre-computed features from ``analyze_sequence()``.  Avoids
        redundant Cys-pair scanning and classification.
    groove_anchor : (c1, c2), optional
        Pre-identified groove anchor pair from the parser, used for
        junction detection without re-scanning.
    """
    if mature_start <= 0 or not sequence:
        return mature_start

    seq = features.seq if features else sequence.upper()
    if mature_start >= len(seq):
        return mature_start
    if (
        features is not None
        and features.sp_shortcut_kind == "exact_sp_prefix30"
        and features.sp_shortcut_state == "sp_present"
        and features.sp_shortcut_estimate > 0
        and abs(mature_start - features.sp_shortcut_estimate) <= 1
    ):
        return features.sp_shortcut_estimate

    group = _refinement_group(species_category)
    window = _REFINEMENT_WINDOW_BY_GROUP[group]

    # Try to find the groove Cys pair and junction position.
    junction_pos = 0
    if groove_anchor:
        junction_pos, _ = _find_junction(seq, groove_anchor[0])
    else:
        pairs = list(features.cys_pairs) if features else find_cys_pairs(seq)
        if pairs:
            annotated = list(features.pair_annotations) if features else [classify_cys_pair(seq, c1, c2) for c1, c2, _sep in pairs]
            groove_pairs = [a for a in annotated if a.domain_type in ("groove", "ambiguous")]
            if groove_pairs:
                groove = min(groove_pairs, key=lambda a: abs(a.c1 - (mature_start + 100)))
                junction_pos, _ = _find_junction(seq, groove.c1)

    # Compute h_region once for all refinement candidates.
    _h_region = features.h_region if features is not None else detect_h_region(seq)

    best_pos = mature_start
    initial_score = _score_refinement_candidate(
        seq,
        mature_start,
        mature_start,
        species_category,
        mhc_class,
        h_region=_h_region,
        features=features,
    )
    best_score = initial_score
    for pos in range(max(5, mature_start - window), min(len(seq) - 3, mature_start + window + 1)):
        s = _score_refinement_candidate(
            seq,
            pos,
            mature_start,
            species_category,
            mhc_class,
            h_region=_h_region,
            features=features,
        )
        if s > best_score + _SP_EARLIER_TIE_MARGIN:
            best_score = s
            best_pos = pos
        elif pos < best_pos and s >= best_score - _SP_EARLIER_TIE_MARGIN:
            best_score = s
            best_pos = pos

    junction_score = 0.0
    if junction_pos > 0:
        junction_score = _score_junction(seq, junction_pos)
    if junction_score >= 6.0:
        for pos in range(max(5, mature_start - 15), min(len(seq) - 3, mature_start + 16)):
            if mature_start - window <= pos <= mature_start + window:
                continue
            s = _score_refinement_candidate(
                seq,
                pos,
                mature_start,
                species_category,
                mhc_class,
                h_region=_h_region,
                features=features,
            )
            alpha1_len = junction_pos - pos
            if 80 <= alpha1_len <= 95:
                s += 2.5
            elif 75 <= alpha1_len <= 125:
                s += 1.0
            else:
                continue
            if s > best_score + _SP_EARLIER_TIE_MARGIN:
                best_score = s
                best_pos = pos
            elif pos < best_pos and s >= best_score - _SP_EARLIER_TIE_MARGIN:
                best_score = s
                best_pos = pos

    if best_pos != mature_start:
        move = abs(best_pos - mature_start)
        if move > 5:
            required_margin = 0.30 + 0.06 * move
            if best_pos > mature_start:
                required_margin += 0.10
            if best_score < initial_score + required_margin:
                return mature_start

    return best_pos


def _is_non_classical_class_i(gene: str, allele: str) -> bool:
    """Detect non-classical MHC class I lineages (teleost L/S/P/H, contaminants)."""
    token = str(gene or allele or "").strip().upper()
    # Strip species prefix
    if "-" in token:
        token = token.rsplit("-", 1)[-1]
    return any(token.startswith(pat) for pat in NON_CLASSICAL_CLASS_I_GENE_PATTERNS)


def _refine_status(result: AlleleRecord) -> AlleleRecord:
    """Post-parse refinement: detect non-classical lineage and short grooves.

    Applied after a successful parse to override status when the gene is
    non-classical or the groove is too short to be functional.
    """
    if not result.ok:
        return result

    flags = list(result.flags)

    # Non-classical class I detection
    if result.mhc_class == "I" and _is_non_classical_class_i(result.gene, result.allele):
        flags.append("non_classical_lineage")
        return replace(result, status="non_classical", flags=_flags_to_tuple(flags))

    grammar = _grammar_for_record(result)
    if grammar is None:
        return result

    # Short groove detection: whichever groove segment(s) this grammar
    # actually owns must still be large enough to be functional.
    for field_name, _role in grammar.groove_segments:
        seg_len = int(getattr(result, f"{field_name}_len", 0) or 0)
        if seg_len > 0 and seg_len < MIN_FUNCTIONAL_GROOVE_HALF_LEN:
            flags.append(f"{field_name}_short({seg_len})")
            return replace(result, status="short", flags=_flags_to_tuple(flags))

    return result


def _flags_to_tuple(flags: Sequence[str]) -> tuple[str, ...]:
    return tuple(str(flag) for flag in flags if flag)


def _slice_or_empty(seq: str, start: int, end: int) -> str:
    lo = max(0, int(start))
    hi = max(lo, min(len(seq), int(end)))
    return seq[lo:hi]


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


def apply_mutations(result: AlleleRecord, mutations: Sequence[object]) -> AlleleRecord:
    """Apply mutations to a AlleleRecord, returning a new result with mutated sequences.

    Mutations are specified in mature protein numbering (1-indexed).
    Accepts any sequence of mutation specs: strings ("K66A"), tuples, or
    mhcgnomes Mutation objects.

    The original groove extraction must have succeeded (result.ok must be True).
    Mutations are applied to the mature sequence reconstructed from domain parts,
    then re-sliced into the same domain boundaries.

    >>> from mhcseqs.domain_parsing import parse_class_i, apply_mutations
    >>> wt = decompose_class_i("G" * 400, allele="test")  # doctest: +SKIP
    """
    if not result.ok:
        raise ValueError(f"Cannot apply mutations to failed result (status={result.status!r})")

    parsed = [_parse_mutation(m) for m in mutations]
    if not parsed:
        return result

    # Reconstruct mature sequence from domain parts
    if result.mhc_class == "I":
        mature = result.groove1 + result.groove2 + result.ig_domain + result.tail
    elif result.chain == "alpha":
        mature = result.groove1 + result.ig_domain + result.tail
    else:  # beta
        mature = result.groove2 + result.ig_domain + result.tail

    # Apply mutations
    mature_list = list(mature)
    mutation_strs: list[str] = []
    for pos, orig, mutant in parsed:
        idx = pos - 1  # convert to 0-indexed
        if idx < 0 or idx >= len(mature_list):
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

    return replace(
        result,
        groove_seq=g1 + g2,
        groove1=g1,
        groove2=g2,
        ig_domain=ig,
        tail=t,
        tail_len=len(t),
        mutations=tuple(mutation_strs),
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


# ---------------------------------------------------------------------------
# Class II chain inference
# ---------------------------------------------------------------------------


def _class_ii_chain_from_name(
    *,
    gene: str,
    allele: str,
) -> Optional[str]:
    gene_token = str(gene or "").strip().upper()
    if not gene_token and allele:
        try:
            gene_token = infer_gene(allele)
        except Exception:
            gene_token = ""
        gene_token = str(gene_token or "").strip().upper()

    if not gene_token:
        return None
    if gene_token.startswith(CLASS_II_ALPHA_GENE_PREFIXES):
        return "alpha"
    if gene_token.startswith(CLASS_II_BETA_GENE_PREFIXES):
        return "beta"
    if gene_token.endswith("A"):
        return "alpha"
    if gene_token.endswith("B"):
        return "beta"
    return None


def _class_i_fragment_result(
    *,
    seq: str,
    allele: str,
    gene: str,
) -> AlleleRecord:
    """Return an alpha1_only or alpha2_only record for class I single-exon fragments.

    Distinguishes which groove half the fragment represents by checking for
    the α2 Ig-fold Cys pair (separation ~55–72):
      - Present → exon 3 / α2 domain → groove2
      - Absent  → exon 2 / α1 domain → groove1
    """
    cleaned = _clean_seq(seq)
    pairs = find_cys_pairs(cleaned)
    has_alpha2_pair = any(55 <= sep <= 72 for _, _, sep in pairs)

    if has_alpha2_pair:
        return _attach_parse_candidate(
            AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class="I",
                chain="alpha",
                sequence=cleaned,
                seq_len=len(cleaned),
                mature_start=0,
                groove_seq=cleaned,
                groove1="",
                groove2=cleaned,
                groove1_len=0,
                groove2_len=len(cleaned),
                status="alpha2_only",
                anchor_type="raw_fragment",
                flags=("single_exon_fragment",),
            )
        )
    return _attach_parse_candidate(
        AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="I",
            chain="alpha",
            sequence=cleaned,
            seq_len=len(cleaned),
            mature_start=0,
            groove_seq=cleaned,
            groove1=cleaned,
            groove2="",
            groove1_len=len(cleaned),
            groove2_len=0,
            status="alpha1_only",
            anchor_type="raw_fragment",
            flags=("single_exon_fragment",),
        )
    )


def _class_ii_fragment_result(
    *,
    seq: str,
    allele: str,
    gene: str,
    chain: str,
) -> AlleleRecord:
    cleaned = _clean_seq(seq)
    if chain == "alpha":
        return _attach_parse_candidate(
            AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class="II",
                chain="alpha",
                sequence=cleaned,
                seq_len=len(cleaned),
                mature_start=0,
                groove_seq=cleaned,
                groove1=cleaned,
                groove2="",
                groove1_len=len(cleaned),
                groove2_len=0,
                status="fragment_fallback",
                anchor_type="raw_fragment",
                flags=("fragment_fallback",),
            )
        )
    return _attach_parse_candidate(
        AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="beta",
            sequence=cleaned,
            seq_len=len(cleaned),
            mature_start=0,
            groove_seq=cleaned,
            groove1="",
            groove2=cleaned,
            groove1_len=0,
            groove2_len=len(cleaned),
            status="fragment_fallback",
            anchor_type="raw_fragment",
            flags=("fragment_fallback",),
        )
    )


def _missing_groove_result(
    *,
    seq_len: int,
    mhc_class: str,
    chain: str,
    sequence: str = "",
    allele: str = "",
    gene: str = "",
    mature_start: int = 0,
    anchor_type: str = "",
    anchor_cys1: Optional[int] = None,
    anchor_cys2: Optional[int] = None,
    flags: Sequence[str] = (),
) -> AlleleRecord:
    """Return a structured failure when no recoverable groove is present."""
    return _attach_parse_candidate(
        AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class=mhc_class,
            chain=chain,
            sequence=sequence,
            seq_len=seq_len,
            mature_start=mature_start,
            status="missing_groove",
            anchor_type=anchor_type,
            anchor_cys1=anchor_cys1,
            anchor_cys2=anchor_cys2,
            flags=_flags_to_tuple(flags),
        )
    )


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


def _prepare_parse_inputs(
    seq: str,
    *,
    grammar: GrammarSpec,
    allele: str,
    gene: str,
    features: Optional[SequenceFeatures] = None,
    use_early_shortcuts: bool = True,
) -> tuple[Optional[AlleleRecord], str, list[tuple[int, int, int]], list[CysPairAnnotation], int]:
    """Return shared parser inputs or an early terminal AlleleRecord."""
    cleaned = features.seq if features else _clean_seq(seq)
    if len(cleaned) < MIN_GROOVE_SOURCE_LEN:
        return (
            _attach_parse_candidate(
                AlleleRecord(
                    allele=allele,
                    gene=gene,
                    mhc_class=grammar.mhc_class,
                    chain=grammar.chain,
                    sequence=cleaned,
                    seq_len=len(cleaned),
                    status="too_short",
                )
            ),
            cleaned,
            [],
            [],
            0,
        )

    pairs = list(features.cys_pairs) if features else find_cys_pairs(cleaned)
    if not pairs:
        if len(cleaned) <= grammar.fragment_max_len:
            if grammar.mhc_class == "I":
                return _class_i_fragment_result(seq=cleaned, allele=allele, gene=gene), cleaned, [], [], 0
            return (
                _class_ii_fragment_result(seq=cleaned, allele=allele, gene=gene, chain=grammar.chain),
                cleaned,
                [],
                [],
                0,
            )
        return (
            _missing_groove_result(
                allele=allele,
                gene=gene,
                mhc_class=grammar.mhc_class,
                chain=grammar.chain,
                sequence=cleaned,
                seq_len=len(cleaned),
                flags=grammar.missing_no_pairs_flags,
            ),
            cleaned,
            [],
            [],
            0,
        )

    annotated = list(features.pair_annotations) if features else [classify_cys_pair(cleaned, c1, c2) for c1, c2, _sep in pairs]
    if features:
        sp_est = _selected_sp_estimate(features, use_early_shortcuts=use_early_shortcuts)
    else:
        h_region = detect_h_region(cleaned)
        if use_early_shortcuts:
            shortcut_est, _shortcut_conf, shortcut_state, shortcut_kind, _shortcut_candidates = fast_sp_triage(cleaned, h_region=h_region)
            if shortcut_kind == "exact_sp_prefix30" and shortcut_state == "sp_present":
                sp_est = shortcut_est
            elif shortcut_kind == "exact_mature10" and shortcut_state == "leaderless":
                sp_est = 0
            else:
                sp_est = infer_signal_peptide(cleaned, h_region=h_region)[0]
        else:
            sp_est = infer_signal_peptide(cleaned, h_region=h_region)[0]
    return None, cleaned, pairs, annotated, sp_est


def _build_class_i_secondary_resolver(
    seq: str,
    annotated: Sequence[CysPairAnnotation],
    features: Optional[SequenceFeatures] = None,
):
    """Return a resolver for the best downstream class-I C-like support pair."""
    ann_scored = sorted(
        [(a, a.ig_score + _cached_topology_support(seq, a.c1, a.c2, "c_like", features)) for a in annotated],
        key=lambda x: x[0].c1,
    )
    if not ann_scored:
        return lambda _ann: None
    suffix_best: list[tuple[CysPairAnnotation, float]] = [ann_scored[-1]] * len(ann_scored)
    for i in range(len(ann_scored) - 2, -1, -1):
        suffix_best[i] = ann_scored[i] if ann_scored[i][1] >= suffix_best[i + 1][1] else suffix_best[i + 1]
    c1_positions = [a.c1 for a, _ in ann_scored]

    def resolve(anchor: CysPairAnnotation) -> Optional[tuple[int, int]]:
        idx = bisect_right(c1_positions, anchor.c2 + 10)
        if idx < len(suffix_best):
            best = suffix_best[idx][0]
            return (best.c1, best.c2)
        return None

    return resolve


def _build_beta_secondary_resolver(
    seq: str,
    annotated: Sequence[CysPairAnnotation],
    features: Optional[SequenceFeatures] = None,
):
    """Return (resolver, strongest β1 annotation) for class-II beta parsing."""
    beta1_ann = (
        max(
            annotated,
            key=lambda a: (a.groove_score + _cached_topology_support(seq, a.c1, a.c2, "g_domain", features)) - (0.01 * a.c1),
        )
        if annotated
        else None
    )

    def resolve(anchor: CysPairAnnotation) -> Optional[tuple[int, int]]:
        if beta1_ann is None:
            return None
        if not _grammar_spec("II", "beta").allow_self_secondary and (anchor.c1, anchor.c2) == (beta1_ann.c1, beta1_ann.c2):
            return None
        return (beta1_ann.c1, beta1_ann.c2)

    return resolve, beta1_ann


def _tail_state_for_candidate_tail(
    seq_len: int,
    tail_start: int,
    tail: str,
    *,
    features: Optional[SequenceFeatures] = None,
) -> str:
    """Return the tail/TM state for a concrete tail segment."""
    if not tail:
        return "tail_absent"
    if features is not None:
        tm_start, tm_end = features.tm_span
        if tm_end > tm_start and tm_start >= tail_start and tm_end <= seq_len:
            return "tm_present"
    tm_start, tm_end = _find_tm_span(tail, raw_start=0)
    if tm_end > tm_start:
        return "tm_present"
    return "tail_no_tm"


def _scaled_contradiction_penalty(
    expected_score: float,
    alternate_score: float,
    *,
    structural_mass: float,
    scale: float,
    label: str,
) -> tuple[float, tuple[str, ...]]:
    """Return a structural contradiction penalty when the wrong topology wins."""
    if alternate_score <= expected_score:
        return 0.0, ()
    gap = alternate_score - expected_score
    divisor = _PRIMARY_CANDIDATE_CONTRADICTION_WEIGHTS["structural_mass_divisor"]
    mass_factor = 1.0 + (max(0.0, structural_mass) / max(divisor, 1.0))
    return -(gap * scale * mass_factor), (label,)


def _score_primary_contradictions(
    seq: str,
    *,
    grammar: GrammarSpec,
    mature_start: int,
    scaffold: ParseScaffold,
    domains: ExtractedDomains,
    features: Optional[SequenceFeatures] = None,
) -> tuple[float, tuple[str, ...]]:
    """Score structural contradictions that should not be rescued additively.

    Strong local groove-length or anchor scores should not compensate for
    topology that points at the wrong domain role, or for a full-length-like
    sequence missing expected support. These penalties scale with the size of
    the structural score so obvious contradictions do not survive on the back
    of one dominant additive term.
    """
    score = 0.0
    evidence: list[str] = []
    structural_mass = max(0.0, scaffold.fixed_groove) + 0.5 * max(0.0, scaffold.fixed_ig) + 0.25 * max(0.0, scaffold.fixed_tail)

    primary_scale = _PRIMARY_CANDIDATE_CONTRADICTION_WEIGHTS["primary_topology_scale"]
    secondary_scale = _PRIMARY_CANDIDATE_CONTRADICTION_WEIGHTS["secondary_topology_scale"]
    support_missing_scale = _PRIMARY_CANDIDATE_CONTRADICTION_WEIGHTS["support_missing_scale"]
    sp_zone_scale = _PRIMARY_CANDIDATE_CONTRADICTION_WEIGHTS["sp_zone_conflict_scale"]

    if grammar.uses_primary_anchor_for_ig:
        delta, ev = _scaled_contradiction_penalty(
            scaffold.anchor_c_topology,
            scaffold.anchor_g_topology,
            structural_mass=structural_mass,
            scale=primary_scale,
            label="primary_anchor_topology_conflict",
        )
    else:
        delta, ev = _scaled_contradiction_penalty(
            scaffold.anchor_g_topology,
            scaffold.anchor_c_topology,
            structural_mass=structural_mass,
            scale=primary_scale,
            label="primary_anchor_topology_conflict",
        )
    score += delta
    evidence.extend(ev)

    if grammar.mhc_class == "I" and scaffold.secondary_c_topology > 0.0:
        delta, ev = _scaled_contradiction_penalty(
            scaffold.secondary_c_topology,
            scaffold.secondary_g_topology,
            structural_mass=structural_mass,
            scale=secondary_scale,
            label="support_anchor_topology_conflict",
        )
        score += delta
        evidence.extend(ev)
    elif grammar.chain == "beta" and scaffold.secondary_g_topology > 0.0:
        delta, ev = _scaled_contradiction_penalty(
            scaffold.secondary_g_topology,
            scaffold.secondary_c_topology,
            structural_mass=structural_mass,
            scale=secondary_scale,
            label="secondary_groove_topology_conflict",
        )
        score += delta
        evidence.extend(ev)

    if not domains.ig_domain and len(seq) > grammar.fragment_max_len:
        score -= max(0.0, structural_mass) * support_missing_scale
        evidence.append("support_missing_scaled_by_structure")

    if features is not None and features.sp_estimate_candidates:
        zone_score = _score_sp_estimate_zone_consistency(
            mature_start,
            features.sp_estimate_candidates,
            features.sp_estimate_confidence,
        )
        if zone_score < 0.0:
            divisor = _PRIMARY_CANDIDATE_CONTRADICTION_WEIGHTS["structural_mass_divisor"]
            mass_factor = 1.0 + (max(0.0, scaffold.fixed_groove) / max(divisor, 1.0))
            score += zone_score * mass_factor * sp_zone_scale
            evidence.append("sp_zone_conflict")

    return score, tuple(evidence)


def _score_primary_candidate_selection(
    seq: str,
    *,
    grammar: GrammarSpec,
    mature_start: int,
    anchor_pair: tuple[int, int],
    secondary_pair: Optional[tuple[int, int]],
    scaffold: ParseScaffold,
    sp_estimate: int,
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Score one full-length primary candidate using state compatibility."""
    domains = _extract_primary_domains(
        seq,
        grammar=grammar,
        mature_start=mature_start,
        anchor_pair=anchor_pair,
        secondary_pair=secondary_pair,
    )
    candidate_type = _full_candidate_type(grammar) if domains.ig_domain else _missing_support_candidate_type(grammar)
    nterm_state = "sp_present" if mature_start > 0 else "leaderless_mature_only"
    support_state = "support_present" if domains.ig_domain else "support_missing"
    tail_start = len(seq) - len(domains.tail)
    tail_state = _tail_state_for_candidate_tail(
        len(seq),
        tail_start,
        domains.tail,
        features=features,
    )
    compatibility_score, _evidence = _candidate_compatibility_score(
        candidate_type=candidate_type,
        nterm_state=nterm_state,
        support_state=support_state,
        tail_state=tail_state,
        seq_len=len(seq),
        mature_start=mature_start,
        seq=seq,
        h_region=scaffold.h_region,
        sp_estimate=sp_estimate,
        features=features,
    )
    contradiction_score, _evidence = _score_primary_contradictions(
        seq,
        grammar=grammar,
        mature_start=mature_start,
        scaffold=scaffold,
        domains=domains,
        features=features,
    )
    return scaffold.score(mature_start) + compatibility_score + contradiction_score


def _primary_selection_key(selection: PrimaryParseSelection) -> tuple[float, float, float, float]:
    """Return a stable ordering key for primary parse candidates."""
    secondary_bonus = 1.0 if selection.secondary_pair is not None else 0.0
    return (
        float(selection.selection_score),
        secondary_bonus,
        float(selection.anchor.groove_score + selection.anchor.ig_score),
        float(-selection.anchor.c1),
    )


def _primary_selection_identity(selection: PrimaryParseSelection) -> tuple[int, int, int, int, int]:
    """Return a deduplication key for primary parse candidates."""
    sec_c1 = selection.secondary_pair[0] if selection.secondary_pair is not None else -1
    sec_c2 = selection.secondary_pair[1] if selection.secondary_pair is not None else -1
    return (
        int(selection.mature_start),
        int(selection.anchor.c1),
        int(selection.anchor.c2),
        int(sec_c1),
        int(sec_c2),
    )


def _primary_selection_is_distinct(
    selection: PrimaryParseSelection,
    kept: Sequence[PrimaryParseSelection],
) -> bool:
    """Return True when a primary candidate adds a new structural explanation."""
    if not kept:
        return True
    selection_has_secondary = selection.secondary_pair is not None
    for other in kept:
        other_has_secondary = other.secondary_pair is not None
        if selection_has_secondary != other_has_secondary:
            continue
        if abs(selection.anchor.c1 - other.anchor.c1) >= 6:
            continue
        if abs(selection.mature_start - other.mature_start) >= 4:
            continue
        return False
    return True


def _record_has_short_flag(record: AlleleRecord) -> bool:
    """Whether a record already looks fragment-like or structurally compressed."""
    return any(flag.endswith("_short") or "_short(" in flag for flag in record.flags)


def _should_expand_primary_competition(
    selections: Sequence[PrimaryParseSelection],
    best_record: Optional[AlleleRecord],
    *,
    grammar: GrammarSpec,
    features: Optional[SequenceFeatures] = None,
    use_early_shortcuts: bool = False,
) -> bool:
    """Whether to broaden competition beyond the top primary explanation.

    Most sequences have a decisive single structural explanation. Expansion
    should only happen when the best primary parse is weak, incomplete, or
    genuinely close to a competitor.
    """
    if best_record is None:
        return True
    if best_record.status != "ok":
        return True
    if _record_has_short_flag(best_record):
        return True

    margin = float("inf")
    if len(selections) >= 2:
        margin = selections[0].selection_score - selections[1].selection_score

    if use_early_shortcuts and features is not None:
        if (
            features.sp_shortcut_kind == "exact_sp_prefix30"
            and features.sp_shortcut_state == "sp_present"
            and best_record.mature_start > 0
            and abs(best_record.mature_start - features.sp_shortcut_estimate) <= 1
            and best_record.parse_candidate.support_state == "support_present"
        ):
            return False
        if features.sp_shortcut_kind == "exact_mature10" and features.sp_shortcut_state == "leaderless" and best_record.mature_start == 0:
            return False
        if (
            features.sp_shortcut_state == "sp_present"
            and features.sp_shortcut_confidence >= 0.90
            and best_record.mature_start > 0
            and abs(best_record.mature_start - features.sp_shortcut_estimate) <= 2
            and margin >= max(4.0, PRIMARY_PARSE_EXPANSION_MARGIN / 2.0)
        ):
            return False
        if (
            features.sp_shortcut_state == "leaderless"
            and features.sp_shortcut_confidence >= 0.90
            and best_record.mature_start == 0
            and margin >= max(4.0, PRIMARY_PARSE_EXPANSION_MARGIN / 2.0)
        ):
            return False

    candidate = best_record.parse_candidate
    if candidate.support_state in {"support_missing", "support_inferred"} and best_record.seq_len > grammar.fragment_max_len:
        return True

    if features is not None:
        if features.sp_estimate_confidence >= PRIMARY_PARSE_LOW_CONFIDENCE and features.sp_estimate > 0 and best_record.mature_start > 0:
            if abs(best_record.mature_start - features.sp_estimate) > 6:
                return True

    if len(selections) < 2:
        return False
    return margin < PRIMARY_PARSE_EXPANSION_MARGIN


def _collect_primary_parse_candidates(
    seq: str,
    *,
    grammar: GrammarSpec,
    annotated: Sequence[CysPairAnnotation],
    sp_estimate: int,
    features: Optional[SequenceFeatures] = None,
    secondary_resolver=None,
) -> list[PrimaryParseSelection]:
    """Enumerate top whole-parse primary candidates under one grammar.

    The parser used to collapse each grammar down to one primary full parse
    before partial, salvage, and missing-evidence candidates could compete.
    Keeping a small top-K set lets structurally distinct anchor hypotheses
    survive into whole-parse ranking, which is the main remaining source of
    generalizable accuracy gains.
    """
    selections: list[PrimaryParseSelection] = []
    resolve_secondary = secondary_resolver or (lambda _ann: None)

    for ann in annotated:
        secondary = resolve_secondary(ann)
        sec_c1 = secondary[0] if secondary else 0
        sec_c2 = secondary[1] if secondary else 0
        scaffold = _build_parse_scaffold(
            seq,
            ann.c1,
            ann.c2,
            grammar.mhc_class,
            grammar.chain,
            sec_c1,
            sec_c2,
            ann.groove_score,
            ann.ig_score,
            features=features,
        )
        ms = _enumerate_mature_starts(scaffold, sp_estimate=sp_estimate)
        if ms < 0:
            continue
        score = _score_primary_candidate_selection(
            seq,
            grammar=grammar,
            mature_start=ms,
            anchor_pair=(ann.c1, ann.c2),
            secondary_pair=secondary,
            scaffold=scaffold,
            sp_estimate=sp_estimate,
            features=features,
        )
        selections.append(
            PrimaryParseSelection(
                mature_start=ms,
                anchor=ann,
                secondary_pair=secondary,
                scaffold=scaffold,
                selection_score=score,
            )
        )

    ranked = sorted(selections, key=_primary_selection_key, reverse=True)
    deduped: list[PrimaryParseSelection] = []
    seen: set[tuple[int, int, int, int, int]] = set()
    for selection in ranked:
        identity = _primary_selection_identity(selection)
        if identity in seen:
            continue
        if not _primary_selection_is_distinct(selection, deduped):
            continue
        seen.add(identity)
        deduped.append(selection)
        if len(deduped) >= PRIMARY_PARSE_CANDIDATE_KEEP:
            break
    return deduped


def _select_best_primary_parse(
    seq: str,
    *,
    grammar: GrammarSpec,
    annotated: Sequence[CysPairAnnotation],
    sp_estimate: int,
    features: Optional[SequenceFeatures] = None,
    secondary_resolver=None,
) -> Optional[PrimaryParseSelection]:
    """Return the best primary parse under one grammar."""
    selections = _collect_primary_parse_candidates(
        seq,
        grammar=grammar,
        annotated=annotated,
        sp_estimate=sp_estimate,
        features=features,
        secondary_resolver=secondary_resolver,
    )
    if not selections:
        return None
    return selections[0]


def _extract_primary_domains(
    seq: str,
    *,
    grammar: GrammarSpec,
    mature_start: int,
    anchor_pair: tuple[int, int],
    secondary_pair: Optional[tuple[int, int]] = None,
) -> ExtractedDomains:
    """Materialize groove / support / tail slices from a grammar-aligned parse."""
    anchor_c1, anchor_c2 = anchor_pair
    groove_boundary = max(mature_start, anchor_c1 - grammar.groove_boundary_from_c1)
    groove_end = groove_boundary if grammar.groove_end_from_c2 == 0 else min(len(seq), anchor_c2 + grammar.groove_end_from_c2)

    if grammar.groove_output == "split_both":
        groove1 = _slice_or_empty(seq, mature_start, groove_boundary)
        groove2 = _slice_or_empty(seq, groove_boundary, groove_end)
        groove_seq = groove1 + groove2
        if secondary_pair:
            ig_end = min(len(seq), secondary_pair[1] + IG_DOMAIN_END_AFTER_CYS2)
            ig_domain = _slice_or_empty(seq, groove_end, ig_end)
            tail = seq[ig_end:]
        else:
            ig_domain = ""
            tail = seq[groove_end:]
        return ExtractedDomains(groove1=groove1, groove2=groove2, groove_seq=groove_seq, ig_domain=ig_domain, tail=tail)

    groove_segment = _slice_or_empty(seq, mature_start, groove_boundary)
    ig_end = min(len(seq), anchor_c2 + IG_DOMAIN_END_AFTER_CYS2)
    ig_domain = _slice_or_empty(seq, groove_boundary, ig_end)
    tail = seq[ig_end:]
    if grammar.groove_output == "groove1":
        return ExtractedDomains(groove1=groove_segment, groove2="", groove_seq=groove_segment, ig_domain=ig_domain, tail=tail)
    return ExtractedDomains(groove1="", groove2=groove_segment, groove_seq=groove_segment, ig_domain=ig_domain, tail=tail)


def _apply_short_flags(grammar: GrammarSpec, domains: ExtractedDomains, flags: list[str]) -> None:
    """Append grammar-specific short-domain flags."""
    if grammar.groove_output == "split_both":
        lengths = (len(domains.groove1), len(domains.groove2))
    elif grammar.groove_output == "groove1":
        lengths = (len(domains.groove1),)
    else:
        lengths = (len(domains.groove2),)
    for (flag_name, min_len), actual_len in zip(grammar.short_flags, lengths):
        if actual_len < min_len:
            flags.append(f"{flag_name}({actual_len})")


def _append_long_sp_flag(flags: list[str], mature_start: int) -> None:
    """Record an SP-length warning when the mature start looks unusually deep."""
    if mature_start > 40:
        flags.append(f"long_sp({mature_start})")


def _domains_match_grammar(grammar: GrammarSpec, domains: ExtractedDomains) -> bool:
    """Return True when the materialized groove segments satisfy the grammar."""
    return all(bool(getattr(domains, field_name)) for field_name, _role in grammar.groove_segments)


def _build_primary_result(
    *,
    allele: str,
    gene: str,
    seq: str,
    grammar: GrammarSpec,
    selection: PrimaryParseSelection,
    sp_estimate: int = 0,
    features: Optional[SequenceFeatures] = None,
    flags: Sequence[str] = (),
    secondary_pair: Optional[tuple[int, int]] = None,
) -> AlleleRecord:
    """Materialize one grammar-driven primary parse into an AlleleRecord."""
    local_flags = list(flags)
    mature_start = selection.mature_start
    anchor_pair = (selection.anchor.c1, selection.anchor.c2)
    _append_long_sp_flag(local_flags, mature_start)
    domains = _extract_primary_domains(
        seq,
        grammar=grammar,
        mature_start=mature_start,
        anchor_pair=anchor_pair,
        secondary_pair=selection.secondary_pair,
    )
    _apply_short_flags(grammar, domains, local_flags)
    if not _domains_match_grammar(grammar, domains):
        return _attach_parse_candidate(
            AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class=grammar.mhc_class,
                chain=grammar.chain,
                sequence=seq,
                seq_len=len(seq),
                mature_start=mature_start,
                status="invalid_boundaries",
                anchor_type=grammar.anchor_type,
                anchor_cys1=anchor_pair[0],
                anchor_cys2=anchor_pair[1],
            ),
            seq=seq,
            h_region=selection.scaffold.h_region,
            sp_estimate=sp_estimate,
            features=features,
        )
    candidate_type = _full_candidate_type(grammar) if domains.ig_domain else _missing_support_candidate_type(grammar)
    if not domains.ig_domain and "missing_support_anchor" not in local_flags:
        local_flags.append("missing_support_anchor")
    components = selection.scaffold.score_components(mature_start)
    record_secondary = secondary_pair if secondary_pair is not None else selection.secondary_pair
    return _attach_parse_candidate(
        AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class=grammar.mhc_class,
            chain=grammar.chain,
            sequence=seq,
            seq_len=len(seq),
            mature_start=mature_start,
            groove_seq=domains.groove_seq,
            groove1=domains.groove1,
            groove2=domains.groove2,
            groove1_len=len(domains.groove1),
            groove2_len=len(domains.groove2),
            ig_domain=domains.ig_domain,
            ig_domain_len=len(domains.ig_domain),
            tail=domains.tail,
            tail_len=len(domains.tail),
            parse_score=components.total,
            sp_subscore=components.sp,
            groove_subscore=components.groove,
            ig_subscore=components.ig,
            tail_subscore=components.tail,
            status="ok",
            anchor_type=grammar.anchor_type,
            anchor_cys1=anchor_pair[0],
            anchor_cys2=anchor_pair[1],
            secondary_cys1=(record_secondary[0] if record_secondary else None),
            secondary_cys2=(record_secondary[1] if record_secondary else None),
            candidate_type=candidate_type,
            flags=_flags_to_tuple(local_flags),
        ),
        seq=seq,
        h_region=selection.scaffold.h_region,
        sp_estimate=sp_estimate,
        features=features,
    )


def _build_class_i_alpha3_salvage_result(
    *,
    seq: str,
    allele: str,
    gene: str,
    alpha3_ann: Optional[CysPairAnnotation],
    features: Optional[SequenceFeatures] = None,
) -> Optional[AlleleRecord]:
    """Build a support-anchored class-I salvage candidate when α2 is missing."""
    if alpha3_ann is None:
        return None
    c1, c2 = alpha3_ann.c1, alpha3_ann.c2
    salvage = _enumerate_class_i_from_alpha3(seq, c1, c2, features=features)
    if salvage is None:
        return None
    mature_start, alpha2_start, alpha3_start, salvage_score = salvage
    half_1 = _slice_or_empty(seq, mature_start, alpha2_start)
    half_2 = _slice_or_empty(seq, alpha2_start, alpha3_start)
    ig_end = min(len(seq), c2 + IG_DOMAIN_END_AFTER_CYS2)
    ig = _slice_or_empty(seq, alpha3_start, ig_end)
    tail = seq[ig_end:]
    flags = ["inferred_from_alpha3"]
    _append_long_sp_flag(flags, mature_start)
    if len(half_1) < 50:
        flags.append(f"alpha1_short({len(half_1)})")
    if len(half_2) < 60:
        flags.append(f"alpha2_short({len(half_2)})")
    return _attach_parse_candidate(
        AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="I",
            chain="alpha",
            sequence=seq,
            seq_len=len(seq),
            mature_start=mature_start,
            groove_seq=half_1 + half_2,
            groove1=half_1,
            groove2=half_2,
            groove1_len=len(half_1),
            groove2_len=len(half_2),
            ig_domain=ig,
            ig_domain_len=len(ig),
            tail=tail,
            tail_len=len(tail),
            parse_score=salvage_score,
            status="inferred_from_alpha3",
            anchor_type="alpha3_cys",
            anchor_cys1=c1,
            anchor_cys2=c2,
            flags=_flags_to_tuple(flags),
        ),
        seq=seq,
        h_region=(features.h_region if features is not None else None),
        sp_estimate=(features.sp_estimate if features is not None else None),
        features=features,
    )


def _build_class_ii_beta1_salvage_result(
    *,
    seq: str,
    allele: str,
    gene: str,
    beta1_ann: Optional[CysPairAnnotation],
    features: Optional[SequenceFeatures] = None,
) -> Optional[AlleleRecord]:
    """Build a β1-only class-II beta candidate when β2 support is missing."""
    if beta1_ann is None:
        return None
    c1, c2 = beta1_ann.c1, beta1_ann.c2
    salvage = _enumerate_class_ii_beta_from_beta1(
        seq,
        c1,
        c2,
        beta1_ann.groove_score,
        features=features,
    )
    if salvage is None:
        return None
    mature_start, groove_end, salvage_score = salvage
    flags = ["beta1_only_fallback"]
    _append_long_sp_flag(flags, mature_start)
    domains = ExtractedDomains(
        groove1="",
        groove2=_slice_or_empty(seq, mature_start, groove_end),
        groove_seq=_slice_or_empty(seq, mature_start, groove_end),
        ig_domain="",
        tail=seq[groove_end:],
    )
    _apply_short_flags(_grammar_spec("II", "beta"), domains, flags)
    if not domains.groove2:
        return None
    components = ParseSubscores(0.0, 0.0, 0.0, 0.0)
    return _attach_parse_candidate(
        AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="beta",
            sequence=seq,
            seq_len=len(seq),
            mature_start=mature_start,
            groove_seq=domains.groove_seq,
            groove1="",
            groove2=domains.groove2,
            groove1_len=0,
            groove2_len=len(domains.groove2),
            ig_domain=domains.ig_domain,
            ig_domain_len=len(domains.ig_domain),
            tail=domains.tail,
            tail_len=len(domains.tail),
            parse_score=salvage_score,
            sp_subscore=components.sp,
            groove_subscore=components.groove,
            ig_subscore=components.ig,
            tail_subscore=components.tail,
            status="beta1_only_fallback",
            anchor_type="beta1_cys",
            anchor_cys1=c1,
            anchor_cys2=c2,
            flags=_flags_to_tuple(flags),
        ),
        seq=seq,
        h_region=(features.h_region if features is not None else None),
        sp_estimate=(features.sp_estimate if features is not None else None),
        features=features,
    )


# ---------------------------------------------------------------------------
# Class I groove parser
# ---------------------------------------------------------------------------


def decompose_class_i(
    seq: str,
    *,
    allele: str = "",
    gene: str = "",
    features: Optional[SequenceFeatures] = None,
    use_early_shortcuts: bool = True,
) -> AlleleRecord:
    """Parse a class-I alpha chain via the shared grammar-driven engine.

    Enumerates a small top-K set of structurally distinct primary parses over
    plausible groove Cys pairs and mature_start positions, then lets those
    full parses compete with fragment, salvage, and missing-groove
    explanations under one candidate selector.

    Pass *features* from ``analyze_sequence()`` to skip redundant
    Cys-pair scanning and fold classification.
    """
    grammar = _grammar_spec("I", "alpha")
    flags: list[str] = []
    early_result, cleaned, _pairs, annotated, sp_est = _prepare_parse_inputs(
        seq,
        grammar=grammar,
        allele=allele,
        gene=gene,
        features=features,
        use_early_shortcuts=use_early_shortcuts,
    )
    if early_result is not None:
        return early_result

    candidates: list[AlleleRecord] = []
    selections = _collect_primary_parse_candidates(
        cleaned,
        grammar=grammar,
        annotated=annotated,
        sp_estimate=sp_est,
        features=features,
        secondary_resolver=_build_class_i_secondary_resolver(cleaned, annotated, features),
    )
    best_record: Optional[AlleleRecord] = None
    expand_primary = True
    if selections:
        best_record = _build_primary_result(
            allele=allele,
            gene=gene,
            seq=cleaned,
            grammar=grammar,
            selection=selections[0],
            sp_estimate=sp_est,
            features=features,
            flags=flags,
        )
        _append_record_candidate(candidates, best_record)
        expand_primary = _should_expand_primary_competition(
            selections,
            best_record,
            grammar=grammar,
            features=features,
            use_early_shortcuts=use_early_shortcuts,
        )
        if expand_primary:
            for selection in selections[1:]:
                record = _build_primary_result(
                    allele=allele,
                    gene=gene,
                    seq=cleaned,
                    grammar=grammar,
                    selection=selection,
                    sp_estimate=sp_est,
                    features=features,
                    flags=flags,
                )
                _append_record_candidate(candidates, record)

    if len(cleaned) <= grammar.fragment_max_len:
        _append_record_candidate(
            candidates,
            _class_i_fragment_result(seq=cleaned, allele=allele, gene=gene),
        )

    alpha3_ann = (
        max(
            annotated,
            key=lambda a: a.ig_score + _cached_topology_support(cleaned, a.c1, a.c2, "c_like", features),
        )
        if annotated
        else None
    )
    if alpha3_ann is not None:
        if best_record is None or best_record.status != "ok" or best_record.parse_candidate.support_state != "support_present":
            _append_record_candidate(
                candidates,
                _build_class_i_alpha3_salvage_result(
                    seq=cleaned,
                    allele=allele,
                    gene=gene,
                    alpha3_ann=alpha3_ann,
                    features=features,
                ),
            )

    missing_flags: tuple[str, ...]
    if alpha3_ann is None:
        missing_flags = ("no_plausible_groove_anchor", "no_plausible_support_anchor")
    else:
        missing_flags = ("no_plausible_groove_anchor", "alpha3_salvage_failed")
    _append_record_candidate(
        candidates,
        _missing_groove_result(
            allele=allele,
            gene=gene,
            mhc_class="I",
            chain="alpha",
            sequence=cleaned,
            seq_len=len(cleaned),
            anchor_type=("alpha3_cys" if alpha3_ann is not None else ""),
            anchor_cys1=(alpha3_ann.c1 if alpha3_ann is not None else None),
            anchor_cys2=(alpha3_ann.c2 if alpha3_ann is not None else None),
            flags=missing_flags,
        ),
    )
    return _select_best_record_candidate(candidates)


def trace_parse_class_i(
    seq: str,
    *,
    allele: str = "",
    gene: str = "",
) -> tuple[AlleleRecord, ParseTrace]:
    """Parse class I with full debug trace of scoring decisions.

    Returns (AlleleRecord, ParseTrace) — the trace contains all Cys pair
    classifications, junction scores, SP candidate scores, and boundary
    motif scores for debugging ground truth disagreements.
    """
    cleaned = _clean_seq(seq)
    trace = ParseTrace(sequence_len=len(cleaned), mhc_class="I")

    # Step 1: find and classify Cys pairs
    pairs = find_cys_pairs(cleaned)
    if pairs:
        annotated = tuple(classify_cys_pair(cleaned, c1, c2) for c1, c2, _sep in pairs)
        trace.cys_pairs = annotated
        for i, ann in enumerate(annotated):
            trace.log.append(
                f"Cys #{i + 1} C{ann.c1}-C{ann.c2} sep={ann.separation} groove={ann.groove_score:.1f} ig={ann.ig_score:.1f} → {ann.domain_type}"
            )
            trace.log.append(f"  evidence: {' '.join(ann.evidence)}")

    # Step 2: find groove pair and junction
    groove_classified = [a for a in trace.cys_pairs if a.domain_type in ("groove", "ambiguous")]
    if groove_classified:
        best_groove = max(groove_classified, key=lambda a: a.groove_score)
        jpos, jscore = _find_junction(cleaned, best_groove.c1)
        trace.junction_pos = jpos
        trace.junction_score = jscore
        trace.log.append(f"Junction search near C{best_groove.c1}-10={best_groove.c1 - 10}: best at {jpos} score={jscore:.1f}")
        # Show top 3 junction candidates
        j_candidates = []
        for p in range(max(0, best_groove.c1 - 15), min(len(cleaned), best_groove.c1 - 5)):
            j_candidates.append((p, _score_junction(cleaned, p)))
        j_candidates.sort(key=lambda x: -x[1])
        for p, s in j_candidates[:3]:
            motif = cleaned[p : p + 6] if p + 6 <= len(cleaned) else cleaned[p:]
            trace.log.append(f"  junction candidate pos={p} motif={motif} score={s:.1f}")

    # Step 3: run the actual parser
    result = decompose_class_i(cleaned, allele=allele, gene=gene)
    trace.selected_sp = result.mature_start

    # Step 4: log domain boundaries
    trace.log.append(f"Parse result: status={result.status} mature_start={result.mature_start}")
    if result.ok:
        trace.log.append(
            f"  [SP 0..{result.mature_start}] "
            f"[groove1 {result.mature_start}..{result.mature_start + result.groove1_len}] "
            f"[groove2 ..+{result.groove2_len}] "
            f"[Ig ..+{result.ig_domain_len}] "
            f"[tail ..+{result.tail_len}]"
        )

    # Step 5: score boundary motifs for the trace
    if result.ok and result.anchor_cys2 is not None:
        boundary_pos = result.anchor_cys2 + CLASS_I_ALPHA2_END_AFTER_CYS2
        bscore = _score_groove_ig_boundary(cleaned, boundary_pos, "I")
        trace.boundary_scores["groove2_ig"] = bscore
        trace.log.append(f"Boundary groove2|Ig @{boundary_pos}: score={bscore:.1f}")

    return result, trace


# ---------------------------------------------------------------------------
# Class II alpha groove parser
# ---------------------------------------------------------------------------


def decompose_class_ii_alpha(
    seq: str,
    *,
    allele: str = "",
    gene: str = "",
    features: Optional[SequenceFeatures] = None,
    use_early_shortcuts: bool = True,
) -> AlleleRecord:
    """Parse a class-II alpha chain via the shared grammar-driven engine.

    Enumerates a small top-K set of primary full parses, then lets them
    compete with fragment and missing-groove explanations.
    """
    grammar = _grammar_spec("II", "alpha")
    flags: list[str] = []
    early_result, cleaned, _pairs, annotated, sp_est = _prepare_parse_inputs(
        seq,
        grammar=grammar,
        allele=allele,
        gene=gene,
        features=features,
        use_early_shortcuts=use_early_shortcuts,
    )
    if early_result is not None:
        return early_result

    candidates: list[AlleleRecord] = []
    selections = _collect_primary_parse_candidates(
        cleaned,
        grammar=grammar,
        annotated=annotated,
        sp_estimate=sp_est,
        features=features,
    )
    best_record = None
    expand_primary = True
    if selections:
        best_record = _build_primary_result(
            allele=allele,
            gene=gene,
            seq=cleaned,
            grammar=grammar,
            selection=selections[0],
            sp_estimate=sp_est,
            features=features,
            flags=flags,
        )
        _append_record_candidate(candidates, best_record)
        expand_primary = _should_expand_primary_competition(
            selections,
            best_record,
            grammar=grammar,
            features=features,
            use_early_shortcuts=use_early_shortcuts,
        )
        if expand_primary:
            for selection in selections[1:]:
                _append_record_candidate(
                    candidates,
                    _build_primary_result(
                        allele=allele,
                        gene=gene,
                        seq=cleaned,
                        grammar=grammar,
                        selection=selection,
                        sp_estimate=sp_est,
                        features=features,
                        flags=flags,
                    ),
                )

    if len(cleaned) <= grammar.fragment_max_len:
        _append_record_candidate(
            candidates,
            _class_ii_fragment_result(seq=cleaned, allele=allele, gene=gene, chain="alpha"),
        )

    _append_record_candidate(
        candidates,
        _missing_groove_result(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="alpha",
            sequence=cleaned,
            seq_len=len(cleaned),
            flags=grammar.missing_no_anchor_flags,
        ),
    )
    return _select_best_record_candidate(candidates)


# ---------------------------------------------------------------------------
# Class II beta groove parser
# ---------------------------------------------------------------------------


def decompose_class_ii_beta(
    seq: str,
    *,
    allele: str = "",
    gene: str = "",
    features: Optional[SequenceFeatures] = None,
    use_early_shortcuts: bool = True,
) -> AlleleRecord:
    """Parse a class-II beta chain via the shared grammar-driven engine.

    Enumerates a small top-K set of primary full parses, then lets them
    compete with β1-only salvage, fragment, and missing-groove
    explanations.
    """
    grammar = _grammar_spec("II", "beta")
    flags: list[str] = []
    early_result, cleaned, _pairs, annotated, sp_est = _prepare_parse_inputs(
        seq,
        grammar=grammar,
        allele=allele,
        gene=gene,
        features=features,
        use_early_shortcuts=use_early_shortcuts,
    )
    if early_result is not None:
        return early_result

    candidates: list[AlleleRecord] = []
    beta_resolver, beta1_ann = _build_beta_secondary_resolver(cleaned, annotated, features)
    selections = _collect_primary_parse_candidates(
        cleaned,
        grammar=grammar,
        annotated=annotated,
        sp_estimate=sp_est,
        features=features,
        secondary_resolver=beta_resolver,
    )
    best_record = None
    expand_primary = True
    if selections:
        best_record = _build_primary_result(
            allele=allele,
            gene=gene,
            seq=cleaned,
            grammar=grammar,
            selection=selections[0],
            sp_estimate=sp_est,
            features=features,
            flags=flags,
            secondary_pair=((beta1_ann.c1, beta1_ann.c2) if beta1_ann is not None else None),
        )
        _append_record_candidate(candidates, best_record)
        expand_primary = _should_expand_primary_competition(
            selections,
            best_record,
            grammar=grammar,
            features=features,
            use_early_shortcuts=use_early_shortcuts,
        )
        if expand_primary:
            for selection in selections[1:]:
                record = _build_primary_result(
                    allele=allele,
                    gene=gene,
                    seq=cleaned,
                    grammar=grammar,
                    selection=selection,
                    sp_estimate=sp_est,
                    features=features,
                    flags=flags,
                    secondary_pair=((beta1_ann.c1, beta1_ann.c2) if beta1_ann is not None else None),
                )
                _append_record_candidate(candidates, record)

    if beta1_ann is not None and (
        best_record is None or best_record.status != "ok" or best_record.parse_candidate.support_state != "support_present"
    ):
        _append_record_candidate(
            candidates,
            _build_class_ii_beta1_salvage_result(
                seq=cleaned,
                allele=allele,
                gene=gene,
                beta1_ann=beta1_ann,
                features=features,
            ),
        )

    if len(cleaned) <= grammar.fragment_max_len:
        _append_record_candidate(
            candidates,
            _class_ii_fragment_result(seq=cleaned, allele=allele, gene=gene, chain="beta"),
        )

    missing_flags = ("missing_support_anchor", "beta1_salvage_failed") if beta1_ann is not None else grammar.missing_no_anchor_flags
    _append_record_candidate(
        candidates,
        _missing_groove_result(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="beta",
            sequence=cleaned,
            seq_len=len(cleaned),
            anchor_type=("beta1_cys" if beta1_ann is not None else ""),
            anchor_cys1=(beta1_ann.c1 if beta1_ann is not None else None),
            anchor_cys2=(beta1_ann.c2 if beta1_ann is not None else None),
            flags=missing_flags,
        ),
    )
    return _select_best_record_candidate(candidates)


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


@lru_cache(maxsize=4096)
def _decompose_domains_cached(
    cleaned_seq: str,
    mhc_class: str,
    chain: str,
    allele: str,
    gene: str,
    use_early_shortcuts: bool,
) -> AlleleRecord:
    """Cache exact-sequence parses for repeated proteins with identical metadata."""
    features = analyze_sequence(cleaned_seq)
    return decompose_domains(
        cleaned_seq,
        mhc_class=mhc_class,
        chain=chain or None,
        allele=allele,
        gene=gene,
        features=features,
        use_early_shortcuts=use_early_shortcuts,
    )


def decompose_domains(
    seq: str,
    *,
    mhc_class: str,
    chain: Optional[str] = None,
    allele: str = "",
    gene: str = "",
    mutations: Sequence[object] = (),
    features: Optional[SequenceFeatures] = None,
    use_early_shortcuts: bool = True,
) -> AlleleRecord:
    """Dispatch groove parsing by class and chain.

    If ``mutations`` is provided, the groove is first extracted from the
    wild-type sequence (preserving Cys-pair detection), then the mutations
    are applied to the result.  Mutations use mature protein numbering
    (1-indexed) and can be strings ("K66A"), tuples, or mhcgnomes Mutation
    objects.
    """
    nc = normalize_mhc_class(mhc_class)
    cleaned_seq = features.seq if features is not None else _clean_seq(seq)
    if features is None and not mutations:
        return _decompose_domains_cached(
            cleaned_seq,
            nc,
            str(chain or ""),
            allele,
            gene,
            use_early_shortcuts,
        )

    gene_token = str(gene or "").strip().upper()
    if not gene_token and allele:
        try:
            gene_token = str(infer_gene(allele) or "").strip().upper()
        except Exception:
            gene_token = ""

    # Filter non-MHC genes that appear in MHC-region datasets.
    # LOC* gene names are NOT filtered here — they can be real MHC with
    # NCBI identifiers.  The pipeline filters LOC* during build instead.
    if gene_token in {g.upper() for g in NON_MHC_GENE_NAMES}:
        return _attach_parse_candidate(
            AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class=nc or mhc_class,
                chain=str(chain or ""),
                sequence=cleaned_seq,
                seq_len=len(cleaned_seq),
                status="non_groove",
                anchor_type="non_mhc_gene",
                flags=("non_mhc_gene",),
            )
        )

    # Infer class from gene name patterns when mhcgnomes can't resolve.
    if nc not in ("I", "II") and gene_token:
        gene_raw = str(gene or "").strip()
        if any(gene_raw.startswith(p) or gene_raw == p for p in GENE_CLASS_I_PATTERNS):
            nc = "I"
        elif any(gene_raw.startswith(p) or gene_raw == p for p in GENE_CLASS_II_ALPHA_PATTERNS):
            nc = "II"
            chain = "alpha"
        elif any(gene_raw.startswith(p) or gene_raw == p for p in GENE_CLASS_II_BETA_PATTERNS):
            nc = "II"
            chain = "beta"

    if gene_token in {"B2M", "BETA-2-MICROGLOBULIN"}:
        chain_token = str(chain or "").strip().lower()
        inferred_chain = ""
        if nc == "I":
            inferred_chain = "alpha"
        elif chain_token in {"a", "alpha", "mhc_a"}:
            inferred_chain = "alpha"
        elif chain_token in {"b", "beta", "mhc_b"}:
            inferred_chain = "beta"
        else:
            inferred_chain = _class_ii_chain_from_name(gene=gene, allele=allele) or ""
        return _attach_parse_candidate(
            AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class=nc,
                chain=inferred_chain,
                sequence=cleaned_seq,
                seq_len=len(cleaned_seq),
                status="non_groove",
                anchor_type="excluded_gene",
                flags=("non_groove_gene",),
            )
        )

    if features is None:
        features = analyze_sequence(cleaned_seq)

    if nc == "I":
        result = decompose_class_i(
            cleaned_seq,
            allele=allele,
            gene=gene,
            features=features,
            use_early_shortcuts=use_early_shortcuts,
        )
        if mutations and result.ok:
            result = apply_mutations(result, mutations)
        result = _refine_status(result)
        if gene_token in NON_GROOVE_GENES and "non_groove_gene" not in result.flags:
            result = replace(result, flags=_flags_to_tuple((*result.flags, "non_groove_gene")))
        return _attach_parse_candidate(result)
    if nc not in ("I", "II"):
        # Unknown class — try all three parsers and pick the best.
        # This handles sequences with missing or unrecognized metadata
        # (common in reptiles, amphibians, and non-model organisms).
        candidates: list[AlleleRecord] = []
        for parser, name in [
            (decompose_class_i, "class_I"),
            (decompose_class_ii_beta, "class_II_beta"),
            (decompose_class_ii_alpha, "class_II_alpha"),
        ]:
            try:
                r = parser(
                    cleaned_seq,
                    allele=allele,
                    gene=gene,
                    features=features,
                    use_early_shortcuts=use_early_shortcuts,
                )
                if r.ok and r.mature_start >= 0:
                    candidates.append(r)
            except Exception:
                pass
        result = _select_best_viable_record(candidates)
        if result is not None:
            if mutations and result.ok:
                result = apply_mutations(result, mutations)
            result = _refine_status(result)
            return _attach_parse_candidate(result)
        # All parsers failed — return a structured failure
        return _attach_parse_candidate(
            AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class=mhc_class or "",
                chain=str(chain or ""),
                sequence=cleaned_seq,
                seq_len=len(cleaned_seq),
                status="missing_groove",
                flags=("unknown_class_all_parsers_failed",),
            )
        )

    chain_token = str(chain or "").strip().lower()
    result: Optional[AlleleRecord] = None
    if chain_token in {"a", "alpha", "mhc_a"}:
        result = decompose_class_ii_alpha(
            cleaned_seq,
            allele=allele,
            gene=gene,
            features=features,
            use_early_shortcuts=use_early_shortcuts,
        )
    elif chain_token in {"b", "beta", "mhc_b"}:
        result = decompose_class_ii_beta(
            cleaned_seq,
            allele=allele,
            gene=gene,
            features=features,
            use_early_shortcuts=use_early_shortcuts,
        )
    elif chain_token == "unknown":
        # Unknown chain for class II — try both and pick best
        alpha_r = decompose_class_ii_alpha(
            cleaned_seq,
            allele=allele,
            gene=gene,
            features=features,
            use_early_shortcuts=use_early_shortcuts,
        )
        beta_r = decompose_class_ii_beta(
            cleaned_seq,
            allele=allele,
            gene=gene,
            features=features,
            use_early_shortcuts=use_early_shortcuts,
        )
        result = _select_best_viable_record([alpha_r, beta_r])
    elif chain_token:
        # Unrecognized chain token — try both as fallback
        alpha_r = decompose_class_ii_alpha(
            cleaned_seq,
            allele=allele,
            gene=gene,
            features=features,
            use_early_shortcuts=use_early_shortcuts,
        )
        beta_r = decompose_class_ii_beta(
            cleaned_seq,
            allele=allele,
            gene=gene,
            features=features,
            use_early_shortcuts=use_early_shortcuts,
        )
        result = _select_best_viable_record([alpha_r, beta_r])
    else:
        # Infer chain from gene name
        name_chain = _class_ii_chain_from_name(gene=gene, allele=allele)
        if name_chain == "alpha":
            result = decompose_class_ii_alpha(
                cleaned_seq,
                allele=allele,
                gene=gene,
                features=features,
                use_early_shortcuts=use_early_shortcuts,
            )
        elif name_chain == "beta":
            result = decompose_class_ii_beta(
                cleaned_seq,
                allele=allele,
                gene=gene,
                features=features,
                use_early_shortcuts=use_early_shortcuts,
            )
        else:
            # Last resort: try both parsers
            alpha_r = decompose_class_ii_alpha(
                cleaned_seq,
                allele=allele,
                gene=gene,
                features=features,
                use_early_shortcuts=use_early_shortcuts,
            )
            beta_r = decompose_class_ii_beta(
                cleaned_seq,
                allele=allele,
                gene=gene,
                features=features,
                use_early_shortcuts=use_early_shortcuts,
            )
            result = _select_best_viable_record([alpha_r, beta_r])
            if result is None:
                result = AlleleRecord(
                    allele=allele,
                    gene=gene,
                    mhc_class="II",
                    chain="",
                    sequence=features.seq,
                    seq_len=features.seq_len,
                    status="chain_inference_failed",
                    anchor_type="chain_inference",
                    flags=_flags_to_tuple(
                        (
                            f"alpha_status={alpha_r.status}",
                            f"beta_status={beta_r.status}",
                        )
                    ),
                )

    if result is None:
        result = AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain=str(chain or ""),
            sequence=features.seq,
            seq_len=features.seq_len,
            status="missing_groove",
            flags=("class_ii_all_parsers_failed",),
        )
    if mutations and result.ok:
        result = apply_mutations(result, mutations)
    result = _refine_status(result)
    if gene_token in NON_GROOVE_GENES and "non_groove_gene" not in result.flags:
        result = replace(result, flags=_flags_to_tuple((*result.flags, "non_groove_gene")))
    return _attach_parse_candidate(result)


def is_class_ii_alpha_gene(gene: str) -> bool:
    """Whether a gene name corresponds to a class II alpha chain.

    NOTE: This function uses a trailing-letter heuristic as a fallback
    (gene ending in "A" → alpha). This can misclassify fish class I
    genes (UBA, UCA, etc.) as class II alpha. Callers should gate on
    mhc_class == "II" before using this function.
    """
    token = str(gene or "").strip().upper()
    return token.startswith(CLASS_II_ALPHA_GENE_PREFIXES) or token.endswith("A")


# Backward compatibility aliases
extract_groove = decompose_domains
parse_class_i = decompose_class_i
parse_class_ii_alpha = decompose_class_ii_alpha
parse_class_ii_beta = decompose_class_ii_beta
