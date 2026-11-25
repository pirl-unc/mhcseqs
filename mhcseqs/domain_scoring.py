"""Candidate generation, scaffold scoring, primary selection, and result assembly."""

from __future__ import annotations

import re
from bisect import bisect_right
from dataclasses import dataclass, replace
from typing import Optional, Sequence

from .alleles import infer_gene
from .domain_features import (
    _cached_between_cys,
    _cached_topology_support,
    _find_tm_span,
    _score_tail,
    classify_cys_pair,
    find_cys_pairs,
    find_junction,
    score_groove_ig_boundary,
    score_junction,
)
from .domain_grammar import (
    CLASS_I_ALPHA2_END_AFTER_CYS2,
    CLASS_II_ALPHA_GENE_PREFIXES,
    CLASS_II_BETA_GENE_PREFIXES,
    IG_DOMAIN_END_AFTER_CYS2,
    MAX_PLAUSIBLE_SP,
    MIN_FUNCTIONAL_GROOVE_HALF_LEN,
    MIN_GROOVE_SOURCE_LEN,
    NON_CLASSICAL_CLASS_I_GENE_PATTERNS,
    PRIMARY_PARSE_CANDIDATE_KEEP,
    PRIMARY_PARSE_EXPANSION_MARGIN,
    PRIMARY_PARSE_LOW_CONFIDENCE,
    GrammarSpec,
)
from .domain_grammar import (
    GROOVE_PRIOR as _GROOVE_PRIOR,
)
from .domain_grammar import (
    GROOVE_SEGMENT_PRIOR as _GROOVE_SEGMENT_PRIOR,
)
from .domain_grammar import (
    IG_PRIOR as _IG_PRIOR,
)
from .domain_grammar import (
    PARSE_CANDIDATE_COMPAT_WEIGHTS as _PARSE_CANDIDATE_COMPAT_WEIGHTS,
)
from .domain_grammar import (
    PRIMARY_CANDIDATE_CONTRADICTION_WEIGHTS as _PRIMARY_CANDIDATE_CONTRADICTION_WEIGHTS,
)
from .domain_grammar import (
    RELATIVE_OFFSET_PRIOR as _RELATIVE_OFFSET_PRIOR,
)
from .domain_grammar import (
    SP_CHARGED as _SP_CHARGED,
)
from .domain_grammar import (
    SP_EARLIER_TIE_MARGIN as _SP_EARLIER_TIE_MARGIN,
)
from .domain_types import (
    AlleleRecord,
    CysPairAnnotation,
    ExtractedDomains,
    ParseCandidate,
    ParseSubscores,
    PrimaryParseSelection,
    SequenceFeatures,
    _clean_seq,
    _flags_to_tuple,
    _grammar_for_record,
    _hydrophobic_fraction,
    _is_truncated_class_i_sp,
    _sigmoid,
    _slice_or_empty,
    _soft_prior_score,
    grammar_spec,
)
from .signal_peptide import (
    _append_long_sp_flag,
    _append_sp_overcapture_flag,
    _score_mature_nterm,
    _score_nterm_lexical_state,
    _score_sp_cleavage,
    _score_sp_estimate_consistency,
    _score_sp_estimate_zone_consistency,
    _score_sp_length,
    _selected_sp_estimate,
    detect_h_region,
    fast_sp_triage,
    infer_signal_peptide,
)


def _score_class_i_groove_lengths(alpha1_len: int, alpha2_len: int) -> tuple[float, float, float]:
    """Return class-I groove1 / groove2 / total length evidence."""
    g1 = _soft_prior_score(alpha1_len, *_GROOVE_SEGMENT_PRIOR[("I", "groove1")])
    g2 = _soft_prior_score(alpha2_len, *_GROOVE_SEGMENT_PRIOR[("I", "groove2")])
    total = _soft_prior_score(alpha1_len + alpha2_len, *_GROOVE_SEGMENT_PRIOR[("I", "total")]) * 0.5
    return g1, g2, total


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
            cleavage_score = _score_sp_cleavage(
                self.seq,
                mature_start,
                h_region=self.h_region,
                features=self.features,
            )
            sp_len_score = _score_sp_length(mature_start)
            mature_nterm_score = _score_mature_nterm(self.seq, mature_start, self.grammar.mhc_class)
            truncated_class_i_sp = self.grammar.mhc_class == "I" and _is_truncated_class_i_sp(
                self.seq,
                mature_start,
                self.groove_boundary,
                cleavage_score=cleavage_score,
                mature_nterm_score=mature_nterm_score,
            )
            sp = (
                cleavage_score
                + (max(0.0, sp_len_score) if truncated_class_i_sp else sp_len_score)
                + mature_nterm_score
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

            if truncated_class_i_sp:
                pass
            elif not starts_met and h_frac < 0.55:
                sp -= 5.0  # no Met + no real h-region = almost certainly no SP
            elif not starts_met and h_frac < 0.65:
                sp -= 3.0  # no Met + weak h-region = very suspicious
            elif not starts_met:
                sp -= 1.0  # no Met but good h-region = suspicious, not impossible
            elif h_frac < 0.55:
                sp -= 2.0  # Met but no real h-region = suspicious
            elif h_frac < 0.65:
                sp -= 0.5  # Met but borderline h-region

            if not cleavage_after_h and h_end > h_start and not truncated_class_i_sp:
                sp -= 1.5  # cleavage before h-region = structurally wrong

        groove = self.fixed_groove
        if self.grammar.length_mode == "class_i_split":
            alpha1_len = self.groove_boundary - mature_start
            if alpha1_len > 0 and self.class_i_alpha2_len > 0:
                g1, g2, total = _score_class_i_groove_lengths(alpha1_len, self.class_i_alpha2_len)
                # SP-stripped mature-only sequences (pos=0, no Met) have
                # truncated α1 because the SP and N-terminal residues are
                # missing from the database deposit.  Cap the α1 length
                # penalty so good Cys pair topology can still produce a
                # viable parse.  This recovers horse Eqca, pig SLA, and
                # similar IPD-MHC deposits.
                if mature_start == 0 and self.seq[:1] != "M" and g1 < -3.0:
                    g1 = -3.0
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
    grammar = grammar_spec(mhc_class, chain)
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
        fixed_groove += score_junction(seq, groove_boundary) * weights.groove_junction
    if weights.groove_boundary and 5 < groove_end < n - 5:
        fixed_groove += score_groove_ig_boundary(seq, groove_end, grammar.mhc_class, grammar.chain) * weights.groove_boundary

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
    score += score_junction(seq, alpha2_start) * 0.8 if alpha2_start + 5 < len(seq) else -8.0
    score += score_groove_ig_boundary(seq, alpha3_start, "I", "alpha") if 5 < alpha3_start < len(seq) - 5 else -2.0
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
    score += score_groove_ig_boundary(seq, groove_end, "II", "beta") if 5 < groove_end < len(seq) - 5 else -1.0
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

    best_pos = -1
    best_score = -999.0

    if search_hi < search_lo:
        # No valid SP range.  For sequences that look SP-stripped (no Met,
        # no h-region), fall through to the pos=0 check below — these are
        # mature-only deposits where α1 is truncated (e.g., horse Eqca,
        # pig SLA in IPD-MHC).  For Met-starting sequences, return -1 so
        # the salvage path can try alternative Cys pairs.
        if scaffold.seq[:1] == "M":
            return -1

    for pos in range(search_lo, search_hi + 1):
        s = scaffold.score(pos)
        if sp_candidates:
            zone_score = _score_sp_estimate_zone_consistency(pos, sp_candidates, sp_confidence)
            # Scale zone weight by h-region quality: when the h-region is
            # very strong (frac ≥ 0.75), trust the zone more heavily.
            h_s, h_e = scaffold.h_region
            if h_e > h_s and scaffold.seq[:1] == "M":
                h_frac = _hydrophobic_fraction(scaffold.seq, h_s, h_e, features=scaffold.features)
                if h_frac >= 0.75:
                    zone_score *= 1.5
            s += zone_score
        elif sp_estimate > 0:
            s += _score_sp_estimate_consistency(pos, sp_estimate, sp_confidence)
        if s > best_score + _SP_EARLIER_TIE_MARGIN:
            best_score = s
            best_pos = pos
        elif best_pos < 0 or (pos < best_pos and s >= best_score - _SP_EARLIER_TIE_MARGIN):
            best_score = s
            best_pos = pos

    # Also check pos=0 (SP-stripped / leaderless).
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
    # For SP-stripped sequences (mature_start=0, no Met), groove1 may be
    # truncated because the N-terminal residues were missing from the
    # database deposit.  Use a lower threshold (40 aa) for these.
    is_sp_stripped = result.mature_start == 0 and not (result.sequence or "").startswith("M")
    min_len = 40 if is_sp_stripped else MIN_FUNCTIONAL_GROOVE_HALF_LEN
    for field_name, _role in grammar.groove_segments:
        seg_len = int(getattr(result, f"{field_name}_len", 0) or 0)
        if seg_len > 0 and seg_len < min_len:
            flags.append(f"{field_name}_short({seg_len})")
            return replace(result, status="short", flags=_flags_to_tuple(flags))

    return result


# ---------------------------------------------------------------------------
# Class II chain inference
# ---------------------------------------------------------------------------


def _normalized_gene_token(gene: str, allele: str = "") -> str:
    """Return an unprefixed uppercase gene token for dispatch checks."""
    raw = str(gene or allele or "").strip()
    if gene and not raw.startswith("~") and "-" in raw:
        body = raw.split("-", 1)[1].split("*", 1)[0]
        if body.upper() in {"H60A", "H60B", "H60C"}:
            return body.upper()
    inferred = infer_gene(raw)
    return str(inferred or raw).strip().upper()


def _class_ii_chain_from_name(
    *,
    gene: str,
    allele: str,
) -> Optional[str]:
    gene_token = _normalized_gene_token(gene, allele)

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


def _build_class_i_junction_salvage_result(
    *,
    seq: str,
    allele: str,
    gene: str,
    features: Optional[SequenceFeatures] = None,
) -> Optional[AlleleRecord]:
    """Recover a Class I α1/α2 split from a conserved junction when an α2 Cys is lost.

    Handles alleles such as ``HLA-A*30:14L`` (C164S), ``Mafa-B*202:01`` (C101S),
    ``BoLA-NC9*001:01``, and ``Chae-G*02:01`` that retain the canonical
    ``[YF]NQS..|GSHT..`` α1/α2 boundary plus at least one surviving α2 Cys
    landmark (Cys101 near ``junction+10`` or Cys164 near ``junction+73``).
    """
    cleaned = features.seq if features else _clean_seq(seq)
    n = len(cleaned)
    if n < 60:
        return None

    best_jpos = -1
    best_jscore = -999.0
    for pos in range(25, n - 18):
        upstream = cleaned[pos - 6 : pos]
        if upstream[0] not in "YF" or upstream[1:4] != "NQS":
            continue
        jscore = score_junction(cleaned, pos)
        if jscore < 8.0:
            continue
        has_cys101 = any(cleaned[p] == "C" and p > 0 and cleaned[p - 1] in "GSA" for p in range(pos + 8, min(n, pos + 14)))
        has_cys164 = any(cleaned[p] == "C" and any(cleaned[p + d] in "WL" for d in (3, 4) if p + d < n) for p in range(pos + 69, min(n, pos + 78)))
        if not (has_cys101 or has_cys164):
            continue
        if jscore > best_jscore:
            best_jscore = jscore
            best_jpos = pos

    if best_jpos < 0:
        return None

    mature_start = 0
    if best_jpos > 94:
        h_region = features.h_region if features is not None else detect_h_region(cleaned)
        sp_pos, _sp_conf = infer_signal_peptide(cleaned, h_region=h_region)
        if sp_pos > 0 and 82 <= (best_jpos - sp_pos) <= 96:
            mature_start = sp_pos
        else:
            # Check for a truncated Class I SP tail upstream of canonical 90-aa α1.
            best_trunc_pos = 0
            best_trunc_score = -999.0
            for cand in range(max(5, best_jpos - 94), min(17, best_jpos - 85)):
                cl = _score_sp_cleavage(cleaned, cand, h_region=h_region, features=features)
                mn = _score_mature_nterm(cleaned, cand, "I")
                if _is_truncated_class_i_sp(cleaned, cand, best_jpos, cleavage_score=cl, mature_nterm_score=mn):
                    total_s = cl + mn
                    if total_s > best_trunc_score:
                        best_trunc_score = total_s
                        best_trunc_pos = cand
            mature_start = best_trunc_pos

    if n <= best_jpos + 96:
        alpha2_end = n
    else:
        alpha2_end = min(n, best_jpos + 93)
        best_bscore = -999.0
        for bpos in range(best_jpos + 88, min(n, best_jpos + 98)):
            bscore = score_groove_ig_boundary(cleaned, bpos, "I", "alpha")
            if bscore > best_bscore:
                best_bscore = bscore
                alpha2_end = bpos

    half_1 = _slice_or_empty(cleaned, mature_start, best_jpos)
    half_2 = _slice_or_empty(cleaned, best_jpos, alpha2_end)
    # Reject runaway junction-inferred alpha1 segments (e.g. chimeric / N-extended proteins)
    if len(half_1) > 120:
        return None
    ig_end = min(n, alpha2_end + 95) if alpha2_end < n else n
    ig = _slice_or_empty(cleaned, alpha2_end, ig_end)
    tail = cleaned[ig_end:]

    flags = ["inferred_from_junction"]
    if not ig:
        flags.append("missing_support_anchor")
    if len(half_1) < 80:
        flags.append(f"alpha1_short({len(half_1)})")
    if len(half_2) < 80:
        flags.append(f"alpha2_short({len(half_2)})")

    # When BOTH groove halves are severely truncated (e.g. 105-aa partial exon2+3
    # fragments with 59+46 aa), downgrade from "ok" to "short".
    status = "short" if (len(half_1) < MIN_FUNCTIONAL_GROOVE_HALF_LEN and len(half_2) < MIN_FUNCTIONAL_GROOVE_HALF_LEN) else "ok"

    return _attach_parse_candidate(
        AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="I",
            chain="alpha",
            sequence=cleaned,
            seq_len=n,
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
            parse_score=best_jscore + 10.0,
            status=status,
            anchor_type="junction",
            flags=_flags_to_tuple(flags),
        ),
        seq=cleaned,
        h_region=(features.h_region if features is not None else None),
        sp_estimate=(features.sp_estimate if features is not None else None),
        features=features,
    )


def _class_i_fragment_result(
    *,
    seq: str,
    allele: str,
    gene: str,
    annotated: Optional[Sequence[CysPairAnnotation]] = None,
    features: Optional[SequenceFeatures] = None,
) -> AlleleRecord:
    """Return an alpha1_only, alpha2_only, partial α1+α2, or missing_groove record.

    Classifies candidate Cys pairs with ``classify_cys_pair`` rather than
    relying on raw separation alone:
      - Strong Ig-only Cys pair (α3 + TM/tail fragment) → ``missing_groove``
        with ``ig_only_fragment`` (never mislabeled ``alpha2_only``).
      - Strong α2 groove Cys pair with an upstream α1/α2 junction → split into
        partial ``groove1``, ``groove2``, and any trailing ``ig_domain``.
      - Strong α2 groove Cys pair starting at α2 → ``alpha2_only``.
      - No α2 groove Cys pair → try junction salvage before ``alpha1_only``.
    """
    cleaned = _clean_seq(seq)
    if annotated is None:
        pairs = find_cys_pairs(cleaned)
        annotated = [classify_cys_pair(cleaned, c1, c2) for c1, c2, _sep in pairs]

    groove_pairs = [ann for ann in annotated if 55 <= ann.separation <= 72 and ann.domain_type == "groove" and ann.groove_score >= 4.0]
    ig_pairs = [ann for ann in annotated if ann.domain_type == "ig" and ann.ig_score >= 4.0 and ann.ig_score > ann.groove_score + 2.0]

    if groove_pairs:
        best_groove = max(groove_pairs, key=lambda a: a.groove_score)
        jpos, jscore = find_junction(cleaned, best_groove.c1)
        if jpos >= 5 and jscore >= 6.0:
            alpha2_end = min(len(cleaned), best_groove.c2 + 20)
            half_1 = _slice_or_empty(cleaned, 0, jpos)
            half_2 = _slice_or_empty(cleaned, jpos, alpha2_end)
            ig = cleaned[alpha2_end:]
            flags = ["truncated_alpha1"]
            if not ig:
                flags.append("missing_support_anchor")
            return _attach_parse_candidate(
                AlleleRecord(
                    allele=allele,
                    gene=gene,
                    mhc_class="I",
                    chain="alpha",
                    sequence=cleaned,
                    seq_len=len(cleaned),
                    mature_start=0,
                    groove_seq=half_1 + half_2,
                    groove1=half_1,
                    groove2=half_2,
                    groove1_len=len(half_1),
                    groove2_len=len(half_2),
                    ig_domain=ig,
                    ig_domain_len=len(ig),
                    parse_score=best_groove.groove_score + jscore,
                    status="ok",
                    anchor_type="alpha2_cys",
                    anchor_cys1=best_groove.c1,
                    anchor_cys2=best_groove.c2,
                    flags=_flags_to_tuple(flags),
                ),
                seq=cleaned,
                h_region=(features.h_region if features is not None else None),
                sp_estimate=(features.sp_estimate if features is not None else None),
                features=features,
            )
        downstream_ig = [ig for ig in ig_pairs if ig.c1 >= best_groove.c2 + 15]
        if downstream_ig:
            best_ig = min(downstream_ig, key=lambda a: a.c1)
            alpha2_end = min(len(cleaned), best_groove.c2 + CLASS_I_ALPHA2_END_AFTER_CYS2)
            ig_end = min(len(cleaned), best_ig.c2 + IG_DOMAIN_END_AFTER_CYS2)
            half_2 = _slice_or_empty(cleaned, 0, alpha2_end)
            ig = _slice_or_empty(cleaned, alpha2_end, ig_end)
            tail = cleaned[ig_end:]
            return _attach_parse_candidate(
                AlleleRecord(
                    allele=allele,
                    gene=gene,
                    mhc_class="I",
                    chain="alpha",
                    sequence=cleaned,
                    seq_len=len(cleaned),
                    mature_start=0,
                    groove_seq=half_2,
                    groove1="",
                    groove2=half_2,
                    groove1_len=0,
                    groove2_len=len(half_2),
                    ig_domain=ig,
                    ig_domain_len=len(ig),
                    tail=tail,
                    tail_len=len(tail),
                    parse_score=best_groove.groove_score + best_ig.ig_score,
                    status="alpha2_only",
                    anchor_type="alpha2_cys",
                    anchor_cys1=best_groove.c1,
                    anchor_cys2=best_groove.c2,
                    secondary_cys1=best_ig.c1,
                    secondary_cys2=best_ig.c2,
                    flags=("missing_alpha1",),
                ),
                seq=cleaned,
                h_region=(features.h_region if features is not None else None),
                sp_estimate=(features.sp_estimate if features is not None else None),
                features=features,
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
                groove1="",
                groove2=cleaned,
                groove1_len=0,
                groove2_len=len(cleaned),
                status="alpha2_only",
                anchor_type="raw_fragment",
                anchor_cys1=best_groove.c1,
                anchor_cys2=best_groove.c2,
                flags=("single_exon_fragment",),
            )
        )

    if ig_pairs:
        best_ig = max(ig_pairs, key=lambda a: a.ig_score)
        ig_end = min(len(cleaned), best_ig.c2 + IG_DOMAIN_END_AFTER_CYS2)
        ig = cleaned[:ig_end]
        tail = cleaned[ig_end:]
        return _attach_parse_candidate(
            AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class="I",
                chain="alpha",
                sequence=cleaned,
                seq_len=len(cleaned),
                mature_start=0,
                ig_domain=ig,
                ig_domain_len=len(ig),
                tail=tail,
                tail_len=len(tail),
                status="missing_groove",
                anchor_type="alpha3_cys",
                anchor_cys1=best_ig.c1,
                anchor_cys2=best_ig.c2,
                flags=("ig_only_fragment",),
            )
        )

    junction_salvage = _build_class_i_junction_salvage_result(seq=cleaned, allele=allele, gene=gene, features=features)
    if junction_salvage is not None:
        return junction_salvage

    # Guard against C-terminal α3 Ig half-fragments (missing α3 Cys203, so no
    # Cys pair exists, e.g. MELVETRP...GTFQKW...CHVHHEGLP) falling through to alpha1_only.
    if re.search(r"G[STN]FQKW[AGSV][ASV]|CHV[QHNYL]H[EDQ]GLP|LVETRP[AST][GSN]D?GTFQ", cleaned) and not re.search(
        r"G[SP]HS[LMF]RYF|VGYVDDT|RFDSDA", cleaned
    ):
        return _attach_parse_candidate(
            AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class="I",
                chain="alpha",
                sequence=cleaned,
                seq_len=len(cleaned),
                mature_start=0,
                ig_domain=cleaned,
                ig_domain_len=len(cleaned),
                status="missing_groove",
                anchor_type="",
                flags=("ig_only_fragment",),
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


def prepare_parse_inputs(
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
        if grammar.mhc_class == "I":
            junction_salvage = _build_class_i_junction_salvage_result(seq=cleaned, allele=allele, gene=gene, features=features)
            if junction_salvage is not None:
                return junction_salvage, cleaned, [], [], 0
        if len(cleaned) <= grammar.fragment_max_len:
            if grammar.mhc_class == "I":
                return (
                    _class_i_fragment_result(seq=cleaned, allele=allele, gene=gene, annotated=(), features=features),
                    cleaned,
                    [],
                    [],
                    0,
                )
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


def build_class_i_secondary_resolver(
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


def build_beta_secondary_resolver(
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
        if not grammar_spec("II", "beta").allow_self_secondary and (anchor.c1, anchor.c2) == (beta1_ann.c1, beta1_ann.c2):
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
    domains = extract_primary_domains(
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


def collect_primary_parse_candidates(
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

    strong_ig_pairs = (
        [a for a in annotated if a.domain_type == "ig" and a.ig_score >= 4.0 and a.ig_score > a.groove_score + 2.0]
        if grammar.mhc_class == "I"
        else []
    )
    for ann in annotated:
        # Class I primary anchor must be the G-ALPHA2 groove disulfide, never a
        # strongly C1-set Ig domain (nor an overlapping Cys pair spanning one).
        if grammar.mhc_class == "I" and any(ann.c1 <= ig.c1 and ig.c2 <= ann.c2 for ig in strong_ig_pairs):
            continue
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


def select_best_primary_parse(
    seq: str,
    *,
    grammar: GrammarSpec,
    annotated: Sequence[CysPairAnnotation],
    sp_estimate: int,
    features: Optional[SequenceFeatures] = None,
    secondary_resolver=None,
) -> Optional[PrimaryParseSelection]:
    """Return the best primary parse under one grammar."""
    selections = collect_primary_parse_candidates(
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


def extract_primary_domains(
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


def _domains_match_grammar(grammar: GrammarSpec, domains: ExtractedDomains) -> bool:
    """Return True when the materialized groove segments satisfy the grammar."""
    return all(bool(getattr(domains, field_name)) for field_name, _role in grammar.groove_segments)


def build_primary_result(
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
    _append_sp_overcapture_flag(local_flags, seq, mature_start)
    domains = extract_primary_domains(
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
            candidate_score=selection.selection_score,
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
    _apply_short_flags(grammar_spec("II", "beta"), domains, flags)
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
