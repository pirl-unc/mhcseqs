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

from dataclasses import replace
from functools import lru_cache
from typing import Optional, Sequence

from .alleles import is_non_mhc_gene, normalize_mhc_class, parse_gene_class
from .domain_features import (
    _apply_requested_mutations,
    _cached_topology_support,
    _normalize_mutations,
    apply_mutations,
    classify_cys_pair,
    classify_domain_fold,
    find_cys_pairs,
    find_junction,
    infer_structural_domains,
    score_cys_flanking_properties,
    score_groove_ig_boundary,
    score_junction,
)
from .domain_grammar import (
    CLASS_I_ALPHA2_END_AFTER_CYS2,
    CLASS_II_ALPHA_GENE_PREFIXES,
    MAX_PLAUSIBLE_SP,
    MIN_FUNCTIONAL_GROOVE_HALF_LEN,
    MIN_GROOVE_SOURCE_LEN,
    NON_GROOVE_GENES,
    PRIMARY_PARSE_CANDIDATE_KEEP,
    GrammarSpec,
)
from .domain_scoring import (
    ParseScaffold,
    _append_record_candidate,
    _attach_parse_candidate,
    _build_class_i_alpha3_salvage_result,
    _build_class_i_junction_salvage_result,
    _build_class_ii_beta1_salvage_result,
    _class_i_fragment_result,
    _class_ii_chain_from_name,
    _class_ii_fragment_result,
    _missing_groove_result,
    _normalized_gene_token,
    _refine_status,
    _select_best_record_candidate,
    _select_best_viable_record,
    _should_expand_primary_competition,
    build_beta_secondary_resolver,
    build_class_i_secondary_resolver,
    build_primary_result,
    collect_primary_parse_candidates,
    extract_primary_domains,
    prepare_parse_inputs,
    select_best_primary_parse,
)
from .domain_types import (
    AlleleRecord,
    BoundaryParseObjective,
    CysPairAnnotation,
    DomainFoldAnnotation,
    ExtractedDomains,
    ParseCandidate,
    ParseSubscores,
    ParseTrace,
    PrimaryParseSelection,
    RawAllele,
    SequenceFeatures,
    StructuralDomain,
    _clean_seq,
    _flags_to_tuple,
    grammar_spec,
)
from .signal_peptide import (
    analyze_sequence,
    detect_h_region,
    detect_sp_internal_met_overcapture,
    estimate_sp_from_h_region,
    fast_sp_triage,
    has_plausible_sp_core,
    infer_signal_peptide,
    infer_signal_peptide_candidates,
    load_sp_boundary_model,
    load_sp_sequence_cue_model,
    refine_signal_peptide,
    refinement_group,
    sp_boundary_excluded,
)

__all__ = [
    "AlleleRecord",
    "BoundaryParseObjective",
    "CLASS_I_ALPHA2_END_AFTER_CYS2",
    "CLASS_II_ALPHA_GENE_PREFIXES",
    "CysPairAnnotation",
    "DomainFoldAnnotation",
    "ExtractedDomains",
    "GrammarSpec",
    "MAX_PLAUSIBLE_SP",
    "MIN_FUNCTIONAL_GROOVE_HALF_LEN",
    "MIN_GROOVE_SOURCE_LEN",
    "NON_GROOVE_GENES",
    "PRIMARY_PARSE_CANDIDATE_KEEP",
    "ParseCandidate",
    "ParseScaffold",
    "ParseSubscores",
    "ParseTrace",
    "PrimaryParseSelection",
    "RawAllele",
    "SequenceFeatures",
    "StructuralDomain",
    "analyze_sequence",
    "apply_mutations",
    "build_beta_secondary_resolver",
    "build_class_i_secondary_resolver",
    "build_primary_result",
    "classify_cys_pair",
    "classify_domain_fold",
    "collect_primary_parse_candidates",
    "decompose_class_i",
    "decompose_class_ii_alpha",
    "decompose_class_ii_beta",
    "decompose_domains",
    "detect_h_region",
    "detect_sp_internal_met_overcapture",
    "estimate_sp_from_h_region",
    "extract_groove",
    "extract_primary_domains",
    "fast_sp_triage",
    "find_cys_pairs",
    "find_junction",
    "grammar_spec",
    "has_plausible_sp_core",
    "infer_signal_peptide",
    "infer_signal_peptide_candidates",
    "infer_structural_domains",
    "is_class_ii_alpha_gene",
    "load_sp_boundary_model",
    "load_sp_sequence_cue_model",
    "parse_class_i",
    "parse_class_ii_alpha",
    "parse_class_ii_beta",
    "prepare_parse_inputs",
    "refine_signal_peptide",
    "refinement_group",
    "score_cys_flanking_properties",
    "score_groove_ig_boundary",
    "score_junction",
    "select_best_primary_parse",
    "sp_boundary_excluded",
    "trace_parse_class_i",
]

# ---------------------------------------------------------------------------
# Class I groove parser
# ---------------------------------------------------------------------------


def decompose_class_i(
    seq: str,
    *,
    allele: str = "",
    gene: str = "",
    mutations: str | Sequence[object] = (),
    features: Optional[SequenceFeatures] = None,
    use_early_shortcuts: bool = True,
) -> AlleleRecord:
    """Parse a class-I alpha chain via the shared grammar-driven engine.

    Enumerates a small top-K set of structurally distinct primary parses over
    plausible groove Cys pairs and mature_start positions, then lets those
    full parses compete with fragment, salvage, and missing-groove
    explanations under one candidate selector.

    Pass *mutations* (e.g. ``"K66A"`` or ``["K66A", "D77S"]``) to apply
    1-indexed mature point mutations to the parsed wild-type record.

    Pass *features* from ``analyze_sequence()`` to skip redundant
    Cys-pair scanning and fold classification.
    """
    grammar = grammar_spec("I", "alpha")
    flags: list[str] = []
    early_result, cleaned, _pairs, annotated, sp_est = prepare_parse_inputs(
        seq,
        grammar=grammar,
        allele=allele,
        gene=gene,
        features=features,
        use_early_shortcuts=use_early_shortcuts,
    )
    if early_result is not None:
        return _apply_requested_mutations(early_result, mutations)

    candidates: list[AlleleRecord] = []
    selections = collect_primary_parse_candidates(
        cleaned,
        grammar=grammar,
        annotated=annotated,
        sp_estimate=sp_est,
        features=features,
        secondary_resolver=build_class_i_secondary_resolver(cleaned, annotated, features),
    )
    best_record: Optional[AlleleRecord] = None
    expand_primary = True
    if selections:
        best_record = build_primary_result(
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
                record = build_primary_result(
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

    nterm_alpha2_grooves = [
        ann for ann in annotated if 55 <= ann.separation <= 72 and ann.domain_type == "groove" and ann.groove_score >= 4.0 and ann.c1 <= 30
    ]
    strong_igs = [ann for ann in annotated if ann.domain_type == "ig" and ann.ig_score >= 4.0 and ann.ig_score > ann.groove_score + 2.0]
    has_nterm_alpha2_plus_ig = any(ig.c1 >= gr.c2 + 15 for gr in nterm_alpha2_grooves for ig in strong_igs)
    if len(cleaned) <= grammar.fragment_max_len or ((best_record is None or best_record.status != "ok") and has_nterm_alpha2_plus_ig):
        _append_record_candidate(
            candidates,
            _class_i_fragment_result(seq=cleaned, allele=allele, gene=gene, annotated=annotated, features=features),
        )

    if best_record is None or best_record.status != "ok":
        _append_record_candidate(
            candidates,
            _build_class_i_junction_salvage_result(seq=cleaned, allele=allele, gene=gene, features=features),
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
    return _apply_requested_mutations(_select_best_record_candidate(candidates), mutations)


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
        jpos, jscore = find_junction(cleaned, best_groove.c1)
        trace.junction_pos = jpos
        trace.junction_score = jscore
        trace.log.append(f"Junction search near C{best_groove.c1}-10={best_groove.c1 - 10}: best at {jpos} score={jscore:.1f}")
        # Show top 3 junction candidates
        j_candidates = []
        for p in range(max(0, best_groove.c1 - 15), min(len(cleaned), best_groove.c1 - 5)):
            j_candidates.append((p, score_junction(cleaned, p)))
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
        bscore = score_groove_ig_boundary(cleaned, boundary_pos, "I")
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
    mutations: str | Sequence[object] = (),
    features: Optional[SequenceFeatures] = None,
    use_early_shortcuts: bool = True,
) -> AlleleRecord:
    """Parse a class-II alpha chain via the shared grammar-driven engine.

    Enumerates a small top-K set of primary full parses, then lets them
    compete with fragment and missing-groove explanations.
    """
    grammar = grammar_spec("II", "alpha")
    flags: list[str] = []
    early_result, cleaned, _pairs, annotated, sp_est = prepare_parse_inputs(
        seq,
        grammar=grammar,
        allele=allele,
        gene=gene,
        features=features,
        use_early_shortcuts=use_early_shortcuts,
    )
    if early_result is not None:
        return _apply_requested_mutations(early_result, mutations)

    candidates: list[AlleleRecord] = []
    selections = collect_primary_parse_candidates(
        cleaned,
        grammar=grammar,
        annotated=annotated,
        sp_estimate=sp_est,
        features=features,
    )
    best_record = None
    expand_primary = True
    if selections:
        best_record = build_primary_result(
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
                    build_primary_result(
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
    return _apply_requested_mutations(_select_best_record_candidate(candidates), mutations)


# ---------------------------------------------------------------------------
# Class II beta groove parser
# ---------------------------------------------------------------------------


def decompose_class_ii_beta(
    seq: str,
    *,
    allele: str = "",
    gene: str = "",
    mutations: str | Sequence[object] = (),
    features: Optional[SequenceFeatures] = None,
    use_early_shortcuts: bool = True,
) -> AlleleRecord:
    """Parse a class-II beta chain via the shared grammar-driven engine.

    Enumerates a small top-K set of primary full parses, then lets them
    compete with β1-only salvage, fragment, and missing-groove
    explanations.
    """
    grammar = grammar_spec("II", "beta")
    flags: list[str] = []
    early_result, cleaned, _pairs, annotated, sp_est = prepare_parse_inputs(
        seq,
        grammar=grammar,
        allele=allele,
        gene=gene,
        features=features,
        use_early_shortcuts=use_early_shortcuts,
    )
    if early_result is not None:
        return _apply_requested_mutations(early_result, mutations)

    candidates: list[AlleleRecord] = []
    beta_resolver, beta1_ann = build_beta_secondary_resolver(cleaned, annotated, features)
    selections = collect_primary_parse_candidates(
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
        best_record = build_primary_result(
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
                record = build_primary_result(
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
    return _apply_requested_mutations(_select_best_record_candidate(candidates), mutations)


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
    species: str,
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
        species=species or None,
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
    species: Optional[str] = None,
    mutations: str | Sequence[object] = (),
    features: Optional[SequenceFeatures] = None,
    use_early_shortcuts: bool = True,
) -> AlleleRecord:
    """Dispatch groove parsing by class and chain.

    If ``mutations`` is provided, the groove is first extracted from the
    wild-type sequence (preserving Cys-pair detection), then the mutations
    are applied to the result.  Mutations use mature protein numbering
    (1-indexed) and can be a single string (``"K66A"``) or a sequence of
    strings, tuples, or mhcgnomes Mutation objects.
    """
    nc = normalize_mhc_class(mhc_class)
    normalized_mutations = _normalize_mutations(mutations)
    cleaned_seq = features.seq if features is not None else _clean_seq(seq)
    if features is None and not normalized_mutations:
        return _decompose_domains_cached(
            cleaned_seq,
            nc,
            str(chain or ""),
            allele,
            gene,
            str(species or ""),
            use_early_shortcuts,
        )

    gene_token = _normalized_gene_token(gene, allele)

    # Filter non-MHC genes (TAP1, CIITA, HM13, etc.).
    # mhcgnomes classification uses the species constraint when present;
    # the helper retains a small local set for unscoped non-MHC labels.
    if is_non_mhc_gene(gene or allele, species=species):
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

    # Infer class from the species-aware mhcgnomes classifier. Ambiguous
    # unscoped gene labels remain unresolved rather than using duplicate local
    # suffix rules that can disagree with the ontology.
    if nc not in ("I", "II") and (gene or allele):
        gc = parse_gene_class(gene or allele, species=species)
        if gc is not None:
            gc_class = gc.get("mhc_class")
            gc_chain = gc.get("chain")
            if gc_class == "I":
                nc = "I"
            elif gc_class == "II":
                nc = "II"
                if gc_chain in ("alpha", "beta"):
                    chain = gc_chain

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
        if normalized_mutations and result.ok:
            result = apply_mutations(result, normalized_mutations)
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
            except ValueError:
                pass
        result = _select_best_viable_record(candidates)
        if result is not None:
            if normalized_mutations and result.ok:
                result = apply_mutations(result, normalized_mutations)
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
    if normalized_mutations and result.ok:
        result = apply_mutations(result, normalized_mutations)
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
