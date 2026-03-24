# SP Structural Parser: Next Steps

Date: 2026-03-24

## Current state

Committed on `feature/sp-pipeline-integration`:
- Property-based SP scoring (universal across species)
- Class I Cys pair classification by flanking motifs (replaces position ranges)
- Junction motif discovery (G..H..Q finds α1/α2 boundary)
- Adaptive-width SP search (±8 base, ±15 with junction)
- Parse trace for debugging
- Class II parsers remain on position-range Cys selection

Accuracy: bird 62.1%, fish 59.6%, human 84.7%, overall 56.6% exact.

## Next: Full candidate parse enumeration

### Problem

The class II motif-based Cys classification (calibrated tables ready in
groove.py) causes regressions when applied directly because it sometimes
selects different pairs than the position-range approach. Without holistic
scoring of the full `[SP]-[groove]-[Ig]-[tail]` decomposition, there's no
way to tell which pair selection produces a more plausible overall parse.

### Architecture

Add `_enumerate_parses()` and `_score_parse()` that:

1. For each Cys pair assignment (motif-based AND position-based candidates):
   a. Compute domain boundaries
   b. Score SP boundary (existing _score_sp_candidate)
   c. Score junction motifs (existing _score_junction)
   d. Score groove/Ig boundary motifs (existing _score_groove_ig_boundary)
   e. **Score domain size plausibility** (new):
      - SP: 14-42 standard → bonus, 0 (stripped) → neutral, >42 → penalty
      - groove1: 75-95 standard → bonus, 70-125 extended → neutral, else → penalty
      - groove2: 85-98 standard → bonus, 80-110 extended → neutral, else → penalty
      - Ig: 75-100 → bonus, else → penalty
      - Total protein = SP + groove + Ig + tail: class I ~340-370, class II ~250-270
   f. **Score structural ordering** (new): domains must appear in expected order
   g. Sum all component scores

2. Pick the candidate with highest total score
3. Convert to AlleleRecord

### Key insight

The domain size scoring prevents impossible parses:
- Can't parse 40aa SP from 90aa fragment (leaves 50aa groove = too short)
- Can't parse 200aa groove (never happens)
- Can't have Ig before groove (wrong order)
- SP + groove + Ig + tail should sum to ≈ full protein length

### Implementation plan

1. Add `DomainSpan` and `CandidateParse` dataclasses (designed in plan, not yet implemented)
2. Add `_score_domain_lengths()` — soft penalties for implausible sizes
3. Add `_enumerate_class_i_parses()` — generates candidates from all valid Cys assignments
4. Add `_enumerate_class_ii_alpha_parses()` and `_enumerate_class_ii_beta_parses()`
5. Add `_score_parse()` — holistic score combining all signals
6. Wire into existing `decompose_class_i/ii` as a new code path
7. Gate on evaluator: no regression for any category

### Class II motif classification (unblocked by enumeration)

With holistic scoring, motif-based Cys selection can compete with position-based:
- Generate candidates from BOTH motif and position pair selections
- Score each holistically
- Best one wins

The calibrated class II flanking tables (already in groove.py):
- β1 groove: E@-1 (66%), K@-3 (52%), N@+3 (81%), Y@+4 (81%)
- β2 Ig: Y@c2-2 (64%), V@c2+2 (74%), H@c2+4 (73%), W@c2+15 (72%)
- α Ig: G/I@c1-1, L@c2+15 (60%)

### Gene-specific constant elimination

With enumeration + holistic scoring + junction discovery:
- Class I: constants already bypassed by motif classification
- Class II alpha: junction discovery doesn't apply (no groove1/groove2 split),
  but domain size scoring can replace the gene-specific Cys1 mature positions
  (DQA=109, DMA=120) by preferring candidates with plausible groove1 lengths
- Class II beta: similar — score groove2 length instead of using DPB=114

## Next: Remaining accuracy gaps

### Human class II (84.7% → target >95%)

39 human misses are almost all class II gene-constant issues:
- 12x DQB alleles with SP=30 (vs 32 for some alleles) — allelic SP variation
- 6x DQA/DMA — evaluator lacks gene context (pipeline has it)
- 4x entries with SP misparse
- Fix: pass gene context to evaluator, add DQB allelic constant

### NHP class I +2 cluster (78.7% → target >85%)

8 NHP class I entries with +2 offset — the mammal short-circuit early-returns
on a valid cleavage residue 2 positions from the correct site. Fix: use scored
approach for NHP (relax mammal short-circuit for NHP only, or for all mammals).

### Bird remaining misses (62.1% → target >75%)

~160 bird misses within ±8 of Cys prediction — the scorer picks nearby
positions with good von Heijne signals but wrong motifs. Many are the
`...GAAA|ELH...` pattern where AAE at -2 competes with ELH at the correct
position. The property-based scorer partially solved this (was 29.2%), but
the remaining cases need the junction confirmation in the wider ±15 range
to fully resolve.

### Fish (59.6% → target >70%)

Fish have 8 divergent alpha1 lineages with varying lengths. The junction
motif (G..H..Q) is less conserved in fish. May need fish-specific junction
motifs or a wider search with weaker junction threshold.

## Next: Class inference from sequence

When MHC class is unknown, infer from Cys pair classification:
- If best groove pair matches class I α2 flanking (G@-1, M@-3) → class I
- If best groove pair matches class II β1 flanking (E@-1, K@-3) → class II beta
- If only Ig-type pairs found → class II alpha
- If ambiguous → try both and pick best holistic parse

This would let the evaluator work without needing to try all three parsers.
