# Issue #15: whole-parse beam competition

## Baseline

- Current branch: `feature/top-k-whole-parse` from v2.6.1.
- The enriched 2,402-record benchmark, deliberately dispatched without class or
  gene metadata, has 405 wrong labeled classes and 584 wrong class/chain
  architectures. Twenty-one class-I records are selected as class-II beta.
- Each grammar already keeps two structurally distinct primary anchor
  candidates and `ParseScaffold.score()` already applies multiplicative gates.
- `_build_primary_result()` discards the selected candidate's multiplicative
  score and materializes only the additive component total. Final record and
  inter-parser ranking therefore do not compare the scores that formed the
  beam.

## Design

1. Carry the complete `PrimaryParseSelection.selection_score` into the
   materialized record's `candidate_score`, while preserving additive
   `parse_score` and component fields for diagnostics.
2. Widen the primary beam to four candidates, as requested by issue #15, and
   retain structurally distinct anchor/start hypotheses.
3. Add regression tests proving that final record selection uses the
   multiplicative/contradiction-aware score and that the beam can retain four
   distinct whole parses.
4. Benchmark classless class/chain dispatch, signal-peptide accuracy, and
   mature-only false positives before accepting the change. Do not ship a
   scoring change that improves one named bucket by silently regressing the
   broader scientific controls.
5. Bump the patch version, run lint and the full test suite, open the PR, wait
   for green CI, merge it, and deploy the release from clean `main`.

## Verification

- [x] Focused beam and score-propagation regressions
- [x] `scripts/evaluate_sp_ground_truth.py` (2,203 parsed, 1,853 exact;
      baseline 2,193 parsed, 1,851 exact)
- [x] `scripts/evaluate_sp_negative_controls.py` (32 false positives;
      baseline 34)
- [x] `./lint.sh`
- [x] `./test.sh` (375 passed)
- [ ] GitHub CI green
- [ ] PyPI release import/version verified

## Follow-up

- Filed #52 for the 241 unresolved/inferred benchmark rows and the evaluator's
  duplicate closest-to-23-aa parser selection heuristic.
