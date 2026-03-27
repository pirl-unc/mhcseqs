# Signal Peptide Scoring Assessment

Date: 2026-03-23

## Executive Summary

Your core diagnosis is substantially correct.

The biggest currently validated failure mode in non-mammalian sequences is
the scorer's bias toward Gly at the `-1` cleavage residue. In birds
specifically, Ala-dominant true cleavage sites are frequently displaced by a
nearby Gly candidate inside the `+-5` refinement window.

That said, the full picture is slightly different from the draft plan:

- `refine_signal_peptide()` is already wired into the pipeline.
- `MAX_PLAUSIBLE_SP = 50` is already implemented in the primary parsers.
- The outlier problem is not fully solved because the class I alpha3 fallback
  can still return implausible `mature_start` values.
- Human/NHP class II errors are still consistent with missing or imperfect
  gene-specific constants and with parser misselection on atypical DM/DO
  molecules.
- The mammal/non-mammal split is itself now part of the problem: the mammal
  branch uses a short-circuiting `+-2` scan, which explains the NHP `+2`
  class I failures.

## Current Code State

Relevant implementation in `mhcseqs/groove.py`:

- `CLASS_I_ALPHA2_CYS1_MATURE_POS = 100`
- `_SP_CLEAVAGE_RESIDUES = AGSC`
- `_SP_CLEAVAGE_STRONG = AG`
- Non-mammal scoring gives:
  - `A/G` at `-1`: `+3.0`
  - `S/C` at `-1`: `+1.5`
  - distance penalty: `-0.5` per residue
- Class I mature start motifs currently include:
  - mammal: `GSH`, `CSH`, `GPH`, `SPH`
  - bird: `AAE`, `ASE`, `ASG`
  - fish: `AVT`, `KHS`
- Missing from the current motif table:
  - `ELH`, `EPH`, `EET`, `AEL`, `TRP`, `VKV`, `VDI`, `QPE`

Relevant implementation in `mhcseqs/pipeline.py`:

- `_signal_peptide_fields()` already calls `refine_signal_peptide()`

Important correction to the earlier notes:

- The older status note claiming SP refinement is "not yet wired into
  pipeline" is stale.

## Local Evaluation Snapshot

From `python scripts/evaluate_sp_ground_truth.py` on the current repo state:

- Parsed: `2324 / 2406` (`96.6%`)
- Exact SP length: `43.2%`
- Within `+-2 aa`: `67.9%`

By species category:

| Category | Exact | <=2 aa |
|---|---:|---:|
| Human | 84.7% | 90.6% |
| NHP | 78.7% | 92.6% |
| Bird | 29.2% | 74.8% |
| Fish | 33.5% | 56.9% |
| Other vertebrate | 37.9% | 61.3% |

These numbers match the broad pattern in your writeup: birds, fish, and
reptile/amphibian sequences remain the main problem.

## What Is Strongly Supported

### 1. Ala-at-`-1` is a real bird failure mode

Bird exact vs missed cases by true `-1` residue:

| `-1` residue | Exact | Miss |
|---|---:|---:|
| `G` | 102 | 35 |
| `A` | 21 | 185 |
| `S` | 0 | 63 |

This is the strongest single signal in the data.

Among bird misses:

- `185 / 298` misses have true `-1 = A`
- `172 / 185` of those Ala misses still have at least one nearby Gly-based
  candidate inside the current refinement window

That is exactly the failure shape you described: the scorer is often finding a
nearby Gly and preferring it over the correct Ala cleavage site.

### 2. Bird mature-start motif coverage is incomplete

Bird start tripeptides among exact calls:

- `EEL`: 29
- `TRP`: 17
- `EET`: 17
- `AEL`: 10
- `KET`: 10

Bird start tripeptides among misses:

- `ELH`: 81
- `EPH`: 56
- `AGG`: 26
- `VKV`: 11

This strongly supports adding `ELH` and `EPH` first. Those are not edge
motifs; they are dominant missed motifs.

### 3. NHP class I really does have a `+2` cluster

NHP delta distribution from the current evaluator:

- `0`: 74
- `+2`: 8
- `-2`: 4

Broken down by parser:

- `7` of the `+2` cases are `class_I`

This is consistent with the current mammalian path:

- if the predicted site already lands on any valid `AGSC` residue, the code
  returns immediately
- otherwise mammals only scan `+1, -1, +2, -2`

So a slightly wrong but still "valid" mammalian prediction can survive
unchallenged.

### 4. Human/NHP class II still has gene-specific constant problems

Current human mismatches greater than `3 aa` include:

- `P01906`: `+5`, `class_II_alpha`
- `P20036`: `-7`, `class_II_beta`
- `Q14954`: `+12`, `class_II_beta`
- `Q99706`: `+6`, `class_II_beta`
- `P43626`: `+12`, `class_II_beta`
- `P43628`: `+12`, `class_II_beta`

This is aligned with the current constant tables:

- class II alpha overrides only exist for `DQA` and `DMA`
- class II beta overrides only exist for `DPB`

So the claim that class II non-classical genes remain under-parameterized is
supported.

## What Is Already Implemented

Two pieces of the draft plan are already present in code:

### 1. Pipeline integration

`refine_signal_peptide()` is already called from the pipeline, so the current
bad bird/fish performance is with refinement enabled, not without it.

### 2. Plausible-SP cap

`MAX_PLAUSIBLE_SP = 50` is already enforced in the primary class I alpha2,
class II alpha, and class II beta anchor paths.

## What Is Not Fully Solved By Current Code

### 1. `MAX_PLAUSIBLE_SP` does not stop all catastrophic class I outliers

Example from the current evaluation:

- `A0A674HN04` (`Taeniopygia guttata`)
- ground truth SP: `23`
- predicted: `171`
- delta: `+148`

Why this still happens:

- class I alpha2 path is capped
- but the class I alpha3 fallback can still return `status = inferred_from_alpha3`
  with `mature_start = 169`

So the cap exists, but not on every fallback path that can emit an SP length.

### 2. The mammal/non-mammal split is now a liability

The mammal branch is intentionally simple and conservative, but it also means:

- wrong predictions that already land on `A/G/S/C` are preserved
- NHP class I `+2` errors survive
- mammals do not benefit from motif competition or distance-aware rescoring

So the proposed unified scored approach is defensible, although I would still
roll it out carefully behind evaluation rather than replacing the current
branch in one shot.

## Overall Assessment Of The Hypothesis

I would summarize your hypothesis like this:

- Correct in substance
- Slightly outdated in implementation details
- Strongest on birds
- Also consistent with fish and reptile/amphibian misses
- Incomplete for catastrophic fallback outliers

The single most compelling validated claim is:

> promoting Ala and adding missing bird start motifs should recover a large
> fraction of current bird mismatches

That is a much stronger conclusion than "possible explanation"; the local data
supports it directly.

## Recommended Priority Order

### Tier 1: Highest confidence, highest leverage

1. Promote non-mammal `A` at `-1`
2. Add bird motifs: `ELH`, `EPH`, `EET`, `AEL`, `TRP`
3. Add fish motifs: `VKV`, `VDI`, `QPE`
4. Add `T` to the valid `-1` residue set or otherwise score it positively

Reason:

- these changes directly target the dominant observed miss patterns
- they do not require taxonomy refactors
- they are easy to validate with the existing evaluator

### Tier 2: Fix obvious structural gaps

1. Apply plausible-SP capping to class I alpha3 fallback as well
2. Revisit human/NHP class II gene-specific constants
3. Audit parser selection for DM/DO-like atypical molecules

Reason:

- these address large deltas and misparser failures
- they improve correctness, not just exact-match rate

### Tier 3: Broader architecture improvements

1. Replace mammalian short-circuiting with scored refinement
2. Use species-tuned distance penalties
3. Widen the non-mammal window from `+-5` to `+-8` where justified
4. Add bird/reptile order-level Cys1 overrides

Reason:

- these are plausible but riskier
- they should come after the quick scoring fixes, not before

## Bottom Line

If the goal is the next highest-yield SP improvement, I would not start with a
full redesign.

I would start with:

1. Ala promotion
2. `ELH` / `EPH` motif support
3. class I alpha3 fallback capping
4. class II constant cleanup

That sequence matches the current evidence best and should produce measurable
gains quickly.
