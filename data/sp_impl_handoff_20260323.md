# SP Improvement Handoff For Claude

Date: 2026-03-23

This is the implementation handoff for the first SP-scoring patch. It is based
on the current repo state plus local evaluation and small controlled
simulations.

Related analysis:

- `data/sp_scoring_assessment_20260323.md`

## Goal

Implement the safest first-pass SP patch that:

- improves birds materially
- improves fish / other vertebrate modestly
- does not perturb Human / NHP yet
- removes the known alpha3-fallback catastrophic outliers

## Key Finding: Current Non-Mammal Scorer Is Often Bypassed

In `mhcseqs/groove.py`, `refine_signal_peptide()` currently returns early if
`seq[mature_start - 1]` is in `_SP_CLEAVAGE_RESIDUES`, before it checks whether
the sequence is mammalian or non-mammalian.

That means many non-mammal entries never reach `_score_nonmammal_candidate()`.

Local counts from `data/sp_ground_truth.csv`:

- birds: `341 / 421` (`81.0%`) are short-circuit eligible
- bird misses blocked from scoring: `232`
- fish: `158 / 367` (`43.1%`) short-circuit eligible
- other vertebrate: `536 / 1072` (`50.0%`) short-circuit eligible

So the first patch should not just tweak weights. It should make sure
non-mammals actually go through the scorer.

## Recommended Patch Scope

Edit only:

- `mhcseqs/groove.py`
- `tests/test_groove.py`

Do not change the mammal branch yet.

## Patch 1: Non-Mammals Should Always Score The Window

File:

- `mhcseqs/groove.py`

Current logic:

- unconditional early return if current `-1` residue is in `AGSC`

Change:

- keep the early-return behavior for mammals only
- for non-mammals, always evaluate the `+-5` scoring window

In practice:

1. In `refine_signal_peptide()`, move the existing early-return check into the
   mammal branch.
2. Leave mammal behavior otherwise unchanged.

Target behavior:

- mammals:
  - if current position is valid, keep it
  - else scan `+1, -1, +2, -2`
- non-mammals:
  - do not short-circuit on `AGSC`
  - always score the full local window

Why:

- this is the lowest-risk way to let the non-mammal scorer fix the exact
  bird/fish/reptile cases it was written for
- it preserves current Human / NHP behavior

## Patch 2: Improve Non-Mammal `-1` Residue Weights

File:

- `mhcseqs/groove.py`

Function:

- `_score_nonmammal_candidate()`

Change only the non-mammal scoring weights. Do not change
`_SP_CLEAVAGE_RESIDUES` globally in this first pass.

Recommended scoring:

- `A`: `+3.0`
- `G`: `+3.0`
- `S`: `+2.0`
- `C`: `+1.5`
- `T`: `+1.0`
- charged or `P`: `-3.0`
- other: `-0.5`

Why not change `_SP_CLEAVAGE_RESIDUES` yet:

- the mammal branch uses that set directly for short-circuiting and `+-2`
  scanning
- adding `T` globally would start changing Human / NHP behavior in the first
  patch, which is not desirable

## Patch 3: Add The Missing Bird Start Motifs

File:

- `mhcseqs/groove.py`

Table:

- `_MATURE_MOTIFS_CLASS_I`

Add:

- `ELH`, bonus `+2.5`
- `EPH`, bonus `+2.5`
- `EET`, bonus `+1.5`
- `AEL`, bonus `+1.5`
- `TRP`, bonus `+1.5`

Why these and not a broader motif expansion yet:

- these are strongly supported by the current bird misses
- they materially improve bird exact rate in local simulation
- fresh fish motif evidence is weaker / noisier than the earlier draft plan

Bird top missed tripeptides from the current repo state:

- `ELH`: `81`
- `EPH`: `56`
- `AGG`: `26`

Bird top exact tripeptides:

- `EEL`: `29`
- `TRP`: `17`
- `EET`: `17`
- `AEL`: `10`

## Patch 4: Apply `MAX_PLAUSIBLE_SP` To Class I Alpha3 Fallback

File:

- `mhcseqs/groove.py`

Function:

- `decompose_class_i()`

Problem:

- primary alpha2 path already enforces `MAX_PLAUSIBLE_SP`
- alpha3 fallback does not
- current evaluator still emits catastrophic outliers via `anchor_type =
  alpha3_cys`

Current examples:

- `A0A674HN04` (`Taeniopygia guttata`): `mature_start = 169`
- `A0A8D2LPA0` (`Varanus komodoensis`): `mature_start = 93`
- `A0A8D2LQE9` (`Varanus komodoensis`): `mature_start = 117`

Implementation:

After computing:

- `mature_start = _infer_mature_start(c1, CLASS_I_ALPHA3_CYS1_MATURE_POS)`

apply the same plausibility guard used in the alpha2 path:

- if `mature_start > MAX_PLAUSIBLE_SP`, return `status = "suspect_anchor"`
  with `anchor_type = "alpha3_cys"`

## What Not To Change In This Patch

Do not do these yet:

- unified scored approach for mammals
- changing `_SP_CLEAVAGE_RESIDUES` globally
- species-specific distance penalties
- wider than `+-5` window
- bird order overrides
- fish motif expansion beyond clearly validated cases
- class II constant changes

Those may be reasonable later, but they should not be bundled into the first
patch.

## Suggested Tests

File:

- `tests/test_groove.py`

Add tests for these behaviors:

### 1. Non-mammal path does not short-circuit on a merely valid residue

Add a direct unit test for `refine_signal_peptide()` using a synthetic sequence
where:

- the initial `mature_start` has a valid but wrong `-1` residue
- a nearby position has a better bird-style motif such as `ELH`

Expected:

- with `species_category="bird"`, refinement moves to the better site
- with `species_category="human"`, behavior remains conservative

### 2. Alpha3 fallback with implausible mature start returns suspect anchor

Add a synthetic class I sequence whose alpha3 fallback yields
`mature_start > MAX_PLAUSIBLE_SP`.

Expected:

- `decompose_class_i()` returns `status == "suspect_anchor"`
- `ok` is `False`

### 3. Bird motif regression test

Add a small direct test of `_score_nonmammal_candidate()` or
`refine_signal_peptide()` demonstrating that a candidate with `ELH` can beat a
nearby weaker alternative in a bird sequence.

If you prefer to avoid testing private helpers, use `refine_signal_peptide()`
only.

## Validation Command

Run:

```bash
python scripts/evaluate_sp_ground_truth.py
```

Also run targeted tests:

```bash
pytest tests/test_groove.py tests/test_pipeline.py
```

## Expected Outcome

From a local simulation of this exact patch shape:

- change set:
  - remove non-mammal short-circuit
  - improve non-mammal `-1` weights
  - add bird motifs `ELH/EPH/EET/AEL/TRP`
  - leave mammal branch unchanged

Expected exact-match movement:

| Category | Baseline exact | Simulated exact |
|---|---:|---:|
| Human | 84.7% | 84.7% |
| NHP | 78.7% | 78.7% |
| Bird | 29.2% | 34.7% |
| Fish | 33.5% | 34.9% |
| Other vertebrate | 37.9% | 42.5% |

Expected `<= 2 aa` movement:

| Category | Baseline <=2 | Simulated <=2 |
|---|---:|---:|
| Human | 90.6% | 90.6% |
| NHP | 92.6% | 92.6% |
| Bird | 74.8% | 76.2% |
| Fish | 56.9% | 59.4% |
| Other vertebrate | 61.3% | 57.9% |

Interpretation:

- bird exact improves materially
- fish improves modestly
- Human/NHP are preserved
- other vertebrate exact improves, though `<=2 aa` may trade off slightly

That tradeoff is acceptable for a first patch only if the evaluator confirms
the real code matches or beats the simulation.

## If You Need To Trim Scope Even Further

Minimum useful patch:

1. remove non-mammal short-circuit
2. add bird motifs `ELH` and `EPH`
3. cap alpha3 fallback by `MAX_PLAUSIBLE_SP`

That is the smallest patch I would still consider worth landing.
