# Signal Peptide Prediction: Failure Analysis and Improvement Plan

Evaluated against `data/sp_ground_truth.csv` (2,406 UniProt entries with
SignalP-annotated or experimentally validated signal peptides, fetched
March 2026).

## Current architecture

```
Full sequence → Cys-pair finder → mature_start estimate
                                        ↓
                            refine_signal_peptide()
                            ├─ mammals:     ±2 simple A/G/S/C scan
                            └─ non-mammals: ±5 scored window (7 signals)
                                        ↓
                                  refined SP length
```

The Cys-pair finder locates conserved disulfide bonds in Ig-fold domains and
back-calculates `mature_start` using hardcoded expected mature positions
(e.g. class I α2 Cys1 at mature position 100).  `refine_signal_peptide()`
then adjusts the cleavage site within a small window.

## Current accuracy

| Category | n | Exact | ±2 | ±5 |
|---|---|---|---|---|
| Human (class I only) | 129 | 99.2% | 100% | 100% |
| Human (all) | 255 | 84.7% | 90.6% | — |
| NHP | 94 | 78.7% | 92.6% | — |
| Bird | 421 | 29.2% | 74.8% | — |
| Fish | 367 | 33.5% | 56.9% | — |
| Other vertebrate | 1072 | 37.9% | 61.3% | — |

---

## Category-by-category failure analysis

### 1. Human class II (39 mismatches out of 260)

**Root cause: Cys-pair constant imprecision for non-classical class II genes.**

The class II parser uses gene-specific Cys1 mature positions:

| Gene | Current constant | Effective SP offset from GT |
|---|---|---|
| DRA (default) | 106 | +6 (P01903: gt=25, pred=31) |
| DQA | 109 | +5 (P01906: gt=23, pred=28) |
| DMA | 120 | +11 (P28068: gt=18, pred=29) |
| DPA | (uses default 106) | -7 (P20036: gt=28, pred=21 via IIb) |
| DOA | (no constant) | +12 (Q14954: gt=21, pred=33) |
| DOB | (no constant) | — |
| DPB | 114 | — |

Key issues:
- **HLA-DM** (DMA/DMB): atypical Ig-fold with intra-domain Cys pair. DMA constant
  120 was set but the IIa parser still lands 11 aa past GT. The predicted site
  `...ESTC|LLDD...` vs GT `...CTGA|GGFV...` — the C at -1 of our predicted site
  is actually the start of the Ig-domain Cys pair, not a cleavage residue.
- **HLA-DO** (DOA/DOB): no gene-specific constant. Parser falls back to default
  and overshoots by 12.
- **HLA-DPA1**: no gene-specific constant for DPA; the evaluation picks the IIb
  parser (closest to typical SP=23) which gives 21 instead of the correct 28.

**Fix**: Add gene-specific Cys1 constants for DPA, DOA, DOB, DMB. For DM/DO,
consider whether the standard Cys-pair approach works at all — these non-classical
molecules may need a separate code path.

### 2. NHP (23 mismatches out of 94)

**Two distinct failure modes:**

**a) Class I: systematic +2 offset (7 entries)**

All share the same pattern: GT cleavage at `...TWA|GSH...` (position 21),
but we predict `...AGS|HSL...` (position 23).

Root cause: `CLASS_I_ALPHA2_CYS1_MATURE_POS = 100` is 2 less than the actual
Cys1-to-mature-start distance in these NHP class I sequences (~102). The Cys
prediction lands at 23, and since `seq[22] = 'S'` is a valid cleavage residue,
the simple ±2 scan stops immediately — it doesn't even try other positions.

Impact: 7/94 NHP entries = 7.4%. All are class I with this specific motif.

**b) Class II: same gene-specific constant issues as human (16 entries)**

DM, DO molecules from gorilla, chimp, macaque show identical offset patterns
to human counterparts.

**Fix for (a)**: Either:
- Adjust `CLASS_I_ALPHA2_CYS1_MATURE_POS` from 100 to 101 (affects all species)
- Or: use the scored approach for mammals too, with `GSH` and `TWA` motifs
  competing — the scorer could recognize that `TWA|GSH` is a better cleavage
  site than `AGS|HSL` because T at -3 is small/aliphatic and W is hydrophobic.

### 3. Birds (302 mismatches out of 421)

**Three distinct failure modes:**

**a) Systematic -2 offset for Gallus/Anas (dominant, ~200 entries)**

Gallus (n=118): 8% exact, 87% ±2, median delta = **-2**
Anas (n=91): 11% exact, 78% ±2, median delta = **-2**

The Cys-pair prediction is close (within 2) but the ±5 scored refinement
picks a position 2 aa before the GT cleavage site.

**Motif evidence** — the -1 residue at the GT cleavage site tells the story:

| -1 residue | Exact | Miss |
|---|---|---|
| G | 103 | 34 |
| A | 14 | **184** |
| S | 0 | **69** |

When the correct -1 is **Gly**, we find it. When it's **Ala** or **Ser**, we
almost always miss — the scorer systematically prefers Gly (score +3.0) over
Ala/Ser (score +1.5). With Cys prediction 2 past the GT site, the scorer
finds a Gly at -3/-4 from GT and scores it higher than the Ala at -1.

**Mature start tripeptide analysis confirms this:**

| Tripeptide | Exact | Miss | Total | Notes |
|---|---|---|---|---|
| ELH | 3 | **81** | 84 | Dominant Gallus/Anas class I start. -1 is A (cleavage after A). |
| EPH | 0 | **64** | 64 | Common Anas class I variant. -1 is A. |
| AGG | 0 | **26** | 26 | -1 is A. |
| EEL | **29** | 2 | 31 | -1 is G. We find these because G scores higher. |
| TRP | **21** | 8 | 29 | -1 is A, but AAE motif bonus helps. |
| KET | **10** | 0 | 10 | -1 is G. |

**Fix**: The scoring function `_score_nonmammal_candidate` gives Gly +3.0
but Ala only +1.5 at -1. For birds, Ala is the dominant cleavage residue.
Solutions:
1. Add `ELH`, `EPH` as bird class I mature start motifs with bonus ≥ 2.0
2. Increase the Ala score at -1 for birds (or globally, since Ala is the
   most common SP cleavage residue in SignalP databases)
3. Add a **"consensus cleavage site"** signal: Ala-X-Ala at -3/-1 (the
   von Heijne AXA motif) is the canonical SP cleavage pattern

**b) 15-aa α1 insertion for passerines/strigiformes (~50 entries)**

Strix (n=31): 6% exact, median delta +3, but ±5 = 100%
Taeniopygia (n=21): 71% exact, but outliers at +148
Apteryx (n=20): 55% exact, avg delta +3.0
Serinus (n=6): 0% exact, avg delta -3.2

These birds have Cys1 at mature position ~110-119 instead of 100. The Cys
prediction overshoots by 10-19 residues. The ±5 window can partially
compensate (Strix gets 100% ±5) but can't fix larger offsets.

**Fix**: Add genus-to-order mapping in species.py, then use order-level
Cys1 mature position overrides:
- Strigiformes (owls): 119
- Passeriformes subset (Taeniopygia, Serinus, Lonchura, Catharus): 115
- Phasianidae with insertion (some Phasianus): 116

**c) Catastrophic outliers (A0A674HN04: delta +148)**

A handful of entries where the Cys-pair finder latches onto a completely
wrong Cys pair deep in the sequence. These are likely non-classical class I
genes or misannotated entries.

**Fix**: Add a MAX_PLAUSIBLE_SP cap (e.g. 50 aa). If predicted mature_start
> 50, flag as suspect rather than emitting a wrong boundary.

### 4. Fish (244 mismatches out of 367)

**Cleavage site motif analysis:**

| -1 residue | Exact | Miss |
|---|---|---|
| A | 87 | **188** |
| S | 32 | **64** |
| G | 20 | 21 |
| T | 2 | **18** |
| P | 0 | **21** |

Same pattern as birds: Ala-dominant cleavage sites are systematically missed.
Additionally, **Thr** (18 misses) and **Pro** (21 misses) appear at -1 in
the ground truth — Pro is never a valid cleavage residue in our model, and
Thr isn't in `_SP_CLEAVAGE_RESIDUES` (`AGSC`).

**+1 residue (mature start):**

| +1 | Exact | Miss |
|---|---|---|
| V | 28 | **106** |
| G | 3 | **57** |
| Q | 5 | **26** |
| E | 56 | 39 |
| D | 14 | 31 |

Val, Gly, and Gln at +1 are strongly associated with misses. These are
fish-specific mature protein starts not in our motif table.

**Per-genus breakdown:**

| Genus | n | Exact | ±2 | Issue |
|---|---|---|---|---|
| Triakis (shark) | 35 | 0% | 0% | Systematic -3.0 delta; shark Cys positions differ |
| Sparus (sea bream) | 22 | 5% | 27% | Cys1 variable, some non-classical |
| Oreochromis (tilapia) | 42 | 21% | 38% | Cichlid-specific α1 length |
| Gadus (cod) | 26 | 50% | 92% | Decent — cod MHC is well-characterized |
| Salmo (salmon) | 36 | 39% | 86% | Decent |
| Acipenser (sturgeon) | 25 | 60% | 64% | Basal fish, Cys1 variable |

**Fix**:
1. Add T (threonine) to `_SP_CLEAVAGE_RESIDUES` — Thr at -1 is valid in
   SignalP for ~8% of fish MHC SPs
2. Add fish-specific mature motifs: VKV, VDI, QVE, etc.
3. For sharks (Triakis, Scyliorhinus): the α2 Cys-pair distance may differ
   from bony fish. Investigate whether cartilaginous fish need separate
   Cys1 constants.
4. Widen the scanning window from ±5 to ±8 for fish, with gentler distance
   penalty (0.3 per residue instead of 0.5)

### 5. Other vertebrate — reptiles/amphibians (698 mismatches out of 1072)

**Hydrophobic core density does NOT discriminate:**
- Exact matches: mean h-region hydrophobicity = 0.638
- Misses: mean h-region hydrophobicity = 0.688

This means signal 4 (hydrophobic core density) is essentially noise for
reptiles/amphibians. Their signal peptides have consistently high h-region
density whether the prediction is right or wrong.

**Delta distribution:**
- Peak at 0 (274 exact), then -3 (53 cases) and +10 (23 cases)
- 21.2% off by >5 aa, 6.8% off by >10 aa
- Asymmetric: negative deltas cluster tightly around -3, while positive
  deltas are spread with a bump at +9/+10

**Fix**: The +9/+10 bump (40 entries) is likely the bird-like α1 insertion
in some squamate reptiles (Varanus, Pogona). Same fix as birds: genus/order
Cys1 overrides. The -3 cluster suggests the scored refinement is consistently
falling 3 short — likely the same Ala-at-minus-1 issue as birds.

---

## Proposed new scoring architecture

### Core idea: Unified scored approach for all species

Replace the mammal/non-mammal split with a single scored approach that works
across all species. The current ±2 simple scan for mammals is fast but misses
cases where the Cys prediction is already at a "valid" residue but wrong
(NHP +2 problem). A unified scorer with species-tuned weights solves both.

### Proposed scanning window

| Category | Current | Proposed | Distance penalty |
|---|---|---|---|
| Mammals | ±2 simple | ±5 scored | 0.8/residue (tight) |
| Birds | ±5 scored | ±8 scored | 0.3/residue (loose) |
| Fish | ±5 scored | ±8 scored | 0.3/residue |
| Reptile/amphibian | ±5 scored | ±8 scored | 0.3/residue |

### Proposed scoring signals (expanded from current 7)

```
Signal 1: -1 cleavage residue (REVISED)
  Current:  G/A → +3.0,  S/C → +1.5
  Proposed: A   → +3.0  (promote — Ala is the most common SP -1 across all clades)
            G   → +3.0
            S   → +2.0  (promote — common in fish/bird)
            C   → +1.5
            T   → +1.0  (NEW — valid in ~8% of fish SPs)
  Charged/P at -1: -3.0 (unchanged)

Signal 2: -3 residue (REVISED)
  Current:  small/aliphatic at -3 → +1.0
  Proposed: Add AXA motif detection:
            If -1 ∈ {A,G,S} AND -3 ∈ {A,V,S,T,I,L} → +2.0 (von Heijne motif)
            This is the canonical eukaryotic SP cleavage pattern.

Signal 3: +1 not proline (unchanged)
  P at +1 → -2.0

Signal 4: Upstream hydrophobic density (REVISED)
  Current:  h-region density in [pos-12:pos-3], threshold 0.4
  Problem:  Doesn't discriminate for reptiles (both hits and misses are ~0.65)
  Proposed: Keep but reduce weight from 2.0 to 1.0 for non-mammals.
            Add minimum absolute density check: if < 0.3, score -2.0
            (catches cases where "SP" is actually just mature protein)

Signal 5: c-region polarity (unchanged)
  Hydrophobic at [pos-3:pos] → -1.5

Signal 6: Distance from Cys prediction (REVISED)
  Current:  -0.5 per residue (all species)
  Proposed: Species-tuned:
            Mammals:    -0.8/residue (tight, Cys is reliable)
            Birds:      -0.3/residue (loose, α1 insertions shift it)
            Fish:       -0.3/residue (loose, variable α1)
            Reptile:    -0.4/residue (moderate)

Signal 7: Mature protein start motif (EXPANDED)
  Current motifs:
    GSH +1.5 (mammalian class I)
    CSH +1.0, GPH +1.0, SPH +0.8 (mammalian variants)
    AAE +1.5, ASE +1.0, ASG +0.8 (avian class I)
    AVT +1.0, KHS +1.0 (fish class I)

  Proposed additions from ground truth analysis:
    # Bird class I (from tripeptide data — these are the top missed motifs)
    ELH +2.5   # dominant Gallus/Anas class I start (84 occurrences, 3 exact)
    EPH +2.5   # common Anas class I variant (64 occurrences, 0 exact)
    EEL +1.5   # already well-predicted but boost for consistency
    EET +1.5   # bird class I variant
    AEL +1.5   # bird class I variant
    TRP +1.5   # bird class I variant

    # Fish class I (from +1 residue data — V/G/Q dominant in misses)
    VKV +1.5   # fish class I start
    VDI +1.0   # fish class I variant
    QPE +1.5   # fish class I variant

    # Reptile/amphibian (needs further analysis, placeholders)
    EEA +1.0
    KHN +1.0

Signal 8 (NEW): Downstream charged residue density
  Mature MHC protein starts typically have 1-2 charged residues (D, E, K, R)
  in positions +1 to +5. This distinguishes real cleavage sites from
  positions still inside the hydrophobic signal peptide.
  Score: count charged in [pos:pos+5] / 5
         if > 0.3 → +1.0
         if < 0.1 → -1.0

Signal 9 (NEW): N-terminal Met check
  If the sequence starts with M and pos < 12, it's likely inside the SP
  n-region, not the cleavage site.
  Score: if pos < 12 → -2.0
```

### Cys1 mature position improvements

**Class I:**
- Bump `CLASS_I_ALPHA2_CYS1_MATURE_POS` from 100 → 101 (fixes NHP +2 bias,
  verified against human class I which stays at 100% ±2)
- Add order-level overrides for birds:
  ```python
  _CLASS_I_CYS1_MATURE_POS_BY_ORDER = {
      "strigiformes": 119,
      "passeriformes_inserted": 115,  # Taeniopygia, Serinus, Lonchura, Catharus
  }
  ```
  Requires genus → order mapping in species.py.

**Class II:**
- Add missing gene constants:
  ```python
  _CLASS_II_ALPHA_CYS1_MATURE_POS_BY_GENE["DPA"] = 104   # needs validation
  _CLASS_II_ALPHA_CYS1_MATURE_POS_BY_GENE["DOA"] = 112   # needs validation
  ```

### MAX_PLAUSIBLE_SP cap

Add a hard cap: if `refined_mature_start > 50`, return 0 (treat as no SP
detected). No vertebrate MHC signal peptide exceeds ~40 aa. This eliminates
catastrophic outliers like the Taeniopygia entry at +148.

---

## Expected impact

| Category | Current exact | Projected exact | Key fix |
|---|---|---|---|
| Human (all) | 84.7% | ~93% | Class II gene constants |
| NHP | 78.7% | ~90% | Cys1 101 + class II constants |
| Bird | 29.2% | ~65% | ELH/EPH motifs + Ala promotion + order overrides |
| Fish | 33.5% | ~55% | Ala promotion + T at -1 + fish motifs + wider window |
| Other vertebrate | 37.9% | ~55% | Ala promotion + wider window + order overrides |

The biggest single improvement is **promoting Ala to +3.0 at -1** and adding
**ELH/EPH** as bird mature motifs — together these address ~250 of the ~800
total mismatches (birds with Ala at -1 whose start tripeptide is ELH/EPH).

## Implementation order

1. **Quick wins** (30 min each, independent):
   - Add MAX_PLAUSIBLE_SP=50 cap
   - Add ELH, EPH, EET, AEL bird motifs
   - Add T to `_SP_CLEAVAGE_RESIDUES`
   - Add class II gene constants (DPA, DOA)

2. **Scoring changes** (needs careful regression testing):
   - Promote Ala to +3.0 at -1 (currently +1.5 via `_SP_CLEAVAGE_RESIDUES`)
   - Promote Ser to +2.0
   - Species-tuned distance penalty
   - AXA motif detection at -3/-1

3. **Structural changes** (needs species.py expansion):
   - Bird order Cys1 overrides
   - Genus → order mapping for ~30 bird genera

4. **Architecture change** (biggest risk, biggest reward):
   - Unified scored approach for mammals
   - Wider ±8 window for non-mammals
   - Signal 8 (downstream charged density) and Signal 9 (Met check)

Each step should be validated against the full ground truth before proceeding
to the next. The evaluation script (`scripts/evaluate_sp_ground_truth.py`)
enables rapid iteration.
