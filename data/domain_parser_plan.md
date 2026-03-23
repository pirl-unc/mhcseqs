# Plan: Domain Parser and Signal Peptide Refinement

## Context

The "groove parser" in mhcseqs decomposes MHC proteins into structural domains:
signal peptide → groove1 (α1) → groove2 (α2) → Ig domain (α3) → tail.

Currently the signal peptide boundary is inferred purely from a positional
heuristic (Cys-pair offset subtraction). This works well for mammals (~96% match
the AxA cleavage motif) but poorly for fish (23% match) and reptiles (30%).

This plan covers:
1. Renaming groove parser → domain parser
2. Extracting ground truth SP annotations from UniProt
3. Evaluating and improving SP boundary prediction
4. Auto-updating README counts

## Phase 1: Rename groove parser → domain parser

### What changes

- `groove.py` → rename functions but keep file name for now
  - `extract_groove()` → `decompose_domains()` (keep old name as alias)
  - `parse_class_i()` → `decompose_class_i()`
  - `parse_class_ii_alpha()` → `decompose_class_ii_alpha()`
  - `parse_class_ii_beta()` → `decompose_class_ii_beta()`
  - `_try_groove_parse()` in pipeline.py → `_try_domain_parse()`
- Keep `groove1`, `groove2` field names in AlleleRecord (they refer to the
  domain, not the parser)
- Update docstrings, comments, variable names

### Why

"Groove parser" implies it only extracts the binding groove. It actually does
full domain decomposition including SP inference, Ig domain, and tail. The
new name reflects what it actually does.

## Phase 2: Extract UniProt SP ground truth

### Approach

Fetch `Signal` feature annotations from UniProt for all reviewed (Swiss-Prot)
MHC proteins. These are experimentally verified or high-confidence SP
annotations with exact cleavage positions.

### Script: `scripts/fetch_sp_ground_truth.py`

1. Query UniProt: `(keyword:KW-0491) AND (reviewed:true) AND (ft_signal:*)`
2. For each entry: extract accession, organism, Signal start/end, sequence
3. Output: `data/sp_ground_truth.csv` with columns:
   - `accession`, `organism`, `sp_length` (UniProt annotation),
     `sequence`, `mhc_class`, `species_category`

### Expected yield

- Human HLA: ~100 reviewed entries with SP (A, B, C, DRA, DRB1, DQA1, etc.)
- Mouse H2: ~50 reviewed entries
- Rat RT1: ~30 reviewed entries
- Cattle/pig/horse: ~20 each
- Chicken: ~5
- Fish: ~0 (all unreviewed)
- Reptiles: ~0

**Limitation**: ground truth is mammal/bird-heavy. Fish and reptile SP
evaluation will need to use our Cys-pair predictions validated against
known sequence features (hydrophobic core, Met start) rather than
UniProt annotations.

## Phase 3: Evaluate SP prediction heuristics

### Current heuristic: Cys-pair offset

```
mature_start = raw_cys_position - expected_mature_cys_position
```

**Performance by clade** (from analysis of 36,974 SPs):

| Category | n | Starts M | -1 = A/G | -1 = A/G/S/T | AxA motif |
|---|---|---|---|---|---|
| human | 27,080 | 100% | 98.6% | 99.7% | ~97% |
| nhp | 6,514 | 100% | 96.2% | 98.1% | ~93% |
| carnivore | 22 | 100% | 100% | 100% | ~100% |
| ungulate | 667 | 100% | 85.2% | 91.8% | ~80% |
| bird | 410 | 100% | 85.9% | 92.0% | ~80% |
| murine | 424 | 100% | 57.3% | 62.3% | ~55% |
| other_mammal | 139 | 100% | 41.0% | 51.1% | ~40% |
| other_vertebrate | 218 | 100% | 30.3% | 45.0% | ~30% |
| fish | 888 | 100% | 23.3% | 43.8% | ~20% |

### Heuristics to evaluate

**Heuristic A: Cys-pair only (current)**
- No sequence analysis
- Pro: completely alignment-free, fast
- Con: off by 1–3 residues for many species, catastrophic for some

**Heuristic B: Cys-pair + AxA cleavage scan (mammals/birds only)**
- After Cys-pair offset, scan ±3 residues for nearest A/G/S at -1 position
- Pro: fixes ~91% of mammal/bird errors
- Con: doesn't help fish/reptiles where -1 residue is diverse
- Implementation: `refine_signal_peptide()` in groove.py (already drafted)

**Heuristic C: Hydrophobic core boundary**
- SP has n-region (charged, 1–5aa) → h-region (hydrophobic, 7–15aa) → c-region (polar, 3–7aa)
- Find the last residue of the hydrophobic core, then scan forward for
  the cleavage site
- Pro: sequence-based, works across clades
- Con: requires tuning hydrophobic threshold, harder to implement

**Heuristic D: Von Heijne weight matrix**
- Classic (-3, -1) weight matrix from von Heijne 1986
- Score each position in a window around the Cys-pair prediction
- Pick the highest-scoring position as the cleavage site
- Pro: well-validated, published, works for eukaryotes broadly
- Con: may need separate matrices for fish vs mammals

**Heuristic E: UniProt SP annotation passthrough**
- For entries with UniProt Signal features, use the annotated SP length directly
- Fall back to heuristic B/C/D for entries without annotations
- Pro: exact for reviewed entries
- Con: only covers mammals/birds in practice

### Evaluation plan

1. Build ground truth from UniProt Signal features (Phase 2)
2. For each heuristic A–E, predict SP length for ground truth entries
3. Compute: exact match %, mean absolute error, ±1 match %, ±3 match %
4. Stratify by species_category
5. Also evaluate: does the SP shift change groove1 boundaries? How much?

### Recommended approach

**Layered**: Use the best available method for each entry:

1. If UniProt has a Signal annotation → use it (exact)
2. If mammal/bird → Cys-pair + AxA scan (heuristic B)
3. If fish/reptile → Cys-pair only (heuristic A), no scan
4. Flag entries where no confident SP can be determined

## Phase 4: Auto-update README counts

### Script: `scripts/update_readme_counts.py` (rewrite)

The current script uses regex replacement and is fragile. Rewrite to:

1. Run after `mhcseqs build`
2. Read counts from built CSVs (requires build to exist)
3. Update the "Current data summary" section in README.md
4. Update species category table
5. Update source table with entry counts

### CI integration

Add to `.github/workflows/build-counts.yml`:

```yaml
- name: Update README counts
  run: python scripts/update_readme_counts.py

- name: Commit if changed
  run: |
    git diff --quiet README.md || (
      git add README.md
      git commit -m "Update README counts [ci]"
      git push
    )
```

This runs on every push to main. If counts change, it auto-commits.

## Implementation order

1. **Phase 4 first** (auto-update counts) — quick win, immediate value
2. **Phase 2** (fetch SP ground truth) — enables evaluation
3. **Phase 3** (evaluate heuristics) — data-driven SP improvement
4. **Phase 1 last** (rename) — cosmetic, can be done anytime

## Files to create/modify

| File | Change |
|---|---|
| `scripts/fetch_sp_ground_truth.py` | NEW — fetch UniProt SP annotations |
| `scripts/evaluate_sp_heuristics.py` | NEW — evaluate A/B/C/D/E on ground truth |
| `scripts/update_readme_counts.py` | REWRITE — auto-update README |
| `data/sp_ground_truth.csv` | NEW — UniProt SP annotations |
| `mhcseqs/groove.py` | ADD refine_signal_peptide(), rename functions |
| `mhcseqs/pipeline.py` | Call refine_signal_peptide() in SP detection |
| `.github/workflows/build-counts.yml` | ADD auto-commit README updates |
| `README.md` | Updated by script |
