# MHC Domain Parsing Status — 2026-03-27

## Overall Numbers

| Metric | Value |
|--------|-------|
| Ground truth sequences | 2,403 |
| Parsed (SP detected + groove decomposed) | 2,197 (91.4%) |
| Unparsed | 206 (8.6%) |
| **SP exact match** | **1,802 / 2,197 (82.0%)** |
| SP within ±1 | 1,911 / 2,197 (87.0%) |
| SP within ±2 | 1,975 / 2,197 (89.9%) |
| SP within ±3 | 2,027 / 2,197 (92.3%) |
| Negative controls: correct zero-SP | 1,951 / 2,155 (90.5%) |
| Negative controls: false positive SPs | 82 / 2,155 (3.8%) |
| Negative controls: abstentions | 122 / 2,155 (5.7%) |
| Tests passing | 318 |
| Evaluator runtime | ~25s for 2,403 sequences |

## Accuracy by MHC class × species

### Class I

| Species | Total | Exact | ≤±2 | ≤±3 |
|---------|-------|-------|-----|-----|
| Human | 130 | **99.2%** | 100.0% | 100.0% |
| NHP | 57 | **100.0%** | 100.0% | 100.0% |
| Bird | 232 | **94.4%** | 98.3% | 99.1% |
| Fish | 275 | **86.9%** | 92.4% | 93.5% |
| Other vertebrate | 645 | **68.2%** | 81.1% | 85.3% |
| Murine | 46 | 67.4% | 71.7% | 71.7% |
| Ungulate | 4 | 100.0% | 100.0% | 100.0% |

### Class II

| Species | Total | Exact | ≤±2 | ≤±3 |
|---------|-------|-------|-----|-----|
| Human | 112 | **88.4%** | 99.1% | 100.0% |
| NHP | 34 | **91.2%** | 100.0% | 100.0% |
| Bird | 150 | **89.3%** | 92.7% | 96.0% |
| Fish | 105 | **81.9%** | 87.6% | 88.6% |
| Other vertebrate | 339 | **86.1%** | 91.7% | 94.7% |
| Murine | 27 | 55.6% | 81.5% | 81.5% |
| Ungulate | 20 | 50.0% | 100.0% | 100.0% |

## Remaining failure modes

### 1. Unparsed (206 sequences)

| Status | Count | Notes |
|--------|-------|-------|
| missing_groove | ~119 | No viable Cys pair anchor; mostly divergent reptile/amphibian class I |
| alpha1_only | ~30 | True single-exon fragments (α1 only) — correctly handled |
| alpha2_only | ~17 | True single-exon fragments (α2 only) — correctly handled |
| fragment_fallback | ~14 | Short sequences — correctly handled |
| short | ~13 | Wrong Cys pair gives groove < 70 aa |
| non_classical | 4 | Danio L-lineage (now parsed but dispatch catches them) |
| other | ~9 | Edge cases |

### 2. Parsed but wrong SP (395 non-exact)

| Mode | Count | Description |
|------|-------|-------------|
| Near miss (±1-2) | 173 | Refinement picks adjacent valid site; +1 bias from preferring S/A over G |
| Overcall, sp_est also wrong | ~93 | Both h-region and parse overshoot; non-mammalian c-region grammar |
| Overcall, sp_est correct but ignored | ~35 | h-region knows the answer but groove-length prior overrides |
| Undercall | ~50 | SP predicted too short |
| Wrong parser selected | ~44 | Class I vs II β confusion in divergent vertebrates |

### 3. False positive SPs (82 controls)

| Species | Count | Notes |
|---------|-------|-------|
| Other vertebrate | 65 | Mostly reptile non-classical (E-S, F10 lineage) mature-only sequences |
| Fish | 9 | |
| Murine | 6 | |
| Bird | 1 | |
| Human | 1 | |

## Key architectural features

### Signal peptide detection
- **h-region detection**: windowed hydrophobicity scan with KD-based end trimming
- **SP grammar**: [n-region] → [h-region ≥8aa, hydro ≥0.50] → [c-region ≥3aa] → [cleavage]
- **Cleavage scoring**: von Heijne -3/-1 rules + property-based exclusion filter
- **Leaderless detection**: Met absence + weak h-region + charged N-terminus
- **Fast triage**: pattern-match top 10 boundary property hexamers (96.3% accuracy when firing)
- **Joint evidence**: Met × h-region quality × cleavage ordering — multiplicative gating

### Domain decomposition
- **Cys pair classification**: Trp41 at c1+14 (92.4% sensitivity, 96.5% specificity for C-like vs G-domain); Ile at c1+14 as tier-3 substitute for fish/reptile
- **Cys-flanking properties**: 3-aa property patterns before/after Cys pairs
- **Groove-length priors**: soft priors for α1/α2/β1 lengths, per-grammar
- **Junction motif**: G..H..Q for class I α1/α2 boundary
- **TM detection**: windowed hydrophobic scan, cached in SequenceFeatures

### Scoring architecture
- **Factored multiplicative scoring**: three factors (SP grammar, domain architecture, completeness) gate the additive total
- **Contradiction = gate**: contradictory evidence in any factor drives score toward 30% of additive
- **Missing = soft penalty**: missing evidence reduces but doesn't kill the parse
- **ParseScaffold**: pair-fixed terms computed once, only SP/length terms in the inner loop

## Performance architecture
- **SequenceFeatures**: all sequence-wide features precomputed once (Cys pairs, fold annotations, topology scores, h-region, TM span, prefix sums)
- **ParseScaffold**: pair-fixed scoring separated from start-dependent scoring
- **Suffix-best downstream support**: O(log P) lookup replaces O(P²) scan for class I α3 selection
- **Cached h-region**: threaded through all SP scoring functions

## Known gaps for future work

### Accuracy — parse competition
1. **Top-K parse competition** — retain multiple candidate parses per parser and compare whole architectures. The largest remaining error bucket is "the wrong whole explanation wins" (72 cases where class II β beats the correct class I parse). A small beam (K=3-5) with multiplicative compatibility scoring would address this directly.
2. **SP as ranked zone** — return top 3-5 cleavage candidates from infer_signal_peptide instead of a single point. 35 overcalls have the correct sp_estimate but it gets overridden by groove-length scoring at a single alternate position.
3. **Explicit fragment states as first-class outcomes** — 119 missing_groove unparsed sequences are mostly structured partial evidence that shouldn't be hard failures. Fragment/partial candidates should compete alongside full parses, with "missing evidence = soft penalty, contradictory evidence = gate."

### Accuracy — SP detection
4. **Murine class II** (55.6% exact) — specific investigation needed; may be allelic SP variation.
5. **Boundary-family models** — residue-class PWMs around -3..+3 for non-mammalian boundary families (e.g., ALA|RIP, LYQ|EKN, SEG|GKE). Each has 3-7 occurrences; residue-class level is the right granularity.
6. **Class I vs II β disambiguation** — multiplicative gating helps but structural disambiguation factors (expected groove length × support geometry × TM consistency) would help more.

### Remaining false positive SPs (82 cases)

The 82 mature-only controls falsely assigned an SP break down as:
- **64 other_vertebrate** (mostly reptile E-S/F10 non-classical genes), 9 fish, 6 murine
- **77/81 don't start with Met** — the Met gate works; these slip through via h-region
- **Root cause:** the h-region detector finds hydrophobic stretches that are internal domain features (groove α-helices), not SP cores. A real SP h-region is followed by a polar c-region and cleavage site — but so is a groove helix followed by a loop. The grammar can't fully distinguish them without knowing whether the sequence is full-length or mature-only.

Possible approaches:
- **Groove-helix fingerprinting:** the first helix of α1/β1 has a different residue composition than an SP h-region (more charged residues, aromatic residues). A trained classifier distinguishing "SP h-region" from "groove helix" could gate false positives.
- **Sequence-length-aware gating:** mature-only sequences are often shorter (no SP, no tail). If the sequence is 200-280 aa and starts without Met, the prior for "this is mature-only" should be much stronger.
- **N-terminal domain compatibility:** if the first 80 aa AFTER the phantom cleavage site look like a well-formed groove domain (good Cys pair, junction motif), that's evidence the protein starts at position 0, not at the phantom SP.

### Coverage
1. **mhcgnomes gaps** — 139 species prefixes (IPD-MHC standard 4-letter codes like Povi, Pagu, Euma) and ~80 gene suffixes (F10, BLB, Q9, E-S, DRA, A-U, DAA, etc.) not recognized by mhcgnomes even with species hint. See `data/mhcgnomes_gaps.md` for the full list with class/chain annotations. These are real IPD-MHC nomenclature, not invented.
2. **DM-specific grammar** — DMA has α1 ≈ 120 aa (vs standard 83-106). DM is cross-vertebrate (frog, caecilian, bird, mouse, human). Gene name "DMA"/"DMB" is conserved and could drive parser parameter selection.
3. **Teleost class I** — α1 domain often HAS a Cys pair (unlike mammals), causing pair selection confusion. The parser assumes α1 has no pair.
4. **Fragment handling** — 47 true fragments correctly identified as alpha1_only/alpha2_only; could be surfaced as partial-domain results with groove sequence for the available half.

### Performance
1. **Speed recovery if top-K widens the beam** — the ParseScaffold + SequenceFeatures infrastructure makes per-candidate scoring cheap once the beam is generated. The inner loop is pair-fixed terms (cached) + start-dependent terms only.
2. **Species-aware dispatch** using gene-name heuristics when mhcgnomes can't resolve — simple suffix matching (F10 → class I, BLB → class II β, DRA → class II α) would correctly dispatch most of the 139 unrecognized species.
