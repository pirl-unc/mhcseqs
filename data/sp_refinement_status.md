# Signal Peptide Refinement: Status and TODO

## What's implemented (this PR)

### `refine_signal_peptide()` in groove.py

Hybrid approach:
- **Mammals**: simple ±2 scan for nearest A/G/S/C at -1 position
- **Non-mammals**: multi-signal scoring over ±5 window combining:
  1. -1 cleavage residue (A/G strong, S/C weaker)
  2. -3 residue (small/aliphatic preferred)
  3. +1 not proline
  4. Upstream hydrophobic density (h-region detection)
  5. c-region polarity (hydrophobic→polar transition)
  6. Distance from Cys-pair prediction (gentle anchor, 0.5 per residue)
  7. Mature protein start motifs (GSH, AAE, AVT, KHS, etc.)

### Performance (class I alpha)

| Group | Cys-pair only | Simple ±2 | Hybrid | n |
|---|---|---|---|---|
| Human | 99.6% | 100.0% | 100.0% | 20,618 |
| Other mammal | 96.2% | 99.7% | 99.7% | 5,620 |
| Non-mammal | 62.9% | 92.2% | 96.5% | 1,037 |

### `parse_allele_name()` species threading

`alleles.py` now accepts optional `species=` parameter, threaded to
`mhcgnomes.parse()` for disambiguation. Auto-detects `species=` vs
`default_species=` for backwards compatibility.

## Not yet wired into pipeline

`refine_signal_peptide()` is implemented and tested but **not yet called
from pipeline.py**. The pipeline still uses raw Cys-pair `mature_start`.
Wiring it in requires passing `species_category` through to the SP
detection points in the pipeline.

## Known failures

### Passerine/Strigiform birds (α1 domain ~15aa longer)

Owls (*Strix*, *Athene*, *Asio*), kiwi (*Apteryx*), thrushes (*Catharus*),
finches (*Serinus*), and some other passerines have Cys1 at mature position
110-119 instead of 100. The α1 domain contains a ~15 residue insertion
compared to galliform birds and mammals.

Affected genera (Cys1 position, n):
- Strix: 119, n=16
- Apteryx: 116, n=10
- Catharus: 115, n=5
- Athene: 115, n=3
- Serinus: 113, n=5
- Phasianus: 116, n=3
- Nipponia: 110-116, n=3
- Ciconia: 110-116, n=6

The hybrid scorer's ±5 window can partially compensate, but a 15-residue
offset exceeds the window. These entries get the wrong SP boundary.

**Fix**: Add order-level Cys1 overrides for Strigiformes (~119) and
affected Passeriformes (~115). Requires knowing the bird order, which
needs either species.py taxonomy expansion or an ete3-based lookup.

### Cypriniform fish (variable α1 length)

Carp (*Cyprinus*), zebrafish (*Danio*), and relatives have Cys1 at
92-130 with high variance. Most are near 101 (fine), but outliers at
92 and 127-130 are likely non-classical class I genes.

**Fix**: Detect and exclude non-classical class I (ZAA, LAA lineages)
before SP refinement. Already partially handled by `non_classical`
groove status.

### Class II not evaluated

All SP analysis was on class I alpha. Class II alpha and beta have
different domain architectures and different SP characteristics.
Gene-specific constants (DQA 109, DMA 120, DPB 114) partially
address this but the AxA scan hasn't been evaluated for class II.

## TODO for next PR

### Pipeline integration
- [ ] Call `refine_signal_peptide()` from pipeline.py at all 3 SP
      detection points (IMGT, IPD, diverse)
- [ ] Thread `species_category` to the SP refinement call
- [ ] Verify build counts don't regress

### Evaluation infrastructure
- [ ] `scripts/fetch_sp_ground_truth.py` — fetch UniProt Signal
      features for all entries (SignalP predictions from TrEMBL)
- [ ] `scripts/evaluate_sp.py` — compare our predictions against
      UniProt ground truth, stratified by clade
- [ ] Store UniProt SP length in diverse CSV as `uniprot_sp_length`

### Bird order constants
- [ ] Add Cys1 overrides for Strigiformes (~119) and affected
      Passeriformes (~115)
- [ ] Requires mapping bird species → order in species.py or via
      NCBI Taxonomy / ete3

### Rename groove → domain
- [ ] `extract_groove()` → `decompose_domains()` (keep alias)
- [ ] `parse_class_i()` → `decompose_class_i()`
- [ ] `_try_groove_parse()` → `_try_domain_parse()`
- [ ] Update all docstrings and comments

### Auto-update README
- [ ] Rewrite `scripts/update_readme_counts.py` for new README format
- [ ] Add CI step to auto-commit count updates on push to main

### Class II SP evaluation
- [ ] Run the same analysis on class II alpha and beta
- [ ] Check if DQA/DMA/DPB gene-specific constants + AxA scan work
