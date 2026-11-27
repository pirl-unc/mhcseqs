# Pan-Class Joint MHC-I / MHC-II IMGT Alignment

`mhcseqs.imgt` projects both MHC Class I (`G-ALPHA1` / `G-ALPHA2`) and MHC Class II (`G-ALPHA` / `G-BETA`) peptide-binding grooves onto a single 3D-homologous 202-slot IMGT G-DOMAIN coordinate system (`IMGT_PAN_CLASS_GROOVE_SLOTS`: 100 `HALF1` slots + 102 `HALF2` slots), following [Lefranc et al. (2005)](https://pubmed.ncbi.nlm.nih.gov/15996672/) and verified via Kabsch 3D superposition of PDB `1HHK` (HLA-A\*02:01) against `1DLH` (HLA-DRA\*01:01 / DRB1\*01:01), `1JK8` (HLA-DQA1\*03:01), and `3LQZ` (HLA-DPB1\*02:01).

This page summarizes our quantitative audit across all **36,693** `status="ok"` vertebrate grooves (`26,588` Class I across 233 species; `2,465` Class II α and `7,637` Class II β across 240 species) and evaluates suitability for a joint pan-class MHC-I / MHC-II binding predictor.

## Verdict for Joint MHC-I / MHC-II Predictors

**Yes — `pan_class_fixed_seq` (202 slots) is structurally sound for joint Class I / Class II transfer**, provided a model follows four design rules:

1. **Use the full 202-slot frame or the 46-slot pan-class pocket union (`IMGT_PAN_CLASS_POCKET_SLOTS`) — never a Class I-only 34-position pseudosequence.**
   - NetMHCpan-34 (Class I) and NetMHCIIpan-34 (Class II, Karosiene et al. 2013: 15 DRA α + 19 DRB1 β) share **22 / 34** exact 3D-aligned IMGT positions (`11 HALF1` + `11 HALF2`).
   - A Class I-only 34-position subset drops **12 critical Class II P1/P4/P7/P9 pocket positions** (`DRA α24/31/52/53` → `HALF1:26,33,56,61`; `DRB1 β26/47/70/77/81/85/89/90` → `HALF2:22,45,66,72A,76,80,84,85`).
   - All **46** positions in the exact union (`22` shared + `12` Class I-only + `12` Class II-only) are occupied columns in `IMGT_PAN_CLASS_GROOVE_SLOTS` and exported as `IMGT_PAN_CLASS_POCKET_SLOTS`.
2. **Keep structural gaps (`'-'`) distinct from missing/truncated residues (`'X'`).**
   - `'-'` encodes real 3D geometry that separates Class I from Class II — most prominently `HALF1:57–60` (`PEYW` closed A-pocket cap in Class I vs `----` open P1 groove end in Class II) and `HALF2:89–92` (`LQRT` in Class I vs `----` in Class II β).
3. **Pass `pan_class_insertion_counts` (or lossless `residues` / `a3m_seq`) alongside `pan_class_fixed_seq`.**
   - Across vertebrates, ~3.6–6.4% of grooves carry internal insertions outside the fixed hotspot slots (e.g. Class I `G-ALPHA2` α2-helix kink after `72A`, Class I `G-ALPHA1` after `43`/`81`, Class II `G-BETA` AB-turn `15–19`).
   - `GrooveAlignment.to_pan_class_features()` exports both the 202-char `seq` and the 202-int `insertion_counts` channel so fixed-grid models never have blind spots on loop-insertion alleles, while sequence/transformer models can also consume the lossless `residues` tokens directly.
4. **Condition the peptide decoder on `mhc_class`.**
   - Even with shared 3D MHC pocket coordinates, Class I peptides are length-constrained inside closed A/F pockets whereas Class II peptides extend freely past P1/P9 through open ends.

## 3D Pocket Equivalence (HLA-A\*02:01 vs HLA-DRA\*01:01 / DRB1\*01:01)

| Pan-Class Slot | Class I (HLA-A\*02:01) | Class II (DRA / DRB1\*01:01) | NetMHCpan / NetMHCIIpan | Pocket Role |
|---|---|---|---|---|
| `HALF1:7` | `Y7` | `αQ9` | Both | A-pocket / P1-P4 rim |
| `HALF1:9` | `F9` | `αE11` | Both | B/C-pocket / P6 floor |
| `HALF1:24` | `A24` | `αF22` | Both | B-pocket / P4 floor |
| `HALF1:26,33` | `G26, F33` | `αF24, αI31` | Class II only | P1/P4 floor |
| `HALF1:56` | `G56` | `αA52` | Class II only | P1 wall |
| `HALF1:57–60` | `P57 E58 Y59 W60` | `----` (gap) | Class I (`Y59`) | **Closed Class I A-cap vs open Class II P1 end** |
| `HALF1:61` | `D61` | `αS53` | Class II only | P1 wall |
| `HALF1:66,67,69` | `K66, V67, A69` | `αG58, αA59, αA61` | Both | B/C-pocket / P4-P6 helix |
| `HALF1:73,74,76` | `T73, H74, V76` | `αV65, αD66, αA68` | Both | C/E/F-pocket / P6-P9 helix |
| `HALF1:80,81` | `T80, L81` | `αI72, αM73` | Both | F-pocket / P9 pocket |
| `HALF1:84` | `Y84` | `αR76` | Class I (`Y84`) | Closed Class I F-cap vs Class II P9 H-bond |
| `HALF2:5,7,9` | `V95, R97, Y99` | `βW9, βL11, βF13` | Both | B/C/F floor / P6-P9 floor |
| `HALF2:11`–`74` | `C101`–`C164` | `βC15`–`βC79` | Structural | Conserved G-domain disulfide (99.5% Class I, 98.1% Class II β) |
| `HALF2:22` | `G112` | `βL26` | Class II only | P4 pocket floor |
| `HALF2:24,26` | `H114, Y116` | `βE28, βC30` | Both | D/E/F floor / P4-P7 floor |
| `HALF2:45` | `W133` | `βY47` | Class II only | P7 pocket floor |
| `HALF2:55` | `T143` | `βD57` | Both | F-pocket / P9 salt bridge (T1D-associated β57) |
| `HALF2:59,61A` | `W147, A150` | `βW61, βQ64` | Class I only | E/F helix wall |
| `HALF2:63,67,70,73` | `V152, L156, Y159, T163` | `βL67, βR71, βA74, βY78` | Both | D/E/A helix / P4-P7 shared epitope |
| `HALF2:66,72A,76` | `Q155, G162, E166` | `βQ70, βT77, βH81` | Class II only | Shared epitope β70 & conserved P1 H-bond βH81 |
| `HALF2:80,81,84,85` | `R170, Y171, N174, G175` | `βV85, βG86, βF89, βT90` | Both (`81`) / II (`80,84,85`) | Class I A-pocket / Class II P1 dimorphism (βG86V) |

## Quantitative Audit Across 36,693 Grooves

- **Shared slot occupancy**: `170 / 202` slots (`81/100 HALF1`, `89/102 HALF2`) are actively occupied (≥5%) in both Class I and Class II.
- **Disulfide conservation**: `HALF2:11`–`HALF2:74` aligns `C`–`C` in **99.45%** of Class I (`26,442 / 26,588`) and **98.13%** of Class II β (`7,494 / 7,637`).
- **Monotonic pan-class hotspot bucketing (`v2.25.0`)**: `_build_fixed_frame()` and `_build_pan_class_frame()` strictly anchor every hotspot slot tuple (`7A`, `19A–19C`, `40–41`, `54A–54C` on `HALF1`; `13A–13C`, `40/40A–40C/41`, `61B–61C` on `HALF2`) to its exact left-flanking IMGT column with uniform left-aligned fill across Class I and Class II (and symmetric `cd_turn` fill on `HALF2:40..41`). Off-boundary insertions stay in `overflow_insertions` and `pan_class_insertion_counts` at their exact anchor column instead of transposing across canonical residues, guaranteeing a strict 100% order-preserving subsequence invariant across all grooves.
- **Overhangs outside fixed frames (preserved losslessly in `residues` and `a3m_seq`)**:
  - `G-ALPHA2` C-overhang: mature 183 (`D`/`E`, α2→α3 linker) on **77.9%** of Class I chains (non-pocket).
  - `G-BETA` N-overhang: mature 1–4 (`GDTR`/`RDSP`/`RATP`, upstream of IMGT 1 = Pro5) on **53.2%** of Class II β chains (non-pocket).
  - `G-ALPHA` C-overhang: mature 83–84 (`T`/`NV`, C-terminal linker past R76) on **98.0%** of Class II α chains; N-overhang (`ED`) on **32.9%**.
- **Reserved slot `HALF1:92A`**: retained per Lefranc et al. (2005) (defined solely for `DMA` C-terminal helix extension); remains `'-'` for standard classical grooves.
- **Non-classical `DMB` note**: classical Class II β chains (`DRB`, `DQB`, `DPB`, `DOB`, `H2-Ab1/Eb1`, `RT1-Bb/Db1`, `BLB`, `DAB`) align cleanly with 0–2 deletions. The non-classical peptide-exchange chaperone `DMB` (which does not bind peptides and has a structurally shifted β1 helix in PDB `1HDM`/`4FQX`) diverges sharply from classical β1 profiles and exhibits alignment gaps across the helix.

## Quick API Example

```python
from mhcseqs.imgt import (
    IMGT_PAN_CLASS_POCKET_SLOTS,
    align_groove_to_imgt,
)

# Class I (single chain)
c1 = align_groove_to_imgt(a0201_g1, a0201_g2, mhc_class="I")
f1 = c1.to_pan_class_features()
assert len(f1["seq"]) == 202
assert len(f1["pocket_seq"]) == 46
assert len(f1["insertion_counts"]) == 202

# Class II (heterodimer)
c2 = align_groove_to_imgt(dra_g1, drb1_g2, mhc_class="II")
f2 = c2.to_pan_class_features()
assert len(f2["seq"]) == 202
assert len(f2["pocket_seq"]) == 46
```
