# Plan: Archosaur MHC Sequence Gathering and Parse Diagnosis

Generated: 2026-03-16

## Current state

| Group | Alleles | Species | Class I | Class II | Key problems |
|---|---|---|---|---|---|
| Crocodilians | 150 | 13 | 68 (85% fragments) | 73 | Only *C. porosus* has full-length; zero Alligatoridae |
| Turtles | 28 | 8 | 7 | 12 | Mostly fragments; 4 sea turtle genera missing |
| Ratites | 25 | 6 | **0** | 21 | No class I at all; no rhea |

### Root cause of failures
1. **85% of archosaur sequences are single-exon fragments** (80–86 aa) — not a pipeline bug, a data availability problem
2. **100% fail mhcgnomes** — species prefixes not in ontology (except Crpo in 3.3.0)
3. **Groove parsing works fine** when given full-length sequences (53/203 = 26% succeed, all full-length)

## Phase 1: Diagnose existing parse failures (no new data)

### 1.1 Classify all 150 non-ok archosaur groove results

For each failed entry, determine:
- Is it a **fragment** (single exon, <120 aa)? → Expected failure, can't fix without more sequence
- Is it a **full-length sequence with wrong class/chain**? → Fix metadata
- Is it a **species-specific Cys offset** (like DQA/DPB)? → Add constant
- Is it a **non-classical gene**? → Label appropriately

**Status**: The failure analysis above shows 126/150 crocodilian failures are fragments. No algorithmic fix needed — we need more data.

### 1.2 Fix known annotation errors
- `Chse-H2-Q9` (snapping turtle) classified as `class=unknown, chain=unknown` — this is a mammalian H2 gene name transferred by automated annotation. Determine correct class/chain from sequence structure.
- Crocodilian `DA-Ex3` entries — gene name indicates single exon submission. Consider tagging with `is_fragment=True` metadata.

### 1.3 Check for archosaur-specific Cys1 offsets
- Run the same UniProt SP validation we did for DQA/DPB, using the few full-length archosaur class I chains (*C. porosus* UA/UB/UC, *Chelonia mydas* class I).
- If offsets differ from mammalian constants (100 for class I, 106/116 for class II), add archosaur-specific constants.

## Phase 2: Gather new sequences from genome assemblies

### 2.1 Priority genomes with MHC annotations

These genomes have NCBI automated annotation that should include predicted MHC proteins (XP_* accessions):

| Species | Assembly | Group | Priority |
|---|---|---|---|
| *Alligator mississippiensis* | GCF_030867095.1 | crocodilian | **HIGH** — zero Alligatoridae currently |
| *Alligator sinensis* | GCF_000455745.1 | crocodilian | HIGH |
| *Gavialis gangeticus* | GCF_001723915.1 | crocodilian | HIGH — zero Gavialidae |
| *Struthio camelus* | GCF_040807025.1 | ratite | **HIGH** — zero ratite class I |
| *Dromaius novaehollandiae* | GCF_016128335.2 | ratite | HIGH |
| *Chelonia mydas* | GCF_015237465.2 | turtle | MEDIUM — already 8 entries |
| *Caretta caretta* | (available 2023) | turtle | MEDIUM |
| *Dermochelys coriacea* | GCF_009764565.3 | turtle | MEDIUM |
| *Rhea americana* | GCA_003343005.1 | ratite | HIGH — zero rhea |

### 2.2 Script: `scripts/gather_refseq_mhc.py`

Write a script that:
1. Queries NCBI Datasets API for each assembly accession
2. Downloads predicted protein sequences from the MHC genomic region
3. Filters for MHC class I alpha, class II alpha/beta, and B2M by:
   - InterPro domain annotations (IPR001039 = MHC class I α1α2)
   - Gene name keywords (MHC, HLA, class I, class II, B2M, beta-2-microglobulin)
   - Sequence length (class I: 300–400 aa, class II: 200–280 aa)
4. Outputs in diverse_mhc_sequences.csv format
5. Assigns species, gene, class, chain metadata from NCBI annotation

### 2.3 Script: `scripts/gather_genbank_amplicons.py`

For species with published MHC amplicon studies (nucleotide only):
1. Fetch nucleotide sequences from GenBank accessions cited in papers
2. Translate to protein (6-frame, identify the frame with MHC-like features)
3. These will be partial sequences (exon 2 or 3 only, ~54–72 aa protein)
4. Add as fragments with `is_fragment=True`

Key accession ranges:
- **Crocodilian class I**: Jaratlerdsiri et al. 2014 (Immunogenetics 66:175–187) — 124 sequences across 20 species
- **Crocodilian class II**: Jaratlerdsiri et al. 2014 (PLoS ONE 9:e87534) — class II β from 20 species, GenBank AF256650–AF256652, FJ886734–FJ886741, AY491421–AY491430
- **Sea turtle class I**: Martin et al. 2022 (R Soc Open Sci) — 116 alleles, GenBank OK135205–OK135305
- **Sea turtle class I+II**: Leighton et al. 2026 (Genome Biol Evol) — 162 class I + 308 class II alleles, BioProject PRJNA1219623
- **Loggerhead class I**: Stiebens et al. 2013 (BMC Evol Biol) — 34 alleles, GenBank KF021627–KF021666
- **Palaeognath class II**: Minias & Babik 2024 (Genome Biol Evol) — 19 species including ratites

## Phase 3: Improve nomenclature handling

### 3.1 Species prefix policy

Current problem: 4-letter prefixes collide (Gaga = chicken AND gharial, Caca = 5+ species).

**Proposed policy:**
- For species with an **established prefix in the literature** (IPD-MHC, published papers): use that prefix
- For species with **no established prefix**: use the full latin binomial as the primary key; do NOT generate a 4-letter prefix
- For species with **colliding prefixes**: use longer prefixes (e.g., 5+5: `ChrysPickt` for *Chrysemys picta* vs `ChrysPicku` for *Chrysolophus pictus*)

### 3.2 mhcgnomes species additions needed

**Tier 1 (crocodilians):**
- *Alligator mississippiensis*, *A. sinensis*
- *Caiman crocodilus*, *C. latirostris*
- *Gavialis gangeticus* — needs non-colliding prefix (NOT Gaga)
- *Mecistops cataphractus*, *Osteolaemus tetraspis* (already in dataset)

**Tier 2 (turtles):**
- *Chelonia mydas*, *Caretta caretta*, *Dermochelys coriacea*
- *Chrysemys picta*, *Chelydra serpentina*, *Pelodiscus sinensis*
- *Gopherus polyphemus*, *Terrapene triunguis*

**Tier 3 (ratites):**
- *Struthio camelus*, *Dromaius novaehollandiae*, *Rhea americana*
- *Apteryx owenii*, *A. mantelli*, *Casuarius casuarius*

### 3.3 Gene definitions needed

| Gene pattern | Class | Chain | Lineage | Used by |
|---|---|---|---|---|
| UA, UB, UC | I | alpha | U (classical) | Crocodilians, turtles |
| DAA, DAB | II | alpha/beta | — | Crocodilians, turtles |
| DBA, DBB | II | alpha/beta | — | Turtles |
| DB01–DB08 | II | beta | — | Multi-species crocodilian |
| DMA, DMB | II | alpha/beta | Non-classical | Turtles |
| MHCIIB | II | beta | Generic | Ratites |
| MHCIA | I | alpha | Generic | Various |
| UAA, UBA | I | alpha | U | Turtles, fish |

## Phase 4: Validation and integration

### 4.1 Re-run build after new data
- Add new sequences to `diverse_mhc_sequences.csv`
- Run `mhcseqs build`
- Compare groove extraction success rates before/after

### 4.2 Archosaur-specific validation
- Verify Cys-pair positions are conserved (α2 C100–C163, α3 C202–C?)
- Check signal peptide lengths against any UniProt annotations
- Verify class II β2 Cys1 position matches mammalian 116

### 4.3 Expected outcomes
- **Crocodilian** class I: expect 3–10 full-length from alligator/gharial genomes
- **Ratite** class I: expect 5–15 from ostrich/emu/rhea genomes — **this fills a critical gap**
- **Turtle**: expect 10–20 additional from sea turtle genomes
- **Total**: ~30–50 new full-length archosaur MHC proteins

## References

- Jaratlerdsiri et al. 2014, Immunogenetics 66:175–187 (crocodilian class I)
- Jaratlerdsiri et al. 2014, PLoS ONE 9:e87534 (crocodilian class II)
- Jaratlerdsiri et al. 2014, PLoS ONE 9:e114631 (crocodilian MHC structure)
- Martin et al. 2022, R Soc Open Sci 9:211190 (sea turtle class I)
- Leighton et al. 2026, Genome Biol Evol 18:evag008 (4-species sea turtle MHC)
- Stiebens et al. 2013, BMC Evol Biol 13:95 (loggerhead class I)
- Minias & Babik 2024, Genome Biol Evol 16:evae211 (palaeognath class II)
- Balakrishnan et al. 2022, Front Genet 13:823686 (avian MHC architecture)
