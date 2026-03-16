# Archosaur MHC Sequence Audit

Generated: 2026-03-16

## Summary

| Group | Alleles | Species | Class I | Class II | B2M | Groove ok | mhcgnomes parsed |
|---|---|---|---|---|---|---|---|
| Crocodilians | 150 | 13 (3 genera) | 68 | 73 | 1 | 23 (15%) | Crpo only |
| Turtles | 28 | 8 (7 genera) | 7 | 12 | 8 | 19 (68%) | 0 |
| Ratites | 25 | 6 (4 genera) | **0** | 21 | 4 | 11 (44%) | 0 |
| **Total** | **203** | **27** | **75** | **106** | **13** | **53 (26%)** | **Crpo only** |

## Crocodilians (150 alleles)

### Coverage
Only Crocodylidae represented. **Missing entire families:**
- **Alligatoridae**: *Alligator mississippiensis* (American alligator), *A. sinensis* (Chinese alligator), *Caiman* spp., *Melanosuchus*, *Paleosuchus* — 0 sequences
- **Gavialidae**: *Gavialis gangeticus* (gharial), *Tomistoma schlegelii* (false gharial) — 0 sequences

*Crocodylus porosus* dominates (70/150 entries) and is the only species with full-length class I alpha chains (UA 408aa, UB 293aa, UC 350aa). All other species have only partial fragments (83–86 aa class I, 86–89 aa class II beta).

### Gene nomenclature problems
- Numbered allele naming: `Crpo-Crpo94`, `Crpo-Crpo103`, ... `Crpo-Crpo182` — not parseable by any nomenclature system
- Multi-species class II: `DB01`–`DB08` without species prefix
- Only `Crpo-UA`, `Crpo-UB`, `Crpo-UC`, `Crpo-DAA`, `Crpo-DAB1`, `Crpo-DAB2` follow standard naming

### Data sources to investigate
- Jaratlerdsiri et al. 2014 (crocodilian MHC class I diversity)
- Jaratlerdsiri et al. 2012 (saltwater croc MHC II)
- NCBI GenBank: search "Alligator MHC" — likely has *A. mississippiensis* and *A. sinensis* sequences from genome projects
- NCBI RefSeq: *Alligator mississippiensis* genome (GCF_000281125.4) has MHC annotations

## Turtles (28 alleles)

### Coverage
Reasonable genus diversity (7 genera) but very few sequences per species.

| Species | Entries | Class I | Class II | Notes |
|---|---|---|---|---|
| *Chelonia mydas* | 8 | 4 | 0 | Green sea turtle, genome available |
| *Chrysemys picta* | 4 | 0 | 4 | Painted turtle, DMA/DMB only (non-classical) |
| *Chelydra serpentina* | 4 | 3 | 0 | One entry has H2-Q9 gene name (annotation error) |
| *Pelodiscus sinensis* | 4 | 0 | 4 | Chinese softshell, genome available |
| *Terrapene triunguis* | 3 | 0 | 3 | Box turtle |
| *Gopherus* spp. | 4 | 0 | 4 | Gopher/Goode's tortoise |
| *Chelonoidis abingdonii* | 1 | 0 | 1 | Lonesome George genome |

**Missing:** *Caretta caretta* (loggerhead), *Dermochelys coriacea* (leatherback), *Eretmochelys imbricata* (hawksbill), *Lepidochelys* — all sea turtles with conservation genomics data available.

### Data sources to investigate
- Genome assemblies: *Chelonia mydas* (GCF_015237465.2), *Chrysemys picta* (GCF_000241765.5), *Gopherus evgoodei* (GCF_007399415.2)
- Qin et al. 2022: green sea turtle MHC characterization
- GenBank: "Testudines MHC" — recent genomic surveys exist

## Ratites (25 alleles)

### Coverage
**Critical gap: zero class I alpha chains for any ratite.** Only class II beta (21) and B2M (4).

| Species | Entries | Class I | Class II beta | B2M |
|---|---|---|---|---|
| *Apteryx owenii* | 14 | 0 | 14 | 0 |
| *Casuarius casuarius* | 3 | 0 | 3 | 0 |
| *Dromaius novaehollandiae* | 3 | 0 | 1 | 2 |
| *Struthio camelus* | 2 | 0 | 2 | 0 |
| *Nothoprocta perdicaria* | 2 | 0 | 0 | 2 |
| *Apteryx mantelli* | 1 | 0 | 1 | 0 |

**Missing:** *Rhea americana* (greater rhea), *Pterocnemia pennata* (lesser rhea), *Tinamus major* (great tinamou).

### Data sources to investigate
- Le Breton et al. 2022: kiwi MHC diversity (*A. mantelli*, *A. owenii*)
- Miller et al. 2011: emu MHC class II
- GenBank: "Struthio camelus MHC" — ostrich genome (GCF_000698965.1) has MHC annotations
- Rhea genome (GCF_000956375.2) — likely has unannotated MHC

## Recommendations for mhcgnomes

### Species prefixes to add (priority order)

**Tier 1 — crocodilians (6 genera):**
| Prefix | Species | Notes |
|---|---|---|
| Almi | *Alligator mississippiensis* | American alligator, genome available |
| Alsi | *Alligator sinensis* | Chinese alligator, genome available |
| Crcr | *Caiman crocodilus* | Spectacled caiman |
| Cala | *Caiman latirostris* | Broad-snouted caiman |
| Gaga | *Gavialis gangeticus* | Gharial — COLLISION with chicken (Gallus gallus) |
| Tosc | *Tomistoma schlegelii* | False gharial |
| Meca | *Mecistops cataphractus* | Slender-snouted croc (already in dataset) |
| Oste | *Osteolaemus tetraspis* | Dwarf croc (already in dataset) |

**Note:** `Gaga` is already taken by *Gallus gallus* (chicken). The gharial would need an alternative prefix — possibly `Ggan` (3+1 split) or `Gavi` (genus-only).

**Tier 2 — turtles (10+ genera):**
| Prefix | Species |
|---|---|
| Chmy | *Chelonia mydas* |
| Chpi | *Chrysemys picta* — COLLISION with *Chrysolophus pictus* (golden pheasant) |
| Chse | *Chelydra serpentina* |
| Pesi | *Pelodiscus sinensis* |
| Tetr | *Terrapene triunguis* |
| Gopo | *Gopherus polyphemus* |
| Caca | *Caretta caretta* — COLLISION with 4 other species |
| Deco | *Dermochelys coriacea* |

**Tier 3 — ratites (5+ genera):**
| Prefix | Species |
|---|---|
| Stca | *Struthio camelus* |
| Drno | *Dromaius novaehollandiae* |
| Caca | *Casuarius casuarius* — COLLISION |
| Apow | *Apteryx owenii* |
| Apma | *Apteryx mantelli* |
| Rham | *Rhea americana* |
| Nope | *Nothoprocta perdicaria* |

### Gene definitions needed
- Reptilian class I: UA, UB, UC (U-lineage, like fish)
- Reptilian class II: DAA, DAB, DBA, DBB (fish-style) plus DMA, DMB
- Generic bird: MHCIIB, MHCIA (non-standard, should map to class II beta / class I alpha)
- B2M for all groups

### Prefix collision policy
The 4-letter prefix space has 26^4 = 456,976 possible codes but species names cluster heavily in the same genus/species initial letters. Known archosaur collisions:
- `Gaga`: *Gavialis gangeticus* vs *Gallus gallus* (chicken)
- `Caca`: *Casuarius casuarius* vs *Caretta caretta* vs *Carassius carassius* vs *Carduelis carduelis* vs *Calidris canutus*
- `Chpi`: *Chrysemys picta* vs *Chrysolophus pictus*

**Recommendation:** mhcgnomes should support latin binomial names as primary keys, with 4-letter prefixes as aliases that resolve unambiguously within the ontology. When a prefix is ambiguous, the full binomial must be used.
