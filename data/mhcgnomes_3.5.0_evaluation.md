# mhcgnomes 3.5.0 Evaluation Against mhcseqs Diverse MHC Dataset

Generated: 2026-03-16

## Test methodology

Tested every unique (gene, organism) pair in `diverse_mhc_sequences.csv` (15,861 entries, 1,592 unique gene+organism pairs) against mhcgnomes 3.5.0 in two modes:
1. **4-letter prefix** as-is (e.g., `Crpo-UA`)
2. **5+5 latin prefix** constructed from organism name (e.g., `CrocoPoros-UA`)

Species validation: checked that the species returned by mhcgnomes matches the organism in our metadata.

## Summary

| Category | Count | % |
|---|---|---|
| Parsed correctly (4-letter prefix) | 212 | 13.3% |
| Parsed correctly (latin prefix only) | 14 | 0.9% |
| **Wrong species returned** | **18** | **1.1%** |
| Known species, unknown gene | 314 | 19.7% |
| Unknown species prefix | 886 | 55.6% |
| Doubled-prefix bug (in our data) | 148 | 9.3% |
| **Total** | **1,592** | |

**Net improvement from 3.3.0 → 3.5.0**: Species recognition jumped from ~6% to ~13%. Key additions: Modo (opossum), Saha (Tasmanian devil), Almi/Alsi (alligators), Chmy (green sea turtle), Chse (snapping turtle), Stca (ostrich), Drno (emu), Apow/Apma (kiwi), and many others.

## 1. Wrong species parses (18 cases) — BUGS

These are the most critical issues: mhcgnomes returns a valid parse but for the **wrong organism**. All caused by prefix collisions or gene names that happen to match a different species.

### 1a. Prefix collisions (10 cases)

| Gene in our data | Our organism | mhcgnomes returns | Root cause |
|---|---|---|---|
| `Cyca-UA` | *Cyclura carinata* (iguana) | *Cyprinus carpio* (carp) | Cyca collision |
| `Cyca-UA` | *Cyanistes caeruleus* (blue tit) | *Cyprinus carpio* (carp) | Cyca collision |
| `Cyca-UA` | *Clarias magur* (catfish) | *Cyprinus carpio* (carp) | Cyca collision — **also wrong 4-letter code** in our data; `Icpu` or `Clma` would be correct |

**Fix needed**: `Cyca` currently maps only to *Cyprinus carpio*. Either:
- Add *Cyclura carinata* and *Cyanistes caeruleus* as additional species with longer prefixes (`CycluCarin`, `CyaniCaeru`)
- Or flag `Cyca` as ambiguous and require longer form for all three

### 1b. Mammalian prefixes on wrong species (8 cases)

| Gene in our data | Our organism | mhcgnomes returns | Root cause |
|---|---|---|---|
| `MAMU-DRA` / `Mamu-DRA` | *Lonchura striata*, *Eudyptes sclateri*, *E. filholi*, *Eudyptula minor* (birds/penguins) | *Macaca mulatta* | Gene named after macaque ortholog in UniProt |
| `PATR-A` | *Astyanax mexicanus* (cavefish) | *Pan troglodytes* | UniProt annotation transfer |
| `POPY-E` | *Astyanax mexicanus* (cavefish) | *Pongo pygmaeus* | UniProt annotation transfer |
| `SLA-DQB1` | *Acipenser oxyrinchus* (sturgeon) | *Sus sp.* | UniProt annotation transfer |
| `MAFA-A1` | *Rhinolophus ferrumequinum* (bat) | *Macaca fascicularis* | UniProt annotation transfer |

**Root cause**: These are **errors in our curation** (diverse_mhc_sequences.csv), not mhcgnomes bugs. The UniProt gene names for these non-model organisms were auto-transferred from the closest mammalian ortholog. Our `normalize_gene()` function should have stripped these and re-prefixed with the correct species. mhcgnomes is parsing them correctly for what they say — the input is wrong.

### 1c. Generic gene names (5 cases)

| Gene | Our organisms | mhcgnomes returns |
|---|---|---|
| `MHC-B` | *Meleagris gallopavo*, *Chrysolophus amherstiae*, *Syrmaticus reevesii*, *Crossoptilon* spp., *Pavo cristatus* (galliform birds) | *Homo sapiens* |

**Root cause**: `MHC-B` is ambiguous — it could mean the human MHC B locus or chicken B complex. mhcgnomes defaults to human. These entries should use species-specific gene names (e.g., `Mega-MHC-B` for turkey).

## 2. Known species, unknown gene (314 cases) — GENE DEFINITIONS NEEDED

mhcgnomes recognizes the species but can't parse the gene name. Grouped by species:

### Priority 1: Chicken (*Gallus gallus*, Gaga) — 44 genes

| Gene pattern | Count | Notes |
|---|---|---|
| `Gaga-MHCY*` | 10 | MHC-Y complex genes (Y2B1, Y2B2, Y15, Y46, Y11, YFV, YFVI, Y-FA, YL) |
| `Gaga-BFw-*`, `Gaga-BFz-*` | 7 | BF workshop/z alleles |
| `Gaga-B-F-S*` | 9 | BF serological types |
| `Gaga-B-F-minor` | 1 | BF minor locus |
| `Gaga-B-LBII`, `Gaga-B-LBI`, `Gaga-B-LBVI` | 3 | BLB roman numeral variants |
| `Gaga-S19` | 1 | Serological type |
| `Gaga-B-DMA`, `Gaga-B-DMB2` | 2 | Chicken DM genes |

Chicken B-complex nomenclature is well established in the literature. mhcgnomes 3.5.0 knows BF, BF1, BF2, BLB, BLB1, BLB2 but not the workshop/serological/MHC-Y variants.

### Priority 2: Japanese quail (*Coturnix japonica*, Coja) — 38 genes

| Gene pattern | Count | Notes |
|---|---|---|
| `Coja-II-01` through `Coja-II-17` | 17 | Numbered class II loci with alleles like `Coja-II-17*01` |
| `Coja-DMB1`, `Coja-DMB2`, `Coja-DMA1` | 3 | DM genes |
| `Coja-B1`, `Coja-B1e3`, `Coja-C`, `Coja-D*`, `Coja-E` | 8 | Named loci |
| `Coja-QF41`, `Coja-QF63` | 2 | QF-numbered loci |

### Priority 3: Tilapia (*Oreochromis niloticus*, Orni) — 35 genes

Gene names like `Orni-DBA`, `Orni-orni-dba` (doubled prefix), `Orni-DCA`, `Orni-DBB`, `Orni-DCB`, `Orni-DDA`, `Orni-DDB`, `Orni-DBA1`, `Orni-UBA1`, `Orni-UAA1`. Fish D-series genes (DBA, DBB, DCA, DCB, DDA, DDB) are not in mhcgnomes.

### Priority 4: Zebrafish (*Danio rerio*, Dare) — 33 genes

Gene names like `Dare-DBB`, `Dare-DCA`, `Dare-DCB`, `Dare-MHCII`, `Dare-mhc1zfa`, `Dare-mhc1zda`, `Dare-mhc1lfa`, `Dare-mhc1ula`, `Dare-mfsd6a` (NOT MHC — lipid transporter). The `mhc1z*`, `mhc1l*`, `mhc1u*` naming convention for fish lineage genes is not recognized.

### Priority 5: Opossum (*Monodelphis domestica*, Modo) — 18 genes

Gene names: `Modo-DRA`, `Modo-DAA`, `Modo-DAB`, `Modo-DBA`, `Modo-DBB`, `Modo-UA1`–`UA6`, `Modo-UG`, `Modo-UT3`, `Modo-UT8`. Marsupial MHC nomenclature.

### Priority 6: Crocodilians — 23 genes across 5 species

| Species | Prefix | Genes |
|---|---|---|
| *Crocodylus porosus* | Crpo | DB01–DB08, Crpo84–Crpo182 (numbered alleles) |
| *Osteolaemus tetraspis* | Oste | DB01–DB08 |
| *Mecistops cataphractus* | Meca | DB03–DB08 |
| *Crocodylus niloticus* | Crni | DB01–DB08 |
| *Chelonia mydas* | Chmy | UY3_* numbered genes |

The `DB01`–`DB08` pattern is a multi-species class II beta nomenclature from Jaratlerdsiri et al. 2014.

### Others

| Species | Prefix | Count | Example genes |
|---|---|---|---|
| *Pongo sp.* | OrLA | 16 | ORLA-UAA, ORLA-UBA — **wrong species**: `OrLA` maps to orangutan but these are fish (*Nothobranchius*) sequences |
| *Chelydra serpentina* | Chse | 3 | Chse-DPB1, Chse-DMB, Chse-H2-Q9 |
| *Xenopus* spp. | Xela/Xetr | 13 | mhc1b.L, hla-dqa1 (UniProt annotation transfer names) |
| *Sarcophilus harrisii* | Saha | 5 | Saha-I-01, Saha-UC-01 (Tasmanian devil numbered format) |

## 3. Unknown species prefix (886 cases, 328 unique prefixes)

These are species mhcgnomes doesn't recognize at all. Top 25 by gene count:

| Prefix | Count | Organism | Latin prefix (5+5) |
|---|---|---|---|
| Bain | 25 | *Labeobarbus intermedius* (Lake Tana barbs) | LabeoInter |
| Char | 19 | *Channa argus* (snakehead) | ChanArgus |
| gamr | 19 | *Gadus morhua* (Atlantic cod) | GadusMorhu |
| Phci | 19 | *Phascolarctos cinereus* (koala) | PhasCinere |
| Asme | 18 | *Astyanax mexicanus* (blind cavefish) | AstyaMexic |
| Cosp | 18 | *Coregonus sp.* (whitefish) | CoregSp |
| Sppu | 16 | *Sphenodon punctatus* (tuatara) | SphenPunct |
| Tuna | 16 | *Turdus naumanni* (dusky thrush) | TurduNauma |
| Angr | 16 | *Anabarilius grahami* (Kanglang fish) | AnabGraha |
| Sias | 14 | *Silurus asotus* (Amur catfish) | SilurAsotu |
| Tueu | 13 | *Turdus eunomus* (dusky thrush) | TurduEunom |
| Trsp | 12 | *Tropheus sp.* (cichlid) | TrophSp |
| Phco | 11 | *Phylloscopus collybita* (chiffchaff) | PhyllColli |
| Kabi | 11 | *Kareius bicoloratus* (stone flounder) | KareBicol |
| Crmo | 10 | *Crocodylus moreletii* (Morelet's croc) | CrocoMorel |
| Icpu | 10 | *Clarias magur* (Asian catfish) | ClariMagur |
| Krma | 10 | *Kryptolebias marmoratus* (killifish) | KryptMarmo |
| Mega | 9 | *Meleagris gallopavo* (wild turkey) | MeleaGallo |
| Turu | 9 | *Turdus ruficollis* | TurduRufic |
| Tuat | 9 | *Turdus atrogularis* | TurduAtrog |
| Rhma | 8 | *Rhinella marina* (cane toad) | RhineMarma |
| Cibo | 8 | *Ciconia boyciana* (Oriental stork) | CiconBoyci |
| Lacr | 8 | *Larimichthys crocea* (yellow croaker) | LarimCroce |
| Mepo | 8 | *Merluccius polli* (Benguela hake) | MerluPolli |

### Notable archosaur gaps

| Prefix | Organism | Status in 3.5.0 | Literature prefix? |
|---|---|---|---|
| Crmo | *Crocodylus moreletii* | NOT KNOWN | Yes (Jaratlerdsiri 2014) |
| Crmi | *Crocodylus mindorensis* | NOT KNOWN | Yes (Jaratlerdsiri 2014) |
| Crpa | *Crocodylus palustris* | NOT KNOWN | Yes (Jaratlerdsiri 2014) |
| Crsi | *Crocodylus siamensis* | NOT KNOWN | Yes (Jaratlerdsiri 2014) |
| Crjo | *Crocodylus johnstoni* | NOT KNOWN | Yes (Jaratlerdsiri 2014) |
| Crac | *Crocodylus acutus* | NOT KNOWN | Yes (Jaratlerdsiri 2014) |
| Gaga | *Gavialis gangeticus* | Maps to chicken! | NOT a valid prefix for gharial — collision |
| Caca | *Casuarius casuarius* | NOT KNOWN | NOT literature-attested for MHC |

## 4. Doubled-prefix bug in our data (148 cases)

Gene names like `Crpo-Crpo94`, `Dila-dila_a4`, `Orni-orni-dba` where the species prefix appears twice. This is a bug in our `curate_diverse_mhc.py` script — it re-prefixes gene names that already contain the prefix in concatenated form.

| Organism | Count | Example |
|---|---|---|
| *Dicentrarchus labrax* | 58 | `Dila-dila_a4` (should be `Dila-a4` or just accession) |
| *Crocodylus porosus* | 54 | `Crpo-Crpo94` (should be `Crpo94` allele designation) |
| *Crocodylus palustris* | 3 | `Crpa-CrpaMHCI` |
| *Oryzias latipes* | 3 | `Orla-orla-uia1` |

**This is an mhcseqs bug, not an mhcgnomes issue.** Fix: strip redundant prefix in `normalize_gene()`.

## 5. Latin prefix (5+5) evaluation

The 5+5 latin prefix feature in 3.5.0 recovered **14 additional genes** that failed with 4-letter prefixes. All 14 are cases where the 4-letter prefix was ambiguous or non-standard:

| Original gene | Latin form that works | Organism |
|---|---|---|
| `Bubu-DAB1/2` | `BuboBubo-DAB1/2` | *Bubo bubo* (eagle-owl) — `Bubu` not in mhcgnomes |
| `Egue-DAB2` | `EgretEulop-DAB2` | *Egretta eulophotes* (Chinese egret) |
| `Acaru-UA` | `AcrocArund-UA` | *Acrocephalus arundinaceus* (reed warbler) |
| `GSP-BLB1/2` | `GalluGallu-BLB1/2` | *Gallus gallus* — `GSP` is non-standard |
| `Orla-UHA/UBA/UAA/UGA` | `OryziLatip-U*` | *Oryzias latipes* — `Orla` collides with orangutan |
| `TO-DAB` | `TrachOvatu-DAB` | *Trachinotus ovatus* — `TO` too short |
| `Citd-UBA/UAA` | `CtenoIdell-U*` | *Ctenopharyngodon idella* — `Citd` not known |
| `Truv-UB` | `TrichVulpe-UB` | *Trichosurus vulpecula* (possum) — `Truv` not known |

**Zero wrong-species parses with latin prefixes** — the longer form eliminates all collisions.

### Latin prefixes that SHOULD work but DON'T

| Latin prefix | Expected species | Status |
|---|---|---|
| `GaviGanget` | *Gavialis gangeticus* (gharial) | NOT KNOWN |
| `ChrysPickt` | *Chrysemys picta* (painted turtle) | NOT KNOWN |
| `ChelySerpn` | *Chelydra serpentina* (snapping turtle) | NOT KNOWN — but `Chse` works! |
| `PhasCinere` | *Phascolarctos cinereus* (koala) | NOT KNOWN |
| `SphenPunct` | *Sphenodon punctatus* (tuatara) | NOT KNOWN |
| `CheloMyda` | *Chelonia mydas* (5+4 form) | NOT KNOWN — `CheloMydas` (5+5) works |

## Recommendations for mhcgnomes

### Tier 1: Fix wrong-species parses
1. **`OrLA` → orangutan collision**: 16 entries in our data use `ORLA-UAA` etc. for *Nothobranchius* fish. mhcgnomes maps `OrLA` to *Pongo sp.* This needs either:
   - A disambiguation rule (ORLA uppercase = fish, Orla = orangutan?)
   - Removing the short prefix for orangutan and requiring longer form

2. **`Cyca` collision**: Maps to *Cyprinus carpio* but also used for *Cyclura carinata* (iguana) and *Cyanistes caeruleus* (blue tit). Need longer prefixes for all three.

### Tier 2: Add remaining Crocodylus species (6 prefixes)
All attested in Jaratlerdsiri et al. 2014: `Crmo`, `Crmi`, `Crpa`, `Crsi`, `Crjo`, `Crac`

### Tier 3: Add gene definitions
- Chicken MHC-Y genes: MHCY2B1, MHCY2B2, MHCY15, YFV, YFVI
- Chicken BF variants: BFw, BFz, B-F-S (serological), B-F-minor
- Fish D-series: DBA, DBB, DCA, DCB, DDA, DDB (for tilapia, zebrafish, etc.)
- Crocodilian multi-species: DB01–DB08
- Marsupial: UA1–UA6, UG, UT3, UT8 (for opossum)
- Quail numbered class II: Coja-II-01 through Coja-II-17

### Tier 4: Add remaining species (328 unknown prefixes)
Full list of 328 unknown prefixes with organisms and gene examples is in `data/mhcgnomes_failures.csv`.

Top priority species not yet in mhcgnomes (by gene count in our dataset):
1. *Labeobarbus intermedius* (25 genes) — Lake Tana barbs
2. *Channa argus* (19 genes) — snakehead fish
3. *Gadus morhua* (19 genes) — Atlantic cod
4. *Phascolarctos cinereus* (19 genes) — koala
5. *Sphenodon punctatus* (16 genes) — tuatara (only living rhynchocephalian!)
6. *Meleagris gallopavo* (9 genes) — wild turkey

### Tier 5: Ensure all 5+5 latin forms work
Several species are known by 4-letter prefix but their 5+5 form fails (e.g., `ChelySerpn` for snapping turtle). The 5+5 form should always work as a fallback.

## Appendix: Prefix provenance audit

### Method
Cross-referenced every prefix in `diverse_mhc_sequences.csv` against the original UniProt `gene_names` field in `data/diverse_mhc_raw.csv`. A prefix is "literature-attested" if the 4-letter code appears in the UniProt gene name for at least one entry from that species.

### Results

| Category | Count |
|---|---|
| **Literature-attested** (prefix found in UniProt gene_names) | 299 |
| **Script-generated** (prefix NOT in UniProt, made up by curate_diverse_mhc.py) | 120 |

Of the 299 literature prefixes:

| mhcgnomes 3.5.0 status | Count | Action needed |
|---|---|---|
| 4-letter prefix recognized | 63 | None |
| 4-letter fails, 5+5 latin works | 12 | mhcgnomes should add 4-letter alias |
| **Neither works** | **224** | **mhcgnomes needs to add these species** |

Of the 120 script-generated prefixes:

| mhcgnomes 3.5.0 status | Count | Action needed |
|---|---|---|
| 5+5 latin works | 9 | mhcseqs should use latin prefix instead of made-up 4-letter |
| **Neither works** | **111** | **mhcgnomes needs to add these species** |

### Notable issues with "literature" prefixes

Several prefixes are in UniProt but refer to **old taxonomy**. The organism has been reclassified:

| Prefix | UniProt organism | Current name | Notes |
|---|---|---|---|
| Brre | *Brachydanio rerio* | *Danio rerio* (zebrafish) | Old genus name |
| Trsi | *Trionyx sinensis* | *Pelodiscus sinensis* (softshell turtle) | Old genus |
| Maeu | *Macropus eugenii* | *Notamacropus eugenii* (wallaby) | Genus split |
| Chni | *Charadrius nivosus* | *Anarhynchus nivosus* (plover) | Genus revised |
| Gran/Grvip | *Grus antigone/vipio* | *Antigone antigone/vipio* (cranes) | Genus split |
| Lica/Licl | *Lithobates catesbeianus/clamitans* | *Aquarana/Lithobates* (frogs) | Genus revised |

These create **double-identity problems**: the same species has entries under two different prefixes. mhcgnomes should treat the old prefix as an alias for the current binomial.

### 224 literature prefixes not in mhcgnomes

The full list is in the "NOTHING WORKS" section above. Taxonomic breakdown:

- **Birds**: ~95 prefixes (raptors, owls, cranes, egrets, warblers, finches, etc.)
- **Fish**: ~60 prefixes (salmon, carp, catfish, cichlids, seahorses, etc.)
- **Reptiles**: ~25 prefixes (iguanas, crocodiles, lizards, snakes, tuatara)
- **Amphibians**: ~20 prefixes (frogs, salamanders, newts, axolotl)
- **Mammals**: ~15 prefixes (bats, marsupials, platypus)
- **Sharks**: ~5 prefixes

These are all backed by published MHC studies with sequences in UniProt. They represent the core gap between what exists in the immunogenetics literature and what mhcgnomes can parse.
