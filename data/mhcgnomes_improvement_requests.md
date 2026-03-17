# mhcgnomes Improvement Requests from mhcseqs v0.9.0

Generated: 2026-03-16
Tested against: mhcgnomes 3.5.0
Dataset: 15,861 diverse MHC entries (1,526 unique gene names)

## Parsing summary

| Result | Count | % |
|---|---|---|
| Parsed correctly | 200 | 13.1% |
| Wrong species returned | 8 | 0.5% |
| Known species, unknown gene | 320 | 21.0% |
| Unknown species | 894 | 58.6% |
| Doubled prefix (our bug) | 148 | 9.7% |

## 1. Wrong-species parses (8 unique, affecting 18 entries)

These are silent data corruption — mhcgnomes returns a confident parse for the wrong organism.

### 1a. `Cyca` collision (3 entries)

```
Cyca-UA  expected: Cyclura carinata (iguana)        got: Cyprinus carpio (carp)
Cyca-UA  expected: Cyanistes caeruleus (blue tit)   got: Cyprinus carpio (carp)
Cyca-UA  expected: Clarias magur (catfish)           got: Cyprinus carpio (carp)
```

**Advice**: `Cyca` is ambiguous across 3+ genera. Either require the long form for all (`CypriCarpi`, `CycluCarin`, `CyaniCaeru`) or restrict `Cyca` to carp and add the others under longer prefixes. The iguana and blue tit gene names come from published UniProt entries, not our fabrication.

### 1b. `MHC-B` defaults to human (6 entries)

```
MHC-B  expected: Meleagris gallopavo (turkey)            got: Homo sapiens
MHC-B  expected: Chrysolophus amherstiae (Lady Amherst)  got: Homo sapiens
MHC-B  expected: Pavo cristatus (Indian peafowl)         got: Homo sapiens
MHC-B  expected: Crossoptilon crossoptilon (pheasant)    got: Homo sapiens
MHC-B  expected: Crossoptilon auritum (pheasant)         got: Homo sapiens
MHC-B  expected: Syrmaticus reevesii (pheasant)          got: Homo sapiens
```

**Advice**: Bare `MHC-B` without a species prefix is fundamentally ambiguous. These are galliform birds using the chicken B-complex naming convention. Consider: (a) not parsing `MHC-B` as human since it's too generic, or (b) only parsing it as human when the context is explicitly human (e.g., `HLA-B` but not `MHC-B`).

### 1c. Mammalian prefixes on wrong species (9 entries)

These are UniProt annotation-transfer artifacts where NCBI named a fish/bird gene after its mammalian ortholog:

```
MAMU-DRA   expected: penguins/munia (5 species)     got: Macaca mulatta
PATR-A     expected: Astyanax mexicanus (cavefish)   got: Pan troglodytes
POPY-E     expected: Astyanax mexicanus (cavefish)   got: Pongo pygmaeus
SLA-DQB1   expected: Acipenser oxyrinchus (sturgeon) got: Sus sp.
MAFA-A1    expected: Rhinolophus ferrumequinum (bat)  got: Macaca fascicularis
```

**Advice**: This is partly our data quality issue (these gene names are wrong in our CSV). But mhcgnomes could help by warning when a well-known mammalian prefix (MAMU, PATR, POPY, SLA, MAFA) appears in an unusual context. No action required from mhcgnomes — we'll fix our curation.

## 2. Known species, unknown gene (320 entries)

mhcgnomes recognizes the species but can't parse the gene name. These are the highest-value additions since no new species definitions are needed.

### Chicken (*Gallus gallus*, Gaga) — 44 genes

mhcgnomes knows BF, BF1, BF2, BLB, BLB1, BLB2 but not:

| Gene pattern | Count | What it is | Example |
|---|---|---|---|
| MHCY2B1, MHCY2B2, MHCY15, MHCY46, MHCY11 | 8 | MHC-Y complex (second MHC in chicken) | `Gaga-MHCY2B2` |
| YFV, YFVI, Y-FA, YL | 5 | Y-F locus class I genes | `Gaga-YFV` |
| BFw-01 through BFw-06 | 6 | BF workshop alleles | `Gaga-BFw-03` |
| BFz-01 through BFz-05 | 5 | BF z-series alleles | `Gaga-BFz-02` |
| B-F-S01 through B-F-S14 | 9 | BF serological types | `Gaga-B-F-S05` |
| B-DMA, B-DMB2 | 2 | Chicken DM genes | `Gaga-B-DMA` |
| B-LBI, B-LBII, B-LBVI | 3 | BLB roman numeral variants | `Gaga-B-LBII` |
| B-F-minor, S19 | 2 | Minor locus, serological | `Gaga-B-F-minor` |

**Advice**: The MHC-Y complex is well-characterized (Kaufman 2000, Hosomichi 2006). Gene definitions for MHCY* and YF* would cover 13 genes. BF workshop/serological types could be handled by a pattern rule: `BFw-NN`, `BFz-NN`, `B-F-SNN`.

### Japanese quail (*Coturnix japonica*, Coja) — 38 genes

| Gene pattern | Count | Example |
|---|---|---|
| II-01 through II-17 (numbered class II) | 17 | `Coja-II-17*01` |
| DMB1, DMB2, DMA1 | 3 | `Coja-DMB1` |
| B1, B1e3, C, D, D1, D2, E | 8 | `Coja-B1` |
| QF41, QF63 | 2 | `Coja-QF41` |

**Advice**: The `Coja-II-NN` format is unique to quail MHC nomenclature (Hosomichi et al. 2006). A pattern rule for `II-\d+` would cover 17 genes. DM genes and single-letter loci (B, C, D, E) are standard.

### Tilapia (*Oreochromis niloticus*, Orni) — 35 genes

Gene names: DBA, DBB, DCA, DCB, DDA, DDB, UAA1, UBA1, UCA, UDA, UEA, UFA, UGA, UHA, UIA.

**Advice**: These are standard fish MHC gene names following the D-series (class II) and U-lineage (class I) conventions. They should be added as gene definitions for all fish species, not just tilapia. The numbered U-lineage loci (UAA through UIA) represent distinct genomic loci characterized by Sato et al. (2006) and Dijkstra et al. (2013).

### Zebrafish (*Danio rerio*, Dare) — 33 genes

Gene names: mhc1zfa, mhc1zda, mhc1laa, mhc1lba, mhc1ula, mhc2dga, mhc2dca, DBB, DCA, DCB.

**Advice**: The `mhc1z*`, `mhc1l*`, `mhc1u*` naming convention is the official zebrafish nomenclature committee system (Bingulac-Popovic et al. 1997, Dirscherl et al. 2014). The letter after `mhc1` indicates the lineage (z=Z, l=L, u=U), and the last two letters are the locus name. These are widely used in the zebrafish community.

### Opossum (*Monodelphis domestica*, Modo) — 18 genes

Gene names: DRA, DAA, DAB, DBA, DBB, UA1-UA6, UG, UT3, UT8, UC, UE, UI, MIC-like, MILL-like.

**Advice**: Marsupial MHC nomenclature from Belov et al. (2006) and Siddle et al. (2009). The UA/UG/UT genes are U-lineage class I, numbered by locus. DRA/DAA/DAB/DBA/DBB are class II using the fish-like D-series convention (marsupials retained the ancestral arrangement).

### Crocodilians (26 genes across Crpo, Oste, Meca, Crni)

Gene names: DB01 through DB08.

**Advice**: Multi-species class II beta nomenclature from Jaratlerdsiri et al. (2014). `DB01`-`DB08` are trans-specific class II beta allele lineages across Crocodylia. They should be defined as class II beta gene patterns for all crocodilian species.

### Other notable gaps

| Species | Prefix | Count | Example genes |
|---|---|---|---|
| *Chelonia mydas* (sea turtle) | Chmy | 8 | UY3_17009 (genome scaffold IDs) |
| *Anas platyrhynchos* (duck) | Anpl | 6 | DRA, U, UEA, UBA |
| *Sarcophilus harrisii* (Tasm. devil) | Saha | 5 | I-01, UC-01 |
| *Chelydra serpentina* (snapping turtle) | Chse | 3 | DPB1, DMB |
| *Trichosurus vulpecula* (possum) | Trvu | 4 | DAA, DRA, DRB2, DRB3 |

## 3. Unknown species (894 entries, 328 unique prefixes)

These species are not in the mhcgnomes ontology at all. Top 25 by impact:

| Prefix | Count | Organism | Literature? | Latin 5+5 |
|---|---|---|---|---|
| Bain | 25 | *Labeobarbus intermedius* | Yes | LabeoInter |
| Char | 19 | *Channa argus* (snakehead) | Yes | ChanArgus |
| gamr | 19 | *Gadus morhua* (Atlantic cod) | Yes | GadusMorhu |
| Phci | 19 | *Phascolarctos cinereus* (koala) | Yes | PhasCinere |
| Asme | 18 | *Astyanax mexicanus* (cavefish) | Yes | AstyaMexic |
| Sppu | 16 | *Sphenodon punctatus* (tuatara) | Yes | SphenPunct |
| Tuna | 16 | *Turdus naumanni* (thrush) | Yes | TurduNauma |
| Angr | 16 | *Anabarilius grahami* | Yes | AnabGraha |
| Sias | 14 | *Silurus asotus* (catfish) | Yes | SilurAsotu |
| Crmo | 10 | *Crocodylus moreletii* | Yes | CrocoMorel |
| Icpu | 10 | *Clarias magur* (catfish) | Yes | ClariMagur |
| Krma | 10 | *Kryptolebias marmoratus* | Yes | KryptMarmo |
| Mega | 9 | *Meleagris gallopavo* (turkey) | Yes | MeleaGallo |
| Rhma | 8 | *Rhinella marina* (cane toad) | Yes | RhineMarma |
| Cibo | 8 | *Ciconia boyciana* (stork) | Yes | CiconBoyci |
| Dila | 62 | *Dicentrarchus labrax* (seabass) | Yes | DicenLabra |
| Epco | 12 | *Epinephelus coioides* (grouper) | Yes | EpineCoioi |
| Hiab | 7 | *Hippocampus abdominalis* (seahorse) | Yes | HippoAbdom |
| Ocle | 9 | *Oceanodroma leucorhoa* (petrel) | Yes | OceanLeuco |
| Napa | 8 | *Nanorana parkeri* (frog) | Yes | NanorParke |
| Haal | 5 | *Haliaeetus albicilla* (eagle) | Yes | HaliaAlbic |
| Amcr | 3 | *Amblyrhynchus cristatus* (iguana) | Yes | AmblyCrist |
| Oran | 3 | *Ornithorhynchus anatinus* (platypus) | Yes | OrnitAnati |
| Sppu | 16 | *Sphenodon punctatus* (tuatara) | Yes | SphenPunct |
| Taac | 2 | *Tachyglossus aculeatus* (echidna) | Yes | TachyAcule |

All 328 prefixes are listed with organisms in the full evaluation at `data/mhcgnomes_3.5.0_evaluation.md`.

**Advice**: These are all backed by published MHC studies with sequences in UniProt. The 4-letter prefixes listed above are literature-attested (they appear in UniProt gene_names fields). Adding the latin 5+5 form for each would be sufficient — mhcseqs can use `GadusMorhu-UA` instead of `gamr-UA` and it will work.

**Priority species for archosaur research**: *Sphenodon punctatus* (tuatara, only living rhynchocephalian), *Crocodylus moreletii* and other Crocodylus spp. (6 species with Jaratlerdsiri 2014 literature prefixes), *Meleagris gallopavo* (turkey, economically important).

## 4. Taxonomy reclassification aliases

Several species have old genus names in UniProt that create double-identity problems:

| Old prefix | Current binomial | Old genus | Entries |
|---|---|---|---|
| Brre | *Danio rerio* (zebrafish) | *Brachydanio* | 1 |
| Trsi | *Pelodiscus sinensis* (softshell turtle) | *Trionyx* | 2 |
| Maeu | *Notamacropus eugenii* (wallaby) | *Macropus* | 30 |
| Gran/Grvip | *Antigone* spp. (cranes) | *Grus* | 8 |
| Chni | *Anarhynchus nivosus* (plover) | *Charadrius* | 46 |
| Lica/Licl | *Aquarana/Lithobates* (frogs) | *Rana* | 28 |

**Advice**: Old prefixes should be accepted as aliases mapping to the current accepted binomial. The old names still appear in UniProt and published literature.

## 5. 5+5 latin prefix gaps

Several species in the mhcgnomes ontology work with 4-letter prefix but NOT with the 5+5 latin form:

| 4-letter (works) | 5+5 (fails) | Species |
|---|---|---|
| Chse | ChelySerpn | *Chelydra serpentina* |
| Gopo | GophePolyp | *Gopherus polyphemus* |
| Tetr | TerraTriun | *Terrapene triunguis* |
| Chab | ChelonAbing | *Chelonoidis abingdonii* |

Also, partial latin forms like `CheloMyda` (5+4) fail even though `CheloMydas` (5+5) works for *Chelonia mydas*.

**Advice**: Every species in the ontology should be addressable by any unambiguous prefix of its latin binomial, not just exactly 5+5. Prefix matching would make the system more robust.

## Summary of asks (prioritized)

1. **Fix 3 prefix collisions** (`Cyca`, `MHC-B` default, `ORLA` → orangutan)
2. **Add ~120 gene definitions** for already-known species (chicken MHC-Y, quail numbered loci, fish D-series, marsupial U-genes, crocodilian DB01-DB08)
3. **Add ~328 species** (all literature-attested, latin 5+5 forms provided above)
4. **Add 6 taxonomy aliases** for reclassified genera
5. **Ensure all 5+5 latin forms work** for existing species + support partial prefix matching
