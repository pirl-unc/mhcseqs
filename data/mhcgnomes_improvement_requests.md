# mhcgnomes Improvement Requests from mhcseqs v0.9.0

Generated: 2026-03-17
Tested against: mhcgnomes 3.8.0
Dataset: 15,861 diverse MHC entries (1,610 unique gene+organism pairs)

## Context

mhcseqs curates MHC sequences from UniProt for species underrepresented in
IMGT/HLA and IPD-MHC. For every entry, we have the organism name from UniProt
metadata. We plan to pass the organism as `species` when calling
mhcgnomes, which eliminates prefix collision issues.

## Parsing summary (mhcgnomes 3.8.0)

Tested 1,610 unique (gene, organism) pairs with both as-is and
bare-gene + `species` strategies.

| Strategy | Count | % |
|---|---|---|
| Parsed correctly as-is | 267 | 16.6% |
| Fixed by `species` | 112 | 7.0% |
| **Still fails** | **1,231** | **76.5%** |

Of the 1,231 remaining failures:

| Root cause | Count | Notes |
|---|---|---|
| Species not in ontology | 758 | Need species added |
| Species known, gene not recognized | 330 | Need gene definitions |
| Doubled prefix (our data bug) | 143 | mhcseqs issue, not mhcgnomes |

## 1. Hierarchical gene definitions (primary recommendation)

Of the 758 unknown-species failures, we tested what would happen if those
species were added. **398 (53%) would immediately parse** because the gene
names (DAB, UA, DRB, etc.) are already known to mhcgnomes for other species.
**360 (47%) would still fail** because the gene names aren't recognized for
any species.

This suggests that MHC gene names should be defined at higher taxonomic
levels and inherited by species:

### Tier 1: Vertebrate-wide gene names

These genes appear across fish, amphibians, reptiles, birds, and mammals.
They should be accepted for **any** vertebrate species in the ontology:

| Gene | Class | Chain | Count in our data | Example species using it |
|---|---|---|---|---|
| UA | I | alpha | 36 | *Agalychnis callidryas* (frog), *Sphenodon punctatus* (tuatara) |
| UAA | I | alpha | 14 | *Andrias davidianus* (salamander), *Clarias batrachus* (catfish) |
| UBA | I | alpha | 6 | *Acipenser sinensis* (sturgeon), *Chroicocephalus scopulinus* (gull) |
| UCA | I | alpha | 2 | *Chroicocephalus scopulinus*, *Clarias batrachus* |
| UDA | I | alpha | 3 | *Chroicocephalus scopulinus*, *Oryzias dancena* (ricefish) |
| DAB | II | beta | 53 | *Acipenser dabryanus* (sturgeon), *Agelaius phoeniceus* (blackbird) |
| DAB1 | II | beta | 42 | *Abramis brama* (bream), *Falco cherrug* (falcon) |
| DAB2 | II | beta | 16 | *Amblyrhynchus cristatus* (iguana), *Ardea alba* (egret) |
| DAA | II | alpha | 15 | *Andrias davidianus*, *Nipponia nippon* (ibis) |
| DRA | II | alpha | 8 | *Laticauda laticaudata* (sea krait), *Naja naja* (cobra) |
| DRB | II | beta | 24 | *Aegypius monachus* (vulture), *Astur gentilis* (goshawk) |
| DRB1 | II | beta | 7 | *Alca torda* (razorbill), *Brachyramphus marmoratus* (murrelet) |
| DMA | II | alpha | 7 | *Ciconia boyciana* (stork), *Gopherus evgoodei* (tortoise) |
| DMB | II | beta | 6 | *Eudyptes chrysocome* (penguin), *Eudyptula minor* (penguin) |
| DMB1 | II | beta | 2 | *Ciconia boyciana*, *Tympanuchus cupido* (grouse) |
| DMB2 | II | beta | 2 | *Ciconia boyciana*, *Tympanuchus cupido* |
| DBA | II | alpha | 3 | *Ciconia boyciana*, *Oceanodroma leucorhoa* (petrel) |
| DBB | II | beta | 2 | *Ciconia boyciana*, *Rhinella marina* (toad) |
| DXB | II | beta | 25 | *Coregonus* sp. (whitefish) — standard fish class II |
| MHCIIB | II | beta | 25 | *Alcedo atthis* (kingfisher), *Ardea cinerea* (heron) |
| B | I | alpha | 5 | *Chrysolophus amherstiae* (pheasant), *Pavo cristatus* (peafowl) |
| BF | I | alpha | 4 | *Lophura nycthemera* (pheasant), *Syrmaticus ellioti* (pheasant) |
| BLB2 | II | beta | 3 | *Lagopus scotica* (grouse), *Lyrurus tetrix* (grouse) |
| F | I | alpha | 5 | *Carduelis carduelis* (goldfinch), *Ophiophagus hannah* (king cobra) |
| A | I | alpha | 2 | *Anser indicus* (goose), *Ophiophagus hannah* |

### Tier 2: Order-level gene definitions

These genes are specific to particular taxonomic orders:

**Crocodylia — DB-series (class II beta)**

| Gene | Count | Species |
|---|---|---|
| DB01 | 4 | *C. porosus*, *C. niloticus*, *O. tetraspis*, *M. cataphractus* |
| DB02 | 4 | same |
| DB03 | 4 | same |
| DB04 | 4 | same |
| DB05 | 4 | same |
| DB06 | 4 | same |
| DB07 | 3 | *C. porosus*, *O. tetraspis*, *M. cataphractus* |
| DB08 | 3 | *C. porosus*, *O. tetraspis*, *M. cataphractus* |

From Jaratlerdsiri et al. 2014. Trans-specific class II beta lineages
across all Crocodylia. Any crocodilian species should accept DB01–DB08.

**Galliformes — BF/BLB workshop and serological variants**

Already partly handled for chicken (Gaga-BF works), but workshop alleles
(BFw-NN), z-series (BFz-NN), and serological types (B-F-SNN) don't parse.
These should work for all galliform species (chicken, quail, turkey,
pheasant, grouse).

### Tier 3: Genus-level gene definitions

**Turdus — PFA alleles (class I)**

| Gene | Count |
|---|---|
| PFA01–PFA30 | 47 entries across *T. naumanni*, *T. eunomus*, *T. ruficollis*, *T. atrogularis* |

Published thrush MHC class I nomenclature. Pattern: `PFA\d+`.

**Coturnix — numbered class II**

| Gene | Count |
|---|---|
| II-01 through II-17 | 17 entries for *Coturnix japonica* |

Hosomichi et al. 2006. Pattern: `II-\d+`. Quail-specific.

## 2. Standard genes that fail for known species (330 failures)

Even with correct `species`, standard gene names fail for many
species mhcgnomes already knows:

```
               UA   UAA  UBA  DAB  DAA  DRA  DRB  DMA  DMB  DBB  DBA
Croc. porosus  OK    -    -    -   OK    -    -    -    -    -    -
Chelonia mydas  -    -    -    -    -    -    -    -    -    -    -
Chelydra serp.  -    -    -    -    -    -    -    -    -    -    -
Struthio cam.   -    -    -    -    -    -    -    -    -    -    -
Apteryx owen.   -    -    -    -    -    -    -    -    -    -    -
Danio rerio     -    -   OK   OK   OK    -    -    -    -   OK    -
Oreochromis n.  -   OK   OK   OK   OK    -    -    -    -   OK   OK
Salmo salar     -    -   OK   OK   OK    -    -    -    -    -    -
Monodelphis d. OK    -    -   OK   OK   OK    -    -    -   OK   OK
Sarcophilus h. OK    -    -   OK    -    -    -    -    -    -    -
Gallus gallus   -    -    -    -    -    -    -   OK    -    -    -
Anas platyrh.   -   OK    -    -    -    -    -    -    -    -    -
Xenopus trop.   -   OK   OK   OK   OK    -    -   OK    -    -    -
```

**Sea turtle, snapping turtle, ostrich, and kiwi accept ZERO standard MHC
gene names** despite being in the species ontology. Croc only accepts 2/11.

With hierarchical gene definitions (Tier 1 above), all species in the
ontology would automatically accept the vertebrate-wide gene names.

## 3. Verify `species` param is authoritative

With the new `species` parameter (replacing `default_species`), the parser
should never return a different species than the one passed. Previously:

```python
mhcgnomes.parse("MHCIIB", default_species="Struthio camelus")
# Returned: Gene(species='Tyto alba', name='DAB') — WRONG
```

The new `species` param should make this impossible by design.

## 4. Unknown species (328 prefixes, 758 gene entries)

These species are not in the ontology. With hierarchical gene definitions,
adding a species would immediately give it access to all vertebrate-wide
and order/genus-level gene names.

Top 15 by gene count, all literature-attested in UniProt:

| Prefix | Count | Organism | Latin 5+5 |
|---|---|---|---|
| Dila | 62 | *Dicentrarchus labrax* (European seabass) | DicenLabra |
| Bain | 25 | *Labeobarbus intermedius* (Lake Tana barbs) | LabeoInter |
| Char | 19 | *Channa argus* (snakehead) | ChannArgus |
| gamr | 19 | *Gadus morhua* (Atlantic cod) | GadusMorhu |
| Phci | 19 | *Phascolarctos cinereus* (koala) | PhascCiner |
| Asme | 18 | *Astyanax mexicanus* (blind cavefish) | AstyaMexic |
| Sppu | 16 | *Sphenodon punctatus* (tuatara) | SphenPunct |
| Tuna | 16 | *Turdus naumanni* (dusky thrush) | TurduNauma |
| Angr | 16 | *Anabarilius grahami* (kanglang fish) | AnabaGraha |
| Sias | 14 | *Silurus asotus* (Amur catfish) | SilurAsotu |
| Epco | 12 | *Epinephelus coioides* (orange-spotted grouper) | EpineCoioi |
| Crmo | 10 | *Crocodylus moreletii* (Morelet's croc) | CrocoMorel |
| Icpu | 10 | *Clarias magur* (Asian catfish) | ClariMagur |
| Mega | 9 | *Meleagris gallopavo* (wild turkey) | MeleaGallo |
| Ocle | 9 | *Oceanodroma leucorhoa* (Leach's petrel) | OceanLeuco |

Notable phylogenetically important species:

| Organism | Prefix | Notes |
|---|---|---|
| *Sphenodon punctatus* (tuatara) | Sppu | Only living rhynchocephalian |
| *Crocodylus moreletii* + 5 spp. | Crmo etc. | Literature prefixes, Jaratlerdsiri 2014 |
| *Meleagris gallopavo* (turkey) | Mega | Economically important poultry |
| *Ornithorhynchus anatinus* (platypus) | Oran | Monotreme |
| *Tachyglossus aculeatus* (echidna) | Taac | Monotreme |
| *Amblyrhynchus cristatus* (marine iguana) | Amcr | Galapagos endemic |

Full list of 328 prefixes with organisms in `data/mhcgnomes_3.5.0_evaluation.md`.

## 5. Entries that will never parse (not mhcgnomes issues)

These bare gene names can't be meaningfully parsed:

| Category | Count | Example | Why |
|---|---|---|---|
| Scaffold/genome gene IDs | 49 | `C0J50_7051`, `EXN66_Car000162` | Not MHC names, NCBI locus tags |
| Bare numbers | 6 | `01`, `94`, `134` | Allele IDs without locus name |
| Single letters | 4 | `r`, `l`, `I` | Too ambiguous |

These need to be handled by mhcseqs (pass full context via `species`
and accept that the gene name is opaque).

## 6. Minor issues

### 5+5 latin prefix gaps

| 4-letter (works) | 5+5 (fails) | Species |
|---|---|---|
| Chse | ChelySerpe | *Chelydra serpentina* |

### Taxonomy alias

| Old prefix | Current binomial | Status in 3.8.0 |
|---|---|---|
| Maeu | *Notamacropus eugenii* (wallaby) | NOT KNOWN |

## Summary of asks (prioritized)

1. **Hierarchical gene definitions** — define standard MHC gene names at
   vertebrate/order/genus level so species inherit them automatically.
   This fixes 398 unknown-species entries immediately when species are added,
   and fixes 330 known-species entries that currently fail on standard genes.
2. **Verify new `species` param** — confirm it's authoritative and can't
   return a different species than passed.
3. **Add 328 species** to the ontology (all literature-attested, latin 5+5
   forms provided).
4. **Add order/genus gene patterns** — croc DB01–DB08, quail II-NN, thrush PFA-NN.
5. **Fix 5+5 gaps and taxonomy alias** (ChelySerpe, Maeu).
