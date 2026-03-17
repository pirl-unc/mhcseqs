# mhcgnomes Improvement Requests from mhcseqs v0.9.0

Generated: 2026-03-17
Tested against: mhcgnomes 3.7.0
Dataset: 15,861 diverse MHC entries (1,526 unique gene names)

## Context

mhcseqs curates MHC sequences from UniProt for species underrepresented in
IMGT/HLA and IPD-MHC. For every entry, we have the organism name from UniProt
metadata. We plan to pass the organism as `default_species` when calling
mhcgnomes, which eliminates prefix collision issues — but requires mhcgnomes
to accept full latin binomials as `default_species` values.

## Parsing summary (mhcgnomes 3.7.0)

| Result | Count | % |
|---|---|---|
| Parsed correctly | 255 | 16.7% |
| Wrong species (prefix collision) | 8 | 0.5% |
| Failed | 1,263 | 82.8% |

Note: 3.7.0 fixed many issues from 3.5.0 — chicken MHC-Y, zebrafish lineage
genes, opossum genes, tilapia DBA, and taxonomy aliases (Brre→Danio,
Trsi→Pelodiscus) now all parse correctly.

## 1. Primary request: accept full latin binomials in `default_species`

This is the single most impactful change. Currently:

```python
# Works (species in ontology):
mhcgnomes.parse("UA", default_species="CrocoPoros")
# → Gene(species='Crocodylus porosus', name='UA')

# Fails (species not in ontology):
mhcgnomes.parse("UA", default_species="Cyclura carinata")
# → ERROR: Could not parse 'UA'
```

**Request**: When `default_species` is a full latin binomial (or long prefix)
that isn't in the ontology, create an ad-hoc `Species` from it rather than
failing. mhcseqs has the organism from UniProt for every entry — if mhcgnomes
accepted it, we could always disambiguate:

```python
# Desired:
mhcgnomes.parse("UA", default_species="Cyclura carinata")
# → Gene(species=Species(name='Cyclura carinata'), name='UA')
```

This solves prefix collisions (`Cyca` = carp vs iguana vs blue tit) without
needing every species in the ontology first. It also makes `MHC-B` on turkeys
and mammalian-prefix-on-wrong-species issues irrelevant, since we'd strip the
prefix and pass the bare gene with the correct `default_species`.

The species should still be required to be in the ontology — if there's a gene
ontology check, the ad-hoc species should support standard MHC gene names
(UA, DRA, DAB, etc.) even without species-specific gene definitions.

## 2. Gene definitions needed for known species

mhcgnomes 3.7.0 recognizes the species but can't parse these gene names.
Grouped by what's needed:

### Quail (*Coturnix japonica*, Coja) — 38 genes

```
Coja-II-17*01  → FAIL   (expected: class II beta allele)
Coja-DMB1      → FAIL   (expected: class II beta, DM gene)
Coja-B1        → FAIL   (expected: class I alpha)
Coja-D2        → FAIL   (expected: class II beta)
```

The `Coja-II-NN` format is quail-specific MHC class II nomenclature from
Hosomichi et al. 2006. A pattern rule for `II-\d+` would cover 17 genes.
DM genes (DMB1, DMB2, DMA1) and single-letter loci (B1, C, D, E) are standard.

### Crocodilian class II beta (DB01–DB08) — 26 genes across 4 species

```
Crpo-DB01  → FAIL   (expected: class II beta)
Oste-DB06  → FAIL   (expected: class II beta)
Meca-DB03  → FAIL   (expected: class II beta)
Crni-DB01  → FAIL   (expected: class II beta)
```

`DB01`–`DB08` are trans-specific class II beta lineages from Jaratlerdsiri
et al. 2014 (Immunogenetics 66:175–187). Should be defined as class II beta
gene patterns for all crocodilian species that mhcgnomes already knows
(Crpo, Crni, Meca, Oste).

### Fish U-lineage numbered loci — scattered

```
Orni-UAA1  → FAIL   (expected: class I alpha, locus UAA1)
Ctid-UAI02 → FAIL   (expected: class I alpha, locus UAI02)
```

Some fish species have numbered U-lineage loci (UAA1, UBA1, UAI01, etc.)
that aren't recognized. The base forms (UAA, UBA) parse for some species
but the numbered variants don't.

### Duck (*Anas platyrhynchos*, Anpl) — 6 genes

```
Anpl-DRA  → FAIL   (expected: class II alpha)
Anpl-UEA  → FAIL   (expected: class I alpha)
Anpl-UBA  → FAIL   (expected: class I alpha)
```

Standard MHC gene names that parse for other species but not duck.

### Tasmanian devil (*Sarcophilus harrisii*, Saha) — 5 genes

```
Saha-I-01   → FAIL   (expected: class I alpha)
Saha-UC-01  → FAIL   (expected: class I alpha)
```

Numbered format unique to Tasmanian devil MHC nomenclature.

### Sea turtle scaffold gene IDs

```
Chmy-UY3_17009  → FAIL   (expected: class I alpha)
```

These are genome assembly scaffold gene identifiers, not standard MHC names.
Probably not worth adding as gene definitions — better handled by the
`default_species` approach where we pass the bare gene name.

## 3. Unknown species (328 unique prefixes)

These species are not in the mhcgnomes ontology at all. With the
`default_species` enhancement, many of these would work immediately for
standard gene names (UA, DRA, DAB, etc.). But adding them to the ontology
would also enable 4-letter prefix parsing.

Top 15 by gene count, all literature-attested in UniProt:

| Prefix | Count | Organism | Latin 5+5 |
|---|---|---|---|
| Dila | 62 | *Dicentrarchus labrax* (European seabass) | DicenLabra |
| Bain | 25 | *Labeobarbus intermedius* (Lake Tana barbs) | LabeoInter |
| Char | 19 | *Channa argus* (snakehead) | ChanArgus |
| gamr | 19 | *Gadus morhua* (Atlantic cod) | GadusMorhu |
| Phci | 19 | *Phascolarctos cinereus* (koala) | PhasCinere |
| Asme | 18 | *Astyanax mexicanus* (blind cavefish) | AstyaMexic |
| Sppu | 16 | *Sphenodon punctatus* (tuatara) | SphenPunct |
| Tuna | 16 | *Turdus naumanni* (dusky thrush) | TurduNauma |
| Angr | 16 | *Anabarilius grahami* (kanglang fish) | AnabGraha |
| Sias | 14 | *Silurus asotus* (Amur catfish) | SilurAsotu |
| Epco | 12 | *Epinephelus coioides* (orange-spotted grouper) | EpineCoioi |
| Crmo | 10 | *Crocodylus moreletii* (Morelet's croc) | CrocoMorel |
| Icpu | 10 | *Clarias magur* (Asian catfish) | ClariMagur |
| Mega | 9 | *Meleagris gallopavo* (wild turkey) | MeleaGallo |
| Ocle | 9 | *Oceanodroma leucorhoa* (Leach's petrel) | OceanLeuco |

Notable archosaurs and other phylogenetically important species:

| Prefix | Organism | Notes |
|---|---|---|
| Sppu | *Sphenodon punctatus* (tuatara) | Only living rhynchocephalian |
| Crmo + 5 others | *Crocodylus* spp. | Literature prefixes from Jaratlerdsiri 2014 |
| Mega | *Meleagris gallopavo* (turkey) | Economically important poultry |
| Oran | *Ornithorhynchus anatinus* (platypus) | Monotreme |
| Taac | *Tachyglossus aculeatus* (echidna) | Monotreme |
| Amcr | *Amblyrhynchus cristatus* (marine iguana) | Galapagos endemic |
| Haal | *Haliaeetus albicilla* (white-tailed eagle) | Raptor |

Full list of 328 prefixes with organisms is in `data/mhcgnomes_3.5.0_evaluation.md`.

## 4. 5+5 latin prefix gaps

Several species mhcgnomes knows by 4-letter prefix don't work with the
5+5 latin form:

| 4-letter (works) | 5+5 (fails) | Species |
|---|---|---|
| Chse | ChelySerpn | *Chelydra serpentina* |
| — | CheloMyda (5+4) | *Chelonia mydas* — only exact 5+5 `CheloMydas` works |

3.7.0 fixed GophePolyp and TerraTriun.

**Suggestion**: Support prefix matching against the full binomial, not just
exact 5+5. Any unambiguous prefix of the latin name should resolve.

## 5. Remaining taxonomy alias

| Old prefix | Current binomial | Status in 3.7.0 |
|---|---|---|
| Maeu | *Notamacropus eugenii* (wallaby) | NOT KNOWN — old genus *Macropus* |

3.7.0 fixed Brre→Danio and Trsi→Pelodiscus. Only Maeu remains.

## Summary of asks (prioritized)

1. **Accept full latin binomials in `default_species`** — our primary request.
   Lets mhcseqs pass organism metadata to disambiguate collisions and handle
   unknown species. Species should still be in ontology, but standard gene
   names (UA, DRA, DAB, B2M, etc.) should work with an ad-hoc species.
2. **Add gene definitions** for quail (Coja-II-NN, DM), crocodilian (DB01–DB08),
   fish numbered U-loci, duck, Tasmanian devil — ~80 genes across known species.
3. **Add ~328 species** to the ontology (all literature-attested, latin 5+5
   forms provided above).
4. **Fix 5+5 latin gaps** — ChelySerpn, partial prefix matching.
5. **Add Maeu taxonomy alias** for *Notamacropus eugenii*.
