# Unknown Species Report for mhcgnomes

Generated: 2026-03-17
Tested against: mhcgnomes 3.9.2
Source: mhcseqs diverse_mhc_sequences.csv (15,861 entries)

## Summary

**630 species** with **1,034 gene entries** fail because mhcgnomes doesn't
know the species. All species have latin binomials from UniProt metadata.

Of the 1,034 entries:
- **388 (38%)** use standard MHC gene names (DAB, UA, DRA, etc.) that would
  parse immediately if the species were added to the ontology
- **446 (43%)** use gene names that mhcgnomes already inherits at the
  Gnathostomata level (B2M, DMA, DMB, F10, IA) — these would also work
  if species were added
- **200 (19%)** use non-standard or species-specific names that need
  gene definitions too

## Gene names that need adding to Gnathostomata root

Two very common gene names appear across many unknown species but aren't
in the mhcgnomes gene ontology:

### F10 — NOT recommended for mhcgnomes (mhcseqs issue)

`F10` is an NCBI genome annotation locus designation extracted from UniProt
protein descriptions like "Class I histocompatibility antigen, F10 alpha
chain-like". The UniProt gene_names for these entries are all `LOC*` locus
tags — `F10` appears only in the protein name, not as an actual gene symbol.

Our curation script (`_infer_gene_from_protein_name`) extracts it from the
protein description and treats it as a gene name. This is an mhcseqs bug:
F10 should be treated as an opaque annotation label, not passed to mhcgnomes.

**No action needed from mhcgnomes.** We will fix this in our curation pipeline.

### IA — NOT recommended for mhcgnomes (naming collision)

`IA` / `Ia` is already an H2 haplotype in mouse (I-region A, class II).
In non-mammalian species (frogs, fish, birds, sharks — 321 entries), `IA`
means "class I A" informally. Adding it as a gene would collide with the
mouse haplotype.

The non-mammalian uses of `IA` are shorthand for a class I alpha locus —
the actual gene is `UA` or `A`. mhcseqs should map `IA` → `UA` before
passing to mhcgnomes.

**No action needed from mhcgnomes.** We will fix this in our curation pipeline.

### Other gene names to add at root level

| Gene | Count | Class | Description |
|---|---|---|---|
| ZAA | 4 | I | Z-lineage class I alpha — legitimate gene, add at Gnathostomata |
| U | 1 | I | U-lineage (short form of UA) — add as alias for UA |
| DDA | 1 | II | D-series class II alpha — add at Gnathostomata |
| BLB | 1 | II | B-LB (chicken-style class II beta, no number) — add at Galliformes |
| UIA | 1 | I | U-lineage I-A locus — add at Actinopterygii |
| B2ML | 3 | I | B2M-like — mhcseqs should map to B2M |
| B2MG | 3 | I | B2M alternate name — mhcseqs should map to B2M |
| MHC | 9 | varies | Too generic — mhcseqs should treat as opaque |
| MHCI | 2 | I | Too generic — mhcseqs should treat as opaque |
| MHCIA | 1 | I | Too generic — mhcseqs should map to UA |
| MHCBETA | 1 | II | Too generic — mhcseqs should map to DAB |
| MHCIIA | 1 | II | Too generic — mhcseqs should map to DAA |

## Top 50 species to add

These are the species with the most gene entries. All have literature-attested
prefixes in UniProt. Adding these 50 species would cover 319 of 1,034 entries.

| Species | Prefix | Entries | Example genes |
|---|---|---|---|
| *Coregonus sp.* | Cosp | 17 | DXB |
| *Tropheus sp.* | Trsp | 12 | DCB2, DDB1, DDB2 |
| *Phasianus colchicus* | Phco | 11 | DAB\*02, DAB\*07, DMB, DRB2 |
| *Kareius bicoloratus* | Kabi | 11 | DAA, DAB, DRB2, DRB8, DRB10 |
| *Labeobarbus intermedius* | Bain | 10 | B2M, DAA\*01, DAB3 |
| *Ciconia boyciana* | Cibo | 8 | DAA, DAB, DBA, DBB, DMA, DMB1, DMB2 |
| *Amblyrhynchus cristatus* | Amcr | 7 | DAB1–4, UB |
| *Rhinella marina* | Rhma | 7 | DAA, DAB, DBA, DBB, DCB, DRA |
| *Scophthalmus maximus* | Scma | 7 | B2M, DAB, DBB1, DBB2 |
| *Monopterus albus* | Mhc | 7 | B2M, DAA, DAB, IA |
| *Rhinolophus ferrumequinum* | Rhfe | 7 | B2M, DOA, DOB, DPA1, DQB |
| *Agelaius phoeniceus* | Agph | 6 | B2M, DAB, DAB1, DAB03 |
| *Oceanodroma leucorhoa* | Ocle | 6 | DAA, DAB1, DAB2, DBA |
| *Eudyptula minor* | Eumi | 6 | DMB, DRA, F10 |
| *Geotrypetes seraphini* | Gese | 5 | B2M, DMA, DMB, DRB1 |
| *Columba livia* | Coli | 5 | B2M, DMB, DQB, F10 |
| *Tympanuchus cupido* | Iibi | 5 | DMA, DMB1, DMB2, I |
| *Egretta garzetta* | Egga | 5 | DAB, DAB1, DAB2, DMB |
| *Astyanax mexicanus* | Asme | 5 | B2M, DRB1, F10, UBA |
| *Sinocyclocheilus grahami* | Sigr | 5 | B2M, B2ML, F10, UBA |
| *Silurus asotus* | Sias | 5 | B2M, UBA, UDA, UEA |
| *Oryzias dancena* | Orda | 5 | UAA, UAA1, UAA3, UAA4 |
| *Oryzias luzonensis* | Orlu | 5 | UAA, UAA1, UAA2, UAA3 |
| *Clarias batrachus* | Clba | 5 | MHCIA, UAA, UBA, UCA |
| *Myotis lucifugus* | Mylu | 5 | B2M, DMA, DMB, DPA1 |
| *Eublepharis macularius* | Euma | 4 | B2M, DMA, DMB, F10 |
| *Pelobates cultripes* | Pecu | 4 | DMA, DMB, F10, IA |
| *Chroicocephalus scopulinus* | Lasc | 4 | UAA, UBA, UCA, UDA |
| *Pipra filicauda* | Pifi | 4 | B2M, DMA, DMB, F10 |
| *Spheniscus mendiculus* | Spme | 4 | DMB, DRB1, F10 |
| *Scatophagus argus* | Scar | 4 | B2M, DXA, DXB, UA |
| *Vombatus ursinus* | Vour | 4 | B2M, DMA, DMB, DRA |
| *Notamacropus* | Nota | 4 | DAB1, DAB2, DNA, DRA |
| *Hipposideros armiger* | Hiar | 4 | B2M, DMA, DMB, DQB |
| *Myotis brandtii* | Mybr | 4 | B2M, DMA, DMB, DQB |
| *Naja naja* | Nana | 3 | B2M, DMA, DRA |
| *Salvator merianae* | Same | 3 | B2M, DMA, DMB |
| *Pogona vitticeps* | Povi | 3 | B2M, DMA, F10 |
| *Ophiophagus hannah* | Opha | 3 | B2M, DRA, F |
| *Gopherus evgoodei* | Goev | 3 | B2M, DMA, F10 |
| *Crotalus adamanteus* | Crad | 3 | B2M, F, F10 |
| *Notechis scutatus* | Nosc | 3 | B2M, DMB, F10 |
| *Pantherophis guttatus* | Pagu | 3 | B2M, DMB, F10 |
| *Rutilus rutilus* | Arbr | 4 | DAB1, DAB3 |
| *Parachondrostoma toxostoma* | Pctn | 4 | DAB1, DAB3 |
| *Pseudonaja textilis* | Pste | 3 | B2M, DRA |
| *Laticauda laticaudata* | Lala | 3 | B2M, DRA |
| *Podarcis muralis* | Pomu | 3 | B2M, DMB, F10 |
| *Lacerta agilis* | Laag | 3 | B2M, F10 |
| *Anolis carolinensis* | Anca | 3 | B2M, DMA, F10 |

## Numbered loci that need gene definitions

These are standard-ish gene names with numbers/suffixes that mhcgnomes
doesn't recognize even for known species:

| Gene | Count | Notes |
|---|---|---|
| DDB1, DDB2 | 10 | Fish class II beta — D-series |
| DCB2 | 2 | Fish class II beta |
| UAA1–UAA4 | 10 | Fish U-lineage numbered loci |
| DRB2, DRB8, DRB10 | 7 | Numbered DRB loci |
| DAB\*02, DAB\*03, etc. | 16 | Allele designations within DAB |
| DMB-0, DMB-1 | 6 | Numbered DMB variants |
| DQB1 | 2 | Standard but not recognized for some species |

## Allele designations (16 entries)

These include allele-field notation (e.g. `DAB*02`, `DAA*0301`) and would
be parseable if the base gene (DAB, DAA) were accepted for the species.
mhcgnomes should strip allele fields from the bare gene before matching.
