# mhcgnomes Next Round: Numbered Loci and Clade Gaps

Tested against: mhcgnomes 3.11.0 (34.6% parse rate, 0 wrong species)
Remaining known-species failures: 269

## 1. Numbered paralogous loci (31 gene names, ~55 entries)

The base gene (DAB, UAA, DMB, etc.) parses correctly for these species,
but the numbered locus variant doesn't. These are **separate paralogous
genes**, not alleles — DAB1 and DAB3 are distinct genomic loci, like how
DRB1 and DRB3 are for human.

They should be added as **independent gene definitions alongside their
base gene** on the same clade nodes. No aliasing — `DAB` and `DAB1` are
both valid, meaning different things (unresolved family vs specific locus).

### D-series class II (place on same nodes as their base gene)

| Gene | Base | Fails | Clades affected |
|---|---|---|---|
| DAB1 | DAB | 12 | birds (4), fish (6), marsupial (1), reptile (1) |
| DAB2 | DAB | 4 | birds (2), fish (1), reptile (1) |
| DAB3 | DAB | 6 | birds (1), fish (4), reptile (1) |
| DAB4 | DAB | 1 | reptile — *Amblyrhynchus cristatus* |
| DAB5 | DAB | 1 | reptile — *Amblyrhynchus cristatus* |
| DAA1 | DAA | 1 | fish — *Ictalurus punctatus* |
| DAA-1 | DAA | 1 | reptile — *Sphenodon punctatus* (hyphen, not number suffix) |
| DAA-2 | DAA | 1 | reptile — *Sphenodon punctatus* |
| DBA1 | DBA | 1 | fish — *Ictalurus punctatus* |
| DBB1–4 | DBB | 4 | fish — *Scophthalmus maximus* |
| DCB2 | DCB | 1 | fish — *Tropheus sp.* |

**Recommendation**: Add DAB1, DAB2, DAB3 at the same clade nodes where
DAB is currently defined (Actinopterygii, Amphibia, Reptilia, Aves,
Marsupialia). DAB4/DAB5 are rare — maybe just on Reptilia. DAA1, DBA1,
DBB1–4, DCB2 on Actinopterygii.

### U-series class I

| Gene | Base | Fails | Clades |
|---|---|---|---|
| UAA1–4 | UAA | 6 | fish — *Oryzias dancena*, *O. luzonensis* |
| UBA2 | UBA | 1 | fish — *Oreochromis niloticus* |
| UBA4 | UBA | 1 | amphibian — *Xenopus tropicalis* |

**Recommendation**: Add UAA1–UAA4, UBA2 on Actinopterygii. UBA4 on Amphibia.

### DM/DR numbered variants

| Gene | Base | Fails | Clades |
|---|---|---|---|
| DRB1 | DRB | 5 | amphibian (2), bird (1), fish (1), turtle (1) |
| DMB1 | DMB | 1 | bird — *Ciconia boyciana* |
| DMB2 | DMB | 1 | bird — *Ciconia boyciana* |
| DMB-0, DMB-1 | DMB | 2 | bird — *Spheniscus mendiculus* |
| DMA1 | DMA | 1 | bird — *Coturnix japonica* |

**Recommendation**: DRB1 is used very broadly (amphibians, birds, fish,
turtles) — add wherever DRB is defined. DMB1/DMB2/DMA1 on Aves.

### Allele-field forms (3 entries)

`DAB1*06`, `DAB2*0102`, `DAB2*04` for *Phasianus colchicus* (pheasant).
These are alleles of DAB1/DAB2 with the `*field` notation. If DAB1 and
DAB2 are defined for pheasant (via Galliformes), these should parse
automatically as alleles.

## 2. Genes missing from specific clade nodes (~45 entries)

These are genes where the base form itself fails — the gene isn't defined
on the clade node for that species.

### UAA — missing from Actinopterygii and Chondrichthyes (9 failures)

Currently UAA works for some fish (Oryzias, Oreochromis, Salmo) but not:
- *Nothobranchius* spp. (5 species) — killifish
- *Iconisemion striatum* — killifish
- *Squalus acanthias* — spiny dogfish (shark)
- *Triakis scyllium* — banded houndshark
- *Gopherus polyphemus* — gopher tortoise

The killifish and sharks suggest UAA isn't on the right ancestor node.
Sharks (Chondrichthyes) are a separate branch from bony fish. Check that
UAA is on a node that covers both.

### UA — missing from Amphibia and some Aves (5 failures)

Works for reptiles and crocs, fails for:
- *Aquarana catesbeiana* (bullfrog) — Amphibia
- *Grus grus* (common crane) — Aves
- *Carassius langsdorfii* — fish
- *Limosa limosa* (godwit) — Aves
- *Spinus spinus* (siskin) — Aves

Suggests UA isn't on Amphibia or isn't inherited by all Aves nodes.

### DAB — 3 species where even the base fails

- *Emberiza cioides* (meadow bunting) — bird
- *Grus grus* (common crane) — bird
- *Hypophthalmichthys molitrix* (silver carp) — fish

These species might not be under the right parent node in the tree.

### DQB, DPB1, DPA1, DRA, DOA, DOB — bats (Chiroptera)

4 bat species (*Rhinolophus ferrumequinum*, *R. macrotis*, *Myotis lucifugus*,
*Myotis brandtii*, *Hipposideros armiger*) fail on standard mammalian class II
genes. These are all Chiroptera — check that bat species inherit from
a node that has DQB/DPB1/DPA1/DOA/DOB.

### B2M — 3 birds

*Leucopsar rothschildi*, *Menura novaehollandiae*, *Rhodinocichla rosea* fail
on B2M which should be on Gnathostomata root. These birds may not be
connected to the tree properly.

### BLB1, BLB2, BLB3 — *Lagopus scotica* (grouse)

BLB is defined for Gallus gallus. Grouse is Galliformes — should inherit.
Check that *Lagopus* is under the Galliformes node.

### Marsupial U-lineage loci (UB, UC, UD, UE, UF, UG, UM)

*Phascolarctos cinereus* (koala) and *Notamacropus eugenii* (wallaby)
fail on UB, UC, UD, UE, UF, UG, UM. These are marsupial-specific
U-lineage class I paralogs from Belov et al. (2006). Add to Marsupialia.

### F — reptiles and fish

*Crotalus adamanteus* (rattlesnake), *Ophiophagus hannah* (king cobra),
*Ictalurus punctatus* (catfish). `F` is a class I gene (like HLA-F).
Check if it's defined broadly enough.

### Shark-specific

*Triakis scyllium* (houndshark): UBA fails.
*Ginglymostoma cirratum* (nurse shark): UAA01, UAA03, UAA05, UFA fail.
Sharks (Chondrichthyes) may need their own clade node with U-series genes.

## 3. Species likely not connected to the tree (3 entries failing on B2M)

If B2M is at Gnathostomata root, every vertebrate should inherit it. These
3 bird species failing on B2M suggests they're not connected:
- *Leucopsar rothschildi* (Bali myna)
- *Menura novaehollandiae* (superb lyrebird)
- *Rhodinocichla rosea* (rosy thrush-tanager)

Quick check: are these species in the ontology at all? If so, are they
under a node that descends from Gnathostomata?

## Summary of changes

| Change | Entries fixed | Effort |
|---|---|---|
| Add DAB1/2/3 at same nodes as DAB | ~22 | Low — copy existing DAB definition |
| Add DRB1 at same nodes as DRB | ~5 | Low |
| Add UAA1–4, UBA2 on Actinopterygii | ~6 | Low |
| Add DMB1/2, DMA1 on Aves | ~4 | Low |
| Add DBB1–4, DCB2, DBA1, DAA1 on Actinopterygii | ~7 | Low |
| Fix UAA coverage (killifish, sharks) | ~9 | Check tree connections |
| Fix UA on Amphibia | ~5 | Check tree connections |
| Fix bat (Chiroptera) class II inheritance | ~8 | Check tree connections |
| Add marsupial U-loci (UB–UM) on Marsupialia | ~7 | Low |
| Fix B2M for 3 disconnected birds | ~3 | Check tree connections |
| Fix BLB on Galliformes for grouse | ~3 | Check tree connections |
| Add DAA-1/DAA-2 style (hyphenated) | ~2 | Low |
| **Total** | **~81** | |
