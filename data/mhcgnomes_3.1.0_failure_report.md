# mhcgnomes 3.1.0 Parse Failures from mhcseqs Diverse MHC Dataset

## Summary

mhcseqs curates 16,208 MHC sequences from UniProt covering taxonomic groups
underrepresented in IMGT/HLA and IPD-MHC (birds, fish, reptiles, amphibians,
sharks, marsupials, monotremes, bats). After upgrading to mhcgnomes 3.1.0:

| Metric | Count |
|--------|-------|
| Total curated sequences | 16,208 |
| Parsed successfully | 7,903 (48.8%) |
| Failed | 8,305 (51.2%) |
| Unique failed gene strings | 1,701 |
| Unique failed species prefixes | 592 |

### Failure categories

| Category | Sequences | Unique prefixes |
|----------|-----------|-----------------|
| Unknown species prefix | 7,086 | 560 |
| Known species, unknown gene | 1,219 | 32 |

## Category 1: Known species, unknown gene (1,219 sequences)

mhcgnomes recognizes the species prefix but cannot parse the gene name.
This is the higher-priority category since only gene definitions need adding.

### Gallus gallus (Gaga) — 287 sequences

Chicken MHC uses the B-locus nomenclature (BF for class I, BLB for class II)
which differs from the standard D-gene system. mhcgnomes knows `Gaga` but
doesn't recognize these gene names.

| Gene | Count | Class | Chain | Notes |
|------|-------|-------|-------|-------|
| Gaga-BF | 128 | I | alpha | B-F locus, class I heavy chain |
| Gaga-B-LB | 37 | II | beta | B-LB locus, class II beta |
| Gaga-B-F | 32 | I | alpha | Alternate hyphenation of BF |
| Gaga-B-LBII | 15 | II | beta | BLB roman numeral variant |
| Gaga-beta | 12 | II | beta | Generic beta chain |
| Gaga-BF12 | 1 | I | alpha | BF with serotype number |
| Gaga-BFIV21 | 1 | I | alpha | BF minor locus |
| Gaga-BLB | 1 | II | beta | No hyphen variant |
| Gaga-B-DMA | 2 | II | alpha | DM alpha |
| Gaga-B-DMB2 | 2 | II | beta | DM beta 2 |
| Gaga-YFV | 3 | I | alpha | Y-F locus (MHC-Y) |
| Gaga-YFVI | 2 | I | alpha | Y-F variant |
| Gaga-Y-FA | 1 | I | alpha | Y-F alternate format |
| Gaga-YL | 1 | I | alpha | Y-L locus |
| Gaga-MHCY2B2 | 3 | II | beta | MHC-Y class II |
| Gaga-MHCY2B1 | 2 | II | beta | MHC-Y class II |
| Gaga-MHCY15 | 2 | I | alpha | MHC-Y class I |
| Gaga-MHCY46 | 1 | I | alpha | MHC-Y class I |
| Gaga-MHCY11 | 1 | I | alpha | MHC-Y class I |
| Gaga-b2m-Z01..Z02 | 2 | I | B2M | Beta-2-microglobulin |
| Gaga-BFw-01..06 | 6 | I | alpha | BF-w (workshop) alleles |
| Gaga-BFz-01..05 | 5 | I | alpha | BF-z alleles |
| Gaga-B-F-S01..S14 | 14 | I | alpha | BF serological types |
| Gaga-B-F-S-b2m01..04 | 4 | I | B2M | B2M from serological study |
| Gaga-S19 | 1 | I | alpha | Serological type |
| Gaga-B-F-minor | 1 | I | alpha | Minor BF locus |
| Gaga-B-LBI | 3 | II | beta | BLB I |
| Gaga-B-LBVI | 1 | II | beta | BLB VI |
| Gaga-B-LB12c | 1 | II | beta | BLB12 variant c |

### Tyto alba (Tyal) — 156 sequences

mhcgnomes knows the species but fails on the genes.

| Gene | Count | Class | Chain |
|------|-------|-------|-------|
| Tyal-UA | 146 | I | alpha |
| Tyal-MhcTyal-DAB1 | 4 | II | beta |
| Tyal-MhcTyal-DAB2 | 2 | II | beta |
| Tyal-MHCIIB | 3 | II | beta |
| Tyal-DRB | 1 | II | beta |

Note: `Tyal-MhcTyal-DAB1` is a double-prefix pattern (see Category 3 below).

### Egretta eulophotes (Egeu) — 141 sequences

| Gene | Count | Class | Chain |
|------|-------|-------|-------|
| Egeu-DAB | 38 | II | beta |
| Egeu-DAB1 | 31 | II | beta |
| Egeu-DAB2 | 21 | II | beta |
| Egeu-DRA | 19 | II | alpha |
| Egeu-DAB3 | 8 | II | beta |
| Egeu-DAB4 | 8 | II | beta |
| Egeu-UAA | 9 | I | alpha |
| Egeu-UBA | 7 | I | alpha |

### Paralichthys olivaceus (Paol) — 86 sequences

mhcgnomes parses `Paol-DAB` and `Paol-DAA` but fails on numbered variants.

| Gene | Count | Class | Chain |
|------|-------|-------|-------|
| Paol-DAB1 | 8 | II | beta |
| Paol-DAB2 | 16 | II | beta |
| Paol-DAB3 | 1 | II | beta |
| Paol-DAB4 | 31 | II | beta |
| Paol-DAB5 | 27 | II | beta |
| Paol-DAB6 | 2 | II | beta |
| Paol-b2m | 1 | I | B2M |

### Oryzias latipes (Orla) — 82 sequences

| Gene | Count | Class | Chain |
|------|-------|-------|-------|
| Orla-DAB | 15 | II | beta |
| Orla-UGA | 13 | I | alpha |
| Orla-DAA | 11 | II | alpha |
| Orla-UAA | 10 | I | alpha |
| Orla-UBA | 9 | I | alpha |
| Orla-UHA | 6 | I | alpha |
| Orla-UIA1..3 | 7 | I | alpha |
| Orla-UDA | 1 | I | alpha |
| Orla-UCA | 1 | I | alpha |
| Orla-DCB | 1 | II | beta |
| Orla-DBB | 1 | II | beta |

### Danio rerio (Dare) — 61 sequences

Zebrafish uses a non-standard nomenclature (`mhc1uma`, `mhc2daa`, etc.).

| Gene | Count | Class | Chain | Notes |
|------|-------|-------|-------|-------|
| Dare-mhc1uma | 4 | I | alpha | MHC class I U-lineage |
| Dare-mhc1lla | 3 | I | alpha | MHC class I L-lineage |
| Dare-mhc1zea | 3 | I | alpha | MHC class I Z-lineage |
| Dare-mhc1zda | 3 | I | alpha | |
| Dare-mhc1lfa | 3 | I | alpha | |
| Dare-mhc2daa | 3 | II | alpha | Class II DAA |
| Dare-mhc2d4b | 2 | II | beta | |
| Dare-DBB | 2 | II | beta | |
| Dare-b2m | 2 | I | B2M | |
| Dare-UCA | 2 | I | alpha | |
| Dare-DAB | 1 | II | beta | |
| Dare-DAB1 | 1 | II | beta | |
| Dare-MHCII | 1 | II | unknown | |
| *(25 more mhc1/mhc2 variants)* | | | | |

### Coturnix japonica (Coja) — 48 sequences

Japanese quail. Has allele-level nomenclature.

| Gene | Count | Class | Chain | Allele examples |
|------|-------|-------|-------|-----------------|
| Coja-DMB1 | 1 | II | beta | |
| Coja-DMB2 | 1 | II | beta | |
| Coja-DMB2e2 | 3 | II | beta | |
| Coja-DMA1 | 1 | II | alpha | |
| Coja-B1 | 1 | I | alpha | |
| Coja-B1e3 | 3 | I | alpha | |
| Coja-B2M | 2 | I | B2M | |
| Coja-II-01..17 | 25 | II | beta | Coja-II-17\*01, Coja-II-13\*01 |
| Coja-C | 1 | I | alpha | |
| Coja-D | 1 | II | beta | |
| Coja-D1 | 1 | II | beta | |
| Coja-D2 | 1 | II | beta | |
| Coja-E | 1 | I | alpha | |
| Coja-QF41 | 1 | I | alpha | |
| Coja-QF63 | 1 | I | alpha | |

### Other known-species failures (< 50 sequences each)

| Prefix | Organism | Seqs | Example genes |
|--------|----------|------|---------------|
| Hymo | Hypophthalmichthys molitrix | 47 | UA, DAB |
| Orni | Oreochromis niloticus | 43 | DBA, UAA1, UBA1 |
| Gogo | Gobio gobio | 40 | DAB1, DAB3 |
| Anpl | Anas platyrhynchos | 31 | U, I, DRA, UEA, UBA |
| Xetr | Xenopus tropicalis | 27 | UAA, b2m, mhc1-uba13.2, mhc2-dab |
| Sphu | Spheniscus humboldti | 25 | MhcSphu_CLSI |
| Bubu | Bufo bufo / Buteo buteo | 22 | Bubu (prefix collision) |
| Xela | Xenopus laevis | 19 | DAB, mhc1-uea.L, mhc2-dab.L |
| Cyse | Cynoglossus semilaevis | 16 | DAB1, DAB2, b2m |
| Ctid | Ctenopharyngodon idella | 12 | UAD03..UAI02 (numbered loci) |
| Satr | Salmo trutta | 12 | DAA, UAA |
| Trvu | Trichosurus vulpecula | 12 | DAA, DBB, DRB2, DRB3 |
| ORLA | Nothobranchius spp. | 11 | UAA (uppercase prefix collision with Orla) |
| Acar | Acrocephalus arundinaceus | 9 | beta, B2M |
| Nini | Nipponia nippon | 7 | DAB, DBA, DAA, DCB, DBB |
| Stal | Strix aluco | 6 | MHCIIB, DRB |
| Grgr | Grus grus | 6 | DAB, UA |
| Chpi | Chrysemys picta bellii | 4 | B2M, DMA, DMB |

## Category 2: Unknown species prefix (7,086 sequences, 560 prefixes)

These are species where mhcgnomes doesn't know the prefix at all.
The full list is in `data/mhcgnomes_failures.csv`. The top 50 by sequence count:

| Prefix | Organism | Group | Seqs | Example genes |
|--------|----------|-------|------|---------------|
| Sthi | Sterna hirundo | bird | 206 | UA, DAB |
| Zhom | Zhangixalus omeimontis | amphibian | 162 | Rhom-beta1 |
| Spma | Spheniscus magellanicus | bird | 157 | DRB1 |
| Phtr | Phylloscopus trochilus | bird | 112 | UA |
| Cyca | Cyprinus carpio | fish | 105 | DAB1, DAB3 |
| Trov | Trachinotus ovatus | fish | 153 | DAA, DAB, UBA |
| Saha | Sarcophilus harrisii | marsupial | 112 | I, DAB, UC |
| Scar | Scatophagus argus | fish | 106 | DXB, DXA, UA |
| Dila | Dicentrarchus labrax | fish | 106 | UA, DAA, DAB |
| Phco | Phylloscopus collybita | bird | 98 | UA, B2M |
| Phci | Phascolarctos cinereus | marsupial | 94 | UA, DAB, UC |
| Ritr | Rissa tridactyla | bird | 87 | DRB1, DRB |
| Amci | Amphilophus citrinellus | fish | 76 | DXB |
| Sqce | Squalius cephalus | fish | 118 | DAB1, DAB3 |
| Scma | Scophthalmus maximus | fish | 73 | DBB1, DAB |
| Spsp | Spinus spinus | bird | 71 | UA |
| Ocle | Oceanodroma leucorhoa | bird | 71 | DAB1, DAB2 |
| Crpo | Crocodylus porosus | reptile | 70 | UA, UB, UC, DAA, B2M |
| Ruru | Rutilus rutilus | fish | 67 | DAB1, DAB3 |
| Napa | Nanorana parkeri | amphibian | 65 | DAB |
| Amme | Ambystoma mexicanum | amphibian | 65 | Mhc, DAB |
| Anda | Andrias davidianus | amphibian | 64 | UAA, DAA, DAB |
| Abbr | Abramis brama | fish | 64 | DAB1, DAB3 |
| RAJA | Rana japonica | amphibian | 60 | UA |
| SAAL | Salvelinus alpinus | fish | 200 | UBA, UGA, UEA |
| Modo | Monodelphis domestica | marsupial | 232 | UT3, UG, UT8 |
| Epco | Epinephelus coioides | fish | 276 | DAB, DBB, DAA |
| Fuat | Fulica atra | bird | 267 | DAB |
| Satr | Salmo trutta | fish | 246 | DAB, DAA |

*(Full list of 560 prefixes in data/mhcgnomes_failures.csv)*

## Category 3: Structural parse issues

These are patterns that may need special handling in the parser.

### Double-prefix gene names (274 sequences)

Some UniProt entries use an old species prefix concatenated with the gene,
then our curation adds the current species prefix, producing a double-prefix.

| Pattern | Count | Organism | Notes |
|---------|-------|----------|-------|
| Zhom-Rhom-beta1 | 162 | Zhangixalus omeimontis | Old genus: Rhacophorus (Rhom) |
| Liya-Raya-beta1 | 33 | Lithobates yavapaiensis | Old genus: Rana (Raya) |
| Pefl-PefuMhc-DXB | 33 | Perca fluviatilis | Old prefix: Pefu |
| Odto-Odto-beta2 | 29 | Odorrana tormota | Same prefix repeated |
| Ocle-MhcOcle-DAA | 9 | Oceanodroma leucorhoa | "Mhc" + prefix embedded |
| Ocle-MhcOcle-DBA | 9 | Oceanodroma leucorhoa | |
| Phci-Phci-mhc | 9 | Phascolarctos cinereus | Same prefix repeated |
| Noeu-Maeu-1..3 | 11 | Notamacropus eugenii | Old prefix: Macropus (Maeu) |
| Tyal-MhcTyal-DAB1..2 | 6 | Tyto alba | "Mhc" + prefix embedded |
| Sqac-MhcSqac-UAA | 2 | Squalus acanthias | "Mhc" + prefix embedded |

### Case-insensitive prefix issues

Some gene names use ALL-CAPS prefixes (e.g., `SAAL`, `RAJA`, `ORLA`, `MODO`)
rather than the standard capitalized form. mhcgnomes should ideally normalize
case during parsing.

| Uppercase prefix | Standard form | Organism | Seqs |
|------------------|---------------|----------|------|
| SAAL | Saal | Salvelinus alpinus | 200 |
| RAJA | Raja | Rana japonica | 60 |
| ORLA | Orla | Nothobranchius spp. | 11 |
| MODO | Modo | Monodelphis domestica | 10 |

### Prefix collisions

Some 4-letter prefixes map to multiple species from different genera:

| Prefix | Species |
|--------|---------|
| Bubu | Bufo bufo, Buteo buteo, Bulweria bulwerii |
| Cyca | Cyprinus carpio, Cyanistes caeruleus, Cyclura carinata |
| Chpi | Chrysolophus pictus, Chrysemys picta |
| Mimi | Miichthys miiuy, Milvus milvus |
| Phco | Phasianus colchicus, Phylloscopus collybita |

## Recommendations

1. **High impact**: Add the 560 unknown species prefixes. The full list with
   organisms and example gene strings is in `data/mhcgnomes_failures.csv`.

2. **Medium impact**: Add gene definitions for the 32 known-species failures,
   especially chicken B-locus genes (BF, BLB and variants) and numbered DAB
   variants (DAB1-DAB6).

3. **Low impact**: Handle double-prefix patterns by stripping redundant
   embedded prefixes (e.g., `Tyal-MhcTyal-DAB1` → `Tyal-DAB1`).

4. **Low impact**: Case-normalize prefixes during parsing (SAAL → Saal).

5. **Consider**: Establish a policy for prefix collisions where the same
   4-letter code maps to different genera.
