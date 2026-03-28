# mhcgnomes Gene Resolution Gaps

Generated from 2403 GT + 2155 control sequences. These gene names appear with
gold-standard class/chain annotations but `mhcgnomes.infer_mhc_class()` returns empty.

## Gene suffix → class/chain rules

These suffixes are consistent across species and should be recognized regardless of
the 4-letter species prefix (e.g., `Povi-F10`, `Geja-F10`, `Orni-F10` are all class I alpha).

### Class I alpha (always)
| Suffix | Count | Notes |
|--------|-------|-------|
| F10 | 344 | Non-classical class I, ubiquitous across reptiles/fish/birds |
| Q9 | 120 | Non-classical class I |
| E-S | 110 | Non-classical class I |
| A-U | 102 | Classical/uncharacterized class I |
| Q10 | 48 | Non-classical class I |
| UAA | 48 | Class I alpha (amphibians) |
| U | 32 | Class I (newts/salamanders) |
| UA | 20 | Class I alpha |
| A-Q | 18 | Class I |
| UB | 18 | Class I alpha |
| I | 12 | Class I |
| D | 10 | Class I |
| I-E | 6 | Class I |
| E | 6 | Class I |
| mr1-like | 8 | MR1-like class I |
| MR1 | 4 | MR1 |
| Mill1, Mill2 | 4 each | Mill family class I |

### Class II alpha (always)
| Suffix | Count | Notes |
|--------|-------|-------|
| DRA | 106 | DR alpha |
| DAA | 90 | DA alpha (amphibians/fish) |
| DPA1 | 18 | DP alpha |
| M | 16 | Class II alpha (reptiles) |
| DR-1 | 14 | DR alpha variant |
| DMA | 8 | DM alpha |
| Aa | 6 | Class II alpha (mouse-style) |
| Ea | 4 | Class II alpha (mouse-style) |

### Class II beta (always)
| Suffix | Count | Notes |
|--------|-------|-------|
| BLB | 120 | B-locus beta (birds) |
| DAB | 30 | DA beta (amphibians/fish) |
| DRB1-4 | 18 | DR beta |
| M1 | 12 | Class II beta (reptiles) |
| DNB | 12 | DN beta (newts) |
| Eb1 | 6 | Class II beta (mouse-style) |
| DRB1-16 | 6 | DR beta |
| hla-drb1 | 6 | DR beta (Xenopus tropicalis naming) |
| DMB | 4 | DM beta |
| DQB1 | 4 | DQ beta |
| DRB5 | 4 | DR beta |
| DRB3 | 4 | DR beta |
| DRB1-10, -12, -15 | 4 each | DR beta |

### Ambiguous
| Suffix | Count | Notes |
|--------|-------|-------|
| A | 14 | Mixed: 8 class I + 6 class II beta |
| B | 6 | Mixed: 4 class I + 2 class II alpha |
| MHCY2B7 | 4 | Unknown (chicken Y-region) |

## Species prefixes (139 total)

4-letter species codes not in mhcgnomes: Anda (Andrias davidianus), Pifi (Pipra filicauda),
Povi (Pogona vitticeps), Pagu (Pantherophis guttatus), Orni (Oreochromis niloticus),
Nosc (Notechis scutatus), Geja (Gekko japonicus), Acda (Acipenser dabryanus),
Limv (Lissotriton vulgaris), Euma (Eublepharis macularius), Safa (Salarias fasciatus),
Dare (Danio rerio MHC-specific names), and 127 more.

## Species unknown to mhcgnomes (top 30 by sequence count)

| Species | Sequences |
|---------|-----------|
| Pogona vitticeps | 132 |
| Andrias davidianus | 128 |
| Pantherophis guttatus | 120 |
| Pipra filicauda | 94 |
| Eublepharis macularius | 90 |
| Crotalus adamanteus | 72 |
| Notechis scutatus | 68 |
| Gekko japonicus | 56 |
| Rhinella marina | 50 |
| Sparus aurata | 46 |
| Lates calcarifer | 42 |
| Salarias fasciatus | 42 |
| Acipenser dabryanus | 40 |
| Python bivittatus | 40 |
| Podarcis lilfordi | 40 |
| Lissotriton vulgaris | 40 |
| Microcaecilia unicolor | 38 |
| Salvator merianae | 36 |
| Scomber scombrus | 32 |
| Strix occidentalis caurina | 31 |
