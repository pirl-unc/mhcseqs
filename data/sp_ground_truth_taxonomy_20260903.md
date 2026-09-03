# Signal-peptide ground-truth taxonomy audit (2026-09-03)

## Scope and source

Issue #63 found that the SP benchmark discarded the source taxonomy clade and
then reconstructed broad categories from handwritten genus/common-name lists.
That approach mislabeled unfamiliar taxa and allowed substring collisions.

All 253 unique taxon IDs in the 2,403-row raw corpus were resolved through the
official UniProt taxonomy endpoint on 2026-09-03. Every response reported
UniProt release `2026_03` (released 2026-09-02). The accession-independent
decision, lineage roots, release, and exact source URL are preserved in
`sp_ground_truth_taxonomy.csv`.

The audit covers all 2,402 rows in the enriched benchmark. The pre-existing
raw-only non-MHC record `Q99004` is tracked separately in issue #66; this
taxonomy-only change does not alter its inclusion status.

## Source-clade coverage

The raw fetcher queries six vertebrate lineage roots. All stored taxa resolve
to exactly one of them:

| Source clade | Root taxon ID | Raw rows | Enriched rows |
| --- | ---: | ---: | ---: |
| Mammalia | 40674 | 497 | 496 |
| Aves | 8782 | 500 | 500 |
| Actinopterygii | 7898 | 500 | 500 |
| Reptilia | 8504 | 470 | 470 |
| Amphibia | 8292 | 359 | 359 |
| Chondrichthyes | 7777 | 77 | 77 |

The exact 500-row Aves and Actinopterygii totals expose an independent missing
pagination problem, tracked in issue #67.

## Reporting-category rules

Non-mammalian reporting categories derive directly from the six source clades:
Aves is `bird`; Actinopterygii and Chondrichthyes are `fish`; Reptilia and
Amphibia are `other_vertebrate`.

Mammalian subcategories use explicit lineage roots while preserving the
package's established category semantics:

- `human`: Homo sapiens (9606)
- `nhp`: Primates (9443), after the human special case
- `murine`: Murinae (39107)
- `ungulate`: Bovidae (9895), Suidae (9821), or Equidae (9789)
- `carnivore`: Canidae (9608) or Felidae (9681)
- `other_mammal`: all remaining Mammalia (40674)

## Migration result

The migration changes only `source_clade` and `species_category`. Class,
chain, gene, signal-peptide boundary, label status, and sequence are unchanged
for every enriched row. Of 2,402 rows, 332 category assignments change:

| Previous category | Lineage-backed category | Rows |
| --- | --- | ---: |
| other_vertebrate | fish | 205 |
| other_vertebrate | bird | 94 |
| other_vertebrate | other_mammal | 21 |
| other_vertebrate | nhp | 8 |
| other_vertebrate | murine | 2 |
| other_vertebrate | carnivore | 1 |
| other_mammal | bird | 1 |

The last row is `Struthio camelus australis`: the old name matcher treated the
ostrich as a mammal because its species epithet matched the camel genus. Known
fish misses such as `Astyanax mexicanus` and `Sinocyclocheilus rhinocerous`
now correctly report as `fish` without adding new genus hints.

Regenerating all accession metadata also exposed loss of archived class/chain
labels for inactive UniProt entries. That separate provenance problem is
tracked in issue #68; the taxonomy migration deliberately preserves the
existing scientific labels.

The checked-in SP boundary model was regenerated from the corrected groups
without changing its 2,402-row training population. Issue #69 tracks aligning
that population with the curated/gold eligibility policy in a separate change.
The lineage-backed model produces 1,825 exact boundaries among 2,186 parsed
curated/gold rows (83.5%), compared with 1,824 (83.4%) before this migration.
