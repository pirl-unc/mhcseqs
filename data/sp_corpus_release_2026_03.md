# SP corpus regeneration: UniProt 2026_03

This report records the corpus-definition and model-population changes shipped
with the release-pinned offline regeneration work in PR #74. It closes the
population audits requested by issues #67, #68, #69, #71, and #73.

## Reproducible source bundle

`fetch_sp_ground_truth.py --refresh` staged and validated one bundle before
publication:

- 55,516 current Vertebrata candidates from UniProt release `2026_03`
- 203 complete latest-version UniSave records for historical accessions absent
  from the current stream
- 1,405 exact taxonomy records covering every entry taxon and every
  classification root
- one manifest containing the release, exact candidate query and fields, row
  counts, byte sizes, and SHA-256 digests for all three artifacts

The candidate query was:

```text
(taxonomy_id:7742) AND (keyword:KW-0491 OR protein_name:histocompatibility OR protein_name:"class I" OR protein_name:"class II")
```

All raw scientific names, enrichment metadata, and taxonomy decisions are now
derived from that bundle. Ordinary raw, enrichment, and taxonomy regeneration
performs no network requests.

## Eligibility funnel

The 55,516 current and 203 archived rows form 55,719 disjoint candidates:

| Decision | Rows |
|---|---:|
| Labelled MHC alpha/beta chain with an exact `SIGNAL 1..N` boundary | 14,721 |
| Unlabelled, complete-N-terminus MHC alpha/beta inference set | 4,483 |
| Fragment | 26,600 |
| Outside the historical 100–500 aa length window | 4,017 |
| Non-MHC-chain query match | 2,872 |
| MHC-like but class/chain unresolved | 1,828 |
| Incomplete N terminus without a signal annotation | 1,198 |
| Invalid or fuzzy signal annotation | 0 |
| **Total** | **55,719** |

Missing `ft_signal` is never used as a negative label. Only complete
N-terminal rows with a usable MHC class/chain enter the separate unlabelled
inference set.

The broad query does contain the contamination described in #73. In
particular, the fixture gate now rejects an IL10RB cytokine-receptor record,
and the local-corpus label fallback was removed after it was shown to admit
invariant chains, class-II transactivators, and minor-histocompatibility
proteins.

## Taxonomy coverage

One shared definition now contains all 11 mutually exclusive source clades and
their category rules. The labelled rows are distributed as follows:

| Source clade | Rows | Category |
|---|---:|---|
| Mammalia | 9,749 | mammal subcategories |
| Aves | 784 | bird |
| Actinopterygii | 3,262 | fish |
| Lepidosauria | 454 | other vertebrate |
| Amphibia | 307 | other vertebrate |
| Chondrichthyes | 71 | fish |
| Testudines | 64 | other vertebrate |
| Crocodylia | 28 | other vertebrate |
| Coelacanthimorpha | 1 | fish |
| Dipnomorpha | 1 | fish |
| Myxini | 0 | fish |

This removes the historical 500-row Aves and Actinopterygii caps and adds the
turtle and crocodilian coverage requested by #71.

## Historical-label audit

Of 2,402 historical enriched rows, 2,238 remain in the new labelled corpus and
12,483 newly eligible rows were added. The 164 removals are explained:

- 151 are fragments (133 formerly gold, 12 curated, and 6 already excluded)
- 10 are already-curated non-MHC records
- 1 is the already-recorded unresolved class-II row `A0A8B9VS73`
- 2 formerly gold non-fragment rows are the non-MHC H60b and H60c proteins

Historical labels that remain source-backed but no longer classify from their
short current names are now explicit, versioned curation decisions. These
include archived accessions such as `A9UM13`, `Q31365`, and `Q8HWB2`, plus the
15 formerly gold non-fragment rows whose pinned short names are now generic or
omit the chain arm. The sole class/chain correction among retained rows is
`Q5TJH0`, changed from class-II alpha to class-II beta because its pinned
protein name and `HLA-DMB` gene both identify the DM beta chain. No gold row
silently degrades to unresolved in the regenerated labelled corpus.

## Model and benchmark population

The old enriched file contained 2,402 rows, but only 2,385 were eligible under
the evaluation policy; the other 16 were excluded non-MHC records and one was
unresolved. The regenerated boundary and sequence-cue models use exactly the
same 14,721 curated/gold rows as evaluation and mature-only control generation:

| Population | Before | After | Change |
|---|---:|---:|---:|
| Enriched rows | 2,402 | 14,721 | +12,319 |
| Curated/gold evaluation rows | 2,385 | 14,721 | +12,336 |
| Boundary-model training rows | 2,402 | 14,721 | +12,319 |
| Mature-only controls | 2,385 | 14,721 | +12,336 |
| Excluded/unresolved model rows | 17 | 0 | -17 |

With the regenerated models, the labelled benchmark parsed 13,571/14,721
(92.2%). Among parsed calls, 11,734/13,571 (86.5%) matched the exact signal
boundary and 12,500/13,571 (92.1%) were within two residues. On 14,721
mature-only controls, 14,034 returned the correct zero-SP result, 558
abstained, and 129 (0.88%) made a false-positive SP call.

Regression tests pin the funnel, curation strata, model population, complete
offline replay, and failed-refresh preservation behavior.
