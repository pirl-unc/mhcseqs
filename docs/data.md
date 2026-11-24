# Datasets and columns

This page describes what a build contains, what each exported column means,
and the separately versioned `mhc-proteins` record layer. Read it before
filtering rows or comparing the two datasets, which are deliberately not the
same population.

## Current data summary

All sources (IMGT/HLA, IPD-MHC, and curated UniProt/GenBank references)
are merged into a single dataset. Representatives by species category and chain type:

| Category | Class I | Class II | Other | Total |
|---|---:|---:|---:|---:|
| human | 18,040 | 8,196 | 25 | 26,261 |
| nhp | 4,475 | 2,414 | 0 | 6,889 |
| murine | 964 | 566 | 73 | 1,603 |
| ungulate | 650 | 1,132 | 6 | 1,788 |
| carnivore | 174 | 330 | 2 | 506 |
| other_mammal | 2,312 | 1,265 | 248 | 3,825 |
| bird | 5,974 | 3,337 | 163 | 9,474 |
| fish | 1,264 | 3,594 | 72 | 4,930 |
| other_vertebrate | 471 | 665 | 28 | 1,164 |
| **total** | **34,324** | **21,499** | **617** | **56,440** |

Covering 558+ species. Groove parse success rate on IMGT/IPD-MHC
entries: 99.3%.

## Key columns

| Column | Description |
|---|---|
| `two_field_allele` | Allele name at two-field resolution |
| `gene` | MHC gene (e.g., A, DRB1, BF, UA) |
| `mhc_class` | I or II |
| `chain` | alpha, beta, or B2M |
| `species` | Latin binomial from source |
| `species_category` | One of the 9 [species categories](#species-categories) |
| `source` | `imgt`, `ipd_mhc`, `uniprot`, or `genbank` |
| `source_id` | Database accession for provenance |
| `groove_status` | See [groove status values](parsing.md#groove-status-values) |
| `is_functional` | True if groove parsed and not null/pseudogene |

## All columns

`mhcseqs build` writes two CSVs. `mhc-seqs-raw.csv` has one row per source
sequence record; `mhc-full-seqs.csv` has one representative row per two-field
allele, with the structural decomposition added.

### `mhc-full-seqs.csv`

| Column | Description |
|---|---|
| `two_field_allele` | Allele name at two-field resolution |
| `representative_allele` | Full-resolution allele chosen as the representative |
| `protein_seq_selection` | How the representative was chosen (`unique`, `nested_longest`, ...) |
| `gene`, `mhc_class`, `chain`, `species`, `species_category`, `species_prefix`, `source`, `source_id` | As in [key columns](#key-columns); `species_prefix` is the full-binomial prefix, e.g. `HomoSapiens` |
| `seq_len`, `sequence` | Length and sequence of the full protein |
| `mature_start` | Offset where the mature protein begins (signal peptide length) |
| `mature_sequence` | Protein with the signal peptide removed |
| `groove1`, `groove2` | The two groove halves (see [structural decomposition](parsing.md)) |
| `groove_seq` | `groove1` followed by `groove2` |
| `groove1_len`, `groove2_len` | Lengths of the groove halves |
| `ig_domain`, `ig_domain_len` | Ig-like support domain and its length |
| `tail`, `tail_len` | Linker, transmembrane and cytoplasmic region, and its length |
| `domain_architecture` | Typed domain grammar, e.g. `signal_peptide>g_alpha1>g_alpha2>c1_alpha3>tail_linker>transmembrane>cytoplasmic_tail` |
| `domain_spans` | 1-based domain coordinates in `sequence`, e.g. `signal_peptide:1-24;g_alpha1:25-114;...` |
| `groove_status` | Parse outcome; see [groove status values](parsing.md#groove-status-values) |
| `groove_flags` | Comma-separated parser flags (empty for a clean parse) |
| `anchor_type` | Which Cys pair anchored the parse (`alpha2_cys`, `alpha3_cys`, `beta1_cys`, `beta2_cys`, `raw_fragment`) |
| `is_null`, `is_questionable`, `is_pseudogene` | Allele status flags |
| `is_functional` | True if groove parsed and not null/pseudogene |

### `mhc-seqs-raw.csv`

| Column | Description |
|---|---|
| `allele_raw` | Allele token or accession from the source record |
| `allele_normalized` | Full-resolution name normalized by mhcgnomes |
| `two_field_allele` | Truncated to two-field resolution |
| `gene`, `mhc_class`, `chain`, `species`, `species_category`, `species_prefix`, `source`, `source_id` | As in [key columns](#key-columns) |
| `seq_len`, `sequence` | Length and sequence of the full protein |
| `has_signal_peptide` | True when the sequence starts with Met and the parsed mature start is plausible |
| `signal_peptide_len` | Mature start offset from the groove parser |
| `signal_peptide_seq` | Residues before the mature start |
| `is_null`, `is_questionable`, `is_pseudogene` | Allele status flags |

## Species categories

`mhcseqs.normalize_species()` assigns a fine-grained class and
`mhcseqs.normalize_mhc_species()` maps it to one of nine `species_category`
values:

| Category | Fine-grained classes |
|---|---|
| `human` | human |
| `nhp` | macaque, chimpanzee, gorilla, orangutan, baboon, other_nhp |
| `murine` | mouse, rat |
| `ungulate` | cattle, pig, horse, sheep, goat |
| `carnivore` | dog, cat |
| `other_mammal` | rabbit, bat, marsupial, monotreme, cetacean, other_mammal |
| `bird` | chicken, other_bird |
| `fish` | salmon, zebrafish, shark, other_fish |
| `other_vertebrate` | reptile, amphibian, other_vertebrate |

## Species prefixes

mhcseqs emits full-binomial species aliases such as `HomoSapiens-A*02:01`
and `MusMusculus-K*b`. Its input registry accepts 476 current, historical, and
external-database prefix assignments—including all 137 designation tokens in
the current 125-organism IPD-MHC taxonomy register—and records an evidence URL
for every one. Colliding short codes require explicit species context. mhcseqs
does not invent abbreviated species prefixes, and mechanically generated
2+2/4+4/5+5 aliases are rejected unless that exact spelling has external
evidence. See the [prefix audit](https://github.com/pirl-unc/mhcseqs/blob/main/data/species_prefix_audit_20260831.md).
The audit includes a versioned UniProtKB/UniSave provenance snapshot for all
239 records behind the 21 historical pairs that previously lacked raw inputs.

## Full UniProt protein records

The independently versioned `mhc-proteins` dataset preserves all 55,719
current and historical records from the release-pinned UniProt Vertebrata MHC
candidate query. It is a source-complete record layer, distinct from both the
56,440 representative IMGT/IPD build summarized above and the smaller
signal-peptide benchmark derived from it.

Install the default data version (`uniprot-2026_03-r2`) with:

```bash
mhcseqs data install mhc-proteins
mhcseqs data path mhc-proteins
```

Code releases and data releases are deliberately independent. Pin a data
version in automated work:

```bash
mhcseqs data install mhc-proteins --version uniprot-2026_03-r2
```

See [Python API](python-api.md#full-uniprot-protein-records) for streaming and
eager loading.

### Column provenance

The columns are intentionally separated by provenance:

- `source_*` plus organism, sequence, lineage, name, and gene columns are
  UniProtKB, UniSave, or taxonomy-snapshot facts. An empty source annotation
  means "not annotated," never "biologically absent."
- `inferred_*` columns contain normalized MHC class, chain, type, gene, and
  candidate-disposition decisions derived from source metadata and explicit
  accession curation.
- `parsed_*`, `mature_sequence`, `groove1`, `groove2`, `ig_domain`, `tail`,
  boundary, span, score, and state columns are outputs from the mhcseqs
  sequence parser. Partial parts can remain populated when `parse_ok` is false;
  inspect `parse_status` before treating them as a complete decomposition.

The broad source query deliberately retains contaminants and ambiguous hits so
that the source population is reproducible. `inferred_disposition` is the MHC
identity decision; `parse_ok` only reports whether the structural parser
produced a usable decomposition and is not an independent identity label.
Known non-MHC genes are excluded; insufficient or conflicting identity
evidence stays `retain_unresolved`. The
[r2 release notes](https://github.com/pirl-unc/mhcseqs/blob/main/data/mhc_protein_dataset_uniprot_2026_03-r2.md)
describe the corrected labels. Historical `r1` remains available unchanged.
These full-dataset identity rules are separate from the stored SP benchmark's
selection policy, so they do not silently change benchmark membership or prevent
its offline enrichment.

### Reproducing the artifact offline

To reproduce the records artifact offline, install the source bundle too:

```bash
mhcseqs data install mhc-proteins --version uniprot-2026_03-r2 --with-sources
python scripts/build_mhc_protein_dataset.py \
  --data-dir ~/.cache/mhcseqs/source-bundles/mhc-proteins/uniprot-2026_03-r2 \
  --label-curation ~/.cache/mhcseqs/source-bundles/mhc-proteins/uniprot-2026_03-r2/sp_ground_truth_label_curation.csv \
  --revision 2
```

The builder publishes records and their manifest together. By default they go
under `<source-root>-generated/mhc-proteins/<data-version>/`, a sibling tree
outside the replaceable source bundle. Custom `--output` and `--manifest-output`
paths must be distinct files in the same dedicated directory outside the source
tree, containing no unrelated files. Source reinstallation therefore cannot
delete generated outputs or their staging directory. Failed builds preserve the previous
pair, and interrupted swaps are recovered on the next access. To reproduce a
release exactly, use its generator and mhcgnomes versions recorded in the
manifest as well as its pinned source bundle.

The release manifest pins source queries, releases, byte counts, SHA-256
digests, model hashes, schema, and output distributions. UniProt-derived data
is redistributed under [CC BY 4.0](https://www.uniprot.org/help/license/) with
attribution to the UniProt Consortium.
