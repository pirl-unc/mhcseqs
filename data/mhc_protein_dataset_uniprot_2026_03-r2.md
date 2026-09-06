# MHC protein dataset: `uniprot-2026_03-r2`

This revision corrects identity labels while preserving all 55,719 records,
sequences, and source annotations from the pinned UniProt `2026_03` bundle:
55,516 current UniProtKB entries and 203 historical UniSave entries. The
73-column schema and source/curation assets are unchanged. Historical r1
assets remain available and are never overwritten.

## Identity corrections

- Known non-MHC genes are checked before permissive protein-name heuristics.
  For example, CIITA is a transactivator, not an MHC class II alpha chain.
  See [UniProt P33076](https://www.uniprot.org/uniprotkb/P33076/entry).
- Missing or vague names stay unresolved instead of being called non-MHC.
  This preserves the uncertainty of A0A0G2JMH6 and Q29967, whose pinned
  sequences exactly match included HLA-DRA and HLA-DQB1 records.
- Mixed MHC/helper-gene synonym lists also stay unresolved; their ordering
  does not decide identity. Explicit accession curation remains authoritative.

| Identity disposition | Records |
|---|---:|
| `include` | 46,019 |
| `exclude_non_mhc` | 1,510 |
| `retain_unresolved` | 8,190 |

`parse_ok` remains a structural decomposition result, not an identity label.
All source rows remain present, including excluded and unresolved candidates.
Signal annotations are optional: 20,053 records have exact source boundaries,
4 have fuzzy source features, and 35,662 have no source signal annotation.

## Install and reproduce

The [data release](https://github.com/pirl-unc/mhcseqs/releases/tag/data-mhc-proteins-uniprot-2026_03-r2)
contains records, their manifest, and all five pinned regeneration inputs.
Every asset's byte count and SHA-256 is pinned in
[`mhc_protein_datasets.json`](../mhcseqs/mhc_protein_datasets.json).

Full-population QA verified 55,719 unique accessions, unchanged source facts
and sequences, no parser exceptions, and exact boundary/part reconstructions.
There are no included records with `non_mhc_gene` parser flags and no identical
taxon/sequence pairs classified as both included and excluded. The parser
produced 44,530 usable decompositions; all other source records are retained.

```bash
mhcseqs data install mhc-proteins --version uniprot-2026_03-r2 --with-sources
mhcseqs data path mhc-proteins --version uniprot-2026_03-r2

python scripts/build_mhc_protein_dataset.py \
  --data-dir ~/.cache/mhcseqs/source-bundles/mhc-proteins/uniprot-2026_03-r2 \
  --label-curation ~/.cache/mhcseqs/source-bundles/mhc-proteins/uniprot-2026_03-r2/sp_ground_truth_label_curation.csv \
  --revision 2
```

Exact regeneration uses the mhcseqs 2.6.12 source and `mhcgnomes==3.41.0`.
The manifest also pins the learned models, curation, schema, and source hashes.
Parser exact/lexical signal-boundary shortcuts remain disabled.

Generated records and their manifest are staged and validated before one
bundle directory is published. Custom output files must share a dedicated
directory with no unrelated files. Failed regeneration preserves the previous
pair; interrupted installation/build swaps recover before requiring network
access. Cache readers use a verified compressed snapshot so concurrent forced
reinstalls cannot invalidate a read.

## License

UniProt-derived data is provided under
[CC BY 4.0](https://www.uniprot.org/help/license/) with attribution to the
UniProt Consortium. The mhcseqs code remains Apache-2.0 licensed.
