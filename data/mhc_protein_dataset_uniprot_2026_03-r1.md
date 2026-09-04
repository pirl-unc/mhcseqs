# MHC protein dataset: `uniprot-2026_03-r1`

This is the first independently versioned `mhc-proteins` data release. It
preserves every current or historically retained row in the UniProt `2026_03`
Vertebrata MHC candidate bundle and adds normalized metadata plus mhcseqs
sequence-parser annotations. It is not the signal-peptide benchmark and is not
tied to the mhcseqs package version.

## Population

| Population | Records |
|---|---:|
| Current UniProtKB records | 55,516 |
| Historical latest-version UniSave records | 203 |
| **Total** | **55,719** |

The release-pinned source query is:

```text
(taxonomy_id:7742) AND (keyword:KW-0491 OR protein_name:histocompatibility OR protein_name:"class I" OR protein_name:"class II")
```

Because that query is intentionally broad, the dataset retains 4,851 records
classified as non-MHC query matches and 4,109 unresolved candidates alongside
46,759 included MHC records. Consumers selecting MHC identities should filter
on `inferred_disposition`; `parse_ok` describes structural decomposition and
is not an identity decision.

## Annotations

The 73-column records schema separates three kinds of information:

- Source columns retain UniProtKB, UniSave, and pinned taxonomy facts. Of the
  55,719 records, 20,053 have exact source signal boundaries, 4 have fuzzy
  source signal annotations, and 35,662 have no source signal annotation.
- Inferred columns normalize MHC class, chain, protein type, gene, label
  status, and record disposition from source metadata and explicit curation.
- Parsed columns contain mhcseqs signal-boundary and domain-part output,
  including groove halves, Ig/C-like support domain, tail, typed architecture,
  spans, scores, anchors, states, flags, and parse status.

The dataset build disables the parser's exact/lexical early shortcuts so the
parsed signal boundary is not copied from a shortcut learned on the same
source population. Source signal annotations remain available independently
in `source_signal_*`.

The parser produced 44,521 usable decompositions and positive parsed signal
boundaries for 24,779 records. Partial sequence parts are also preserved for
264 `short` results where `parse_ok` is false. All 55,719 source records remain
present even when identity inference or parsing abstains.

## Integrity

| Asset | Bytes | SHA-256 |
|---|---:|---|
| `mhc-proteins-uniprot-2026_03-r1.csv.gz` | 9,992,921 | `0ff706494d7d371e6bcf107c2e8ffc55bd47903befeaf8304ea7afa61c2b1830` |
| `mhc-proteins-uniprot-2026_03-r1.manifest.json` | 6,247 | `cf9fd6e6a058609628b63373a0505f9afaeaff8036333faaee89a176d200433d` |

Release QA confirmed:

- 55,719 unique accessions and no current/archive overlap
- no source sequence-length mismatches
- no parser exceptions
- no parsed signal-prefix/mature-sequence reconstruction mismatches
- no complete part/mature-sequence reconstruction mismatches
- deterministic sorted CSV and gzip output

The downloadable release also contains the current, deleted, and taxonomy
source artifacts, their manifest, and the exact label-curation file. Each is
pinned by byte count and SHA-256 in `mhcseqs/mhc_protein_datasets.json`.

## Install and reproduce

```bash
mhcseqs data install mhc-proteins --version uniprot-2026_03-r1

# Also install the complete offline-regeneration inputs.
mhcseqs data install mhc-proteins --version uniprot-2026_03-r1 --with-sources

python scripts/build_mhc_protein_dataset.py \
  --data-dir ~/.cache/mhcseqs/source-bundles/mhc-proteins/uniprot-2026_03-r1 \
  --label-curation ~/.cache/mhcseqs/source-bundles/mhc-proteins/uniprot-2026_03-r1/sp_ground_truth_label_curation.csv
```

The records manifest pins the source query and release, source artifact hashes,
schema, curation hash, learned-model hashes, generator version, coordinate
conventions, and output distributions.

## License

The redistributed UniProt-derived data is provided under
[CC BY 4.0](https://www.uniprot.org/help/license/) with attribution to the
UniProt Consortium. The mhcseqs code remains Apache-2.0 licensed.
