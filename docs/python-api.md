# Python API

This page covers the functions mhcseqs exports: looking up one allele and
reading its `AlleleRecord`, loading the whole build as a table, the allele and
species helpers, streaming the `mhc-proteins` records, and which module each
lower-level function lives in.

## Look up one allele

```python
import mhcseqs

# Build the database (downloads to ~/.cache/mhcseqs/, only needed once)
paths = mhcseqs.build()  # BuildPaths dataclass

# Look up any allele -> AlleleRecord with everything
r = mhcseqs.lookup("HLA-A*02:01")
r.allele            # "HomoSapiens-A*02:01" (two-field)
r.full_allele       # "HomoSapiens-A*02:01:01:100" (full resolution)
r.sequence          # full protein (with signal peptide)
r.mature_sequence   # signal peptide removed (computed property)
r.mature_start      # signal peptide length (24 for HLA-A*02:01)
r.groove1           # α1 domain
r.groove2           # α2 domain
r.ig_domain         # α3 Ig-fold
r.tail              # TM + cytoplasmic
r.domains           # typed domain spans
r.domain_architecture
r.domain_spans
r.species_category  # "human"

# Apply mutations (IEDB-style, e.g. "K66A")
m = mhcseqs.lookup("HLA-A*02:01", mutations=["K66A", "D77S"])
```

`r.allele` and `r.full_allele` carry the full-binomial species prefix that
mhcseqs emits (`HomoSapiens-A*02:01`), not the short IMGT form. See
[datasets and columns](data.md#species-prefixes) for why.

## Load the whole build

```python
import mhcseqs

# As a DataFrame (full sequence + groove decomposition + metadata)
df = mhcseqs.load_sequences_dataframe()

# Or as a list of dicts (no pandas dependency)
rows = mhcseqs.load_sequences_dict()
```

## Allele and species utilities

```python
>>> import mhcseqs
>>> mhcseqs.normalize_allele_name("A0201")
'HomoSapiens-A*02:01'
>>> mhcseqs.infer_gene("HLA-A*02:01")
'A'
>>> mhcseqs.infer_mhc_class("HLA-A*02:01")
'I'
>>> mhcseqs.normalize_species("Macaca mulatta")
'macaque'
>>> mhcseqs.normalize_mhc_species("macaque")
'nhp'
>>> mhcseqs.get_latin_name("macaque")
'Macaca mulatta'
>>> mhcseqs.get_canonical_prefix("macaque")
'MacacaMulatta'
```

`normalize_species` gives the fine-grained class and `normalize_mhc_species`
the coarse category used in the `species_category` column.

## Full UniProt protein records

Records from the independently versioned `mhc-proteins` dataset can be
streamed without pandas or loaded eagerly. The requested data version is
downloaded once, checksum-verified, and then reused from the local cache.
Cache readers capture a verified compressed snapshot (about 10 MB), releasing
the version lock before parsing or yielding rows:

```python
import mhcseqs

rows = mhcseqs.iter_mhc_protein_records(version="uniprot-2026_03-r2")
first = next(rows)

all_rows = mhcseqs.load_mhc_protein_records(version="uniprot-2026_03-r2")
df = mhcseqs.load_mhc_protein_dataframe(version="uniprot-2026_03-r2")
```

For read-only/shared deployments, preinstall the cache and retain its sibling
`.<version>.lock` files. Reads need only read access to those files; creating
an installation or recovering an interrupted swap still requires write access.
Listing a missing version requires no writes and reports it as not installed.

## Module reference

### `mhcseqs.domain_parsing`

Core groove extraction.

- `extract_groove(seq, mhc_class, ...)` — main dispatcher
- `parse_class_i(seq, ...)` — class I alpha chain parser
- `parse_class_ii_alpha(seq, ...)` — class II alpha chain parser
- `parse_class_ii_beta(seq, ...)` — class II beta chain parser
- `find_cys_pairs(seq)` — locate Ig-fold Cys-Cys pairs
- `AlleleRecord` — dataclass with all parsed fields

### `mhcseqs.alleles`

Allele name parsing via mhcgnomes.

- `parse_allele_name(name)` — parse with mhcgnomes
- `normalize_allele_name(name)` — canonical two-field form
- `infer_gene(allele)` — extract gene (e.g. `"A"`, `"DRB1"`)
- `infer_mhc_class(allele)` — `"I"` or `"II"`

### `mhcseqs.species`

Species normalization.

- `normalize_species(raw)` — fine-grained class
- `normalize_mhc_species(raw)` — coarse category
- `get_latin_name(raw)` — canonical Latin name
- `get_canonical_prefix(raw)` — MHC naming prefix
- `LATIN_NAMES`, `CANONICAL_MHC_PREFIXES` — lookup dicts

### `mhcseqs.pipeline`

Two-step build pipeline.

- `build_raw_index(fasta_paths, out_csv)` — step 1: parse FASTA into the raw CSV
- `build_full_seqs(raw_csv, out_csv)` — step 2: select representatives and extract grooves

### `mhcseqs.validate`

Post-build sanity checks.

- `validate_build(raw_csv, full_csv)` — run all checks
- `format_validation_report(warnings, stats)` — human-readable report

### `mhcseqs.download`

FASTA source downloading.

- `download_all(dest_dir)` — download IMGT/HLA and IPD-MHC
- `SOURCES` — dict of source URLs and metadata
