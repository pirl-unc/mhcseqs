# mhcseqs

Self-contained pipeline for downloading, curating, and extracting binding
grooves from MHC (Major Histocompatibility Complex) protein sequences.

## Install

```bash
pip install mhcseqs
```

For development:

```bash
git clone https://github.com/pirl-unc/mhcseqs.git
cd mhcseqs
./develop.sh          # uv pip install -e ".[dev]"
./test.sh             # pytest
./lint.sh             # ruff
```

## Quick start

### CLI

```bash
# Build all output CSVs (writes to ~/.cache/mhcseqs/)
mhcseqs build

# Build to a specific directory instead
mhcseqs build --output-dir output/

# Look up a specific allele
mhcseqs lookup "HLA-A*02:01"

# Check version
mhcseqs --version
```

### Python API

```python
import mhcseqs

# Build the database (downloads to ~/.cache/mhcseqs/, only needed once)
paths = mhcseqs.build()  # BuildPaths dataclass

# Look up any allele → AlleleRecord with everything
r = mhcseqs.lookup("HLA-A*02:01")
r.sequence          # full protein (with signal peptide)
r.mature_sequence   # signal peptide removed (computed property)
r.mature_start      # signal peptide length (24 for HLA-A*02:01)
r.groove1           # α1 domain
r.groove2           # α2 domain
r.ig_domain         # α3 Ig-fold
r.tail              # TM + cytoplasmic
r.species_category  # "human"

# Apply mutations (IEDB-style, e.g. "K66A")
m = mhcseqs.lookup("HLA-A*02:01", mutations=["K66A", "D77S"])
```

### Load as a DataFrame

```python
import mhcseqs

# As a DataFrame (full sequence + groove decomposition + metadata)
df = mhcseqs.load_sequences_dataframe()

# Or as a list of dicts (no pandas dependency)
rows = mhcseqs.load_sequences_dict()
```

## Current data summary

Three sources are merged into a single dataset:

| Source | Entries | Species | Notes |
|---|---:|---|---|
| IMGT/HLA | 44,630 | Human | Downloaded at build time |
| IPD-MHC | 12,380 | Non-human mammals, birds, fish | Downloaded at build time |
| UniProt | 20,566 | 500+ species | Curated diverse MHC, B2M, H-2 references (shipped in package) |
| **Total raw** | **77,576** | | |
| **After merge/dedup** | **55,658** | | One representative per two-field allele |
| **Groove OK** | **54,121** | | 97.2% of representatives |

### By species category

| Category | Count |
|---|---:|
| human | 25,364 |
| bird | 9,312 |
| nhp | 7,125 |
| fish | 4,859 |
| ungulate | 1,768 |
| murine | 1,587 |
| other_vertebrate | 1,137 |
| other_mammal | 943 |
| carnivore | 484 |

Species categories: human, nhp (non-human primates), murine (mice, rats, rodents),
ungulate (cattle, pig, horse, sheep, goat), carnivore (dog, cat),
other_mammal (marsupials, monotremes, bats, cetaceans, rabbit),
bird, fish, other_vertebrate (reptiles, amphibians).

## Data directory

By default, `mhcseqs build` downloads FASTA files and writes output CSVs to
`~/.cache/mhcseqs/` (override with `$MHCSEQS_DATA` or `--output-dir`).

```
~/.cache/mhcseqs/
├── fasta/                     # Downloaded FASTA source files
│   ├── hla_prot.fasta         # IMGT/HLA (human)
│   └── ipd_mhc_prot.fasta    # IPD-MHC (non-human)
├── mhc-seqs-raw.csv           # Every protein entry from all sources
├── mhc-full-seqs.csv          # One representative per two-field allele (with grooves)
├── mhc-merge-report.txt       # Deduplication decisions
└── mhc-validation-report.txt  # Sanity checks
```

## Structural decomposition

Each protein chain is decomposed into four contiguous regions:

| Column | Class I alpha | Class II alpha | Class II beta |
|---|---|---|---|
| `groove1` | α1 domain (~90 aa) | α1 domain (~83 aa) | — |
| `groove2` | α2 domain (~93 aa) | — | β1 domain (~93 aa) |
| `ig_domain` | α3 Ig-fold (~95 aa) | α2 Ig-fold (~95 aa) | β2 Ig-fold (~95 aa) |
| `tail` | TM + cytoplasmic | TM + cytoplasmic | TM + cytoplasmic |

For a class I chain: `mature_protein = groove1 + groove2 + ig_domain + tail`

## Groove extraction algorithm

The groove parser is **alignment-free** — it uses conserved Cys-Cys disulfide
pairs in Ig-fold domains as structural landmarks to slice domain boundaries
without multiple sequence alignment.

### Signal peptide inference

Signal peptide length is inferred from the Cys pair position, not from
sequence motifs. The conserved Ig-fold Cys has a known position in the
mature protein (e.g., position 100 for class I α2). The offset between the
raw sequence position and the expected mature position gives the signal
peptide length:

```
mature_start = raw_cys_position - expected_mature_cys_position
```

Gene-specific constants account for groove domain length variation:
DQA (109), DMA (120), DPB (114). All others use the defaults (class I: 100,
class II alpha: 106, class II beta: 116).

### Groove status values

| Status | Meaning |
|---|---|
| `ok` | Full decomposition: groove1 + groove2 + ig_domain + tail |
| `alpha1_only` | Single-exon class I fragment — α1 domain (no Cys pair) |
| `alpha2_only` | Single-exon class I fragment — α2 domain (has Cys pair) |
| `beta1_only_fallback` | Class II beta with β1 pair only (no β2 Ig pair) |
| `fragment_fallback` | Short fragment used as raw groove sequence |
| `inferred_from_alpha3` | Groove boundaries estimated from α3 Cys pair |
| `not_applicable` | Non-groove gene (B2M, MICA, MICB, HFE, MR1) |
| `non_classical` | Non-classical MHC lineage (fish L/S/P/H) |
| `short` | Groove half < 70 aa — unlikely to be functional |
| `suspect_anchor` | Cys mutation produced implausible mature_start |

See [groove.py](mhcseqs/groove.py) module docstring for detailed algorithm
documentation with ASCII structural diagrams.

## Key columns

| Column | Description |
|---|---|
| `two_field_allele` | Allele name at two-field resolution |
| `gene` | MHC gene (e.g., A, DRB1, BF, UA) |
| `mhc_class` | I or II |
| `chain` | alpha, beta, or B2M |
| `species` | Latin binomial from source |
| `species_category` | One of 9 categories above |
| `source` | `imgt`, `ipd_mhc`, or `uniprot` |
| `source_id` | Database accession for provenance |
| `groove_status` | See table above |
| `is_functional` | True if groove parsed and not null/pseudogene |

## Dependencies

- **Python 3.10+**
- **[mhcgnomes](https://github.com/pirl-unc/mhcgnomes) >= 3.14.0** — allele name parsing with species-directed disambiguation

No alignment tools, BLAST, or structure databases are required.

## License

Apache 2.0
