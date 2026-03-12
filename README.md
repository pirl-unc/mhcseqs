# mhcseqs

Self-contained pipeline for downloading, curating, and extracting binding
grooves from MHC (Major Histocompatibility Complex) protein sequences.

## Install

```bash
pip install mhcseqs
```

For development:

```bash
git clone https://github.com/openvax/mhcseqs.git
cd mhcseqs
./develop.sh          # uv pip install -e ".[dev]"
./test.sh             # pytest
./lint.sh             # ruff
```

## Quick start

### CLI

```bash
# Build all output CSVs from upstream databases
mhcseqs build

# Build to a specific directory
mhcseqs build --output-dir output/

# Look up a specific allele
mhcseqs lookup "HLA-A*02:01"

# Check version
mhcseqs --version
```

### Python API

```python
import mhcseqs

# Build the database (downloads FASTA sources, only needed once)
mhcseqs.build()

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

### Load all sequences as a DataFrame

```python
import mhcseqs, pandas as pd

mhcseqs.build()  # once — downloads + builds CSVs

# All binding grooves with structural decomposition (35K rows)
df = pd.read_csv("mhc-binding-grooves.csv")

# Or full protein sequences (one per two-field allele)
df = pd.read_csv("mhc-full-seqs.csv")

# Without pandas
rows = mhcseqs.load_grooves()     # list[dict]
rows = mhcseqs.load_full_seqs()   # list[dict]
```

## Output files

| File | Description |
|---|---|
| `mhc-seqs-raw.csv` | Every protein entry from both databases |
| `mhc-full-seqs.csv` | One representative mature protein per two-field allele group |
| `mhc-binding-grooves.csv` | Extracted binding groove + Ig domain + tail for each representative |

## Current data summary

Species category x MHC class counts (from `mhc-binding-grooves.csv`):

| Category | Class I | Class II | Total |
|---|---:|---:|---:|
| human | 17,462 | 7,878 | 25,340 |
| nhp | 4,639 | 2,466 | 7,105 |
| murine | 49 | 9 | 58 |
| other_mammal | 807 | 1,564 | 2,371 |
| bird | 28 | 0 | 28 |
| fish | 90 | 85 | 175 |
| **total** | **23,075** | **12,002** | **35,077** |

Groove parse success rate: 99.6% (146 failures out of 35,077 entries).

## Structural decomposition

Each protein chain is decomposed into four contiguous regions:

| Column | Class I alpha | Class II alpha | Class II beta |
|---|---|---|---|
| `groove1` | α1 domain (~90 aa) | α1 domain (~83 aa) | — |
| `groove2` | α2 domain (~93 aa) | — | β1 domain (~93 aa) |
| `ig_domain` | α3 Ig-fold (~95 aa) | α2 Ig-fold (~95 aa) | β2 Ig-fold (~95 aa) |
| `tail` | TM + cytoplasmic | TM + cytoplasmic | TM + cytoplasmic |

For a class I chain: `mature_protein = groove1 + groove2 + ig_domain + tail`

## Key columns

All three CSVs share: `gene`, `mhc_class`, `chain`, `species`,
`species_category`, `species_prefix`, `source`.

`species_category` is one of: `human`, `nhp`, `murine`, `other_mammal`,
`bird`, `fish`, `other_vertebrate`.

## Data sources

| Source | Species | URL |
|---|---|---|
| IMGT/HLA | Human | `https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/hla_prot.fasta` |
| IPD-MHC | Non-human | `https://raw.githubusercontent.com/ANHIG/IPDMHC/Latest/MHC_prot.fasta` |
| UniProt | B2M references | Curated (shipped in `mhcseqs/b2m_sequences.csv`) |

## Species prefixes

| Species | Latin name | MHC prefix |
|---|---|---|
| human | *Homo sapiens* | HLA |
| macaque | *Macaca mulatta* | Mamu |
| chimpanzee | *Pan troglodytes* | Patr |
| gorilla | *Gorilla gorilla* | Gogo |
| mouse | *Mus musculus* | H2 |
| rat | *Rattus norvegicus* | RT1 |
| cattle | *Bos taurus* | BoLA |
| pig | *Sus scrofa* | SLA |
| horse | *Equus caballus* | ELA |
| sheep | *Ovis aries* | OLA |
| dog | *Canis lupus familiaris* | DLA |
| cat | *Felis catus* | FLA |
| chicken | *Gallus gallus* | Gaga |
| salmon | *Salmo salar* | Sasa |
| zebrafish | *Danio rerio* | Dare |

## Groove extraction algorithm

The groove parser is **alignment-free** — it uses conserved Cys-Cys disulfide
pairs in Ig-fold domains as structural landmarks to slice domain boundaries
without multiple sequence alignment.

See [groove.py](mhcseqs/groove.py) module docstring for detailed algorithm
documentation with ASCII structural diagrams.

## Dependencies

- **Python 3.10+**
- **[mhcgnomes](https://github.com/openvax/mhcgnomes)** — allele name parsing

No alignment tools, BLAST, or structure databases are required.

## Repository structure

```
mhcseqs/
├── mhcseqs/
│   ├── __init__.py        # Public API
│   ├── __main__.py        # CLI entry point
│   ├── version.py         # Package version
│   ├── download.py        # FASTA source downloading
│   ├── species.py         # Species taxonomy (29-class → 7-class)
│   ├── alleles.py         # Allele name parsing (mhcgnomes wrapper)
│   ├── groove.py          # Binding groove extraction + mutation support
│   ├── imgt.py            # IMGT G-DOMAIN position numbering
│   ├── pipeline.py        # Three-step build pipeline
│   ├── validate.py        # Post-build validation
│   └── b2m_sequences.csv  # Reference B2M sequences
├── tests/                 # pytest test suite
├── data/
├── build.py               # Convenience shim
├── pyproject.toml         # Package metadata
├── develop.sh             # Install in dev mode
├── lint.sh                # Run ruff
├── test.sh                # Run pytest
└── deploy.sh              # Build + publish to PyPI
```

## License

Apache 2.0
