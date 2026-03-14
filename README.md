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

Pre-built CSVs are also attached to each
[GitHub release](https://github.com/openvax/mhcseqs/releases).

## Output files

| File | Description |
|---|---|
| `mhc-seqs-raw.csv` | Every protein entry from all sources |
| `mhc-full-seqs.csv` | One representative per two-field allele: full sequence, groove decomposition, and metadata |

## Current data summary

### Built sequences (from `mhc-full-seqs.csv`)

Species category x MHC class counts from IMGT/HLA and IPD-MHC:

| Category | Class I | Class II | Total |
|---|---:|---:|---:|
| human | 17,462 | 7,878 | 25,340 |
| nhp | 4,639 | 2,486 | 7,125 |
| murine | 59 | 29 | 88 |
| ungulate | 638 | 1,128 | 1,766 |
| carnivore | 166 | 318 | 484 |
| cetacean | 3 | 98 | 101 |
| bird | 28 | 0 | 28 |
| fish | 90 | 85 | 175 |
| **total** | **23,085** | **12,022** | **35,107** |

Groove parse success rate: 99.6% (146 failures out of 35,107 entries).

### Diverse MHC sequences (shipped with package)

16,208 curated UniProt sequences from underrepresented taxonomic groups,
fetched with `scripts/fetch_diverse_mhc.py` and curated with
`scripts/curate_diverse_mhc.py`:

| Source group | Sequences | Class I | Class II |
|---|---:|---:|---:|
| bird (non-chicken) | 8,869 | 5,564 | 3,305 |
| bony fish | 4,747 | 1,313 | 3,434 |
| amphibian | 847 | 305 | 542 |
| chicken | 589 | 274 | 315 |
| marsupial | 579 | 432 | 147 |
| bat | 181 | 22 | 159 |
| reptile (lepidosauria) | 144 | 93 | 51 |
| reptile (crocodylia) | 150 | 133 | 17 |
| shark/ray | 63 | 52 | 11 |
| reptile (testudines) | 28 | 14 | 14 |
| monotreme | 11 | 5 | 6 |
| **total** | **16,208** | **8,207** | **8,001** |

Covering 614 unique species prefixes across 11 taxonomic groups.

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

Both CSVs share: `gene`, `mhc_class`, `chain`, `species`,
`species_category`, `species_prefix`, `source`, `source_id`.

`source` is one of: `imgt`, `ipd_mhc`, `uniprot_curated`, `uniprot_reference`.

`source_id` is the database accession for provenance tracking (e.g.,
`HLA00001` for IMGT, `NHP00001` for IPD-MHC, `P01901` for UniProt).

`species_category` is one of: `human`, `nhp`, `murine`, `ungulate`,
`carnivore`, `cetacean`, `other_mammal`, `bird`, `fish`, `other_vertebrate`.

## Data sources

| Source | `source` value | Species | Data |
|---|---|---|---|
| IMGT/HLA | `imgt` | Human | Downloaded at build time |
| IPD-MHC | `ipd_mhc` | Non-human | Downloaded at build time |
| UniProt | `uniprot_reference` | Multi-species | B2M references (shipped in package) |
| UniProt | `uniprot_curated` | Mouse | 30 H-2 alleles (shipped in package) |
| UniProt | `uniprot_diverse` | 614 species | 16,208 diverse MHC sequences (shipped in package) |

IMGT/HLA and IPD-MHC FASTA files are downloaded on first `build` and cached.
The UniProt curated CSVs (`b2m_sequences.csv`, `mouse_h2_sequences.csv`,
`diverse_mhc_sequences.csv`) ship inside the `mhcseqs` package — no download needed.

To refresh the diverse MHC dataset from UniProt:

```bash
python scripts/fetch_diverse_mhc.py    # Download raw data → data/diverse_mhc_raw.csv
python scripts/curate_diverse_mhc.py   # Curate → mhcseqs/diverse_mhc_sequences.csv
```

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
- **[mhcgnomes](https://github.com/openvax/mhcgnomes) >= 3.1.0** — allele name parsing

No alignment tools, BLAST, or structure databases are required.

## Repository structure

```
mhcseqs/
├── mhcseqs/
│   ├── __init__.py        # Public API
│   ├── __main__.py        # CLI entry point
│   ├── version.py         # Package version
│   ├── download.py        # FASTA source downloading
│   ├── species.py         # Species taxonomy (29-class → 10-class)
│   ├── alleles.py         # Allele name parsing (mhcgnomes wrapper)
│   ├── groove.py          # Binding groove extraction + mutation support
│   ├── imgt.py            # IMGT G-DOMAIN position numbering
│   ├── pipeline.py        # Two-step build pipeline
│   ├── validate.py        # Post-build validation
│   ├── b2m_sequences.csv                # Reference B2M sequences (UniProt)
│   ├── mouse_h2_sequences.csv          # Mouse H-2 sequences (UniProt)
│   └── diverse_mhc_sequences.csv      # 16k diverse MHC sequences (UniProt)
├── tests/                 # pytest test suite
├── scripts/
│   ├── fetch_diverse_mhc.py      # Download diverse MHC from UniProt
│   ├── curate_diverse_mhc.py     # Curate into shipped CSV
│   └── validate_signal_peptides.py
├── data/                          # Intermediate data (not shipped)
├── build.py               # Convenience shim
├── pyproject.toml         # Package metadata
├── develop.sh             # Install in dev mode
├── lint.sh                # Run ruff
├── test.sh                # Run pytest
└── deploy.sh              # Build + publish to PyPI
```

## License

Apache 2.0
