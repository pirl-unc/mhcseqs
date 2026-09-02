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

# Inspect cached downloads and built CSVs (path, size, age)
mhcseqs data list

# Re-download source FASTAs (IMGT/HLA + IPD-MHC publish from a rolling "Latest")
mhcseqs data refresh
mhcseqs build --force-download   # equivalently, rebuild with fresh sources

# Delete cached files
mhcseqs data clear               # source FASTAs
mhcseqs data clear --built       # also remove built CSVs/reports
mhcseqs data clear --built-only  # only built CSVs/reports

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
r.domains           # typed domain spans
r.domain_architecture
r.domain_spans
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

All sources (IMGT/HLA, IPD-MHC, and curated UniProt/GenBank references)
are merged into a single dataset. Categorized class I/II representatives:

| Category | Class I | Class II | Total |
|---|---:|---:|---:|
| human | 18,041 | 8,196 | 26,237 |
| nhp | 4,475 | 2,414 | 6,889 |
| murine | 967 | 566 | 1,533 |
| ungulate | 654 | 1,132 | 1,786 |
| carnivore | 176 | 330 | 506 |
| other_mammal | 549 | 449 | 998 |
| bird | 5,975 | 3,337 | 9,312 |
| fish | 1,265 | 3,594 | 4,859 |
| other_vertebrate | 471 | 665 | 1,136 |
| **total** | **32,573** | **20,683** | **53,256** |

Covering 558+ species. Groove parse success rate on IMGT/IPD-MHC
entries: 99.3%.

## Structural decomposition

The parser materializes an explicit domain grammar:

| Chain | Grammar |
|---|---|
| Class I alpha | `signal_peptide? -> g_alpha1 -> g_alpha2 -> c1_alpha3 -> transmembrane? -> cytoplasmic_tail?` |
| Class II alpha | `signal_peptide? -> g_alpha1 -> c1_alpha2 -> transmembrane? -> cytoplasmic_tail?` |
| Class II beta | `signal_peptide? -> g_beta1 -> c1_beta2 -> transmembrane? -> cytoplasmic_tail?` |

The exported contiguous sequence fields are:

| Column | Class I alpha | Class II alpha | Class II beta |
|---|---|---|---|
| `groove1` | α1 domain (~80-95 aa typical) | α1 domain (~75-95 aa typical) | — |
| `groove2` | α2 domain (~80-100 aa typical) | — | β1 domain (~70-100 aa typical) |
| `ig_domain` | α3 C-like support domain | α2 C-like support domain | β2 C-like support domain |
| `tail` | linker + TM + cytoplasmic tail | linker + TM + cytoplasmic tail | linker + TM + cytoplasmic tail |

`domain_architecture` and `domain_spans` expose the typed domain grammar
directly, for example:

- class I: `signal_peptide>g_alpha1>g_alpha2>c1_alpha3>tail_linker>transmembrane>cytoplasmic_tail`
- class II beta: `signal_peptide>g_beta1>c1_beta2>tail_linker>transmembrane>cytoplasmic_tail`

## How Parsing Works

The parser is **alignment-free** and **holistic**. It does not rely on one
absolute Cys position to define the mature start.

For each sequence it:

1. Enumerates all plausible Cys-Cys pairs in the Ig/C-like separation range.
2. Scores each pair as a candidate G-domain or C-like anchor using fold-topology
   evidence, especially the Trp41-like signal around `c1+14`.
3. Enumerates candidate SP boundaries and whole domain parses, including partial
   parses when only fragment evidence is available.
4. Chooses the best full parse using factored multiplicative scoring:
   three structural claims (SP grammar, domain architecture, completeness) each
   produce a [0,1] factor.  Contradictory evidence in any factor gates the score
   down multiplicatively, while missing evidence is a softer penalty.

The strongest evidence types are:

- SP cleavage grammar: hydrophobic h-region, short c-region, von Heijne `-3/-1`
  compatibility, exclusion of impossible `-3/-1` property pairs, and mild `+1`
  mature-sequence penalties.
- Domain-fold grammar: canonical G-domain versus C-like disulfide topology,
  including the IMGT-style `Cys11-Cys74` G-domain signature and the
  `Cys23/Trp41/Cys104` C-like grammar.
- Class-specific groove boundaries:
  - class I α1/α2 junction motifs
  - class I α2 -> α3 boundary motifs
  - class II α1 -> α2 and β1 -> β2 boundary motifs
- Soft priors on groove/support-domain lengths and TM support downstream.

The parser handles:

- full-length proteins with or without signal peptides
- SP-stripped deposits (`mature_start = 0`)
- common fragments:
  - class I exon 2 only -> `alpha1_only`
  - class I exon 3 only -> `alpha2_only`
  - class II exon 2-like fragments -> `fragment_fallback`
- low-evidence salvage:
  - class I from α3 C-like support only -> `inferred_from_alpha3`
  - class II beta from β1 groove pair only -> `beta1_only_fallback`
- true groove absence / insufficient structural evidence -> `missing_groove`

### Groove Status Values

These are the important parser-facing statuses:

| Status | Meaning |
|---|---|
| `ok` | Full decomposition from the main structural grammar |
| `alpha1_only` | Class I fragment consistent with α1 / exon 2 only |
| `alpha2_only` | Class I fragment consistent with α2 / exon 3 only |
| `fragment_fallback` | Short fragment retained as the observable groove half |
| `inferred_from_alpha3` | Class I salvage parse using a downstream α3 C-like anchor |
| `beta1_only_fallback` | Class II beta salvage parse using only the β1 groove pair |
| `missing_groove` | No recoverable groove architecture from the available evidence |
| `non_classical` | Non-classical class-I lineage flagged post-parse |
| `short` | Groove half too short to look functionally peptide-binding |

Pipeline-only statuses can still appear in CSV outputs:

| Status | Meaning |
|---|---|
| `not_applicable` | Row intentionally excluded from groove functionality, mainly B2M in build outputs |

## Literature Basis

The parser is built around conserved sequence grammar from the MHC literature:

- MHC domain organization is more conserved than short local motifs across
  vertebrates: [Primordial Linkage of β2-Microglobulin to the MHC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3805034/)
- IMGT domain numbering and the G-domain versus C-like disulfide grammar:
  [PMC3913909](https://pmc.ncbi.nlm.nih.gov/articles/PMC3913909/)
- Classical class-I domain layout and landmarks:
  [PMC2434379](https://pmc.ncbi.nlm.nih.gov/articles/PMC2434379/)
- Salmonid class-II alpha/beta cysteine topology and lineage-specific extra
  cysteines: [PMC2386828](https://pmc.ncbi.nlm.nih.gov/articles/PMC2386828/)
- Teleost class-II evolutionary divergence while retaining the same modular
  architecture: [PMC4219347](https://pmc.ncbi.nlm.nih.gov/articles/PMC4219347/)

Signal-peptide logic follows the standard SPase grammar:

- von Heijne `-3/-1` cleavage rule:
  [How signal sequences maintain cleavage specificity](https://pubmed.ncbi.nlm.nih.gov/6423828/)
- flanking `-2/+1` effects:
  [Flanking signal and mature peptide residues influence signal peptide cleavage](https://pmc.ncbi.nlm.nih.gov/articles/PMC2638155/)
- strong penalty for `+1 Pro`:
  [PubMed 1544500](https://pubmed.ncbi.nlm.nih.gov/1544500/)
- h-region / c-region structural context:
  [Structure of the human signal peptidase complex reveals the determinants for signal peptide cleavage](https://www.sciencedirect.com/science/article/pii/S1097276521006006)

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
- **[mhcgnomes](https://github.com/pirl-unc/mhcgnomes) >= 3.41.0** — allele parsing, species provenance, NHP taxonomy, and species-directed gene classification

mhcseqs emits full-binomial species aliases such as `HomoSapiens-A*02:01`
and `MusMusculus-K*b`. Its input registry accepts 476 current, historical, and
external-database prefix assignments—including all 137 designation tokens in
the current 125-organism IPD-MHC taxonomy register—and records an evidence URL
for every one. Colliding short codes require explicit species context. mhcseqs
does not invent abbreviated species prefixes, and mechanically generated
2+2/4+4/5+5 aliases are rejected unless that exact spelling has external
evidence. See the [prefix audit](data/species_prefix_audit_20260831.md).
The audit includes a versioned UniProtKB/UniSave provenance snapshot for all
239 records behind the 21 historical pairs that previously lacked raw inputs.

No alignment tools, BLAST, or structure databases are required.

## License

Apache 2.0
