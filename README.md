[![CI](https://github.com/pirl-unc/mhcseqs/actions/workflows/ci.yml/badge.svg)](https://github.com/pirl-unc/mhcseqs/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/mhcseqs.svg?maxAge=1000)](https://pypi.org/project/mhcseqs/)
[![Python versions](https://img.shields.io/pypi/pyversions/mhcseqs.svg)](https://pypi.org/project/mhcseqs/)
[![Docs](https://img.shields.io/badge/docs-pirl--unc.github.io%2Fmhcseqs-blue)](https://pirl-unc.github.io/mhcseqs/)

# mhcseqs

Self-contained pipeline for downloading, curating, and extracting binding
grooves from MHC (Major Histocompatibility Complex) protein sequences. Given
an allele name it returns the full protein sequence together with an explicit
structural decomposition: signal peptide, groove domains, Ig-fold support
domain, and tail.

In the openvax stack it is the sequence layer beside
[mhcgnomes](https://github.com/pirl-unc/mhcgnomes), which it uses for allele
parsing, species provenance and gene classification. No alignment tools,
BLAST, or structure databases are required.

## Install

```bash
pip install mhcseqs
```

Requires Python 3.11 or later. The only required dependency is
[mhcgnomes](https://github.com/pirl-unc/mhcgnomes) >= 3.41.0. The DataFrame
loaders also need pandas 3, which `pip install "mhcseqs[pandas]"` installs.

## Quick start

Build the database once (downloads IMGT/HLA and IPD-MHC sources into
`~/.cache/mhcseqs/`), then look an allele up:

```python
import mhcseqs

mhcseqs.build()  # only needed once; returns a BuildPaths dataclass

r = mhcseqs.lookup("HLA-A*02:01")
print(r.allele, r.species_category, r.mhc_class, r.chain)
print("mature starts at", r.mature_start)
print("groove1", r.groove1)
print("architecture", r.domain_architecture)
```

```
HomoSapiens-A*02:01 human I alpha
mature starts at 24
groove1 GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEA
architecture signal_peptide>g_alpha1>g_alpha2>c1_alpha3>tail_linker>transmembrane>cytoplasmic_tail
```

`lookup` also takes IEDB-style mutations
(`mhcseqs.lookup("HLA-A*02:01", mutations=["K66A", "D77S"])`), and
`load_sequences_dataframe()` / `load_sequences_dict()` load the whole build at
once. See the [Python API guide](https://pirl-unc.github.io/mhcseqs/python-api/).

## Command line

```bash
mhcseqs build                    # download sources and build CSVs
mhcseqs lookup "HLA-A*02:01"     # print one allele's record
mhcseqs data list                # cached downloads and built CSVs
mhcseqs data install mhc-proteins  # versioned UniProt record dataset
```

`mhcseqs <command> --help` lists every option; the
[command-line guide](https://pirl-unc.github.io/mhcseqs/cli/) shows the output
of each.

## Documentation

Full documentation: [pirl-unc.github.io/mhcseqs](https://pirl-unc.github.io/mhcseqs/).

- [Command line](https://pirl-unc.github.io/mhcseqs/cli/): `build`, `lookup`, and the `data` subcommands
- [Python API](https://pirl-unc.github.io/mhcseqs/python-api/): `AlleleRecord` fields, DataFrame loaders, allele and species utilities, module reference
- [Structural decomposition](https://pirl-unc.github.io/mhcseqs/parsing/): the domain grammar, how the alignment-free parser scores a decomposition, groove status codes, literature basis
- [Datasets and columns](https://pirl-unc.github.io/mhcseqs/data/): what a build contains, column meanings, and the versioned `mhc-proteins` record layer

## Development

```bash
git clone https://github.com/pirl-unc/mhcseqs.git
cd mhcseqs
./develop.sh          # uv pip install -e ".[dev]"
./test.sh             # pytest
./lint.sh             # ruff check and ruff format --check
```

To work on the documentation site:

```bash
pip install -e ".[docs]"
mkdocs serve
mkdocs build --strict
```

## License

The mhcseqs code is licensed under
[Apache 2.0](https://github.com/pirl-unc/mhcseqs/blob/main/LICENSE). The independently
distributed UniProt-derived `mhc-proteins` data is licensed under
[CC BY 4.0](https://www.uniprot.org/help/license/) with attribution to the
UniProt Consortium.
