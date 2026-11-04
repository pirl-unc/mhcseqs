# mhcseqs

mhcseqs is a self-contained pipeline for downloading, curating, and extracting
binding grooves from MHC (Major Histocompatibility Complex) protein sequences.
Given an allele name it returns the full protein sequence together with an
explicit structural decomposition: signal peptide, groove domains, Ig-fold
support domain, and tail.

In the openvax stack it is the sequence layer beside
[mhcgnomes](https://github.com/pirl-unc/mhcgnomes), which it uses for allele
parsing, species provenance and gene classification.

No alignment tools, BLAST, or structure databases are required.

## Install

```bash
pip install mhcseqs
```

Requires Python 3.11 or later. The only required dependency is
`mhcgnomes >= 3.41.0`. The DataFrame loaders also need pandas 3, which
`pip install "mhcseqs[pandas]"` installs.

## Quick start

Build the database once (downloads IMGT/HLA and IPD-MHC sources into
`~/.cache/mhcseqs/`), then look an allele up:

```bash
mhcseqs build
mhcseqs lookup "HLA-A*02:01"
```

The same thing from Python:

```python
import mhcseqs

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

`mhcseqs.build()` returns a `BuildPaths` dataclass, and
`mhcseqs.load_sequences_dataframe()` loads the whole table at once.

## Where to go next

Roughly in reading order:

1. [Command line](cli.md): `build`, `lookup`, and the `data` subcommands for
   cached downloads and versioned datasets.
2. [Python API](python-api.md): `lookup` and the `AlleleRecord` fields, the
   DataFrame loaders, the allele and species utilities, and the module
   reference.
3. [Structural decomposition](parsing.md): the domain grammar, how the
   alignment-free parser decides a decomposition, the groove status codes, and
   the literature the grammar comes from.
4. [Datasets and columns](data.md): what the build contains, what each column
   means, and the independently versioned `mhc-proteins` record layer.

## Project links

- Repository: <https://github.com/pirl-unc/mhcseqs>
- PyPI: <https://pypi.org/project/mhcseqs/>
- README: <https://github.com/pirl-unc/mhcseqs/blob/main/README.md>
