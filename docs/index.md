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

Requires Python 3.11 or later. Runtime dependencies are
[mhcgnomes](https://github.com/pirl-unc/mhcgnomes) >= 3.92.0, `pandas >= 3, < 4`,
`datacache`, `biotite`, and `argcomplete`.

## Quick start

Install the prebuilt release catalog in seconds (or run `mhcseqs build` to
rebuild locally from IMGT/HLA and IPD-MHC FASTAs into `~/.cache/mhcseqs/`),
then look an allele up:

```bash
mhcseqs data install catalog
mhcseqs lookup "HLA-A*02:01"
```

The same thing from Python:

```python
import mhcseqs

mhcseqs.install_catalog()  # fast prebuilt release download; or mhcseqs.build()

r = mhcseqs.lookup("HLA-A*02:01")
print(r.allele, r.species_category, r.mhc_class, r.chain)
print("mature starts at", r.mature_start)
print("groove1", r.groove1)
print("architecture", r.domain_architecture)

# Fast cached O(1) index for batch lookups and pan-class IMGT features
catalog = mhcseqs.load_catalog()
feats = catalog.pan_class_features("HLA-DRB1*01:01")  # auto-pairs with DRA*01:01
```

```
HomoSapiens-A*02:01 human I alpha
mature starts at 24
groove1 GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEA
architecture signal_peptide>g_alpha1>g_alpha2>c1_alpha3>tail_linker>transmembrane>cytoplasmic_tail
```

`mhcseqs.install_catalog()` returns a `CatalogPaths` dataclass,
`mhcseqs.load_catalog()` loads a cached O(1) `AlleleCatalog`, and
`mhcseqs.load_sequences_dataframe()` loads the whole table at once.

## Where to go next

Roughly in reading order:

1. [Command line](cli.md): `build`, `lookup`, and the `data` subcommands for
   prebuilt catalog installs, cached downloads, and versioned datasets.
2. [Python API](python-api.md): `lookup` and `AlleleRecord`, `AlleleCatalog`
   batch/heterodimer/pan-class helpers, DataFrame loaders, allele and species
   utilities, and the module reference.
3. [Structural decomposition](parsing.md): the domain grammar, how the
   alignment-free parser decides a decomposition, the groove status codes, and
   the literature the grammar comes from.
4. [Pan-class IMGT alignment](pan-class-alignment.md): the unified 202-slot
   MHC-I / MHC-II IMGT frame, 46-slot NetMHCpan/NetMHCIIpan pocket union,
   insertion tracking, and joint predictor design.
5. [Datasets and columns](data.md): what the build contains, what each column
   means, and the independently versioned `mhc-proteins` record layer.

## Project links

- Repository: <https://github.com/pirl-unc/mhcseqs>
- PyPI: <https://pypi.org/project/mhcseqs/>
- README: <https://github.com/pirl-unc/mhcseqs/blob/main/README.md>
