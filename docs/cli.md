# Command line

`mhcseqs` is one executable with three subcommands: `build` downloads the
sources and writes the CSVs, `lookup` reads one allele back out of them, and
`data` manages the cache and the versioned datasets. Run
`mhcseqs <command> --help` for every option.

| Command | Does |
|---|---|
| `mhcseqs build` | Download sources and build CSVs |
| `mhcseqs lookup` | Look up an allele from built CSVs |
| `mhcseqs data list` | List cached source and built data files |
| `mhcseqs data install` | Install an immutable versioned dataset |
| `mhcseqs data path` | Print the path to an installed versioned dataset |
| `mhcseqs data refresh` | Re-download source FASTA(s) |
| `mhcseqs data clear` | Delete cached source FASTA(s) |

## Build the database

```bash
# Build all output CSVs (writes to ~/.cache/mhcseqs/; set MHCSEQS_DATA to use another
# data directory for downloads, built CSVs and versioned datasets)
mhcseqs build

# Build to a specific directory instead
mhcseqs build --output-dir output/

# Re-download source FASTAs (IMGT/HLA and IPD-MHC publish from a rolling "Latest")
mhcseqs build --force-download
```

`--data-dir` sets where the downloaded FASTA files live, separately from
`--output-dir` for the built CSVs.

## Look up an allele

```bash
mhcseqs lookup "HLA-A*02:01"
mhcseqs lookup "HLA-A*02:01" --format fasta
mhcseqs lookup "HLA-A*02:01" --format json --all --output a0201.json
```

Supported `--format` values are `text` (default), `csv`, `tsv`, `json`, and
`fasta`. Every format emits full protein and mature sequences (no truncation).
Pass `--all` to emit every matching row instead of only the first match per
CSV, and `--output PATH` (`-o`) to write output to a file. Missing matches
report an error on `stderr` and exit with status 1.

`mhcseqs lookup` uses the same CSV search order as `mhcseqs.lookup()`:
`--output-dir` (when given), then the current working directory, then the
resolved `mhcseqs` data directory. Non-human alleles work the same way:
`mhcseqs lookup "Mamu-A1*001:01"`.

## Inspect and manage cached data

```bash
mhcseqs data list
mhcseqs data list --details  # full paths and source URLs
```

```
Data directory: /Users/you/.cache/mhcseqs

Source FASTA files
FILE                SIZE     UPDATED  STATUS
hla_prot.fasta      14.0 MB  61d ago  present
ipd_mhc_prot.fasta  3.6 MB   61d ago  present

Built CSVs and reports
FILE                       SIZE     UPDATED  STATUS
mhc-seqs-raw.csv           30.6 MB  61d ago  present
mhc-full-seqs.csv          58.6 MB  61d ago  present
mhc-merge-report.txt       16.4 KB  61d ago  present
mhc-validation-report.txt  72.8 KB  61d ago  present

Versioned datasets
DATASET       VERSION             STATUS
mhc-proteins  uniprot-2026_03-r1  not installed
mhc-proteins  uniprot-2026_03-r2  not installed

Use 'mhcseqs data list --details' for full paths and source URLs.
```

```bash
# Re-download source FASTAs without rebuilding
mhcseqs data refresh

# Delete cached files
mhcseqs data clear               # source FASTAs
mhcseqs data clear --built       # also remove built CSVs/reports
mhcseqs data clear --built-only  # only built CSVs/reports
mhcseqs data clear --source imgt_hla   # just one source
```

## Versioned datasets

```bash
mhcseqs data install mhc-proteins
mhcseqs data path mhc-proteins
mhcseqs data install mhc-proteins --version uniprot-2026_03-r2
```

Code releases and data releases are deliberately independent, so pin a data
version in automated work. See [datasets and columns](data.md) for what
`mhc-proteins` contains.

## Tab completion

`mhcseqs` supports shell tab completion via `argcomplete`. Enable it once for
every OpenVax command-line tool:

```bash
activate-global-python-argcomplete
```

Or enable it for `mhcseqs` alone:

```bash
eval "$(register-python-argcomplete mhcseqs)"
```

## Version

```bash
mhcseqs --version
```
