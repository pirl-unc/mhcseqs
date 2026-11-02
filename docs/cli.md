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
# Build all output CSVs (writes to ~/.cache/mhcseqs/)
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
```

```
--- mhc-full-seqs.csv ---
  two_field_allele: HomoSapiens-A*02:01
  representative_allele: HomoSapiens-A*02:01:01:100
  protein_seq_selection: nested_longest
  gene: A
  mhc_class: I
  chain: alpha
  species: Homo sapiens
  species_category: human
  species_prefix: HomoSapiens
  source: imgt
  source_id: HLA24674
  seq_len: 365
  sequence: MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVG... (365 aa)
  mature_start: 24
  mature_sequence: GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAP... (341 aa)
  groove1: GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEA
  groove2: GSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTD
  groove1_len: 90
  groove2_len: 93
  ig_domain: APKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSS
  ig_domain_len: 95
  tail: QPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSLTACKV
  tail_len: 63
  domain_architecture: signal_peptide>g_alpha1>g_alpha2>c1_alpha3>tail_linker>transmembrane>cytoplasmic_tail
  domain_spans: signal_peptide:1-24;g_alpha1:25-114;g_alpha2:115-207;c1_alpha3:208-302;tail_linker:303-307;transmembrane:308-331;cytoplasmic_tail:332-365
  groove_status: ok
  groove_flags:
```

The `groove_seq` line is omitted above for width; the command prints it
between `groove2` and `groove1_len`. Non-human alleles work the same way:
`mhcseqs lookup "Mamu-A1*001:01"`. Use `--output-dir` when the CSVs are not in
`~/.cache/mhcseqs`.

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

## Version

```bash
mhcseqs --version
```
