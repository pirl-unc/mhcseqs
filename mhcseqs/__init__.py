"""MHC sequence curation and binding groove extraction.

Quick start::

    import mhcseqs

    # Build the database (downloads to ~/.cache/mhcseqs/, only needed once)
    paths = mhcseqs.build()

    # Look up any allele → AlleleRecord
    r = mhcseqs.lookup("HLA-A*02:01")
    r.groove1            # α1 domain
    r.groove2            # α2 domain
    r.domains            # typed domain grammar spans
    r.domain_architecture
    r.mature_sequence    # signal peptide removed
    r.sequence           # full protein (with signal peptide)

    # Look up with mutations (IEDB-style)
    m = mhcseqs.lookup("HLA-A*02:01", mutations=["K66A"])

    # Load all sequences as a DataFrame
    df = mhcseqs.load_sequences_dataframe()
"""

from dataclasses import dataclass
from typing import Sequence

from .alleles import (
    infer_gene,
    infer_mhc_class,
    infer_species,
    normalize_allele_name,
    normalize_mhc_class,
    parse_allele_name,
)
from .datafiles import (
    FULL_SEQS_CSV,
    MERGE_REPORT,
    RAW_CSV,
    VALIDATION_REPORT,
    build_query_candidates,
    fasta_dir,
    find_csv,
    match_allele_row,
    resolve_data_dir,
)
from .domain_parsing import (
    NON_GROOVE_GENES,
    AlleleRecord,
    RawAllele,
    SequenceFeatures,
    StructuralDomain,
    analyze_sequence,
    apply_mutations,
    decompose_class_i,
    decompose_class_ii_alpha,
    decompose_class_ii_beta,
    decompose_domains,
    extract_groove,
    find_cys_pairs,
    infer_structural_domains,
    parse_class_i,
    parse_class_ii_alpha,
    parse_class_ii_beta,
)
from .download import SOURCES, download_all
from .imgt import (
    CONSERVED_CYS_POSITIONS,
    GALPHA2_GAP_POSITIONS,
    GALPHA2_POSITIONS,
    HELIX_INSERTIONS,
    STRUCTURAL_ELEMENTS,
    imgt_to_mature,
    imgt_to_mature_class_i,
    mature_to_imgt,
    mature_to_imgt_class_i,
    structural_element,
)
from .mhc_protein_dataset import (
    ProteinDatasetError,
    ProteinDatasetPaths,
    ProteinSourceBundlePaths,
    available_mhc_protein_dataset_versions,
    default_mhc_protein_dataset_version,
    install_mhc_protein_dataset,
    install_mhc_protein_source_bundle,
    iter_mhc_protein_records,
    load_mhc_protein_dataframe,
    load_mhc_protein_records,
    mhc_protein_dataset_paths,
    validate_mhc_protein_dataset,
)
from .pipeline import (
    FULL_FIELDS,
    RAW_FIELDS,
    build_full_seqs,
    build_raw_index,
)
from .species import (
    LATIN_NAMES,
    MHC_SPECIES_CATEGORIES,
    full_species_name_alias,
    get_canonical_prefix,
    get_latin_name,
    normalize_mhc_species,
    normalize_species,
)
from .version import __version__


@dataclass(frozen=True)
class BuildPaths:
    """Paths to CSV files produced by :func:`build`.

    Attributes
    ----------
    raw : str
        Every protein entry from all sources (``mhc-seqs-raw.csv``).
    sequences : str
        One representative per two-field allele group, with full sequence,
        groove decomposition, and metadata (``mhc-full-seqs.csv``).
    merge_report : str
        Human-readable report of the merge/dedup step
        (``mhc-merge-report.txt``).
    validation_report : str
        Human-readable validation sanity checks
        (``mhc-validation-report.txt``).
    """

    raw: str
    sequences: str
    merge_report: str
    validation_report: str


def default_data_dir() -> str:
    """Return the default directory for downloads and built CSVs.

    Uses ``$MHCSEQS_DATA`` if set to a non-empty value, otherwise
    ``~/.cache/mhcseqs`` (resolved by :func:`mhcseqs.datafiles.resolve_data_dir`).
    The directory is created on first use.
    """
    d = resolve_data_dir()
    d.mkdir(parents=True, exist_ok=True)
    return str(d)


def build(
    output_dir: str | None = None,
    data_dir: str | None = None,
    force_download: bool = False,
) -> BuildPaths:
    """Run the full build pipeline: download, parse, extract grooves.

    Parameters
    ----------
    output_dir : str, optional
        Directory for the built CSVs and reports. Defaults to
        :func:`default_data_dir` (``$MHCSEQS_DATA`` or the resolved cache dir).
    data_dir : str, optional
        Directory for the downloaded IMGT/HLA and IPD-MHC FASTA files.
        Defaults to the ``fasta/`` subdirectory of :func:`default_data_dir`,
        matching ``mhcseqs build --output-dir``.
    force_download : bool, optional
        Re-fetch the source FASTA files even when a cached copy exists.

    Returns
    -------
    BuildPaths
        Paths to the four generated files.

    Raises
    ------
    OSError
        If a source download fails (``urllib.error.URLError`` is a subclass)
        or an output file cannot be written.
    """
    from pathlib import Path

    from .validate import format_validation_report, validate_build

    dd_root = Path(default_data_dir())
    out = dd_root if output_dir is None else Path(output_dir)
    dd = fasta_dir(dd_root) if data_dir is None else Path(data_dir)

    paths = download_all(dd, force=force_download)
    fasta_inputs = [
        (paths["imgt_hla"], SOURCES["imgt_hla"]["label"]),
        (paths["ipd_mhc"], SOURCES["ipd_mhc"]["label"]),
    ]

    raw_csv = out / RAW_CSV
    build_raw_index(fasta_inputs, raw_csv)

    full_csv = out / FULL_SEQS_CSV
    report_path = out / MERGE_REPORT
    build_full_seqs(raw_csv, full_csv, report_path=report_path)

    warnings, stats = validate_build(raw_csv, full_csv)
    validation_path = out / VALIDATION_REPORT
    with open(validation_path, "w", encoding="utf-8") as f:
        f.write(format_validation_report(warnings, stats) + "\n")

    return BuildPaths(
        raw=str(raw_csv),
        sequences=str(full_csv),
        merge_report=str(report_path),
        validation_report=str(validation_path),
    )


def load_raw_dict(path: str | None = None) -> list[dict]:
    """Load the raw sequence CSV as a list of dicts.

    Contains every individual protein entry from all sources,
    including duplicates and fragments.

    Parameters
    ----------
    path : str, optional
        CSV to read. By default ``mhc-seqs-raw.csv`` is located via
        :func:`find_csv` (current directory, then resolved data directory,
        then package parent directory).

    Returns
    -------
    list of dict
        One ``{column: value}`` dict per row; all values are strings.

    Raises
    ------
    FileNotFoundError
        If *path* is not given and no built CSV is found.
    """
    import csv as _csv

    p = path or find_csv(RAW_CSV)
    with open(p, "r", encoding="utf-8") as f:
        return list(_csv.DictReader(f))


def load_raw_dataframe(path: str | None = None):
    """Load the raw sequence CSV as a pandas DataFrame.

    Parameters
    ----------
    path : str, optional
        CSV to read. Located as in :func:`load_raw_dict` when omitted.

    Returns
    -------
    pandas.DataFrame
        One row per raw protein entry.

    Raises
    ------
    ImportError
        If pandas is not installed (``pip install "mhcseqs[pandas]"``).
    FileNotFoundError
        If *path* is not given and no built CSV is found.
    """
    import pandas as pd

    p = path or find_csv(RAW_CSV)
    return pd.read_csv(p)


def load_sequences_dict(path: str | None = None) -> list[dict]:
    """Load the full-sequence representatives CSV as a list of dicts.

    One row per two-field allele with the best available full-length
    protein sequence, groove decomposition, species metadata, and
    signal peptide annotation.

    Parameters
    ----------
    path : str, optional
        CSV to read. By default ``mhc-full-seqs.csv`` is located via
        :func:`find_csv` (current directory, then resolved data directory,
        then package parent directory).

    Returns
    -------
    list of dict
        One ``{column: value}`` dict per row; all values are strings.

    Raises
    ------
    FileNotFoundError
        If *path* is not given and no built CSV is found.
    """
    import csv as _csv

    p = path or find_csv(FULL_SEQS_CSV)
    with open(p, "r", encoding="utf-8") as f:
        return list(_csv.DictReader(f))


def load_sequences_dataframe(path: str | None = None):
    """Load the full-sequence representatives CSV as a pandas DataFrame.

    One row per two-field allele with full protein, mature sequence,
    groove decomposition, species metadata, and allele status flags.

    Parameters
    ----------
    path : str, optional
        CSV to read. Located as in :func:`load_sequences_dict` when omitted.

    Returns
    -------
    pandas.DataFrame
        One row per two-field allele.

    Raises
    ------
    ImportError
        If pandas is not installed (``pip install "mhcseqs[pandas]"``).
    FileNotFoundError
        If *path* is not given and no built CSV is found.
    """
    import pandas as pd

    p = path or find_csv(FULL_SEQS_CSV)
    return pd.read_csv(p)


def lookup(
    allele: str,
    *,
    mutations: Sequence = (),
    search_dir: str | None = None,
) -> AlleleRecord:
    """Look up a single allele and return a fully parsed AlleleRecord.

    Finds the allele in the built full-sequences CSV, then runs
    ``decompose_domains()`` on it to produce a complete record with all
    fields populated (including anchor Cys positions).

    Optionally applies mutations (same format as ``decompose_domains()``).

    Parameters
    ----------
    allele : str
        Allele name. A row matches when its ``two_field_allele`` or
        ``representative_allele`` equals the query exactly, ignoring case, or
        after :func:`normalize_allele_name` (so ``"HLA-A*02:01"`` finds
        ``"HomoSapiens-A*02:01"``). The first matching row is used.
    mutations : sequence of str, optional
        Point mutations such as ``"K66A"`` to apply to the looked-up sequence.
    search_dir : str, optional
        Directory searched first for ``mhc-full-seqs.csv``, before the current
        directory, :func:`default_data_dir`, and the package's parent directory.

    Returns
    -------
    AlleleRecord
        Domain decomposition of the representative sequence, with the CSV's
        species, source, and full sequence fields filled in.

    Raises
    ------
    FileNotFoundError
        If the full-seqs CSV hasn't been built yet.
    KeyError
        If no matching allele is found.

    Examples
    --------
    >>> r = mhcseqs.lookup("HLA-A*02:01")           # doctest: +SKIP
    >>> r.groove1[:10]                               # doctest: +SKIP
    'GSHSMRYFFT'
    >>> m = mhcseqs.lookup("HLA-A*02:01",            # doctest: +SKIP
    ...                    mutations=["K66A"])
    >>> m.groove1[65]                                # doctest: +SKIP
    'A'
    """
    import csv as _csv
    from dataclasses import replace

    query = allele.strip()
    query_candidates = build_query_candidates(query)

    p = find_csv(FULL_SEQS_CSV, search_dir=search_dir)
    hit: dict | None = None
    with open(p, "r", encoding="utf-8") as f:
        for row in _csv.DictReader(f):
            if match_allele_row(row, query, query_candidates=query_candidates):
                hit = row
                break

    if hit is None:
        raise KeyError(f"Allele {query!r} not found in {p}")

    seq = hit.get("sequence", "")
    mhc_class = hit.get("mhc_class", "")
    chain = hit.get("chain", "")
    gene = hit.get("gene", "")
    allele_name = hit.get("two_field_allele", "")
    full_allele_name = hit.get("representative_allele", "")

    result = decompose_domains(
        seq,
        mhc_class=mhc_class,
        chain=chain,
        allele=allele_name,
        gene=gene,
        species=hit.get("species", ""),
        mutations=mutations,
    )

    return replace(
        result,
        full_allele=full_allele_name,
        species=hit.get("species", ""),
        species_category=hit.get("species_category", ""),
        species_prefix=hit.get("species_prefix", ""),
        source=hit.get("source", ""),
        source_id=hit.get("source_id", ""),
        sequence=seq,
        seq_len=len(seq),
    )


__all__ = [
    "__version__",
    "BuildPaths",
    "default_data_dir",
    "find_csv",
    "build",
    "load_raw_dict",
    "load_raw_dataframe",
    "load_sequences_dict",
    "load_sequences_dataframe",
    "available_mhc_protein_dataset_versions",
    "default_mhc_protein_dataset_version",
    "install_mhc_protein_dataset",
    "install_mhc_protein_source_bundle",
    "iter_mhc_protein_records",
    "load_mhc_protein_records",
    "load_mhc_protein_dataframe",
    "mhc_protein_dataset_paths",
    "validate_mhc_protein_dataset",
    "ProteinDatasetError",
    "ProteinDatasetPaths",
    "ProteinSourceBundlePaths",
    "lookup",
    "build_raw_index",
    "build_full_seqs",
    "apply_mutations",
    "decompose_domains",
    "decompose_class_i",
    "decompose_class_ii_alpha",
    "decompose_class_ii_beta",
    "extract_groove",
    "find_cys_pairs",
    "parse_class_i",
    "parse_class_ii_alpha",
    "parse_class_ii_beta",
    "analyze_sequence",
    "infer_structural_domains",
    "RawAllele",
    "AlleleRecord",
    "SequenceFeatures",
    "StructuralDomain",
    "NON_GROOVE_GENES",
    "parse_allele_name",
    "normalize_allele_name",
    "infer_gene",
    "infer_mhc_class",
    "infer_species",
    "normalize_mhc_class",
    "normalize_species",
    "normalize_mhc_species",
    "get_latin_name",
    "full_species_name_alias",
    "get_canonical_prefix",
    "MHC_SPECIES_CATEGORIES",
    "LATIN_NAMES",
    "download_all",
    "SOURCES",
    "RAW_FIELDS",
    "FULL_FIELDS",
    "mature_to_imgt",
    "mature_to_imgt_class_i",
    "imgt_to_mature",
    "imgt_to_mature_class_i",
    "structural_element",
    "STRUCTURAL_ELEMENTS",
    "GALPHA2_GAP_POSITIONS",
    "GALPHA2_POSITIONS",
    "HELIX_INSERTIONS",
    "CONSERVED_CYS_POSITIONS",
]
