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

from __future__ import annotations

import importlib
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Sequence

from .version import __version__

if TYPE_CHECKING:
    from .alleles import (
        infer_gene,
        infer_mhc_class,
        infer_species,
        normalize_allele_name,
        normalize_mhc_class,
        parse_allele_name,
    )
    from .catalog import (
        AlleleCatalog,
        CatalogEntry,
        CatalogError,
        CatalogNotInstalledError,
        CatalogPaths,
        clear_catalog_cache,
        install_catalog,
        load_catalog,
        validate_catalog,
    )
    from .datafiles import find_csv
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
        FIXED_FRAME_ALPHABET,
        GALPHA1_POSITIONS,
        GALPHA2_GAP_POSITIONS,
        GALPHA2_POSITIONS,
        GALPHA_GAP_POSITIONS,
        GALPHA_POSITIONS,
        GBETA_GAP_POSITIONS,
        GBETA_POSITIONS,
        HELIX_INSERTIONS,
        HLA_A0201_GALPHA1,
        HLA_A0201_GALPHA2,
        HLA_DRA0101_GALPHA,
        HLA_DRB1_0101_GBETA,
        IMGT_CLASS_I_GROOVE_SLOTS,
        IMGT_CLASS_I_PROFILE_PATH,
        IMGT_CLASS_II_GROOVE_SLOTS,
        IMGT_CLASS_II_PROFILE_PATH,
        IMGT_GALPHA1_SLOTS,
        IMGT_GALPHA2_SLOTS,
        IMGT_GALPHA_SLOTS,
        IMGT_GBETA_SLOTS,
        IMGT_PAN_CLASS_GROOVE_SLOTS,
        IMGT_PAN_HALF1_SLOTS,
        IMGT_PAN_HALF2_SLOTS,
        STRUCTURAL_ELEMENTS,
        AlignedResidue,
        DomainAlignment,
        GrooveAlignment,
        GrooveDeletion,
        GrooveInsertion,
        align_domain_to_imgt,
        align_groove_to_imgt,
        imgt_to_mature,
        imgt_to_mature_class_i,
        map_class_ii_groove_positions,
        map_groove_positions,
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
    from .pipeline import FULL_FIELDS, RAW_FIELDS, build_full_seqs, build_raw_index
    from .species import (
        LATIN_NAMES,
        MHC_SPECIES_CATEGORIES,
        full_species_name_alias,
        get_canonical_prefix,
        get_latin_name,
        get_ncbi_taxid,
        get_species_clade,
        normalize_mhc_species,
        normalize_species,
    )

_LAZY_EXPORTS: dict[str, tuple[str, str]] = {
    # alleles
    "infer_gene": (".alleles", "infer_gene"),
    "infer_mhc_class": (".alleles", "infer_mhc_class"),
    "infer_species": (".alleles", "infer_species"),
    "normalize_allele_name": (".alleles", "normalize_allele_name"),
    "normalize_mhc_class": (".alleles", "normalize_mhc_class"),
    "parse_allele_name": (".alleles", "parse_allele_name"),
    # catalog
    "AlleleCatalog": (".catalog", "AlleleCatalog"),
    "CatalogEntry": (".catalog", "CatalogEntry"),
    "CatalogError": (".catalog", "CatalogError"),
    "CatalogNotInstalledError": (".catalog", "CatalogNotInstalledError"),
    "CatalogPaths": (".catalog", "CatalogPaths"),
    "clear_catalog_cache": (".catalog", "clear_catalog_cache"),
    "install_catalog": (".catalog", "install_catalog"),
    "load_catalog": (".catalog", "load_catalog"),
    "validate_catalog": (".catalog", "validate_catalog"),
    # datafiles
    "FULL_SEQS_CSV": (".datafiles", "FULL_SEQS_CSV"),
    "MERGE_REPORT": (".datafiles", "MERGE_REPORT"),
    "RAW_CSV": (".datafiles", "RAW_CSV"),
    "VALIDATION_REPORT": (".datafiles", "VALIDATION_REPORT"),
    "build_query_candidates": (".datafiles", "build_query_candidates"),
    "fasta_dir": (".datafiles", "fasta_dir"),
    "find_csv": (".datafiles", "find_csv"),
    "match_allele_row": (".datafiles", "match_allele_row"),
    "resolve_data_dir": (".datafiles", "resolve_data_dir"),
    # domain_parsing
    "NON_GROOVE_GENES": (".domain_parsing", "NON_GROOVE_GENES"),
    "AlleleRecord": (".domain_parsing", "AlleleRecord"),
    "RawAllele": (".domain_parsing", "RawAllele"),
    "SequenceFeatures": (".domain_parsing", "SequenceFeatures"),
    "StructuralDomain": (".domain_parsing", "StructuralDomain"),
    "analyze_sequence": (".domain_parsing", "analyze_sequence"),
    "apply_mutations": (".domain_parsing", "apply_mutations"),
    "decompose_class_i": (".domain_parsing", "decompose_class_i"),
    "decompose_class_ii_alpha": (".domain_parsing", "decompose_class_ii_alpha"),
    "decompose_class_ii_beta": (".domain_parsing", "decompose_class_ii_beta"),
    "decompose_domains": (".domain_parsing", "decompose_domains"),
    "extract_groove": (".domain_parsing", "extract_groove"),
    "find_cys_pairs": (".domain_parsing", "find_cys_pairs"),
    "infer_structural_domains": (".domain_parsing", "infer_structural_domains"),
    "parse_class_i": (".domain_parsing", "parse_class_i"),
    "parse_class_ii_alpha": (".domain_parsing", "parse_class_ii_alpha"),
    "parse_class_ii_beta": (".domain_parsing", "parse_class_ii_beta"),
    # download
    "SOURCES": (".download", "SOURCES"),
    "download_all": (".download", "download_all"),
    # imgt
    "CONSERVED_CYS_POSITIONS": (".imgt", "CONSERVED_CYS_POSITIONS"),
    "FIXED_FRAME_ALPHABET": (".imgt", "FIXED_FRAME_ALPHABET"),
    "GALPHA1_POSITIONS": (".imgt", "GALPHA1_POSITIONS"),
    "GALPHA2_GAP_POSITIONS": (".imgt", "GALPHA2_GAP_POSITIONS"),
    "GALPHA2_POSITIONS": (".imgt", "GALPHA2_POSITIONS"),
    "GALPHA_GAP_POSITIONS": (".imgt", "GALPHA_GAP_POSITIONS"),
    "GALPHA_POSITIONS": (".imgt", "GALPHA_POSITIONS"),
    "GBETA_GAP_POSITIONS": (".imgt", "GBETA_GAP_POSITIONS"),
    "GBETA_POSITIONS": (".imgt", "GBETA_POSITIONS"),
    "HELIX_INSERTIONS": (".imgt", "HELIX_INSERTIONS"),
    "HLA_A0201_GALPHA1": (".imgt", "HLA_A0201_GALPHA1"),
    "HLA_A0201_GALPHA2": (".imgt", "HLA_A0201_GALPHA2"),
    "HLA_DRA0101_GALPHA": (".imgt", "HLA_DRA0101_GALPHA"),
    "HLA_DRB1_0101_GBETA": (".imgt", "HLA_DRB1_0101_GBETA"),
    "IMGT_CLASS_I_GROOVE_SLOTS": (".imgt", "IMGT_CLASS_I_GROOVE_SLOTS"),
    "IMGT_CLASS_I_PROFILE_PATH": (".imgt", "IMGT_CLASS_I_PROFILE_PATH"),
    "IMGT_CLASS_II_GROOVE_SLOTS": (".imgt", "IMGT_CLASS_II_GROOVE_SLOTS"),
    "IMGT_CLASS_II_PROFILE_PATH": (".imgt", "IMGT_CLASS_II_PROFILE_PATH"),
    "IMGT_GALPHA1_SLOTS": (".imgt", "IMGT_GALPHA1_SLOTS"),
    "IMGT_GALPHA2_SLOTS": (".imgt", "IMGT_GALPHA2_SLOTS"),
    "IMGT_GALPHA_SLOTS": (".imgt", "IMGT_GALPHA_SLOTS"),
    "IMGT_GBETA_SLOTS": (".imgt", "IMGT_GBETA_SLOTS"),
    "IMGT_PAN_CLASS_GROOVE_SLOTS": (".imgt", "IMGT_PAN_CLASS_GROOVE_SLOTS"),
    "IMGT_PAN_HALF1_SLOTS": (".imgt", "IMGT_PAN_HALF1_SLOTS"),
    "IMGT_PAN_HALF2_SLOTS": (".imgt", "IMGT_PAN_HALF2_SLOTS"),
    "STRUCTURAL_ELEMENTS": (".imgt", "STRUCTURAL_ELEMENTS"),
    "AlignedResidue": (".imgt", "AlignedResidue"),
    "DomainAlignment": (".imgt", "DomainAlignment"),
    "GrooveAlignment": (".imgt", "GrooveAlignment"),
    "GrooveDeletion": (".imgt", "GrooveDeletion"),
    "GrooveInsertion": (".imgt", "GrooveInsertion"),
    "align_domain_to_imgt": (".imgt", "align_domain_to_imgt"),
    "align_groove_to_imgt": (".imgt", "align_groove_to_imgt"),
    "imgt_to_mature": (".imgt", "imgt_to_mature"),
    "imgt_to_mature_class_i": (".imgt", "imgt_to_mature_class_i"),
    "map_class_ii_groove_positions": (".imgt", "map_class_ii_groove_positions"),
    "map_groove_positions": (".imgt", "map_groove_positions"),
    "mature_to_imgt": (".imgt", "mature_to_imgt"),
    "mature_to_imgt_class_i": (".imgt", "mature_to_imgt_class_i"),
    "structural_element": (".imgt", "structural_element"),
    # mhc_protein_dataset
    "ProteinDatasetError": (".mhc_protein_dataset", "ProteinDatasetError"),
    "ProteinDatasetPaths": (".mhc_protein_dataset", "ProteinDatasetPaths"),
    "ProteinSourceBundlePaths": (".mhc_protein_dataset", "ProteinSourceBundlePaths"),
    "available_mhc_protein_dataset_versions": (".mhc_protein_dataset", "available_mhc_protein_dataset_versions"),
    "default_mhc_protein_dataset_version": (".mhc_protein_dataset", "default_mhc_protein_dataset_version"),
    "install_mhc_protein_dataset": (".mhc_protein_dataset", "install_mhc_protein_dataset"),
    "install_mhc_protein_source_bundle": (".mhc_protein_dataset", "install_mhc_protein_source_bundle"),
    "iter_mhc_protein_records": (".mhc_protein_dataset", "iter_mhc_protein_records"),
    "load_mhc_protein_dataframe": (".mhc_protein_dataset", "load_mhc_protein_dataframe"),
    "load_mhc_protein_records": (".mhc_protein_dataset", "load_mhc_protein_records"),
    "mhc_protein_dataset_paths": (".mhc_protein_dataset", "mhc_protein_dataset_paths"),
    "validate_mhc_protein_dataset": (".mhc_protein_dataset", "validate_mhc_protein_dataset"),
    # pipeline
    "FULL_FIELDS": (".pipeline", "FULL_FIELDS"),
    "RAW_FIELDS": (".pipeline", "RAW_FIELDS"),
    "build_full_seqs": (".pipeline", "build_full_seqs"),
    "build_raw_index": (".pipeline", "build_raw_index"),
    # species
    "LATIN_NAMES": (".species", "LATIN_NAMES"),
    "MHC_SPECIES_CATEGORIES": (".species", "MHC_SPECIES_CATEGORIES"),
    "full_species_name_alias": (".species", "full_species_name_alias"),
    "get_canonical_prefix": (".species", "get_canonical_prefix"),
    "get_latin_name": (".species", "get_latin_name"),
    "get_ncbi_taxid": (".species", "get_ncbi_taxid"),
    "get_species_clade": (".species", "get_species_clade"),
    "normalize_mhc_species": (".species", "normalize_mhc_species"),
    "normalize_species": (".species", "normalize_species"),
}

_LAZY_SUBMODULES: frozenset[str] = frozenset(
    {
        "alleles",
        "catalog",
        "cli_help",
        "datafiles",
        "domain_features",
        "domain_grammar",
        "domain_parsing",
        "domain_scoring",
        "domain_types",
        "download",
        "imgt",
        "mhc_protein_dataset",
        "pipeline",
        "signal_peptide",
        "species",
        "validate",
        "version",
    }
)


def _resolve(name: str) -> Any:
    """Resolve a lazily exported attribute, honoring monkeypatches in globals."""
    if name in globals():
        return globals()[name]
    return __getattr__(name)


def __getattr__(name: str) -> Any:
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        module = importlib.import_module(module_name, __name__)
        value = getattr(module, attr_name)
        globals()[name] = value
        return value
    if name in _LAZY_SUBMODULES:
        module = importlib.import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__) | set(_LAZY_SUBMODULES))


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
    d = _resolve("resolve_data_dir")()
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
    dd = _resolve("fasta_dir")(dd_root) if data_dir is None else Path(data_dir)

    paths = _resolve("download_all")(dd, force=force_download)
    sources = _resolve("SOURCES")
    fasta_inputs = [
        (paths["imgt_hla"], sources["imgt_hla"]["label"]),
        (paths["ipd_mhc"], sources["ipd_mhc"]["label"]),
    ]

    raw_csv = out / _resolve("RAW_CSV")
    _resolve("build_raw_index")(fasta_inputs, raw_csv)

    full_csv = out / _resolve("FULL_SEQS_CSV")
    report_path = out / _resolve("MERGE_REPORT")
    _resolve("build_full_seqs")(raw_csv, full_csv, report_path=report_path)

    warnings, stats = validate_build(raw_csv, full_csv)
    validation_path = out / _resolve("VALIDATION_REPORT")
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

    p = path or _resolve("find_csv")(_resolve("RAW_CSV"))
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
    FileNotFoundError
        If *path* is not given and no built CSV is found.
    """
    import pandas as pd

    p = path or _resolve("find_csv")(_resolve("RAW_CSV"))
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

    p = path or _resolve("find_csv")(_resolve("FULL_SEQS_CSV"))
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
    FileNotFoundError
        If *path* is not given and no built CSV is found.
    """
    import pandas as pd

    p = path or _resolve("find_csv")(_resolve("FULL_SEQS_CSV"))
    return pd.read_csv(p)


def lookup(
    allele: str,
    *,
    mutations: str | Sequence = (),
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
    from dataclasses import replace

    query = allele.strip()
    p = _resolve("find_csv")(_resolve("FULL_SEQS_CSV"), search_dir=search_dir)
    catalog = _resolve("load_catalog")(p)
    entry = catalog.get(query)
    if entry is None:
        raise KeyError(f"Allele {query!r} not found in {p}")

    result = _resolve("decompose_domains")(
        entry.sequence,
        mhc_class=entry.mhc_class,
        chain=entry.chain,
        allele=entry.two_field_allele,
        gene=entry.gene,
        species=entry.species,
        mutations=mutations,
    )

    final_seq = result.sequence if result.mutations and result.sequence else entry.sequence
    return replace(
        result,
        full_allele=entry.representative_allele,
        species=entry.species,
        species_category=entry.species_category,
        species_prefix=entry.species_prefix,
        source=entry.source,
        source_id=entry.source_id,
        sequence=final_seq,
        seq_len=len(final_seq),
    )


__all__ = [
    "__version__",
    "BuildPaths",
    "AlleleCatalog",
    "CatalogEntry",
    "CatalogError",
    "CatalogNotInstalledError",
    "CatalogPaths",
    "install_catalog",
    "load_catalog",
    "validate_catalog",
    "clear_catalog_cache",
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
    "get_species_clade",
    "get_ncbi_taxid",
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
    "align_domain_to_imgt",
    "align_groove_to_imgt",
    "map_groove_positions",
    "map_class_ii_groove_positions",
    "AlignedResidue",
    "DomainAlignment",
    "GrooveAlignment",
    "GrooveInsertion",
    "GrooveDeletion",
    "FIXED_FRAME_ALPHABET",
    "IMGT_GALPHA1_SLOTS",
    "IMGT_GALPHA2_SLOTS",
    "IMGT_CLASS_I_GROOVE_SLOTS",
    "IMGT_CLASS_I_PROFILE_PATH",
    "IMGT_GALPHA_SLOTS",
    "IMGT_GBETA_SLOTS",
    "IMGT_CLASS_II_GROOVE_SLOTS",
    "IMGT_CLASS_II_PROFILE_PATH",
    "IMGT_PAN_HALF1_SLOTS",
    "IMGT_PAN_HALF2_SLOTS",
    "IMGT_PAN_CLASS_GROOVE_SLOTS",
    "STRUCTURAL_ELEMENTS",
    "GALPHA1_POSITIONS",
    "GALPHA2_GAP_POSITIONS",
    "GALPHA2_POSITIONS",
    "GALPHA_GAP_POSITIONS",
    "GALPHA_POSITIONS",
    "GBETA_GAP_POSITIONS",
    "GBETA_POSITIONS",
    "HELIX_INSERTIONS",
    "HLA_A0201_GALPHA1",
    "HLA_A0201_GALPHA2",
    "HLA_DRA0101_GALPHA",
    "HLA_DRB1_0101_GBETA",
    "CONSERVED_CYS_POSITIONS",
]
