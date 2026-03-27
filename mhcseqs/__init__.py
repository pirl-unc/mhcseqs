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
from .download import SOURCES, download_all
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
from .pipeline import (
    FULL_FIELDS,
    RAW_FIELDS,
    build_full_seqs,
    build_raw_index,
)
from .species import (
    CANONICAL_MHC_PREFIXES,
    LATIN_NAMES,
    MHC_SPECIES_CATEGORIES,
    get_canonical_prefix,
    get_latin_name,
    normalize_mhc_species,
    normalize_species,
)
from .version import __version__


@dataclass(frozen=True)
class BuildPaths:
    """Paths to CSV files produced by :func:`build`.

    Attributes:
        raw: Every protein entry from all sources (``mhc-seqs-raw.csv``).
        sequences: One representative per two-field allele group, with full
            sequence, groove decomposition, and metadata (``mhc-full-seqs.csv``).
        merge_report: Human-readable report of the merge/dedup step.
        validation_report: Human-readable validation sanity checks.
    """

    raw: str
    sequences: str
    merge_report: str
    validation_report: str


def default_data_dir() -> str:
    """Return the default directory for downloads and built CSVs.

    Uses ``$MHCSEQS_DATA`` if set, otherwise ``~/.cache/mhcseqs``.
    The directory is created on first use.
    """
    import os
    from pathlib import Path

    d = Path(os.environ.get("MHCSEQS_DATA", Path.home() / ".cache" / "mhcseqs"))
    d.mkdir(parents=True, exist_ok=True)
    return str(d)


def build(
    output_dir: str | None = None,
    data_dir: str | None = None,
) -> BuildPaths:
    """Run the full build pipeline: download, parse, extract grooves.

    If *output_dir* is ``None`` (default), CSVs are written to
    :func:`default_data_dir` (``~/.cache/mhcseqs``).
    FASTA downloads are stored in a ``fasta/`` subdirectory of the data dir.

    Returns a :class:`BuildPaths` with the paths to the generated files.
    """
    from pathlib import Path

    from .validate import format_validation_report, validate_build

    if output_dir is None:
        out = Path(default_data_dir())
    else:
        out = Path(output_dir)
    if data_dir is None:
        dd = out / "fasta"
    else:
        dd = Path(data_dir)

    paths = download_all(dd)
    fasta_inputs = [
        (paths["imgt_hla"], SOURCES["imgt_hla"]["label"]),
        (paths["ipd_mhc"], SOURCES["ipd_mhc"]["label"]),
    ]

    raw_csv = out / "mhc-seqs-raw.csv"
    build_raw_index(fasta_inputs, raw_csv)

    full_csv = out / "mhc-full-seqs.csv"
    report_path = out / "mhc-merge-report.txt"
    build_full_seqs(raw_csv, full_csv, report_path=report_path)

    warnings, stats = validate_build(raw_csv, full_csv)
    validation_path = out / "mhc-validation-report.txt"
    with open(validation_path, "w", encoding="utf-8") as f:
        f.write(format_validation_report(warnings, stats) + "\n")

    return BuildPaths(
        raw=str(raw_csv),
        sequences=str(full_csv),
        merge_report=str(report_path),
        validation_report=str(validation_path),
    )


def _find_csv(name: str, search_dir: str | None = None) -> str:
    """Locate a built CSV file.

    Search order: *search_dir* (if given), current directory,
    :func:`default_data_dir`, package parent directory.
    """
    from pathlib import Path

    candidates = []
    if search_dir:
        candidates.append(Path(search_dir) / name)
    candidates.extend(
        [
            Path(".") / name,
            Path(default_data_dir()) / name,
            Path(__file__).resolve().parent.parent / name,
        ]
    )
    for p in candidates:
        if p.exists():
            return str(p)
    raise FileNotFoundError(f"{name} not found. Run mhcseqs.build() or 'mhcseqs build' first.")


def load_raw_dict(path: str | None = None) -> list[dict]:
    """Load the raw sequence CSV as a list of dicts.

    Contains every individual protein entry from all sources,
    including duplicates and fragments.
    """
    import csv as _csv

    p = path or _find_csv("mhc-seqs-raw.csv")
    with open(p, "r", encoding="utf-8") as f:
        return list(_csv.DictReader(f))


def load_raw_dataframe(path: str | None = None):
    """Load the raw sequence CSV as a pandas DataFrame."""
    import pandas as pd

    p = path or _find_csv("mhc-seqs-raw.csv")
    return pd.read_csv(p)


def load_sequences_dict(path: str | None = None) -> list[dict]:
    """Load the full-sequence representatives CSV as a list of dicts.

    One row per two-field allele with the best available full-length
    protein sequence, groove decomposition, species metadata, and
    signal peptide annotation.
    """
    import csv as _csv

    p = path or _find_csv("mhc-full-seqs.csv")
    with open(p, "r", encoding="utf-8") as f:
        return list(_csv.DictReader(f))


def load_sequences_dataframe(path: str | None = None):
    """Load the full-sequence representatives CSV as a pandas DataFrame.

    One row per two-field allele with full protein, mature sequence,
    groove decomposition, species metadata, and allele status flags.
    """
    import pandas as pd

    p = path or _find_csv("mhc-full-seqs.csv")
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

    >>> r = mhcseqs.lookup("HLA-A*02:01")           # doctest: +SKIP
    >>> r.groove1[:10]                               # doctest: +SKIP
    'GSHSMRYFFT'
    >>> m = mhcseqs.lookup("HLA-A*02:01",            # doctest: +SKIP
    ...                    mutations=["K66A"])
    >>> m.groove1[65]                                # doctest: +SKIP
    'A'

    Raises:
        FileNotFoundError: If the full-seqs CSV hasn't been built yet.
        KeyError: If no matching allele is found.
    """
    import csv as _csv
    from dataclasses import replace

    query = allele.strip()
    query_upper = query.upper()
    try:
        query_normalized = normalize_allele_name(query)
    except Exception:
        query_normalized = query

    def _match(row: dict) -> bool:
        for f in ("two_field_allele", "representative_allele"):
            val = row.get(f, "").strip()
            if val and (val == query or val.upper() == query_upper or val == query_normalized):
                return True
        return False

    p = _find_csv("mhc-full-seqs.csv", search_dir=search_dir)
    hit: dict | None = None
    with open(p, "r", encoding="utf-8") as f:
        for row in _csv.DictReader(f):
            if _match(row):
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
    "build",
    "load_raw_dict",
    "load_raw_dataframe",
    "load_sequences_dict",
    "load_sequences_dataframe",
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
    "RawAllele",
    "AlleleRecord",
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
    "get_canonical_prefix",
    "MHC_SPECIES_CATEGORIES",
    "LATIN_NAMES",
    "CANONICAL_MHC_PREFIXES",
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
