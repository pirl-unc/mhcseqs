"""MHC sequence curation and binding groove extraction.

Quick start::

    # Build all CSVs from upstream databases
    import mhcseqs
    paths = mhcseqs.build()

    # Load results as DataFrames (requires pandas)
    import pandas as pd
    df = pd.read_csv(paths["grooves"])

    # Or use the convenience loaders (returns list[dict], no pandas needed)
    rows = mhcseqs.load_grooves()
"""

from .alleles import (
    infer_gene,
    infer_mhc_class,
    infer_species,
    normalize_allele_name,
    normalize_mhc_class,
    parse_allele_name,
)
from .download import SOURCES, download_all
from .groove import (
    NON_GROOVE_GENES,
    GrooveResult,
    extract_groove,
    find_cys_pairs,
    parse_class_i,
    parse_class_ii_alpha,
    parse_class_ii_beta,
)
from .pipeline import (
    FULL_FIELDS,
    GROOVE_FIELDS,
    RAW_FIELDS,
    build_binding_grooves,
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


def build(
    output_dir: str = ".",
    data_dir: str | None = None,
) -> dict[str, str]:
    """Run the full build pipeline: download, parse, extract grooves.

    Returns a dict mapping output names to file paths::

        {
            "raw": "path/to/mhc-seqs-raw.csv",
            "full_seqs": "path/to/mhc-full-seqs.csv",
            "grooves": "path/to/mhc-binding-grooves.csv",
            "merge_report": "path/to/mhc-merge-report.txt",
            "validation_report": "path/to/mhc-validation-report.txt",
        }
    """
    from pathlib import Path

    from .validate import format_validation_report, validate_build

    out = Path(output_dir)
    if data_dir is None:
        dd = Path(__file__).resolve().parent.parent / "data" / "fasta"
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

    groove_csv = out / "mhc-binding-grooves.csv"
    build_binding_grooves(full_csv, groove_csv)

    warnings, stats = validate_build(raw_csv, full_csv, groove_csv)
    validation_path = out / "mhc-validation-report.txt"
    with open(validation_path, "w", encoding="utf-8") as f:
        f.write(format_validation_report(warnings, stats) + "\n")

    return {
        "raw": str(raw_csv),
        "full_seqs": str(full_csv),
        "grooves": str(groove_csv),
        "merge_report": str(report_path),
        "validation_report": str(validation_path),
    }


def _find_csv(name: str, search_dir: str | None = None) -> str:
    """Locate a built CSV file."""
    from pathlib import Path

    candidates = []
    if search_dir:
        candidates.append(Path(search_dir) / name)
    candidates.extend(
        [
            Path(".") / name,
            Path(__file__).resolve().parent.parent / name,
        ]
    )
    for p in candidates:
        if p.exists():
            return str(p)
    raise FileNotFoundError(f"{name} not found. Run mhcseqs.build() or 'mhcseqs build' first.")


def load_raw(path: str | None = None) -> list[dict]:
    """Load the raw sequence CSV as a list of dicts."""
    import csv as _csv

    p = path or _find_csv("mhc-seqs-raw.csv")
    with open(p, "r", encoding="utf-8") as f:
        return list(_csv.DictReader(f))


def load_full_seqs(path: str | None = None) -> list[dict]:
    """Load the full-sequence representatives CSV as a list of dicts."""
    import csv as _csv

    p = path or _find_csv("mhc-full-seqs.csv")
    with open(p, "r", encoding="utf-8") as f:
        return list(_csv.DictReader(f))


def load_grooves(path: str | None = None) -> list[dict]:
    """Load the binding grooves CSV as a list of dicts."""
    import csv as _csv

    p = path or _find_csv("mhc-binding-grooves.csv")
    with open(p, "r", encoding="utf-8") as f:
        return list(_csv.DictReader(f))


__all__ = [
    "__version__",
    "build",
    "load_raw",
    "load_full_seqs",
    "load_grooves",
    "build_raw_index",
    "build_full_seqs",
    "build_binding_grooves",
    "extract_groove",
    "parse_class_i",
    "parse_class_ii_alpha",
    "parse_class_ii_beta",
    "find_cys_pairs",
    "GrooveResult",
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
    "GROOVE_FIELDS",
]
