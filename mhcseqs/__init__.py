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
    AlleleRecord,
    RawAllele,
    apply_mutations,
    extract_groove,
    find_cys_pairs,
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


def lookup(allele: str, *, search_dir: str | None = None) -> AlleleRecord:
    """Look up a single allele and return an AlleleRecord.

    Merges data from both the full-sequences and binding-grooves CSVs,
    giving you everything: full sequence, signal peptide boundary,
    groove domains, Ig domain, tail, species info, etc.

    >>> result = mhcseqs.lookup("HLA-A*02:01")  # doctest: +SKIP
    >>> result.groove1[:10]                      # doctest: +SKIP
    'GSHSMRYFFT'
    >>> result.mature_sequence[:10]              # doctest: +SKIP
    'GSHSMRYFFT'

    Raises:
        FileNotFoundError: If the CSVs haven't been built yet.
        KeyError: If no matching allele is found.
    """
    import csv as _csv

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

    merged: dict = {}

    # Try full-seqs first (has sequence, mature_start, etc.)
    try:
        p = _find_csv("mhc-full-seqs.csv", search_dir=search_dir)
        with open(p, "r", encoding="utf-8") as f:
            for row in _csv.DictReader(f):
                if _match(row):
                    merged.update(row)
                    break
    except FileNotFoundError:
        pass

    # Then grooves (has groove1, groove2, ig_domain, tail, etc.)
    try:
        p = _find_csv("mhc-binding-grooves.csv", search_dir=search_dir)
        with open(p, "r", encoding="utf-8") as f:
            for row in _csv.DictReader(f):
                if _match(row):
                    merged.update(row)
                    break
    except FileNotFoundError:
        pass

    if not merged:
        try:
            _find_csv("mhc-binding-grooves.csv", search_dir=search_dir)
        except FileNotFoundError:
            raise FileNotFoundError("No built CSVs found. Run mhcseqs.build() or 'mhcseqs build' first.") from None
        raise KeyError(f"Allele {query!r} not found")

    return _dict_to_allele_record(merged)


def _dict_to_allele_record(d: dict) -> AlleleRecord:
    """Convert a merged CSV row dict to an AlleleRecord."""

    def _int(key: str, default: int = 0) -> int:
        v = d.get(key, "")
        if not v:
            return default
        try:
            return int(v)
        except (ValueError, TypeError):
            return default

    def _opt_int(key: str) -> int | None:
        v = d.get(key, "")
        if not v:
            return None
        try:
            return int(v)
        except (ValueError, TypeError):
            return None

    return AlleleRecord(
        allele=d.get("two_field_allele", d.get("allele", "")),
        representative_allele=d.get("representative_allele", ""),
        gene=d.get("gene", ""),
        mhc_class=d.get("mhc_class", ""),
        chain=d.get("chain", ""),
        species=d.get("species", ""),
        species_category=d.get("species_category", ""),
        species_prefix=d.get("species_prefix", ""),
        source=d.get("source", ""),
        sequence=d.get("sequence", ""),
        seq_len=_int("seq_len"),
        mature_start=_int("mature_start"),
        groove_seq=d.get("groove_seq", ""),
        groove1=d.get("groove1", ""),
        groove2=d.get("groove2", ""),
        groove1_len=_int("groove1_len"),
        groove2_len=_int("groove2_len"),
        ig_domain=d.get("ig_domain", ""),
        ig_domain_len=_int("ig_domain_len"),
        tail=d.get("tail", ""),
        tail_len=_int("tail_len"),
        status=d.get("groove_status", d.get("status", "ok")),
        anchor_type=d.get("anchor_type", ""),
        anchor_cys1=_opt_int("anchor_cys1"),
        anchor_cys2=_opt_int("anchor_cys2"),
        secondary_cys1=_opt_int("secondary_cys1"),
        secondary_cys2=_opt_int("secondary_cys2"),
    )


__all__ = [
    "__version__",
    "build",
    "load_raw",
    "load_full_seqs",
    "load_grooves",
    "lookup",
    "build_raw_index",
    "build_full_seqs",
    "build_binding_grooves",
    "apply_mutations",
    "extract_groove",
    "parse_class_i",
    "parse_class_ii_alpha",
    "parse_class_ii_beta",
    "find_cys_pairs",
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
    "GROOVE_FIELDS",
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
