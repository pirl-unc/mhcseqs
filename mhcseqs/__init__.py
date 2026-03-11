"""MHC sequence curation and binding groove extraction.

Quick start::

    import mhcseqs

    # Build the database (downloads FASTA sources, only needed once)
    mhcseqs.build()

    # Look up any allele → AlleleRecord
    r = mhcseqs.lookup("HLA-A*02:01")
    r.groove1            # α1 domain
    r.groove2            # α2 domain
    r.mature_sequence    # signal peptide removed
    r.sequence           # full protein (with signal peptide)

    # Look up with mutations (IEDB-style)
    m = mhcseqs.lookup("HLA-A*02:01", mutations=["K66A"])

    # Extract groove from a raw sequence (no build needed)
    r = mhcseqs.extract_groove(seq, mhc_class="I")
"""

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


def lookup(
    allele: str,
    *,
    mutations: Sequence = (),
    search_dir: str | None = None,
) -> AlleleRecord:
    """Look up a single allele and return a fully parsed AlleleRecord.

    Finds the allele in the built full-sequences CSV, then runs
    ``extract_groove()`` on it to produce a complete record with all
    fields populated (including anchor Cys positions).

    Optionally applies mutations (same format as ``extract_groove()``).

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

    result = extract_groove(
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
        sequence=seq,
        seq_len=len(seq),
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
