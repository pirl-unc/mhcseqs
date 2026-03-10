"""CLI entry point for mhcseqs.

Usage:
    mhcseqs build [--output-dir DIR] [--data-dir DIR]
    mhcseqs lookup ALLELE [--data-dir DIR]
    python -m mhcseqs build
    python -m mhcseqs lookup "HLA-A*02:01"
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from .download import SOURCES, download_all
from .pipeline import build_binding_grooves, build_full_seqs, build_raw_index
from .validate import format_validation_report, validate_build
from .version import __version__


def _default_data_dir() -> Path:
    return Path(__file__).resolve().parent.parent / "data" / "fasta"


def _default_output_dir() -> Path:
    return Path(__file__).resolve().parent.parent


def cmd_build(args):
    out_dir: Path = args.output_dir
    data_dir: Path = args.data_dir

    print("=" * 60)
    print("Step 1/5: Downloading FASTA source files")
    print("=" * 60)
    paths = download_all(data_dir)
    print()

    fasta_inputs = [
        (paths["imgt_hla"], SOURCES["imgt_hla"]["label"]),
        (paths["ipd_mhc"], SOURCES["ipd_mhc"]["label"]),
    ]

    raw_csv = out_dir / "mhc-seqs-raw.csv"
    print("=" * 60)
    print("Step 2/5: Building raw sequence index")
    print("=" * 60)
    raw_stats = build_raw_index(fasta_inputs, raw_csv)
    print(f"  Raw index: {json.dumps(raw_stats, indent=2)}")
    print(f"  -> {raw_csv}")
    print()

    full_csv = out_dir / "mhc-full-seqs.csv"
    report_path = out_dir / "mhc-merge-report.txt"
    print("=" * 60)
    print("Step 3/5: Selecting two-field representatives")
    print("=" * 60)
    full_stats = build_full_seqs(raw_csv, full_csv, report_path=report_path)
    print(f"  Full seqs: {json.dumps(full_stats, indent=2)}")
    print(f"  -> {full_csv}")
    print(f"  -> {report_path}")
    print()

    groove_csv = out_dir / "mhc-binding-grooves.csv"
    print("=" * 60)
    print("Step 4/5: Extracting binding grooves")
    print("=" * 60)
    groove_stats = build_binding_grooves(full_csv, groove_csv)
    print(f"  Grooves: {json.dumps(groove_stats, indent=2)}")
    print(f"  -> {groove_csv}")
    print()

    print("=" * 60)
    print("Step 5/5: Validation sanity checks")
    print("=" * 60)
    warnings, val_stats = validate_build(raw_csv, full_csv, groove_csv)
    report = format_validation_report(warnings, val_stats)
    print(report)
    validation_path = out_dir / "mhc-validation-report.txt"
    with open(validation_path, "w", encoding="utf-8") as f:
        f.write(report + "\n")
    print(f"  -> {validation_path}")
    print()

    print("=" * 60)
    print("Done!")
    print(f"  {raw_csv}")
    print(f"  {full_csv}")
    print(f"  {groove_csv}")
    print(f"  {report_path}")
    print(f"  {validation_path}")


def cmd_lookup(args):
    data_dir: Path = args.data_dir

    # Find the CSVs — look in data_dir's parent (the repo root) by default
    search_dirs = [data_dir.parent, data_dir, Path(".")]
    csv_files = {}
    for name in ("mhc-seqs-raw.csv", "mhc-full-seqs.csv", "mhc-binding-grooves.csv"):
        for d in search_dirs:
            p = d / name
            if p.exists():
                csv_files[name] = p
                break

    if not csv_files:
        print("No built CSVs found. Run 'mhcseqs build' first.", file=sys.stderr)
        sys.exit(1)

    query = args.allele.strip()
    query_upper = query.upper()

    # Normalize query for matching
    try:
        from .alleles import normalize_allele_name

        query_normalized = normalize_allele_name(query)
    except Exception:
        query_normalized = query

    def _matches(row: dict) -> bool:
        for field in (
            "two_field_allele",
            "representative_allele",
            "allele_normalized",
            "allele_raw",
        ):
            val = row.get(field, "").strip()
            if val and (val == query or val.upper() == query_upper or val == query_normalized):
                return True
        return False

    found_any = False
    for csv_name, csv_path in sorted(csv_files.items()):
        with open(csv_path, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if _matches(row):
                    if not found_any:
                        found_any = True
                    print(f"\n--- {csv_name} ---")
                    for k, v in row.items():
                        if k == "sequence":
                            print(f"  {k}: {v[:50]}... ({len(v)} aa)")
                        elif k == "mature_sequence":
                            print(f"  {k}: {v[:50]}... ({len(v)} aa)")
                        else:
                            print(f"  {k}: {v}")
                    break  # one match per CSV is enough

    if not found_any:
        print(f"No match found for '{query}'")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        prog="mhcseqs",
        description="MHC sequence curation and binding groove extraction",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    subparsers = parser.add_subparsers(dest="command")

    # build
    build_parser = subparsers.add_parser("build", help="Download sources and build CSVs")
    build_parser.add_argument(
        "--output-dir",
        type=Path,
        default=_default_output_dir(),
        help="Directory for output CSVs",
    )
    build_parser.add_argument(
        "--data-dir",
        type=Path,
        default=_default_data_dir(),
        help="Directory for downloaded FASTA files",
    )

    # lookup
    lookup_parser = subparsers.add_parser("lookup", help="Look up an allele from built CSVs")
    lookup_parser.add_argument("allele", help="Allele name to look up (e.g. HLA-A*02:01)")
    lookup_parser.add_argument(
        "--data-dir",
        type=Path,
        default=_default_data_dir(),
        help="Directory containing built CSVs",
    )

    args = parser.parse_args()

    if args.command == "build":
        cmd_build(args)
    elif args.command == "lookup":
        cmd_lookup(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
