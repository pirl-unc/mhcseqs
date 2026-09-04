"""CLI entry point for mhcseqs.

Usage:
    mhcseqs build [--output-dir DIR] [--data-dir DIR] [--force-download]
    mhcseqs lookup ALLELE [--output-dir DIR]
    mhcseqs data list
    mhcseqs data install mhc-proteins [--version VERSION] [--with-sources]
    mhcseqs data refresh [--source imgt_hla|ipd_mhc]
    mhcseqs data clear [--source ...] [--built | --built-only]
    python -m mhcseqs build
    python -m mhcseqs lookup "HLA-A*02:01"
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from . import default_data_dir
from .datafiles import clear as clear_data
from .datafiles import fasta_dir, inventory
from .download import SOURCES, download_all
from .mhc_protein_dataset import (
    DATASET_NAME as MHC_PROTEIN_DATASET_NAME,
)
from .mhc_protein_dataset import (
    ProteinDatasetError,
    available_mhc_protein_dataset_versions,
    default_mhc_protein_dataset_version,
    install_mhc_protein_dataset,
    install_mhc_protein_source_bundle,
    mhc_protein_dataset_paths,
    validate_mhc_protein_dataset,
)
from .pipeline import build_full_seqs, build_raw_index
from .validate import format_validation_report, validate_build
from .version import __version__


def cmd_build(args):
    dd = Path(default_data_dir())
    out_dir: Path = args.output_dir or dd
    data_dir: Path = args.data_dir or (dd / "fasta")

    print("=" * 60)
    print("Step 1/4: Downloading FASTA source files")
    print("=" * 60)
    paths = download_all(data_dir, force=args.force_download)
    print()

    fasta_inputs = [
        (paths["imgt_hla"], SOURCES["imgt_hla"]["label"]),
        (paths["ipd_mhc"], SOURCES["ipd_mhc"]["label"]),
    ]

    raw_csv = out_dir / "mhc-seqs-raw.csv"
    print("=" * 60)
    print("Step 2/4: Building raw sequence index")
    print("=" * 60)
    raw_stats = build_raw_index(fasta_inputs, raw_csv)
    print(f"  Raw index: {json.dumps(raw_stats, indent=2)}")
    print(f"  -> {raw_csv}")
    print()

    full_csv = out_dir / "mhc-full-seqs.csv"
    report_path = out_dir / "mhc-merge-report.txt"
    print("=" * 60)
    print("Step 3/4: Selecting two-field representatives + groove extraction")
    print("=" * 60)
    full_stats = build_full_seqs(raw_csv, full_csv, report_path=report_path)
    print(f"  Full seqs: {json.dumps(full_stats, indent=2)}")
    print(f"  -> {full_csv}")
    print(f"  -> {report_path}")
    print()

    print("=" * 60)
    print("Step 4/4: Validation sanity checks")
    print("=" * 60)
    warnings, val_stats = validate_build(raw_csv, full_csv)
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
    print(f"  {report_path}")
    print(f"  {validation_path}")


def cmd_lookup(args):
    dd = Path(default_data_dir())
    out_dir: Path = args.output_dir or dd

    # Find the CSVs — look in output dir, default data dir, and CWD
    search_dirs = [out_dir, dd, Path(".")]
    csv_files = {}
    for name in ("mhc-seqs-raw.csv", "mhc-full-seqs.csv"):
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


def _human_size(n: int) -> str:
    size = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} GB"


def _human_age(mtime: float) -> str:
    import time

    secs = max(0.0, time.time() - mtime)
    days = secs / 86400
    if days >= 1:
        return f"{days:.0f}d ago"
    hours = secs / 3600
    if hours >= 1:
        return f"{hours:.0f}h ago"
    return f"{secs / 60:.0f}m ago"


def cmd_data(args):
    dd = Path(default_data_dir())
    action = args.data_command

    if action == "list":
        print(f"Data directory: {dd}")
        print(f"{'FILE':<26} {'KIND':<7} {'SIZE':>10} {'AGE':>9}  STATUS")
        for item in inventory(dd):
            if item.exists:
                status = item.url or "built"
                print(f"{item.path.name:<26} {item.kind:<7} {_human_size(item.size):>10} {_human_age(item.mtime):>9}  {status}")
            else:
                print(f"{item.path.name:<26} {item.kind:<7} {'—':>10} {'—':>9}  (not present)")
        print("\nVersioned datasets:")
        for version in available_mhc_protein_dataset_versions():
            paths = mhc_protein_dataset_paths(version, data_dir=dd)
            if paths.records.parent.exists():
                try:
                    validate_mhc_protein_dataset(version, data_dir=dd)
                    status = _human_size(paths.records.stat().st_size)
                except ProteinDatasetError:
                    status = "invalid; reinstall with --force"
            else:
                status = "not installed"
            print(f"  {MHC_PROTEIN_DATASET_NAME} {version}: {status}")
        return

    if action == "install":
        paths = install_mhc_protein_dataset(args.version, data_dir=dd, force=args.force)
        print(f"Installed {MHC_PROTEIN_DATASET_NAME} {paths.version}")
        print(f"  records: {paths.records}")
        print(f"  manifest: {paths.manifest}")
        if args.with_sources:
            sources = install_mhc_protein_source_bundle(args.version, data_dir=dd, force=args.force)
            print(f"  source bundle: {sources.root}")
            print(f"  label curation: {sources.label_curation}")
        return

    if action == "path":
        try:
            paths = validate_mhc_protein_dataset(args.version, data_dir=dd)
        except ProteinDatasetError as exc:
            version = args.version or default_mhc_protein_dataset_version()
            print(
                f"{MHC_PROTEIN_DATASET_NAME} {version} is missing or corrupt: {exc}\n"
                f"Reinstall it with 'mhcseqs data install {MHC_PROTEIN_DATASET_NAME} --version {version} --force'.",
                file=sys.stderr,
            )
            sys.exit(1)
        print(paths.records)
        return

    if action == "refresh":
        from .download import download_fasta

        keys = [args.source] if args.source else list(SOURCES)
        print(f"Refreshing source FASTA(s): {', '.join(keys)}")
        for k in keys:
            download_fasta(k, fasta_dir(dd), force=True)
        print("Done.")
        return

    if action == "clear":
        # Default: clear source FASTAs. --built also/only clears built CSVs.
        keys = [args.source] if args.source else None
        removed = clear_data(
            dd,
            sources=not args.built_only,
            built=args.built or args.built_only,
            source_keys=keys,
        )
        if removed:
            for p in removed:
                print(f"Removed {p}")
            print(f"Cleared {len(removed)} file(s).")
        else:
            print("Nothing to clear.")
        return


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
        default=None,
        help="Directory for output CSVs (default: ~/.cache/mhcseqs)",
    )
    build_parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Directory for downloaded FASTA files (default: ~/.cache/mhcseqs/fasta)",
    )
    build_parser.add_argument(
        "--force-download",
        action="store_true",
        help="Re-download source FASTA files even when a cached copy exists",
    )

    # lookup
    lookup_parser = subparsers.add_parser("lookup", help="Look up an allele from built CSVs")
    lookup_parser.add_argument("allele", help="Allele name to look up (e.g. HLA-A*02:01)")
    lookup_parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory where built CSVs are located (default: ~/.cache/mhcseqs)",
    )

    # data — manage cached source FASTAs and built CSVs
    data_parser = subparsers.add_parser("data", help="Manage cached downloads and built CSVs")
    data_sub = data_parser.add_subparsers(dest="data_command")
    data_sub.add_parser("list", help="List cached source and built data files")
    install_parser = data_sub.add_parser("install", help="Install an immutable versioned dataset")
    install_parser.add_argument("dataset", choices=[MHC_PROTEIN_DATASET_NAME])
    install_parser.add_argument("--version", choices=available_mhc_protein_dataset_versions())
    install_parser.add_argument("--force", action="store_true", help="Re-download and replace this exact cached version")
    install_parser.add_argument(
        "--with-sources",
        action="store_true",
        help="Also install the pinned source and curation artifacts needed for offline regeneration",
    )
    path_parser = data_sub.add_parser("path", help="Print the path to an installed versioned dataset")
    path_parser.add_argument("dataset", choices=[MHC_PROTEIN_DATASET_NAME])
    path_parser.add_argument("--version", choices=available_mhc_protein_dataset_versions())
    refresh_parser = data_sub.add_parser("refresh", help="Re-download source FASTA(s)")
    refresh_parser.add_argument(
        "--source",
        choices=list(SOURCES),
        default=None,
        help="Only refresh this source (default: all)",
    )
    clear_parser = data_sub.add_parser("clear", help="Delete cached source FASTA(s)")
    clear_parser.add_argument(
        "--source",
        choices=list(SOURCES),
        default=None,
        help="Only clear this source FASTA (default: all sources)",
    )
    clear_parser.add_argument(
        "--built",
        action="store_true",
        help="Also delete built CSVs/reports (in addition to source FASTAs)",
    )
    clear_parser.add_argument(
        "--built-only",
        action="store_true",
        help="Delete only built CSVs/reports, keep source FASTAs",
    )

    args = parser.parse_args()

    if args.command == "build":
        cmd_build(args)
    elif args.command == "lookup":
        cmd_lookup(args)
    elif args.command == "data":
        if not args.data_command:
            data_parser.print_help()
            sys.exit(1)
        cmd_data(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
