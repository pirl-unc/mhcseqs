"""CLI entry point for mhcseqs.

Usage:
    mhcseqs build [--output-dir DIR] [--data-dir DIR] [--force-download]
    mhcseqs lookup ALLELE [--output-dir DIR] [--format {text,csv,tsv,json,fasta}] [--all] [--output PATH]
    mhcseqs data list [--details]
    mhcseqs data install mhc-proteins [--version VERSION] [--with-sources] [--force]
    mhcseqs data path mhc-proteins [--version VERSION]
    mhcseqs data refresh [--source imgt_hla|ipd_mhc]
    mhcseqs data clear [--source ...] [--built | --built-only]
    python -m mhcseqs build
    python -m mhcseqs lookup "HLA-A*02:01"
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import sys
from contextlib import contextmanager
from pathlib import Path

from . import default_data_dir
from . import default_data_dir as _default_data_dir_impl
from .datafiles import (
    FULL_SEQS_CSV,
    MERGE_REPORT,
    RAW_CSV,
    VALIDATION_REPORT,
    build_query_candidates,
    fasta_dir,
    find_csv,
    inventory,
    match_allele_row,
    resolve_data_dir,
    resolve_data_location,
)
from .datafiles import clear as clear_data
from .download import SOURCES, download_all
from .mhc_protein_dataset import (
    DATASET_NAME as MHC_PROTEIN_DATASET_NAME,
)
from .mhc_protein_dataset import (
    ProteinDatasetError,
    ProteinDatasetNotInstalledError,
    available_mhc_protein_dataset_versions,
    default_mhc_protein_dataset_version,
    install_mhc_protein_dataset,
    install_mhc_protein_source_bundle,
    validate_mhc_protein_dataset,
)
from .pipeline import build_full_seqs, build_raw_index
from .validate import format_validation_report, validate_build
from .version import __version__


def cmd_build(args):
    dd = Path(default_data_dir())
    out_dir: Path = args.output_dir or dd
    data_dir: Path = args.data_dir or fasta_dir(dd)

    print("=" * 60)
    print("Step 1/4: Downloading FASTA source files")
    print("=" * 60)
    paths = download_all(data_dir, force=args.force_download)
    print()

    fasta_inputs = [
        (paths["imgt_hla"], SOURCES["imgt_hla"]["label"]),
        (paths["ipd_mhc"], SOURCES["ipd_mhc"]["label"]),
    ]

    raw_csv = out_dir / RAW_CSV
    print("=" * 60)
    print("Step 2/4: Building raw sequence index")
    print("=" * 60)
    raw_stats = build_raw_index(fasta_inputs, raw_csv)
    print(f"  Raw index: {json.dumps(raw_stats, indent=2)}")
    print(f"  -> {raw_csv}")
    print()

    full_csv = out_dir / FULL_SEQS_CSV
    report_path = out_dir / MERGE_REPORT
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
    validation_path = out_dir / VALIDATION_REPORT
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


LOOKUP_FORMATS = ("text", "csv", "tsv", "json", "fasta")


def format_lookup_results(
    matches: list[tuple[str, dict[str, str]]],
    fmt: str = "text",
) -> str:
    """Format matched lookup rows in the requested output format.

    All formats emit full sequences without truncation.
    """
    if fmt == "text":
        blocks = []
        for csv_name, row in matches:
            lines = [f"--- {csv_name} ---"]
            for k, v in row.items():
                lines.append(f"  {k}: {v}")
            blocks.append("\n".join(lines))
        return "\n\n".join(blocks) + "\n"

    if fmt == "json":
        payload = [{"csv": csv_name, **row} for csv_name, row in matches]
        return json.dumps(payload, indent=2) + "\n"

    if fmt in ("csv", "tsv"):
        buf = io.StringIO()
        fieldnames = ["csv"]
        seen_fields = {"csv"}
        for _, row in matches:
            for key in row:
                if key not in seen_fields:
                    seen_fields.add(key)
                    fieldnames.append(key)
        delimiter = "\t" if fmt == "tsv" else ","
        writer = csv.DictWriter(
            buf,
            fieldnames=fieldnames,
            delimiter=delimiter,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for csv_name, row in matches:
            writer.writerow({"csv": csv_name, **row})
        return buf.getvalue()

    if fmt == "fasta":
        records = []
        for csv_name, row in matches:
            allele_id = (
                row.get("representative_allele") or row.get("allele_normalized") or row.get("two_field_allele") or row.get("allele_raw") or "unknown"
            )
            meta = [f"csv={csv_name}"]
            for key in ("gene", "mhc_class", "chain", "species"):
                val = (row.get(key) or "").strip()
                if val:
                    meta.append(f"{key}={val}")
            seq = (row.get("sequence") or "").strip()
            records.append(f">{allele_id} {' '.join(meta)}\n{seq}")
        return "\n".join(records) + "\n"

    raise ValueError(f"Unsupported lookup format: {fmt!r}")


def cmd_lookup(args):
    dd = _cli_data_dir(create=False)

    # Locate built CSVs using the same search order as mhcseqs.lookup()
    csv_files: dict[str, Path] = {}
    for name in (FULL_SEQS_CSV, RAW_CSV):
        try:
            csv_files[name] = find_csv(name, search_dir=args.output_dir, data_dir=dd)
        except FileNotFoundError:
            continue

    if not csv_files:
        print("No built CSVs found. Run 'mhcseqs build' first.", file=sys.stderr)
        sys.exit(1)

    query = args.allele.strip()
    query_candidates = build_query_candidates(query)

    matches: list[tuple[str, dict[str, str]]] = []
    for csv_name, csv_path in csv_files.items():
        with open(csv_path, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if match_allele_row(row, query, query_candidates=query_candidates):
                    matches.append((csv_name, row))
                    if not getattr(args, "all", False):
                        break

    if not matches:
        print(f"No match found for '{query}'", file=sys.stderr)
        sys.exit(1)

    rendered = format_lookup_results(matches, fmt=getattr(args, "format", "text"))
    output_path: Path | None = getattr(args, "output", None)
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(rendered, encoding="utf-8")
    else:
        sys.stdout.write(rendered)


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


def _print_file_group(title: str, items, *, details: bool) -> None:
    print(f"\n{title}")
    rows = [
        (
            item.path.name,
            _human_size(item.size) if item.exists else "—",
            _human_age(item.mtime) if item.exists else "—",
            "present" if item.exists else "missing",
        )
        for item in items
    ]
    headers = ("FILE", "SIZE", "UPDATED", "STATUS")
    widths = [max(len(headers[i]), *(len(row[i]) for row in rows)) for i in range(4)]
    print("  ".join(header.ljust(widths[i]) for i, header in enumerate(headers)).rstrip())
    for row in rows:
        print("  ".join(value.ljust(widths[i]) for i, value in enumerate(row)).rstrip())
    if details:
        for item in items:
            print(f"  {item.path.name}: {item.path}")
            if item.url:
                print(f"    Source: {item.url}")


@contextmanager
def _offline_override(offline: bool):
    if not offline:
        yield
        return
    previous = os.environ.get("DATACACHE_OFFLINE")
    os.environ["DATACACHE_OFFLINE"] = "1"
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop("DATACACHE_OFFLINE", None)
        else:
            os.environ["DATACACHE_OFFLINE"] = previous


def _dir_stats(path: Path) -> tuple[bool, int, int]:
    if not path.exists():
        return False, 0, 0
    total_bytes = 0
    file_count = 0
    for root, _dirs, files in os.walk(path, followlinks=False):
        for name in files:
            try:
                st = os.lstat(os.path.join(root, name))
            except OSError:
                continue
            total_bytes += st.st_size
            file_count += 1
    return True, total_bytes, file_count


def _cli_data_dir(*, create: bool) -> Path:
    if default_data_dir is not _default_data_dir_impl:
        return Path(default_data_dir())
    if create:
        return Path(default_data_dir())
    return resolve_data_dir()


def cmd_data(args):
    import datacache

    action = args.data_command
    read_only = action in ("dirs", "list", "verify", "path")
    dd = _cli_data_dir(create=not read_only)
    json_output = getattr(args, "json", False)

    if action == "dirs":
        loc = resolve_data_location()
        # Honor monkeypatched default_data_dir() in tests if it differs from loc.path
        effective_path = dd
        exists, total_bytes, file_count = _dir_stats(effective_path)
        offline = datacache.offline_source()
        entries = []
        if effective_path.is_dir():
            try:
                children = sorted(os.scandir(effective_path), key=lambda e: e.name)
            except OSError:
                children = []
            for entry in children:
                if not entry.is_dir(follow_symlinks=False):
                    continue
                _, sub_bytes, sub_files = _dir_stats(Path(entry.path))
                entries.append(
                    dict(
                        name=entry.name,
                        path=entry.path,
                        size_bytes=sub_bytes,
                        file_count=sub_files,
                    )
                )
        payload = dict(
            data_dir=str(effective_path),
            origin=loc.origin,
            detail=loc.detail,
            exists=exists,
            size_bytes=total_bytes,
            file_count=file_count,
            offline=offline is not None,
            offline_source=offline,
            subdirectories=entries,
        )
        if json_output:
            print(json.dumps(payload, indent=2))
            return
        print(f"Data directory: {effective_path}")
        print(f"Origin:         {loc.detail} ({loc.origin})")
        if exists:
            print(f"Size:           {datacache.format_bytes(total_bytes)} ({file_count} file(s))")
        else:
            print("Size:           not created yet")
        if offline:
            print(f"Offline:        enabled ({offline})")
        if entries:
            print()
            for entry in entries:
                print(f"  {entry['name']:<20} {datacache.format_bytes(entry['size_bytes']):>10}  ({entry['file_count']} file(s))")
        return

    if action == "verify":
        checks = []
        ok = True
        for item in inventory(dd):
            if not item.exists:
                continue
            status = "ok" if item.size > 0 else "corrupt"
            if status != "ok":
                ok = False
            checks.append(
                dict(
                    kind=item.kind,
                    name=item.key,
                    path=str(item.path),
                    status=status,
                    size_bytes=item.size,
                )
            )
        for version in available_mhc_protein_dataset_versions():
            try:
                paths = validate_mhc_protein_dataset(version, data_dir=dd)
                checks.append(
                    dict(
                        kind="dataset",
                        name=f"{MHC_PROTEIN_DATASET_NAME}:{version}",
                        path=str(paths.records),
                        status="verified",
                    )
                )
            except ProteinDatasetNotInstalledError:
                continue
            except ProteinDatasetError as exc:
                ok = False
                checks.append(
                    dict(
                        kind="dataset",
                        name=f"{MHC_PROTEIN_DATASET_NAME}:{version}",
                        path=str(dd),
                        status="corrupt",
                        error=str(exc),
                    )
                )
        payload = dict(data_dir=str(dd), ok=ok, checks=checks)
        if json_output:
            print(json.dumps(payload, indent=2))
        else:
            if not checks:
                print(f"No downloaded files found under {dd}.")
            else:
                for check in checks:
                    extra = f" ({check['error']})" if check.get("error") else ""
                    print(f"  [{check['status']:<8}] {check['kind']}:{check['name']}{extra}")
        if not ok:
            sys.exit(1)
        return

    if action == "list":
        items = inventory(dd)
        dataset_rows = []
        for version in available_mhc_protein_dataset_versions():
            try:
                paths = validate_mhc_protein_dataset(version, data_dir=dd)
                dataset_rows.append(
                    dict(
                        dataset=MHC_PROTEIN_DATASET_NAME,
                        version=version,
                        status="installed; verified",
                        records=str(paths.records),
                        manifest=str(paths.manifest),
                    )
                )
            except ProteinDatasetNotInstalledError:
                dataset_rows.append(
                    dict(
                        dataset=MHC_PROTEIN_DATASET_NAME,
                        version=version,
                        status="not installed",
                    )
                )
            except ProteinDatasetError:
                dataset_rows.append(
                    dict(
                        dataset=MHC_PROTEIN_DATASET_NAME,
                        version=version,
                        status="invalid; reinstall with --force",
                    )
                )
        if json_output:
            print(
                json.dumps(
                    dict(
                        data_dir=str(dd),
                        files=[
                            dict(
                                key=item.key,
                                kind=item.kind,
                                path=str(item.path),
                                exists=item.exists,
                                size_bytes=item.size,
                                mtime=item.mtime,
                                url=item.url,
                            )
                            for item in items
                        ],
                        datasets=dataset_rows,
                    ),
                    indent=2,
                )
            )
            return
        print(f"Data directory: {dd}")
        _print_file_group("Source FASTA files", [item for item in items if item.kind == "source"], details=args.details)
        _print_file_group("Built CSVs and reports", [item for item in items if item.kind == "built"], details=args.details)
        print("\nVersioned datasets")
        print("DATASET       VERSION             STATUS")
        for row in dataset_rows:
            print(f"{row['dataset']:<13} {row['version']:<19} {row['status']}")
            if args.details and row["status"] == "installed; verified":
                print(f"  Records: {row['records']}")
                print(f"  Manifest: {row['manifest']}")
        if not args.details:
            print("\nUse 'mhcseqs data list --details' for full paths and source URLs.")
        return

    if action in ("install", "fetch"):
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
        if json_output:
            print(
                json.dumps(
                    dict(
                        dataset=MHC_PROTEIN_DATASET_NAME,
                        version=paths.version,
                        records=str(paths.records),
                        manifest=str(paths.manifest),
                    ),
                    indent=2,
                )
            )
        else:
            print(paths.records)
        return

    if action == "refresh":
        from .download import download_fasta

        keys = [args.source] if args.source else list(SOURCES)
        print(f"Refreshing source FASTA(s): {', '.join(keys)}", file=sys.stderr)
        for k in keys:
            download_fasta(k, fasta_dir(dd), force=True)
        print("Done.", file=sys.stderr)
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


def main(argv: list[str] | None = None) -> None:
    """Run the ``mhcseqs`` command line.

    Parameters
    ----------
    argv : list of str, optional
        Arguments after the program name. Defaults to ``sys.argv[1:]``.

    Returns
    -------
    None

    Raises
    ------
    SystemExit
        With status 0 after ``--help``/``--version``, 2 on invalid
        arguments, and 1 when no command (or no ``data`` subcommand) is
        given, ``lookup`` finds no built CSV or no match, or ``data path``
        finds the dataset missing or corrupt. Other failures (for example a
        failed download) propagate as exceptions.
    """
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
        help="Directory for output CSVs (default: resolved mhcseqs data dir)",
    )
    build_parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Directory for downloaded FASTA files (default: <data dir>/fasta)",
    )
    build_parser.add_argument(
        "--force-download",
        action="store_true",
        help="Re-download source FASTA files even when a cached copy exists",
    )
    build_parser.add_argument(
        "--offline",
        action="store_true",
        help="Forbid remote downloads (equivalent to OPENVAX_OFFLINE=1)",
    )

    # lookup
    lookup_parser = subparsers.add_parser("lookup", help="Look up an allele from built CSVs")
    lookup_parser.add_argument("allele", help="Allele name to look up (e.g. HLA-A*02:01)")
    lookup_parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory searched first for built CSVs (before CWD and resolved mhcseqs data dir)",
    )
    lookup_parser.add_argument(
        "--format",
        choices=LOOKUP_FORMATS,
        default="text",
        help="Output format (default: text). All formats emit full sequences.",
    )
    lookup_parser.add_argument(
        "--all",
        action="store_true",
        help="Emit every matching row instead of only the first match per CSV",
    )
    lookup_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Write lookup output to PATH instead of stdout",
    )

    # data — manage cached source FASTAs and built CSVs
    data_parser = subparsers.add_parser("data", help="Manage cached downloads and built CSVs")
    data_sub = data_parser.add_subparsers(dest="data_command")
    list_parser = data_sub.add_parser("list", help="List cached source and built data files")
    list_parser.add_argument("--details", action="store_true", help="Show full paths and source URLs")
    list_parser.add_argument("--json", action="store_true", help="Print machine-readable JSON")
    dirs_parser = data_sub.add_parser("dirs", help="Show where mhcseqs data lives, why it resolved there, and disk usage")
    dirs_parser.add_argument("--json", action="store_true", help="Print machine-readable JSON")
    verify_parser = data_sub.add_parser("verify", help="Verify cached source files and installed datasets offline")
    verify_parser.add_argument("--json", action="store_true", help="Print machine-readable JSON")
    install_parser = data_sub.add_parser("install", aliases=["fetch"], help="Install an immutable versioned dataset")
    install_parser.add_argument("dataset", choices=[MHC_PROTEIN_DATASET_NAME])
    install_parser.add_argument("--version", choices=available_mhc_protein_dataset_versions())
    install_parser.add_argument("--force", action="store_true", help="Re-download and replace this exact cached version")
    install_parser.add_argument(
        "--with-sources",
        action="store_true",
        help="Also install the pinned source and curation artifacts needed for offline regeneration",
    )
    install_parser.add_argument(
        "--offline",
        action="store_true",
        help="Forbid remote downloads (equivalent to OPENVAX_OFFLINE=1)",
    )
    path_parser = data_sub.add_parser("path", help="Print the path to an installed versioned dataset")
    path_parser.add_argument("dataset", choices=[MHC_PROTEIN_DATASET_NAME])
    path_parser.add_argument("--version", choices=available_mhc_protein_dataset_versions())
    path_parser.add_argument("--json", action="store_true", help="Print machine-readable JSON")
    refresh_parser = data_sub.add_parser("refresh", help="Re-download source FASTA(s)")
    refresh_parser.add_argument(
        "--source",
        choices=list(SOURCES),
        default=None,
        help="Only refresh this source (default: all)",
    )
    refresh_parser.add_argument(
        "--offline",
        action="store_true",
        help="Forbid remote downloads (equivalent to OPENVAX_OFFLINE=1)",
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

    args = parser.parse_args(argv)

    import datacache

    with _offline_override(getattr(args, "offline", False)):
        try:
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
        except datacache.OfflineError as exc:
            print(f"error: {exc}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
