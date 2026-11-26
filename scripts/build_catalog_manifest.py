"""Generate gzipped release CSVs and mhc-catalog-manifest.json for Auto Release."""

from __future__ import annotations

import argparse
from pathlib import Path

from mhcseqs.catalog import build_catalog_manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Build compressed release CSVs and catalog manifest.")
    parser.add_argument("output_dir", type=Path, help="Directory containing built mhc-full-seqs.csv and companions.")
    parser.add_argument("--version", default=None, help="Release version (default: mhcseqs.__version__).")
    args = parser.parse_args(argv)
    manifest = build_catalog_manifest(args.output_dir, version=args.version, compress=True)
    print(f"Wrote {manifest}")


if __name__ == "__main__":
    main()
