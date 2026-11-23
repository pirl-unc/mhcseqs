"""Download IMGT/HLA and IPD-MHC FASTA files."""

from __future__ import annotations

import sys
from pathlib import Path

import datacache

# Official FASTA sources
SOURCES = {
    "imgt_hla": {
        "url": "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/hla_prot.fasta",
        "filename": "hla_prot.fasta",
        "label": "imgt",
        "description": "IMGT/HLA human HLA protein sequences (all alleles)",
    },
    "ipd_mhc": {
        "url": "https://raw.githubusercontent.com/ANHIG/IPDMHC/Latest/MHC_prot.fasta",
        "filename": "ipd_mhc_prot.fasta",
        "label": "ipd_mhc",
        "description": "IPD-MHC non-human MHC protein sequences (all species)",
    },
}


def download_fasta(key: str, dest_dir: Path, *, force: bool = False) -> Path:
    """Download a FASTA source file through datacache. Returns local path.

    A cached copy is reused unless *force* is set, in which case it is always
    re-downloaded (used by ``mhcseqs data refresh`` / ``build --force-download``).
    """
    info = SOURCES[key]
    dest = dest_dir / info["filename"]
    if not force and dest.exists() and dest.stat().st_size > 0:
        print(f"  {info['filename']} already downloaded ({dest.stat().st_size:,} bytes)", file=sys.stderr)
        return dest
    print(f"  Downloading {info['filename']} from {info['url']} ...", file=sys.stderr)
    datacache.fetch_file(
        info["url"],
        destination=dest,
        force=force,
        raw=True,
        timeout=300,
        show_progress=datacache.should_show_progress("auto"),
        record_provenance=True,
    )
    print(f"  Saved {dest.stat().st_size:,} bytes", file=sys.stderr)
    return dest


def download_all(dest_dir: Path, *, force: bool = False) -> dict[str, Path]:
    """Download all FASTA sources. Returns {key: local_path}."""
    paths = {}
    for key in SOURCES:
        paths[key] = download_fasta(key, dest_dir, force=force)
    return paths
