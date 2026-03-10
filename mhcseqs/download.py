"""Download IMGT/HLA and IPD-MHC FASTA files."""

from __future__ import annotations

import urllib.request
from pathlib import Path

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


def download_fasta(key: str, dest_dir: Path) -> Path:
    """Download a FASTA source file if not already present. Returns local path."""
    info = SOURCES[key]
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / info["filename"]
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  {info['filename']} already downloaded ({dest.stat().st_size:,} bytes)")
        return dest
    print(f"  Downloading {info['filename']} from {info['url']} ...")
    urllib.request.urlretrieve(info["url"], dest)
    print(f"  Saved {dest.stat().st_size:,} bytes")
    return dest


def download_all(dest_dir: Path) -> dict[str, Path]:
    """Download all FASTA sources. Returns {key: local_path}."""
    paths = {}
    for key in SOURCES:
        paths[key] = download_fasta(key, dest_dir)
    return paths
