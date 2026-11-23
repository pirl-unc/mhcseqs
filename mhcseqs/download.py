"""Download IMGT/HLA, IPD-MHC, and dataset assets with timeouts and retries."""

from __future__ import annotations

import os
import shutil
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

from .version import __version__

DEFAULT_TIMEOUT_SECONDS = 300
DEFAULT_RETRIES = 3
DEFAULT_RETRY_BACKOFF_SECONDS = 1.0
DEFAULT_USER_AGENT = f"mhcseqs/{__version__} (https://github.com/pirl-unc/mhcseqs)"
_CHUNK_SIZE = 1024 * 1024

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


def download_url(
    url: str,
    dest: Path,
    *,
    timeout: float = DEFAULT_TIMEOUT_SECONDS,
    retries: int = DEFAULT_RETRIES,
    backoff_seconds: float = DEFAULT_RETRY_BACKOFF_SECONDS,
    user_agent: str = DEFAULT_USER_AGENT,
) -> Path:
    """Download *url* atomically to *dest* with a socket timeout and transient retries."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    attempts = max(1, int(retries))
    last_error: BaseException | None = None

    for attempt in range(attempts):
        fd, tmp_path = tempfile.mkstemp(dir=dest.parent, suffix=".tmp")
        os.close(fd)
        try:
            request = urllib.request.Request(url, headers={"User-Agent": user_agent})
            with urllib.request.urlopen(request, timeout=timeout) as response, open(tmp_path, "wb") as handle:
                shutil.copyfileobj(response, handle, length=_CHUNK_SIZE)
            os.replace(tmp_path, dest)
            return dest
        except urllib.error.HTTPError as exc:
            _unlink_quietly(tmp_path)
            # Fail fast on permanent client errors (4xx other than 408/429).
            if 400 <= exc.code < 500 and exc.code not in (408, 429):
                raise
            last_error = exc
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            _unlink_quietly(tmp_path)
            last_error = exc
        except BaseException:
            _unlink_quietly(tmp_path)
            raise

        if attempt + 1 < attempts:
            time.sleep(backoff_seconds * (2**attempt))

    assert last_error is not None
    raise last_error


def _unlink_quietly(path: str | Path) -> None:
    try:
        os.unlink(path)
    except OSError:
        pass


def download_fasta(
    key: str,
    dest_dir: Path,
    *,
    force: bool = False,
    timeout: float = DEFAULT_TIMEOUT_SECONDS,
    retries: int = DEFAULT_RETRIES,
) -> Path:
    """Download a FASTA source file. Returns local path.

    A cached copy is reused unless *force* is set, in which case it is always
    re-downloaded (used by ``mhcseqs data refresh`` / ``build --force-download``).
    """
    info = SOURCES[key]
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / info["filename"]
    if not force and dest.exists() and dest.stat().st_size > 0:
        print(f"  {info['filename']} already downloaded ({dest.stat().st_size:,} bytes)")
        return dest
    print(f"  Downloading {info['filename']} from {info['url']} ...")
    download_url(info["url"], dest, timeout=timeout, retries=retries)
    print(f"  Saved {dest.stat().st_size:,} bytes")
    return dest


def download_all(dest_dir: Path, *, force: bool = False) -> dict[str, Path]:
    """Download all FASTA sources. Returns {key: local_path}."""
    paths = {}
    for key in SOURCES:
        paths[key] = download_fasta(key, dest_dir, force=force)
    return paths
