"""On-disk layout of the mhcseqs data directory: where files live and what is there.

This module is the one place that decides where mhcseqs reads and writes
cached data, and it backs the ``mhcseqs data`` CLI commands.  The data
directory is ``$MHCSEQS_DATA`` when set, otherwise ``~/.cache/mhcseqs``
(see :func:`resolve_data_dir`).  It holds:

* **source** FASTA files downloaded from IMGT/HLA and IPD-MHC, under ``fasta/``
* **built** CSVs and reports produced by ``mhcseqs build``
* versioned datasets and source bundles (managed by
  :mod:`mhcseqs.mhc_protein_dataset`)

These functions are pure (no printing) so they can be unit-tested; the CLI
layer in ``__main__`` handles formatting.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from .download import SOURCES

# Environment variable that overrides the default data directory.
DATA_DIR_ENV_VAR = "MHCSEQS_DATA"

# Built artifacts produced by the build pipeline, relative to the output dir.
RAW_CSV = "mhc-seqs-raw.csv"
FULL_SEQS_CSV = "mhc-full-seqs.csv"
MERGE_REPORT = "mhc-merge-report.txt"
VALIDATION_REPORT = "mhc-validation-report.txt"
BUILT_FILES = (
    RAW_CSV,
    FULL_SEQS_CSV,
    MERGE_REPORT,
    VALIDATION_REPORT,
)


def resolve_data_dir(data_dir: str | Path | None = None) -> Path:
    """Return the mhcseqs data directory, without creating it.

    Parameters
    ----------
    data_dir : str or Path, optional
        Explicit directory chosen by the caller. When given, it is returned
        as a :class:`~pathlib.Path` and the environment is ignored.

    Returns
    -------
    Path
        *data_dir* if given; otherwise the value of ``$MHCSEQS_DATA`` when it
        is set to a non-empty string; otherwise ``~/.cache/mhcseqs``.
    """
    if data_dir is not None:
        return Path(data_dir)
    configured = os.environ.get(DATA_DIR_ENV_VAR)
    if configured:
        return Path(configured)
    return Path.home() / ".cache" / "mhcseqs"


@dataclass
class DataFile:
    """One tracked data file (a source FASTA or a built artifact)."""

    key: str  # SOURCES key for sources, filename for built artifacts
    kind: str  # "source" or "built"
    path: Path
    exists: bool
    size: int  # bytes (0 if missing)
    mtime: float  # epoch seconds (0.0 if missing)
    url: str = ""  # upstream URL for sources, "" for built artifacts


def fasta_dir(data_dir: Path) -> Path:
    """Directory holding downloaded source FASTA files."""
    return data_dir / "fasta"


def _stat(path: Path, key: str, kind: str, url: str = "") -> DataFile:
    if path.exists():
        st = path.stat()
        return DataFile(key, kind, path, True, st.st_size, st.st_mtime, url)
    return DataFile(key, kind, path, False, 0, 0.0, url)


def inventory(data_dir: Path) -> list[DataFile]:
    """Return the status of every tracked source and built data file."""
    fdir = fasta_dir(data_dir)
    items = [_stat(fdir / info["filename"], key, "source", info["url"]) for key, info in SOURCES.items()]
    items += [_stat(data_dir / name, name, "built", "") for name in BUILT_FILES]
    return items


def clear(
    data_dir: Path,
    *,
    sources: bool = True,
    built: bool = False,
    source_keys: list[str] | None = None,
) -> list[Path]:
    """Delete cached data files; return the paths actually removed.

    Removes source FASTA files when *sources* is True (restricted to
    *source_keys* if given) and built CSVs/reports when *built* is True.
    """
    removed: list[Path] = []
    for item in inventory(data_dir):
        if item.kind == "source":
            if not sources:
                continue
            if source_keys and item.key not in source_keys:
                continue
        elif item.kind == "built" and not built:
            continue
        if item.exists:
            item.path.unlink()
            removed.append(item.path)
    return removed
