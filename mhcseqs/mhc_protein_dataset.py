"""Install and load immutable, independently versioned MHC protein data.

The code package and the data have separate versions. Dataset assets live in
GitHub data releases and are cached under ``$MHCSEQS_DATA/datasets`` (or
``~/.cache/mhcseqs/datasets``). Every downloaded byte count and SHA-256 is
checked before a complete staged directory is published into the cache.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import os
import shutil
import tempfile
import urllib.request
from contextlib import contextmanager
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Iterator

DATASET_NAME = "mhc-proteins"
_REGISTRY_FILE = "mhc_protein_datasets.json"
_READ_SIZE = 1 << 20


class ProteinDatasetError(RuntimeError):
    """Invalid version, failed transfer, or corrupt cached protein data."""


class ProteinDatasetNotInstalledError(ProteinDatasetError):
    """The requested version has no cached installation to validate."""


@dataclass(frozen=True)
class ProteinDatasetPaths:
    """Paths for one installed records artifact and its provenance manifest."""

    version: str
    records: Path
    manifest: Path


@dataclass(frozen=True)
class ProteinSourceBundlePaths:
    """Paths needed to reproduce one records artifact without network access."""

    version: str
    root: Path
    uniprot: Path
    label_curation: Path


def _registry() -> dict[str, object]:
    resource = files("mhcseqs").joinpath(_REGISTRY_FILE)
    with resource.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def available_mhc_protein_dataset_versions() -> tuple[str, ...]:
    """Return all immutable data versions known to this code release."""
    versions = _registry()["versions"]
    return tuple(sorted(versions))


def default_mhc_protein_dataset_version() -> str:
    """Return the data version used when callers do not pin one explicitly."""
    return str(_registry()["default_version"])


def _version_spec(version: str | None) -> tuple[str, dict[str, object]]:
    registry = _registry()
    resolved = version or str(registry["default_version"])
    try:
        spec = registry["versions"][resolved]
    except KeyError:
        available = ", ".join(sorted(registry["versions"]))
        raise ProteinDatasetError(f"Unknown {DATASET_NAME} version {resolved!r}; available: {available}") from None
    return resolved, spec


def _data_root(data_dir: str | Path | None) -> Path:
    if data_dir is not None:
        return Path(data_dir)
    configured = os.environ.get("MHCSEQS_DATA")
    return Path(configured) if configured else Path.home() / ".cache" / "mhcseqs"


def _version_dir(version: str, data_dir: str | Path | None) -> Path:
    return _data_root(data_dir) / "datasets" / DATASET_NAME / version


def mhc_protein_dataset_paths(
    version: str | None = None,
    *,
    data_dir: str | Path | None = None,
) -> ProteinDatasetPaths:
    """Return the expected paths for a data version without downloading it."""
    resolved, spec = _version_spec(version)
    root = _version_dir(resolved, data_dir)
    return ProteinDatasetPaths(
        version=resolved,
        records=root / spec["records"]["filename"],
        manifest=root / spec["records_manifest"]["filename"],
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_READ_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_file(path: Path, metadata: dict[str, object]) -> None:
    try:
        if not path.is_file():
            raise ProteinDatasetError(f"Missing dataset asset: {path}")
        actual_bytes = path.stat().st_size
    except OSError as exc:
        raise ProteinDatasetError(f"Cannot read dataset asset {path}: {exc}") from exc
    expected_bytes = int(metadata["bytes"])
    if actual_bytes != expected_bytes:
        raise ProteinDatasetError(f"Wrong byte count for {path.name}: expected {expected_bytes}, found {actual_bytes}")
    try:
        actual_sha256 = _sha256(path)
    except OSError as exc:
        raise ProteinDatasetError(f"Cannot read dataset asset {path}: {exc}") from exc
    if actual_sha256 != metadata["sha256"]:
        raise ProteinDatasetError(f"SHA-256 mismatch for {path.name}: expected {metadata['sha256']}, found {actual_sha256}")


def _download(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "mhcseqs-data/1.0"})
    with urllib.request.urlopen(request, timeout=300) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=_READ_SIZE)


def _download_asset(base_url: str, metadata: dict[str, object], destination: Path) -> None:
    filename = str(metadata["filename"])
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        _download(f"{base_url}/{filename}", destination)
        _validate_file(destination, metadata)
    except Exception as exc:
        destination.unlink(missing_ok=True)
        if isinstance(exc, ProteinDatasetError):
            raise
        raise ProteinDatasetError(f"Failed to download {filename}: {exc}") from exc


@contextmanager
def _publication_lock(destination: Path):
    """Coordinate readers, recovery, and publication for one data version."""
    lock_path = destination.with_name(f".{destination.name}.lock")
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        handle = lock_path.open("a+b")
    except OSError as exc:
        raise ProteinDatasetError(f"Cannot lock dataset directory {destination}: {exc}") from exc
    with handle:
        if os.name == "nt":
            import msvcrt

            handle.seek(0, os.SEEK_END)
            if handle.tell() == 0:
                handle.write(b"\0")
                handle.flush()
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_LOCK, 1)
            try:
                yield
            finally:
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl

            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _recover_previous_directory(destination: Path) -> None:
    """Restore an interrupted swap before cache checks or downloads, under lock.

    Callers must validate the restored contents before using them. If both
    directories exist, publication completed; keep the backup until a verified
    replacement is ready to be published.
    """
    backup = destination.with_name(f".{destination.name}.previous")
    if backup.exists() and not destination.exists():
        os.replace(backup, destination)


def _publish_staged_directory_locked(staged: Path, destination: Path, *, force: bool) -> bool:
    """Publish a verified directory while the caller holds its version lock."""
    _recover_previous_directory(destination)
    if not force:
        if destination.exists():
            return False
        os.replace(staged, destination)
        return True

    backup = destination.with_name(f".{destination.name}.previous")
    if backup.exists():
        shutil.rmtree(backup)
    if not destination.exists():
        os.replace(staged, destination)
        return True
    os.replace(destination, backup)
    try:
        os.replace(staged, destination)
    except BaseException:
        os.replace(backup, destination)
        raise
    shutil.rmtree(backup)
    return True


def _validate_records(paths: ProteinDatasetPaths, spec: dict[str, object]) -> None:
    if not paths.records.parent.exists():
        raise ProteinDatasetNotInstalledError(f"Dataset is not installed: {paths.version}")
    _validate_file(paths.records, spec["records"])
    _validate_file(paths.manifest, spec["records_manifest"])
    try:
        manifest = json.loads(paths.manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ProteinDatasetError(f"Invalid records manifest {paths.manifest}: {exc}") from exc
    if manifest.get("dataset") != DATASET_NAME or manifest.get("dataset_version") != paths.version:
        raise ProteinDatasetError(f"Records manifest does not describe {DATASET_NAME} {paths.version}")
    records = manifest.get("records", {})
    if records.get("filename") != paths.records.name or records.get("sha256") != spec["records"]["sha256"]:
        raise ProteinDatasetError(f"Records manifest and registry disagree for {paths.version}")


def validate_mhc_protein_dataset(
    version: str | None = None,
    *,
    data_dir: str | Path | None = None,
) -> ProteinDatasetPaths:
    """Validate a cached records/manifest pair and return its paths."""
    resolved, spec = _version_spec(version)
    paths = mhc_protein_dataset_paths(resolved, data_dir=data_dir)
    with _publication_lock(paths.records.parent):
        _recover_previous_directory(paths.records.parent)
        _validate_records(paths, spec)
    return paths


def install_mhc_protein_dataset(
    version: str | None = None,
    *,
    data_dir: str | Path | None = None,
    force: bool = False,
) -> ProteinDatasetPaths:
    """Download, verify, and atomically install one records data version."""
    resolved, spec = _version_spec(version)
    paths = mhc_protein_dataset_paths(resolved, data_dir=data_dir)
    with _publication_lock(paths.records.parent):
        _recover_previous_directory(paths.records.parent)
        if paths.records.parent.exists() and not force:
            _validate_records(paths, spec)
            return paths

    parent = paths.records.parent.parent
    parent.mkdir(parents=True, exist_ok=True)
    staged = Path(tempfile.mkdtemp(dir=parent, prefix=f".{resolved}-"))
    try:
        staged_paths = ProteinDatasetPaths(
            resolved,
            staged / paths.records.name,
            staged / paths.manifest.name,
        )
        base_url = str(spec["release_url"])
        _download_asset(base_url, spec["records"], staged_paths.records)
        _download_asset(base_url, spec["records_manifest"], staged_paths.manifest)
        _validate_records(staged_paths, spec)
        with _publication_lock(paths.records.parent):
            _publish_staged_directory_locked(staged, paths.records.parent, force=force)
            _validate_records(paths, spec)
    finally:
        if staged.exists():
            shutil.rmtree(staged)
    return paths


def install_mhc_protein_source_bundle(
    version: str | None = None,
    *,
    data_dir: str | Path | None = None,
    force: bool = False,
) -> ProteinSourceBundlePaths:
    """Install every source artifact needed for offline data regeneration."""
    resolved, spec = _version_spec(version)
    destination = _data_root(data_dir) / "source-bundles" / DATASET_NAME / resolved
    curation_name = "sp_ground_truth_label_curation.csv"
    result = ProteinSourceBundlePaths(resolved, destination, destination / "uniprot", destination / curation_name)
    assets = spec["source_bundle"]

    def validate(root: Path) -> None:
        for filename, metadata in assets.items():
            relative = Path(filename) if filename == curation_name else Path("uniprot") / filename
            _validate_file(root / relative, metadata)

    with _publication_lock(destination):
        _recover_previous_directory(destination)
        if destination.exists() and not force:
            validate(destination)
            return result

    destination.parent.mkdir(parents=True, exist_ok=True)
    staged = Path(tempfile.mkdtemp(dir=destination.parent, prefix=f".{resolved}-"))
    try:
        for filename, metadata in assets.items():
            relative = Path(filename) if filename == curation_name else Path("uniprot") / filename
            asset_metadata = {"filename": filename, **metadata}
            _download_asset(str(spec["release_url"]), asset_metadata, staged / relative)
        validate(staged)
        with _publication_lock(destination):
            _publish_staged_directory_locked(staged, destination, force=force)
            validate(destination)
    finally:
        if staged.exists():
            shutil.rmtree(staged)
    return result


@contextmanager
def _open_mhc_protein_records(path, *, version, data_dir):
    if path is not None:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            yield handle
        return
    paths = install_mhc_protein_dataset(version, data_dir=data_dir)
    _, spec = _version_spec(paths.version)
    # Capture the compressed bytes (~10 MB for 55k rows) under the same lock as
    # validation. Release before parsing/yielding so slow or nested iterators
    # cannot block other readers or installers, including on Windows.
    with _publication_lock(paths.records.parent):
        _recover_previous_directory(paths.records.parent)
        _validate_records(paths, spec)
        try:
            compressed = paths.records.read_bytes()
        except OSError as exc:
            raise ProteinDatasetError(f"Cannot read dataset asset {paths.records}: {exc}") from exc
    with gzip.open(io.BytesIO(compressed), "rt", encoding="utf-8", newline="") as handle:
        yield handle


def iter_mhc_protein_records(
    path: str | Path | None = None,
    *,
    version: str | None = None,
    data_dir: str | Path | None = None,
) -> Iterator[dict[str, str]]:
    """Yield records, installing if needed and taking a verified cache snapshot.

    Only compressed bytes are buffered; rows are parsed lazily. Explicit paths
    are caller-managed and stream directly without cache validation or locks.
    """
    with _open_mhc_protein_records(path, version=version, data_dir=data_dir) as handle:
        yield from csv.DictReader(handle)


def load_mhc_protein_records(
    path: str | Path | None = None,
    *,
    version: str | None = None,
    data_dir: str | Path | None = None,
) -> list[dict[str, str]]:
    """Load all records as dictionaries, installing the data when absent."""
    return list(iter_mhc_protein_records(path, version=version, data_dir=data_dir))


def load_mhc_protein_dataframe(
    path: str | Path | None = None,
    *,
    version: str | None = None,
    data_dir: str | Path | None = None,
):
    """Load all records as a pandas DataFrame (pandas is an optional dependency)."""
    import pandas as pd

    with _open_mhc_protein_records(path, version=version, data_dir=data_dir) as handle:
        return pd.read_csv(handle, keep_default_na=False)
