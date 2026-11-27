"""Prebuilt release catalog installer and fast cached in-memory allele index.

Two capabilities live here:

1. **Prebuilt release catalog installation** (:func:`install_catalog`,
   :func:`validate_catalog`): downloads version-matched ``mhc-full-seqs.csv``
   (and optionally ``mhc-seqs-raw.csv`` and build reports) published on GitHub
   releases into the standard ``mhcseqs`` data directory via ``datacache`` —
   turning fresh setup from a ~15-minute local build into a few-second
   download.
2. **Fast cached in-memory indexing** (:class:`AlleleCatalog`,
   :class:`CatalogEntry`, :func:`load_catalog`): loads ``mhc-full-seqs.csv``
   once per ``(path, mtime, size)`` into compact slotted records with O(1)
   exact, normalized, genus-bridged, and ``source_id`` lookup.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
import threading
import urllib.error
from collections import OrderedDict
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable, Iterator, Mapping, Sequence

from .datafiles import (
    FULL_SEQS_CSV,
    MERGE_REPORT,
    RAW_CSV,
    VALIDATION_REPORT,
    build_query_candidates,
    find_csv,
    resolve_data_dir,
)
from .mhc_protein_dataset import publication_lock
from .version import __version__

if TYPE_CHECKING:
    from .domain_types import AlleleRecord
    from .imgt import GrooveAlignment

CATALOG_DATASET_NAME = "catalog"
CATALOG_MANIFEST_FILENAME = "mhc-catalog-manifest.json"
RELEASE_BASE_URL = "https://github.com/pirl-unc/mhcseqs/releases/download"

_READ_CHUNK_BYTES = 1 << 20
_MAX_CACHED_CATALOGS = 8
_REQUIRED_FULL_CSV_COLUMNS = frozenset({"two_field_allele", "sequence"})
_REQUIRED_RAW_CSV_COLUMNS = frozenset({"allele_normalized", "sequence"})

_CATALOG_CACHE_LOCK = threading.RLock()
_CATALOG_CACHE: OrderedDict[tuple[Path, int, int], AlleleCatalog] = OrderedDict()


class CatalogError(RuntimeError):
    """Failed catalog download, missing artifact, or corrupt cached catalog."""


class CatalogNotInstalledError(CatalogError):
    """The built/installed catalog CSV is not present on disk."""


@dataclass(frozen=True)
class CatalogPaths:
    """Paths to an installed or built ``mhcseqs`` catalog."""

    version: str | None
    sequences: Path
    raw: Path | None = None
    merge_report: Path | None = None
    validation_report: Path | None = None
    manifest: Path | None = None


def _normalize_release_version(version: str | None) -> str:
    raw = (__version__ if version is None else str(version)).strip()
    if raw.startswith(("v", "V")):
        raw = raw[1:].strip()
    if not raw:
        raise CatalogError("Catalog release version must be non-empty.")
    return raw


def _sha256_file(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    total = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_READ_CHUNK_BYTES), b""):
            total += len(chunk)
            digest.update(chunk)
    return total, digest.hexdigest()


def _validate_csv_header(path: Path, required_columns: frozenset[str]) -> None:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader, None)
    except OSError as exc:
        raise CatalogError(f"Cannot read CSV {path}: {exc}") from exc
    if not header:
        raise CatalogError(f"CSV file is empty: {path}")
    missing = required_columns - set(header)
    if missing:
        raise CatalogError(f"CSV {path.name} is missing required column(s): {', '.join(sorted(missing))}")


def _verify_file_metadata(path: Path, metadata: Mapping[str, Any]) -> None:
    if not path.is_file():
        raise CatalogError(f"Missing catalog file: {path}")
    actual_bytes, actual_sha = _sha256_file(path)
    expected_bytes = int(metadata["bytes"])
    if actual_bytes != expected_bytes:
        raise CatalogError(f"Wrong byte count for {path.name}: expected {expected_bytes}, found {actual_bytes}")
    expected_sha = str(metadata["sha256"])
    if actual_sha != expected_sha:
        raise CatalogError(f"SHA-256 mismatch for {path.name}: expected {expected_sha}, found {actual_sha}")


def _fetch_url(url: str, destination: Path) -> None:
    import datacache

    datacache.fetch_file(
        url,
        destination=destination,
        force=True,
        raw=True,
        timeout=300,
        show_progress=datacache.should_show_progress("auto"),
    )


def _is_not_found_error(exc: BaseException) -> bool:
    cur: BaseException | None = exc
    seen: set[int] = set()
    while cur is not None and id(cur) not in seen:
        seen.add(id(cur))
        if isinstance(cur, urllib.error.HTTPError) and cur.code == 404:
            return True
        if "404" in str(cur):
            return True
        cur = cur.__cause__ or cur.__context__
    return False


def _gunzip_file(src_gz: Path, dst_path: Path) -> None:
    with gzip.open(src_gz, "rb") as src, dst_path.open("wb") as dst:
        shutil.copyfileobj(src, dst, length=_READ_CHUNK_BYTES)


def build_catalog_manifest(
    output_dir: str | Path,
    version: str | None = None,
    *,
    compress: bool = True,
) -> Path:
    """Create ``mhc-catalog-manifest.json`` (and optional ``.csv.gz`` archives) for a release build.

    Used by the release workflow to attach integrity metadata and compressed
    archives alongside the uncompressed release CSVs.
    """
    out = Path(output_dir)
    resolved_version = _normalize_release_version(version)
    files_meta: dict[str, dict[str, Any]] = {}

    for filename in (FULL_SEQS_CSV, RAW_CSV, MERGE_REPORT, VALIDATION_REPORT):
        path = out / filename
        if not path.is_file():
            continue
        size, sha = _sha256_file(path)
        entry: dict[str, Any] = {
            "filename": filename,
            "bytes": size,
            "sha256": sha,
        }
        if compress and filename.endswith(".csv"):
            gz_name = f"{filename}.gz"
            gz_path = out / gz_name
            with path.open("rb") as src, gzip.open(gz_path, "wb", compresslevel=9) as dst:
                shutil.copyfileobj(src, dst, length=_READ_CHUNK_BYTES)
            gz_size, gz_sha = _sha256_file(gz_path)
            entry["gzip"] = {
                "filename": gz_name,
                "bytes": gz_size,
                "sha256": gz_sha,
            }
        files_meta[filename] = entry

    if FULL_SEQS_CSV not in files_meta:
        raise CatalogError(f"Cannot build catalog manifest: missing {out / FULL_SEQS_CSV}")

    manifest_payload = {
        "schema_version": 1,
        "dataset": CATALOG_DATASET_NAME,
        "version": resolved_version,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "files": files_meta,
    }
    manifest_path = out / CATALOG_MANIFEST_FILENAME
    manifest_path.write_text(json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_path


def _install_asset(
    base_url: str,
    filename: str,
    dest_path: Path,
    *,
    meta: Mapping[str, Any] | None,
    required_columns: frozenset[str] | None,
    tmp_dir: Path,
) -> dict[str, Any]:
    staged = tmp_dir / filename
    gzip_meta = meta.get("gzip") if isinstance(meta, Mapping) else None
    if isinstance(gzip_meta, Mapping) and gzip_meta.get("filename"):
        gz_name = str(gzip_meta["filename"])
        staged_gz = tmp_dir / gz_name
        _fetch_url(f"{base_url}/{gz_name}", staged_gz)
        _verify_file_metadata(staged_gz, gzip_meta)
        _gunzip_file(staged_gz, staged)
        staged_gz.unlink(missing_ok=True)
    else:
        _fetch_url(f"{base_url}/{filename}", staged)

    if isinstance(meta, Mapping) and "bytes" in meta and "sha256" in meta:
        _verify_file_metadata(staged, meta)
    if required_columns is not None:
        _validate_csv_header(staged, required_columns)

    size, sha = _sha256_file(staged)
    os.replace(staged, dest_path)
    return {"filename": filename, "bytes": size, "sha256": sha}


def install_catalog(
    version: str | None = None,
    *,
    data_dir: str | Path | None = None,
    include_raw: bool = False,
    include_reports: bool = True,
    force: bool = False,
    base_url: str | None = None,
) -> CatalogPaths:
    """Download prebuilt catalog CSVs from a GitHub release into the mhcseqs cache.

    Parameters
    ----------
    version : str, optional
        Package release version to install (e.g. ``"2.15.0"`` or ``"v2.15.0"``).
        Defaults to the installed :data:`mhcseqs.__version__`.
    data_dir : str or Path, optional
        Target directory for ``mhc-full-seqs.csv`` and companion files.
        Defaults to :func:`mhcseqs.datafiles.resolve_data_dir`.
    include_raw : bool, optional
        Also download ``mhc-seqs-raw.csv``. Default: False.
    include_reports : bool, optional
        Also download ``mhc-merge-report.txt`` and
        ``mhc-validation-report.txt`` when available. Default: True.
    force : bool, optional
        Re-download and overwrite existing files even if already present.
    base_url : str, optional
        Override the release asset base URL (primarily for mirrors and tests).

    Returns
    -------
    CatalogPaths
        Paths to the installed catalog files.
    """
    import datacache

    resolved_version = _normalize_release_version(version)
    dest_dir = resolve_data_dir(data_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    release_url = base_url.rstrip("/") if base_url else f"{RELEASE_BASE_URL}/v{resolved_version}"
    lock_target = dest_dir / CATALOG_DATASET_NAME

    with publication_lock(lock_target):
        manifest_path = dest_dir / CATALOG_MANIFEST_FILENAME
        existing_manifest: dict[str, Any] | None = None
        if manifest_path.is_file():
            try:
                existing_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                existing_manifest = None

        full_dest = dest_dir / FULL_SEQS_CSV
        raw_dest = dest_dir / RAW_CSV
        merge_dest = dest_dir / MERGE_REPORT
        val_dest = dest_dir / VALIDATION_REPORT

        if not force and existing_manifest and existing_manifest.get("version") == resolved_version:
            installed_files = existing_manifest.get("files", {})
            full_ok = full_dest.is_file() and FULL_SEQS_CSV in installed_files
            raw_ok = (not include_raw) or (raw_dest.is_file() and RAW_CSV in installed_files)
            if full_ok and raw_ok:
                try:
                    _verify_file_metadata(full_dest, installed_files[FULL_SEQS_CSV])
                    if include_raw:
                        _verify_file_metadata(raw_dest, installed_files[RAW_CSV])
                    return CatalogPaths(
                        version=resolved_version,
                        sequences=full_dest,
                        raw=raw_dest if raw_dest.is_file() else None,
                        merge_report=merge_dest if merge_dest.is_file() else None,
                        validation_report=val_dest if val_dest.is_file() else None,
                        manifest=manifest_path,
                    )
                except CatalogError:
                    pass  # fall through and re-install corrupt files

        tmp_root = Path(tempfile.mkdtemp(prefix=".catalog-install-", dir=dest_dir))
        try:
            remote_manifest: dict[str, Any] | None = None
            staged_manifest = tmp_root / CATALOG_MANIFEST_FILENAME
            try:
                _fetch_url(f"{release_url}/{CATALOG_MANIFEST_FILENAME}", staged_manifest)
                remote_manifest = json.loads(staged_manifest.read_text(encoding="utf-8"))
            except datacache.OfflineError:
                raise
            except Exception as exc:
                if not _is_not_found_error(exc):
                    raise CatalogError(f"Failed to fetch catalog manifest from {release_url}: {exc}") from exc

            remote_files = remote_manifest.get("files", {}) if isinstance(remote_manifest, dict) else {}
            recorded_files: dict[str, dict[str, Any]] = {}
            if existing_manifest and existing_manifest.get("version") == resolved_version:
                for fname, fmeta in existing_manifest.get("files", {}).items():
                    if (dest_dir / fname).is_file():
                        recorded_files[fname] = dict(fmeta)

            targets: list[tuple[str, Path, frozenset[str] | None, bool]] = [
                (FULL_SEQS_CSV, full_dest, _REQUIRED_FULL_CSV_COLUMNS, True),
            ]
            if include_raw:
                targets.append((RAW_CSV, raw_dest, _REQUIRED_RAW_CSV_COLUMNS, True))
            if include_reports:
                targets.append((MERGE_REPORT, merge_dest, None, False))
                targets.append((VALIDATION_REPORT, val_dest, None, False))

            for filename, target_path, req_cols, required in targets:
                if not force and target_path.is_file() and filename in recorded_files:
                    try:
                        _verify_file_metadata(target_path, recorded_files[filename])
                        continue
                    except CatalogError:
                        pass
                try:
                    recorded_files[filename] = _install_asset(
                        release_url,
                        filename,
                        target_path,
                        meta=remote_files.get(filename),
                        required_columns=req_cols,
                        tmp_dir=tmp_root,
                    )
                except datacache.OfflineError:
                    raise
                except Exception as exc:
                    if not required and _is_not_found_error(exc):
                        continue
                    if isinstance(exc, CatalogError):
                        raise
                    raise CatalogError(f"Failed to install {filename} for mhcseqs v{resolved_version} from {release_url}: {exc}") from exc

            local_manifest = {
                "schema_version": 1,
                "dataset": CATALOG_DATASET_NAME,
                "version": resolved_version,
                "source_url": release_url,
                "installed_at": datetime.now(timezone.utc).isoformat(),
                "files": recorded_files,
            }
            staged_local_manifest = tmp_root / "local-manifest.json"
            staged_local_manifest.write_text(json.dumps(local_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            os.replace(staged_local_manifest, manifest_path)
        finally:
            shutil.rmtree(tmp_root, ignore_errors=True)

    return CatalogPaths(
        version=resolved_version,
        sequences=full_dest,
        raw=raw_dest if raw_dest.is_file() else None,
        merge_report=merge_dest if merge_dest.is_file() else None,
        validation_report=val_dest if val_dest.is_file() else None,
        manifest=manifest_path,
    )


def validate_catalog(
    *,
    data_dir: str | Path | None = None,
    search_dir: str | Path | None = None,
) -> CatalogPaths:
    """Validate an installed or locally built catalog offline and return its paths."""
    if search_dir is None and data_dir is not None:
        candidate = resolve_data_dir(data_dir) / FULL_SEQS_CSV
        if not candidate.is_file():
            raise CatalogNotInstalledError("Catalog is not installed. Run 'mhcseqs data install catalog' or 'mhcseqs build' first.")
        full_path = candidate.resolve()
    else:
        try:
            full_path = find_csv(FULL_SEQS_CSV, search_dir=search_dir, data_dir=data_dir).resolve()
        except FileNotFoundError as exc:
            raise CatalogNotInstalledError("Catalog is not installed. Run 'mhcseqs data install catalog' or 'mhcseqs build' first.") from exc

    parent = full_path.parent
    manifest_path = parent / CATALOG_MANIFEST_FILENAME
    version: str | None = None

    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise CatalogError(f"Invalid catalog manifest {manifest_path}: {exc}") from exc
        version = str(manifest.get("version")) if manifest.get("version") else None
        files_meta = manifest.get("files", {})
        if FULL_SEQS_CSV in files_meta:
            _verify_file_metadata(full_path, files_meta[FULL_SEQS_CSV])
        for fname, fmeta in files_meta.items():
            fpath = parent / fname
            if fpath.is_file():
                _verify_file_metadata(fpath, fmeta)
    else:
        _validate_csv_header(full_path, _REQUIRED_FULL_CSV_COLUMNS)

    raw_path = parent / RAW_CSV
    merge_path = parent / MERGE_REPORT
    val_path = parent / VALIDATION_REPORT
    return CatalogPaths(
        version=version,
        sequences=full_path,
        raw=raw_path if raw_path.is_file() else None,
        merge_report=merge_path if merge_path.is_file() else None,
        validation_report=val_path if val_path.is_file() else None,
        manifest=manifest_path if manifest_path.is_file() else None,
    )


# ---------------------------------------------------------------------------
# Fast cached in-memory catalog
# ---------------------------------------------------------------------------


def _parse_int(value: str) -> int:
    if not value:
        return 0
    try:
        return int(value)
    except ValueError:
        return 0


def _parse_bool(value: str) -> bool:
    return value.strip().lower() == "true" if value else False


@dataclass(frozen=True, slots=True)
class CatalogEntry:
    """Compact, typed record for one row of ``mhc-full-seqs.csv``."""

    two_field_allele: str
    representative_allele: str = ""
    protein_seq_selection: str = ""
    gene: str = ""
    mhc_class: str = ""
    chain: str = ""
    species: str = ""
    species_category: str = ""
    species_prefix: str = ""
    source: str = ""
    source_id: str = ""
    seq_len: int = 0
    sequence: str = ""
    mature_start: int = 0
    mature_sequence: str = ""
    groove1: str = ""
    groove2: str = ""
    groove_seq: str = ""
    groove1_len: int = 0
    groove2_len: int = 0
    ig_domain: str = ""
    ig_domain_len: int = 0
    tail: str = ""
    tail_len: int = 0
    domain_architecture: str = ""
    domain_spans: str = ""
    groove_status: str = ""
    groove_flags: tuple[str, ...] = ()
    anchor_type: str = ""
    is_null: bool = False
    is_questionable: bool = False
    is_pseudogene: bool = False
    is_functional: bool = False

    @property
    def allele(self) -> str:
        """Primary canonical allele identifier."""
        return self.two_field_allele or self.representative_allele

    @property
    def is_complete_groove(self) -> bool:
        """Whether this entry has a complete ``status == 'ok'`` binding groove."""
        if self.groove_status != "ok":
            return False
        if self.mhc_class == "I":
            return bool(self.groove1 and self.groove2)
        if self.mhc_class == "II":
            return bool(self.groove_seq)
        return False

    def to_dict(self) -> dict[str, str]:
        """Serialize back to the exact ``mhc-full-seqs.csv`` string-dict schema."""
        return {
            "two_field_allele": self.two_field_allele,
            "representative_allele": self.representative_allele,
            "protein_seq_selection": self.protein_seq_selection,
            "gene": self.gene,
            "mhc_class": self.mhc_class,
            "chain": self.chain,
            "species": self.species,
            "species_category": self.species_category,
            "species_prefix": self.species_prefix,
            "source": self.source,
            "source_id": self.source_id,
            "seq_len": str(self.seq_len),
            "sequence": self.sequence,
            "mature_start": str(self.mature_start),
            "mature_sequence": self.mature_sequence,
            "groove1": self.groove1,
            "groove2": self.groove2,
            "groove_seq": self.groove_seq,
            "groove1_len": str(self.groove1_len),
            "groove2_len": str(self.groove2_len),
            "ig_domain": self.ig_domain,
            "ig_domain_len": str(self.ig_domain_len),
            "tail": self.tail,
            "tail_len": str(self.tail_len),
            "domain_architecture": self.domain_architecture,
            "domain_spans": self.domain_spans,
            "groove_status": self.groove_status,
            "groove_flags": ",".join(self.groove_flags),
            "anchor_type": self.anchor_type,
            "is_null": str(self.is_null),
            "is_questionable": str(self.is_questionable),
            "is_pseudogene": str(self.is_pseudogene),
            "is_functional": str(self.is_functional),
        }

    def as_allele_record(self, *, mutations: str | Sequence[str] = ()) -> AlleleRecord:
        """Run :func:`mhcseqs.decompose_domains` on this entry (with optional mutations)."""
        from .domain_parsing import decompose_domains

        result = decompose_domains(
            self.sequence,
            mhc_class=self.mhc_class,
            chain=self.chain,
            allele=self.two_field_allele,
            gene=self.gene,
            species=self.species,
            mutations=mutations,
        )
        final_seq = result.sequence if result.mutations and result.sequence else self.sequence
        return replace(
            result,
            full_allele=self.representative_allele,
            species=self.species,
            species_category=self.species_category,
            species_prefix=self.species_prefix,
            source=self.source,
            source_id=self.source_id,
            sequence=final_seq,
            seq_len=len(final_seq),
        )

    def align_imgt(self) -> GrooveAlignment:
        """Align this entry's binding groove onto the IMGT G-DOMAIN & pan-class frames."""
        from .imgt import align_groove_to_imgt

        return align_groove_to_imgt(self.groove1, self.groove2, mhc_class=self.mhc_class)


def _entry_from_mapping(row: Mapping[str, Any]) -> CatalogEntry:
    def _s(key: str) -> str:
        val = row.get(key)
        return "" if val is None else str(val).strip()

    seq = _s("sequence")
    flags_raw = _s("groove_flags")
    flags = tuple(tok.strip() for tok in flags_raw.split(",") if tok.strip()) if flags_raw else ()
    seq_len_raw = _s("seq_len")
    g1 = _s("groove1")
    g2 = _s("groove2")
    ig = _s("ig_domain")
    tail = _s("tail")
    return CatalogEntry(
        two_field_allele=_s("two_field_allele"),
        representative_allele=_s("representative_allele"),
        protein_seq_selection=_s("protein_seq_selection"),
        gene=_s("gene"),
        mhc_class=_s("mhc_class"),
        chain=_s("chain"),
        species=_s("species"),
        species_category=_s("species_category"),
        species_prefix=_s("species_prefix"),
        source=_s("source"),
        source_id=_s("source_id"),
        seq_len=_parse_int(seq_len_raw) if seq_len_raw else len(seq),
        sequence=seq,
        mature_start=_parse_int(_s("mature_start")),
        mature_sequence=_s("mature_sequence"),
        groove1=g1,
        groove2=g2,
        groove_seq=_s("groove_seq"),
        groove1_len=_parse_int(_s("groove1_len")) if _s("groove1_len") else len(g1),
        groove2_len=_parse_int(_s("groove2_len")) if _s("groove2_len") else len(g2),
        ig_domain=ig,
        ig_domain_len=_parse_int(_s("ig_domain_len")) if _s("ig_domain_len") else len(ig),
        tail=tail,
        tail_len=_parse_int(_s("tail_len")) if _s("tail_len") else len(tail),
        domain_architecture=_s("domain_architecture"),
        domain_spans=_s("domain_spans"),
        groove_status=_s("groove_status"),
        groove_flags=flags,
        anchor_type=_s("anchor_type"),
        is_null=_parse_bool(_s("is_null")),
        is_questionable=_parse_bool(_s("is_questionable")),
        is_pseudogene=_parse_bool(_s("is_pseudogene")),
        is_functional=_parse_bool(_s("is_functional")),
    )


class AlleleCatalog:
    """In-memory indexed view of ``mhc-full-seqs.csv`` with O(1) allele lookup.

    Supports exact canonical keys, case-insensitive short aliases (``HLA-A*02:01``,
    ``A0201``, ``H2-Kb``), genus/species bridging (``FLA-E*003:01`` <->
    ``FelisCatus-E*003:01``), ``AccessionedAllele`` identifiers
    (``Saha-DAB*uniprot:B3F369``), and bare ``source_id`` database accessions
    (``B3F369``, ``P01898``).
    """

    def __init__(
        self,
        entries: Sequence[CatalogEntry],
        *,
        path: Path | None = None,
        extra_lookup_tokens: Mapping[int, tuple[str, ...]] | None = None,
    ) -> None:
        self._entries: tuple[CatalogEntry, ...] = tuple(entries)
        self.path: Path | None = path

        # Map casefolded lookup keys -> CSV-ordered row indices.
        key_index: dict[str, list[int]] = {}
        primary_index: dict[str, int] = {}

        for idx, entry in enumerate(self._entries):
            seen_for_row: set[str] = set()
            tfa_key = entry.two_field_allele.casefold()
            if tfa_key:
                primary_index.setdefault(tfa_key, idx)
            raw_tokens: list[str] = [
                entry.two_field_allele,
                entry.representative_allele,
                entry.source_id,
            ]
            # Bridge pre-v2.18.0 CSV rows that used legacy gene-less
            # "<Species>*<db>:<acc> [class I|II]" formatting so lookups
            # normalized by mhcgnomes >= 3.92.0 ("<Species>-<db>:<acc>") still
            # hit cached older CSVs in O(1) time.
            tfa = entry.two_field_allele
            if "*" in tfa and ":" in tfa:
                left, right = tfa.split("*", 1)
                if left and "-" not in left and ":" in right:
                    raw_tokens.append(f"{left}-{right.split()[0]}")
            # Index mhcgnomes-canonicalized aliases for curated short-prefix mouse H2 / rat RT1
            # rows (e.g. "H2-Ab1*b" -> "h2-ab*b", "H2-Eb1*k" -> "h2-eb*k") so Class II heterodimer
            # queries like "H2-IAb" and "H2-IEk" resolve in O(1) time.
            if tfa.startswith(("H2-", "RT1-")):
                raw_tokens.extend(build_query_candidates(tfa))
            if extra_lookup_tokens and idx in extra_lookup_tokens:
                raw_tokens.extend(extra_lookup_tokens[idx])
            for token in raw_tokens:
                if not token:
                    continue
                k = token.casefold()
                if k in seen_for_row:
                    continue
                seen_for_row.add(k)
                bucket = key_index.get(k)
                if bucket is None:
                    key_index[k] = [idx]
                else:
                    bucket.append(idx)

        self._primary_index: dict[str, int] = primary_index
        self._key_index: dict[str, tuple[int, ...]] = {k: tuple(v) for k, v in key_index.items()}

    def __len__(self) -> int:
        return len(self._entries)

    def __repr__(self) -> str:
        if self.path is not None:
            return f"AlleleCatalog({len(self._entries):,} entries, path={str(self.path)!r})"
        return f"AlleleCatalog({len(self._entries):,} entries)"

    def __iter__(self) -> Iterator[CatalogEntry]:
        return iter(self._entries)

    def __contains__(self, allele: object) -> bool:
        if not isinstance(allele, str):
            return False
        return self.get(allele) is not None

    @property
    def entries(self) -> tuple[CatalogEntry, ...]:
        """All catalog entries in file order."""
        return self._entries

    @classmethod
    def from_rows(
        cls,
        rows: Iterable[Mapping[str, Any]],
        *,
        path: Path | None = None,
    ) -> AlleleCatalog:
        """Build an :class:`AlleleCatalog` from an iterable of row dicts."""
        entries: list[CatalogEntry] = []
        extra_tokens: dict[int, tuple[str, ...]] = {}
        for idx, row in enumerate(rows):
            entries.append(_entry_from_mapping(row))
            extras = tuple(str(row.get(f) or "").strip() for f in ("allele_normalized", "allele_raw") if str(row.get(f) or "").strip())
            if extras:
                extra_tokens[idx] = extras
        return cls(entries, path=path, extra_lookup_tokens=extra_tokens)

    @classmethod
    def from_csv(cls, path: str | Path) -> AlleleCatalog:
        """Load an :class:`AlleleCatalog` directly from ``mhc-full-seqs.csv``."""
        csv_path = Path(path)
        entries: list[CatalogEntry] = []
        extra_tokens: dict[int, tuple[str, ...]] = {}

        with csv_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader, None)
            if not header:
                return cls((), path=csv_path)

            col = {name: i for i, name in enumerate(header)}
            i_tfa = col.get("two_field_allele", -1)
            i_rep = col.get("representative_allele", -1)
            i_sel = col.get("protein_seq_selection", -1)
            i_gene = col.get("gene", -1)
            i_cls = col.get("mhc_class", -1)
            i_ch = col.get("chain", -1)
            i_sp = col.get("species", -1)
            i_spcat = col.get("species_category", -1)
            i_sppref = col.get("species_prefix", -1)
            i_src = col.get("source", -1)
            i_srcid = col.get("source_id", -1)
            i_slen = col.get("seq_len", -1)
            i_seq = col.get("sequence", -1)
            i_ms = col.get("mature_start", -1)
            i_mseq = col.get("mature_sequence", -1)
            i_g1 = col.get("groove1", -1)
            i_g2 = col.get("groove2", -1)
            i_gseq = col.get("groove_seq", -1)
            i_g1l = col.get("groove1_len", -1)
            i_g2l = col.get("groove2_len", -1)
            i_ig = col.get("ig_domain", -1)
            i_igl = col.get("ig_domain_len", -1)
            i_tail = col.get("tail", -1)
            i_taill = col.get("tail_len", -1)
            i_arch = col.get("domain_architecture", -1)
            i_spans = col.get("domain_spans", -1)
            i_gstat = col.get("groove_status", -1)
            i_gflg = col.get("groove_flags", -1)
            i_anch = col.get("anchor_type", -1)
            i_null = col.get("is_null", -1)
            i_quest = col.get("is_questionable", -1)
            i_pseudo = col.get("is_pseudogene", -1)
            i_func = col.get("is_functional", -1)
            i_anorm = col.get("allele_normalized", -1)
            i_araw = col.get("allele_raw", -1)

            append = entries.append
            for idx, row in enumerate(reader):
                n = len(row)

                def _get(pos: int) -> str:
                    return row[pos].strip() if 0 <= pos < n else ""

                seq = _get(i_seq)
                g1 = _get(i_g1)
                g2 = _get(i_g2)
                ig = _get(i_ig)
                tail = _get(i_tail)
                slen_s = _get(i_slen)
                g1l_s = _get(i_g1l)
                g2l_s = _get(i_g2l)
                igl_s = _get(i_igl)
                taill_s = _get(i_taill)
                flags_s = _get(i_gflg)
                flags = tuple(tok.strip() for tok in flags_s.split(",") if tok.strip()) if flags_s else ()

                append(
                    CatalogEntry(
                        two_field_allele=_get(i_tfa),
                        representative_allele=_get(i_rep),
                        protein_seq_selection=_get(i_sel),
                        gene=_get(i_gene),
                        mhc_class=_get(i_cls),
                        chain=_get(i_ch),
                        species=_get(i_sp),
                        species_category=_get(i_spcat),
                        species_prefix=_get(i_sppref),
                        source=_get(i_src),
                        source_id=_get(i_srcid),
                        seq_len=_parse_int(slen_s) if slen_s else len(seq),
                        sequence=seq,
                        mature_start=_parse_int(_get(i_ms)),
                        mature_sequence=_get(i_mseq),
                        groove1=g1,
                        groove2=g2,
                        groove_seq=_get(i_gseq),
                        groove1_len=_parse_int(g1l_s) if g1l_s else len(g1),
                        groove2_len=_parse_int(g2l_s) if g2l_s else len(g2),
                        ig_domain=ig,
                        ig_domain_len=_parse_int(igl_s) if igl_s else len(ig),
                        tail=tail,
                        tail_len=_parse_int(taill_s) if taill_s else len(tail),
                        domain_architecture=_get(i_arch),
                        domain_spans=_get(i_spans),
                        groove_status=_get(i_gstat),
                        groove_flags=flags,
                        anchor_type=_get(i_anch),
                        is_null=_parse_bool(_get(i_null)),
                        is_questionable=_parse_bool(_get(i_quest)),
                        is_pseudogene=_parse_bool(_get(i_pseudo)),
                        is_functional=_parse_bool(_get(i_func)),
                    )
                )
                if i_anorm >= 0 or i_araw >= 0:
                    extras = tuple(t for t in (_get(i_anorm), _get(i_araw)) if t)
                    if extras:
                        extra_tokens[idx] = extras

        return cls(entries, path=csv_path, extra_lookup_tokens=extra_tokens)

    def _matching_indices(self, query: str) -> list[int]:
        token = (query or "").strip()
        if not token:
            return []
        candidates = build_query_candidates(token)
        if not candidates:
            return []
        matched: set[int] = set()
        key_index = self._key_index
        for cand in candidates:
            hits = key_index.get(cand)
            if hits:
                matched.update(hits)
        return sorted(matched)

    def get(self, allele: str) -> CatalogEntry | None:
        """Return the first matching :class:`CatalogEntry` in O(1) time, or ``None``."""
        token = (allele or "").strip()
        if not token:
            return None
        key = token.casefold()
        # Ultra-fast exact canonical match without invoking mhcgnomes.
        primary_idx = self._primary_index.get(key)
        if primary_idx is not None:
            return self._entries[primary_idx]
        indices = self._matching_indices(token)
        if not indices:
            return None
        return self._entries[indices[0]]

    def get_all(self, allele: str) -> list[CatalogEntry]:
        """Return all matching :class:`CatalogEntry` rows in CSV order."""
        return [self._entries[idx] for idx in self._matching_indices(allele)]

    def lookup(
        self,
        allele: str,
        *,
        mutations: str | Sequence[str] = (),
    ) -> AlleleRecord:
        """Look up *allele* and return a fully decomposed :class:`AlleleRecord`."""
        token = (allele or "").strip()
        entry = self.get(token)
        if entry is None:
            import mhcgnomes

            parsed = mhcgnomes.parse(token, raise_on_error=False)
            if isinstance(parsed, mhcgnomes.Pair):
                raise KeyError(
                    f"Allele {token!r} is a Class II heterodimer pair; "
                    f"use catalog.resolve_heterodimer({token!r}) or catalog.pan_class_features({token!r})."
                )
            loc = f" in {self.path}" if self.path is not None else ""
            raise KeyError(f"Allele {token!r} not found{loc}")
        return entry.as_allele_record(mutations=mutations)

    def lookup_many(
        self,
        alleles: Iterable[str],
        *,
        missing: str = "omit",
    ) -> dict[str, CatalogEntry | None]:
        """Resolve a batch of allele queries into ``{query: CatalogEntry}``.

        Parameters
        ----------
        alleles : iterable of str
            Allele names, aliases, ``AccessionedAllele`` identifiers, or
            database accessions.
        missing : {"omit", "none", "raise"}, optional
            How to handle unresolved queries:
            - ``"omit"`` (default): omit missing alleles from the returned dict.
            - ``"none"``: include missing alleles mapped to ``None``.
            - ``"raise"``: raise :class:`KeyError` listing the missing alleles.
        """
        if missing not in ("omit", "none", "raise"):
            raise ValueError(f"Unsupported missing policy: {missing!r}. Expected 'omit', 'none', or 'raise'.")
        resolved: dict[str, CatalogEntry | None] = {}
        missing_keys: list[str] = []
        for raw in alleles:
            query = str(raw or "").strip()
            if not query:
                continue
            if query in resolved:
                continue
            hit = self.get(query)
            if hit is not None:
                resolved[query] = hit
            elif missing == "none":
                resolved[query] = None
            elif missing == "raise":
                missing_keys.append(query)
        if missing_keys:
            raise KeyError(f"Alleles not found in catalog: {missing_keys[:10]!r} ({len(missing_keys)} missing)")
        return resolved

    def filter(
        self,
        *,
        mhc_class: str | None = None,
        chain: str | None = None,
        species: str | None = None,
        species_category: str | None = None,
        source: str | None = None,
        gene: str | None = None,
        functional_only: bool = False,
        complete_only: bool = False,
    ) -> list[CatalogEntry]:
        """Filter catalog entries by metadata and groove quality flags."""
        mc_norm = mhc_class.strip().upper() if mhc_class else None
        chain_norm = chain.strip().lower() if chain else None
        sp_norm = species.strip().casefold() if species else None
        cat_norm = species_category.strip().casefold() if species_category else None
        src_norm = source.strip().casefold() if source else None
        gene_norm = gene.strip().casefold() if gene else None

        out: list[CatalogEntry] = []
        for entry in self._entries:
            if mc_norm is not None and entry.mhc_class.upper() != mc_norm:
                continue
            if chain_norm is not None and entry.chain.lower() != chain_norm:
                continue
            if sp_norm is not None and entry.species.casefold() != sp_norm:
                continue
            if cat_norm is not None and entry.species_category.casefold() != cat_norm:
                continue
            if src_norm is not None and entry.source.casefold() != src_norm:
                continue
            if gene_norm is not None and entry.gene.casefold() != gene_norm:
                continue
            if functional_only and not entry.is_functional:
                continue
            if complete_only and not entry.is_complete_groove:
                continue
            out.append(entry)
        return out

    def to_dicts(self) -> list[dict[str, str]]:
        """Return all entries as ``mhc-full-seqs.csv`` row dicts."""
        return [entry.to_dict() for entry in self._entries]

    def to_dataframe(self):
        """Return all entries as a :class:`pandas.DataFrame`."""
        import pandas as pd

        return pd.DataFrame(self.to_dicts())

    def _find_default_dra_alpha(self, beta_entry: CatalogEntry) -> CatalogEntry | None:
        """Locate the canonical monomorphic DRA alpha entry for a DRB beta entry's species."""
        prefix = (beta_entry.species_prefix or "").strip()
        if not prefix:
            return None
        for candidate in (f"{prefix}-DRA*01:01", f"{prefix}-DRA*01:02"):
            hit = self.get(candidate)
            if hit is not None and hit.mhc_class == "II" and hit.chain == "alpha" and hit.groove1:
                return hit
        # Fall back to any complete DRA alpha entry for this species prefix
        prefix_cf = prefix.casefold()
        for entry in self._entries:
            if (
                entry.species_prefix.casefold() == prefix_cf
                and entry.mhc_class == "II"
                and entry.chain == "alpha"
                and entry.gene.upper() in ("DRA", "DRA1")
                and entry.groove1
            ):
                return entry
        return None

    def resolve_heterodimer(
        self,
        alpha_or_pair: str,
        beta: str | None = None,
        *,
        allow_default_dra: bool = True,
    ) -> tuple[CatalogEntry, CatalogEntry]:
        """Resolve a Class II alpha/beta heterodimer into ``(alpha_entry, beta_entry)``.

        Accepts:
        - Explicit ``(alpha, beta)`` allele names, e.g. ``("HLA-DQA1*01:02", "HLA-DQB1*06:02")``
        - A slash-paired heterodimer string parsed via ``mhcgnomes``, e.g. ``"HLA-DQA1*01:02/DQB1*06:02"``
        - A single monomorphic-alpha ``DRB`` beta allele (e.g. ``"HLA-DRB1*01:01"``) when
          ``allow_default_dra=True``, pairing it with that species' canonical ``DRA*01:01`` entry.
        """
        import mhcgnomes

        alpha_query: str | None = None
        beta_query: str | None = None

        if beta is not None:
            alpha_query = str(alpha_or_pair or "").strip()
            beta_query = str(beta or "").strip()
        else:
            raw = str(alpha_or_pair or "").strip()
            if not raw:
                raise KeyError("Empty heterodimer query")
            parsed = mhcgnomes.parse(raw, raise_on_error=False)
            if isinstance(parsed, mhcgnomes.Pair) and parsed.alpha is not None and parsed.beta is not None:
                alpha_query = parsed.alpha.to_string()
                beta_query = parsed.beta.to_string()
            elif "/" in raw:
                parts = [p.strip() for p in raw.split("/", 1)]
                if len(parts) == 2 and parts[0] and parts[1]:
                    alpha_query, beta_query = parts[0], parts[1]
            if alpha_query is None:
                # Single-allele query: allow default DRA pairing for DRB beta chains only
                entry = self.get(raw)
                if entry is None:
                    raise KeyError(f"Class II allele {raw!r} not found in catalog")
                if entry.mhc_class != "II" or entry.chain != "beta":
                    raise ValueError(
                        f"Cannot resolve single allele {raw!r} (class={entry.mhc_class!r}, chain={entry.chain!r}) "
                        "as a Class II heterodimer without an explicit partner."
                    )
                if not allow_default_dra or not entry.gene.upper().startswith("DRB"):
                    raise ValueError(f"Class II beta allele {raw!r} (gene={entry.gene!r}) requires an explicit alpha chain.")
                dra = self._find_default_dra_alpha(entry)
                if dra is None:
                    raise KeyError(f"No canonical DRA alpha entry found for {entry.species_prefix!r} to pair with {raw!r}")
                return dra, entry

        alpha_entry = self.get(alpha_query)
        if alpha_entry is None:
            raise KeyError(f"Class II alpha allele {alpha_query!r} not found in catalog")
        beta_entry = self.get(beta_query or "")
        if beta_entry is None:
            raise KeyError(f"Class II beta allele {beta_query!r} not found in catalog")

        if alpha_entry.mhc_class != "II" or alpha_entry.chain != "alpha":
            raise ValueError(f"Expected Class II alpha chain for {alpha_query!r}, got class={alpha_entry.mhc_class!r}, chain={alpha_entry.chain!r}")
        if beta_entry.mhc_class != "II" or beta_entry.chain != "beta":
            raise ValueError(f"Expected Class II beta chain for {beta_query!r}, got class={beta_entry.mhc_class!r}, chain={beta_entry.chain!r}")
        return alpha_entry, beta_entry

    def align_pan_class(
        self,
        allele: str,
        beta: str | None = None,
        *,
        allow_default_dra: bool = True,
    ) -> GrooveAlignment:
        """Align a Class I allele or Class II heterodimer onto the unified 202-slot pan-class IMGT frame."""
        import mhcgnomes

        from .imgt import align_groove_to_imgt

        if beta is not None:
            alpha_entry, beta_entry = self.resolve_heterodimer(allele, beta, allow_default_dra=allow_default_dra)
            return align_groove_to_imgt(alpha_entry.groove1, beta_entry.groove2, mhc_class="II")

        raw = str(allele or "").strip()
        if not raw:
            raise KeyError("Empty allele query")

        parsed = mhcgnomes.parse(raw, raise_on_error=False)
        if isinstance(parsed, mhcgnomes.Pair) or "/" in raw:
            alpha_entry, beta_entry = self.resolve_heterodimer(raw, allow_default_dra=allow_default_dra)
            return align_groove_to_imgt(alpha_entry.groove1, beta_entry.groove2, mhc_class="II")

        entry = self.get(raw)
        if entry is None:
            loc = f" in {self.path}" if self.path is not None else ""
            raise KeyError(f"Allele {raw!r} not found{loc}")

        if entry.mhc_class == "I":
            return align_groove_to_imgt(entry.groove1, entry.groove2, mhc_class="I")
        if entry.mhc_class == "II" and entry.chain == "beta" and allow_default_dra and entry.gene.upper().startswith("DRB"):
            dra = self._find_default_dra_alpha(entry)
            if dra is not None:
                return align_groove_to_imgt(dra.groove1, entry.groove2, mhc_class="II")
        return align_groove_to_imgt(entry.groove1, entry.groove2, mhc_class=entry.mhc_class or "II")

    def pan_class_features(
        self,
        allele: str,
        beta: str | None = None,
        *,
        allow_default_dra: bool = True,
    ) -> dict[str, Any]:
        """Return predictor-ready 202-slot pan-class IMGT & 46-slot pocket features for an allele or Class II pair."""
        return self.align_pan_class(allele, beta, allow_default_dra=allow_default_dra).to_pan_class_features()


def clear_catalog_cache() -> None:
    """Clear the in-memory :func:`load_catalog` cache."""
    with _CATALOG_CACHE_LOCK:
        _CATALOG_CACHE.clear()


def load_catalog(
    path: str | Path | None = None,
    *,
    search_dir: str | Path | None = None,
    data_dir: str | Path | None = None,
    auto_install: bool = False,
) -> AlleleCatalog:
    """Load and cache an :class:`AlleleCatalog` for ``mhc-full-seqs.csv``.

    Results are cached in memory by ``(resolved_path, mtime_ns, size)``, so
    repeated calls are instantaneous and automatically refresh if the CSV on
    disk is rebuilt or replaced.

    Parameters
    ----------
    path : str or Path, optional
        Explicit path to ``mhc-full-seqs.csv``. When omitted, located via
        :func:`mhcseqs.datafiles.find_csv`.
    search_dir : str or Path, optional
        Directory searched first before CWD and the resolved data directory.
    data_dir : str or Path, optional
        Override the default ``mhcseqs`` data directory.
    auto_install : bool, optional
        If True and no local ``mhc-full-seqs.csv`` is found, automatically
        download the prebuilt catalog for the installed package version via
        :func:`install_catalog`. Default: False.
    """
    if path is not None:
        csv_path = Path(path)
    elif search_dir is None and data_dir is not None:
        candidate = resolve_data_dir(data_dir) / FULL_SEQS_CSV
        if candidate.is_file():
            csv_path = candidate
        elif auto_install:
            installed = install_catalog(data_dir=data_dir)
            csv_path = installed.sequences
        else:
            raise FileNotFoundError(
                f"{candidate} not found. Run 'mhcseqs data install catalog' (prebuilt release) or mhcseqs.build() / 'mhcseqs build' first."
            )
    else:
        try:
            csv_path = find_csv(FULL_SEQS_CSV, search_dir=search_dir, data_dir=data_dir)
        except FileNotFoundError:
            if not auto_install:
                raise
            installed = install_catalog(data_dir=data_dir)
            csv_path = installed.sequences

    resolved = csv_path.resolve()
    st = resolved.stat()
    cache_key = (resolved, st.st_mtime_ns, st.st_size)

    with _CATALOG_CACHE_LOCK:
        cached = _CATALOG_CACHE.get(cache_key)
        if cached is not None:
            _CATALOG_CACHE.move_to_end(cache_key)
            return cached

    catalog = AlleleCatalog.from_csv(resolved)
    with _CATALOG_CACHE_LOCK:
        _CATALOG_CACHE[cache_key] = catalog
        _CATALOG_CACHE.move_to_end(cache_key)
        while len(_CATALOG_CACHE) > _MAX_CACHED_CATALOGS:
            _CATALOG_CACHE.popitem(last=False)
    return catalog
