"""On-disk layout of the mhcseqs data directory: where files live and what is there.

The build, the CLI, the versioned-dataset installer and the maintainer
scripts that use the data directory find it through
:func:`resolve_data_dir`: ``$MHCSEQS_DATA`` when set to a non-empty value,
otherwise ``~/.cache/mhcseqs``.  This module also names the built files and
backs the ``mhcseqs data`` CLI commands.  The data directory holds:

* **source** FASTA files downloaded from IMGT/HLA and IPD-MHC, under ``fasta/``
* **built** CSVs and reports produced by ``mhcseqs build``
* versioned datasets and source bundles (managed by
  :mod:`mhcseqs.mhc_protein_dataset`)

These functions are pure (no printing) so they can be unit-tested; the CLI
layer in ``__main__`` handles formatting.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .download import SOURCES

# Environment variable that overrides the default data directory.
DATA_DIR_ENV_VAR = "MHCSEQS_DATA"
SHARED_CACHE_DIR_ENV_VAR = "OPENVAX_DATA_CACHE"

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


def resolve_data_location(data_dir: str | Path | None = None):
    """Resolve the mhcseqs data directory and why it was chosen, without creating it.

    Precedence: explicit *data_dir* -> ``$MHCSEQS_DATA`` ->
    ``$OPENVAX_DATA_CACHE/mhcseqs`` -> populated shared platform cache ->
    populated legacy ``~/.cache/mhcseqs`` -> default shared platform cache.
    """
    import datacache

    return datacache.resolve_cache_dir(
        "mhcseqs",
        DATA_DIR_ENV_VAR,
        override=data_dir,
        legacy=(Path.home() / ".cache" / "mhcseqs",),
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
        Resolved data directory path.
    """
    return Path(resolve_data_location(data_dir).path)


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


def find_csv(
    name: str,
    search_dir: str | Path | None = None,
    *,
    data_dir: str | Path | None = None,
) -> Path:
    """Locate a built CSV file using the unified search order.

    Search order:
    1. *search_dir* (when explicitly provided)
    2. Current working directory (``Path(".")``)
    3. Resolved data directory (*data_dir* if given, else :func:`resolve_data_dir`)
    4. Package parent directory

    Raises
    ------
    FileNotFoundError
        If *name* does not exist in any candidate directory.
    """
    candidates: list[Path] = []
    if search_dir is not None:
        candidates.append(Path(search_dir) / name)
    resolved_data = Path(data_dir) if data_dir is not None else resolve_data_dir()
    candidates.extend(
        [
            Path(".") / name,
            resolved_data / name,
            Path(__file__).resolve().parent.parent / name,
        ]
    )
    seen: set[Path] = set()
    for p in candidates:
        resolved = p.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if p.exists():
            return p
    raise FileNotFoundError(f"{name} not found. Run mhcseqs.build() or 'mhcseqs build' first.")


ALLELE_LOOKUP_FIELDS = (
    "two_field_allele",
    "representative_allele",
    "allele_normalized",
    "allele_raw",
    "source_id",
)


def build_query_candidates(query: str) -> set[str]:
    """Return casefolded candidate strings for matching an allele query against CSV rows."""
    import mhcgnomes

    from .alleles import canonical_allele_name, parse_allele_name
    from .species import full_species_name_alias

    token = query.strip()
    if not token:
        return set()
    candidates = {token.casefold()}
    parsed = parse_allele_name(token)
    if parsed is not None:
        restricted = parsed.restrict_allele_fields(2)
        for obj in (parsed, restricted):
            candidates.add(canonical_allele_name(obj).casefold())
            candidates.add(str(obj.to_string()).casefold())

        # Bridge genus-level prefixes (e.g. FLA -> Felis sp., BoLA -> Bos sp.)
        # and member species (e.g. Felis catus, Bos taurus) using mhcgnomes'
        # species hierarchy so lookup("FLA-E*003:01") finds FelisCatus-E*003:01
        # and lookup("Bota-1*007:01") finds BoLA-1*007:01 (as well as legacy
        # BosSp-1*007:01 rows).
        species = getattr(parsed, "species", None)
        if species is not None and species.name:
            related_prefixes: set[str] = set()
            if species.name.endswith(" sp.") and species.mhc_prefix:
                for sp_obj in mhcgnomes.species.find_matching_species_objects(species.mhc_prefix):
                    if sp_obj.name:
                        alias = full_species_name_alias(sp_obj.name)
                        if alias:
                            related_prefixes.add(alias)
                    if sp_obj.mhc_prefix:
                        related_prefixes.add(sp_obj.mhc_prefix)
            parent = getattr(species, "parent", None)
            if parent is not None and getattr(parent, "name", None) and parent.name.endswith(" sp."):
                if parent.mhc_prefix:
                    related_prefixes.add(parent.mhc_prefix)
                parent_alias = full_species_name_alias(parent.name)
                if parent_alias:
                    related_prefixes.add(parent_alias)
            for rel_prefix in related_prefixes:
                for obj in (parsed, restricted):
                    body = str(obj.to_string(include_species=False))
                    if body:
                        candidates.add(f"{rel_prefix}-{body}".casefold())
    return candidates


def match_allele_row(
    row: dict[str, str],
    query: str,
    *,
    query_candidates: set[str] | None = None,
) -> bool:
    """Return whether a raw or full-seqs CSV row matches *query*."""
    candidates = query_candidates if query_candidates is not None else build_query_candidates(query)
    if not candidates:
        return False
    for field in ALLELE_LOOKUP_FIELDS:
        val = (row.get(field) or "").strip()
        if val and val.casefold() in candidates:
            return True
    return False


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
