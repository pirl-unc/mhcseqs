"""Gates keeping runtime data files inside the installed package.

The learned signal-peptide models used to live in a sibling ``data/`` directory
outside ``mhcseqs/``. Editable installs and repo checkouts found them, so tests
and CI passed, but wheels shipped without them and every installed user silently
fell back to no model at all.

Every loader for these files fails soft -- ``{}`` for the models, ``[]`` for the
GenBank SahaI references -- so a file missing from the wheel produces degraded
output rather than an error. The gate therefore has to run in both directions:
every declared entry must exist, and every shipped data file must be declared.
"""

import csv
import json
import sys
from pathlib import Path

import pytest

import mhcseqs
from mhcseqs.domain_grammar import SP_BOUNDARY_MODEL_PATH, SP_SEQUENCE_CUE_MODEL_PATH
from mhcseqs.domain_parsing import _load_sp_boundary_model, _load_sp_sequence_cue_model
from mhcseqs.pipeline import _GENBANK_SAHA_CSV, _load_genbank_saha_references

PACKAGE_DIR = Path(mhcseqs.__file__).resolve().parent
PYPROJECT = PACKAGE_DIR.parent / "pyproject.toml"

# Files mhcseqs opens at runtime, each behind a fail-soft loader.
RUNTIME_DATA_PATHS = (SP_BOUNDARY_MODEL_PATH, SP_SEQUENCE_CUE_MODEL_PATH, _GENBANK_SAHA_CSV)


def _package_data_entries() -> set[str]:
    """Read the declared package-data filenames, or skip if unreadable here."""
    if sys.version_info >= (3, 11):
        import tomllib
    else:  # pragma: no cover - exercised only on Python 3.10
        tomllib = pytest.importorskip("tomli", reason="reading pyproject.toml needs tomllib (3.11+) or tomli")
    if not PYPROJECT.exists():
        pytest.skip("pyproject.toml is not available next to an installed package")
    with PYPROJECT.open("rb") as handle:
        config = tomllib.load(handle)
    return set(config["tool"]["setuptools"]["package-data"]["mhcseqs"])


def _shipped_data_files() -> set[str]:
    """Every non-Python file directly inside the package directory."""
    return {p.name for p in PACKAGE_DIR.iterdir() if p.is_file() and p.suffix not in {".py", ".pyc", ".typed"}}


def test_runtime_data_files_exist_and_are_non_empty():
    for path in RUNTIME_DATA_PATHS:
        assert path.exists() and path.stat().st_size > 0, f"{path.name} is missing or empty"


def test_every_package_data_entry_exists():
    for name in _package_data_entries():
        assert (PACKAGE_DIR / name).exists(), f"package-data lists {name}, which is not in mhcseqs/"


def test_every_shipped_data_file_is_declared():
    """The direction that matters: an undeclared file is silently absent from wheels."""
    undeclared = _shipped_data_files() - _package_data_entries()
    assert not undeclared, f"data files in mhcseqs/ missing from [tool.setuptools.package-data]: {sorted(undeclared)}"


def test_runtime_data_files_are_declared():
    declared = _package_data_entries()
    for path in RUNTIME_DATA_PATHS:
        assert path.name in declared, f"{path.name} is read at runtime but not declared as package-data"


def test_fail_soft_loaders_actually_loaded_something():
    """Each loader swallows a missing file, so assert on loaded content, not existence."""
    assert _load_sp_boundary_model().get("groups"), "SP boundary model did not load; scoring would degrade to 0.0"
    assert _load_sp_sequence_cue_model(), "SP sequence-cue model did not load; cue scores would degrade to 0.0"
    assert _load_genbank_saha_references(), "GenBank SahaI references did not load; devil alleles would be dropped"


def test_shipped_data_files_are_parseable():
    """A truncated or corrupt file passes an existence check but breaks fail-soft loaders."""
    for name in sorted(_shipped_data_files()):
        path = PACKAGE_DIR / name
        if path.suffix == ".json":
            with path.open("r", encoding="utf-8") as handle:
                assert json.load(handle), f"{name} parsed as empty JSON"
        elif path.suffix == ".csv":
            with path.open("r", encoding="utf-8") as handle:
                assert next(csv.reader(handle), None), f"{name} has no header row"
