"""Gates keeping runtime model files inside the installed package.

The learned signal-peptide models used to live in a sibling ``data/`` directory
outside ``mhcseqs/``. Editable installs and repo checkouts found them, so tests
and CI passed, but wheels shipped without them and every installed user silently
fell back to no model at all. These tests fail if that regresses.
"""

import sys
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - exercised only on Python 3.10
    import tomli as tomllib

import mhcseqs
from mhcseqs.domain_grammar import SP_BOUNDARY_MODEL_PATH, SP_SEQUENCE_CUE_MODEL_PATH
from mhcseqs.domain_parsing import _load_sp_boundary_model

PACKAGE_DIR = Path(mhcseqs.__file__).resolve().parent
PYPROJECT = PACKAGE_DIR.parent / "pyproject.toml"
RUNTIME_MODEL_PATHS = (SP_BOUNDARY_MODEL_PATH, SP_SEQUENCE_CUE_MODEL_PATH)


def _package_data_entries() -> set[str]:
    with PYPROJECT.open("rb") as handle:
        config = tomllib.load(handle)
    return set(config["tool"]["setuptools"]["package-data"]["mhcseqs"])


def test_runtime_models_live_inside_the_package():
    for path in RUNTIME_MODEL_PATHS:
        assert path.resolve().parent == PACKAGE_DIR, f"{path.name} must sit in mhcseqs/ to ship in the wheel"
        assert path.exists() and path.stat().st_size > 0


def test_runtime_models_are_declared_as_package_data():
    declared = _package_data_entries()
    for path in RUNTIME_MODEL_PATHS:
        assert path.name in declared, f"{path.name} is missing from [tool.setuptools.package-data]"


def test_every_package_data_entry_exists():
    for name in _package_data_entries():
        assert (PACKAGE_DIR / name).exists(), f"package-data lists {name}, which is not in mhcseqs/"


def test_boundary_model_loads_non_empty():
    model = _load_sp_boundary_model()
    assert model.get("groups"), "SP boundary model failed to load; scoring would silently degrade to 0.0"
