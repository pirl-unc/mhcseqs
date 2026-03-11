import csv
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

import mhcseqs
from mhcseqs.__main__ import _default_data_dir, _default_output_dir, main


def test_default_data_dir():
    d = _default_data_dir()
    assert isinstance(d, Path)
    assert d.name == "fasta"


def test_default_output_dir():
    d = _default_output_dir()
    assert isinstance(d, Path)


def test_version_flag(capsys):
    with patch.object(sys, "argv", ["mhcseqs", "--version"]):
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "mhcseqs" in captured.out
    assert "." in captured.out  # version string contains a dot


def test_no_args_prints_help(capsys):
    with patch.object(sys, "argv", ["mhcseqs"]):
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1


def test_lookup_no_csvs(capsys, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with patch.object(sys, "argv", ["mhcseqs", "lookup", "HLA-A*02:01", "--output-dir", str(tmp_path)]):
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert "No built CSVs found" in captured.err


# ---------------------------------------------------------------------------
# lookup() API tests
# ---------------------------------------------------------------------------


def _write_grooves_csv(path, rows):
    """Write a minimal grooves CSV for testing."""
    from mhcseqs.pipeline import GROOVE_FIELDS

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=GROOVE_FIELDS)
        writer.writeheader()
        for row in rows:
            full = {k: "" for k in GROOVE_FIELDS}
            full.update(row)
            writer.writerow(full)


def test_lookup_finds_allele(tmp_path, monkeypatch):
    csv_path = tmp_path / "mhc-binding-grooves.csv"
    _write_grooves_csv(csv_path, [
        {"two_field_allele": "HLA-A*02:01", "groove1": "AAA", "groove2": "BBB"},
    ])
    monkeypatch.chdir(tmp_path)
    rec = mhcseqs.lookup("HLA-A*02:01")
    assert isinstance(rec, mhcseqs.AlleleRecord)
    assert rec.groove1 == "AAA"
    assert rec.groove2 == "BBB"


def test_lookup_case_insensitive(tmp_path, monkeypatch):
    csv_path = tmp_path / "mhc-binding-grooves.csv"
    _write_grooves_csv(csv_path, [
        {"two_field_allele": "HLA-A*02:01", "groove1": "XYZ"},
    ])
    monkeypatch.chdir(tmp_path)
    rec = mhcseqs.lookup("hla-a*02:01")
    assert rec.groove1 == "XYZ"


def test_lookup_not_found(tmp_path, monkeypatch):
    csv_path = tmp_path / "mhc-binding-grooves.csv"
    _write_grooves_csv(csv_path, [
        {"two_field_allele": "HLA-A*02:01"},
    ])
    # Write an empty full-seqs CSV so it doesn't fall through to the repo root
    from mhcseqs.pipeline import FULL_FIELDS

    full_path = tmp_path / "mhc-full-seqs.csv"
    with open(full_path, "w", newline="", encoding="utf-8") as f:
        csv.DictWriter(f, fieldnames=FULL_FIELDS).writeheader()
    monkeypatch.chdir(tmp_path)
    with pytest.raises(KeyError, match="HLA-B"):
        mhcseqs.lookup("HLA-B*07:02")


def test_lookup_no_csv(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(mhcseqs, "_find_csv", lambda name, **kw: (_ for _ in ()).throw(
        FileNotFoundError(f"{name} not found. Run mhcseqs.build() or 'mhcseqs build' first.")
    ))
    with pytest.raises(FileNotFoundError, match="build"):
        mhcseqs.lookup("HLA-A*02:01")
