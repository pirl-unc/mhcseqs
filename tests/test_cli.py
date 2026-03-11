import sys
from pathlib import Path
from unittest.mock import patch

import pytest

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
