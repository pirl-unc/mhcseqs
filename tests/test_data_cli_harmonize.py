"""Tests for mhcseqs OPENVAX_DATA_CACHE, data dirs/verify/list --json, and --offline."""

import json
import os

import pytest

from mhcseqs.__main__ import main


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    for key in (
        "MHCSEQS_DATA",
        "OPENVAX_DATA_CACHE",
        "DATACACHE_OFFLINE",
        "OPENVAX_OFFLINE",
    ):
        monkeypatch.delenv(key, raising=False)


def test_data_dirs_and_list_json_readonly(tmp_path, monkeypatch, capsys):
    shared = tmp_path / "shared"
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(shared))
    monkeypatch.setenv("OPENVAX_OFFLINE", "1")

    main(["data", "dirs", "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert payload["data_dir"] == str(shared / "mhcseqs")
    assert payload["origin"] == "shared_env"
    assert payload["detail"] == "OPENVAX_DATA_CACHE"
    assert payload["exists"] is False
    assert payload["offline"] is True
    assert not shared.exists()

    main(["data", "list", "--json"])
    listed = json.loads(capsys.readouterr().out)
    assert listed["data_dir"] == str(shared / "mhcseqs")
    assert len(listed["files"]) >= 2
    assert not shared.exists()


def test_data_verify_detects_empty_source(tmp_path, monkeypatch, capsys):
    data_dir = tmp_path / "mhcseqs"
    monkeypatch.setenv("MHCSEQS_DATA", str(data_dir))
    fasta = data_dir / "fasta" / "hla_prot.fasta"
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">HLA-A*02:01\nMAVMAP\n")

    main(["data", "verify", "--json"])
    report = json.loads(capsys.readouterr().out)
    assert report["ok"] is True
    assert report["checks"][0]["status"] == "ok"

    fasta.write_text("")
    with pytest.raises(SystemExit) as exc:
        main(["data", "verify", "--json"])
    assert exc.value.code == 1
    bad = json.loads(capsys.readouterr().out)
    assert bad["ok"] is False
    assert bad["checks"][0]["status"] == "corrupt"


def test_data_refresh_offline_exits_cleanly_and_restores_env(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("MHCSEQS_DATA", str(tmp_path / "mhcseqs"))
    with pytest.raises(SystemExit) as exc:
        main(["data", "refresh", "--source", "imgt_hla", "--offline"])
    assert exc.value.code == 1
    assert "Offline mode is enabled" in capsys.readouterr().err
    assert "DATACACHE_OFFLINE" not in os.environ
