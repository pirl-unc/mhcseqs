import csv
import sys
from types import SimpleNamespace
from unittest.mock import patch

import pytest

import mhcseqs
from mhcseqs import default_data_dir
from mhcseqs.__main__ import main


def test_default_data_dir():
    d = default_data_dir()
    assert isinstance(d, str)
    assert "mhcseqs" in d


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


def test_data_list_reports_versioned_dataset(capsys, tmp_path, monkeypatch):
    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    with patch.object(sys, "argv", ["mhcseqs", "data", "list"]):
        main()
    captured = capsys.readouterr()
    assert "mhc-proteins uniprot-2026_03-r1: not installed" in captured.out


def test_data_install_dispatches_version_and_sources(capsys, tmp_path, monkeypatch):
    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    records = tmp_path / "records.csv.gz"
    manifest = tmp_path / "manifest.json"
    source_root = tmp_path / "source-bundle"
    curation = source_root / "sp_ground_truth_label_curation.csv"
    seen = []

    def install(version, *, data_dir, force):
        seen.append(("records", version, data_dir, force))
        return SimpleNamespace(version=version, records=records, manifest=manifest)

    def install_sources(version, *, data_dir, force):
        seen.append(("sources", version, data_dir, force))
        return SimpleNamespace(root=source_root, label_curation=curation)

    monkeypatch.setattr("mhcseqs.__main__.install_mhc_protein_dataset", install)
    monkeypatch.setattr("mhcseqs.__main__.install_mhc_protein_source_bundle", install_sources)
    with patch.object(
        sys,
        "argv",
        ["mhcseqs", "data", "install", "mhc-proteins", "--version", "uniprot-2026_03-r1", "--with-sources"],
    ):
        main()

    assert seen == [
        ("records", "uniprot-2026_03-r1", tmp_path, False),
        ("sources", "uniprot-2026_03-r1", tmp_path, False),
    ]
    assert str(records) in capsys.readouterr().out


def test_data_path_prints_only_validated_cache(capsys, tmp_path, monkeypatch):
    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    records = tmp_path / "records.csv.gz"
    seen = []

    def validate(version, *, data_dir):
        seen.append((version, data_dir))
        return SimpleNamespace(version="uniprot-2026_03-r1", records=records)

    monkeypatch.setattr("mhcseqs.__main__.validate_mhc_protein_dataset", validate)
    with patch.object(sys, "argv", ["mhcseqs", "data", "path", "mhc-proteins"]):
        main()

    assert seen == [(None, tmp_path)]
    assert capsys.readouterr().out.strip() == str(records)


def test_data_path_rejects_corrupt_cache(capsys, tmp_path, monkeypatch):
    from mhcseqs.mhc_protein_dataset import mhc_protein_dataset_paths

    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    paths = mhc_protein_dataset_paths(data_dir=tmp_path)
    paths.records.parent.mkdir(parents=True)
    paths.records.write_bytes(b"truncated")

    with patch.object(sys, "argv", ["mhcseqs", "data", "path", "mhc-proteins"]):
        with pytest.raises(SystemExit) as exc_info:
            main()

    assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "missing or corrupt" in captured.err
    assert "--force" in captured.err


def test_data_list_marks_corrupt_version_invalid(capsys, tmp_path, monkeypatch):
    from mhcseqs.mhc_protein_dataset import mhc_protein_dataset_paths

    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    paths = mhc_protein_dataset_paths(data_dir=tmp_path)
    paths.records.parent.mkdir(parents=True)
    paths.records.write_bytes(b"truncated")

    with patch.object(sys, "argv", ["mhcseqs", "data", "list"]):
        main()

    assert f"mhc-proteins {paths.version}: invalid; reinstall with --force" in capsys.readouterr().out


def test_lookup_no_csvs(capsys, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    with patch.object(sys, "argv", ["mhcseqs", "lookup", "HLA-A*02:01", "--output-dir", str(tmp_path)]):
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert "No built CSVs found" in captured.err


# ---------------------------------------------------------------------------
# lookup() API tests
# ---------------------------------------------------------------------------

# HLA-A*02:01 full sequence (with signal peptide)
_A0201_SEQ = (
    "MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLG"
    "TLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRA"
    "YLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGD"
    "GTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSLTACKV"
)


def _write_full_seqs_csv(path, rows):
    """Write a minimal full-seqs CSV for testing."""
    from mhcseqs.pipeline import FULL_FIELDS

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FULL_FIELDS)
        writer.writeheader()
        for row in rows:
            full = {k: "" for k in FULL_FIELDS}
            full.update(row)
            writer.writerow(full)


def test_lookup_finds_allele(tmp_path, monkeypatch):
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [
            {
                "two_field_allele": "HLA-A*02:01",
                "representative_allele": "HLA-A*02:01:01:01",
                "gene": "A",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Homo sapiens",
                "species_category": "human",
                "sequence": _A0201_SEQ,
                "seq_len": str(len(_A0201_SEQ)),
            },
        ],
    )
    monkeypatch.chdir(tmp_path)
    rec = mhcseqs.lookup("HLA-A*02:01")
    assert isinstance(rec, mhcseqs.AlleleRecord)
    assert rec.ok
    assert rec.allele == "HLA-A*02:01"
    assert rec.full_allele == "HLA-A*02:01:01:01"
    assert rec.species_category == "human"
    assert len(rec.groove1) == 90
    assert len(rec.groove2) == 93
    assert rec.anchor_cys1 is not None  # populated by extract_groove


def test_lookup_case_insensitive(tmp_path, monkeypatch):
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [
            {
                "two_field_allele": "HLA-A*02:01",
                "mhc_class": "I",
                "chain": "alpha",
                "sequence": _A0201_SEQ,
            },
        ],
    )
    monkeypatch.chdir(tmp_path)
    rec = mhcseqs.lookup("hla-a*02:01")
    assert rec.ok
    assert len(rec.groove1) == 90


def test_lookup_with_mutations(tmp_path, monkeypatch):
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [
            {
                "two_field_allele": "HLA-A*02:01",
                "mhc_class": "I",
                "chain": "alpha",
                "sequence": _A0201_SEQ,
            },
        ],
    )
    monkeypatch.chdir(tmp_path)
    rec = mhcseqs.lookup("HLA-A*02:01", mutations=["K66A"])
    assert rec.groove1[65] == "A"
    assert rec.mutations == ("K66A",)


def test_lookup_not_found(tmp_path, monkeypatch):
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [
            {"two_field_allele": "HLA-A*02:01", "mhc_class": "I", "sequence": _A0201_SEQ},
        ],
    )
    monkeypatch.chdir(tmp_path)
    with pytest.raises(KeyError, match="HLA-B"):
        mhcseqs.lookup("HLA-B*07:02")


def test_lookup_no_csv(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        mhcseqs,
        "_find_csv",
        lambda name, **kw: (_ for _ in ()).throw(FileNotFoundError(f"{name} not found. Run mhcseqs.build() or 'mhcseqs build' first.")),
    )
    with pytest.raises(FileNotFoundError, match="build"):
        mhcseqs.lookup("HLA-A*02:01")
