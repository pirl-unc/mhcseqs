"""Readable, offline inventory output for cached mhcseqs data."""

from argparse import Namespace

from mhcseqs import __main__ as cli
from mhcseqs.datafiles import fasta_dir
from mhcseqs.download import SOURCES


def test_data_list_groups_files_and_keeps_provenance_in_details(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(cli, "default_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr(cli, "available_mhc_protein_dataset_versions", lambda: ())
    source = next(iter(SOURCES.values()))
    source_dir = fasta_dir(tmp_path)
    source_dir.mkdir()
    (source_dir / source["filename"]).write_text(">example\nACDEFG\n")
    (tmp_path / "mhc-seqs-raw.csv").write_text("allele,sequence\nA,ACDEFG\n")

    cli.cmd_data(Namespace(data_command="list", details=False))
    output = capsys.readouterr().out
    assert "Source FASTA files" in output
    assert "Built CSVs and reports" in output
    assert "Versioned datasets" in output
    assert "FILE" in output and "STATUS" in output
    assert source["filename"] in output
    assert "present" in output and "missing" in output
    assert source["url"] not in output

    cli.cmd_data(Namespace(data_command="list", details=True))
    details = capsys.readouterr().out
    assert str(source_dir / source["filename"]) in details
    assert source["url"] in details
