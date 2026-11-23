from pathlib import Path
from unittest.mock import patch

from mhcseqs.download import SOURCES, download_all, download_fasta


def test_sources_structure():
    assert "imgt_hla" in SOURCES
    assert "ipd_mhc" in SOURCES
    for key, info in SOURCES.items():
        assert "url" in info
        assert "filename" in info
        assert "label" in info
        assert info["url"].startswith("https://")


def test_download_fasta_skips_existing(tmp_path):
    # Create a fake existing file
    fasta_file = tmp_path / "hla_prot.fasta"
    fasta_file.write_text(">test\nACDEFG\n")

    with patch("mhcseqs.download.datacache.fetch_file") as mock_fetch:
        result = download_fasta("imgt_hla", tmp_path)
        mock_fetch.assert_not_called()
        assert result == fasta_file


def test_download_fasta_downloads_new(tmp_path):
    def fake_fetch(url, *, destination, **kwargs):
        Path(destination).parent.mkdir(parents=True, exist_ok=True)
        Path(destination).write_text(">test\nACDEFGH\n")
        return str(destination)

    with patch("mhcseqs.download.datacache.fetch_file", side_effect=fake_fetch):
        result = download_fasta("imgt_hla", tmp_path)
        assert result.exists()
        assert result.stat().st_size > 0
        assert result.name == "hla_prot.fasta"


def test_download_fasta_cleans_up_on_failure(tmp_path):
    def fail_fetch(url, *, destination, **kwargs):
        raise ConnectionError("network failure")

    with patch("mhcseqs.download.datacache.fetch_file", side_effect=fail_fetch):
        try:
            download_fasta("imgt_hla", tmp_path)
        except ConnectionError:
            pass

    # The destination file should NOT exist (atomic download)
    dest = tmp_path / "hla_prot.fasta"
    assert not dest.exists()
    # No temp files should remain
    tmp_files = list(tmp_path.glob("*.tmp"))
    assert len(tmp_files) == 0


def test_download_all_calls_each_source(tmp_path):
    call_log = []

    def fake_fetch(url, *, destination, **kwargs):
        call_log.append(url)
        Path(destination).parent.mkdir(parents=True, exist_ok=True)
        Path(destination).write_text(">test\nACDEFGH\n")
        return str(destination)

    with patch("mhcseqs.download.datacache.fetch_file", side_effect=fake_fetch):
        paths = download_all(tmp_path)

    assert "imgt_hla" in paths
    assert "ipd_mhc" in paths
    assert len(call_log) == 2
