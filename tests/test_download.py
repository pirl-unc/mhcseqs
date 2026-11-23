"""Tests for mhcseqs.download timeout, retry, User-Agent, and atomic writes (#107)."""

from __future__ import annotations

import io
import urllib.error
from pathlib import Path

import pytest

from mhcseqs import download


class _FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
        return False


def test_download_url_sets_timeout_user_agent_and_writes_atomically(tmp_path: Path, monkeypatch):
    seen = {}

    def fake_urlopen(request, timeout=None):
        seen["url"] = request.full_url
        seen["user_agent"] = request.get_header("User-agent")
        seen["timeout"] = timeout
        return _FakeResponse(b">seq1\nAAAA\n")

    monkeypatch.setattr(download.urllib.request, "urlopen", fake_urlopen)
    dest = tmp_path / "out.fasta"
    result = download.download_url("https://example.org/a.fasta", dest, timeout=42)
    assert result == dest
    assert dest.read_bytes() == b">seq1\nAAAA\n"
    assert seen["url"] == "https://example.org/a.fasta"
    assert seen["timeout"] == 42
    assert "pirl-unc/mhcseqs" in seen["user_agent"]
    assert list(tmp_path.glob("*.tmp")) == []


def test_download_url_retries_transient_error_then_succeeds(tmp_path: Path, monkeypatch):
    calls = 0

    def fake_urlopen(request, timeout=None):
        nonlocal calls
        calls += 1
        if calls < 2:
            raise urllib.error.URLError("transient network blip")
        return _FakeResponse(b"ok")

    monkeypatch.setattr(download.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(download.time, "sleep", lambda _s: None)
    dest = tmp_path / "retry.fasta"
    download.download_url("https://example.org/b.fasta", dest, retries=3)
    assert calls == 2
    assert dest.read_bytes() == b"ok"


def test_download_url_fails_fast_on_404_and_leaves_no_partial_file(tmp_path: Path, monkeypatch):
    calls = 0

    def fake_urlopen(request, timeout=None):
        nonlocal calls
        calls += 1
        raise urllib.error.HTTPError(request.full_url, 404, "Not Found", {}, None)

    monkeypatch.setattr(download.urllib.request, "urlopen", fake_urlopen)
    dest = tmp_path / "missing.fasta"
    with pytest.raises(urllib.error.HTTPError):
        download.download_url("https://example.org/404.fasta", dest, retries=3)
    assert calls == 1
    assert not dest.exists()
    assert list(tmp_path.glob("*.tmp")) == []
