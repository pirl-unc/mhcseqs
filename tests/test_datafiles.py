"""Tests for the cached-data inventory/clear helpers behind `mhcseqs data`."""

from mhcseqs.datafiles import BUILT_FILES, clear, fasta_dir, inventory
from mhcseqs.download import SOURCES


def _populate(data_dir):
    """Create stub source FASTAs and built CSVs under *data_dir*."""
    fdir = fasta_dir(data_dir)
    fdir.mkdir(parents=True, exist_ok=True)
    for info in SOURCES.values():
        (fdir / info["filename"]).write_text(">x\nACDEFG\n")
    for name in BUILT_FILES:
        (data_dir / name).write_text("col\nval\n")


def test_inventory_empty(tmp_path):
    items = inventory(tmp_path)
    assert len(items) == len(SOURCES) + len(BUILT_FILES)
    assert all(not it.exists for it in items)
    assert sum(it.kind == "source" for it in items) == len(SOURCES)
    assert sum(it.kind == "built" for it in items) == len(BUILT_FILES)


def test_inventory_reports_present_files(tmp_path):
    _populate(tmp_path)
    by_kind = {"source": [], "built": []}
    for it in inventory(tmp_path):
        by_kind[it.kind].append(it)
        assert it.exists and it.size > 0 and it.mtime > 0
    # sources carry their upstream URL; built artifacts don't
    assert all(it.url for it in by_kind["source"])
    assert all(it.url == "" for it in by_kind["built"])


def test_clear_sources_only_keeps_built(tmp_path):
    _populate(tmp_path)
    removed = clear(tmp_path)  # default: sources only
    assert len(removed) == len(SOURCES)
    remaining = {it.path.name for it in inventory(tmp_path) if it.exists}
    assert remaining == set(BUILT_FILES)


def test_clear_built_only(tmp_path):
    _populate(tmp_path)
    removed = clear(tmp_path, sources=False, built=True)
    assert len(removed) == len(BUILT_FILES)
    remaining = {it.path.name for it in inventory(tmp_path) if it.exists}
    assert remaining == {info["filename"] for info in SOURCES.values()}


def test_clear_single_source(tmp_path):
    _populate(tmp_path)
    key = next(iter(SOURCES))
    removed = clear(tmp_path, source_keys=[key])
    assert removed == [fasta_dir(tmp_path) / SOURCES[key]["filename"]]
    # the other source FASTA survives
    others = [k for k in SOURCES if k != key]
    survivor = fasta_dir(tmp_path) / SOURCES[others[0]]["filename"]
    assert survivor.exists()


def test_clear_everything(tmp_path):
    _populate(tmp_path)
    removed = clear(tmp_path, sources=True, built=True)
    assert len(removed) == len(SOURCES) + len(BUILT_FILES)
    assert not any(it.exists for it in inventory(tmp_path))


def test_clear_nothing_present(tmp_path):
    assert clear(tmp_path, sources=True, built=True) == []
