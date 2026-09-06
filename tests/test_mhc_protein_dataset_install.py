import csv
import gzip
import hashlib
import json
import threading
from concurrent.futures import ThreadPoolExecutor, TimeoutError
from pathlib import Path
from types import SimpleNamespace

import pytest

from mhcseqs import mhc_protein_dataset as datasets


def _metadata(filename, content):
    return {
        "filename": filename,
        "bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def _registry(records, manifest, sources=None):
    version = "uniprot-test-r1"
    return {
        "dataset": "mhc-proteins",
        "default_version": version,
        "versions": {
            version: {
                "release_url": "https://example.test/release",
                "records": _metadata("records.csv.gz", records),
                "records_manifest": _metadata("records.manifest.json", manifest),
                "source_bundle": sources or {},
            }
        },
    }


def _record_bytes():
    import io

    raw = io.BytesIO()
    with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
        with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
            writer = csv.DictWriter(text, fieldnames=["accession", "sequence"])
            writer.writeheader()
            writer.writerow({"accession": "TEST123", "sequence": "MAAA"})
    return raw.getvalue()


def _asset_pair():
    records = _record_bytes()
    manifest = json.dumps(
        {
            "dataset": "mhc-proteins",
            "dataset_version": "uniprot-test-r1",
            "records": {
                "filename": "records.csv.gz",
                "sha256": hashlib.sha256(records).hexdigest(),
            },
        }
    ).encode()
    return records, manifest


def _fake_download(monkeypatch, files, fail_name=None):
    def download(url, destination):
        name = Path(url).name
        if name == fail_name:
            raise OSError("simulated transfer failure")
        destination.write_bytes(files[name])

    monkeypatch.setattr(datasets, "_download", download)


def test_install_validates_and_publishes_complete_pair(monkeypatch, tmp_path):
    records, manifest = _asset_pair()
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest))
    _fake_download(monkeypatch, {"records.csv.gz": records, "records.manifest.json": manifest})

    paths = datasets.install_mhc_protein_dataset(data_dir=tmp_path)
    assert paths.records.read_bytes() == records
    assert paths.manifest.read_bytes() == manifest
    assert datasets.load_mhc_protein_records(data_dir=tmp_path) == [{"accession": "TEST123", "sequence": "MAAA"}]
    assert not list(paths.records.parent.parent.glob(".uniprot-test-r1-*"))


def test_concurrent_first_installs_validate_and_reuse_winner(monkeypatch, tmp_path):
    records, manifest = _asset_pair()
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest))
    files = {"records.csv.gz": records, "records.manifest.json": manifest}
    manifests_downloaded = threading.Barrier(2)

    def download(url, destination):
        name = Path(url).name
        destination.write_bytes(files[name])
        if name == "records.manifest.json":
            manifests_downloaded.wait(timeout=5)

    monkeypatch.setattr(datasets, "_download", download)
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(datasets.install_mhc_protein_dataset, data_dir=tmp_path) for _ in range(2)]
        installed = [future.result(timeout=10) for future in futures]

    assert installed[0] == installed[1]
    assert installed[0] == datasets.validate_mhc_protein_dataset(data_dir=tmp_path)
    assert not list(installed[0].records.parent.parent.glob(".uniprot-test-r1-*"))


def test_concurrent_forced_installs_serialize_publication(monkeypatch, tmp_path):
    records, manifest = _asset_pair()
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest))
    files = {"records.csv.gz": records, "records.manifest.json": manifest}
    _fake_download(monkeypatch, files)
    paths = datasets.install_mhc_protein_dataset(data_dir=tmp_path)

    manifests_downloaded = threading.Barrier(2)

    def download(url, destination):
        name = Path(url).name
        destination.write_bytes(files[name])
        if name == "records.manifest.json":
            manifests_downloaded.wait(timeout=5)

    monkeypatch.setattr(datasets, "_download", download)
    real_replace = datasets.os.replace
    backup = paths.records.parent.with_name(f".{paths.version}.previous")
    delayed = threading.Event()
    delay_once = threading.Event()

    def replace(source, destination):
        real_replace(source, destination)
        if Path(destination) == backup and not delayed.is_set():
            delayed.set()
            delay_once.wait(timeout=0.1)

    monkeypatch.setattr(datasets.os, "replace", replace)
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(datasets.install_mhc_protein_dataset, data_dir=tmp_path, force=True) for _ in range(2)]
        installed = [future.result(timeout=10) for future in futures]

    assert installed == [paths, paths]
    assert datasets.validate_mhc_protein_dataset(data_dir=tmp_path) == paths
    assert not backup.exists()


def test_force_install_recovers_interrupted_swap(monkeypatch, tmp_path):
    records, manifest = _asset_pair()
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest))
    _fake_download(monkeypatch, {"records.csv.gz": records, "records.manifest.json": manifest})
    paths = datasets.install_mhc_protein_dataset(data_dir=tmp_path)
    backup = paths.records.parent.with_name(f".{paths.version}.previous")
    datasets.os.replace(paths.records.parent, backup)

    assert datasets.install_mhc_protein_dataset(data_dir=tmp_path, force=True) == paths
    assert datasets.validate_mhc_protein_dataset(data_dir=tmp_path) == paths
    assert not backup.exists()


@pytest.fixture(params=["records", "sources"])
def installed_bundle(request, monkeypatch, tmp_path):
    records, manifest = _asset_pair()
    source_files = {"uniprot_artifacts.json": b"{}\n", "sp_ground_truth_label_curation.csv": b"accession,disposition\n"}
    metadata = {name: _metadata(name, content) for name, content in source_files.items()}
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest, metadata))
    files = {"records.csv.gz": records, "records.manifest.json": manifest, **source_files}
    _fake_download(monkeypatch, files)
    installer = datasets.install_mhc_protein_dataset if request.param == "records" else datasets.install_mhc_protein_source_bundle

    def install(**kwargs):
        return installer(data_dir=tmp_path, **kwargs)

    paths = install()
    root = paths.records.parent if request.param == "records" else paths.root
    return SimpleNamespace(kind=request.param, paths=paths, root=root, install=install, data_dir=tmp_path)


def test_normal_install_recovers_interrupted_swap_without_network(monkeypatch, installed_bundle):
    bundle = installed_bundle
    backup = bundle.root.with_name(f".{bundle.root.name}.previous")
    datasets.os.replace(bundle.root, backup)

    def deny_download(*_args, **_kwargs):
        pytest.fail("A complete previous installation must be recovered offline")

    monkeypatch.setattr(datasets, "_download", deny_download)
    assert bundle.install() == bundle.paths
    assert not backup.exists()
    if bundle.kind == "records":
        assert datasets.load_mhc_protein_records(data_dir=bundle.data_dir)[0]["accession"] == "TEST123"


def test_recovered_backup_is_validated_before_reuse(monkeypatch, installed_bundle):
    bundle = installed_bundle
    asset = bundle.paths.records if bundle.kind == "records" else bundle.paths.label_curation
    asset.write_bytes(b"corrupt")
    backup = bundle.root.with_name(f".{bundle.root.name}.previous")
    datasets.os.replace(bundle.root, backup)

    def deny_download(*_args, **_kwargs):
        pytest.fail("A corrupt cache should require an explicit forced reinstall")

    monkeypatch.setattr(datasets, "_download", deny_download)
    with pytest.raises(datasets.ProteinDatasetError, match="byte count"):
        bundle.install()


def test_force_recovers_previous_installation_before_failed_download(monkeypatch, installed_bundle):
    bundle = installed_bundle
    backup = bundle.root.with_name(f".{bundle.root.name}.previous")
    datasets.os.replace(bundle.root, backup)

    def fail_download(*_args, **_kwargs):
        assert bundle.root.is_dir() and not backup.exists()
        raise OSError("offline")

    monkeypatch.setattr(datasets, "_download", fail_download)
    with pytest.raises(datasets.ProteinDatasetError, match="offline"):
        bundle.install(force=True)
    assert bundle.install() == bundle.paths


@pytest.mark.parametrize("reader_kind", ["reuse", "validate", "records", "dataframe"])
def test_readers_wait_for_forced_publication(monkeypatch, installed_bundle, reader_kind):
    bundle = installed_bundle
    if bundle.kind == "sources" and reader_kind != "reuse":
        pytest.skip("Source bundles are validated by their installer")
    if reader_kind == "dataframe":
        pytest.importorskip("pandas")
    readers = {
        "reuse": bundle.install,
        "validate": lambda: datasets.validate_mhc_protein_dataset(data_dir=bundle.data_dir),
        "records": lambda: datasets.load_mhc_protein_records(data_dir=bundle.data_dir),
        "dataframe": lambda: datasets.load_mhc_protein_dataframe(data_dir=bundle.data_dir),
    }
    backup = bundle.root.with_name(f".{bundle.root.name}.previous")
    moved, resume, reading = threading.Event(), threading.Event(), threading.Event()
    real_replace = datasets.os.replace

    def replace(source, destination):
        real_replace(source, destination)
        if Path(destination) == backup:
            moved.set()
            assert resume.wait(timeout=5)

    def read():
        reading.set()
        return readers[reader_kind]()

    monkeypatch.setattr(datasets.os, "replace", replace)
    with ThreadPoolExecutor(max_workers=2) as pool:
        writer = pool.submit(bundle.install, force=True)
        try:
            assert moved.wait(timeout=5)
            reader = pool.submit(read)
            assert reading.wait(timeout=5)
            # Deliberately keep the destination absent until the reader has
            # attempted access; it must wait rather than fail or download.
            with pytest.raises(TimeoutError):
                reader.result(timeout=0.1)
        finally:
            resume.set()
        assert writer.result(timeout=5) == bundle.paths
        result = reader.result(timeout=5)
    if reader_kind == "records":
        assert result == [{"accession": "TEST123", "sequence": "MAAA"}]
    elif reader_kind == "dataframe":
        assert list(result["accession"]) == ["TEST123"]
    else:
        assert result == bundle.paths


def test_iterator_releases_lock_before_yielding_records(monkeypatch, tmp_path):
    records, manifest = _asset_pair()
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest))
    _fake_download(monkeypatch, {"records.csv.gz": records, "records.manifest.json": manifest})
    iterator = datasets.iter_mhc_protein_records(data_dir=tmp_path)
    assert next(iterator)["accession"] == "TEST123"
    installing = threading.Event()

    def force_install():
        installing.set()
        return datasets.install_mhc_protein_dataset(data_dir=tmp_path, force=True)

    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(force_install)
        try:
            assert installing.wait(timeout=5)
            future.result(timeout=5)
            # Both a nested reader and a writer can run while the first
            # iterator is suspended. It retains a stable compressed snapshot.
            assert datasets.load_mhc_protein_records(data_dir=tmp_path)[0]["accession"] == "TEST123"
            assert list(iterator) == []
        finally:
            iterator.close()


def test_validation_recovers_previous_pair_without_download(monkeypatch, tmp_path):
    records, manifest = _asset_pair()
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest))
    _fake_download(monkeypatch, {"records.csv.gz": records, "records.manifest.json": manifest})
    paths = datasets.install_mhc_protein_dataset(data_dir=tmp_path)
    datasets.os.replace(paths.records.parent, paths.records.parent.with_name(f".{paths.version}.previous"))
    assert datasets.validate_mhc_protein_dataset(data_dir=tmp_path) == paths


@pytest.mark.parametrize("suffix", [".csv", ".csv.gz"])
def test_explicit_dataframe_paths_preserve_pandas_compression_detection(monkeypatch, tmp_path, suffix):
    pytest.importorskip("pandas")
    path = tmp_path / f"records{suffix}"
    compressed = _record_bytes()
    path.write_bytes(compressed if suffix.endswith(".gz") else gzip.decompress(compressed))

    def unexpected_install(*_args, **_kwargs):
        pytest.fail("Explicit paths must remain caller-managed")

    monkeypatch.setattr(datasets, "install_mhc_protein_dataset", unexpected_install)
    frame = datasets.load_mhc_protein_dataframe(path)
    assert frame.to_dict("records") == [{"accession": "TEST123", "sequence": "MAAA"}]


def test_failed_force_install_preserves_complete_cached_version(monkeypatch, tmp_path):
    records, manifest = _asset_pair()
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest))
    files = {"records.csv.gz": records, "records.manifest.json": manifest}
    _fake_download(monkeypatch, files)
    paths = datasets.install_mhc_protein_dataset(data_dir=tmp_path)
    before = (paths.records.read_bytes(), paths.manifest.read_bytes())

    _fake_download(monkeypatch, files, fail_name="records.manifest.json")
    with pytest.raises(datasets.ProteinDatasetError, match="simulated transfer failure"):
        datasets.install_mhc_protein_dataset(data_dir=tmp_path, force=True)

    assert (paths.records.read_bytes(), paths.manifest.read_bytes()) == before


def test_checksum_failure_never_publishes_partial_version(monkeypatch, tmp_path):
    records, manifest = _asset_pair()
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest))
    _fake_download(
        monkeypatch,
        {"records.csv.gz": b"CORRUPT", "records.manifest.json": manifest},
    )

    with pytest.raises(datasets.ProteinDatasetError, match="byte count"):
        datasets.install_mhc_protein_dataset(data_dir=tmp_path)

    paths = datasets.mhc_protein_dataset_paths(data_dir=tmp_path)
    assert not paths.records.parent.exists()


def test_unreadable_cache_is_reported_as_dataset_error(monkeypatch, tmp_path):
    records, manifest = _asset_pair()
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest))
    _fake_download(monkeypatch, {"records.csv.gz": records, "records.manifest.json": manifest})
    datasets.install_mhc_protein_dataset(data_dir=tmp_path)

    def fail_read(_path):
        raise OSError("simulated read failure")

    monkeypatch.setattr(datasets, "_sha256", fail_read)
    with pytest.raises(datasets.ProteinDatasetError, match="Cannot read dataset asset"):
        datasets.validate_mhc_protein_dataset(data_dir=tmp_path)


def test_source_bundle_is_installed_in_offline_builder_layout(monkeypatch, tmp_path):
    records, manifest = _asset_pair()
    source_files = {
        "uniprot_artifacts.json": b"{}\n",
        "uniprot_mhc_vertebrata.tsv.gz": b"CURRENT",
        "uniprot_deleted_entries.tsv.gz": b"DELETED",
        "uniprot_taxonomy.tsv.gz": b"TAXONOMY",
        "sp_ground_truth_label_curation.csv": b"accession,disposition\n",
    }
    source_metadata = {
        name: {key: value for key, value in _metadata(name, content).items() if key != "filename"} for name, content in source_files.items()
    }
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest, source_metadata))
    _fake_download(monkeypatch, source_files)

    paths = datasets.install_mhc_protein_source_bundle(data_dir=tmp_path)
    assert (paths.uniprot / "uniprot_artifacts.json").read_bytes() == b"{}\n"
    assert paths.label_curation.read_bytes() == b"accession,disposition\n"


def test_unknown_data_version_lists_available_versions(monkeypatch):
    records, manifest = _asset_pair()
    monkeypatch.setattr(datasets, "_registry", lambda: _registry(records, manifest))
    with pytest.raises(datasets.ProteinDatasetError, match="uniprot-test-r1"):
        datasets.mhc_protein_dataset_paths("missing")


def test_shipped_registry_pins_every_asset():
    registry = datasets._registry()
    assert registry["default_version"] in registry["versions"]
    for spec in registry["versions"].values():
        assets = [spec["records"], spec["records_manifest"], *spec["source_bundle"].values()]
        for metadata in assets:
            assert int(metadata["bytes"]) > 0
            assert len(metadata["sha256"]) == 64
            int(metadata["sha256"], 16)
