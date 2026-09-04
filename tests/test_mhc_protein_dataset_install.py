import csv
import gzip
import hashlib
import json
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

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
