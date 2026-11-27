"""Tests for prebuilt release catalog installation and fast cached AlleleCatalog."""

from __future__ import annotations

import csv
import json
import os
import urllib.error
from pathlib import Path

import pytest

import mhcseqs
from mhcseqs.__main__ import main
from mhcseqs.catalog import (
    CATALOG_MANIFEST_FILENAME,
    AlleleCatalog,
    CatalogEntry,
    CatalogError,
    CatalogNotInstalledError,
    build_catalog_manifest,
    clear_catalog_cache,
    install_catalog,
    load_catalog,
    validate_catalog,
)
from mhcseqs.pipeline import FULL_FIELDS, RAW_FIELDS

_A0201_SEQ = (
    "MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLG"
    "TLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRA"
    "YLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGD"
    "GTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSLTACKV"
)


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            full = {k: "" for k in fieldnames}
            full.update(row)
            writer.writerow(full)
    return path


def _sample_full_rows() -> list[dict[str, str]]:
    return [
        {
            "two_field_allele": "HomoSapiens-A*02:01",
            "representative_allele": "HomoSapiens-A*02:01:01:01",
            "gene": "A",
            "mhc_class": "I",
            "chain": "alpha",
            "species": "Homo sapiens",
            "species_category": "human",
            "species_prefix": "HomoSapiens",
            "source": "imgt",
            "source_id": "HLA00005",
            "seq_len": str(len(_A0201_SEQ)),
            "sequence": _A0201_SEQ,
            "mature_start": "24",
            "mature_sequence": _A0201_SEQ[24:],
            "groove1": _A0201_SEQ[24:114],
            "groove2": _A0201_SEQ[114:207],
            "groove_seq": _A0201_SEQ[24:207],
            "groove1_len": "90",
            "groove2_len": "93",
            "groove_status": "ok",
            "is_functional": "True",
        },
        {
            "two_field_allele": "FelisCatus-E*003:01",
            "representative_allele": "FelisCatus-E*003:01",
            "gene": "E",
            "mhc_class": "I",
            "chain": "alpha",
            "species": "Felis catus",
            "species_category": "carnivore",
            "species_prefix": "FelisCatus",
            "source": "uniprot",
            "source_id": "Q95482",
            "seq_len": str(len(_A0201_SEQ)),
            "sequence": _A0201_SEQ,
            "groove1": _A0201_SEQ[24:114],
            "groove2": _A0201_SEQ[114:207],
            "groove_seq": _A0201_SEQ[24:207],
            "groove_status": "ok",
            "is_functional": "True",
        },
        {
            "two_field_allele": "SarcophilusHarrisii-uniprot:G3V1L5",
            "representative_allele": "SarcophilusHarrisii-uniprot:G3V1L5",
            "gene": "~ref:SarcophilusHarrisii|G3V1L5:MHC-I",
            "mhc_class": "I",
            "chain": "alpha",
            "species": "Sarcophilus harrisii",
            "species_category": "other_mammal",
            "species_prefix": "SarcophilusHarrisii",
            "source": "uniprot",
            "source_id": "G3V1L5",
            "seq_len": str(len(_A0201_SEQ)),
            "sequence": _A0201_SEQ,
            "groove1": _A0201_SEQ[24:114],
            "groove2": _A0201_SEQ[114:207],
            "groove_seq": _A0201_SEQ[24:207],
            "groove_status": "ok",
            "is_functional": "True",
        },
        {
            "two_field_allele": "HomoSapiens-DRB1*01:01",
            "representative_allele": "HomoSapiens-DRB1*01:01:01",
            "gene": "DRB1",
            "mhc_class": "II",
            "chain": "beta",
            "species": "Homo sapiens",
            "species_category": "human",
            "species_prefix": "HomoSapiens",
            "source": "imgt",
            "source_id": "HLA00664",
            "seq_len": "266",
            "sequence": "M" * 266,
            "groove2": "A" * 94,
            "groove_seq": "A" * 94,
            "groove_status": "beta1_only_fallback",
            "is_functional": "True",
        },
    ]


def test_allele_catalog_fast_lookup_bridging_and_batch(tmp_path):
    csv_path = _write_csv(tmp_path / "mhc-full-seqs.csv", FULL_FIELDS, _sample_full_rows())
    catalog = AlleleCatalog.from_csv(csv_path)
    assert len(catalog) == 4
    assert "HLA-A*02:01" in catalog
    assert "FLA-E*003:01" in catalog
    assert "G3V1L5" in catalog
    assert "HLA-B*99:99" not in catalog

    # Exact, short alias, 4-field, genus-bridged, accessioned, and bare source_id lookups
    assert catalog.get("HomoSapiens-A*02:01").source_id == "HLA00005"
    assert catalog.get("HLA-A*02:01").source_id == "HLA00005"
    assert catalog.get("A0201").source_id == "HLA00005"
    assert catalog.get("HLA-A*02:01:01:01").source_id == "HLA00005"
    assert catalog.get("FLA-E*003:01").two_field_allele == "FelisCatus-E*003:01"
    assert catalog.get("Q95482").two_field_allele == "FelisCatus-E*003:01"
    assert catalog.get("Saha-uniprot:G3V1L5").source_id == "G3V1L5"
    assert catalog.get("Saha-classI*uniprot:G3V1L5").source_id == "G3V1L5"
    assert catalog.get("Saha*uniprot:G3V1L5 class I").source_id == "G3V1L5"
    assert catalog.get("G3V1L5").two_field_allele == "SarcophilusHarrisii-uniprot:G3V1L5"

    # mhcgnomes >= 3.95.0 class-qualified rows ("-classI*uniprot:...") bridge to class-agnostic queries
    v395_rows = [dict(r) for r in _sample_full_rows()]
    v395_rows[2]["two_field_allele"] = "SarcophilusHarrisii-classI*uniprot:G3V1L5"
    v395_rows[2]["representative_allele"] = "SarcophilusHarrisii-classI*uniprot:G3V1L5"
    v395_catalog = AlleleCatalog.from_rows(v395_rows)
    assert v395_catalog.get("Saha-classI*uniprot:G3V1L5").source_id == "G3V1L5"
    assert v395_catalog.get("Saha-uniprot:G3V1L5").source_id == "G3V1L5"
    assert v395_catalog.get("Saha*uniprot:G3V1L5 class I").source_id == "G3V1L5"

    # Legacy pre-v2.18.0 CSV rows with "*uniprot:... class I" still resolve seamlessly
    legacy_rows = [dict(r) for r in _sample_full_rows()]
    legacy_rows[2]["two_field_allele"] = "SarcophilusHarrisii*uniprot:G3V1L5 class I"
    legacy_rows[2]["representative_allele"] = "SarcophilusHarrisii*uniprot:G3V1L5 class I"
    legacy_catalog = AlleleCatalog.from_rows(legacy_rows)
    assert legacy_catalog.get("Saha-uniprot:G3V1L5").source_id == "G3V1L5"
    assert legacy_catalog.get("Saha-classI*uniprot:G3V1L5").source_id == "G3V1L5"
    assert legacy_catalog.get("Saha*uniprot:G3V1L5 class I").source_id == "G3V1L5"

    # Full AlleleRecord decomposition with mutations
    rec = catalog.lookup("HLA-A*02:01", mutations=["K66A"])
    assert rec.ok
    assert rec.groove1[65] == "A"
    assert rec.source_id == "HLA00005"
    with pytest.raises(KeyError):
        catalog.lookup("HLA-B*99:99")

    # Batch lookup_many policies
    batch_omit = catalog.lookup_many(["HLA-A*02:01", "FLA-E*003:01", "G3V1L5", "NotAnAllele"])
    assert set(batch_omit) == {"HLA-A*02:01", "FLA-E*003:01", "G3V1L5"}
    batch_none = catalog.lookup_many(["HLA-A*02:01", "NotAnAllele"], missing="none")
    assert batch_none["NotAnAllele"] is None
    with pytest.raises(KeyError, match="NotAnAllele"):
        catalog.lookup_many(["HLA-A*02:01", "NotAnAllele"], missing="raise")

    # Filtering & conversions
    c1_complete = catalog.filter(mhc_class="I", functional_only=True, complete_only=True)
    assert len(c1_complete) == 3
    c2_complete = catalog.filter(mhc_class="II", complete_only=True)
    assert len(c2_complete) == 0  # DRB1 is beta1_only_fallback
    c2_func = catalog.filter(mhc_class="II", functional_only=True)
    assert len(c2_func) == 1
    assert isinstance(c2_func[0], CatalogEntry)
    assert len(catalog.to_dicts()) == 4
    df = catalog.to_dataframe()
    assert len(df) == 4
    assert set(df["two_field_allele"]) == {r["two_field_allele"] for r in _sample_full_rows()}


def test_load_catalog_caches_by_mtime_and_invalidates_on_update(tmp_path):
    clear_catalog_cache()
    csv_path = _write_csv(tmp_path / "mhc-full-seqs.csv", FULL_FIELDS, _sample_full_rows()[:2])
    cat1 = load_catalog(csv_path)
    cat2 = load_catalog(csv_path)
    assert cat1 is cat2
    assert len(cat1) == 2

    # Rewrite with more rows and bump mtime so cache refreshes automatically
    _write_csv(csv_path, FULL_FIELDS, _sample_full_rows())
    st = csv_path.stat()
    os.utime(csv_path, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
    cat3 = load_catalog(csv_path)
    assert cat3 is not cat1
    assert len(cat3) == 4


def test_install_catalog_with_manifest_gzip_and_checksum_verification(tmp_path, monkeypatch):
    # Stage a fake GitHub release directory with compressed CSVs + manifest
    release_dir = tmp_path / "release"
    _write_csv(release_dir / "mhc-full-seqs.csv", FULL_FIELDS, _sample_full_rows())
    _write_csv(
        release_dir / "mhc-seqs-raw.csv",
        RAW_FIELDS,
        [{"allele_normalized": "HomoSapiens-A*02:01", "two_field_allele": "HomoSapiens-A*02:01", "sequence": _A0201_SEQ}],
    )
    (release_dir / "mhc-merge-report.txt").write_text("Merge report\n", encoding="utf-8")
    (release_dir / "mhc-validation-report.txt").write_text("Validation report\n", encoding="utf-8")
    manifest_src = build_catalog_manifest(release_dir, version="2.15.0", compress=True)
    assert manifest_src.is_file()
    assert (release_dir / "mhc-full-seqs.csv.gz").is_file()
    assert (release_dir / "mhc-seqs-raw.csv.gz").is_file()

    fetched_urls: list[str] = []

    def fake_fetch(url: str, destination: Path) -> None:
        fetched_urls.append(url)
        filename = url.rsplit("/", 1)[-1]
        src = release_dir / filename
        if not src.is_file():
            raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)
        destination.write_bytes(src.read_bytes())

    monkeypatch.setattr("mhcseqs.catalog._fetch_url", fake_fetch)

    cache_dir = tmp_path / "cache"
    paths = install_catalog("2.15.0", data_dir=cache_dir, include_raw=True)
    assert paths.version == "2.15.0"
    assert paths.sequences == cache_dir / "mhc-full-seqs.csv"
    assert paths.raw == cache_dir / "mhc-seqs-raw.csv"
    assert paths.manifest == cache_dir / CATALOG_MANIFEST_FILENAME
    # Used gzipped CSVs, not raw uncompressed CSVs!
    assert any(u.endswith("/mhc-full-seqs.csv.gz") for u in fetched_urls)
    assert any(u.endswith("/mhc-seqs-raw.csv.gz") for u in fetched_urls)
    assert not any(u.endswith("/mhc-full-seqs.csv") for u in fetched_urls)

    # Validate offline & load catalog
    validated = validate_catalog(data_dir=cache_dir)
    assert validated.version == "2.15.0"
    cat = load_catalog(data_dir=cache_dir)
    assert len(cat) == 4

    # Idempotent re-install makes zero network calls when already verified
    fetched_urls.clear()
    install_catalog("2.15.0", data_dir=cache_dir, include_raw=True)
    assert fetched_urls == []

    # Corrupt file is caught by validate_catalog()
    paths.sequences.write_text("corrupted", encoding="utf-8")
    with pytest.raises(CatalogError, match="Wrong byte count|SHA-256"):
        validate_catalog(data_dir=cache_dir)


def test_install_catalog_falls_back_when_release_has_no_manifest(tmp_path, monkeypatch):
    # Legacy release (like v2.14.0) with only uncompressed CSVs and no manifest
    release_dir = tmp_path / "legacy-release"
    _write_csv(release_dir / "mhc-full-seqs.csv", FULL_FIELDS, _sample_full_rows())
    (release_dir / "mhc-validation-report.txt").write_text("ok\n", encoding="utf-8")

    def fake_fetch(url: str, destination: Path) -> None:
        filename = url.rsplit("/", 1)[-1]
        src = release_dir / filename
        if not src.is_file():
            raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)
        destination.write_bytes(src.read_bytes())

    monkeypatch.setattr("mhcseqs.catalog._fetch_url", fake_fetch)

    cache_dir = tmp_path / "cache"
    paths = install_catalog("v2.14.0", data_dir=cache_dir)
    assert paths.version == "2.14.0"
    assert paths.sequences.is_file()
    assert paths.raw is None
    assert paths.merge_report is None
    assert paths.validation_report == cache_dir / "mhc-validation-report.txt"
    assert paths.manifest is not None and paths.manifest.is_file()
    assert load_catalog(data_dir=cache_dir).get("HLA-A*02:01") is not None


def test_load_catalog_auto_install_and_validate_missing(tmp_path, monkeypatch):
    cache_dir = tmp_path / "empty-cache"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("mhcseqs.datafiles.__file__", str(tmp_path / "pkg" / "datafiles.py"))
    with pytest.raises(CatalogNotInstalledError):
        validate_catalog(data_dir=cache_dir)

    installed_csv = _write_csv(tmp_path / "installed" / "mhc-full-seqs.csv", FULL_FIELDS, _sample_full_rows())
    called = []

    def fake_install(version=None, *, data_dir=None, **kwargs):
        called.append(data_dir)
        return mhcseqs.CatalogPaths(version="2.15.0", sequences=installed_csv)

    monkeypatch.setattr("mhcseqs.catalog.install_catalog", fake_install)
    cat = load_catalog(data_dir=cache_dir, auto_install=True)
    assert called == [cache_dir]
    assert len(cat) == 4


def test_cli_data_install_and_path_catalog(capsys, tmp_path, monkeypatch):
    monkeypatch.setenv("MHCSEQS_DATA", str(tmp_path))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("mhcseqs.datafiles.__file__", str(tmp_path / "pkg" / "datafiles.py"))

    # Missing catalog -> data path catalog exits 1
    with pytest.raises(SystemExit) as exc_info:
        main(["data", "path", "catalog"])
    assert exc_info.value.code == 1
    assert "mhcseqs data install catalog" in capsys.readouterr().err

    # Mock install_catalog via CLI dispatch
    seq_csv = _write_csv(tmp_path / "mhc-full-seqs.csv", FULL_FIELDS, _sample_full_rows())
    seen = []

    def fake_install(version, *, data_dir, include_raw, force):
        seen.append((version, data_dir, include_raw, force))
        return mhcseqs.CatalogPaths(version=version or "2.15.0", sequences=seq_csv)

    monkeypatch.setattr("mhcseqs.__main__.install_catalog", fake_install)
    main(["data", "install", "catalog", "--version", "2.15.0", "--include-raw", "--force"])
    assert seen == [("2.15.0", tmp_path, True, True)]
    assert str(seq_csv) in capsys.readouterr().out

    # Now data path catalog succeeds (plain and JSON)
    main(["data", "path", "catalog"])
    assert capsys.readouterr().out.strip() == str(seq_csv)

    main(["data", "path", "catalog", "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert payload["dataset"] == "catalog"
    assert payload["sequences"] == str(seq_csv)


def test_catalog_resolve_pair_and_pan_class_features():
    """Issue #167 / #173: AlleleCatalog resolves Class II pairs (incl. DRA default) and exports pan-class features."""
    from mhcseqs.imgt import (
        HLA_A0201_G_ALPHA1,
        HLA_A0201_G_ALPHA2,
        HLA_DRA0101_G_ALPHA,
        HLA_DRB1_0101_G_BETA,
        IMGT_PAN_CLASS_GROOVE_SLOTS,
        IMGT_PAN_CLASS_POCKET_SLOTS,
    )

    rows = [
        {
            "two_field_allele": "HomoSapiens-A*02:01",
            "representative_allele": "HomoSapiens-A*02:01:01:01",
            "gene": "A",
            "mhc_class": "I",
            "chain": "alpha",
            "species": "Homo sapiens",
            "species_category": "human",
            "species_prefix": "HomoSapiens",
            "source": "imgt",
            "source_id": "HLA00005",
            "sequence": _A0201_SEQ,
            "groove1": HLA_A0201_G_ALPHA1,
            "groove2": HLA_A0201_G_ALPHA2,
            "groove_seq": HLA_A0201_G_ALPHA1 + HLA_A0201_G_ALPHA2,
            "groove_status": "ok",
            "is_functional": "True",
        },
        {
            "two_field_allele": "HomoSapiens-DRA*01:01",
            "representative_allele": "HomoSapiens-DRA*01:01:01:01",
            "gene": "DRA",
            "mhc_class": "II",
            "chain": "alpha",
            "species": "Homo sapiens",
            "species_category": "human",
            "species_prefix": "HomoSapiens",
            "source": "imgt",
            "source_id": "HLA00662",
            "sequence": HLA_DRA0101_G_ALPHA,
            "groove1": HLA_DRA0101_G_ALPHA,
            "groove2": "",
            "groove_seq": HLA_DRA0101_G_ALPHA,
            "groove_status": "ok",
            "is_functional": "True",
        },
        {
            "two_field_allele": "HomoSapiens-DRB1*01:01",
            "representative_allele": "HomoSapiens-DRB1*01:01:01",
            "gene": "DRB1",
            "mhc_class": "II",
            "chain": "beta",
            "species": "Homo sapiens",
            "species_category": "human",
            "species_prefix": "HomoSapiens",
            "source": "imgt",
            "source_id": "HLA00664",
            "sequence": HLA_DRB1_0101_G_BETA,
            "groove1": "",
            "groove2": HLA_DRB1_0101_G_BETA,
            "groove_seq": HLA_DRB1_0101_G_BETA,
            "groove_status": "ok",
            "is_functional": "True",
        },
        {
            "two_field_allele": "HomoSapiens-DQB1*06:02",
            "representative_allele": "HomoSapiens-DQB1*06:02:01",
            "gene": "DQB1",
            "mhc_class": "II",
            "chain": "beta",
            "species": "Homo sapiens",
            "species_category": "human",
            "species_prefix": "HomoSapiens",
            "source": "imgt",
            "source_id": "HLA00646",
            "sequence": HLA_DRB1_0101_G_BETA,
            "groove1": "",
            "groove2": HLA_DRB1_0101_G_BETA,
            "groove_seq": HLA_DRB1_0101_G_BETA,
            "groove_status": "ok",
            "is_functional": "True",
        },
    ]
    catalog = AlleleCatalog.from_rows(rows)

    # Explicit (alpha, beta), slash-pair string, and default DRA pairing for DRB1
    a1, b1 = catalog.resolve_pair("HLA-DRA*01:01", "HLA-DRB1*01:01")
    a2, b2 = catalog.resolve_pair("HLA-DRA*01:01/DRB1*01:01")
    a3, b3 = catalog.resolve_pair("HLA-DRB1*01:01")
    assert (a1.two_field_allele, b1.two_field_allele) == ("HomoSapiens-DRA*01:01", "HomoSapiens-DRB1*01:01")
    assert (a2.two_field_allele, b2.two_field_allele) == (a1.two_field_allele, b1.two_field_allele)
    assert (a3.two_field_allele, b3.two_field_allele) == (a1.two_field_allele, b1.two_field_allele)

    # Non-DRB polymorphic beta chain must NOT auto-pair with DRA
    with pytest.raises(ValueError, match="explicit alpha chain"):
        catalog.resolve_pair("HLA-DQB1*06:02")
    with pytest.raises(ValueError, match="explicit alpha chain"):
        catalog.resolve_pair("HLA-DRB1*01:01", allow_default_dra=False)
    # Wrong chain order raises ValueError
    with pytest.raises(ValueError, match="Expected Class II alpha"):
        catalog.resolve_pair("HLA-DRB1*01:01", "HLA-DRA*01:01")

    # Unified pan-class features for Class I and Class II
    f1 = catalog.pan_class_features("HLA-A*02:01")
    f2 = catalog.pan_class_features("HLA-DRB1*01:01")
    f2_pair = catalog.pan_class_features("HLA-DRA*01:01/DRB1*01:01")
    assert f1["mhc_class"] == "I"
    assert f2["mhc_class"] == "II"
    assert f2["seq"] == f2_pair["seq"]
    assert f2["pocket_seq"] == f2_pair["pocket_seq"]
    assert len(f1["seq"]) == len(IMGT_PAN_CLASS_GROOVE_SLOTS) == 202
    assert len(f2["seq"]) == len(IMGT_PAN_CLASS_GROOVE_SLOTS) == 202
    assert len(f1["pocket_seq"]) == len(IMGT_PAN_CLASS_POCKET_SLOTS) == 46
    assert len(f2["pocket_seq"]) == len(IMGT_PAN_CLASS_POCKET_SLOTS) == 46

    # CatalogEntry.align_imgt() convenience method
    entry_a = catalog.get("HLA-A*02:01")
    assert entry_a is not None
    assert entry_a.align_imgt().pan_class_fixed_seq == f1["seq"]


def test_catalog_ux_h2_heterodimer_bf2_padding_repr_and_pair_hint():
    """Issue #171: H2-IAb/IEk, 2<->3-digit BF2*21:01 padding, repr, and pair hint."""
    from mhcseqs.imgt import HLA_A0201_G_ALPHA1, HLA_A0201_G_ALPHA2, HLA_DRA0101_G_ALPHA, HLA_DRB1_0101_G_BETA

    rows = [
        {
            "two_field_allele": "H2-Aa*b",
            "representative_allele": "H2-Aa*b",
            "gene": "Aa",
            "mhc_class": "II",
            "chain": "alpha",
            "species": "Mus musculus",
            "species_category": "murine",
            "species_prefix": "H2",
            "source": "uniprot",
            "source_id": "P04228",
            "sequence": HLA_DRA0101_G_ALPHA,
            "groove1": HLA_DRA0101_G_ALPHA,
            "groove2": "",
            "groove_seq": HLA_DRA0101_G_ALPHA,
            "groove_status": "ok",
            "is_functional": "True",
        },
        {
            "two_field_allele": "H2-Ab1*b",
            "representative_allele": "H2-Ab1*b",
            "gene": "Ab1",
            "mhc_class": "II",
            "chain": "beta",
            "species": "Mus musculus",
            "species_category": "murine",
            "species_prefix": "H2",
            "source": "uniprot",
            "source_id": "P14483",
            "sequence": HLA_DRB1_0101_G_BETA,
            "groove1": "",
            "groove2": HLA_DRB1_0101_G_BETA,
            "groove_seq": HLA_DRB1_0101_G_BETA,
            "groove_status": "ok",
            "is_functional": "True",
        },
        {
            "two_field_allele": "GallusGallus-BF2*021:01",
            "representative_allele": "GallusGallus-BF2*021:01",
            "gene": "BF2",
            "mhc_class": "I",
            "chain": "alpha",
            "species": "Gallus gallus",
            "species_category": "bird",
            "species_prefix": "GallusGallus",
            "source": "ipd_mhc",
            "source_id": "GAGA001",
            "sequence": _A0201_SEQ,
            "groove1": HLA_A0201_G_ALPHA1,
            "groove2": HLA_A0201_G_ALPHA2,
            "groove_seq": HLA_A0201_G_ALPHA1 + HLA_A0201_G_ALPHA2,
            "groove_status": "ok",
            "is_functional": "True",
        },
    ]
    cat = AlleleCatalog.from_rows(rows)
    assert repr(cat) == "AlleleCatalog(num_entries=3)"
    assert str(cat) == repr(cat)

    # Mouse H2-Ab*b / H2-IAb / H2-Aa*b/Ab1*b resolve seamlessly
    assert cat.get("H2-Ab*b") is not None
    assert cat.get("H2-Ab*b").two_field_allele == "H2-Ab1*b"
    a_ent, b_ent = cat.resolve_pair("H2-IAb")
    assert (a_ent.two_field_allele, b_ent.two_field_allele) == ("H2-Aa*b", "H2-Ab1*b")
    a_ent2, b_ent2 = cat.resolve_pair("H2-Aa*b/Ab1*b")
    assert (a_ent2.two_field_allele, b_ent2.two_field_allele) == ("H2-Aa*b", "H2-Ab1*b")

    # Helpful KeyError hint when a Class II Pair is passed to single-allele lookup()
    with pytest.raises(KeyError, match="resolve_pair"):
        cat.lookup("H2-IAb")

    # 2-digit <-> 3-digit zero-padded & colon-less chicken BF2 queries hit GallusGallus-BF2*021:01 (mhcgnomes >= 3.96.0)
    for q in ("BF2*21:01", "BF2*2101", "Gaga-BF2*21:01", "Gaga-BF2*2101", "GallusGallus-BF2*21:01", "Gaga-BF2*021:01"):
        hit = cat.get(q)
        assert hit is not None, f"Expected hit for {q}"
        assert hit.two_field_allele == "GallusGallus-BF2*021:01"


def test_catalog_inline_mhcgnomes_mutations_for_alleles_and_pairs():
    """Issue #175: get, lookup, resolve_pair, and pan_class_features apply mhcgnomes inline mutations."""
    from mhcseqs.domain_parsing import decompose_domains
    from mhcseqs.imgt import HLA_A0201_G_ALPHA1, HLA_A0201_G_ALPHA2, HLA_DRA0101_G_ALPHA, HLA_DRB1_0101_G_BETA

    # Full-length DRB1*01:01 sequence (SP + mature GDTR... beta1 + beta2 + TM, UniProt P04229)
    drb1_full = (
        "MVCLKLPGGSCMTALTVTLMVLSSPLALA"
        + "GDTR"
        + HLA_DRB1_0101_G_BETA
        + "RVQPKVTVYPSKTQPLQHHNLLVCSVSGFYPGSIEVRWFRNGQEEKAGVVSTGLIQNGDWTFQTLVMLETVPRSGEVYTCQVEHPSVTSPLTVEW"
        + "RARSESAQSKMLSGVGGFVLGLLFLGAGLFIYFRNQKGHSGLQPTGFLS"
    )
    dra_full = (
        "MAISGVPVLGFFIIAVLMSAQESWA"
        + "IKEEH"
        + HLA_DRA0101_G_ALPHA
        + "SPVTVLTNSPVELREPNVLICFIDKFTPPVVNVTWLRNGKPVTTGVSETVFLPREDHLFRKFHYLPFLPSTEDVYDCRVEHWGLDEPLLKHWEF"
        + "DAPSPLPETTENVVCALGLTVGLVGIIIGTIFIIKGLRKSNAAERRGPL"
    )
    dra_rec = decompose_domains(dra_full, mhc_class="II", chain="alpha", allele="HomoSapiens-DRA*01:01")
    drb1_rec = decompose_domains(drb1_full, mhc_class="II", chain="beta", allele="HomoSapiens-DRB1*01:01")
    rows = [
        {
            "two_field_allele": "HomoSapiens-A*02:01",
            "representative_allele": "HomoSapiens-A*02:01:01:01",
            "gene": "A",
            "mhc_class": "I",
            "chain": "alpha",
            "species": "Homo sapiens",
            "species_category": "human",
            "species_prefix": "HomoSapiens",
            "source": "imgt",
            "source_id": "HLA00005",
            "sequence": _A0201_SEQ,
            "groove1": HLA_A0201_G_ALPHA1,
            "groove2": HLA_A0201_G_ALPHA2,
            "groove_seq": HLA_A0201_G_ALPHA1 + HLA_A0201_G_ALPHA2,
            "groove_status": "ok",
            "is_functional": "True",
        },
        {
            "two_field_allele": "HomoSapiens-DRA*01:01",
            "representative_allele": "HomoSapiens-DRA*01:01:01:01",
            "gene": "DRA",
            "mhc_class": "II",
            "chain": "alpha",
            "species": "Homo sapiens",
            "species_category": "human",
            "species_prefix": "HomoSapiens",
            "source": "imgt",
            "source_id": "HLA00662",
            "sequence": dra_full,
            "mature_start": str(dra_rec.mature_start),
            "groove1": dra_rec.groove1,
            "groove2": "",
            "groove_seq": dra_rec.groove_seq,
            "groove_status": "ok",
            "is_functional": "True",
        },
        {
            "two_field_allele": "HomoSapiens-DRB1*01:01",
            "representative_allele": "HomoSapiens-DRB1*01:01:01:01",
            "gene": "DRB1",
            "mhc_class": "II",
            "chain": "beta",
            "species": "Homo sapiens",
            "species_category": "human",
            "species_prefix": "HomoSapiens",
            "source": "imgt",
            "source_id": "HLA00664",
            "sequence": drb1_full,
            "mature_start": str(drb1_rec.mature_start),
            "groove1": "",
            "groove2": drb1_rec.groove2,
            "groove_seq": drb1_rec.groove_seq,
            "groove_status": "ok",
            "is_functional": "True",
        },
    ]
    cat = AlleleCatalog.from_rows(rows)

    # Class I inline mutant via get, lookup, and pan_class_features
    entry_mut = cat.get("HLA-A*02:01 K66A mutant")
    assert entry_mut is not None and entry_mut.groove1[65] == "A"
    rec_mut = cat.lookup("HLA-A*02:01 K66A mutant")
    assert rec_mut.groove1[65] == "A"
    f_wt_i = cat.pan_class_features("HLA-A*02:01")
    f_mut_i = cat.pan_class_features("HLA-A*02:01 K66A mutant")
    assert f_wt_i["seq"] != f_mut_i["seq"]

    # Class II Pair inline mutant (885 quantitative BA rows in presto!)
    a_mut, b_mut = cat.resolve_pair("HLA-DRA*01:01/DRB1*01:01 C30S mutant")
    assert a_mut.groove1 == dra_rec.groove1
    assert drb1_rec.groove2[29] == "C"
    assert b_mut.groove2[29] == "S"

    # Single DRB1 inline mutant via get and resolve_pair (auto-pairs with DRA*01:01)
    drb1_mut = cat.get("HLA-DRB1*01:01 C30S mutant")
    assert drb1_mut is not None and drb1_mut.groove2 == b_mut.groove2
    a_mut2, b_mut2 = cat.resolve_pair("HLA-DRB1*01:01 C30S mutant")
    assert b_mut2.groove2 == b_mut.groove2

    # Pan-class features reflect the exact C30S substitution
    f_wt_ii = cat.pan_class_features("HLA-DRA*01:01/DRB1*01:01")
    f_mut_ii = cat.pan_class_features("HLA-DRA*01:01/DRB1*01:01 C30S mutant")
    assert f_wt_ii["seq"] != f_mut_ii["seq"]
    diffs = [(i, a, b) for i, (a, b) in enumerate(zip(f_wt_ii["seq"], f_mut_ii["seq"])) if a != b]
    assert len(diffs) == 1 and diffs[0][1:] == ("C", "S")

    # Wrong wild-type residue raises ValueError instead of silently corrupting
    with pytest.raises(ValueError):
        cat.resolve_pair("HLA-DRA*01:01/DRB1*01:01 A30S mutant")
