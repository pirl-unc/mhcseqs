import io
import json
import urllib.request
from pathlib import Path

import pytest

from scripts import (
    audit_sp_ground_truth_taxonomy,
    enrich_sp_ground_truth,
    fetch_sp_ground_truth,
    sp_corpus_artifacts,
)
from scripts.sp_ground_truth_taxonomy import CATEGORY_LINEAGE_RULES, SOURCE_CLADES

artifacts = sp_corpus_artifacts


UNISAVE_TEXT = """ID   TEST_ENTRY              Unreviewed;       120 AA.
AC   TEST123;
DT   01-JAN-2020, integrated into UniProtKB/TrEMBL.
DT   02-JAN-2020, sequence version 1.
DT   03-JAN-2020, entry version 2.
DE   SubName: Full=MHC class II antigen alpha chain {ECO:0000313|EMBL:ABC.1};
GN   Name=DRA {ECO:0000313|EMBL:ABC.1};
OS   Example vertebrate.
OC   Eukaryota; Metazoa; Chordata; Vertebrata; Mammalia.
OX   NCBI_TaxID=12345;
PE   3: Inferred from homology;
KW   Immunity; MHC II; Signal.
FT   SIGNAL          1..20
FT                   /evidence="ECO:0000256|SAM:SignalP"
FT   CHAIN           21..120
FT                   /id="PRO_123"
SQ   SEQUENCE   120 AA;  12000 MW;  ABC CRC64;
     MAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
     AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
//
"""

METADATA = {
    "accession": "TEST123",
    "entryVersion": 2,
    "lastRelease": "2026_01/2026_01",
    "lastReleaseDate": "28-Jan-2026",
}

TAXONOMY = {
    "12345": {
        "Taxon ID": "12345",
        "Scientific name": "Example vertebrate",
        "Common name": "",
        "Taxonomic lineage (Ids)": "7742, 40674",
        "Taxonomic lineage": "Vertebrata, Mammalia",
        "Taxonomy release": "2026_03",
        "Source URL": "https://rest.uniprot.org/taxonomy/12345",
        "Retrieved on": "2026-09-03",
    }
}


class Response(io.BytesIO):
    def __init__(self, payload: bytes, release: str = "2026_03"):
        super().__init__(payload)
        self.headers = {"X-UniProt-Release": release}


def test_parse_unisave_record_preserves_boundary_sequence_and_provenance():
    row = artifacts.parse_unisave_record(
        UNISAVE_TEXT,
        metadata=METADATA,
        source_url="https://rest.uniprot.org/unisave/TEST123?format=txt&versions=2",
        source_clade="Mammalia",
    )

    assert row["Entry"] == "TEST123"
    assert row["Organism"] == "Example vertebrate"
    assert row["Organism (ID)"] == "12345"
    assert row["Gene Names"] == "DRA"
    assert row["Length"] == "120"
    assert len(row["Sequence"]) == 120
    assert row["Signal peptide"].startswith("SIGNAL 1..20")
    assert row["Chain"].startswith("CHAIN 21..120")
    assert row["Source clade"] == "Mammalia"
    assert row["Archive release"] == "2026_01"
    assert row["Archive release date"] == "28-Jan-2026"
    assert row["Date of last modification"] == "2020-01-03"
    assert row["Entry version"] == "2"


@pytest.mark.parametrize("flag", ["Fragment", "Fragments"])
def test_parse_unisave_record_reads_fragment_flag_from_de_metadata(flag):
    text = UNISAVE_TEXT.replace(
        "DE   SubName:",
        f"DE   Flags: {flag};\nDE   SubName:",
    )
    row = artifacts.parse_unisave_record(text, metadata=METADATA, source_url="https://example.test/archive")
    assert row["Fragment"] == "fragment"


def test_download_deleted_artifact_fetches_selected_flat_file(monkeypatch, tmp_path):
    requested_urls = []

    def fake_request(url, *, timeout):
        requested_urls.append((url, timeout))
        if url.endswith("?format=json"):
            return Response(json.dumps({"results": [METADATA]}).encode())
        if url.endswith("?format=txt&versions=2"):
            return Response(UNISAVE_TEXT.encode())
        raise AssertionError(url)

    monkeypatch.setattr(artifacts, "_request", fake_request)
    artifact = artifacts.download_deleted_artifact(
        ["TEST123"],
        tmp_path,
        source_clades={"TEST123": "Mammalia"},
    )

    assert artifact.rows == 1
    assert artifact.release == "2026_03"
    assert [row["Entry"] for row in artifacts.read_artifact(artifact.path)] == ["TEST123"]
    assert [url for url, _timeout in requested_urls] == [
        "https://rest.uniprot.org/unisave/TEST123?format=json",
        "https://rest.uniprot.org/unisave/TEST123?format=txt&versions=2",
    ]


def test_download_deleted_artifact_rejects_missing_history(monkeypatch, tmp_path):
    def fake_request(_url, *, timeout):
        return Response(b'{"results": []}')

    monkeypatch.setattr(artifacts, "_request", fake_request)
    with pytest.raises(ValueError, match="UniSave contains no versions for MISSING"):
        artifacts.download_deleted_artifact(["MISSING"], tmp_path)
    assert not (artifacts.artifact_dir(tmp_path) / artifacts.DELETED_ARTIFACT).exists()


def test_every_exhaustive_clade_is_reachable_from_artifact_lineage():
    assert artifacts.CORPUS_CLADES is SOURCE_CLADES
    for name, taxon_id in artifacts.CORPUS_CLADES:
        row = {
            "Entry": name,
            "Organism (ID)": "999999",
            "Taxonomic lineage (Ids)": f"7742, {taxon_id}",
        }
        assert artifacts.source_clade_from_artifact_row(row) == name


def test_artifact_row_to_ground_truth_uses_exact_signal_and_stored_clade():
    row = artifacts.parse_unisave_record(
        UNISAVE_TEXT,
        metadata=METADATA,
        source_url="https://example.test/archive",
        source_clade="Mammalia",
    )
    result = fetch_sp_ground_truth.artifact_row_to_ground_truth(row, TAXONOMY)

    assert result == {
        "accession": "TEST123",
        "organism": "Example vertebrate",
        "taxon_id": "12345",
        "source_clade": "Mammalia",
        "sp_length": "20",
        "reviewed": "N",
        "sequence": "M" + "A" * 119,
    }


def test_raw_organism_uses_snapshot_scientific_name_not_tsv_display_label():
    row = artifacts.parse_unisave_record(
        UNISAVE_TEXT,
        metadata=METADATA,
        source_url="https://example.test/archive",
        source_clade="Mammalia",
    )
    row["Organism"] = "Example vertebrate (Display label)"
    assert fetch_sp_ground_truth.artifact_row_to_ground_truth(row, TAXONOMY)["organism"] == "Example vertebrate"


@pytest.mark.parametrize("feature", ["", "SIGNAL <1..20", "SIGNAL 1..?20"])
def test_artifact_row_to_ground_truth_rejects_missing_or_fuzzy_boundary(feature):
    row = {
        "Entry": "TEST123",
        "Organism": "Example vertebrate",
        "Organism (ID)": "12345",
        "Source clade": "Mammalia",
        "Protein names": "MHC class II antigen alpha chain",
        "Gene Names": "DRA",
        "Reviewed": "unreviewed",
        "Sequence": "M" + "A" * 119,
        "Signal peptide": feature,
    }
    assert fetch_sp_ground_truth.artifact_row_to_ground_truth(row, TAXONOMY) is None


def test_corpus_derivation_separates_labels_inference_and_exclusions():
    labelled = artifacts.parse_unisave_record(
        UNISAVE_TEXT,
        metadata=METADATA,
        source_url="https://example.test/archive",
        source_clade="Mammalia",
    )
    unlabelled = {**labelled, "Entry": "UNLABELLED", "Signal peptide": ""}
    fragment = {**labelled, "Entry": "FRAGMENT", "Fragment": "fragment"}
    non_mhc = {
        **labelled,
        "Entry": "NOTMHC",
        "Protein names": "Interleukin-10 receptor subunit beta (Cytokine receptor class-II member 4)",
        "Gene Names": "IL10RB",
    }

    gt_rows, unlabelled_rows, stats = fetch_sp_ground_truth.derive_corpus_rows(
        [labelled, unlabelled, fragment, non_mhc],
        TAXONOMY,
        curation={"NOTMHC": {"disposition": "exclude_non_mhc"}},
    )

    assert [row["accession"] for row in gt_rows] == ["TEST123"]
    assert [row["accession"] for row in unlabelled_rows] == ["UNLABELLED"]
    assert stats == {
        "candidates": 4,
        "labelled": 1,
        "unlabelled": 1,
        "excluded_fragment": 1,
        "excluded_non_mhc": 1,
    }


def test_failed_refresh_preserves_entire_published_bundle(monkeypatch, tmp_path):
    published = artifacts.artifact_dir(tmp_path)
    published.mkdir(parents=True)
    old = {}
    for name in (*artifacts.REQUIRED_ARTIFACTS, artifacts.RELEASE_MANIFEST):
        old[name] = f"old-{name}".encode()
        (published / name).write_bytes(old[name])

    def fake_download_corpus(staged_data_dir):
        path = artifacts.artifact_dir(staged_data_dir) / artifacts.CORPUS_ARTIFACT
        path.parent.mkdir(parents=True)
        path.write_bytes(b"new-corpus")
        return artifacts.Artifact(artifacts.CORPUS_ARTIFACT, path, "2026_03", 1)

    monkeypatch.setattr(fetch_sp_ground_truth, "download_corpus_artifact", fake_download_corpus)
    monkeypatch.setattr(fetch_sp_ground_truth, "load_corpus_artifact", lambda _data_dir: [{"Entry": "CURRENT"}])

    def fail_deleted(_accessions, _data_dir, **_kwargs):
        raise TimeoutError("simulated UniSave failure")

    monkeypatch.setattr(fetch_sp_ground_truth, "download_deleted_artifact", fail_deleted)
    with pytest.raises(TimeoutError, match="simulated UniSave failure"):
        fetch_sp_ground_truth.refresh_artifacts(
            [{"accession": "DELETED", "source_clade": "Aves"}],
            tmp_path,
        )

    assert {name: (published / name).read_bytes() for name in old} == old


def _taxonomy_rows() -> list[dict[str, str]]:
    rows = []
    names = {taxon_id: name for name, taxon_id in SOURCE_CLADES}
    required = {taxon_id for _name, taxon_id in SOURCE_CLADES}
    required.update(root for _category, roots in CATEGORY_LINEAGE_RULES for root in roots)
    required.add(12345)
    for taxon_id in sorted(required):
        lineage_ids = "7742, 40674" if taxon_id == 12345 else ""
        rows.append(
            {
                "Taxon ID": str(taxon_id),
                "Scientific name": "Example vertebrate" if taxon_id == 12345 else names.get(taxon_id, f"Taxon {taxon_id}"),
                "Common name": "",
                "Taxonomic lineage (Ids)": lineage_ids,
                "Taxonomic lineage": "Vertebrata, Mammalia" if taxon_id == 12345 else "",
                "Taxonomy release": "2026_03",
                "Source URL": f"https://rest.uniprot.org/taxonomy/{taxon_id}",
                "Retrieved on": "2026-09-03",
            }
        )
    return rows


def _write_fixture_bundle(data_dir: Path) -> None:
    parsed = artifacts.parse_unisave_record(UNISAVE_TEXT, metadata=METADATA, source_url="https://example.test/current")
    current = {column: parsed.get(column, "") for column in artifacts.ARTIFACT_COLUMNS}
    current["Taxonomic lineage (Ids)"] = "7742, 40674"
    bundle_dir = artifacts.artifact_dir(data_dir)
    artifacts._write_gzip_tsv(bundle_dir / artifacts.CORPUS_ARTIFACT, [current], artifacts.ARTIFACT_COLUMNS)
    artifacts._write_gzip_tsv(bundle_dir / artifacts.DELETED_ARTIFACT, [], artifacts.DELETED_COLUMNS)
    taxonomy_rows = _taxonomy_rows()
    artifacts._write_gzip_tsv(bundle_dir / artifacts.TAXONOMY_ARTIFACT, taxonomy_rows, artifacts.TAXONOMY_COLUMNS)
    artifact_records = [
        artifacts.Artifact(artifacts.CORPUS_ARTIFACT, bundle_dir / artifacts.CORPUS_ARTIFACT, "2026_03", 1),
        artifacts.Artifact(artifacts.DELETED_ARTIFACT, bundle_dir / artifacts.DELETED_ARTIFACT, "2026_03", 0),
        artifacts.Artifact(artifacts.TAXONOMY_ARTIFACT, bundle_dir / artifacts.TAXONOMY_ARTIFACT, "2026_03", len(taxonomy_rows)),
    ]
    artifacts.write_manifest(artifact_records, data_dir)
    artifacts.validate_artifact_bundle(data_dir)


def test_full_regeneration_replays_without_network(monkeypatch, tmp_path):
    _write_fixture_bundle(tmp_path)
    raw = tmp_path / "raw.csv"
    unlabelled = tmp_path / "unlabelled.csv"
    enriched = tmp_path / "enriched.csv"
    controls = tmp_path / "controls.csv"
    audit = tmp_path / "taxonomy.csv"
    missing_curation = tmp_path / "no-curation.csv"

    def deny_network(*_args, **_kwargs):
        raise AssertionError("offline regeneration attempted network access")

    monkeypatch.setattr(urllib.request, "urlopen", deny_network)
    fetch_sp_ground_truth.main(
        [
            "--data-dir",
            str(tmp_path),
            "--output",
            str(raw),
            "--unlabelled-output",
            str(unlabelled),
            "--label-curation",
            str(missing_curation),
        ]
    )
    enrich_sp_ground_truth.main(
        [
            "--data-dir",
            str(tmp_path),
            "--raw",
            str(raw),
            "--output",
            str(enriched),
            "--controls",
            str(controls),
            "--label-curation",
            str(missing_curation),
        ]
    )
    audit_sp_ground_truth_taxonomy.main(
        [
            "--data-dir",
            str(tmp_path),
            "--raw",
            str(raw),
            "--enriched",
            str(enriched),
            "--controls",
            str(controls),
            "--output",
            str(audit),
        ]
    )

    first_run = {path: path.read_bytes() for path in (raw, unlabelled, enriched, controls, audit)}
    fetch_sp_ground_truth.main(
        [
            "--data-dir",
            str(tmp_path),
            "--output",
            str(raw),
            "--unlabelled-output",
            str(unlabelled),
            "--label-curation",
            str(missing_curation),
        ]
    )
    enrich_sp_ground_truth.main(
        [
            "--data-dir",
            str(tmp_path),
            "--raw",
            str(raw),
            "--output",
            str(enriched),
            "--controls",
            str(controls),
            "--label-curation",
            str(missing_curation),
        ]
    )
    audit_sp_ground_truth_taxonomy.main(
        [
            "--data-dir",
            str(tmp_path),
            "--raw",
            str(raw),
            "--enriched",
            str(enriched),
            "--controls",
            str(controls),
            "--output",
            str(audit),
        ]
    )

    assert {path: path.read_bytes() for path in first_run} == first_run
    assert "Example vertebrate" in raw.read_text()
    assert "uniprot_corpus_artifact" in enriched.read_text()
    assert "2026_03" in audit.read_text()
    manifest = artifacts.read_manifest(tmp_path)
    assert manifest["provenance"]["current"]["query"] == artifacts.corpus_query()


def test_refresh_flag_is_supported():
    args = fetch_sp_ground_truth._parse_args(["--refresh"])
    assert args.refresh is True
