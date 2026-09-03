import io
import json

import pytest

from scripts import fetch_sp_ground_truth, sp_corpus_artifacts

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
    assert row["Date of last modification"] == "2026-01-28"
    assert row["Entry version"] == "2"


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
    result = fetch_sp_ground_truth.artifact_row_to_ground_truth(row)

    assert result == {
        "accession": "TEST123",
        "organism": "Example vertebrate",
        "taxon_id": "12345",
        "source_clade": "Mammalia",
        "sp_length": "20",
        "reviewed": "N",
        "sequence": "M" + "A" * 119,
    }


@pytest.mark.parametrize("feature", ["", "SIGNAL <1..20", "SIGNAL 1..?20"])
def test_artifact_row_to_ground_truth_rejects_missing_or_fuzzy_boundary(feature):
    row = {
        "Entry": "TEST123",
        "Organism": "Example vertebrate",
        "Organism (ID)": "12345",
        "Source clade": "Mammalia",
        "Reviewed": "unreviewed",
        "Sequence": "M" + "A" * 119,
        "Signal peptide": feature,
    }
    assert fetch_sp_ground_truth.artifact_row_to_ground_truth(row) is None


def test_refresh_archives_every_preserved_accession_missing_from_stream(monkeypatch, tmp_path):
    corpus_path = tmp_path / "corpus.tsv.gz"
    corpus_path.write_bytes(b"corpus")
    deleted_path = tmp_path / "deleted.tsv.gz"
    deleted_path.write_bytes(b"deleted")
    corpus = artifacts.Artifact(artifacts.CORPUS_ARTIFACT, corpus_path, "2026_03", 2)
    deleted = artifacts.Artifact(artifacts.DELETED_ARTIFACT, deleted_path, "2026_03", 1)
    captured = {}

    monkeypatch.setattr(fetch_sp_ground_truth, "download_corpus_artifact", lambda _data_dir: corpus)
    monkeypatch.setattr(
        fetch_sp_ground_truth,
        "load_corpus_artifact",
        lambda _data_dir: [{"Entry": "CURRENT"}, {"Entry": "ALSO_CURRENT"}],
    )

    def fake_download_deleted(accessions, _data_dir, *, source_clades):
        captured["accessions"] = accessions
        captured["source_clades"] = source_clades
        return deleted

    monkeypatch.setattr(fetch_sp_ground_truth, "download_deleted_artifact", fake_download_deleted)
    monkeypatch.setattr(fetch_sp_ground_truth, "write_manifest", lambda _artifacts, _data_dir: tmp_path / "manifest.json")

    fetch_sp_ground_truth.refresh_artifacts(
        [
            {"accession": "CURRENT", "source_clade": "Mammalia"},
            {"accession": "DELETED", "source_clade": "Aves"},
        ],
        tmp_path,
    )

    assert captured == {
        "accessions": ["DELETED"],
        "source_clades": {"DELETED": "Aves"},
    }


def test_refresh_flag_is_supported():
    args = fetch_sp_ground_truth._parse_args(["--refresh"])
    assert args.refresh is True
