import csv
import gzip
import json
import urllib.request

from mhcseqs.domain_parsing import AlleleRecord
from scripts import build_mhc_protein_dataset as dataset
from tests.test_sp_corpus_artifacts import _write_fixture_bundle


def _artifact(**overrides):
    row = {
        "Entry": "TEST123",
        "Entry Name": "TEST_ENTRY",
        "Organism": "Example vertebrate (display label)",
        "Organism (ID)": "12345",
        "Taxonomic lineage (Ids)": "7742, 40674",
        "Taxonomic lineage": "Vertebrata, Mammalia",
        "Protein names": "MHC class I antigen alpha chain",
        "Gene Names": "HLA-A",
        "Reviewed": "reviewed",
        "Protein existence": "Evidence at protein level",
        "Fragment": "",
        "Length": "120",
        "Sequence": "M" + "A" * 119,
        "Signal peptide": 'SIGNAL 1..20; /evidence="ECO:0000256"',
        "Chain": 'CHAIN 21..120; /id="PRO_123"',
        "Keywords": "Immunity; MHC I",
        "Date of last modification": "2026-08-01",
        "Entry version": "7",
        "Source clade": "",
        "Archive release": "",
        "Archive release date": "",
        "Archive source URL": "",
    }
    row.update(overrides)
    return row


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

CURATION = {
    "TEST123": {
        "accession": "TEST123",
        "disposition": "include",
        "mhc_class": "I",
        "chain": "alpha",
        "label_status": "curated",
    }
}


def test_prepare_record_separates_nullable_source_and_inferred_annotations():
    row = dataset.prepare_record(
        _artifact(),
        source_release="2026_03",
        version="uniprot-2026_03-r1",
        taxonomy=TAXONOMY,
        curation=CURATION,
    )

    assert row["dataset_version"] == "uniprot-2026_03-r1"
    assert row["organism"] == "Example vertebrate"
    assert row["taxonomic_lineage"] == "Vertebrata, Mammalia"
    assert row["source_signal_feature"].startswith("SIGNAL 1..20")
    assert row["source_signal_status"] == "exact"
    assert (row["source_signal_start"], row["source_signal_end"]) == ("0", "20")
    assert row["inferred_protein_type"] == "mhc_i_alpha"
    assert row["inferred_disposition"] == "include"
    assert row["_parse_mhc_class"] == "I"
    assert row["_parse_use_early_shortcuts"] == "false"


def test_prepare_record_preserves_records_without_signal_annotations():
    row = dataset.prepare_record(
        _artifact(**{"Signal peptide": "", "Fragment": "fragment"}),
        source_release="2026_03",
        version="uniprot-2026_03-r1",
        taxonomy=TAXONOMY,
        curation=CURATION,
    )
    assert row["source_signal_status"] == "not_annotated"
    assert row["source_signal_start"] == row["source_signal_end"] == ""
    assert row["source_fragment"] == "true"
    assert row["sequence"] == "M" + "A" * 119


def test_source_signal_fuzzy_coordinates_remain_nullable():
    fields = dataset._source_signal_fields("SIGNAL <1..>20")
    assert fields["source_signal_status"] == "fuzzy"
    assert fields["source_signal_start"] == fields["source_signal_end"] == ""


def test_parser_fields_expose_boundary_parts_and_status_separately():
    record = AlleleRecord(
        allele="TEST123",
        mhc_class="I",
        chain="alpha",
        sequence="M" * 10 + "A" * 90,
        seq_len=100,
        mature_start=10,
        groove1="A" * 30,
        groove2="A" * 30,
        groove_seq="A" * 60,
        groove1_len=30,
        groove2_len=30,
        ig_domain="A" * 20,
        ig_domain_len=20,
        tail="A" * 10,
        tail_len=10,
        parse_score=4.5,
        status="ok",
        anchor_type="alpha2_cys",
        anchor_cys1=42,
        anchor_cys2=105,
        flags=("example",),
    )
    fields = dataset._parser_fields(record)
    assert fields["parse_ok"] == "true"
    assert fields["parsed_protein_type"] == "mhc_i_alpha"
    assert fields["parsed_signal_status"] == "present"
    assert fields["parsed_signal_end"] == "10"
    assert fields["parsed_signal_sequence"] == "M" * 10
    assert fields["groove1_len"] == fields["groove2_len"] == "30"
    assert fields["parse_flags"] == "example"


def test_parser_fields_preserve_partial_parts_when_parse_is_not_final():
    record = AlleleRecord(
        allele="TEST123",
        mhc_class="II",
        chain="alpha",
        sequence="M" * 5 + "A" * 25,
        seq_len=30,
        mature_start=5,
        groove1="A" * 25,
        groove_seq="A" * 25,
        groove1_len=25,
        status="short",
    )
    fields = dataset._parser_fields(record)
    assert fields["parse_ok"] == "false"
    assert fields["parsed_signal_end"] == "5"
    assert fields["mature_sequence"] == "A" * 25
    assert fields["groove1"] == "A" * 25
    assert fields["groove1_len"] == "25"


def test_failed_parse_without_boundary_does_not_claim_no_signal_peptide():
    record = AlleleRecord(sequence="MAAA", seq_len=4, status="too_short")
    fields = dataset._parser_fields(record)
    assert fields["parse_ok"] == "false"
    assert fields["parsed_signal_status"] == "unresolved"
    assert fields["parsed_signal_end"] == ""
    assert fields["mature_sequence"] == ""


def test_parser_error_keeps_the_source_record(monkeypatch):
    def fail_parse(*_args, **_kwargs):
        raise RuntimeError("simulated parser failure")

    monkeypatch.setattr(dataset, "decompose_domains", fail_parse)
    rows = dataset.build_records(
        [_artifact()],
        source_release="2026_03",
        version="uniprot-2026_03-r1",
        taxonomy=TAXONOMY,
        curation=CURATION,
    )
    assert len(rows) == 1
    assert rows[0]["accession"] == "TEST123"
    assert rows[0]["parse_status"] == "parser_error"
    assert rows[0]["parse_error"] == "RuntimeError: simulated parser failure"
    assert rows[0]["sequence"] == "M" + "A" * 119


def test_gzip_dataset_output_is_byte_deterministic(tmp_path):
    row = {field: "" for field in dataset.PROTEIN_FIELDS}
    row.update({"dataset_version": "uniprot-2026_03-r1", "accession": "TEST123", "sequence": "MAAA"})
    first = tmp_path / "first.csv.gz"
    second = tmp_path / "second.csv.gz"
    dataset._write_gzip_csv(first, [row])
    dataset._write_gzip_csv(second, [row])

    assert first.read_bytes() == second.read_bytes()
    with gzip.open(first, "rt", encoding="utf-8") as handle:
        assert list(csv.DictReader(handle))[0]["accession"] == "TEST123"


def test_dataset_manifest_records_source_models_schema_and_counts(tmp_path):
    row = {field: "" for field in dataset.PROTEIN_FIELDS}
    row.update(
        {
            "dataset_version": "uniprot-2026_03-r1",
            "accession": "TEST123",
            "source_kind": "uniprotkb",
            "source_signal_status": "not_annotated",
            "inferred_disposition": "include",
            "parse_status": "ok",
        }
    )
    output = tmp_path / "records.csv.gz"
    dataset._write_gzip_csv(output, [row])
    curation = tmp_path / "curation.csv"
    curation.write_text("accession,disposition\n")
    source_manifest = {
        "schema_version": 2,
        "release": "2026_03",
        "artifacts": {"source.tsv.gz": {"sha256": "a" * 64}},
        "provenance": {"current": {"query": "example"}},
    }
    manifest_path = tmp_path / "manifest.json"
    payload = dataset.write_dataset_manifest(
        manifest_path,
        output=output,
        rows=[row],
        source_manifest=source_manifest,
        curation_path=curation,
        version="uniprot-2026_03-r1",
        use_early_shortcuts=True,
    )

    assert json.loads(manifest_path.read_text()) == payload
    assert payload["source"]["release"] == "2026_03"
    assert payload["license"]["identifier"] == "CC-BY-4.0"
    assert payload["records"]["rows"] == 1
    assert payload["records"]["columns"] == list(dataset.PROTEIN_FIELDS)
    assert len(payload["records"]["sha256"]) == 64


def test_dataset_regenerates_from_stored_bundle_without_network(monkeypatch, tmp_path):
    _write_fixture_bundle(tmp_path)
    output = tmp_path / "records.csv.gz"
    manifest = tmp_path / "records.manifest.json"
    curation = tmp_path / "curation.csv"
    curation.write_text("accession,disposition,mhc_class,chain,label_status,reason\n")

    def deny_network(*_args, **_kwargs):
        raise AssertionError("dataset regeneration attempted network access")

    monkeypatch.setattr(urllib.request, "urlopen", deny_network)
    dataset.main(
        [
            "--data-dir",
            str(tmp_path),
            "--output",
            str(output),
            "--manifest-output",
            str(manifest),
            "--label-curation",
            str(curation),
            "--workers",
            "1",
        ]
    )

    with gzip.open(output, "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["accession"] for row in rows] == ["TEST123"]
    assert rows[0]["source_signal_status"] == "exact"
    assert rows[0]["sequence"] == "M" + "A" * 119
    assert json.loads(manifest.read_text())["records"]["rows"] == 1
