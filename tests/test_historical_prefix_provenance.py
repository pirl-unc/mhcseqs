import csv
from collections import Counter
from pathlib import Path

from scripts.recover_historical_prefix_provenance import (
    FOREIGN_NOMENCLATURE,
    GENERIC_LABEL,
    PAIR_DECISIONS,
    SOURCE_PREFIX,
    parse_unisave_record,
)

ROOT = Path(__file__).resolve().parent.parent
SNAPSHOT = ROOT / "data" / "historical_prefix_provenance_20260831.csv"
AUDIT = ROOT / "data" / "species_prefix_audit_pre_full_binomial.csv"
CORPUS = ROOT / "mhcseqs" / "diverse_mhc_sequences.csv"
REGISTRY = ROOT / "mhcseqs" / "mhc_prefix_aliases.csv"


def _read(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_snapshot_resolves_every_historical_pair_with_versioned_source_labels():
    rows = _read(SNAPSHOT)
    assert len(rows) == 239
    assert {(row["species"], row["legacy_prefix"]) for row in rows} == set(PAIR_DECISIONS)
    assert Counter(row["label_classification"] for row in rows) == {
        SOURCE_PREFIX: 170,
        GENERIC_LABEL: 68,
        FOREIGN_NOMENCLATURE: 1,
    }
    assert Counter(row["source_kind"] for row in rows) == {
        "uniprotkb_current": 200,
        "unisave_archived": 39,
    }
    assert all(row["source_gene_names"] for row in rows)
    assert all(row["uniprot_entry_version"].isdigit() for row in rows)
    assert all(row["source_url"].startswith("https://rest.uniprot.org/") for row in rows)
    assert {row["retrieved_on"] for row in rows} == {"2026-08-31"}


def test_recovered_labels_are_packaged_without_changing_canonical_genes():
    snapshot = {row["uniprot_accession"]: row for row in _read(SNAPSHOT)}
    corpus = {row["uniprot_accession"]: row for row in _read(CORPUS)}
    assert snapshot.keys() <= corpus.keys()
    for accession, source in snapshot.items():
        row = corpus[accession]
        assert row["raw_gene_label"] == source["source_gene_names"]
        assert not row["gene"].startswith(source["legacy_prefix"] + "-")


def test_audit_resolves_prefixes_without_promoting_artifacts():
    audit = {(row["species"], row["legacy_prefix"]): row for row in _read(AUDIT)}
    assert all(row["status"] != "unverified_historical" for row in audit.values())
    for pair, decision in PAIR_DECISIONS.items():
        row = audit[pair]
        if decision.classification == SOURCE_PREFIX:
            assert row["status"] == "external_uniprot_label"
            assert row["source_attested_row_count"] == row["row_count"]
        elif decision.classification == GENERIC_LABEL:
            assert row["status"] == "source_label_artifact"
            assert row["source_attested_row_count"] == "0"
        else:
            assert row["status"] == "source_foreign_nomenclature_alias"
            assert row["source_attested_row_count"] == "0"


def test_registry_contains_recovered_prefixes_but_not_generic_tokens():
    registry = {(row["species"], row["prefix"]) for row in _read(REGISTRY)}
    for pair, decision in PAIR_DECISIONS.items():
        if decision.classification == SOURCE_PREFIX:
            assert pair in registry
        else:
            assert pair not in registry


def test_parse_unisave_record_preserves_primary_and_synonym_gene_names():
    text = """ID   EXAMPLE_RAT             Unreviewed;       100 AA.
AC   A0EXAMPLE;
DT   01-JAN-2024, integrated into UniProtKB/TrEMBL.
DT   01-JAN-2024, sequence version 1.
DT   10-JUN-2026, entry version 3.
DE   RecName: Full=Ig-like domain-containing protein {ECO:0000259|PROSITE:PS50835};
GN   Name=RT1-M10-ps4 {ECO:0000313|RGD:1};
GN   Synonyms=H2-M10.3 {ECO:0000313|Ensembl:ENSP1};
//
"""
    assert parse_unisave_record(text) == {
        "Entry": "A0EXAMPLE",
        "Gene Names": "RT1-M10-ps4 H2-M10.3",
        "Protein names": "Ig-like domain-containing protein",
        "Date of last modification": "10-JUN-2026",
        "Entry version": "3",
    }
