import csv
import json
from collections import Counter

import pytest

from mhcseqs.domain_parsing import analyze_sequence, decompose_domains
from scripts.enrich_sp_ground_truth import (
    GT_ENRICHED_CSV,
    NEGATIVE_CONTROL_CSV,
    _apply_label_curation,
    _build_controls,
    _gene_from_protein_name,
    _load_label_curation,
)
from scripts.evaluate_sp_ground_truth import (
    GT_RAW_CSV,
    _parser_name_for_dispatch,
    _row_benchmark_stratum,
    _row_dispatch_metadata,
    _species_category,
    _try_parse,
    predict_sp_for_row,
)
from scripts.sp_ground_truth_taxonomy import (
    SOURCE_CLADES,
    category_from_lineage,
    load_taxonomy_audit,
    source_clade_from_lineage,
)
from scripts.train_sp_boundary_model import OUT_JSON, train


def test_species_category_uses_audited_taxon_lineage():
    assert _species_category("Sparus aurata", "8175") == "fish"
    assert _species_category("Astyanax mexicanus", "7994") == "fish"
    assert _species_category("Sinocyclocheilus rhinocerous", "307959") == "fish"
    assert _species_category("Athene cunicularia", "194338") == "bird"


def test_species_category_uses_exact_taxon_id_when_available():
    assert _species_category("Unknown vertebrate", "9606") == "human"


def test_species_category_uses_source_clade_without_genus_hints():
    assert _species_category("Salvator merianae", "96440") == "other_vertebrate"
    assert _species_category("Tor putitora", "", "Actinopterygii", audit={}) == "fish"
    assert _species_category("Unknown bird", "", "Aves", audit={}) == "bird"
    assert _species_category("Unknown squamate", "", "Lepidosauria", audit={}) == "other_vertebrate"
    with pytest.raises(ValueError, match="Unknown SP source clade"):
        _species_category("Unknown vertebrate", "", "made-up", audit={})


def test_lineage_rules_preserve_established_mammal_categories():
    assert category_from_lineage(9606, {9443, 40674}) == ("human", 9606)
    assert category_from_lineage(10089, {39107, 40674}) == ("murine", 39107)
    assert category_from_lineage(9913, {9895, 40674}) == ("ungulate", 9895)
    assert category_from_lineage(9615, {9608, 40674}) == ("carnivore", 9608)
    assert category_from_lineage(9739, {40674}) == ("other_mammal", 40674)
    assert category_from_lineage(9796, {9788, 40674}) == ("ungulate", 9788)
    assert source_clade_from_lineage(7994, {7898}) == ("Actinopterygii", 7898)
    assert source_clade_from_lineage(8508, {8504}) == ("Lepidosauria", 8504)


def test_human_subspecies_resolve_through_the_lineage_not_an_exact_id():
    """Homo sapiens subspecies must stay human rather than falling back to nhp."""
    assert category_from_lineage(63221, {9606, 9443, 40674}) == ("human", 9606)


def test_source_clade_roots_are_mutually_exclusive():
    """A taxon may only ever resolve to one benchmark fetch clade."""
    roots = [taxon_id for _name, taxon_id in SOURCE_CLADES]
    assert len(set(roots)) == len(roots)
    for name, taxon_id in SOURCE_CLADES:
        assert source_clade_from_lineage(taxon_id, set()) == (name, taxon_id)


def test_taxonomy_audit_covers_every_raw_and_enriched_row():
    taxonomy = load_taxonomy_audit()
    with GT_RAW_CSV.open() as handle:
        raw_rows = list(csv.DictReader(handle))
    with GT_ENRICHED_CSV.open() as handle:
        enriched_rows = list(csv.DictReader(handle))

    assert len(taxonomy) == 253
    assert {row["taxonomy_release"] for row in taxonomy.values()} == {"2026_03"}
    assert {row["source_clade"] for row in taxonomy.values()} == {name for name, _taxon_id in SOURCE_CLADES}
    assert len(raw_rows) == 2403
    assert len(enriched_rows) == 2402
    assert {row["accession"] for row in enriched_rows} < {row["accession"] for row in raw_rows}

    expected_clade_counts = {
        "Mammalia": 497,
        "Aves": 500,
        "Actinopterygii": 500,
        "Lepidosauria": 470,
        "Amphibia": 359,
        "Chondrichthyes": 77,
    }
    assert {clade: sum(row["source_clade"] == clade for row in raw_rows) for clade, _taxon_id in SOURCE_CLADES} == expected_clade_counts

    for row in raw_rows:
        assert row["source_clade"] == taxonomy[row["taxon_id"]]["source_clade"]
    for row in enriched_rows:
        decision = taxonomy[row["taxon_id"]]
        assert row["source_clade"] == decision["source_clade"]
        assert row["species_category"] == decision["species_category"]

    assert Counter(row["species_category"] for row in enriched_rows) == {
        "other_vertebrate": 829,
        "fish": 577,
        "bird": 500,
        "human": 257,
        "nhp": 102,
        "murine": 79,
        "other_mammal": 36,
        "ungulate": 17,
        "carnivore": 5,
    }


def test_stored_boundary_model_matches_lineage_backed_training_groups():
    with OUT_JSON.open() as handle:
        stored = json.load(handle)
    assert {group: payload["n_positive"] for group, payload in stored["groups"].items()} == {
        "bird": 500,
        "fish": 577,
        "global": 2402,
        "mammal": 496,
        "other_vertebrate": 829,
    }
    assert stored == train()


def test_row_dispatch_metadata_reads_enriched_fields():
    row = {"mhc_class": "II", "chain": "beta", "gene": "DRB1"}
    assert _row_dispatch_metadata(row) == ("II", "beta", "DRB1")


def test_parser_name_for_dispatch_maps_class_and_chain():
    assert _parser_name_for_dispatch("I", "alpha") == "class_I"
    assert _parser_name_for_dispatch("II", "alpha") == "class_II_alpha"
    assert _parser_name_for_dispatch("II", "beta") == "class_II_beta"


def test_label_curation_audits_every_incomplete_row():
    curation = _load_label_curation()
    assert len(curation) == 264
    assert sum(row["label_status"] == "curated" for row in curation.values()) == 247
    assert sum(row["label_status"] == "excluded_non_mhc" for row in curation.values()) == 16
    assert sum(row["label_status"] == "unresolved" for row in curation.values()) == 1
    assert curation["A0A8D0CD71"]["mhc_class"] == "II"
    assert curation["A0A8D0CD71"]["chain"] == "alpha"
    assert curation["A0A8D0CD71"]["source_url"].endswith("A0A8D0CD71?format=txt&versions=21")
    assert {"Q08334", "Q61190", "A0A401NJB8"} <= {accession for accession, row in curation.items() if row["label_status"] == "excluded_non_mhc"}
    assert all(row["evidence"] and row["source_version"] and row["source_url"] and row["audited_on"] for row in curation.values())


def test_enriched_ground_truth_contains_every_curation_decision():
    curation = _load_label_curation()
    with GT_ENRICHED_CSV.open() as handle:
        enriched = {row["accession"]: row for row in csv.DictReader(handle)}
    assert set(curation) <= set(enriched)
    for accession, decision in curation.items():
        row = enriched[accession]
        assert row["label_status"] == decision["label_status"]
        assert row["mhc_class"] == decision["mhc_class"]
        assert row["chain"] == decision["chain"]

    assert Counter(row["label_status"] for row in enriched.values()) == {
        "gold": 2138,
        "curated": 247,
        "excluded_non_mhc": 16,
        "unresolved": 1,
    }


def test_controls_only_use_complete_curated_or_gold_labels():
    with GT_ENRICHED_CSV.open() as handle:
        rows = list(csv.DictReader(handle))
    controls = _build_controls(rows)
    assert len(controls) == 2385
    eligible = {row["accession"] for row in rows if _row_benchmark_stratum(row) == "curated/gold"}
    assert {row["source_accession"] for row in controls} <= eligible
    with NEGATIVE_CONTROL_CSV.open() as handle:
        stored_controls = list(csv.DictReader(handle))
    assert [row["control_id"] for row in stored_controls] == [row["control_id"] for row in controls]
    taxonomy = load_taxonomy_audit()
    for row in stored_controls:
        decision = taxonomy[row["taxon_id"]]
        assert row["source_clade"] == decision["source_clade"]
        assert row["species_category"] == decision["species_category"]


def test_apply_label_curation_uses_source_backed_decision():
    decision = _load_label_curation()["A0A8D0CD71"]
    row = {
        "accession": "A0A8D0CD71",
        "mhc_class": "",
        "chain": "",
        "label_status": "unresolved",
    }
    result = _apply_label_curation(row, {"A0A8D0CD71": decision})
    assert _row_dispatch_metadata(result)[:2] == ("II", "alpha")
    assert _row_benchmark_stratum(result) == "curated/gold"
    assert "IPR001003" in decision["evidence"]
    assert "UniProtKB entry" in decision["source_version"]


def test_excluded_non_mhc_rows_do_not_enter_benchmark_denominators():
    row = _apply_label_curation(
        {
            "accession": "Q08334",
            "mhc_class": "",
            "chain": "",
            "label_status": "unresolved",
        },
        _load_label_curation(),
    )
    assert _row_benchmark_stratum(row) == "excluded"
    row["sequence"] = "NOT_A_SEQUENCE"
    assert predict_sp_for_row(row)["dispatch_mode"] == "excluded"


def test_classless_evaluator_delegates_to_production_whole_parse():
    with GT_ENRICHED_CSV.open() as handle:
        row = next(row for row in csv.DictReader(handle) if row["accession"] == "A0A8D0CD71")
    features = analyze_sequence(row["sequence"])
    production = decompose_domains(row["sequence"], mhc_class="", features=features)
    mature_start, parser_name = _try_parse(row["sequence"], features=features)
    assert mature_start == production.mature_start
    assert parser_name == _parser_name_for_dispatch(production.mhc_class, production.chain)


def test_source_indeterminate_row_is_reported_as_inferred():
    with GT_ENRICHED_CSV.open() as handle:
        row = next(row for row in csv.DictReader(handle) if row["accession"] == "A0A8B9VS73")
    assert _row_benchmark_stratum(row) == "unresolved/inferred"
    assert predict_sp_for_row(row)["dispatch_mode"] == "inferred"


def test_gene_from_protein_name_handles_class_ii_arm_tokens():
    gene, status = _gene_from_protein_name(
        "HLA class II histocompatibility antigen, DR alpha chain-like",
        "Naja naja",
        "II",
        "alpha",
    )
    assert gene == "NajaNaja-DRA"
    assert status == "protein_name_token"


def test_gene_from_protein_name_canonicalizes_hla_parenthetical_names():
    gene, status = _gene_from_protein_name(
        "HLA class II histocompatibility antigen, DR alpha chain (MHC class II antigen HLA-DRA)",
        "Homo sapiens",
        "II",
        "alpha",
    )
    assert gene == "HomoSapiens-DRA"
    assert status == "protein_name_parenthetical"
