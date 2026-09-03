import csv
import json
from collections import Counter

import pytest

from mhcseqs.domain_parsing import analyze_sequence, decompose_domains
from scripts.audit_sp_ground_truth_taxonomy import _apply_audit
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
from scripts.fetch_sp_ground_truth import UNLABELLED_OUTPUT
from scripts.sp_ground_truth_taxonomy import (
    CATEGORY_LINEAGE_RULES,
    MAMMAL_CATEGORIES,
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
    assert _species_category("Unknown turtle", "", "Testudines", audit={}) == "other_vertebrate"
    assert _species_category("Unknown lungfish", "", "Dipnomorpha", audit={}) == "fish"
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
    assert source_clade_from_lineage(8469, {8459}) == ("Testudines", 8459)
    assert source_clade_from_lineage(13397, {1294634}) == ("Crocodylia", 1294634)


def test_human_subspecies_resolve_through_the_lineage_not_an_exact_id():
    """Homo sapiens subspecies must stay human rather than falling back to nhp."""
    assert category_from_lineage(63221, {9606, 9443, 40674}) == ("human", 9606)


def test_category_rule_order_is_specific_to_general():
    """Rule order decides the answer, so pin it: a later broad root hides earlier ones."""
    assert [category for category, _roots in CATEGORY_LINEAGE_RULES] == [
        "human",
        "nhp",
        "murine",
        "ungulate",
        "carnivore",
        "other_mammal",
        "bird",
        "fish",
        "other_vertebrate",
    ]
    # Mammalia matches every mammal, so no mammal rule may follow it.
    categories = [category for category, _roots in CATEGORY_LINEAGE_RULES]
    after_other_mammal = set(categories[categories.index("other_mammal") + 1 :])
    assert not (after_other_mammal & MAMMAL_CATEGORIES)


def test_every_category_rule_is_reachable():
    """A rule shadowed by an earlier root would be silent dead code."""
    for category, roots in CATEGORY_LINEAGE_RULES:
        for root in roots:
            assert category_from_lineage(root, {root}) == (category, root)


def test_missing_taxonomy_audit_fails_loudly(tmp_path):
    """The audit is authoritative; absence must not degrade to name matching."""
    with pytest.raises(FileNotFoundError, match="SP taxonomy audit is missing"):
        load_taxonomy_audit(tmp_path / "absent.csv")


def test_apply_audit_rejects_source_clade_drift():
    """A stored clade disagreeing with the lineage is the signal, not noise."""
    audit = {"9606": {"source_clade": "Mammalia", "species_category": "human"}}
    rows = [{"taxon_id": "9606", "source_clade": "Aves", "species_category": "bird"}]
    with pytest.raises(ValueError, match="Source clade drift"):
        _apply_audit(rows, audit)


def test_apply_audit_fills_fields_without_touching_labels():
    audit = {"9606": {"source_clade": "Mammalia", "species_category": "human"}}
    raw_row = {"taxon_id": "9606", "organism": "Homo sapiens", "sp_length": "21", "source_clade": ""}
    enriched_row = {"taxon_id": "9606", "organism": "Homo sapiens", "source_clade": "", "species_category": "stale", "mhc_class": "I"}
    _apply_audit([raw_row], audit)
    _apply_audit([enriched_row], audit)
    assert raw_row["source_clade"] == "Mammalia"
    assert "species_category" not in raw_row  # raw rows carry no category column
    assert raw_row["organism"] == "Homo sapiens" and raw_row["sp_length"] == "21"
    assert enriched_row["species_category"] == "human"
    assert enriched_row["mhc_class"] == "I"


def test_apply_audit_requires_every_taxon():
    with pytest.raises(ValueError, match="Missing taxonomy audit"):
        _apply_audit([{"taxon_id": "9999", "source_clade": ""}], {})


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

    assert len(taxonomy) == 533
    assert {row["taxonomy_release"] for row in taxonomy.values()} == {"2026_03"}
    assert {row["source_clade"] for row in taxonomy.values()} <= {name for name, _taxon_id in SOURCE_CLADES}
    assert len(raw_rows) == 14721
    assert len(enriched_rows) == 14721
    assert {row["accession"] for row in enriched_rows} == {row["accession"] for row in raw_rows}

    expected_clade_counts = {
        "Mammalia": 9749,
        "Aves": 784,
        "Actinopterygii": 3262,
        "Lepidosauria": 454,
        "Amphibia": 307,
        "Chondrichthyes": 71,
        "Testudines": 64,
        "Crocodylia": 28,
        "Coelacanthimorpha": 1,
        "Dipnomorpha": 1,
        "Myxini": 0,
    }
    assert {clade: sum(row["source_clade"] == clade for row in raw_rows) for clade, _taxon_id in SOURCE_CLADES} == expected_clade_counts

    for row in raw_rows:
        assert row["source_clade"] == taxonomy[row["taxon_id"]]["source_clade"]
    for row in enriched_rows:
        decision = taxonomy[row["taxon_id"]]
        assert row["source_clade"] == decision["source_clade"]
        assert row["species_category"] == decision["species_category"]

    assert Counter(row["species_category"] for row in enriched_rows) == {
        "other_vertebrate": 853,
        "fish": 3335,
        "bird": 784,
        "human": 847,
        "nhp": 5121,
        "murine": 330,
        "other_mammal": 1804,
        "ungulate": 1371,
        "carnivore": 276,
    }


def test_unlabelled_inference_set_is_complete_eligible_and_disjoint():
    with GT_RAW_CSV.open() as handle:
        labelled_accessions = {row["accession"] for row in csv.DictReader(handle)}
    with UNLABELLED_OUTPUT.open() as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 4483
    assert labelled_accessions.isdisjoint(row["accession"] for row in rows)
    assert all(row["sequence"].startswith("M") for row in rows)
    assert {(row["mhc_class"], row["chain"]) for row in rows} <= {
        ("I", "alpha"),
        ("II", "alpha"),
        ("II", "beta"),
    }
    assert Counter(row["source_clade"] for row in rows) == {
        "Mammalia": 1779,
        "Aves": 257,
        "Actinopterygii": 1778,
        "Lepidosauria": 312,
        "Amphibia": 162,
        "Chondrichthyes": 26,
        "Testudines": 102,
        "Crocodylia": 64,
        "Coelacanthimorpha": 1,
        "Myxini": 2,
    }


def test_stored_boundary_model_matches_lineage_backed_training_groups():
    with OUT_JSON.open() as handle:
        stored = json.load(handle)
    assert {group: payload["n_positive"] for group, payload in stored["groups"].items()} == {
        "bird": 784,
        "fish": 3335,
        "global": 14721,
        "mammal": 9749,
        "other_vertebrate": 853,
    }
    assert stored["training"] == {
        "source": "sp_ground_truth_enriched.csv",
        "policy": "label_status in {gold,curated} and MHC class/chain in {I alpha,II alpha,II beta}",
        "source_rows": 14721,
        "eligible_rows": 14721,
        "used_rows": 14721,
        "excluded_rows": 0,
        "unresolved_rows": 0,
    }
    assert stored == train()


def test_boundary_model_training_excludes_non_mhc_and_unresolved_rows(tmp_path):
    path = tmp_path / "training.csv"
    fieldnames = ["accession", "sequence", "sp_length", "species_category", "mhc_class", "chain", "label_status"]
    rows = [
        {
            "accession": "ELIGIBLE",
            "sequence": "M" + "A" * 79,
            "sp_length": "20",
            "species_category": "human",
            "mhc_class": "I",
            "chain": "alpha",
            "label_status": "gold",
        },
        {
            "accession": "CONTAMINANT",
            "sequence": "M" + "A" * 79,
            "sp_length": "20",
            "species_category": "human",
            "mhc_class": "",
            "chain": "",
            "label_status": "excluded_non_mhc",
        },
        {
            "accession": "UNRESOLVED",
            "sequence": "M" + "A" * 79,
            "sp_length": "20",
            "species_category": "human",
            "mhc_class": "II",
            "chain": "unknown",
            "label_status": "unresolved",
        },
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    training = train(path)["training"]
    assert training["source_rows"] == 3
    assert training["eligible_rows"] == 1
    assert training["used_rows"] == 1
    assert training["excluded_rows"] == 1
    assert training["unresolved_rows"] == 1


def test_row_dispatch_metadata_reads_enriched_fields():
    row = {"mhc_class": "II", "chain": "beta", "gene": "DRB1"}
    assert _row_dispatch_metadata(row) == ("II", "beta", "DRB1")


def test_parser_name_for_dispatch_maps_class_and_chain():
    assert _parser_name_for_dispatch("I", "alpha") == "class_I"
    assert _parser_name_for_dispatch("II", "alpha") == "class_II_alpha"
    assert _parser_name_for_dispatch("II", "beta") == "class_II_beta"


def test_label_curation_audits_every_incomplete_row():
    curation = _load_label_curation()
    assert len(curation) == 283
    assert sum(row["label_status"] == "curated" for row in curation.values()) == 266
    assert sum(row["label_status"] == "excluded_non_mhc" for row in curation.values()) == 16
    assert sum(row["label_status"] == "unresolved" for row in curation.values()) == 1
    assert curation["A0A8D0CD71"]["mhc_class"] == "II"
    assert curation["A0A8D0CD71"]["chain"] == "alpha"
    assert curation["A0A8D0CD71"]["source_url"].endswith("A0A8D0CD71?format=txt&versions=21")
    assert curation["Q5TJH0"]["chain"] == "beta"
    assert {"A9UM13", "Q31365", "Q8HWB2"} <= set(curation)
    assert {"Q08334", "Q61190", "A0A401NJB8"} <= {accession for accession, row in curation.items() if row["label_status"] == "excluded_non_mhc"}
    assert all(row["evidence"] and row["source_version"] and row["source_url"] and row["audited_on"] for row in curation.values())


def test_enriched_ground_truth_contains_every_curation_decision():
    curation = _load_label_curation()
    with GT_ENRICHED_CSV.open() as handle:
        enriched = {row["accession"]: row for row in csv.DictReader(handle)}
    for accession in curation.keys() & enriched.keys():
        decision = curation[accession]
        row = enriched[accession]
        assert row["label_status"] == decision["label_status"]
        assert row["mhc_class"] == decision["mhc_class"]
        assert row["chain"] == decision["chain"]

    assert {accession for accession, row in curation.items() if row["disposition"] != "include"}.isdisjoint(enriched)
    assert {accession for accession, row in curation.items() if row["disposition"] == "include"} - set(enriched) == {
        "A0A093GBE9",
        "A0A094N606",
        "A0A2G9R3I4",
        "A0A8J6EIP1",
        "Q0P3R9",
        "Q561L0",
        "Q56A59",
        "Q5PPV7",
        "Q5U4N3",
        "Q6GP56",
        "Q7SYW8",
        "Q7SZ32",
    }
    assert Counter(row["label_status"] for row in enriched.values()) == {
        "gold": 14467,
        "curated": 254,
    }
    assert Counter(row["metadata_source"] for row in enriched.values()) == {
        "uniprot_corpus_artifact": 14531,
        "unisave_artifact": 190,
    }


def test_controls_only_use_complete_curated_or_gold_labels():
    with GT_ENRICHED_CSV.open() as handle:
        rows = list(csv.DictReader(handle))
    controls = _build_controls(rows)
    assert len(controls) == 14721
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


def test_source_indeterminate_row_is_not_in_the_labelled_corpus():
    with GT_ENRICHED_CSV.open() as handle:
        accessions = {row["accession"] for row in csv.DictReader(handle)}
    assert "A0A8B9VS73" not in accessions


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
