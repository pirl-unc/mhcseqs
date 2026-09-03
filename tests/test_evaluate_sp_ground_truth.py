import csv

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
    _parser_name_for_dispatch,
    _row_benchmark_stratum,
    _row_dispatch_metadata,
    _species_category,
    _try_parse,
    predict_sp_for_row,
)


def test_species_category_uses_word_boundaries_for_fallback_hints():
    assert _species_category("Sparus aurata", "8175") == "fish"


def test_species_category_uses_exact_taxon_id_when_available():
    assert _species_category("Unknown vertebrate", "9606") == "human"


def test_species_category_does_not_match_tor_inside_salvator():
    assert _species_category("Salvator merianae", "96440") == "other_vertebrate"
    assert _species_category("Tor putitora", "") == "fish"


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


def test_controls_only_use_complete_curated_or_gold_labels():
    with GT_ENRICHED_CSV.open() as handle:
        rows = list(csv.DictReader(handle))
    controls = _build_controls(rows)
    eligible = {row["accession"] for row in rows if _row_benchmark_stratum(row) == "curated/gold"}
    assert {row["source_accession"] for row in controls} <= eligible
    with NEGATIVE_CONTROL_CSV.open() as handle:
        stored_controls = list(csv.DictReader(handle))
    assert [row["control_id"] for row in stored_controls] == [row["control_id"] for row in controls]


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
