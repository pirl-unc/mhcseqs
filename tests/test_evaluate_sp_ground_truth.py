from scripts.enrich_sp_ground_truth import _gene_from_protein_name
from scripts.evaluate_sp_ground_truth import (
    _parser_name_for_dispatch,
    _row_dispatch_metadata,
    _species_category,
)


def test_species_category_uses_word_boundaries_for_fallback_hints():
    assert _species_category("Sparus aurata", "8175") == "fish"


def test_species_category_uses_exact_taxon_id_when_available():
    assert _species_category("Unknown vertebrate", "9606") == "human"


def test_row_dispatch_metadata_reads_enriched_fields():
    row = {"mhc_class": "II", "chain": "beta", "gene": "DRB1"}
    assert _row_dispatch_metadata(row) == ("II", "beta", "DRB1")


def test_parser_name_for_dispatch_maps_class_and_chain():
    assert _parser_name_for_dispatch("I", "alpha") == "class_I"
    assert _parser_name_for_dispatch("II", "alpha") == "class_II_alpha"
    assert _parser_name_for_dispatch("II", "beta") == "class_II_beta"


def test_gene_from_protein_name_handles_class_ii_arm_tokens():
    gene, status = _gene_from_protein_name(
        "HLA class II histocompatibility antigen, DR alpha chain-like",
        "Naja naja",
        "II",
        "alpha",
    )
    assert gene == "Nana-DRA"
    assert status == "protein_name_token"


def test_gene_from_protein_name_preserves_hla_parenthetical_names():
    gene, status = _gene_from_protein_name(
        "HLA class II histocompatibility antigen, DR alpha chain (MHC class II antigen HLA-DRA)",
        "Homo sapiens",
        "II",
        "alpha",
    )
    assert gene == "HLA-DRA"
    assert status == "protein_name_parenthetical"
