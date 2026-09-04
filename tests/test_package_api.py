from mhcseqs import (
    available_mhc_protein_dataset_versions,
    decompose_class_i,
    decompose_class_ii_alpha,
    decompose_class_ii_beta,
    decompose_domains,
    extract_groove,
    parse_class_i,
    parse_class_ii_alpha,
    parse_class_ii_beta,
    validate_mhc_protein_dataset,
)


def test_package_root_legacy_groove_exports_are_preserved():
    assert extract_groove is decompose_domains
    assert parse_class_i is decompose_class_i
    assert parse_class_ii_alpha is decompose_class_ii_alpha
    assert parse_class_ii_beta is decompose_class_ii_beta


def test_package_root_exports_versioned_mhc_protein_data():
    assert available_mhc_protein_dataset_versions() == ("uniprot-2026_03-r1",)
    assert callable(validate_mhc_protein_dataset)
