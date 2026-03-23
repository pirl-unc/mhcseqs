from mhcseqs import (
    decompose_class_i,
    decompose_class_ii_alpha,
    decompose_class_ii_beta,
    decompose_domains,
    extract_groove,
    parse_class_i,
    parse_class_ii_alpha,
    parse_class_ii_beta,
)


def test_package_root_legacy_groove_exports_are_preserved():
    assert extract_groove is decompose_domains
    assert parse_class_i is decompose_class_i
    assert parse_class_ii_alpha is decompose_class_ii_alpha
    assert parse_class_ii_beta is decompose_class_ii_beta
