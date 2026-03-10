from mhcseqs.alleles import (
    allele_suffix_flags,
    normalize_mhc_class,
)


def test_normalize_mhc_class():
    assert normalize_mhc_class("I") == "I"
    assert normalize_mhc_class("II") == "II"
    assert normalize_mhc_class("class I") == "I"
    assert normalize_mhc_class("class II") == "II"
    assert normalize_mhc_class(None) is None


def test_normalize_mhc_class_default():
    assert normalize_mhc_class(None, default="I") == "I"
    assert normalize_mhc_class("", default="II") == "II"


def test_allele_suffix_flags_null():
    flags = allele_suffix_flags("HLA-A*02:01N")
    assert flags["is_null"] is True
    assert flags["is_questionable"] is False
    assert flags["is_pseudogene"] is False


def test_allele_suffix_flags_questionable():
    flags = allele_suffix_flags("HLA-A*02:01Q")
    assert flags["is_questionable"] is True


def test_allele_suffix_flags_normal():
    flags = allele_suffix_flags("HLA-A*02:01")
    assert flags["is_null"] is False
    assert flags["is_questionable"] is False
    assert flags["is_pseudogene"] is False
