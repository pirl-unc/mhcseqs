from mhcseqs.alleles import (
    allele_suffix_flags,
    infer_gene,
    infer_mhc_class,
    infer_species,
    normalize_allele_name,
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


def test_allele_suffix_flags_rano_haplotype_not_null():
    """Rano-A1*n is a haplotype name, not a null allele."""
    flags = allele_suffix_flags("Rano-A1*n")
    assert flags["is_null"] is False
    assert flags["is_questionable"] is False
    assert flags["is_pseudogene"] is False


def test_allele_suffix_flags_rano_haplotype_q_not_questionable():
    """Rano-A1*q is a haplotype name, not a questionable allele."""
    flags = allele_suffix_flags("Rano-A*q")
    assert flags["is_questionable"] is False


def test_allele_suffix_flags_h2_haplotype_not_null():
    """H2-D*b is a haplotype name, not a suffix-bearing allele."""
    flags = allele_suffix_flags("H2-D*b")
    assert flags["is_null"] is False
    assert flags["is_questionable"] is False


def test_normalize_allele_name_full():
    assert normalize_allele_name("HLA-A*02:01") == "HLA-A*02:01"


def test_normalize_allele_name_compact():
    result = normalize_allele_name("A0201")
    assert "02" in result
    assert "01" in result


def test_infer_gene_hla_a():
    result = infer_gene("HLA-A*02:01")
    assert result == "A"


def test_infer_gene_drb1():
    result = infer_gene("HLA-DRB1*01:01")
    assert result == "DRB1"


def test_infer_mhc_class_i():
    assert infer_mhc_class("HLA-A*02:01") == "I"
    assert infer_mhc_class("HLA-B*07:02") == "I"
    assert infer_mhc_class("HLA-C*04:01") == "I"


def test_infer_mhc_class_ii():
    assert infer_mhc_class("HLA-DRB1*01:01") == "II"
    assert infer_mhc_class("HLA-DQA1*01:01") == "II"


def test_infer_mhc_class_none():
    assert infer_mhc_class(None) is None
    assert infer_mhc_class("") is None


def test_infer_species_human():
    result = infer_species("HLA-A*02:01")
    assert result == "human"


def test_infer_species_returns_7_class():
    from mhcseqs.species import MHC_SPECIES_CATEGORIES

    result = infer_species("HLA-A*02:01")
    assert result in MHC_SPECIES_CATEGORIES


def test_normalize_mhc_class_aliases():
    assert normalize_mhc_class("MHC-I") == "I"
    assert normalize_mhc_class("MHC-II") == "II"
    assert normalize_mhc_class("CLASSI") == "I"
    assert normalize_mhc_class("CLASSII") == "II"
    assert normalize_mhc_class("IA") == "I"
    assert normalize_mhc_class("IIA") == "II"
