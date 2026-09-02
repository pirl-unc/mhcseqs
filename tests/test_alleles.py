from mhcseqs.alleles import (
    _coerce_allele_name,
    allele_suffix_flags,
    infer_gene,
    infer_mhc_class,
    infer_species,
    is_non_mhc_gene,
    normalize_allele_name,
    normalize_mhc_class,
    parse_allele_name,
    parse_gene_class,
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
    assert normalize_allele_name("HLA-A*02:01") == "HomoSapiens-A*02:01"


def test_normalize_allele_name_compact():
    result = normalize_allele_name("A0201")
    assert result == "HomoSapiens-A*02:01"


def test_bare_hla_allele_uses_mhcgnomes_default_species():
    parsed = parse_allele_name("A*02:01")
    assert parsed.to_string() == "HLA-A*02:01"
    assert parsed.species.name == "Homo sapiens"
    assert parsed.species_source == "default"


def test_bare_hla_allele_can_require_explicit_species():
    assert parse_allele_name("A*02:01", require_explicit_species=True) is None


def test_normalize_bare_hla_allele():
    assert normalize_allele_name("A*02:01") == "HomoSapiens-A*02:01"


def test_generated_short_alias_is_not_accepted_without_external_evidence():
    assert parse_allele_name("CyanCaer-DAB1") is None
    assert parse_gene_class("CyanCaer-DAB1") is None
    assert normalize_allele_name("CyanistesCaeruleus-DAB1") == "CyanistesCaeruleus-DAB1"


def test_source_attested_long_tail_prefix_normalizes_to_full_binomial():
    parsed = parse_allele_name("Acsi-UA", species="Acipenser sinensis")
    assert parsed.species.name == "Acipenser sinensis"
    assert parsed.gene.name == "UA"


def test_source_attested_concatenated_prefix_is_supported():
    parsed = parse_allele_name("XimuDXB", species="Xiphophorus multilineatus")
    assert parsed.species.name == "Xiphophorus multilineatus"
    assert parsed.gene.name == "DXB"


def test_historical_nomenclature_prefixes_are_supported():
    assert parse_allele_name("HL-A-A*02:01").species.name == "Homo sapiens"
    assert parse_allele_name("H-2-Kb").species.name == "Mus musculus"
    assert parse_allele_name("RhLA-A1*001:01").species.name == "Macaca mulatta"
    assert parse_allele_name("GPLA-DRA").species.name == "Cavia porcellus"
    assert parse_allele_name("XL-A-DAB").species.name == "Xenopus laevis"
    assert parse_allele_name("B-BF2", species="Gallus gallus").species.name == "Gallus gallus"


def test_historical_hyphenated_h2_survives_ambiguous_modern_coercion():
    assert parse_allele_name("H2-Kb") is None
    assert parse_allele_name("H-2-Kb").species.name == "Mus musculus"


def test_colliding_attested_prefix_requires_species_context():
    assert parse_allele_name("Gogo-DAB") is None
    parsed = parse_allele_name("Gogo-DAB", species="Gobio gobio")
    assert parsed.species.name == "Gobio gobio"

    assert parse_allele_name("RAJA-UA") is None
    assert parse_allele_name("RAJA-UA", species="Rana japonica").species.name == "Rana japonica"


def test_parse_allele_name_preserves_strict_species_constraint():
    assert parse_allele_name("HLA-A*02:01", species="Bos taurus") is None
    assert parse_allele_name("BoLA-DRB3*01:01", species="Homo sapiens") is None


def test_parse_allele_name_rejects_non_molecule_results():
    assert parse_allele_name("MHC class II") is None
    assert parse_allele_name("HLA-DR15") is None
    assert parse_allele_name("H2-b") is None


def test_parse_allele_name_exposes_and_can_require_species_provenance():
    explicit = parse_allele_name("Gaga-BLB2*02")
    inferred = parse_allele_name("BLB2*02")
    assert explicit.species_source == "explicit"
    assert inferred.species_source == "inferred"
    assert parse_allele_name("BLB2*02", require_explicit_species=True) is None


def test_parse_gene_class_is_species_aware():
    assert parse_gene_class("BLB") is None
    result = parse_gene_class("BLB", species="Falco peregrinus")
    assert result == {
        "mhc_class": "II",
        "chain": "beta",
        "non_mhc": False,
        "source": "heuristic_suffix",
    }


def test_local_non_mhc_names_override_allele_like_mhcgnomes_parse():
    assert is_non_mhc_gene("Kdm5d", species="Mus musculus")
    assert is_non_mhc_gene("Daxx", species="Mus musculus")


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


def test_normalize_mhc_class_mhcgnomes_subclasses():
    """mhcgnomes returns 'Ia', 'Ib', 'IIa' — subclass designations, not mouse I-A."""
    assert normalize_mhc_class("Ia") == "I"
    assert normalize_mhc_class("Ib") == "I"
    assert normalize_mhc_class("IIa") == "II"
    assert normalize_mhc_class("IIb") == "II"


def test_coerce_allele_name_hla_short_forms():
    """Known HLA genes are coerced to allele notation."""
    assert _coerce_allele_name("A2") == "HLA-A*02"
    assert _coerce_allele_name("B7") == "HLA-B*07"
    assert _coerce_allele_name("C4") == "HLA-C*04"
    assert _coerce_allele_name("E1") == "HLA-E*01"
    assert _coerce_allele_name("G1") == "HLA-G*01"
    assert _coerce_allele_name("HLA-A2") == "HLA-A*02"
    assert _coerce_allele_name("A0201") == "HLA-A*0201"


def test_coerce_allele_name_rejects_non_hla():
    """Non-HLA gene-like tokens must NOT be coerced to HLA format."""
    assert _coerce_allele_name("RT1") == "RT1"
    assert _coerce_allele_name("SLA1") == "SLA1"
    assert _coerce_allele_name("K2") == "K2"
    assert _coerce_allele_name("D3") == "D3"
    assert _coerce_allele_name("DRB1") == "DRB1"


def test_coerce_allele_name_h2_passthrough():
    """Mouse H-2 alleles pass through correctly."""
    assert _coerce_allele_name("H2-Kb") == "H2-Kb"
    assert _coerce_allele_name("H-2-Kb") == "H2-Kb"
    assert _coerce_allele_name("H-2Kb") == "H2-Kb"
