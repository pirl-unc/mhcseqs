import mhcgnomes
import pytest

from mhcseqs.alleles import (
    allele_suffix_flags,
    build_accessioned_allele,
    canonical_allele_name,
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


def test_allele_suffix_flags_gene_level_pseudogenes():
    """mhcgnomes marks gene-level pseudogenes even without a Ps suffix (#95)."""
    for name in ("HLA-H*01:01:01:01", "HLA-DRB9*01:01", "HomoSapiens-H*01:01:01:01"):
        flags = allele_suffix_flags(name)
        assert flags == {"is_null": False, "is_questionable": False, "is_pseudogene": True}

    j_flags = allele_suffix_flags("HLA-J*01:01:01:02N")
    assert j_flags == {"is_null": True, "is_questionable": False, "is_pseudogene": True}

    ps_flags = allele_suffix_flags("Patr-A*05:01Ps")
    assert ps_flags == {"is_null": False, "is_questionable": False, "is_pseudogene": True}


def test_allele_suffix_flags_accepts_parsed_mhcgnomes_object():
    parsed = parse_allele_name("HLA-H*01:01:01:01")
    assert allele_suffix_flags(parsed) == {"is_null": False, "is_questionable": False, "is_pseudogene": True}


def test_allele_suffix_flags_unparseable_returns_false():
    assert allele_suffix_flags("B2M_human") == {"is_null": False, "is_questionable": False, "is_pseudogene": False}
    assert allele_suffix_flags(None) == {"is_null": False, "is_questionable": False, "is_pseudogene": False}


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
    assert parse_allele_name("H2-Kb").species.name == "Mus musculus"
    assert parse_allele_name("H-2Kb").species.name == "Mus musculus"
    assert parse_allele_name("H_2-Kb").species.name == "Mus musculus"
    assert parse_allele_name("HLA_A*02:01").species.name == "Homo sapiens"
    assert parse_allele_name("RhLA-A1*001:01").species.name == "Macaca sp."
    assert parse_allele_name("GPLA-DRA").species.name == "Cavia porcellus"
    assert parse_allele_name("XL-A-DAB").species.name == "Xenopus laevis"
    chicken_bf2 = parse_allele_name("B-F2", species="Gallus gallus")
    assert chicken_bf2.species.name == "Gallus gallus"
    assert chicken_bf2.name == "BF2"


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


def test_is_non_mhc_gene_unwraps_synthetic_ref_and_predicted_tags():
    """Regression test for #135: ~ref:Species|Acc:Gene tags and RFX symbols must be caught."""
    assert is_non_mhc_gene("~ref:OctodonDegus|A0A6P6E5A6:Ciita", species="Octodon degus")
    assert is_non_mhc_gene("~ref:ChinchillaLanigera|A0A8C2VVP9:CIITA", species="Chinchilla lanigera")
    assert is_non_mhc_gene("~ref:CastorCanadensis|A0A250XZQ4:Hm13", species="Castor canadensis")
    assert is_non_mhc_gene("~ref:DipodomysOrdii|A0A1S3GLH2:Rfx1", species="Dipodomys ordii")
    assert is_non_mhc_gene("~ref:RhinolophusFerrumequinum|A0A671DWX3:PRR3", species="Rhinolophus ferrumequinum")
    assert is_non_mhc_gene("~ref:MusMusculus|P04186:Cfb", species="Mus musculus")
    assert is_non_mhc_gene("~ref:MusMusculus|P50171:Hsd17b8", species="Mus musculus")
    assert is_non_mhc_gene("~ref:RattusNorvegicus|A6KHQ1:H13_predicted", species="Rattus norvegicus")
    assert is_non_mhc_gene("~ref:DipodomysOrdii|A0A1S3FNJ3:Hmha1", species="Dipodomys ordii")
    assert is_non_mhc_gene("HMHA1")
    assert is_non_mhc_gene("RFX5")
    assert not is_non_mhc_gene("~ref:AnasPlatyrhynchos|A0A1B0P060:UAA", species="Anas platyrhynchos")
    assert not is_non_mhc_gene("HLA-A", species="Homo sapiens")


def test_parse_gene_class_enforces_iia_alpha_chain():
    """Regression test for #137 / mhcgnomes#286: Class IIa genes like DAA-1 must have chain='alpha'."""
    result = parse_gene_class("SphenodonPunctatus-DAA-1", species="Sphenodon punctatus")
    assert result is not None
    assert result["mhc_class"] == "II"
    assert result["chain"] == "alpha"


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


def test_short_hla_serotypes_are_not_coerced_to_allele_groups():
    """Serotypes such as A2, HLA-A2, B15, B27 must not become allele groups (#95)."""
    for name in ("A2", "HLA-A2", "B7", "B15", "B27"):
        assert parse_allele_name(name) is None
        serotype = parse_allele_name(name, allow_serotypes=True)
        assert isinstance(serotype, mhcgnomes.Serotype)
        assert serotype.species.name == "Homo sapiens"


def test_infer_gene_returns_none_without_string_fallback():
    """Unparseable strings return None rather than splitting on '-'/'*' (#95)."""
    assert infer_gene("NotASpecies-FakeGene*99:99") is None
    assert infer_gene("") is None
    assert infer_gene(None) is None


def test_parse_allele_name_propagates_unexpected_errors(monkeypatch):
    """Unexpected mhcgnomes bugs must not be silently swallowed (#97)."""

    def boom(*_args, **_kwargs):
        raise RuntimeError("simulated bug")

    monkeypatch.setattr(mhcgnomes, "parse", boom)
    with pytest.raises(RuntimeError, match="simulated bug"):
        parse_allele_name("HLA-A*02:01")


def test_parse_gene_class_propagates_unexpected_errors(monkeypatch):
    def boom(*_args, **_kwargs):
        raise RuntimeError("simulated bug")

    monkeypatch.setattr(mhcgnomes, "parse_gene_class", boom)
    with pytest.raises(RuntimeError, match="simulated bug"):
        parse_gene_class("HLA-A*02:01")


def test_canonical_allele_name_uses_full_binomial_species_alias():
    parsed = parse_allele_name("HLA-A*02:01:01:02L")
    assert canonical_allele_name(parsed) == "HomoSapiens-A*02:01:01:02L"
    assert canonical_allele_name(parsed, allele_fields=2) == "HomoSapiens-A*02:01L"
    assert canonical_allele_name(parse_allele_name("HLA-DRB1"), allele_fields=2) == "HomoSapiens-DRB1"


def test_genus_level_species_emit_roundtrip_mhc_prefix():
    """Issues #92 and #130: Sus sp./Bos sp./Canis sp./Felis sp. emit SLA/BoLA/DLA/FLA and round-trip."""
    cases = [
        ("SLA-1*01:01:01", "SLA-1*01:01"),
        ("BoLA-1*007:01:01:01", "BoLA-1*007:01"),
        ("DLA-12*001:01:02", "DLA-12*001:01"),
        ("FLA-E*003:01", "FLA-E*003:01"),
    ]
    for raw, expected in cases:
        norm = normalize_allele_name(raw)
        assert norm == expected
        assert parse_allele_name(norm) is not None
        assert normalize_allele_name(norm) == expected


def test_canonical_species_prefix_roundtrips_for_every_mhcgnomes_species():
    """Every species in mhcgnomes must produce a canonical prefix that round-trips cleanly (#92, #130)."""
    from mhcseqs.alleles import canonical_species_prefix

    seen: set[str] = set()
    for sp in mhcgnomes.species.species_name_to_species_object.values():
        if sp.name in seen:
            continue
        seen.add(sp.name)
        prefix = canonical_species_prefix(sp)
        assert prefix, f"Empty canonical_species_prefix for {sp.name!r}"
        rt = mhcgnomes.Species.get(prefix)
        assert rt is not None and rt.name == sp.name, f"Prefix {prefix!r} for {sp.name!r} round-tripped to {rt!r}"
        fake_allele = f"{prefix}-A*01:01"
        parsed = parse_allele_name(fake_allele)
        if parsed is not None:
            norm = canonical_allele_name(parsed, allele_fields=2)
            assert parse_allele_name(norm) is not None, f"Canonical allele {norm!r} failed round-trip"
            assert normalize_allele_name(norm) == norm, f"normalize_allele_name not idempotent for {norm!r}"


def test_accessioned_alleles_parse_normalize_and_roundtrip():
    """Issue #150 / mhcgnomes#309 / mhcgnomes#326: AccessionedAllele parsing, canonicalization, and metadata."""
    cases = [
        ("Saha-UA*uniprot:G3V1L5", "SarcophilusHarrisii-UA*uniprot:G3V1L5", "UA", "I", "other_mammal"),
        ("TriakisScyllium-UAA*uniprot:Q863D3", "TriakisScyllium-UAA*uniprot:Q863D3", "UAA", "I", "fish"),
        ("Saha-uniprot:G3V1L5", "SarcophilusHarrisii-uniprot:G3V1L5", None, None, "other_mammal"),
        ("Saha-classI*uniprot:G3V1L5", "SarcophilusHarrisii-classI*uniprot:G3V1L5", None, "I", "other_mammal"),
        ("Saha*uniprot:G3V1L5 class I", "SarcophilusHarrisii-classI*uniprot:G3V1L5", None, "I", "other_mammal"),
        ("Saha*uniprot:G3V1L5 class II", "SarcophilusHarrisii-classII*uniprot:G3V1L5", None, "II", "other_mammal"),
        ("BoLA*uniprot:P12345 class II", "BoLA-classII*uniprot:P12345", None, "II", "ungulate"),
    ]
    for raw, expected, expected_gene, expected_class, expected_species_bucket in cases:
        parsed = parse_allele_name(raw)
        assert isinstance(parsed, mhcgnomes.AccessionedAllele)
        assert normalize_allele_name(raw) == expected
        assert " " not in expected
        assert normalize_allele_name(expected) == expected
        assert infer_gene(expected) == expected_gene
        assert infer_mhc_class(raw) == expected_class
        if expected_gene is not None:
            assert infer_mhc_class(expected) == expected_class
        assert infer_species(expected) == expected_species_bucket
        assert allele_suffix_flags(parsed) == {"is_null": False, "is_questionable": False, "is_pseudogene": False}

    # Bare UniProt accessions and mechanically truncated derived prefixes are rejected.
    for bare in ("G3V1L5", "H2LR16", "P01898", "TriaScyl-uniprot:Q863D3", "TriaScyl*uniprot:Q863D3 class I"):
        assert parse_allele_name(bare) is None


def test_build_accessioned_allele_handles_gene_and_placeholder_fallbacks():
    with_gene = build_accessioned_allele(
        species="Triakis scyllium",
        accession_source="uniprot",
        accession="Q863D3",
        gene="TriakisScyllium-UAA",
        mhc_class="I",
    )
    assert isinstance(with_gene, mhcgnomes.AccessionedAllele)
    assert canonical_allele_name(with_gene) == "TriakisScyllium-UAA*uniprot:Q863D3"
    assert with_gene.is_class1

    from_allele_label = build_accessioned_allele(
        species="Xenopus laevis",
        accession_source="uniprot",
        accession="Q9TPA0",
        gene="XenopusLaevis-UAA*l",
        mhc_class="I",
    )
    assert isinstance(from_allele_label, mhcgnomes.AccessionedAllele)
    assert canonical_allele_name(from_allele_label) == "XenopusLaevis-UAA*uniprot:Q9TPA0"

    placeholder = build_accessioned_allele(
        species="Sarcophilus harrisii",
        accession_source="uniprot",
        accession="G3V1L5",
        gene="~ref:SarcophilusHarrisii|G3V1L5:MHC-I",
        mhc_class="I",
    )
    assert isinstance(placeholder, mhcgnomes.AccessionedAllele)
    assert canonical_allele_name(placeholder) == "SarcophilusHarrisii-classI*uniprot:G3V1L5"
    assert placeholder.is_class1

    empty_gene = build_accessioned_allele(
        species="Cricetulus griseus",
        accession_source="uniprot",
        accession="A0A061HWC8",
        gene="",
        mhc_class="I",
    )
    assert isinstance(empty_gene, mhcgnomes.AccessionedAllele)
    assert canonical_allele_name(empty_gene) == "CricetulusGriseus-classI*uniprot:A0A061HWC8"
    assert empty_gene.is_class1
