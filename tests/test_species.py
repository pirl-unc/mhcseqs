from mhcseqs.species import (
    ESTABLISHED_MHC_PREFIX_SOURCES,
    ESTABLISHED_MHC_PREFIXES,
    MHC_SPECIES_CATEGORIES,
    extract_latin_binomial,
    find_mhc_prefix_aliases,
    full_species_name_alias,
    get_canonical_prefix,
    get_established_mhc_prefix,
    get_established_mhc_prefixes,
    get_latin_name,
    get_mhc_prefix_registry,
    normalize_mhc_species,
    normalize_species,
)


def test_categories_include_fish():
    assert "fish" in MHC_SPECIES_CATEGORIES


def test_normalize_species_human():
    assert normalize_species("Homo sapiens") == "human"
    assert normalize_species("human") == "human"


def test_normalize_species_macaque():
    assert normalize_species("Macaca mulatta") == "macaque"
    assert normalize_species("rhesus") == "macaque"


def test_normalize_species_mouse():
    assert normalize_species("Mus musculus") == "mouse"
    assert normalize_species("mouse") == "mouse"


def test_normalize_species_salmon():
    assert normalize_species("Salmo salar") == "salmon"


def test_normalize_species_none():
    assert normalize_species(None) is None
    assert normalize_species("") is None


def test_mhc_species_human():
    assert normalize_mhc_species("Homo sapiens") == "human"


def test_mhc_species_nhp():
    assert normalize_mhc_species("Macaca mulatta") == "nhp"
    assert normalize_mhc_species("Pan troglodytes") == "nhp"


def test_mhc_species_murine():
    assert normalize_mhc_species("Mus musculus") == "murine"
    assert normalize_mhc_species("Rattus") == "murine"


def test_mhc_species_fish():
    assert normalize_mhc_species("Salmo salar") == "fish"
    assert normalize_mhc_species("Danio") == "fish"


def test_mhc_species_bird():
    assert normalize_mhc_species("Gallus") == "bird"


def test_mhc_species_ungulate():
    assert normalize_mhc_species("Bos taurus") == "ungulate"
    assert normalize_mhc_species("Sus scrofa") == "ungulate"
    assert normalize_mhc_species("Equus caballus") == "ungulate"


def test_mhc_species_carnivore():
    assert normalize_mhc_species("Canis lupus familiaris") == "carnivore"
    assert normalize_mhc_species("Felis catus") == "carnivore"


def test_mhc_species_cetacean_merged_to_other_mammal():
    assert normalize_mhc_species("Tursiops truncatus") == "other_mammal"
    assert normalize_mhc_species("Balaenoptera musculus") == "other_mammal"


def test_alouatta_is_nhp_not_pig():
    assert normalize_species("Alouatta pigra") == "other_nhp"
    assert normalize_mhc_species("Alouatta pigra") == "nhp"


def test_latin_name():
    assert get_latin_name("human") == "Homo sapiens"
    assert get_latin_name("Homo sapiens") == "Homo sapiens"
    assert get_latin_name("mouse") == "Mus musculus"


def test_extract_latin_binomial():
    croc = "Crocodylus porosus (Saltwater crocodile) (Estuarine crocodile)"
    assert extract_latin_binomial(croc) == "Crocodylus porosus"
    assert extract_latin_binomial("Homo sapiens (Human)") == "Homo sapiens"
    assert extract_latin_binomial("Mus musculus") == "Mus musculus"
    assert extract_latin_binomial("Gallus gallus (Chicken)") == "Gallus gallus"
    assert extract_latin_binomial("Cyanistes caeruleus caeruleus") == "Cyanistes caeruleus"
    assert extract_latin_binomial("") == ""
    assert extract_latin_binomial(None) == ""


def test_canonical_prefix():
    assert get_canonical_prefix("human") == "HomoSapiens"
    assert get_canonical_prefix("Macaca mulatta") == "MacacaMulatta"
    assert get_canonical_prefix("mouse") == "MusMusculus"
    assert get_canonical_prefix("Bos taurus") == "BosTaurus"


def test_full_species_name_alias_does_not_invent_short_codes():
    assert full_species_name_alias("Cyanistes caeruleus") == "CyanistesCaeruleus"
    assert full_species_name_alias("Macaca fascicularis") == "MacacaFascicularis"
    assert full_species_name_alias("Bos indicus") == "BosIndicus"
    assert full_species_name_alias("Coregonus sp.") == "CoregonusSp"
    assert full_species_name_alias("Notamacropus") == "Notamacropus"


def test_established_prefixes_are_separate_source_aliases():
    assert get_established_mhc_prefix("human") == "HLA"
    assert get_established_mhc_prefix("Macaca mulatta") == "Mamu"
    assert get_established_mhc_prefix("Macaca fascicularis") == "Mafa"
    assert get_established_mhc_prefix("Macaca nemestrina") == "Mane"
    assert get_established_mhc_prefix("Bos taurus") == "BoLA"
    assert get_established_mhc_prefix("Cyanistes caeruleus") == "Cyca"


def test_current_and_historical_prefixes_have_explicit_provenance():
    assert get_established_mhc_prefixes("Homo sapiens")[:2] == ("HLA", "HL-A")
    assert {"Mamu", "RhLA"}.issubset(get_established_mhc_prefixes("Macaca mulatta"))
    assert {"H2", "H-2"}.issubset(get_established_mhc_prefixes("Mus musculus"))
    assert {"RT1", "Ag-B", "H-1"}.issubset(get_established_mhc_prefixes("Rattus norvegicus"))
    assert {"Gaga", "B"}.issubset(get_established_mhc_prefixes("Gallus gallus"))
    assert {"GPLA", "GPL-A"}.issubset(get_established_mhc_prefixes("Cavia porcellus"))
    assert {"XLA", "XL-A"}.issubset(get_established_mhc_prefixes("Xenopus laevis"))
    assert "Livea" in get_established_mhc_prefixes("Litoria verreauxii")
    assert "Frcoe" in get_established_mhc_prefixes("Fringilla coelebs")


def test_registry_covers_ipd_and_external_database_prefixes():
    registry = get_mhc_prefix_registry()
    assert sum(entry.status == "ipd_current" for entry in registry) == 137
    assert sum(entry.status == "external_database" for entry in registry) == 298
    assert all(entry.evidence.startswith(("https://", "http://")) for entry in registry)

    for entry in registry:
        prefix, body, matches = find_mhc_prefix_aliases(
            f"{entry.prefix}-DRA",
            species=entry.species,
        )
        assert prefix
        assert body == "DRA"
        assert entry in matches


def test_colliding_published_prefixes_require_species_context():
    _, _, gogo = find_mhc_prefix_aliases("Gogo-DAB")
    assert {entry.species for entry in gogo} == {"Gobio gobio", "Gorilla gorilla"}

    _, _, fish = find_mhc_prefix_aliases("Gogo-DAB", species="Gobio gobio")
    assert {entry.species for entry in fish} == {"Gobio gobio"}


def test_every_established_short_prefix_has_external_evidence():
    assert set(ESTABLISHED_MHC_PREFIXES.values()) == set(ESTABLISHED_MHC_PREFIX_SOURCES)
    assert all(url.startswith("https://") for url in ESTABLISHED_MHC_PREFIX_SOURCES.values())
