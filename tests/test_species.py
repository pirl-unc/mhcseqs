from mhcseqs.species import (
    MHC_SPECIES_CATEGORIES,
    get_canonical_prefix,
    get_latin_name,
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


def test_mhc_species_other_mammal():
    assert normalize_mhc_species("Bos taurus") == "other_mammal"
    assert normalize_mhc_species("Sus scrofa") == "other_mammal"


def test_latin_name():
    assert get_latin_name("human") == "Homo sapiens"
    assert get_latin_name("Homo sapiens") == "Homo sapiens"
    assert get_latin_name("mouse") == "Mus musculus"


def test_canonical_prefix():
    assert get_canonical_prefix("human") == "HLA"
    assert get_canonical_prefix("Macaca mulatta") == "Mamu"
    assert get_canonical_prefix("mouse") == "H2"
    assert get_canonical_prefix("Bos taurus") == "BoLA"
