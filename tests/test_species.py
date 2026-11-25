from mhcseqs.species import (
    MHC_SPECIES_CATEGORIES,
    extract_latin_binomial,
    full_species_name_alias,
    get_canonical_prefix,
    get_latin_name,
    get_ncbi_taxid,
    get_species_clade,
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


def test_species_clade_and_ncbi_taxid_from_mhcgnomes():
    assert get_species_clade("Homo sapiens") == "human"
    assert get_ncbi_taxid("Homo sapiens") == 9606
    assert get_species_clade("Canis lupus familiaris") == "carnivore"
    assert get_ncbi_taxid("Canis lupus familiaris") == 9615
    assert get_species_clade("Canis lupus") == "carnivore"
    assert get_ncbi_taxid("Canis lupus") == 9612
    assert get_species_clade("Polyodon spathula (Squalus spathula)") == "bony_fish"
    assert get_species_clade("Squalus acanthias") == "shark_ray"
    assert get_species_clade("Pteropus alecto") == "bat"
    assert get_species_clade("Trichosurus vulpecula") == "marsupial"
    assert get_species_clade("Ornithorhynchus anatinus") == "monotreme"
    assert get_species_clade("Crocodylus porosus (Saltwater crocodile)") == "reptile"
    assert get_species_clade("Triturus cristatus") == "amphibian"
    assert get_ncbi_taxid("Triturus cristatus") == 8323
    assert get_species_clade("Bubulcus ibis (Cattle egret) (Ardea ibis)") == "bird"
    assert get_species_clade("Struthio camelus (Common ostrich)") == "bird"
    assert get_species_clade(None) is None
    assert get_ncbi_taxid(None) is None


def test_normalize_mhc_species_avoids_common_name_substring_false_positives():
    """Authoritative mhcgnomes clade prevents regex traps like Cattle egret -> ungulate."""
    assert normalize_mhc_species("Bubulcus ibis (Cattle egret) (Ardea ibis)") == "bird"
    assert normalize_mhc_species("Struthio camelus (Common ostrich)") == "bird"
    assert normalize_mhc_species("Polyodon spathula (Squalus spathula)") == "fish"
    assert normalize_mhc_species("Acipenser sinensis (Chinese sturgeon)") == "fish"
    assert normalize_mhc_species("Triturus cristatus") == "other_vertebrate"


def test_build_imgt_profile_classify_clade_uses_mhcgnomes():
    from scripts.build_imgt_profile import classify_clade

    assert classify_clade({"species": "Polyodon spathula (Squalus spathula)", "species_category": "fish"}) == "bony_fish"
    assert classify_clade({"species": "Squalus acanthias", "species_category": "fish"}) == "shark_ray"
    assert classify_clade({"species": "Triturus cristatus", "species_category": "other_vertebrate"}) == "amphibian"
    assert classify_clade({"species": "Ornithorhynchus anatinus", "species_category": "other_mammal"}) == "monotreme"
    assert classify_clade({"species": "Pteropus alecto", "species_category": "other_mammal"}) == "bat"
    assert classify_clade({"species": "Bubulcus ibis (Cattle egret)", "species_category": "bird"}) == "bird"
