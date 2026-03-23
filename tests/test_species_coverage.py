"""Tests that expected species are present in the diverse MHC dataset."""

import csv
from pathlib import Path

import pytest

DIVERSE_CSV = Path(__file__).resolve().parent.parent / "mhcseqs" / "diverse_mhc_sequences.csv"
B2M_CSV = Path(__file__).resolve().parent.parent / "mhcseqs" / "b2m_sequences.csv"


def _diverse_organisms():
    """Load all unique organisms from the diverse CSV."""
    organisms = set()
    with open(DIVERSE_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            org = row.get("organism", "")
            if org:
                organisms.add(org.split("(")[0].strip())
    return organisms


def _diverse_species_genes():
    """Load (organism, gene) pairs from the diverse CSV."""
    pairs = set()
    with open(DIVERSE_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            org = row.get("organism", "").split("(")[0].strip()
            gene = row.get("gene", "")
            if org and gene:
                pairs.add((org, gene))
    return pairs


@pytest.fixture(scope="module")
def organisms():
    return _diverse_organisms()


@pytest.fixture(scope="module")
def species_genes():
    return _diverse_species_genes()


# --- Mammals ---


class TestMammalCoverage:
    def test_rabbit(self, organisms):
        assert "Oryctolagus cuniculus" in organisms

    def test_mouse(self, organisms):
        assert "Mus musculus" in organisms

    def test_rat(self, organisms):
        assert "Rattus norvegicus" in organisms

    def test_guinea_pig(self, organisms):
        assert "Cavia porcellus" in organisms

    def test_golden_hamster(self, organisms):
        assert "Mesocricetus auratus" in organisms

    def test_chinese_hamster(self, organisms):
        assert "Cricetulus griseus" in organisms

    def test_naked_mole_rat(self, organisms):
        assert "Heterocephalus glaber" in organisms

    def test_opossum(self, organisms):
        assert "Monodelphis domestica" in organisms

    def test_koala(self, organisms):
        assert "Phascolarctos cinereus" in organisms

    def test_tasmanian_devil(self, organisms):
        assert "Sarcophilus harrisii" in organisms

    def test_platypus(self, organisms):
        assert "Ornithorhynchus anatinus" in organisms

    def test_vampire_bat(self, organisms):
        assert "Desmodus rotundus" in organisms

    def test_beaver(self, organisms):
        assert "Castor canadensis" in organisms

    def test_bank_vole(self, organisms):
        assert "Myodes glareolus" in organisms

    def test_deer_mouse(self, organisms):
        assert "Peromyscus maniculatus" in organisms

    def test_prairie_vole(self, organisms):
        assert "Microtus ochrogaster" in organisms


# --- Birds ---


class TestBirdCoverage:
    def test_chicken(self, organisms):
        assert "Gallus gallus" in organisms

    def test_ostrich(self, organisms):
        assert "Struthio camelus" in organisms

    def test_emu(self, organisms):
        assert "Dromaius novaehollandiae" in organisms

    def test_kiwi(self, organisms):
        assert "Apteryx owenii" in organisms

    def test_barn_owl(self, organisms):
        assert "Tyto alba" in organisms

    def test_quail(self, organisms):
        assert "Coturnix japonica" in organisms


# --- Reptiles ---


class TestReptileCoverage:
    def test_saltwater_croc(self, organisms):
        assert "Crocodylus porosus" in organisms

    def test_green_sea_turtle(self, organisms):
        assert "Chelonia mydas" in organisms

    def test_tuatara(self, organisms):
        assert "Sphenodon punctatus" in organisms

    def test_marine_iguana(self, organisms):
        assert "Amblyrhynchus cristatus" in organisms


# --- Amphibians ---


class TestAmphibianCoverage:
    def test_axolotl(self, organisms):
        assert "Ambystoma mexicanum" in organisms

    def test_xenopus_tropicalis(self, organisms):
        assert "Xenopus tropicalis" in organisms

    def test_xenopus_laevis(self, organisms):
        assert "Xenopus laevis" in organisms


# --- Fish ---


class TestFishCoverage:
    def test_zebrafish(self, organisms):
        assert "Danio rerio" in organisms

    def test_salmon(self, organisms):
        assert any("Salmo" in o or "Oncorhynchus" in o for o in organisms)

    def test_tilapia(self, organisms):
        assert "Oreochromis niloticus" in organisms

    def test_stickleback(self, organisms):
        assert "Gasterosteus aculeatus" in organisms

    def test_nurse_shark(self, organisms):
        assert "Ginglymostoma cirratum" in organisms


# --- Gene coverage for key species ---


class TestRabbitGenes:
    def test_rabbit_class_i(self, species_genes):
        assert any(g.startswith("RLA-A") for o, g in species_genes if o == "Oryctolagus cuniculus")

    def test_rabbit_class_ii(self, species_genes):
        rabbit_genes = {g for o, g in species_genes if o == "Oryctolagus cuniculus"}
        assert any("DRB" in g for g in rabbit_genes)


class TestMouseGenes:
    def test_mouse_present(self, species_genes):
        mouse_genes = {g for o, g in species_genes if "Mus musculus" in o}
        assert len(mouse_genes) > 0


class TestCrocGenes:
    def test_croc_ua(self, species_genes):
        assert ("Crocodylus porosus", "Crpo-UA") in species_genes

    def test_croc_class_ii(self, species_genes):
        croc_genes = {g for o, g in species_genes if o == "Crocodylus porosus"}
        assert any("DAB" in g or "DB0" in g for g in croc_genes)
