"""
Tests for mhcgnomes species= param enforcement.

The species= param should be authoritative: the parsed result must have
exactly the requested species. If the gene name embeds a different species
(e.g., HLA-DRB8 implies human), passing species="Kareius bicoloratus"
should raise an error or return None, never silently return Homo sapiens.

Generated from mhcseqs diverse dataset analysis (2026-03-17).
All 61 cases below return the wrong species with mhcgnomes 3.9.0.

To run:
    pytest data/mhcgnomes_species_param_tests.py -v
"""

import mhcgnomes
import pytest

requires_species_param = pytest.mark.skipif(
    "species" not in __import__("inspect").signature(mhcgnomes.parse).parameters,
    reason=f"mhcgnomes {mhcgnomes.__version__} does not support species= param (needs >=3.9.0)",
)
pytestmark = requires_species_param


def parse_with_species(gene: str, species: str) -> str | None:
    """Parse gene with species= param, return species name or None on failure."""
    try:
        r = mhcgnomes.parse(gene, species=species)
        tp = type(r).__name__
        if tp in ("Gene", "Allele", "AlleleWithoutGene"):
            return getattr(getattr(r, "species", None), "name", None)
    except Exception:
        pass
    return None


def assert_species_match_or_fail(gene: str, species: str):
    """Assert that parse either returns the requested species or fails cleanly."""
    result_species = parse_with_species(gene, species)
    if result_species is not None:
        assert result_species == species or species.lower() in result_species.lower(), (
            f"parse({gene!r}, species={species!r}) returned {result_species!r} "
            f"instead of {species!r} — species= param was ignored"
        )


# ---------------------------------------------------------------------------
# MHCIIB: resolves to Tyto alba for 26 different bird (and fish) species.
# MHCIIB should either be a Gnathostomata-level gene (accepted for any
# species) or should fail for species that don't define it.
# ---------------------------------------------------------------------------


class TestMHCIIBSpeciesEnforcement:
    @pytest.mark.parametrize(
        "species",
        [
            "Milvus milvus",
            "Tachymarptis melba",
            "Limosa limosa",
            "Pagophila eburnea",
            "Calidris pugnax",
            "Ciconia ciconia",
            "Ardea cinerea",
            "Alcedo atthis",
            "Tauraco hartlaubi",
            "Dendrocopos major",
            "Bulweria bulwerii",
            "Procellaria aequinoctialis",
            "Calonectris diomedea",
            "Halobaena caerulea",
            "Oceanodroma castro",
            "Oceanodroma monteiroi",
            "Pachyptila belcheri",
            "Pachyptila desolata",
            "Pagodroma nivea",
            "Diomedea exulans",
            "Bubulcus ibis",
            "Columba livia",
            "Rallus aquaticus",
            "Jynx torquilla",
            "Philetairus socius",
            "Odontobutis potamophilus",
        ],
    )
    def test_mhciib_does_not_return_tyto_alba(self, species):
        """MHCIIB with species={species} must not return Tyto alba."""
        assert_species_match_or_fail("MHCIIB", species)


# ---------------------------------------------------------------------------
# H2-* genes: UniProt annotation transfer gives mouse gene names to
# fish/bird species. species= param must not return Mus musculus.
# ---------------------------------------------------------------------------


class TestH2SpeciesEnforcement:
    @pytest.mark.parametrize(
        "gene,species",
        [
            ("H2-EB1", "Lonchura striata"),
            ("H2-EB1", "Astyanax mexicanus"),
            ("H2-D1", "Astyanax mexicanus"),
            ("H2-D1", "Merluccius polli"),
            ("H2-DMb1", "Liparis tanakae"),
            ("H2-DMB1", "Nibea albiflora"),
            ("H2-Q10", "Astyanax mexicanus"),
            ("H2-L", "Astyanax mexicanus"),
        ],
    )
    def test_h2_gene_does_not_return_mouse(self, gene, species):
        """H2 gene with non-mouse species= must not return Mus musculus."""
        assert_species_match_or_fail(gene, species)


# ---------------------------------------------------------------------------
# Chicken B-complex genes (BLB1, BLB2, BF2): defined for Gallus gallus
# but used by other galliform birds. species= must not return chicken.
# ---------------------------------------------------------------------------


class TestChickenBComplexSpeciesEnforcement:
    @pytest.mark.parametrize(
        "gene,species",
        [
            ("BLB2", "Lagopus scotica"),
            ("BLB2", "Numida meleagris"),
            ("BLB2", "Lyrurus tetrix"),
            ("BLB1", "Lagopus scotica"),
            ("BF2", "Numida meleagris"),
        ],
    )
    def test_blb_bf_does_not_return_chicken(self, gene, species):
        """B-complex gene with non-chicken galliform must not return Gallus gallus."""
        assert_species_match_or_fail(gene, species)


# ---------------------------------------------------------------------------
# DDB1: defined for Coturnix japonica (quail), but used by Tropheus
# (cichlid fish). species= must not return quail.
# ---------------------------------------------------------------------------


class TestDDB1SpeciesEnforcement:
    @pytest.mark.parametrize(
        "species",
        [
            "Tropheus sp. T1378",
            "Tropheus sp. T1002",
            "Tropheus sp. T1373",
            "Tropheus sp. T1368",
            "Tropheus sp. T1001",
        ],
    )
    def test_ddb1_does_not_return_quail(self, species):
        """DDB1 with Tropheus species= must not return Coturnix japonica."""
        assert_species_match_or_fail("DDB1", species)


# ---------------------------------------------------------------------------
# UAA1: defined for Oreochromis niloticus (tilapia), but used by Oryzias
# ricefish species. species= must not return tilapia.
# ---------------------------------------------------------------------------


class TestUAA1SpeciesEnforcement:
    @pytest.mark.parametrize(
        "species",
        [
            "Oryzias dancena",
            "Oryzias luzonensis",
            "Labeobarbus acutirostris",
        ],
    )
    def test_uaa1_does_not_return_tilapia(self, species):
        """UAA1 with non-tilapia species= must not return Oreochromis niloticus."""
        assert_species_match_or_fail("UAA1", species)


# ---------------------------------------------------------------------------
# ZE alleles: defined for Cyprinus carpio (carp), but used by
# Labeobarbus intermedius. species= must not return carp.
# ---------------------------------------------------------------------------


class TestZESpeciesEnforcement:
    @pytest.mark.parametrize(
        "gene",
        ["ZE*0301", "ZE*0102", "ZE*0101", "ZE*0501", "ZE*0401", "ZE*0201", "ZE"],
    )
    def test_ze_does_not_return_carp(self, gene):
        """ZE allele with Labeobarbus species= must not return Cyprinus carpio."""
        assert_species_match_or_fail(gene, "Labeobarbus intermedius")


# ---------------------------------------------------------------------------
# Other single cases where species= is ignored.
# ---------------------------------------------------------------------------


class TestOtherSpeciesEnforcement:
    def test_hladrb1_does_not_return_human(self):
        """HLADRB1 with gecko species= must not return Homo sapiens."""
        assert_species_match_or_fail("HLADRB1", "Sphaerodactylus townsendi")

    def test_dda_does_not_return_tilapia(self):
        """DDA with toad species= must not return Oreochromis niloticus."""
        assert_species_match_or_fail("DDA", "Rhinella marina")

    def test_dxa_does_not_return_carp(self):
        """DXA with scat species= must not return Cyprinus carpio."""
        assert_species_match_or_fail("DXA", "Scatophagus argus")

    def test_dbb1_does_not_return_quail(self):
        """DBB1 with turbot species= must not return Coturnix japonica."""
        assert_species_match_or_fail("DBB1", "Scophthalmus maximus")

    def test_c1_does_not_return_quail(self):
        """C1 with killifish species= must not return Coturnix japonica."""
        assert_species_match_or_fail("C1", "Kryptolebias marmoratus")

    def test_c2_does_not_return_quail(self):
        """C2 with killifish species= must not return Coturnix japonica."""
        assert_species_match_or_fail("C2", "Kryptolebias marmoratus")

    def test_drb8_does_not_return_human(self):
        """DRB8 with stone flounder species= must not return Homo sapiens."""
        assert_species_match_or_fail("DRB8", "Kareius bicoloratus")


# ---------------------------------------------------------------------------
# Positive tests: species= should work for these known-good cases.
# ---------------------------------------------------------------------------


class TestSpeciesParamPositive:
    """These should parse successfully with the correct species."""

    @pytest.mark.parametrize(
        "gene,species",
        [
            ("A*02:01", "Homo sapiens"),
            ("K", "Mus musculus"),
            ("BF", "Gallus gallus"),
            ("UA", "Crocodylus porosus"),
            ("DB01", "Crocodylus porosus"),
            ("DRA", "Monodelphis domestica"),
            ("DBA", "Oreochromis niloticus"),
            ("DMA", "Sphenodon punctatus"),
        ],
    )
    def test_positive_parse(self, gene, species):
        r = mhcgnomes.parse(gene, species=species)
        assert r is not None, f"parse({gene!r}, species={species!r}) returned None"
        sp = getattr(getattr(r, "species", None), "name", None)
        assert sp is not None, f"parse({gene!r}, species={species!r}) has no species"
        assert sp == species, f"Expected {species!r}, got {sp!r}"
