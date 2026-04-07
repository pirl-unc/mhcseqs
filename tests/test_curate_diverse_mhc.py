import csv
from pathlib import Path

from scripts.curate_diverse_mhc import curate_row, normalize_gene, resolve_gene_annotation

ROOT = Path(__file__).resolve().parent.parent


def _load_mouse_class_i_sequence() -> str:
    with open(ROOT / "mhcseqs" / "mouse_h2_sequences.csv", "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["mhc_class"] == "I":
                return row["sequence"]
    raise AssertionError("expected a mouse class I sequence")


def _load_b2m_sequence() -> str:
    with open(ROOT / "mhcseqs" / "b2m_sequences.csv", "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            return row["sequence"]
    raise AssertionError("expected a B2M sequence")


def test_resolve_gene_annotation_demotes_opaque_species_numbering():
    assert normalize_gene("Crpo94", "", "Crpo") == ""
    assert normalize_gene("Crpo-Crpo94", "", "Crpo") == ""
    assert resolve_gene_annotation("Crpo94", "", "Crpo") == ("", "Crpo94", "opaque_unassigned")
    assert resolve_gene_annotation("Crpo-Crpo94", "", "Crpo") == ("", "Crpo94", "opaque_unassigned")


def test_normalize_gene_fixes_lowercase_and_transferred_prefixes():
    assert normalize_gene("orni-dba", "", "Orni") == "Orni-DBA"
    assert normalize_gene("hla-dqa1", "", "Xetr") == "Xetr-DQA1"
    assert normalize_gene("HLADRB1 K3G42_000618", "", "Spto") == "Spto-DRB1"
    # Transferred primate/model-organism prefixes on unrelated species
    assert normalize_gene("PATR-A", "", "Asme") == "Asme-A"
    assert normalize_gene("POPY-E", "", "Asme") == "Asme-E"
    assert normalize_gene("MAMU-DRA", "", "Lost") == "Lost-DRA"
    assert normalize_gene("Mamu-DRA", "", "Eufi") == "Eufi-DRA"
    assert normalize_gene("MAFA-A1", "", "Rhfe") == "Rhfe-A1"
    assert normalize_gene("SLA-DQB1", "", "Acox") == "Acox-DQB1"
    # GOGO on Gobio gobio is legitimate (prefix collision), not an artifact
    assert normalize_gene("Gogo-DAB1", "", "Gogo") == "Gogo-DAB1"
    # Legitimate uses on the actual source species are preserved
    assert normalize_gene("MAMU-DRA", "", "Mamu") == "Mamu-DRA"
    assert normalize_gene("PATR-A", "", "Patr") == "Patr-A"


def test_normalize_gene_preserves_literature_prefixes():
    assert normalize_gene("PochUA", "", "Zhch") == "Poch-UA"
    assert normalize_gene("Hyam_DAB1", "", "Hyam") == "Hyam-DAB1"


def test_curate_row_rescues_structurally_valid_class_i_without_gene():
    row = {
        "uniprot_accession": "TEST001",
        "gene_names": "",
        "protein_name": "MHC class I antigen",
        "organism": "Example species",
        "organism_id": "1",
        "is_fragment": "False",
        "source_group": "amphibian",
        "sequence": _load_mouse_class_i_sequence(),
    }
    curated, stats = curate_row(row)
    assert curated is not None
    assert curated["gene"] == ""
    assert curated["raw_gene_label"] == ""
    assert curated["gene_status"] == "missing"
    assert curated["mhc_class"] == "I"
    assert curated["chain"] == "alpha"
    assert stats["rescued_no_gene"] == 1


def test_curate_row_infers_b2m_from_protein_name():
    row = {
        "uniprot_accession": "TEST002",
        "gene_names": "",
        "protein_name": "Beta-2-microglobulin",
        "organism": "Example species",
        "organism_id": "1",
        "is_fragment": "False",
        "source_group": "amphibian",
        "sequence": _load_b2m_sequence(),
    }
    curated, stats = curate_row(row)
    assert curated is not None
    assert curated["gene"] == "B2M"
    assert curated["raw_gene_label"] == ""
    assert curated["gene_status"] == "inferred"
    assert curated["chain"] == "B2M"
    assert stats["kept"] == 1


def test_curate_row_preserves_opaque_label_without_promoting_it_to_gene():
    row = {
        "uniprot_accession": "TEST004",
        "gene_names": "Crpo94",
        "protein_name": "MHC class I antigen",
        "organism": "Crocodylus porosus",
        "organism_id": "8502",
        "is_fragment": "True",
        "source_group": "reptile_crocodylia",
        "sequence": "YFYTGVSDPGLDVPHFTAVGYVDDQQILCYNSEMRSPEPRGARVRAPSAHTCGTGNQVLAVLAVGISIQPGPLLYRYNQSQT",
    }
    curated, stats = curate_row(row)
    assert curated is not None
    assert curated["gene"] == ""
    assert curated["raw_gene_label"] == "Crpo94"
    assert curated["gene_status"] == "opaque_unassigned"
    assert stats["rescued_no_gene"] == 1
    assert stats["opaque_gene_label"] == 1


def test_curate_row_still_rejects_unstructured_no_gene_sequence():
    row = {
        "uniprot_accession": "TEST003",
        "gene_names": "",
        "protein_name": "MHC class I antigen",
        "organism": "Example species",
        "organism_id": "1",
        "is_fragment": "False",
        "source_group": "amphibian",
        "sequence": "A" * 120,
    }
    curated, stats = curate_row(row)
    assert curated is None
    assert stats["no_gene"] == 1
