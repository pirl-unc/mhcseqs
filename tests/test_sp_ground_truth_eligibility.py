import json
from pathlib import Path

import pytest

from scripts import sp_ground_truth_eligibility as eligibility


@pytest.fixture
def source_rows():
    # Verbatim fields from the pinned 2026_03 source bundle, not live lookups.
    fixture = Path(__file__).parent / "fixtures" / "mhc_protein_identity.json"
    return {row["Entry"]: row for row in json.loads(fixture.read_text())["rows"]}


@pytest.mark.parametrize("accession", ["P33076", "A0A087X276", "A0A061ILP8"])
def test_known_non_mhc_genes_override_mhc_words_in_protein_name(source_rows, accession):
    label = eligibility.resolve_mhc_label(source_rows[accession], {})
    assert label.disposition == "exclude_non_mhc"
    assert label.mhc_class == label.chain == ""
    assert not label.eligible


@pytest.mark.parametrize("unknown,known", [("A0A0G2JMH6", "P01903"), ("Q29967", "A0A0G2JHZ7")])
def test_vague_names_do_not_exclude_exact_mhc_sequence_matches(source_rows, unknown, known):
    candidate, reference = source_rows[unknown], source_rows[known]
    assert (candidate["Organism (ID)"], candidate["Sequence"]) == (reference["Organism (ID)"], reference["Sequence"])
    assert eligibility.resolve_mhc_label(reference, {}).eligible
    label = eligibility.resolve_mhc_label(candidate, {})
    assert label.disposition == "retain_unresolved"
    assert label.label_status == "unresolved"
    # Keywords alone cannot resolve the alpha/beta chain; do not guess it.
    assert not label.eligible


@pytest.mark.parametrize("name", ["", "Uncharacterized protein", "Ig-like domain-containing protein"])
def test_missing_identity_evidence_stays_unresolved(name):
    label = eligibility.resolve_mhc_label({"Protein names": name}, {})
    assert label.disposition == "retain_unresolved"


def test_each_gene_token_is_checked_with_its_source_species(monkeypatch):
    seen = []

    def non_mhc(gene, *, species):
        seen.append((gene, species))
        return gene == "CIITA"

    monkeypatch.setattr(eligibility, "is_non_mhc_gene", non_mhc)
    label = eligibility.resolve_mhc_label({"Gene Names": "opaque;alias,CIITA", "Organism": "Mus musculus"}, {})
    assert label.disposition == "exclude_non_mhc"
    assert seen == [(token, "Mus musculus") for token in ["opaque", "alias", "CIITA"]]


@pytest.mark.parametrize("disposition", ["include", "exclude_non_mhc", "retain_unresolved"])
def test_explicit_accession_curation_still_takes_precedence(source_rows, disposition):
    decision = {"disposition": disposition, "mhc_class": "I", "chain": "alpha", "label_status": "curated"}
    label = eligibility.resolve_mhc_label(source_rows["P33076"], {"P33076": decision})
    assert label.disposition == disposition
    assert label.label_status == "curated"
