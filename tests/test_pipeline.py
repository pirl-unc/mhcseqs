from mhcseqs.pipeline import (
    FULL_FIELDS,
    FUNCTIONAL_GROOVE_STATUSES,
    GROOVE_FIELDS,
    RAW_FIELDS,
    _candidate_tokens,
    _infer_chain,
    _load_b2m_references,
    _looks_like_nucleotide,
)


def test_looks_like_nucleotide_dna():
    assert _looks_like_nucleotide("ACGTACGTACGT") is True


def test_looks_like_nucleotide_rna():
    assert _looks_like_nucleotide("ACGUACGUACGU") is True


def test_looks_like_nucleotide_protein():
    assert _looks_like_nucleotide("MKVLWAALLVTFLAGCQA") is False


def test_looks_like_nucleotide_ambiguous():
    # Pure ACGT with ambiguity codes
    assert _looks_like_nucleotide("ACGTNWSMKRY") is True


def test_looks_like_nucleotide_empty():
    assert _looks_like_nucleotide("") is False
    assert _looks_like_nucleotide(None) is False


def test_looks_like_nucleotide_short_protein():
    # Short sequence with some overlap chars (A, C, G) but also protein-only chars
    assert _looks_like_nucleotide("ACDEFGHIKLM") is False


def test_candidate_tokens_hla():
    tokens = _candidate_tokens("HLA-A*02:01:01:01 123 prot")
    assert tokens[0] == "HLA-A*02:01:01:01"


def test_candidate_tokens_pipe():
    tokens = _candidate_tokens("A*02:01 | some description | extra")
    assert "A*02:01" in tokens


def test_candidate_tokens_empty():
    assert _candidate_tokens("") == []


def test_candidate_tokens_h2():
    tokens = _candidate_tokens("H-2-Kb 200 prot")
    # H-2 prefix should be scored higher
    assert tokens[0] == "H-2-Kb"


def test_infer_chain_class_i():
    assert _infer_chain("A", "I") == "alpha"
    assert _infer_chain("B", "I") == "alpha"
    assert _infer_chain("C", "I") == "alpha"


def test_infer_chain_class_ii_alpha():
    assert _infer_chain("DRA", "II") == "alpha"
    assert _infer_chain("DQA1", "II") == "alpha"
    assert _infer_chain("DPA1", "II") == "alpha"


def test_infer_chain_class_ii_beta():
    assert _infer_chain("DRB1", "II") == "beta"
    assert _infer_chain("DQB1", "II") == "beta"
    assert _infer_chain("DPB1", "II") == "beta"


def test_infer_chain_unknown():
    assert _infer_chain("", "") == ""


def test_load_b2m_references():
    rows = _load_b2m_references()
    assert len(rows) > 0
    first = rows[0]
    assert "allele_normalized" in first
    assert "sequence" in first
    assert "gene" in first
    assert first["gene"] == "B2M"
    assert first["chain"] == "B2M"
    assert first["mhc_class"] == "I"


def test_load_b2m_has_species():
    rows = _load_b2m_references()
    species_keys = {r["allele_normalized"].replace("B2M_", "") for r in rows}
    assert "human" in species_keys
    assert "mouse" in species_keys


def test_raw_fields():
    assert isinstance(RAW_FIELDS, list)
    assert len(RAW_FIELDS) > 0
    assert "allele_raw" in RAW_FIELDS
    assert "sequence" in RAW_FIELDS
    assert "mhc_class" in RAW_FIELDS


def test_full_fields():
    assert isinstance(FULL_FIELDS, list)
    assert "two_field_allele" in FULL_FIELDS
    assert "representative_allele" in FULL_FIELDS
    assert "mature_sequence" in FULL_FIELDS
    assert "groove_status" in FULL_FIELDS
    assert "is_functional" in FULL_FIELDS


def test_groove_fields():
    assert isinstance(GROOVE_FIELDS, list)
    assert "groove1" in GROOVE_FIELDS
    assert "groove2" in GROOVE_FIELDS
    assert "groove_seq" in GROOVE_FIELDS
    assert "ig_domain" in GROOVE_FIELDS
    assert "anchor_type" in GROOVE_FIELDS


def test_functional_groove_statuses():
    assert "ok" in FUNCTIONAL_GROOVE_STATUSES
    assert "alpha3_fallback" in FUNCTIONAL_GROOVE_STATUSES
    assert "beta1_only_fallback" in FUNCTIONAL_GROOVE_STATUSES
    assert "fragment_fallback" in FUNCTIONAL_GROOVE_STATUSES
    assert "too_short" not in FUNCTIONAL_GROOVE_STATUSES
    assert "no_cys_pairs" not in FUNCTIONAL_GROOVE_STATUSES


def test_full_fields_does_not_include_raw_only_columns():
    """FULL_FIELDS should not include allele_raw, allele_normalized, or signal peptide columns."""
    assert "allele_raw" not in FULL_FIELDS
    assert "allele_normalized" not in FULL_FIELDS
    assert "has_signal_peptide" not in FULL_FIELDS
    assert "signal_peptide_seq" not in FULL_FIELDS
