import csv

from mhcseqs.domain_parsing import AlleleRecord
from mhcseqs.pipeline import (
    FULL_FIELDS,
    FUNCTIONAL_GROOVE_STATUSES,
    RAW_FIELDS,
    _candidate_tokens,
    _emit_full_row,
    _extract_source_id,
    _infer_chain,
    _load_b2m_references,
    _load_diverse_mhc_references,
    _load_mouse_h2_references,
    _looks_like_nucleotide,
    build_raw_index,
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
    # Groove columns merged into full fields
    assert "groove1" in FULL_FIELDS
    assert "groove2" in FULL_FIELDS
    assert "groove_seq" in FULL_FIELDS
    assert "ig_domain" in FULL_FIELDS
    assert "anchor_type" in FULL_FIELDS
    assert "domain_architecture" in FULL_FIELDS
    assert "domain_spans" in FULL_FIELDS


def test_functional_groove_statuses():
    assert "ok" in FUNCTIONAL_GROOVE_STATUSES
    assert "inferred_from_alpha3" in FUNCTIONAL_GROOVE_STATUSES
    assert "alpha1_only" in FUNCTIONAL_GROOVE_STATUSES
    assert "alpha2_only" in FUNCTIONAL_GROOVE_STATUSES
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


def test_source_id_in_all_field_lists():
    assert "source_id" in RAW_FIELDS
    assert "source_id" in FULL_FIELDS


def test_extract_source_id_imgt():
    assert _extract_source_id("HLA:HLA00001 A*01:01:01:01 365 bp", "imgt") == "HLA00001"
    assert _extract_source_id("HLA:HLA02169 A*01:01:01:02N 200 bp", "imgt") == "HLA02169"


def test_extract_source_id_ipd():
    assert _extract_source_id("IPD-MHC:NHP00001 Aona-DQA1*27:01 73 bp", "ipd_mhc") == "NHP00001"


def test_extract_source_id_unknown():
    assert _extract_source_id("some random header", "imgt") == ""
    assert _extract_source_id("some random header", "other") == ""


def test_load_mouse_h2_references():
    rows = _load_mouse_h2_references()
    assert len(rows) == 30
    genes = {r["gene"] for r in rows}
    assert "K" in genes
    assert "D" in genes
    assert "L" in genes
    assert "Aa" in genes
    assert "Ab1" in genes
    assert "Ea" in genes
    assert "Eb1" in genes


def test_mouse_h2_class_i():
    rows = _load_mouse_h2_references()
    class_i = [r for r in rows if r["mhc_class"] == "I"]
    assert len(class_i) == 10
    assert all(r["chain"] == "alpha" for r in class_i)
    alleles = {r["allele_normalized"] for r in class_i}
    assert "H2-K*b" in alleles
    assert "H2-D*b" in alleles
    assert "H2-L*d" in alleles


def test_mouse_h2_class_ii():
    rows = _load_mouse_h2_references()
    class_ii = [r for r in rows if r["mhc_class"] == "II"]
    assert len(class_ii) == 20
    alpha = [r for r in class_ii if r["chain"] == "alpha"]
    beta = [r for r in class_ii if r["chain"] == "beta"]
    assert len(alpha) == 10
    assert len(beta) == 10


def test_mouse_h2_has_source_id():
    rows = _load_mouse_h2_references()
    assert all(r["source_id"] for r in rows)
    # Spot check known accessions
    by_allele = {r["allele_normalized"]: r for r in rows}
    assert by_allele["H2-K*b"]["source_id"] == "P01901"
    assert by_allele["H2-D*b"]["source_id"] == "P01899"
    assert by_allele["H2-Ab1*b"]["source_id"] == "P14483"


def test_mouse_h2_metadata():
    rows = _load_mouse_h2_references()
    for r in rows:
        assert r["species"] == "Mus musculus"
        assert r["species_category"] == "murine"
        assert r["species_prefix"] == "H2"
        assert r["source"] == "uniprot"


def test_b2m_has_source_id():
    rows = _load_b2m_references()
    for r in rows:
        assert r.get("source_id"), f"B2M entry {r['allele_normalized']} missing source_id"
    human = [r for r in rows if "human" in r["allele_normalized"]]
    assert human[0]["source_id"] == "P61769"


def test_load_diverse_references_signal_peptide_metadata_is_consistent():
    rows = _load_diverse_mhc_references()
    assert rows
    for row in rows:
        sp_len = int(row["signal_peptide_len"] or "0")
        has_sp = row["has_signal_peptide"] == "True"
        expected_has_sp = sp_len >= 15 and row["sequence"][:1].upper() == "M"
        assert has_sp == expected_has_sp
        if has_sp:
            assert row["signal_peptide_seq"] == row["sequence"][:sp_len]
        else:
            assert row["signal_peptide_seq"] == ""


def test_build_raw_index_refines_signal_peptide_before_serializing(monkeypatch, tmp_path):
    seq = "M" + "A" * 80
    fasta_path = tmp_path / "synthetic.fasta"
    out_csv = tmp_path / "raw.csv"
    fasta_path.write_text(">synthetic\n" + seq + "\n", encoding="utf-8")

    monkeypatch.setattr(
        "mhcseqs.pipeline._resolve_header_allele",
        lambda header: ("HLA-A*02:01", "A", "I", "Homo sapiens", "HLA-A*02:01"),
    )
    monkeypatch.setattr("mhcseqs.pipeline._extract_source_id", lambda header, source_label: "SYN001")
    monkeypatch.setattr(
        "mhcseqs.pipeline._try_domain_parse",
        lambda seq, *, mhc_class, gene, allele, chain="", features=None: AlleleRecord(status="ok", mature_start=27),
    )
    monkeypatch.setattr(
        "mhcseqs.pipeline.refine_signal_peptide",
        lambda seq, mature_start, species_category, mhc_class="", features=None: 29,
    )
    monkeypatch.setattr("mhcseqs.pipeline._load_b2m_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline._load_mouse_h2_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline._load_diverse_mhc_references", lambda: [])

    build_raw_index([(fasta_path, "synthetic")], out_csv)

    with open(out_csv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    assert len(rows) == 1
    row = rows[0]
    assert row["has_signal_peptide"] == "True"
    assert row["signal_peptide_len"] == "29"
    assert row["signal_peptide_seq"] == seq[:29]


def test_emit_full_row_realigns_domain_fields_after_refinement(monkeypatch):
    seq = "M" * 40 + "A" * 260
    representative = {
        "allele_normalized": "HLA-A*02:01",
        "gene": "A",
        "mhc_class": "I",
        "chain": "alpha",
        "species": "Homo sapiens",
        "species_category": "bird",
        "species_prefix": "HLA",
        "source": "test",
        "source_id": "SYN001",
    }
    groove = AlleleRecord(
        status="ok",
        mhc_class="I",
        chain="alpha",
        mature_start=27,
        groove1=seq[27:117],
        groove2=seq[117:210],
        groove_seq=seq[27:210],
        groove1_len=90,
        groove2_len=93,
        ig_domain=seq[210:300],
        ig_domain_len=90,
        tail=seq[300:],
        tail_len=len(seq) - 300,
    )

    monkeypatch.setattr(
        "mhcseqs.pipeline.refine_signal_peptide",
        lambda seq, mature_start, species_category, mhc_class="": 29,
    )

    row = _emit_full_row("HLA-A*02:01", representative, "unique", seq, groove)

    assert row["mature_start"] == "29"
    reconstructed = row["groove1"] + row["groove2"] + row["ig_domain"] + row["tail"]
    assert row["mature_sequence"] == reconstructed
    assert row["groove1"] == seq[29:117]
    assert row["domain_architecture"].startswith("signal_peptide>g_alpha1>g_alpha2>c1_alpha3")
    assert row["domain_spans"].startswith("signal_peptide:1-29;g_alpha1:30-117")
