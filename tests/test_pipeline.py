import csv

import pytest

from mhcseqs.domain_parsing import AlleleRecord
from mhcseqs.pipeline import (
    FULL_FIELDS,
    FUNCTIONAL_GROOVE_STATUSES,
    RAW_FIELDS,
    ResolvedHeaderAllele,
    _emit_full_row,
    _extract_source_id,
    _infer_chain,
    _load_b2m_references,
    _load_curated_uniprot_allele_map,
    _load_diverse_mhc_references,
    _load_genbank_saha_references,
    _load_mouse_h2_references,
    _looks_like_nucleotide,
    _parse_diverse_gene_class,
    build_raw_index,
    header_tokens,
    parse_domain_record,
    resolve_header_allele,
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


def test_header_tokens_splits_delimiters():
    assert header_tokens("HLA-A*02:01:01:01 123 prot") == ["HLA-A*02:01:01:01", "123", "prot"]
    assert "A*02:01" in header_tokens("A*02:01 | some description | extra")
    assert header_tokens("") == []


def test_resolve_header_allele_picks_allele_over_bare_gene_word():
    """Header resolution delegates to mhcgnomes without prefix hacks (#95)."""
    resolved = resolve_header_allele("a protein H-2-Kb 200 prot")
    assert resolved.allele_normalized == "MusMusculus-K*b"
    assert resolved.two_field_allele == "MusMusculus-K*b"
    assert resolved.gene == "K"
    assert resolved.species == "Mus musculus"
    assert resolved.species_prefix == "MusMusculus"
    assert resolved.allele_raw == "H-2-Kb"


def test_resolve_header_allele_marks_gene_level_pseudogenes():
    resolved = resolve_header_allele("HLA:HLA00001 H*01:01:01:01 365 bp")
    assert resolved.allele_normalized == "HomoSapiens-H*01:01:01:01"
    assert resolved.is_pseudogene is True
    assert resolved.is_null is False


def test_resolve_header_allele_raises_value_error_when_unparseable():
    with pytest.raises(ValueError, match="Failed to parse allele"):
        resolve_header_allele("unknown_header 123 prot")


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


def test_load_b2m_preserves_curated_species_metadata():
    rows = _load_b2m_references()
    sturgeon = next(r for r in rows if r["source_id"] == "Q9PRF8")
    assert sturgeon["species"] == "Acipenser baerii"
    assert sturgeon["species_category"] == "fish"
    assert sturgeon["species_prefix"] == "AcipenserBaerii"
    assert all(r["species_category"] for r in rows)


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


def test_resolve_official_imgt_header_with_bare_hla_allele():
    result = resolve_header_allele("HLA:HLA00001 A*01:01:01:01 365 bp")
    assert result.allele_normalized == "HomoSapiens-A*01:01:01:01"
    assert result.two_field_allele == "HomoSapiens-A*01:01"
    assert result.gene == "A"
    assert result.mhc_class == "I"
    assert result.species == "Homo sapiens"
    assert result.species_prefix == "HomoSapiens"
    assert result.allele_raw == "A*01:01:01:01"
    assert result.is_null is False
    assert result.is_questionable is False
    assert result.is_pseudogene is False


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


@pytest.fixture(scope="module")
def diverse_mhc_reference_rows():
    return _load_diverse_mhc_references()


def test_diverse_references_fall_back_to_curated_source_category(diverse_mhc_reference_rows):
    rows = diverse_mhc_reference_rows
    bank_vole = next(r for r in rows if r["source_id"] == "A0AAW0H2V6")
    assert bank_vole["species"] == "Myodes glareolus"
    assert bank_vole["species_category"] == "other_mammal"
    assert bank_vole["species_prefix"] == "MyodesGlareolus"
    assert all(r["species_category"] for r in rows)
    assert all(r["species_prefix"] for r in rows)


def test_diverse_references_include_felis_catus_and_curated_alleles(diverse_mhc_reference_rows):
    """Issues #117 and #119: FLA Class I heavy chains and named non-IPD reference alleles."""
    rows = diverse_mhc_reference_rows
    by_source_id = {r["source_id"]: r for r in rows}

    # Every curated UniProt accession must be present and mapped to a canonical two-field name
    curated_map = _load_curated_uniprot_allele_map()
    assert len(curated_map) >= 20
    for acc in curated_map:
        assert acc in by_source_id, f"Curated accession {acc} missing from diverse_mhc_sequences.csv"
        row = by_source_id[acc]
        assert row["two_field_allele"] != acc
        assert row["species_prefix"]

    # Spot-check key benchmark alleles from #117 and #119
    expected = {
        "Q95IT0": ("TrichosurusVulpecula-UB*01:01", "TrichosurusVulpecula", "other_mammal"),
        "A0A125R585": ("PteropusAlecto-N*01:01", "PteropusAlecto", "other_mammal"),
        "Q6I7L2": ("AnasPlatyrhynchos-UAA*01", "AnasPlatyrhynchos", "bird"),
        "A0A6M3RI17": ("AnasPlatyrhynchos-UAA*76", "AnasPlatyrhynchos", "bird"),
        "Q95482": ("FelisCatus-E*003:01", "FelisCatus", "carnivore"),
        "C6ZK69": ("FelisCatus-K*007:01", "FelisCatus", "carnivore"),
        "C6ZK72": ("FelisCatus-E*018:01", "FelisCatus", "carnivore"),
        "A0A7H1D3E9": ("FelisCatus-E*001:01", "FelisCatus", "carnivore"),
        "A0A516V164": ("FelisCatus-E*005:01", "FelisCatus", "carnivore"),
        "A0A7H1D3F5": ("FelisCatus-K*012:02", "FelisCatus", "carnivore"),
    }
    for acc, (allele_canon, sp_prefix, sp_cat) in expected.items():
        row = by_source_id[acc]
        assert row["two_field_allele"] == allele_canon
        assert row["species_prefix"] == sp_prefix
        assert row["species_category"] == sp_cat
        parsed = parse_domain_record(
            row["sequence"],
            mhc_class=row["mhc_class"],
            gene=row["gene"],
            allele=row["two_field_allele"],
            chain=row["chain"],
            species=row["species"],
        )
        assert parsed is not None and parsed.ok and parsed.groove1 and parsed.groove2

    # Felis catus Class I heavy-chain coverage (#117)
    fla_ci = [r for r in rows if r["species_prefix"] == "FelisCatus" and r["mhc_class"] == "I" and r["chain"] == "alpha"]
    assert len(fla_ci) >= 15


def test_genbank_saha_references_use_canonical_names():
    """Issue #119: genbank_saha_sequences.csv alleles normalize to SarcophilusHarrisii-I*..."""
    rows = _load_genbank_saha_references()
    assert len(rows) >= 60
    by_raw = {r["allele_raw"]: r for r in rows}
    assert by_raw["SahaI*27"]["two_field_allele"] == "SarcophilusHarrisii-I*27"
    assert by_raw["SahaI*27"]["species_prefix"] == "SarcophilusHarrisii"
    assert by_raw["SahaI*27"]["gene"] == "I"
    assert all(r["species_prefix"] == "SarcophilusHarrisii" for r in rows)


def test_load_diverse_references_signal_peptide_metadata_is_consistent(diverse_mhc_reference_rows):
    rows = diverse_mhc_reference_rows
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
        "mhcseqs.pipeline.resolve_header_allele",
        lambda header: ResolvedHeaderAllele(
            parsed=None,
            allele_raw="HLA-A*02:01",
            allele_normalized="HomoSapiens-A*02:01",
            two_field_allele="HomoSapiens-A*02:01",
            gene="A",
            mhc_class="I",
            species="Homo sapiens",
            species_prefix="HomoSapiens",
            is_null=False,
            is_questionable=False,
            is_pseudogene=False,
        ),
    )
    monkeypatch.setattr("mhcseqs.pipeline._extract_source_id", lambda header, source_label: "SYN001")
    monkeypatch.setattr(
        "mhcseqs.pipeline.parse_domain_record",
        lambda seq, *, mhc_class, gene, allele, chain="", species="", features=None, stats=None: AlleleRecord(status="ok", mature_start=27),
    )
    monkeypatch.setattr(
        "mhcseqs.pipeline.refine_signal_peptide",
        lambda seq, mature_start, species_category, mhc_class="", features=None: 29,
    )
    monkeypatch.setattr("mhcseqs.pipeline._load_b2m_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline._load_mouse_h2_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline._load_diverse_mhc_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline._load_genbank_saha_references", lambda: [])

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


def test_emit_full_row_marks_gene_level_pseudogene_not_functional():
    seq = "M" * 40 + "A" * 260
    representative = {
        "allele_normalized": "HomoSapiens-H*01:01:01:01",
        "gene": "H",
        "mhc_class": "I",
        "chain": "alpha",
        "species": "Homo sapiens",
        "species_category": "human",
        "species_prefix": "HomoSapiens",
        "source": "imgt",
        "source_id": "HLA00001",
    }
    groove = AlleleRecord(
        status="ok",
        mhc_class="I",
        chain="alpha",
        mature_start=24,
        groove1=seq[24:114],
        groove2=seq[114:207],
        groove_seq=seq[24:207],
        groove1_len=90,
        groove2_len=93,
    )
    row = _emit_full_row("HomoSapiens-H*01:01", representative, "unique", seq, groove)
    assert row["is_pseudogene"] == "True"
    assert row["is_functional"] == "False"


def test_parse_domain_record_counts_value_error_and_propagates_bugs(monkeypatch):
    from collections import Counter

    stats = Counter()

    def bad_value(*_args, **_kwargs):
        raise ValueError("bad sequence")

    monkeypatch.setattr("mhcseqs.pipeline.decompose_domains", bad_value)
    assert parse_domain_record("MAVMAP" * 20, mhc_class="I", gene="A", allele="HLA-A*02:01", stats=stats) is None
    assert stats["domain_parse_errors"] == 1

    def bug(*_args, **_kwargs):
        raise RuntimeError("unexpected bug")

    monkeypatch.setattr("mhcseqs.pipeline.decompose_domains", bug)
    with pytest.raises(RuntimeError, match="unexpected bug"):
        parse_domain_record("MAVMAP" * 20, mhc_class="I", gene="A", allele="HLA-A*02:01", stats=stats)


def test_diverse_references_filter_non_mhc_gene_rows():
    """Regression test for #135: non-MHC bystander rows must not leak from diverse_mhc_sequences.csv."""
    rows = _load_diverse_mhc_references()
    leaked = []
    banned = {"CIITA", "HM13", "H13", "RFX1", "RFX2", "RFX5", "RFXAP", "RFXANK", "PRR3", "HSD17B8", "CFB"}
    for row in rows:
        gene = (row.get("gene") or "").strip()
        tail = gene.rsplit(":", 1)[-1].strip() if gene.startswith("~") and ":" in gene else gene
        base = tail[: -len("_predicted")] if tail.lower().endswith("_predicted") else tail
        if base.upper() in banned:
            leaked.append((row["source_id"], gene))
    assert leaked == []
    # Sanity check that legitimate curated rows remain present
    accessions = {row["source_id"] for row in rows}
    assert "A0A125R585" in accessions  # Ptal-N*01:01

    # Regression checks for #137: Q31152 is Class I alpha (not B2M), Sphenodon DAA-1 is II alpha,
    # and ~ortho:/~ref: Class II genes override wrong UniProt mhc_class="I".
    by_acc = {row["source_id"]: row for row in rows}
    assert by_acc["Q31152"]["mhc_class"] == "I"
    assert by_acc["Q31152"]["chain"] == "alpha"
    assert by_acc["A0A0F7GWZ7"]["mhc_class"] == "II"
    assert by_acc["A0A0F7GWZ7"]["chain"] == "alpha"
    assert by_acc["A0A250YE21"]["mhc_class"] == "II"
    assert by_acc["A0A250YE21"]["chain"] == "beta"
    assert by_acc["P04224"]["mhc_class"] == "II"
    assert by_acc["P04224"]["chain"] == "alpha"
    assert by_acc["P01915"]["mhc_class"] == "II"
    assert by_acc["P01915"]["chain"] == "beta"
    assert by_acc["Q9TQA7"]["mhc_class"] == "II"
    assert by_acc["Q9TQA7"]["chain"] == "beta"
    # Issue #141: curated mouse H2-D*q / H2-L*q mappings for UniProt Q31145 / Q31151
    assert by_acc["Q31145"]["two_field_allele"] == "MusMusculus-D*q"
    assert by_acc["Q31145"]["gene"] == "D"
    assert by_acc["Q31151"]["two_field_allele"] == "MusMusculus-L*q"
    assert by_acc["Q31151"]["gene"] == "L"


def test_parse_diverse_gene_class_unwraps_ortho_and_ref_safely():
    """Issue #137: unwrap ~ortho:/~ref: tags without heuristic false positives on *-like tokens."""
    assert _parse_diverse_gene_class("~ortho:CastorCanadensis|H2:EB1", species="Castor canadensis")["mhc_class"] == "II"
    assert _parse_diverse_gene_class("~ortho:CastorCanadensis|H2:EB1", species="Castor canadensis")["chain"] == "beta"
    assert _parse_diverse_gene_class("~ref:MusMusculus|Q9XRN9:H2-Ea", species="Mus musculus")["mhc_class"] == "II"
    assert _parse_diverse_gene_class("~ref:MusMusculus|Q9XRN9:H2-Ea", species="Mus musculus")["chain"] == "alpha"
    assert _parse_diverse_gene_class("~ref:MusMusculus|X1:Q5k-like", species="Mus musculus") is None
    assert _parse_diverse_gene_class("~ref:MusMusculus|X2:T3-like", species="Mus musculus") is None


def test_build_raw_index_filters_non_mhc_tap_genes_and_sanitizes_question_marks(tmp_path, monkeypatch):
    """Issue #137: TAP1/TAP2 in IMGT/IPD FASTAs are excluded and '?' residues normalize to 'X'."""
    fasta = tmp_path / "test.fasta"
    fasta.write_text(
        ">HLA:HLA99001 TAP1*01:01 748 aa\n"
        + "M" * 100
        + "\n>HLA:HLA99002 TAP2*01:01 686 aa\n"
        + "M" * 100
        + "\n>IPD-MHC:BOV00001 BoLA-2*030:01 180 aa\n"
        + "GSHSMRYF?TSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEA"
        + "GSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQR\n",
        encoding="utf-8",
    )
    monkeypatch.setattr("mhcseqs.pipeline._load_b2m_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline._load_mouse_h2_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline._load_diverse_mhc_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline._load_genbank_saha_references", lambda: [])

    out_csv = tmp_path / "raw.csv"
    stats = build_raw_index([(fasta, "ipd_mhc")], out_csv)
    assert stats["skipped_non_mhc"] == 2
    assert stats["parsed"] == 1
    rows = list(csv.DictReader(out_csv.open(encoding="utf-8")))
    assert len(rows) == 1
    assert rows[0]["two_field_allele"] == "BoLA-2*030:01"
    assert "?" not in rows[0]["sequence"]
    assert rows[0]["sequence"][8] == "X"


def test_b2m_sequences_are_unique_per_species_and_reasonable_length():
    """Issue #137: b2m_sequences.csv has one canonical entry per species and no >200 aa non-B2M proteins."""
    rows = _load_b2m_references()
    alleles = [r["two_field_allele"] for r in rows]
    assert len(alleles) == len(set(alleles))
    assert all(len(r["sequence"]) <= 200 for r in rows)
    by_allele = {r["two_field_allele"]: r for r in rows}
    assert len(by_allele["B2M_oryzias_latipes"]["sequence"]) == 116
    assert len(by_allele["B2M_sphenodon_punctatus"]["sequence"]) == 119
    assert len(by_allele["B2M_chrysemys_picta_bellii"]["sequence"]) == 119
