import csv

import pytest

from mhcseqs.domain_parsing import AlleleRecord
from mhcseqs.pipeline import (
    FULL_FIELDS,
    FUNCTIONAL_GROOVE_STATUSES,
    RAW_FIELDS,
    ResolvedHeaderAllele,
    build_raw_index,
    emit_full_row,
    extract_source_id,
    header_tokens,
    infer_chain,
    load_b2m_references,
    load_curated_uniprot_allele_map,
    load_diverse_mhc_references,
    load_genbank_saha_references,
    load_mouse_h2_references,
    load_rat_rt1_references,
    looks_like_nucleotide,
    parse_diverse_gene_class,
    parse_domain_record,
    resolve_header_allele,
)


def test_looks_like_nucleotide_dna():
    assert looks_like_nucleotide("ACGTACGTACGT") is True


def test_looks_like_nucleotide_rna():
    assert looks_like_nucleotide("ACGUACGUACGU") is True


def test_looks_like_nucleotide_protein():
    assert looks_like_nucleotide("MKVLWAALLVTFLAGCQA") is False


def test_looks_like_nucleotide_ambiguous():
    # Pure ACGT with ambiguity codes
    assert looks_like_nucleotide("ACGTNWSMKRY") is True


def test_looks_like_nucleotide_empty():
    assert looks_like_nucleotide("") is False
    assert looks_like_nucleotide(None) is False


def test_looks_like_nucleotide_short_protein():
    # Short sequence with some overlap chars (A, C, G) but also protein-only chars
    assert looks_like_nucleotide("ACDEFGHIKLM") is False


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
    assert infer_chain("A", "I") == "alpha"
    assert infer_chain("B", "I") == "alpha"
    assert infer_chain("C", "I") == "alpha"


def test_infer_chain_class_ii_alpha():
    assert infer_chain("DRA", "II") == "alpha"
    assert infer_chain("DQA1", "II") == "alpha"
    assert infer_chain("DPA1", "II") == "alpha"


def test_infer_chain_class_ii_beta():
    assert infer_chain("DRB1", "II") == "beta"
    assert infer_chain("DQB1", "II") == "beta"
    assert infer_chain("DPB1", "II") == "beta"


def test_infer_chain_unknown():
    assert infer_chain("", "") == ""


def test_load_b2m_references():
    rows = load_b2m_references()
    assert len(rows) > 0
    first = rows[0]
    assert "allele_normalized" in first
    assert "sequence" in first
    assert "gene" in first
    assert first["gene"] == "B2M"
    assert first["chain"] == "B2M"
    assert first["mhc_class"] == "I"


def test_load_b2m_has_species():
    rows = load_b2m_references()
    species_keys = {r["allele_normalized"].replace("B2M_", "") for r in rows}
    assert "human" in species_keys
    assert "mouse" in species_keys


def test_load_b2m_preserves_curated_species_metadata():
    rows = load_b2m_references()
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
    assert extract_source_id("HLA:HLA00001 A*01:01:01:01 365 bp", "imgt") == "HLA00001"
    assert extract_source_id("HLA:HLA02169 A*01:01:01:02N 200 bp", "imgt") == "HLA02169"


def test_extract_source_id_ipd():
    assert extract_source_id("IPD-MHC:NHP00001 Aona-DQA1*27:01 73 bp", "ipd_mhc") == "NHP00001"


def test_extract_source_id_unknown():
    assert extract_source_id("some random header", "imgt") == ""
    assert extract_source_id("some random header", "other") == ""


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
    rows = load_mouse_h2_references()
    assert len(rows) == 43
    genes = {r["gene"] for r in rows}
    assert {"K", "D", "L", "T23", "Aa", "Ab1", "Ea", "Eb1"} <= genes


def test_mouse_h2_class_i():
    rows = load_mouse_h2_references()
    class_i = [r for r in rows if r["mhc_class"] == "I"]
    assert len(class_i) == 15
    assert all(r["chain"] == "alpha" for r in class_i)
    alleles = {r["allele_normalized"] for r in class_i}
    assert {"H2-K*b", "H2-K*bm8", "H2-K*wm7", "H2-K*s", "H2-D*b", "H2-L*d", "H2-T23*a", "H2-T23*b"} <= alleles


def test_mouse_h2_class_ii():
    rows = load_mouse_h2_references()
    class_ii = [r for r in rows if r["mhc_class"] == "II"]
    assert len(class_ii) == 28
    alpha = [r for r in class_ii if r["chain"] == "alpha"]
    beta = [r for r in class_ii if r["chain"] == "beta"]
    assert len(alpha) == 11
    assert len(beta) == 17
    alleles = {r["allele_normalized"] for r in class_ii}
    assert {"H2-Ab1*g7", "H2-Ab1*u", "H2-Ea*k", "H2-Eb1*d", "H2-Eb1*g7", "H2-Eb1*p", "H2-Eb1*q", "H2-Eb1*r"} <= alleles


def test_load_rat_rt1_references():
    rows = load_rat_rt1_references()
    assert len(rows) == 39
    by_allele = {r["allele_normalized"]: r for r in rows}
    for expected in (
        "RT1-A1*a",
        "RT1-A1*l",
        "RT1-A1*u",
        "RT1-A1*c",
        "RT1-A2*c",
        "RT1-A1*n",
        "RT1-A2*n",
        "RT1-Ba*l",
        "RT1-Bb*l",
        "RT1-Da*l",
        "RT1-Db1*l",
        "RT1-Ba*u",
        "RT1-Bb*u",
        "RT1-Da*u",
        "RT1-Db1*u",
        "RT1-Ba*n",
        "RT1-Bb*n",
        "RT1-Da*n",
        "RT1-Db1*n",
    ):
        assert expected in by_allele
        assert by_allele[expected]["species"] == "Rattus norvegicus"
        assert by_allele[expected]["species_prefix"] == "RT1"
        assert by_allele[expected]["source_id"]


def test_mouse_h2_has_source_id():
    rows = load_mouse_h2_references()
    assert all(r["source_id"] for r in rows)
    # Spot check known accessions
    by_allele = {r["allele_normalized"]: r for r in rows}
    assert by_allele["H2-K*b"]["source_id"] == "P01901"
    assert by_allele["H2-D*b"]["source_id"] == "P01899"
    assert by_allele["H2-Ab1*b"]["source_id"] == "P14483"


def test_mouse_h2_metadata():
    rows = load_mouse_h2_references()
    for r in rows:
        assert r["species"] == "Mus musculus"
        assert r["species_category"] == "murine"
        assert r["species_prefix"] == "H2"
        assert r["source"] == "uniprot"


def test_b2m_has_source_id():
    rows = load_b2m_references()
    for r in rows:
        assert r.get("source_id"), f"B2M entry {r['allele_normalized']} missing source_id"
    human = [r for r in rows if "human" in r["allele_normalized"]]
    assert human[0]["source_id"] == "P61769"


@pytest.fixture(scope="module")
def diverse_mhc_reference_rows():
    return load_diverse_mhc_references()


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
    curated_map = load_curated_uniprot_allele_map()
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
        "P01898": ("MusMusculus-Q10", "MusMusculus", "murine"),
        "P06339": ("MusMusculus-T23*b", "MusMusculus", "murine"),
        "P14429": ("MusMusculus-Q7", "MusMusculus", "murine"),
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


def test_diverse_references_emit_roundtrippable_accessioned_alleles(diverse_mhc_reference_rows):
    """Issue #150: every UniProt row emits a unique, mhcgnomes-parseable identifier."""
    import mhcgnomes

    from mhcseqs.alleles import parse_allele_name

    rows = diverse_mhc_reference_rows
    by_source_id = {r["source_id"]: r for r in rows}

    # Spot-check gene-backed and gene-less AccessionedAllele rows
    shark = by_source_id["O46817"]
    assert shark["two_field_allele"] == "TriakisScyllium-UAA*uniprot:O46817"
    assert shark["allele_raw"] == "O46817"
    assert shark["source_id"] == "O46817"

    medaka = by_source_id["H2LR16"]
    assert medaka["two_field_allele"] == "OryziasLatipes-UBA*uniprot:H2LR16"
    parsed_medaka = parse_allele_name(medaka["two_field_allele"])
    assert isinstance(parsed_medaka, mhcgnomes.AccessionedAllele)
    assert parsed_medaka.species.name == "Oryzias latipes"
    assert parsed_medaka.is_class1

    hamster = by_source_id["A0A061HWC8"]
    assert hamster["two_field_allele"] == "CricetulusGriseus-classI*uniprot:A0A061HWC8"
    assert hamster["mhc_class"] == "I"
    parsed_hamster = parse_allele_name(hamster["two_field_allele"])
    assert isinstance(parsed_hamster, mhcgnomes.AccessionedAllele)
    assert parsed_hamster.species.name == "Cricetulus griseus"
    assert parsed_hamster.gene is None

    # Zero bare accessions, zero whitespace, zero collisions, and 100% mhcgnomes round-trip parseability.
    seen_alleles: set[str] = set()
    for row in rows:
        allele = row["two_field_allele"]
        acc = row["source_id"]
        assert " " not in allele, f"Whitespace leaked into two_field_allele: {allele!r}"
        assert allele != acc, f"Bare UniProt accession leaked into two_field_allele: {acc}"
        assert allele not in seen_alleles, f"Duplicate two_field_allele emitted: {allele}"
        seen_alleles.add(allele)
        parsed = parse_allele_name(allele)
        assert parsed is not None, f"two_field_allele failed parse_allele_name(): {allele!r} (source_id={acc})"


def test_genbank_saha_references_use_canonical_names():
    """Issue #119: genbank_saha_sequences.csv alleles normalize to SarcophilusHarrisii-I*..."""
    rows = load_genbank_saha_references()
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
    monkeypatch.setattr("mhcseqs.pipeline.extract_source_id", lambda header, source_label: "SYN001")
    monkeypatch.setattr(
        "mhcseqs.pipeline.parse_domain_record",
        lambda seq, *, mhc_class, gene, allele, chain="", species="", features=None, stats=None: AlleleRecord(status="ok", mature_start=27),
    )
    monkeypatch.setattr(
        "mhcseqs.pipeline.refine_signal_peptide",
        lambda seq, mature_start, species_category, mhc_class="", features=None: 29,
    )
    monkeypatch.setattr("mhcseqs.pipeline.load_b2m_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline.load_mouse_h2_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline.load_rat_rt1_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline.load_diverse_mhc_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline.load_genbank_saha_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline.load_curated_genbank_references", lambda: [])

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

    row = emit_full_row("HLA-A*02:01", representative, "unique", seq, groove)

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
    row = emit_full_row("HomoSapiens-H*01:01", representative, "unique", seq, groove)
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
    rows = load_diverse_mhc_references()
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
    # Issue #143: reclassified mislabeled & unknown-class rodent/rabbit MHC chains
    assert (by_acc["A0A091CXW4"]["mhc_class"], by_acc["A0A091CXW4"]["chain"]) == ("II", "beta")
    assert (by_acc["A0A061IKF1"]["mhc_class"], by_acc["A0A061IKF1"]["chain"]) == ("II", "alpha")
    assert (by_acc["Q9XRP1"]["mhc_class"], by_acc["Q9XRP1"]["chain"]) == ("II", "alpha")
    assert (by_acc["Q30847"]["mhc_class"], by_acc["Q30847"]["chain"]) == ("II", "alpha")
    assert (by_acc["A0A8C6QEU6"]["mhc_class"], by_acc["A0A8C6QEU6"]["chain"]) == ("II", "beta")
    assert (by_acc["B0BNJ8"]["mhc_class"], by_acc["B0BNJ8"]["chain"]) == ("I", "alpha")


def test_parse_diverse_gene_class_unwraps_ortho_and_ref_safely():
    """Issue #137: unwrap ~ortho:/~ref: tags without heuristic false positives on *-like tokens."""
    assert parse_diverse_gene_class("~ortho:CastorCanadensis|H2:EB1", species="Castor canadensis")["mhc_class"] == "II"
    assert parse_diverse_gene_class("~ortho:CastorCanadensis|H2:EB1", species="Castor canadensis")["chain"] == "beta"
    assert parse_diverse_gene_class("~ref:MusMusculus|Q9XRN9:H2-Ea", species="Mus musculus")["mhc_class"] == "II"
    assert parse_diverse_gene_class("~ref:MusMusculus|Q9XRN9:H2-Ea", species="Mus musculus")["chain"] == "alpha"
    assert parse_diverse_gene_class("~ref:MusMusculus|X1:Q5k-like", species="Mus musculus") is None
    assert parse_diverse_gene_class("~ref:MusMusculus|X2:T3-like", species="Mus musculus") is None


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
    monkeypatch.setattr("mhcseqs.pipeline.load_b2m_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline.load_mouse_h2_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline.load_rat_rt1_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline.load_diverse_mhc_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline.load_genbank_saha_references", lambda: [])
    monkeypatch.setattr("mhcseqs.pipeline.load_curated_genbank_references", lambda: [])

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
    rows = load_b2m_references()
    alleles = [r["two_field_allele"] for r in rows]
    assert len(alleles) == len(set(alleles))
    assert all(len(r["sequence"]) <= 200 for r in rows)
    by_allele = {r["two_field_allele"]: r for r in rows}
    assert len(by_allele["B2M_oryzias_latipes"]["sequence"]) == 116
    assert len(by_allele["B2M_sphenodon_punctatus"]["sequence"]) == 119
    assert len(by_allele["B2M_chrysemys_picta_bellii"]["sequence"]) == 119


def test_diverse_mhc_csv_has_no_bundled_non_mhc_or_unknown_class_rows():
    """Issue #165: diverse_mhc_sequences.csv itself contains no non-MHC, obsolete, or unknown-class rows."""
    from pathlib import Path

    from mhcseqs.alleles import is_non_mhc_gene
    from mhcseqs.domain_grammar import NON_MHC_ACCESSIONS

    csv_path = Path(__file__).resolve().parent.parent / "mhcseqs" / "diverse_mhc_sequences.csv"
    with open(csv_path, encoding="utf-8", newline="") as f:
        raw_rows = list(csv.DictReader(f))

    bad_non_mhc = [r["uniprot_accession"] for r in raw_rows if is_non_mhc_gene(r["gene"]) or r["uniprot_accession"] in NON_MHC_ACCESSIONS]
    assert bad_non_mhc == []

    bad_unknown = [r["uniprot_accession"] for r in raw_rows if r["mhc_class"] in ("", "unknown") or r["chain"] in ("", "unknown")]
    assert bad_unknown == []

    by_acc = {r["uniprot_accession"]: r for r in raw_rows}
    # Verified non-MHC / hidden H13 SPP / chimeric / inactive removals
    for removed in (
        "Q1LUQ4",
        "B0UYT5",
        "A0A091DH87",
        "A0A5E4CGV5",
        "A0AA41TDD6",
        "A0AA41MEA0",
        "A0AA41MGX6",
        "A0A4Z2F9A0",
        "A0ACB7F4S7",
        "A0A8J4T2F2",
        "F1R0K2",
        "A0A8D2IDC0",
    ):
        assert removed not in by_acc

    # Verified InterPro reclassifications
    assert (by_acc["A0A8J0V155"]["mhc_class"], by_acc["A0A8J0V155"]["chain"]) == ("I", "alpha")
    assert (by_acc["A0A8T1RYK3"]["mhc_class"], by_acc["A0A8T1RYK3"]["chain"]) == ("I", "alpha")
    assert (by_acc["Q5U3R9"]["mhc_class"], by_acc["Q5U3R9"]["chain"]) == ("II", "alpha")
    assert (by_acc["A0A8C8W1F8"]["mhc_class"], by_acc["A0A8C8W1F8"]["chain"]) == ("II", "beta")


def test_classify_group_assembles_short_class_ii_overlap():
    """Issue #166: short Class II exon-2 pairs with >=80% exact overlap assemble cleanly."""
    from mhcseqs.pipeline import _classify_group

    # 88-aa Class II DPB1 exon-2 groove split into two shifted overlapping fragments (87 + 85 aa)
    full_exon2 = "ENYLYQGRQECYAFNGTQRFLERYIYNREELVRFDSDVGEFRAVTELGRPDEEYWNSQKDILEEKRAVPDRVCRHNYELDEAVTLQRR"
    frag_a = full_exon2[1:]  # 87 aa (missing leading E)
    frag_b = full_exon2[:-3]  # 85 aa (missing trailing QRR)
    members = [
        {
            "allele_normalized": "MacacaThibetana-DPB1*01:01:02",
            "two_field_allele": "MacacaThibetana-DPB1*01:01",
            "sequence": frag_a,
            "mhc_class": "II",
            "gene": "DPB1",
            "chain": "beta",
            "species": "Macaca thibetana",
            "signal_peptide_len": "0",
        },
        {
            "allele_normalized": "MacacaThibetana-DPB1*01:01:01",
            "two_field_allele": "MacacaThibetana-DPB1*01:01",
            "sequence": frag_b,
            "mhc_class": "II",
            "gene": "DPB1",
            "chain": "beta",
            "species": "Macaca thibetana",
            "signal_peptide_len": "0",
        },
    ]
    anchor, assembled, policy, _details = _classify_group("MacacaThibetana-DPB1*01:01", members)
    assert policy == "assembled_overlap"
    assert anchor is not None
    assert assembled == full_exon2


def test_classify_group_resolves_groove_nested_isoforms_and_fragments():
    """Issue #166: tail isoforms + trimmed fragments sharing nested grooves resolve as groove_equivalent."""
    from mhcseqs.pipeline import _classify_group

    a1 = "GSHSMRYFYTAVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASPRGEPRAPWVEQEGPEYWDRETQKYKRQAQTDRVSLRNLRGYYNQSEA"
    a2 = "GSHTLQRMYGCDLGPDGRLLRGYDQSAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQWRAYLEGECVEWLRRYLENGKETLQRA"
    a3 = "DPPKTHVTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPEPLTLRW"
    full_1 = "MRVMAPRTLILLLSGALALTETWA" + a1 + a2 + a3 + "EPSSQPTIPIVGIVAGLAVLAVLAVLGAVVAVVMCRRKSSGGKGGSCSQAASSNSAQGSDESLIACKA"
    full_2 = "MRVMAPRTLILLLSGALALTETWA" + a1 + a2 + a3 + "EPSSQPTIPIVGIVAGLAVLAVLAVLGAVVAVVMCRRKSSGGKGGSCSQAASSNSAQGSDESLIAZZZ"
    frag = a1 + a2[:-2]  # exon 2+3 fragment missing last 2 aa of alpha2
    members = [
        {
            "allele_normalized": "HomoSapiens-C*03:03:01:01",
            "two_field_allele": "HomoSapiens-C*03:03",
            "sequence": full_1,
            "mhc_class": "I",
            "gene": "C",
            "chain": "alpha",
            "species": "Homo sapiens",
            "signal_peptide_len": "24",
        },
        {
            "allele_normalized": "HomoSapiens-C*03:03:01:32",
            "two_field_allele": "HomoSapiens-C*03:03",
            "sequence": full_2,
            "mhc_class": "I",
            "gene": "C",
            "chain": "alpha",
            "species": "Homo sapiens",
            "signal_peptide_len": "24",
        },
        {
            "allele_normalized": "HomoSapiens-C*03:03:02",
            "two_field_allele": "HomoSapiens-C*03:03",
            "sequence": frag,
            "mhc_class": "I",
            "gene": "C",
            "chain": "alpha",
            "species": "Homo sapiens",
            "signal_peptide_len": "0",
        },
    ]
    rep, assembled, policy, _details = _classify_group("HomoSapiens-C*03:03", members)
    assert policy == "groove_equivalent"
    assert assembled is None
    assert rep is not None
    assert len(rep["sequence"]) == len(full_1)


def test_emit_full_row_preserves_truncated_class_i_leader_stub():
    """Issue #166: emit_full_row must not revert validated non-Met Class I SP stubs back to mature_start=0."""
    from mhcseqs.domain_parsing import decompose_domains
    from mhcseqs.pipeline import emit_full_row

    # Non-Met partial signal peptide stub (12 aa) upstream of canonical GSHSMRYF alpha1
    a1 = "GSHSMRYFYTSMSRPGRGQPRFIAVGYVDDTQFVRFDSDAASPRMEPRAPWIEQEGPEYWDGETRNMKAHSQTDRANLGTLRGYYNQSEA"
    a2 = "GSHTIQMMYGCDVGSDGRLLRGYEQDAYDGKDYIALNEDLRSWTAADMAAQNTQRKWEAARVAEQLRAYLEGECVEWLRRYLENGKETLQR"
    a3 = "ADPPKTHVTHHPVSDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPEPLTLRW"
    seq = "LSGVLALTQTRA" + a1 + a2 + a3
    groove = decompose_domains(seq, mhc_class="I", chain="alpha", allele="MacacaFascicularis-A1*015:03")
    assert groove.ok
    assert groove.mature_start == 12
    assert groove.groove1_len == 90

    rep = {
        "allele_normalized": "MacacaFascicularis-A1*015:03",
        "gene": "A1",
        "mhc_class": "I",
        "chain": "alpha",
        "species": "Macaca fascicularis",
        "species_category": "nhp",
        "species_prefix": "MacacaFascicularis",
        "source": "ipd_mhc",
        "source_id": "Mafa00001",
        "is_null": "False",
        "is_questionable": "False",
        "is_pseudogene": "False",
    }
    row = emit_full_row("MacacaFascicularis-A1*015:03", rep, "unique", seq, groove)
    assert row["mature_start"] == "12"
    assert row["groove1_len"] == "90"
    assert row["groove1"].startswith("GSHSMRYF")


def test_try_assemble_overlap_single_domain_91aa():
    """Issue #183: complementary Class I fragments overlapping by one full ~91-aa alpha2 exon assemble."""
    from mhcseqs.pipeline import _classify_group

    a1 = "GSHSMRYFYTSMSRPGRGQPRFIAVGYVDDTQFVRFDSDAASPRMEPRAPWIEQEGPEYWDGETRNMKAHSQTDRANLGTLRGYYNQSEA"
    a2 = "GSHTIQMMYGCDVGSDGRLLRGYEQDAYDGKDYIALNEDLRSWTAADMAAQNTQRKWEAARVAEQLRAYLEGECVEWLRRYLENGKETLQR"  # 91 aa
    a3 = "ADPPKTHVTHHPVSDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPEPLTLRW"
    tail = "EPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSLTA"
    members = [
        {
            "allele_normalized": "EquusCaballus-1*010:01",
            "two_field_allele": "EquusCaballus-1*010:01",
            "sequence": a2 + a3 + tail,
            "mhc_class": "I",
            "gene": "1",
            "chain": "alpha",
            "species": "Equus caballus",
            "signal_peptide_len": "0",
        },
        {
            "allele_normalized": "EquusCaballus-1*010:01:01",
            "two_field_allele": "EquusCaballus-1*010:01",
            "sequence": a1 + a2,
            "mhc_class": "I",
            "gene": "1",
            "chain": "alpha",
            "species": "Equus caballus",
            "signal_peptide_len": "0",
        },
    ]
    rep, assembled, policy, _details = _classify_group("EquusCaballus-1*010:01", members)
    assert policy == "assembled_overlap"
    assert rep is not None
    assert assembled == a1 + a2 + a3 + tail


def test_classify_group_discards_missing_groove_member():
    """Issue #183: out-of-frame missing_groove deposits do not block valid groove-bearing alleles."""
    from mhcseqs.pipeline import _classify_group

    valid_beta = (
        "RGSRMAALIVMLMALSPPLAWARKIQPHFLEYTKKECRFSNGTERVRFLDRYFHNGEETLRFDSDVGEFRAV"
        "TELGRPDAEYWNSQKDFLEQRRAAVDTYCRHNYGVGESFTVQRRVEPKVTVYPAKTQPLQHHNLLVCSVNGF"
        "YPGSIEVRWFRNG"
    )
    out_of_frame = "HISWSILRKSVVSPTGRSGCGSWTDTSIMEKRPCASTAIGASTERWPSWGGRTPSTGTARRTGRPGPSWCGWTGTATWWRVSVCSG"
    members = [
        {
            "allele_normalized": "OvisAries-DRB1*01:03",
            "two_field_allele": "OvisAries-DRB1*01:03",
            "sequence": valid_beta,
            "mhc_class": "II",
            "gene": "DRB1",
            "chain": "beta",
            "species": "Ovis aries",
            "signal_peptide_len": "22",
        },
        {
            "allele_normalized": "OvisAries-DRB1*01:03:02",
            "two_field_allele": "OvisAries-DRB1*01:03",
            "sequence": out_of_frame,
            "mhc_class": "II",
            "gene": "DRB1",
            "chain": "beta",
            "species": "Ovis aries",
            "signal_peptide_len": "0",
        },
    ]
    rep, assembled, policy, _details = _classify_group("OvisAries-DRB1*01:03", members)
    assert policy == "groove_equivalent"
    assert assembled is None
    assert rep is not None
    assert rep["allele_normalized"] == "OvisAries-DRB1*01:03"


def test_classify_group_imgt_canonical_core_equivalence_vs_true_conflict():
    """Issue #183: N-overhang indels and trailing-X overhangs resolve via IMGT core; coding diffs stay conflict."""
    from mhcseqs.pipeline import _classify_group

    sp = "MMVLQVSAAPRTVALTALLMVLLTSVVQG"
    dpb_core = "NYLFQGRQECYAFNGTQRFLERYIYNREEFVRFDSDVGEFRAVTELGRPDEEYWNSQKDILEEERAVPDRMCRHNYELGGPMTLQRRV"
    dpb_rest = (
        "QPRVNVSPSKKGPLQHHNLLVCHVTDFYPGSIQVRWFLNGQEETAGVVSTNLIRNGDWTFQILVMLEMTPQQGDVYTCQVEHTSLDSPVTVEWKAQSDSARSKTLTGAGGFVLGLIICGVGIFMHRRSKKVQRGSA"
    )
    # Overhang difference upstream of IMGT G-BETA pos 1 (RATPAE vs RATPE) -> groove_equivalent
    m1 = sp + "RATPAE" + dpb_core + dpb_rest
    m2 = sp + "RATPE" + dpb_core + dpb_rest
    members_equiv = [
        {
            "allele_normalized": "PanTroglodytes-DPB1*01:13:01",
            "two_field_allele": "PanTroglodytes-DPB1*01:13",
            "sequence": m1,
            "mhc_class": "II",
            "gene": "DPB1",
            "chain": "beta",
            "species": "Pan troglodytes",
            "signal_peptide_len": "29",
        },
        {
            "allele_normalized": "PanTroglodytes-DPB1*01:13:02",
            "two_field_allele": "PanTroglodytes-DPB1*01:13",
            "sequence": m2,
            "mhc_class": "II",
            "gene": "DPB1",
            "chain": "beta",
            "species": "Pan troglodytes",
            "signal_peptide_len": "29",
        },
    ]
    rep, _, policy, _ = _classify_group("PanTroglodytes-DPB1*01:13", members_equiv)
    assert policy == "groove_equivalent"
    assert rep is not None

    # Genuine substitution inside canonical IMGT G-BETA core -> stays conflict
    m_diff = sp + "RATPAE" + dpb_core[:25] + ("W" if dpb_core[25] != "W" else "A") + dpb_core[26:] + dpb_rest
    members_conflict = [
        members_equiv[0],
        dict(members_equiv[1], sequence=m_diff),
    ]
    rep_c, _, policy_c, _ = _classify_group("PanTroglodytes-DPB1*01:13", members_conflict)
    assert rep_c is None
    assert policy_c.startswith("conflict")


def test_emit_full_row_protects_class_ii_beta_mature_start_from_sp_reset(monkeypatch):
    """Issue #184: emit_full_row preserves valid Class II beta mature_start when SP refinement resets to 0."""
    from mhcseqs.domain_parsing import decompose_domains
    from mhcseqs.pipeline import emit_full_row

    sp = "MVCLKLPGGSCMTALTVTLMVLSSPLALA"
    b1 = "GDTRPRFLWQLKFECHFFNGTERVRLLERCIYNQEESVRFDSDVGEYRAVTELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVGESFTVQRR"
    b2 = "VEPKVTVYPSKTQPLQHHNLLVCSVSGFYPGSIEVRWFRNGQEEKAGVVSTGLIQNGDWTFQTLVMLETVPRSGEVYTCQVEHPSVTSPLTVEW"
    seq = sp + b1 + b2
    groove = decompose_domains(seq, mhc_class="II", chain="beta", allele="HomoSapiens-DRB1*01:01")
    assert groove.ok
    assert groove.mature_start == 29
    assert 88 <= groove.groove2_len <= 96

    # Simulate refine_signal_peptide() failing has_plausible_sp_core() and returning 0
    monkeypatch.setattr("mhcseqs.pipeline.refine_signal_peptide", lambda *a, **kw: 0)
    rep = {
        "allele_normalized": "HomoSapiens-DRB1*01:01",
        "gene": "DRB1",
        "mhc_class": "II",
        "chain": "beta",
        "species": "Homo sapiens",
        "species_category": "human",
        "species_prefix": "HomoSapiens",
        "source": "imgt",
        "source_id": "HLA00001",
        "is_null": "False",
        "is_questionable": "False",
        "is_pseudogene": "False",
    }
    row = emit_full_row("HomoSapiens-DRB1*01:01", rep, "unique", seq, groove)
    assert int(row["mature_start"]) == 29
    assert int(row["groove2_len"]) == groove.groove2_len


def test_diverse_mhc_csv_issue_184_reclassification_and_removal():
    """Issue #184: A0AA41SZB3 is Class II alpha and obsolete A0A8B7UNK5 is removed."""
    from pathlib import Path

    csv_path = Path(__file__).resolve().parent.parent / "mhcseqs" / "diverse_mhc_sequences.csv"
    with open(csv_path, encoding="utf-8", newline="") as f:
        by_acc = {r["uniprot_accession"]: r for r in csv.DictReader(f)}
    assert "A0A8B7UNK5" not in by_acc
    assert (by_acc["A0AA41SZB3"]["mhc_class"], by_acc["A0AA41SZB3"]["chain"]) == ("II", "alpha")


def test_issue_186_curated_sla_and_h2_alleles():
    """Issue #186: 5 curated SLA/H2 Class I alleles load cleanly with canonical grooves and pseudosequences."""
    from mhcseqs.domain_parsing import decompose_class_i
    from mhcseqs.imgt import align_groove_to_imgt
    from mhcseqs.pipeline import load_curated_genbank_references, load_mouse_h2_references

    by_allele = {r["two_field_allele"]: r for r in load_mouse_h2_references() + load_curated_genbank_references()}
    # NetMHCpan 34 mature positions (1-indexed)
    pos34 = [
        7,
        9,
        24,
        45,
        59,
        62,
        63,
        66,
        67,
        69,
        70,
        73,
        74,
        76,
        77,
        80,
        81,
        84,
        95,
        97,
        99,
        114,
        116,
        118,
        143,
        147,
        150,
        152,
        156,
        158,
        159,
        163,
        167,
        171,
    ]
    expected_p34 = {
        "SLA-3*02:02": "YYAEYREISDNTYVNNLYFNYHYYTWSWESYLSY",
        "SLA-1*03:02": "YYIMYREISETNYGSNLYFNYHFYTWAEDNYTWY",
        "SLA-2*01:06": "YSAMYRENVGSTDVNTLYISYSDYTWAVMSYLGY",
        "H2-K*wm7": "YVEYYGQIAGNSSVDTAYIRRQVYTWAADAYAWY",
        "H2-K*bm8": "YVSYYREKAGNSFVDTLYIVSQYYTWAELAYTWY",
    }
    for allele, exp_p34 in expected_p34.items():
        assert allele in by_allele, f"{allele} missing from curated references"
        row = by_allele[allele]
        dp = decompose_class_i(row["sequence"], allele=allele, gene=row["gene"])
        assert dp.status == "ok", f"{allele}: unexpected status {dp.status}"
        assert dp.groove1_len == 90
        assert dp.groove2_len == 93
        aln = align_groove_to_imgt(dp.groove1, dp.groove2, mhc_class="I")
        assert aln.alpha1 is not None and aln.alpha2 is not None
        mature_map = {**aln.alpha1.mature_to_aa, **aln.alpha2.mature_to_aa}
        p34 = "".join(mature_map.get(p, "X") for p in pos34)
        assert p34 == exp_p34, f"{allele}: {p34} != {exp_p34}"
        assert "X" not in p34
