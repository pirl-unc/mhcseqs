from mhcseqs.groove import (
    NON_GROOVE_GENES,
    AlleleRecord,
    extract_groove,
    find_cys_pairs,
    parse_class_i,
)

# HLA-A*02:01 mature sequence (from UniProt P01892, no signal peptide)
HLA_A0201_MATURE = (
    "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLG"
    "TLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRA"
    "YLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGD"
    "GTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSLTACKV"
)


def test_find_cys_pairs_basic():
    seq = "A" * 20 + "C" + "A" * 60 + "C" + "A" * 20
    pairs = find_cys_pairs(seq)
    assert len(pairs) == 1
    assert pairs[0] == (20, 81, 61)


def test_find_cys_pairs_no_pairs():
    assert find_cys_pairs("AAAAAAA") == []


def test_parse_class_i_a0201():
    result = parse_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    assert result.ok
    assert result.status == "ok"
    assert result.mhc_class == "I"
    assert result.chain == "alpha"
    assert 85 <= result.groove1_len <= 95
    assert 88 <= result.groove2_len <= 98
    assert result.groove_seq == result.groove1 + result.groove2
    # Ig domain and tail should be populated for full-length sequences
    assert result.ig_domain_len > 0
    assert result.tail_len > 0


def test_parse_class_i_too_short():
    result = parse_class_i("ACDEFGH" * 5, allele="short")
    assert not result.ok
    assert result.status == "too_short"


def test_parse_class_i_no_cys():
    result = parse_class_i("A" * 200, allele="no_cys")
    assert not result.ok
    assert result.status == "no_cys_pairs"


def test_extract_groove_class_i():
    result = extract_groove(HLA_A0201_MATURE, mhc_class="I", allele="HLA-A*02:01")
    assert result.ok
    assert result.mhc_class == "I"


def test_non_groove_genes():
    assert "MICA" in NON_GROOVE_GENES
    assert "MICB" in NON_GROOVE_GENES
    assert "HFE" in NON_GROOVE_GENES


def test_groove_result_fields():
    r = AlleleRecord()
    assert r.groove1 == ""
    assert r.groove2 == ""
    assert r.ig_domain == ""
    assert r.tail == ""
    assert r.ig_domain_len == 0
    assert r.tail_len == 0
