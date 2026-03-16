from mhcseqs.groove import (
    MAX_PLAUSIBLE_SP,
    NON_GROOVE_GENES,
    AlleleRecord,
    extract_groove,
    find_cys_pairs,
    parse_class_i,
    parse_class_ii_alpha,
    parse_class_ii_beta,
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


def test_parse_class_i_no_cys_short():
    """Short class I with no Cys pairs gets fragment_fallback."""
    result = parse_class_i("A" * 200, allele="no_cys")
    assert result.ok
    assert result.status == "fragment_fallback"


def test_parse_class_i_no_cys_long():
    """Long class I with no Cys pairs still gets no_cys_pairs failure."""
    result = parse_class_i("A" * 300, allele="no_cys")
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


# ---------------------------------------------------------------------------
# N-terminal truncation detection
# ---------------------------------------------------------------------------


def test_class_i_truncated_fragment():
    """A class I fragment (< 200 aa) with no valid Cys pair should get
    fragment_fallback, not a hard failure."""
    # Take just the groove region (first 183 aa) but kill the Cys residues
    # so no valid pair is found — simulates a truncated fragment without anchors
    fragment = HLA_A0201_MATURE[:180].replace("C", "A")
    assert len(fragment) <= 200
    result = parse_class_i(fragment, allele="HLA-A*02:01-frag")
    assert result.status == "fragment_fallback"
    assert result.ok
    assert result.groove1 == fragment


def test_class_i_moderate_truncation_still_parses():
    """A class I sequence missing only the signal peptide (SP-stripped)
    should still parse correctly with status=ok."""
    result = parse_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    assert result.ok
    assert result.status == "ok"
    assert result.groove1_len >= 80


# ---------------------------------------------------------------------------
# Suspect anchor detection (Cys mutation)
# ---------------------------------------------------------------------------


def test_class_i_suspect_anchor():
    """Mutating both Cys of the alpha2 Ig-fold pair should trigger suspect_anchor
    or a failure status, not a silently wrong parse."""
    # Replace the alpha2 Cys pair with alanines
    seq = HLA_A0201_MATURE.replace("C", "A", 2)  # kill first two Cys
    result = parse_class_i(seq, allele="mutant")
    # Should not silently succeed with a wrong mature_start
    if result.status == "ok":
        # If it does parse, mature_start must be plausible
        assert result.mature_start <= MAX_PLAUSIBLE_SP


def test_class_ii_alpha_suspect_anchor():
    """Class II alpha with implausible mature_start should get suspect_anchor."""
    # Class II alpha CYS1_RAW range is [40, 160], mature_pos = 106.
    # Place Cys1 at position 158 → mature_start = 158 - 106 = 52 > MAX_PLAUSIBLE_SP
    seq = "A" * 158 + "C" + "A" * 55 + "C" + "A" * 50
    result = parse_class_ii_alpha(seq, allele="synth-alpha")
    assert result.status == "suspect_anchor"
    assert not result.ok


def test_class_ii_beta_suspect_anchor():
    """Class II beta with implausible mature_start should get suspect_anchor."""
    # Class II beta2 CYS1_RAW range is [100, 180], mature_pos = 116.
    # Place Cys1 at position 178 → mature_start = 178 - 116 = 62 > MAX_PLAUSIBLE_SP
    seq = "A" * 178 + "C" + "A" * 55 + "C" + "A" * 50
    result = parse_class_ii_beta(seq, allele="synth-beta")
    assert result.status == "suspect_anchor"
    assert not result.ok
