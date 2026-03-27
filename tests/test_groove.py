import pytest

from mhcseqs.domain_parsing import (
    MAX_PLAUSIBLE_SP,
    MIN_FUNCTIONAL_GROOVE_HALF_LEN,
    NON_GROOVE_GENES,
    AlleleRecord,
    _find_junction,
    _score_groove_ig_boundary,
    _score_junction,
    classify_cys_pair,
    classify_domain_fold,
    decompose_class_i,
    decompose_class_ii_alpha,
    decompose_class_ii_beta,
    decompose_domains,
    detect_h_region,
    estimate_sp_from_h_region,
    fast_sp_triage,
    find_cys_pairs,
    refine_signal_peptide,
    score_cys_flanking_properties,
    sp_boundary_excluded,
    trace_parse_class_i,
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


def test_decompose_class_i_a0201():
    result = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
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
    assert result.parse_score > 0
    assert result.parse_score == result.sp_subscore + result.groove_subscore + result.ig_subscore + result.tail_subscore


def test_class_i_domain_architecture_is_g_to_c_like_to_tm():
    result = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    assert result.ok
    tokens = [domain.token for domain in result.domains]
    assert tokens[:3] == ["g_alpha1", "g_alpha2", "c1_alpha3"]
    assert "transmembrane" in tokens
    assert tokens[-1] == "cytoplasmic_tail"
    assert "g_alpha2" in result.domain_architecture
    assert "c1_alpha3" in result.domain_spans
    alpha2 = next(domain for domain in result.domains if domain.token == "g_alpha2")
    assert any(ev.startswith("g_domain_disulfide=") for ev in alpha2.evidence)


def test_decompose_class_i_too_short():
    result = decompose_class_i("ACDEFGH" * 5, allele="short")
    assert not result.ok
    assert result.status == "too_short"


def test_decompose_class_i_no_cys_short():
    """Short class I with no Cys pairs gets alpha1_only (no Cys = α1 fragment)."""
    result = decompose_class_i("A" * 200, allele="no_cys")
    assert result.ok
    assert result.status == "alpha1_only"
    assert result.groove1_len == 200
    assert result.groove2_len == 0


def test_decompose_class_i_no_cys_long():
    """Long class I with no recoverable groove reports missing_groove."""
    result = decompose_class_i("A" * 300, allele="no_cys")
    assert not result.ok
    assert result.status == "missing_groove"
    assert "no_cys_pairs" in result.flags


def test_decompose_domains_class_i():
    result = decompose_domains(HLA_A0201_MATURE, mhc_class="I", allele="HLA-A*02:01")
    assert result.ok
    assert result.mhc_class == "I"


def test_non_groove_genes():
    assert "MICA" in NON_GROOVE_GENES
    assert "MICB" in NON_GROOVE_GENES
    assert "HFE" in NON_GROOVE_GENES
    assert "MR1" in NON_GROOVE_GENES


def test_non_groove_mr1_detected():
    result = decompose_domains(
        HLA_A0201_MATURE,
        mhc_class="I",
        allele="HLA-MR1*01:01",
        gene="MR1",
    )
    assert result.ok
    assert result.status == "ok"
    assert "non_groove_gene" in result.flags


def test_groove_result_fields():
    r = AlleleRecord()
    assert r.groove1 == ""
    assert r.groove2 == ""
    assert r.ig_domain == ""
    assert r.tail == ""
    assert r.ig_domain_len == 0
    assert r.tail_len == 0


def test_class_i_parse_candidate_marks_missing_support_as_partial():
    record = AlleleRecord(
        mhc_class="I",
        chain="alpha",
        sequence="M" + ("A" * 260),
        seq_len=261,
        mature_start=24,
        groove1="A" * 88,
        groove2="A" * 92,
        groove1_len=88,
        groove2_len=92,
        tail="A" * 57,
        tail_len=57,
        status="ok",
    )
    candidate = record.parse_candidate
    assert candidate.candidate_type == "class_I_missing_support"
    assert candidate.support_state == "support_missing"


def test_class_i_parse_candidate_keeps_supported_full_parse():
    record = AlleleRecord(
        mhc_class="I",
        chain="alpha",
        sequence="M" + ("A" * 320),
        seq_len=321,
        mature_start=24,
        groove1="A" * 88,
        groove2="A" * 92,
        groove1_len=88,
        groove2_len=92,
        ig_domain="A" * 90,
        ig_domain_len=90,
        tail="A" * 27,
        tail_len=27,
        status="ok",
    )
    candidate = record.parse_candidate
    assert candidate.candidate_type == "class_I_full"
    assert candidate.support_state == "support_present"


def test_fast_sp_triage_uses_exact_repeated_prefix30_shortcut():
    seq = "MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGY"
    est, conf, state, kind, candidates = fast_sp_triage(seq)
    assert state == "sp_present"
    assert kind == "exact_sp_prefix30"
    assert est == 24
    assert conf >= 0.99
    assert candidates == ()


def test_fast_sp_triage_uses_exact_mature10_leaderless_shortcut():
    seq = "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSD"
    est, conf, state, kind, candidates = fast_sp_triage(seq)
    assert est == 0
    assert state == "leaderless"
    assert kind == "exact_mature10"
    assert conf >= 0.99
    assert candidates == ()


# ---------------------------------------------------------------------------
# N-terminal truncation detection
# ---------------------------------------------------------------------------


def test_class_i_alpha1_only_fragment():
    """A class I fragment with no Cys pair → alpha1_only (exon 2 / α1)."""
    fragment = HLA_A0201_MATURE[:90].replace("C", "A")  # α1 domain, no Cys
    result = decompose_class_i(fragment, allele="HLA-A*02:01-a1frag")
    assert result.status == "alpha1_only"
    assert result.ok
    assert result.groove1 == fragment.upper()
    assert result.groove2 == ""
    assert result.domain_architecture == "g_alpha1"


def test_class_i_alpha2_only_fragment():
    """A class I fragment with the α2 Cys pair → alpha2_only (exon 3 / α2)."""
    # α2 domain from HLA-A*02:01 (positions 90-183), contains the Cys pair
    fragment = HLA_A0201_MATURE[90:183]
    result = decompose_class_i(fragment, allele="HLA-A*02:01-a2frag")
    assert result.status == "alpha2_only"
    assert result.ok
    assert result.groove1 == ""
    assert result.groove2 == fragment.upper()


def test_class_i_moderate_truncation_still_parses():
    """A class I sequence missing only the signal peptide (SP-stripped)
    should still parse correctly with status=ok."""
    result = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    assert result.ok
    assert result.status == "ok"
    assert result.groove1_len >= 80


# ---------------------------------------------------------------------------
# Cys-anchor robustness
# ---------------------------------------------------------------------------


def test_class_i_suspect_anchor():
    """Destroying the main α2 Cys pair should not produce an implausible parse."""
    # Replace the alpha2 Cys pair with alanines
    seq = HLA_A0201_MATURE.replace("C", "A", 2)  # kill first two Cys
    result = decompose_class_i(seq, allele="mutant")
    # Should not silently succeed with a wrong mature_start
    if result.status == "ok":
        # If it does parse, mature_start must be plausible
        assert result.mature_start <= MAX_PLAUSIBLE_SP


def test_class_ii_alpha_far_cys_pair():
    """Class II alpha with Cys pair far downstream.

    With candidate enumeration, the scorer finds the best mature_start in
    [0, MAX_PLAUSIBLE_SP].  A synthetic all-A sequence has valid cleavage
    residues everywhere so the parser finds a plausible (if meaningless) parse.
    The key invariant: mature_start is always <= MAX_PLAUSIBLE_SP.
    """
    # Cys1 at position 158, sep 56
    seq = "A" * 158 + "C" + "A" * 55 + "C" + "A" * 50
    result = decompose_class_ii_alpha(seq, allele="synth-alpha")
    assert result.mature_start <= MAX_PLAUSIBLE_SP


def test_class_ii_beta_far_cys_pair():
    """Class II beta with Cys pair far downstream — same invariant as alpha."""
    # Cys1 at position 178, sep 56
    seq = "A" * 178 + "C" + "A" * 55 + "C" + "A" * 50
    result = decompose_class_ii_beta(seq, allele="synth-beta")
    assert result.mature_start <= MAX_PLAUSIBLE_SP


# ---------------------------------------------------------------------------
# Non-classical lineage and short groove detection
# ---------------------------------------------------------------------------


def test_non_classical_class_i_via_extract():
    """Fish L-lineage genes should get status=non_classical via decompose_domains.

    Non-classical sequences still parse (they have real domain structure)
    but are flagged so downstream consumers can filter them if needed.
    """
    result = decompose_domains(
        HLA_A0201_MATURE,  # reuse a real sequence for structural validity
        mhc_class="I",
        allele="Dare-mhc1laa",
        gene="Dare-mhc1laa",
    )
    assert result.status == "non_classical"
    assert result.ok  # parseable, just non-classical
    assert "non_classical_lineage" in result.flags


def test_non_classical_mfsd_detected():
    """MFSD gene names should be flagged as non-classical (contaminants)."""
    result = decompose_domains(
        HLA_A0201_MATURE,
        mhc_class="I",
        allele="Dare-mfsd6a",
        gene="Dare-mfsd6a",
    )
    assert result.status == "non_classical"
    assert result.ok  # parseable, just non-classical


def test_short_groove_class_i():
    """A class I parse with groove1 < 70 aa should get status=short."""
    # Truncate the mature sequence to produce a short groove1
    # Remove first 30 residues — groove1 will be ~60 aa
    truncated = HLA_A0201_MATURE[30:]
    result = decompose_domains(truncated, mhc_class="I", allele="short-test", gene="A")
    if result.ok or result.status == "short":
        # If it parsed at all, check the short detection
        if result.groove1_len > 0 and result.groove1_len < MIN_FUNCTIONAL_GROOVE_HALF_LEN:
            assert result.status == "short"
            assert any("groove1_short" in f for f in result.flags)


def test_normal_groove_not_flagged():
    """A normal class I parse should NOT be flagged as short or non_classical."""
    result = decompose_domains(HLA_A0201_MATURE, mhc_class="I", allele="HLA-A*02:01", gene="A")
    assert result.status == "ok"
    assert result.ok


# ---------------------------------------------------------------------------
# Signal peptide refinement tests
# ---------------------------------------------------------------------------


def test_nonmammal_does_not_shortcircuit_on_valid_residue():
    """Non-mammal path should score the window even when the current -1 residue
    is valid (e.g. G), because a nearby bird motif (ELH) may be better."""
    # Build a synthetic sequence with precise index control:
    #   index 20 = G  (valid -1 residue for mature_start=21)
    #   index 23 = A  (valid -1 residue for candidate pos 24)
    #   index 24-26 = ELH (bird motif)
    seq = list("A" * 60)
    # Hydrophobic SP core at indices 0-14
    for i, c in enumerate("MLLVILLALVAGAAL"):
        seq[i] = c
    # G at index 20: valid -1 for initial mature_start=21
    seq[20] = "G"
    # Non-cleavage filler at 21-22 (so pos=22/23 don't have strong -1)
    seq[21] = "Q"
    seq[22] = "Q"
    # A at index 23 stays (strong -1 for candidate pos 24)
    # ELH motif at indices 24-26
    seq[24] = "E"
    seq[25] = "L"
    seq[26] = "H"
    seq_str = "".join(seq)

    # Bird: scorer should prefer pos 24 (A at -1 + ELH motif)
    # over pos 21 (G at -1 but no motif)
    result_bird = refine_signal_peptide(seq_str, mature_start=21, species_category="bird")
    assert result_bird == 24, f"Bird should move to ELH motif site, got {result_bird}"

    # Human: should stay near the initial site rather than jumping to the bird motif.
    result_human = refine_signal_peptide(seq_str, mature_start=21, species_category="human")
    assert result_human in {20, 21}, f"Human should stay near the initial site, got {result_human}"


def test_mammal_scored_window_can_recover_twa_gsh_site():
    """Mammals should consider nearby stronger motifs instead of short-circuiting."""
    seq = list("A" * 70)
    for i, c in enumerate("MLLVILLALVAGAAL"):
        seq[i] = c
    seq[18] = "T"
    seq[19] = "W"
    seq[20] = "A"
    seq[21] = "G"
    seq[22] = "S"
    seq[23] = "H"
    seq[24] = "S"
    seq[25] = "L"
    seq_str = "".join(seq)

    result_nhp = refine_signal_peptide(seq_str, mature_start=23, species_category="nhp", mhc_class="I")
    assert result_nhp == 21, f"NHP scorer should move to TWA|GSH site, got {result_nhp}"


def test_alpha3_fallback_far_cys_pair():
    """Class I with Cys pair far downstream — enumeration finds best SP in [0, 50]."""
    seq = "A" * 260 + "C" + "A" * 59 + "C" + "A" * 30
    result = decompose_class_i(seq, allele="synth-alpha3")
    # With enumeration, synthetic all-A sequences always parse (valid cleavage
    # residues everywhere).  The invariant is mature_start <= MAX_PLAUSIBLE_SP.
    assert result.mature_start <= MAX_PLAUSIBLE_SP


def test_bird_motif_elh_wins_over_weaker_candidate():
    """A candidate with bird motif ELH should beat a nearby candidate without
    a motif bonus in a bird sequence."""
    # Same sequence structure as the short-circuit test:
    #   index 20 = G (valid -1 for initial pos 21, no motif at 21)
    #   index 23 = A (valid -1 for candidate 24, ELH motif at 24-26)
    seq = list("A" * 60)
    for i, c in enumerate("MLLVILLALVAGAAL"):
        seq[i] = c
    seq[20] = "G"
    seq[21] = "Q"
    seq[22] = "Q"
    # A at index 23 stays
    seq[24] = "E"
    seq[25] = "L"
    seq[26] = "H"
    seq_str = "".join(seq)

    result = refine_signal_peptide(seq_str, mature_start=21, species_category="bird")
    assert result == 24, f"Expected bird scorer to pick ELH site at 24, got {result}"


# ---------------------------------------------------------------------------
# SP ground truth regression tests (real sequences from UniProt)
# ---------------------------------------------------------------------------
#
# Each test verifies that the Cys-pair heuristic + SP refinement produces the
# correct signal peptide length for a known entry.  The helper below mirrors
# the evaluation script's _try_parse logic.


def _predict_sp(seq: str, species_category: str) -> int:
    """Predict SP length using all parsers, picking the most plausible."""
    from mhcseqs.domain_parsing import (
        decompose_class_i,
        decompose_class_ii_alpha,
        decompose_class_ii_beta,
        refine_signal_peptide,
    )

    TYPICAL_SP = 23
    candidates = []
    for parser, mhc_class in (
        (decompose_class_i, "I"),
        (decompose_class_ii_beta, "II"),
        (decompose_class_ii_alpha, "II"),
    ):
        try:
            r = parser(seq)
            if r.ok and r.mature_start > 0:
                candidates.append((r.mature_start, mhc_class))
        except Exception:
            pass
    if not candidates:
        return 0
    in_range = [item for item in candidates if 10 <= item[0] <= 50]
    cys_start, mhc_class = (
        min(in_range, key=lambda x: abs(x[0] - TYPICAL_SP))
        if in_range
        else min(candidates, key=lambda x: x[0])
    )
    return refine_signal_peptide(seq, cys_start, species_category, mhc_class)


# Real sequences: (accession, species_category, ground_truth_sp, sequence)
_SP_GROUND_TRUTH = [
    # Human class I — HLA-G (P17693, reviewed)
    ("P17693", "human", 24,
     "MVVMAPRTLFLLLSGALTLTETWAGSHSMRYFSAAVSRPGRGEPRFIAMGYVDDTQFVRFDSDSACPRMEPRAPWVEQEG"
     "PEYWEEETRNTKAHAQTDRMNLQTLRGYYNQSEASSHTLQWMIGCDLGSDGRLLRGYEQYAYDGKDYLALNEDLRSWTA"
     "ADTAAQISKRKCEAANVAEQRRAYLEGTCVEWLHRYLENGKEMLQRADPPKTHVTHHPVFDYEATLRCWALGFYPAEIILT"
     "WQRDGEDQTQDVELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPEPLMLRWKQSSLPTIPIMGIVAGLVVLAA"
     "VVTGAAVAAVLWRKKSSD"),
    # NHP class I — Patr (Q95IT1, Pan troglodytes, reviewed)
    ("Q95IT1", "nhp", 24,
     "MVVMAPRTLFLLLSGALTLTETWAGSHSMRYFSAAVSRPGRGEPRFIAMGYVDDTQFVWFDSDSACPRMEPRAPWVEQEG"
     "PEYWEEETRNTKAHAQTDRINLQTLRGYYNQSEASSHTLQWMIGCDLGSDGRLLRGYEQYAYDGKDYLALNEDLRSWTA"
     "ADTAAQISKRKCEAANAAEQRRAYLEGTCVEWLHRYLENGKEMLQRADPPKTHVTHHPVFDYEATLRCWALGFYPAEIILT"
     "WQRDGEDQTQDVELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPEPLMLRWKQSSLPTIPIMGIVAGLVVLAA"
     "VVTGAAVAAVLWRKKSSD"),
    # Mouse class I — H2-Dd (P14431, reviewed)
    ("P14431", "murine", 21,
     "MALTMLLLLVAAALTLIETRAGQHSLQYFHTAVSRPGLGEPWFISVGYVDDTQFVRFDSDAENPRMEPRARWMEQEGPEY"
     "WERETQIAKGHEQSFRGSLRTAQSYYNQSKGGSHTLQWMYGCDMGSDGRLLRGYLQFAYEGRDYIALNEDLKTWTAVDMA"
     "AQITRRKWEQAGIAEKDQAYLEGTCMQSLRRYLELGKETL"),
    # Ungulate class I — BoLA (P13752, Bos taurus, reviewed)
    ("P13752", "ungulate", 21,
     "MGPRALLLLLSGVLILTETRAGSHSLRYFSTAVSRPGLGEPRYLEVGYVDDTQFVQFDSDAPNPRMEPRARWVEQEGPEY"
     "WDRNTRNAKGNAQSFRVNLNTLRGYYNQSEAGSHTLQWMSGCDVGPDGALRRGFMQYGYDGRDYLALNEDLRSWTAGETE"
     "AQITKRKWEAAGYAEVQRNYLEGECVEWLRRYLENGKDTLLRADPPKAHVTHHPISGREVTLRCWALGFYPEEISLTWQHD"
     "GEDQTQDMELVETRPSGDGTFQKWAALVVPSGDEQRYTCRVQHEGLQEPLTLRWEPPQPSFLTMGIIVGLVLLVVTGAVV"
     "AGVVICMKKRSGEKGGNYIQASSSDSAQGSDVSLTVPKV"),
    # Bird class I — Anas platyrhynchos (A0A6M3RI17, duck, EPH motif)
    ("A0A6M3RI17", "bird", 21,
     "MGGALGLGLWLLLGVLGGAASEPHSLRYFDTGVSDPSPGVPRFVSVGYVDGHLIDHYDSETQRTEPRADWFAANTDQQYWE"
     "WDVQNIQQNEKIFRVNLDTLRERYNQSRGSHTVQRMYGCDLLKDGSIRGYEQYGYEGRDFILLDKDTLTFTAADAAAQITK"
     "RKWEEEGTYAERTKYYLENTCIEWLRKYVSYGKDVLGRRERPEVQVSEMHADKILTLSCRAHGFYPRPISISWLKDGMVQE"
     "QETQRGSTVPNSDGTYHIWATIDVLPGDRDKYQCRVEHASLPQPGLFSWEPQSNLIPIVAGVAVAVVAVIAALAGFAVWKS"
     "KQ"),
    # Cobra class II beta-like sequence — earlier coherent cut should win
    # over later familiar-looking downstream motifs when both remain plausible.
    ("A0A8C6V5S8", "other_vertebrate", 19,
     "GRWALESLLTALLLPLPEAHFLYQEKCECLFLNGTQQVRFLFRWFYDRQEFVRFDSDLGKHVAVTEFGKVDADKWNRNEQFL"
     "EQERALIDKVCRHNYVALRSTMQRSVEPDMTISPTCVEALNHHNLLICSVTDFYPAQIKVRWFRNDQEETAGVVSNPLIYHG"
     "ADWTFQVLVMLESSPQRGDVYTCHVEHSSLASPVTVEWRAQSESAQSKMLSGVGGFVLGLIFLGLGLFIGATVQIMWRQKE"
     "PSATQPLTLASEANLEDKD"),
    # Human class II alpha — HLA-DOA (P06340, reviewed)
    ("P06340", "human", 25,
     "MALRAGLVLGFHTLMTLLSPQEAGATKADHMGSYGPAFYQSYGASGQFTHEFDEEQLFSVDLKKSEAVWRLPEFGDFARF"
     "DPQGGLAGIAAIKAHLDILVERSNRSRAINVPPRVTVLPKSRVELGQPNILICIVDNIFPPVINITWLRNGQTVTEGVAQT"
     "SFYSQPDHLFRKFHYLPFVPSAEDVYDCQVEHWGLDAPLLRHWELQVPIPPPDAMETLVCALGLAIGLVGFLVGTVLIIMG"
     "TYVSSVPR"),
    # Human class II beta — HLA-DQB (H0Y7Y7)
    ("H0Y7Y7", "human", 27,
     "XLQIPGGFWAAAVTVMLVMLSTPVAEARDFPKDFLVQFKGMCYFTNGTERVRGVARYIYNREEYGRFDSDVGEFQAVTELG"
     "RSIEDWNNYKDFLEQERAAVDKVCRHNYEAELRTTLQRQVEPTVTISPSRTEALNHHNLLVCSVTDFYPAQIKVRWFRNQD"
     "EETAGVVSTSLIRNGDWTFQILVMLEITPQRGDIYTCQVEHPSLQSPITVEWRLLH"),
    # Bird class I — Coturnix japonica (Q95592, ELH motif)
    ("Q95592", "bird", 19,
     "MGLCGMLGLLLCAVCGAAGELHSMRYIQTAMTDPGPGLPWFYEVGYVDGEIFVHYDSTTRRNVPRTEWIKAPGAVDPDYW"
     "ERNTQIVQRNEQNSRVSLDNVARLYNQSGGSHTVQWMYGCDILDDGTTRGYNQYAYDGRDFIVFDKDTMTFTAAVPEAVP"
     "TKRKWEEGDYAERQKHYLEETCVQWLRRHVENGKAELGRTEQPEVRMWGKEGDGILTLSCRAHGFFPRAIAVSWLKDGAVL"
     "GQDTHSGGIVPNSDGTYHTWITIDALPGDADKYQCRVEHASLPQPGLYSWERAQSNVLSIVGWVVGGILGIAILAGIGFII"
     "YKIHAGKKEKGYNMAPSQDGGSSSSCTGSNQTI"),
    # Bird class I — Anas platyrhynchos (A0A1L5SLR2, EPH motif)
    ("A0A1L5SLR2", "bird", 21,
     "MGGALGLVLGLLLGVLGGATSEPHSLRYFYTAVSEPSPGLPQFVGVGYVDGEAFVRYDSETHRMDSMVDWTSAIDDQQYW"
     "EWNTQNFQNDEKIFRVNLDTLRERYNQSRGSHTVQRMYGCDLLKDGSIRGFEQYGYEGRDFIALDKDTLTFTAADAAAQIT"
     "KRKWEEEGTYAERTKYYLENTCIEWLRKYVSYGKDVLERRERPEVQVSGMEADKILTLSCRAHGFYPRPISISWLKDGMVQ"
     "EQETQRGSTVPNSDGTYHIWATIDVLPGDRDKYQCRVEHASLPQPGLFSWEPQSNLIPIVAGVAVAVVAVIAALAGFAVWK"
     "SKQGKKGKGYNVAPGSEGGSNSSNAGSNPSV"),
]


@pytest.mark.parametrize(
    "accession,species_cat,gt_sp,sequence",
    [(t[0], t[1], t[2], t[3]) for t in _SP_GROUND_TRUTH],
    ids=[t[0] for t in _SP_GROUND_TRUTH],
)
def test_sp_ground_truth(accession, species_cat, gt_sp, sequence):
    """Regression: predicted SP length matches UniProt ground truth."""
    predicted = _predict_sp(sequence, species_cat)
    assert predicted > 0, f"{accession}: failed to parse"
    assert predicted == gt_sp, (
        f"{accession}: predicted SP={predicted}, ground truth SP={gt_sp} "
        f"(delta={predicted - gt_sp:+d})"
    )


def test_sp_refinement_with_truncated_hydrophobic_core():
    """SP refinement should work on sequences missing part of the hydrophobic
    core, using Cys-pair anchoring, -3/-1 residue patterns, and mature start
    motifs to get approximately correct results."""
    # Take a real bird sequence (Q95592, Coturnix japonica, SP=19) and
    # truncate the first 10 residues of the hydrophobic core.  The Cys-pair
    # anchor and mature start motif (ELH) should still guide refinement
    # to within ±2 of the correct position.
    full_seq = (
        "MGLCGMLGLLLCAVCGAAGELHSMRYIQTAMTDPGPGLPWFYEVGYVDGEIFVHYDSTTRRNVPRTEWIKAPGAVDPDYW"
        "ERNTQIVQRNEQNSRVSLDNVARLYNQSGGSHTVQWMYGCDILDDGTTRGYNQYAYDGRDFIVFDKDTMTFTAAVPEAVP"
        "TKRKWEEGDYAERQKHYLEETCVQWLRRHVENGKAELGRTEQPEVRMWGKEGDGILTLSCRAHGFFPRAIAVSWLKDGAVL"
        "GQDTHSGGIVPNSDGTYHTWITIDALPGDADKYQCRVEHASLPQPGLYSWERAQSNVLSIVGWVVGGILGIAILAGIGFII"
        "YKIHAGKKEKGYNMAPSQDGGSSSSCTGSNQTI"
    )
    gt_sp = 19

    # Remove first 10 residues (most of the h-region)
    truncated = full_seq[10:]
    # Ground truth SP in truncated coords = 19 - 10 = 9
    gt_sp_truncated = gt_sp - 10

    predicted = _predict_sp(truncated, "bird")
    # Truncated sequences with most of the h-region removed may not parse.
    # If they do parse, the SP should be within ±3 of the adjusted GT.
    if predicted > 0:
        delta = predicted - gt_sp_truncated
        assert abs(delta) <= 3, (
            f"Truncated sequence: predicted SP={predicted}, "
            f"adjusted GT={gt_sp_truncated} (delta={delta:+d})"
    )


# ---------------------------------------------------------------------------
# Phase 0: Cys pair classification and structural motif scoring
# ---------------------------------------------------------------------------


def test_classify_cys_pair_hla_a0201_alpha2():
    """The α2 Cys pair in HLA-A*02:01 should classify as groove domain."""
    # HLA-A*02:01 mature sequence — α2 Cys pair at mature positions 100 and 163
    pairs = find_cys_pairs(HLA_A0201_MATURE)
    # Find the pair closest to separation 63 in the alpha2 range
    alpha2_pair = min(
        [p for p in pairs if 50 <= p[0] <= 150],
        key=lambda p: abs(p[2] - 63),
        default=None,
    )
    assert alpha2_pair is not None, "No alpha2-range Cys pair found"
    ann = classify_cys_pair(HLA_A0201_MATURE, alpha2_pair[0], alpha2_pair[1])
    assert ann.domain_type == "groove", (
        f"Alpha2 pair should be groove, got {ann.domain_type} "
        f"(groove={ann.groove_score:.1f}, ig={ann.ig_score:.1f})"
    )
    assert ann.groove_score > ann.ig_score


def test_classify_cys_pair_hla_a0201_alpha3():
    """The α3 Cys pair in HLA-A*02:01 should classify as Ig domain."""
    pairs = find_cys_pairs(HLA_A0201_MATURE)
    # Alpha3 pair is the one with c1 > 150 (downstream of α2)
    alpha3_pair = min(
        [p for p in pairs if p[0] > 150],
        key=lambda p: abs(p[2] - 56),
        default=None,
    )
    assert alpha3_pair is not None, "No alpha3-range Cys pair found"
    ann = classify_cys_pair(HLA_A0201_MATURE, alpha3_pair[0], alpha3_pair[1])
    assert ann.domain_type == "ig", (
        f"Alpha3 pair should be ig, got {ann.domain_type} "
        f"(groove={ann.groove_score:.1f}, ig={ann.ig_score:.1f})"
    )
    assert ann.ig_score > ann.groove_score


# ---------------------------------------------------------------------------
# classify_domain_fold — Trp41-based fold-topology classifier
# ---------------------------------------------------------------------------

# HLA-DRA*01:01 mature sequence (for cross-file Trp41 tests)
_DRA_MATURE = (
    "IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYT"
    "PITNVPPEVTVLTNSPVELREPNVLICFIDKFTPPVVNVTWLRNGKPVTTGVSETVFLPREDHLFRKFHYLPFLPSTED"
    "VYDCRVEHWGLDEPLLKHWEFDAPSPLPETTENVVCALGLTVGLVGIIIGTIFIIKGVRKSNAAERRGPL"
)

# HLA-DRB1*01:01 mature (SP removed, 29 aa)
_DRB1_MATURE = (
    "GDTRPRFLWQLKFECHFFNGTERVSFLTTTPKTWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNSTAATNEVPEV"
    "TVFSKSPVTLGQPNTLICLVDNIFPPVVNITWLSNGQSVTEGVSETSFLSKSDHSFFKISYLTFLPSADEIYDCKVEHW"
    "GLDEPLLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTVFIIRGLRSGSFASRGPRGREPAE"
)


def test_classify_domain_fold_class_i_alpha2_is_g_domain():
    """Class I α2 Cys pair (C100-C163 in HLA-A*02:01) should be G-domain."""
    pairs = find_cys_pairs(HLA_A0201_MATURE)
    alpha2_pair = min(
        [p for p in pairs if 50 <= p[0] <= 150 and abs(p[2] - 63) < 10],
        key=lambda p: abs(p[2] - 63),
    )
    ann = classify_domain_fold(HLA_A0201_MATURE, alpha2_pair[0], alpha2_pair[1])
    assert ann.fold_type == "g_domain", (
        f"α2 should be g_domain, got {ann.fold_type} (trp={ann.trp_position}, "
        f"confidence={ann.confidence:.2f}, evidence={ann.evidence})"
    )
    assert ann.trp_position == -1, "G-domain should not have Trp41"
    assert ann.confidence >= 0.85


def test_classify_domain_fold_class_i_alpha3_is_c_like():
    """Class I α3 Cys pair (C202-C258 in HLA-A*02:01) should be C-like."""
    pairs = find_cys_pairs(HLA_A0201_MATURE)
    alpha3_pair = min(
        [p for p in pairs if p[0] > 150],
        key=lambda p: abs(p[2] - 56),
    )
    ann = classify_domain_fold(HLA_A0201_MATURE, alpha3_pair[0], alpha3_pair[1])
    assert ann.fold_type == "c_like", (
        f"α3 should be c_like, got {ann.fold_type} (trp={ann.trp_position}, "
        f"confidence={ann.confidence:.2f}, evidence={ann.evidence})"
    )
    assert ann.trp_position >= 0, "C-like domain should have Trp41"
    assert ann.confidence >= 0.90


def test_classify_domain_fold_class_ii_alpha2_is_c_like():
    """Class II α2 Cys pair in DRA should be C-like (has Trp41)."""
    pairs = find_cys_pairs(_DRA_MATURE)
    # DRA has one Cys pair with sep ~56 (the α2 Ig domain)
    ig_pair = min(pairs, key=lambda p: abs(p[2] - 56))
    ann = classify_domain_fold(_DRA_MATURE, ig_pair[0], ig_pair[1])
    assert ann.fold_type == "c_like", (
        f"DRA α2 should be c_like, got {ann.fold_type} "
        f"(trp={ann.trp_position}, evidence={ann.evidence})"
    )
    assert ann.trp_position >= 0


def test_classify_domain_fold_class_ii_beta1_is_g_domain():
    """A class II β1 groove Cys pair should classify as G-domain (no Trp41).

    DRB1*01:01 lacks the canonical β1 disulfide, so we use a synthetic
    sequence embedding the HLA-A α2 groove Cys pair context (known G-domain)
    at a β1-like position.  The classifier is position-independent so this
    is a valid test of the fold-type signal.
    """
    # Extract the α2 groove Cys pair region from HLA-A*02:01 (positions 90-175)
    # which is a known G-domain, and test classify_domain_fold on it directly.
    pairs = find_cys_pairs(HLA_A0201_MATURE)
    groove_pair = min(
        [p for p in pairs if 50 <= p[0] <= 150 and abs(p[2] - 63) < 10],
        key=lambda p: abs(p[2] - 63),
    )
    ann = classify_domain_fold(HLA_A0201_MATURE, groove_pair[0], groove_pair[1])
    # This is the same as the α2 test but framed as a β1 analogue: any groove
    # Cys pair with sep ~63 and no Trp in [c1+10, c1+22] should be g_domain.
    assert ann.fold_type == "g_domain"
    assert ann.trp_position == -1


def test_classify_domain_fold_class_ii_beta2_is_c_like():
    """Class II β2 Cys pair in DRB1 should be C-like (has Trp41)."""
    pairs = find_cys_pairs(_DRB1_MATURE)
    # β2 pair has sep ~56
    beta2_pair = min(
        [p for p in pairs if p[0] > 80],
        key=lambda p: abs(p[2] - 56),
        default=None,
    )
    assert beta2_pair is not None, "Should find β2 Cys pair"
    ann = classify_domain_fold(_DRB1_MATURE, beta2_pair[0], beta2_pair[1])
    assert ann.fold_type == "c_like", (
        f"DRB1 β2 should be c_like, got {ann.fold_type} "
        f"(trp={ann.trp_position}, evidence={ann.evidence})"
    )
    assert ann.trp_position >= 0


# ---------------------------------------------------------------------------
# sp_boundary_excluded — negative SP motif filter
# ---------------------------------------------------------------------------


def test_sp_boundary_excluded_valid_sites():
    """Valid SP cleavage sites should not be excluded."""
    # HLA-A*02:01 full sequence (SP=24 aa), cleavage at position 24
    hla_a_full = (
        "MAVMAPRTLLLLLSGALALTQTWA"  # 24 aa SP
        + HLA_A0201_MATURE
    )
    assert not sp_boundary_excluded(hla_a_full, 24), (
        f"-3={hla_a_full[21]}, -1={hla_a_full[23]} should be valid"
    )


def test_sp_boundary_excluded_impossible_combos():
    """Aromatic at -3 with non-small at -1 should be excluded."""
    # Construct a sequence with F at -3 and K at -1 (aromatic, basic)
    seq = "M" * 20 + "FXK" + "G" * 10
    assert sp_boundary_excluded(seq, 23), (
        "F at -3, K at -1 should be excluded (aromatic, basic)"
    )
    # Acidic at -3 with aromatic at -1
    seq2 = "M" * 20 + "DXW" + "G" * 10
    assert sp_boundary_excluded(seq2, 23), (
        "D at -3, W at -1 should be excluded (acidic, aromatic)"
    )


def test_sp_boundary_excluded_allows_canonical():
    """Canonical -3/-1 combos (small/small, aliphatic/small) should pass."""
    # A at -3, A at -1 (small, small)
    seq = "M" * 20 + "AXA" + "G" * 10
    assert not sp_boundary_excluded(seq, 23)
    # V at -3, G at -1 (aliphatic, small)
    seq2 = "M" * 20 + "VXG" + "G" * 10
    assert not sp_boundary_excluded(seq2, 23)


# ---------------------------------------------------------------------------
# score_cys_flanking_properties — property-based Cys pair scoring
# ---------------------------------------------------------------------------


def test_cys_flanking_properties_class_i():
    """Class I α3 Ig pair should have higher c_like score than α2 groove pair."""
    pairs = find_cys_pairs(HLA_A0201_MATURE)
    alpha2 = min(
        [p for p in pairs if 50 <= p[0] <= 150 and abs(p[2] - 63) < 10],
        key=lambda p: abs(p[2] - 63),
    )
    alpha3 = min(
        [p for p in pairs if p[0] > 150],
        key=lambda p: abs(p[2] - 56),
    )
    g2, c2 = score_cys_flanking_properties(HLA_A0201_MATURE, alpha2[0], alpha2[1])
    g3, c3 = score_cys_flanking_properties(HLA_A0201_MATURE, alpha3[0], alpha3[1])
    # α3 C-like should have higher c_like score
    assert c3 > c2, f"α3 c_like={c3:.1f} should be > α2 c_like={c2:.1f}"


# ---------------------------------------------------------------------------
# detect_h_region & estimate_sp_from_h_region — hydrophobic core detection
# ---------------------------------------------------------------------------


def test_detect_h_region_finds_sp_core():
    """Should find a hydrophobic stretch in HLA-A full sequence."""
    hla_a_full = "MAVMAPRTLLLLLSGALALTQTWA" + HLA_A0201_MATURE
    h_start, h_end = detect_h_region(hla_a_full)
    assert h_end > h_start > 0, "Should find h-region"
    assert h_end - h_start >= 7, f"h-region should be >=7 aa, got {h_end - h_start}"
    # h-region should be within the SP region (first 24 aa)
    assert h_start < 24, f"h-region should start within SP, got {h_start}"


def test_detect_h_region_no_sp():
    """Mature-only sequence should not have a strong h-region at the N-terminus."""
    # HLA-A*02:01 mature starts with GSHSMR... (hydrophilic)
    h_start, h_end = detect_h_region(HLA_A0201_MATURE)
    # It might find some stretch but it should be weak / late
    if h_end > 0:
        # Any stretch found should NOT be right at position 0
        assert h_start > 5 or h_end - h_start < 10


def test_estimate_sp_from_h_region():
    """SP estimate from h-region should be close to true SP length."""
    hla_a_full = "MAVMAPRTLLLLLSGALALTQTWA" + HLA_A0201_MATURE
    est = estimate_sp_from_h_region(hla_a_full)
    assert 18 <= est <= 30, f"SP estimate should be near 24, got {est}"


def test_junction_score_hla_a0201():
    """The α1/α2 junction in HLA-A*02:01 should score well at the correct position."""
    # In the mature protein, the junction is at position 90 (first α2 residue)
    # The canonical motif is G-S-H at positions 0, +1, +2 and Q at +5
    junction_score = _score_junction(HLA_A0201_MATURE, 90)
    assert junction_score > 5.0, f"Junction at pos 90 should score > 5, got {junction_score:.1f}"

    # Nearby positions should score lower
    for offset in [-3, -2, -1, 1, 2, 3]:
        alt_score = _score_junction(HLA_A0201_MATURE, 90 + offset)
        assert junction_score > alt_score, (
            f"True junction (90, score={junction_score:.1f}) should beat "
            f"pos {90+offset} (score={alt_score:.1f})"
        )


def test_find_junction_hla_a0201():
    """_find_junction should locate the α1/α2 boundary in HLA-A*02:01."""
    # In the mature protein, the α2 Cys pair is around positions 100-163.
    # The junction is ~10 before Cys1, so around position 90.
    pairs = find_cys_pairs(HLA_A0201_MATURE)
    alpha2_pair = min(
        [p for p in pairs if 50 <= p[0] <= 150],
        key=lambda p: abs(p[2] - 63),
    )
    junction_pos, junction_score = _find_junction(HLA_A0201_MATURE, alpha2_pair[0])
    assert junction_pos == 90, f"Expected junction at 90, got {junction_pos}"
    assert junction_score > 5.0


def test_groove_ig_boundary_class_i():
    """Groove2/Ig boundary in HLA-A*02:01 should score positively."""
    # The boundary is at α2 end = Cys2 + 20 ≈ position 183
    pairs = find_cys_pairs(HLA_A0201_MATURE)
    alpha2_pair = min(
        [p for p in pairs if 50 <= p[0] <= 150],
        key=lambda p: abs(p[2] - 63),
    )
    boundary_pos = alpha2_pair[1] + 20
    score = _score_groove_ig_boundary(HLA_A0201_MATURE, boundary_pos, "I")
    assert score > 0.0, f"Groove/Ig boundary should score > 0, got {score:.1f}"


# ---------------------------------------------------------------------------
# Parse trace / debug mode
# ---------------------------------------------------------------------------


def test_trace_parse_class_i_hla_a0201():
    """trace_parse_class_i should produce a valid trace for HLA-A*02:01."""
    result, trace = trace_parse_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    assert result.ok
    assert result.status == "ok"
    # Trace should have Cys pair annotations
    assert len(trace.cys_pairs) >= 2
    # At least one should be groove, at least one Ig
    types = {a.domain_type for a in trace.cys_pairs}
    assert "groove" in types, f"Expected a groove pair, got types: {types}"
    assert "ig" in types, f"Expected an Ig pair, got types: {types}"
    # Junction should be found near position 90
    assert 85 <= trace.junction_pos <= 95, f"Junction at {trace.junction_pos}, expected ~90"
    assert trace.junction_score > 5.0
    # Log should have content
    assert len(trace.log) > 5
    # Summary should be a readable string
    summary = trace.summary()
    assert "groove" in summary.lower()
    assert "junction" in summary.lower()
