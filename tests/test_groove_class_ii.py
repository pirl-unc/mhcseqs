from mhcseqs.domain_parsing import (
    decompose_class_ii_alpha,
    decompose_class_ii_beta,
    decompose_domains,
    is_class_ii_alpha_gene,
)

# HLA-DRA*01:01 mature sequence (from UniProt P01903, no signal peptide)
HLA_DRA0101_MATURE = (
    "IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYT"
    "PITNVPPEVTVLTNSPVELREPNVLICFIDKFTPPVVNVTWLRNGKPVTTGVSETVFLPREDHLFRKFHYLPFLPSTED"
    "VYDCRVEHWGLDEPLLKHWEFDAPSPLPETTENVVCALGLTVGLVGIIIGTIFIIKGVRKSNAAERRGPL"
)

# HLA-DRB1*01:01 full sequence (from UniProt P01911, including signal peptide)
# The parser expects full sequences from FASTA files; the β2 Cys pair
# positional range requires SP context to be in the expected raw range.
HLA_DRB10101_FULL = (
    "MVCLKLPGGSCMTALTVTLMVLSSPLALA"
    "GDTRPRFLWQLKFECHFFNGTERVSFLTTTPKTWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNSTAATNEVPEV"
    "TVFSKSPVTLGQPNTLICLVDNIFPPVVNITWLSNGQSVTEGVSETSFLSKSDHSFFKISYLTFLPSADEIYDCKVEHW"
    "GLDEPLLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTVFIIRGLRSGSFASRGPRGREPAE"
)

# HLA-DQB2 full sequence (UniProt P05538). This was a live regression where the
# parser incorrectly reused the β1 groove pair as the β2 anchor and returned
# invalid_boundaries instead of the correct β1→β2 decomposition.
HLA_DQB2_FULL = (
    "MSWKMALQIPGGFWAAAVTVMLVMLSTPVAEARDFPKDFLVQFKGMCYFTNGTERVRGVARYIYNREEYGRFDSDVGEFQAVTELGRSI"
    "EDWNNYKDFLEQERAAVDKVCRHNYEAELRTTLQRQVEPTVTISPSRTEALNHHNLLVCSVTDFYPAQIKVRWFRNDQEETAGVVSTSL"
    "IRNGDWTFQILVMLEITPQRGDIYTCQVEHPSLQSPITVEWRAQSESAQSKMLSGIGGFVLGLIFLGLGLIIRHRGQKGPRGPPPAGLLH"
)

# Chicken BLB1 full sequence (UniProt A0A1D5PLI8). This exercises the same β1/β2
# distinction in a non-mammalian chain with extra downstream cysteines.
GAGA_BLB1_FULL = (
    "MGSGRVPAAGAVLVALLALGARPAAGTRPSAFFQWTFKAECHYLNGTERVRYLVRYVYNRQEYAHFDSDVGKHVADTPLGEPQAEYWNSN"
    "AEILENRMNEVDTYCRHNYGVVESFTVQRSVEPKVRVSALQSGSLPETDRLACYVTGFYPPEIEVKWFLNGREETERVVSTDVMQNGDWT"
    "YQVLVVLETVPRRGDSYVCRVEHASLRQPISQAWEPPADAGRSKLDAELAAAPPSRCTRTPRSPGRRLGSPSGSTPPVPPLCCRIARSVP"
    "SPRSRDGVALSAYPRGPSVVRQQGRGAGGLCSSIPFSTRRWFGFFNQIYSFVFA"
)


def test_decompose_class_ii_alpha_dra0101():
    result = decompose_class_ii_alpha(HLA_DRA0101_MATURE, allele="HLA-DRA*01:01")
    assert result.ok
    assert result.status == "ok"
    assert result.mhc_class == "II"
    assert result.chain == "alpha"
    assert result.groove1_len > 0
    assert result.groove2 == ""
    assert result.groove2_len == 0
    assert result.groove_seq == result.groove1
    assert result.ig_domain_len > 0
    assert result.tail_len > 0
    # α1 domain is ~83 aa
    assert 70 <= result.groove1_len <= 95


def test_class_ii_alpha_domain_architecture_is_g_to_c_like_to_tm():
    result = decompose_class_ii_alpha(HLA_DRA0101_MATURE, allele="HLA-DRA*01:01")
    assert result.ok
    tokens = [domain.token for domain in result.domains]
    assert tokens[:2] == ["g_alpha1", "c1_alpha2"]
    assert "transmembrane" in tokens
    assert tokens[-1] == "cytoplasmic_tail"
    support = next(domain for domain in result.domains if domain.token == "c1_alpha2")
    assert any(ev.startswith("c_like_disulfide=") for ev in support.evidence)


def test_decompose_class_ii_beta_drb10101():
    result = decompose_class_ii_beta(HLA_DRB10101_FULL, allele="HLA-DRB1*01:01")
    assert result.ok
    assert result.status == "ok"
    assert result.mhc_class == "II"
    assert result.chain == "beta"
    assert result.groove2_len > 0
    assert result.groove1 == ""
    assert result.groove1_len == 0
    assert result.groove_seq == result.groove2
    assert result.ig_domain_len > 0
    assert result.tail_len > 0
    # β1 domain is ~73-93 aa depending on SP accuracy. The current holistic
    # scorer lands near the true SP=29, so the groove can be slightly shorter
    # than the old heuristic slices while still being structurally correct.
    assert 65 <= result.groove2_len <= 105


def test_class_ii_beta_domain_architecture_has_signal_and_tm():
    result = decompose_class_ii_beta(HLA_DRB10101_FULL, allele="HLA-DRB1*01:01")
    assert result.ok
    tokens = [domain.token for domain in result.domains]
    assert tokens[:3] == ["signal_peptide", "g_beta1", "c1_beta2"]
    assert "transmembrane" in tokens
    groove = next(domain for domain in result.domains if domain.token == "g_beta1")
    assert "peptide_binding_module" in groove.evidence
    assert any(ev in {"g_domain_anchor", "bounded_by_c_like_anchor"} for ev in groove.evidence)
    assert result.domain_spans.startswith("signal_peptide:")


def test_decompose_class_ii_beta_dqb2_regression():
    result = decompose_class_ii_beta(HLA_DQB2_FULL, allele="HLA-DQB2", gene="DQB2")
    assert result.ok
    assert result.status == "ok"
    assert abs(result.mature_start - 32) <= 2
    assert 70 <= result.groove2_len <= 100


def test_decompose_class_ii_beta_chicken_blb1_regression():
    result = decompose_class_ii_beta(GAGA_BLB1_FULL, allele="Gaga-BLB1", gene="BLB1")
    assert result.ok
    assert result.status == "ok"
    assert abs(result.mature_start - 26) <= 3
    assert 70 <= result.groove2_len <= 105


def test_decompose_class_ii_alpha_too_short():
    result = decompose_class_ii_alpha("ACDEFGH" * 5, allele="short")
    assert not result.ok
    assert result.status == "too_short"


def test_decompose_class_ii_beta_too_short():
    result = decompose_class_ii_beta("ACDEFGH" * 5, allele="short")
    assert not result.ok
    assert result.status == "too_short"


def test_decompose_class_ii_alpha_fragment():
    # A short sequence should trigger fragment fallback
    fragment = HLA_DRA0101_MATURE[:85]
    result = decompose_class_ii_alpha(fragment, allele="fragment")
    assert result.ok
    assert result.status == "fragment_fallback"
    assert result.groove1 == fragment.upper()


def test_decompose_class_ii_beta_fragment():
    fragment = HLA_DRB10101_FULL[29:124]  # ~95 aa from within the mature region
    result = decompose_class_ii_beta(fragment, allele="fragment")
    assert result.ok
    assert result.status == "fragment_fallback"
    assert result.groove2 == fragment.upper()


def test_decompose_domains_class_ii_alpha():
    result = decompose_domains(
        HLA_DRA0101_MATURE,
        mhc_class="II",
        chain="alpha",
        allele="HLA-DRA*01:01",
        gene="DRA",
    )
    assert result.ok
    assert result.mhc_class == "II"
    assert result.chain == "alpha"


def test_decompose_domains_class_ii_beta():
    result = decompose_domains(
        HLA_DRB10101_FULL,
        mhc_class="II",
        chain="beta",
        allele="HLA-DRB1*01:01",
        gene="DRB1",
    )
    assert result.ok
    assert result.mhc_class == "II"
    assert result.chain == "beta"


def test_decompose_domains_class_ii_infer_chain_alpha():
    # Without explicit chain, gene name "DRA" should infer alpha
    result = decompose_domains(
        HLA_DRA0101_MATURE,
        mhc_class="II",
        gene="DRA",
        allele="HLA-DRA*01:01",
    )
    assert result.ok
    assert result.chain == "alpha"


def test_decompose_domains_class_ii_infer_chain_beta():
    # Without explicit chain, gene name "DRB1" should infer beta
    result = decompose_domains(
        HLA_DRB10101_FULL,
        mhc_class="II",
        gene="DRB1",
        allele="HLA-DRB1*01:01",
    )
    assert result.ok
    assert result.chain == "beta"


def test_is_class_ii_alpha_gene():
    assert is_class_ii_alpha_gene("DRA") is True
    assert is_class_ii_alpha_gene("DQA1") is True
    assert is_class_ii_alpha_gene("DPA1") is True
    assert is_class_ii_alpha_gene("DMA") is True
    assert is_class_ii_alpha_gene("DOA") is True


def test_is_class_ii_alpha_gene_beta():
    assert is_class_ii_alpha_gene("DRB1") is False
    assert is_class_ii_alpha_gene("DQB1") is False
    assert is_class_ii_alpha_gene("DPB1") is False


# ---------------------------------------------------------------------------
# Gene-specific Cys1 positions
# ---------------------------------------------------------------------------


def test_dqa1_mature_start_with_sp():
    """DQA1 with signal peptide should still land near the annotated boundary."""
    # HLA-DQA1*01:01:01:01 full sequence (UniProt P01909)
    hla_dqa1_full = (
        "MILNKALLLGALALTTVMSPCGGEDIVADHVASCGVNLYQFYGPSGQYTHEFDGDEQFYVDLERKETAWRWPEFSKFGGFDPQGALR"
        "NIATQKHNLNIVIKRSNSTAATNEVPEVTVFSKSPVTLGQPNILICFIDKFTPPVVNVTWLRNGKPVTTGVSETVFLPREDHLFRK"
        "FHYLPFLPSTDDYDCRVEHWGLDQPLLKHWEAQEPIQMPETPENVVACLQNLMKLAQINRLNKEDPA"
    )
    result = decompose_class_ii_alpha(hla_dqa1_full, allele="HLA-DQA1*01:01", gene="DQA1")
    assert result.ok
    # The current parser does not use gene-specific absolute Cys-position
    # constants. It should still converge near the true SP = 23 from the
    # structural grammar and SP evidence.
    assert abs(result.mature_start - 23) <= 4, f"mature_start={result.mature_start}, expected ~23"
    assert result.groove1_len >= 78  # DQA α1 is ~86; enumeration may get ~82-87


def test_class_ii_decomposition_complete():
    """Verify that groove + ig_domain + tail covers the inferred mature sequence for class II."""
    # Alpha chain
    alpha = decompose_class_ii_alpha(HLA_DRA0101_MATURE, allele="HLA-DRA*01:01")
    assert alpha.ok
    reconstructed_alpha = alpha.groove1 + alpha.ig_domain + alpha.tail
    mature_alpha = HLA_DRA0101_MATURE[alpha.mature_start :].upper()
    assert reconstructed_alpha == mature_alpha

    # Beta chain (full sequence with SP)
    beta = decompose_class_ii_beta(HLA_DRB10101_FULL, allele="HLA-DRB1*01:01")
    assert beta.ok
    reconstructed_beta = beta.groove2 + beta.ig_domain + beta.tail
    inferred_mature = HLA_DRB10101_FULL[beta.mature_start :].upper()
    assert reconstructed_beta == inferred_mature
