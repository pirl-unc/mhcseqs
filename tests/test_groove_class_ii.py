from mhcseqs.groove import (
    extract_groove,
    is_class_ii_alpha_gene,
    parse_class_ii_alpha,
    parse_class_ii_beta,
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


def test_parse_class_ii_alpha_dra0101():
    result = parse_class_ii_alpha(HLA_DRA0101_MATURE, allele="HLA-DRA*01:01")
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


def test_parse_class_ii_beta_drb10101():
    result = parse_class_ii_beta(HLA_DRB10101_FULL, allele="HLA-DRB1*01:01")
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
    # β1 domain is ~93 aa
    assert 80 <= result.groove2_len <= 105


def test_parse_class_ii_alpha_too_short():
    result = parse_class_ii_alpha("ACDEFGH" * 5, allele="short")
    assert not result.ok
    assert result.status == "too_short"


def test_parse_class_ii_beta_too_short():
    result = parse_class_ii_beta("ACDEFGH" * 5, allele="short")
    assert not result.ok
    assert result.status == "too_short"


def test_parse_class_ii_alpha_fragment():
    # A short sequence should trigger fragment fallback
    fragment = HLA_DRA0101_MATURE[:85]
    result = parse_class_ii_alpha(fragment, allele="fragment")
    assert result.ok
    assert result.status == "fragment_fallback"
    assert result.groove1 == fragment.upper()


def test_parse_class_ii_beta_fragment():
    fragment = HLA_DRB10101_FULL[29:124]  # ~95 aa from within the mature region
    result = parse_class_ii_beta(fragment, allele="fragment")
    assert result.ok
    assert result.status == "fragment_fallback"
    assert result.groove2 == fragment.upper()


def test_extract_groove_class_ii_alpha():
    result = extract_groove(
        HLA_DRA0101_MATURE,
        mhc_class="II",
        chain="alpha",
        allele="HLA-DRA*01:01",
        gene="DRA",
    )
    assert result.ok
    assert result.mhc_class == "II"
    assert result.chain == "alpha"


def test_extract_groove_class_ii_beta():
    result = extract_groove(
        HLA_DRB10101_FULL,
        mhc_class="II",
        chain="beta",
        allele="HLA-DRB1*01:01",
        gene="DRB1",
    )
    assert result.ok
    assert result.mhc_class == "II"
    assert result.chain == "beta"


def test_extract_groove_class_ii_infer_chain_alpha():
    # Without explicit chain, gene name "DRA" should infer alpha
    result = extract_groove(
        HLA_DRA0101_MATURE,
        mhc_class="II",
        gene="DRA",
        allele="HLA-DRA*01:01",
    )
    assert result.ok
    assert result.chain == "alpha"


def test_extract_groove_class_ii_infer_chain_beta():
    # Without explicit chain, gene name "DRB1" should infer beta
    result = extract_groove(
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


def test_class_ii_decomposition_complete():
    """Verify that groove + ig_domain + tail covers the inferred mature sequence for class II."""
    # Alpha chain
    alpha = parse_class_ii_alpha(HLA_DRA0101_MATURE, allele="HLA-DRA*01:01")
    assert alpha.ok
    reconstructed_alpha = alpha.groove1 + alpha.ig_domain + alpha.tail
    mature_alpha = HLA_DRA0101_MATURE[alpha.mature_start :].upper()
    assert reconstructed_alpha == mature_alpha

    # Beta chain (full sequence with SP)
    beta = parse_class_ii_beta(HLA_DRB10101_FULL, allele="HLA-DRB1*01:01")
    assert beta.ok
    reconstructed_beta = beta.groove2 + beta.ig_domain + beta.tail
    inferred_mature = HLA_DRB10101_FULL[beta.mature_start :].upper()
    assert reconstructed_beta == inferred_mature
