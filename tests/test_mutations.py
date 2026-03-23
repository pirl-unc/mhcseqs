"""Tests for mutant MHC groove extraction."""

import pytest

from mhcseqs.groove import apply_mutations, decompose_class_i, decompose_domains

# HLA-A*02:01 mature sequence (UniProt P01892, 341 aa)
HLA_A0201_MATURE = (
    "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLG"
    "TLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRA"
    "YLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGD"
    "GTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSLTACKV"
)

# H-2-Kb mature sequence (UniProt P01901, SP=21 aa removed, 348 aa)
H2KB_MATURE = (
    "GPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEGPEYWERETQKAKGNEQSFRVDLRT"
    "LLGYYNQSKGGSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAERLRAYL"
    "EGTCVEWLRRYLKNGNATLLRTDSPKAHVTHHSRPEDKVTLRCWALGFYPADITLTWQLNGEELIQDMELVETRPAGDGT"
    "FQKWASVVVPLGKEQYYTCHVYHQGLPEPLTLRWEPPPSTVSNMATVAVLVVLGAAIVTGAVVAFVMKMRRRNTGGKGG"
    "DYALAPGSQTSDLSLPDCKVMVHDPHSLA"
)

# HLA-DRB1*01:01 full sequence (from UniProt P01911, including signal peptide)
HLA_DRB10101_FULL = (
    "MVCLKLPGGSCMTALTVTLMVLSSPLALA"
    "GDTRPRFLWQLKFECHFFNGTERVSFLTTTPKTWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNSTAATNEVPEV"
    "TVFSKSPVTLGQPNTLICLVDNIFPPVVNITWLSNGQSVTEGVSETSFLSKSDHSFFKISYLTFLPSADEIYDCKVEHW"
    "GLDEPLLKHWEPEIPAPMSELTETVVCALGLSVGLVGIVVGTVFIIRGLRSGSFASRGPRGREPAE"
)


# ---------------------------------------------------------------------------
# apply_mutations with string format ("K66A")
# ---------------------------------------------------------------------------


def test_apply_mutation_string_class_i():
    """HLA-A*02:01 K66A: lysine at mature pos 66 -> alanine."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    assert wt.ok
    assert wt.groove1[65] == "K"  # mature pos 66, 0-indexed = 65

    mut = apply_mutations(wt, ["K66A"])
    assert mut.groove1[65] == "A"
    assert mut.mutations == ("K66A",)
    # Rest of groove unchanged
    assert mut.groove1[:65] == wt.groove1[:65]
    assert mut.groove1[66:] == wt.groove1[66:]
    assert mut.groove2 == wt.groove2


def test_apply_mutation_d77s():
    """HLA-A*02:01 D77S: aspartate at mature pos 77 -> serine."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    assert wt.groove1[76] == "D"

    mut = apply_mutations(wt, ["D77S"])
    assert mut.groove1[76] == "S"
    assert mut.mutations == ("D77S",)


def test_apply_mutation_y84a():
    """HLA-A*02:01 Y84A: tyrosine at mature pos 84 -> alanine."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    assert wt.groove1[83] == "Y"

    mut = apply_mutations(wt, ["Y84A"])
    assert mut.groove1[83] == "A"


def test_apply_mutation_in_groove2():
    """HLA-A*02:01 Y116A: tyrosine at mature pos 116 (alpha2 domain) -> alanine."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    # Mature pos 116 is at index 115, which is in groove2 at offset 115 - 90 = 25
    assert wt.groove2[25] == "Y"

    mut = apply_mutations(wt, ["Y116A"])
    assert mut.groove2[25] == "A"
    assert mut.groove1 == wt.groove1  # groove1 unaffected


def test_apply_multiple_mutations():
    """Apply K66A and D77S simultaneously."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")

    mut = apply_mutations(wt, ["K66A", "D77S"])
    assert mut.groove1[65] == "A"  # K66A
    assert mut.groove1[76] == "S"  # D77S
    assert len(mut.mutations) == 2
    assert "K66A" in mut.mutations
    assert "D77S" in mut.mutations


# ---------------------------------------------------------------------------
# H-2-Kb mutations
# ---------------------------------------------------------------------------


def test_h2kb_k66a():
    """H-2-Kb K66A."""
    wt = decompose_class_i(H2KB_MATURE, allele="H-2-Kb")
    assert wt.ok
    assert wt.groove1[65] == "K"

    mut = apply_mutations(wt, ["K66A"])
    assert mut.groove1[65] == "A"


def test_h2kb_d77s():
    """H-2-Kb D77S: aspartate at mature pos 77 -> serine (pocket F)."""
    wt = decompose_class_i(H2KB_MATURE, allele="H-2-Kb")
    assert wt.ok
    assert wt.groove1[76] == "D"

    mut = apply_mutations(wt, ["D77S"])
    assert mut.groove1[76] == "S"


# ---------------------------------------------------------------------------
# Tuple format
# ---------------------------------------------------------------------------


def test_apply_mutation_tuple_3():
    """Mutation as (pos, original, mutant) tuple."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    mut = apply_mutations(wt, [(66, "K", "A")])
    assert mut.groove1[65] == "A"
    assert mut.mutations == ("K66A",)


def test_apply_mutation_tuple_2():
    """Mutation as (pos, mutant) tuple — no original AA check."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    mut = apply_mutations(wt, [(66, "A")])
    assert mut.groove1[65] == "A"
    assert mut.mutations == ("K66A",)  # original discovered from sequence


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


def test_wrong_original_aa_raises():
    """Providing wrong original AA should raise ValueError."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    with pytest.raises(ValueError, match="expected R.*found K"):
        apply_mutations(wt, ["R66A"])  # position 66 is K, not R


def test_out_of_range_raises():
    """Mutation position beyond mature sequence should raise ValueError."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    with pytest.raises(ValueError, match="out of range"):
        apply_mutations(wt, ["K999A"])


def test_failed_result_raises():
    """Cannot apply mutations to a failed AlleleRecord."""
    from mhcseqs.groove import AlleleRecord

    bad = AlleleRecord(status="too_short")
    with pytest.raises(ValueError, match="failed result"):
        apply_mutations(bad, ["K66A"])


def test_bad_mutation_string_raises():
    """Unparseable mutation string should raise ValueError."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    with pytest.raises(ValueError, match="Cannot parse"):
        apply_mutations(wt, ["66A"])  # missing original


# ---------------------------------------------------------------------------
# decompose_domains with mutations parameter
# ---------------------------------------------------------------------------


def test_decompose_domains_with_mutations():
    """decompose_domains should accept mutations parameter."""
    result = decompose_domains(
        HLA_A0201_MATURE,
        mhc_class="I",
        allele="HLA-A*02:01",
        mutations=["K66A"],
    )
    assert result.ok
    assert result.groove1[65] == "A"
    assert result.mutations == ("K66A",)


def test_decompose_domains_empty_mutations():
    """Empty mutations list should return wild-type result."""
    wt = decompose_domains(HLA_A0201_MATURE, mhc_class="I", allele="HLA-A*02:01")
    mut = decompose_domains(HLA_A0201_MATURE, mhc_class="I", allele="HLA-A*02:01", mutations=[])
    assert wt.groove_seq == mut.groove_seq
    assert mut.mutations == ()


# ---------------------------------------------------------------------------
# Class II mutations
# ---------------------------------------------------------------------------


def test_class_ii_beta_mutation():
    """Apply mutation to class II beta chain groove."""
    wt = decompose_domains(
        HLA_DRB10101_FULL,
        mhc_class="II",
        chain="beta",
        allele="HLA-DRB1*01:01",
        gene="DRB1",
    )
    assert wt.ok
    # Pick a known residue in the beta1 groove
    pos1_aa = wt.groove2[0]  # mature position 1 of beta chain
    mut = decompose_domains(
        HLA_DRB10101_FULL,
        mhc_class="II",
        chain="beta",
        allele="HLA-DRB1*01:01",
        gene="DRB1",
        mutations=[(1, "X")],  # mutate first residue to X
    )
    assert mut.ok
    assert mut.groove2[0] == "X"
    assert mut.mutations == (f"{pos1_aa}1X",)


# ---------------------------------------------------------------------------
# Structural integrity after mutation
# ---------------------------------------------------------------------------


def test_mutation_preserves_domain_lengths():
    """Mutation should not change domain lengths."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    mut = apply_mutations(wt, ["K66A", "D77S", "Y84A"])

    assert mut.groove1_len == wt.groove1_len
    assert mut.groove2_len == wt.groove2_len
    assert mut.ig_domain_len == wt.ig_domain_len
    assert mut.tail_len == wt.tail_len
    assert mut.mature_start == wt.mature_start
    assert mut.status == wt.status


def test_mutation_roundtrip_reconstruction():
    """Mutated domain parts should reconstruct the mutated mature sequence."""
    wt = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    mut = apply_mutations(wt, ["K66A"])

    # Reconstruct mature from parts
    reconstructed = mut.groove1 + mut.groove2 + mut.ig_domain + mut.tail
    expected = list(HLA_A0201_MATURE)
    expected[65] = "A"
    assert reconstructed == "".join(expected)


# ---------------------------------------------------------------------------
# Known IEDB mutant alleles (comprehensive spot-checks)
# ---------------------------------------------------------------------------

# These are real IEDB mutant MHC molecules. Each entry:
# (allele, sequence_source, mutations, check_pos_0indexed, expected_aa)
IEDB_SPOT_CHECKS = [
    # HLA-A*02:01 single-site mutants
    ("HLA-A*02:01", HLA_A0201_MATURE, "K66A", 65, "A"),
    ("HLA-A*02:01", HLA_A0201_MATURE, "D77S", 76, "S"),
    ("HLA-A*02:01", HLA_A0201_MATURE, "T80A", 79, "A"),
    ("HLA-A*02:01", HLA_A0201_MATURE, "Y84A", 83, "A"),
    ("HLA-A*02:01", HLA_A0201_MATURE, "Y84C", 83, "C"),
    ("HLA-A*02:01", HLA_A0201_MATURE, "W167A", 166, "A"),  # mature pos 167, 0-indexed = 166
    # H-2-Kb single-site mutants
    ("H-2-Kb", H2KB_MATURE, "K66A", 65, "A"),
    ("H-2-Kb", H2KB_MATURE, "D77S", 76, "S"),
    ("H-2-Kb", H2KB_MATURE, "Y84A", 83, "A"),
    ("H-2-Kb", H2KB_MATURE, "T80A", 79, "A"),
]


@pytest.mark.parametrize(
    "allele, seq, mutation_str, check_idx, expected_aa",
    IEDB_SPOT_CHECKS,
    ids=[f"{row[0]}_{row[2]}" for row in IEDB_SPOT_CHECKS],
)
def test_iedb_mutant(allele, seq, mutation_str, check_idx, expected_aa):
    """Verify IEDB mutant alleles produce the expected mutated groove."""
    result = decompose_domains(seq, mhc_class="I", allele=allele, mutations=[mutation_str])
    assert result.ok
    # Figure out which domain the check_idx falls in
    if check_idx < result.groove1_len:
        assert result.groove1[check_idx] == expected_aa
    else:
        g2_idx = check_idx - result.groove1_len
        assert result.groove2[g2_idx] == expected_aa
