"""Tests for IMGT G-DOMAIN numbering and mutant MHC validation."""

import pytest

from mhcseqs.domain_parsing import decompose_class_i
from mhcseqs.imgt import (
    CONSERVED_CYS_POSITIONS,
    FIXED_FRAME_ALPHABET,
    G_ALPHA1_POSITIONS,
    G_ALPHA2_POSITIONS,
    G_ALPHA_POSITIONS,
    G_BETA_POSITIONS,
    HELIX_INSERTIONS,
    HLA_A0201_G_ALPHA1,
    HLA_A0201_G_ALPHA2,
    HLA_DRA0101_G_ALPHA,
    HLA_DRB1_0101_G_BETA,
    IMGT_CLASS_I_GROOVE_SLOTS,
    IMGT_CLASS_II_GROOVE_SLOTS,
    IMGT_G_ALPHA1_SLOTS,
    IMGT_G_ALPHA2_SLOTS,
    IMGT_G_ALPHA_SLOTS,
    IMGT_G_BETA_SLOTS,
    IMGT_PAN_CLASS_GROOVE_SLOTS,
    IMGT_PAN_CLASS_HALF1_SLOTS,
    IMGT_PAN_CLASS_HALF2_SLOTS,
    IMGT_PAN_CLASS_NETMHCIIPAN_SLOTS,
    IMGT_PAN_CLASS_NETMHCPAN_SLOTS,
    IMGT_PAN_CLASS_POCKET_SLOTS,
    AlignedResidue,
    DomainAlignment,
    GrooveAlignment,
    align_domain_to_imgt,
    align_groove_to_imgt,
    imgt_to_mature,
    imgt_to_mature_class_i,
    load_class_i_profile,
    load_class_ii_profile,
    map_class_ii_groove_positions,
    map_groove_positions,
    mature_to_imgt,
    mature_to_imgt_class_i,
    structural_element,
)

# ---------------------------------------------------------------------------
# Reference sequences (from UniProt, signal peptide removed)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# G-ALPHA1 mapping (trivial 1:1)
# ---------------------------------------------------------------------------


def test_galpha1_mapping():
    domain, pos = mature_to_imgt_class_i(1)
    assert domain == "G-ALPHA1"
    assert pos == "1"


def test_galpha1_mapping_pos90():
    domain, pos = mature_to_imgt_class_i(90)
    assert domain == "G-ALPHA1"
    assert pos == "90"


def test_galpha1_roundtrip():
    for i in range(1, 91):
        domain, imgt = mature_to_imgt_class_i(i)
        assert domain == "G-ALPHA1"
        assert imgt_to_mature_class_i("G-ALPHA1", imgt) == i


# ---------------------------------------------------------------------------
# G-ALPHA2 mapping
# ---------------------------------------------------------------------------


def test_galpha2_starts_at_91():
    domain, pos = mature_to_imgt_class_i(91)
    assert domain == "G-ALPHA2"
    assert pos == "1"


def test_galpha2_conserved_cys_11():
    """IMGT CYS-11 in G-ALPHA2 = mature position 101 (conserved disulfide)."""
    assert imgt_to_mature_class_i("G-ALPHA2", "11") == 101
    # Verify it's actually a Cys in HLA-A*02:01
    assert HLA_A0201_MATURE[100] == "C"


def test_galpha2_conserved_cys_74():
    """IMGT CYS-74 in G-ALPHA2 = mature position 164 (conserved disulfide)."""
    assert imgt_to_mature_class_i("G-ALPHA2", "74") == 164
    assert HLA_A0201_MATURE[163] == "C"


def test_galpha2_gap_positions():
    """IMGT positions 40 and 41 are gaps in G-ALPHA2 (CD-TURN)."""
    assert imgt_to_mature_class_i("G-ALPHA2", "40") is None
    assert imgt_to_mature_class_i("G-ALPHA2", "41") is None


def test_galpha2_insertion_61a():
    """IMGT 61A is a helix insertion mapping to mature 150."""
    domain, pos = mature_to_imgt_class_i(150)
    assert domain == "G-ALPHA2"
    assert pos == "61A"
    assert imgt_to_mature_class_i("G-ALPHA2", "61A") == 150


def test_galpha2_insertion_72a():
    """IMGT 72A is a helix insertion mapping to mature 162."""
    domain, pos = mature_to_imgt_class_i(162)
    assert domain == "G-ALPHA2"
    assert pos == "72A"
    assert imgt_to_mature_class_i("G-ALPHA2", "72A") == 162


def test_galpha2_last_position():
    """IMGT 92 in G-ALPHA2 = mature 182 (last position of the domain)."""
    domain, pos = mature_to_imgt_class_i(182)
    assert domain == "G-ALPHA2"
    assert pos == "92"


def test_galpha2_roundtrip():
    """All non-gap G-ALPHA2 positions round-trip correctly."""
    for imgt_pos in G_ALPHA2_POSITIONS:
        mature = imgt_to_mature_class_i("G-ALPHA2", imgt_pos)
        assert mature is not None
        domain, back = mature_to_imgt_class_i(mature)
        assert domain == "G-ALPHA2"
        assert back == imgt_pos


def test_galpha2_total_positions():
    """G-ALPHA2 should have 92 occupied positions (for HLA-A*02:01)."""
    assert len(G_ALPHA2_POSITIONS) == 92


# ---------------------------------------------------------------------------
# C-LIKE domain
# ---------------------------------------------------------------------------


def test_clike_starts_at_183():
    domain, pos = mature_to_imgt_class_i(183)
    assert domain == "C-LIKE"
    assert pos == "1"


# ---------------------------------------------------------------------------
# Structural element lookup
# ---------------------------------------------------------------------------


def test_structural_element_strand():
    assert structural_element("1") == "A-STRAND"
    assert structural_element("14") == "A-STRAND"
    assert structural_element("18") == "B-STRAND"
    assert structural_element("31") == "C-STRAND"
    assert structural_element("42") == "D-STRAND"


def test_structural_element_turns():
    assert structural_element("15") == "AB-TURN"
    assert structural_element("29") == "BC-TURN"
    assert structural_element("39") == "CD-TURN"


def test_structural_element_helix():
    assert structural_element("50") == "HELIX"
    assert structural_element("61A") == "HELIX"
    assert structural_element("72A") == "HELIX"
    assert structural_element("92") == "HELIX"


def test_structural_element_insertions():
    assert structural_element("1.2") == "A-STRAND"
    assert structural_element("1.1") == "A-STRAND"
    assert structural_element("49.1") == "D-STRAND"


# ---------------------------------------------------------------------------
# Class II approximate mapping
# ---------------------------------------------------------------------------


def test_class_ii_alpha():
    domain, pos = mature_to_imgt(10, mhc_class="II", chain="alpha")
    assert domain == "G-ALPHA"
    assert pos == "10"


def test_class_ii_beta():
    domain, pos = mature_to_imgt(50, mhc_class="II", chain="beta")
    assert domain == "G-BETA"
    assert pos == "50"


def test_class_ii_requires_chain():
    with pytest.raises(ValueError, match="chain"):
        mature_to_imgt(10, mhc_class="II")


# ---------------------------------------------------------------------------
# High-level mature_to_imgt dispatching
# ---------------------------------------------------------------------------


def test_mature_to_imgt_class_i():
    domain, pos = mature_to_imgt(66, mhc_class="I")
    assert domain == "G-ALPHA1"
    assert pos == "66"


def test_mature_to_imgt_invalid():
    with pytest.raises(ValueError):
        mature_to_imgt(0, mhc_class="I")


# ---------------------------------------------------------------------------
# Validate known mutant MHC reference residues
# ---------------------------------------------------------------------------
# These are well-characterized MHC mutations from IEDB and crystallography
# literature.  The residue positions use mature protein numbering (1-indexed).
# We verify that:
#   1) The IMGT position mapping is correct
#   2) The reference residue at that position matches the actual sequence

# (allele, mature_pos, expected_aa, imgt_domain, imgt_pos, description)
KNOWN_CLASS_I_MUTANTS = [
    # --- HLA-A*02:01 groove pocket residues ---
    ("HLA-A*02:01", 9, "F", "G-ALPHA1", "9", "pocket B floor"),
    ("HLA-A*02:01", 24, "A", "G-ALPHA1", "24", "pocket A"),
    ("HLA-A*02:01", 45, "M", "G-ALPHA1", "45", "pocket B"),
    ("HLA-A*02:01", 59, "Y", "G-ALPHA1", "59", "pocket A"),
    ("HLA-A*02:01", 63, "E", "G-ALPHA1", "63", "pocket A/B"),
    ("HLA-A*02:01", 66, "K", "G-ALPHA1", "66", "pocket A/B, K66A mutant"),
    ("HLA-A*02:01", 67, "V", "G-ALPHA1", "67", "pocket B"),
    ("HLA-A*02:01", 70, "H", "G-ALPHA1", "70", "pocket A/B"),
    ("HLA-A*02:01", 74, "H", "G-ALPHA1", "74", "pocket A/helix"),
    ("HLA-A*02:01", 77, "D", "G-ALPHA1", "77", "pocket F, D77S mutant"),
    ("HLA-A*02:01", 80, "T", "G-ALPHA1", "80", "pocket F, T80A mutant"),
    ("HLA-A*02:01", 81, "L", "G-ALPHA1", "81", "pocket F"),
    ("HLA-A*02:01", 84, "Y", "G-ALPHA1", "84", "pocket F floor, Y84A mutant"),
    # --- HLA-A*02:01 alpha-2 domain ---
    ("HLA-A*02:01", 95, "V", "G-ALPHA2", "5", "pocket B floor"),
    ("HLA-A*02:01", 97, "R", "G-ALPHA2", "7", "pocket B/C"),
    ("HLA-A*02:01", 99, "Y", "G-ALPHA2", "9", "pocket B/C floor"),
    ("HLA-A*02:01", 114, "H", "G-ALPHA2", "24", "pocket B"),
    ("HLA-A*02:01", 116, "Y", "G-ALPHA2", "26", "pocket F specificity"),
    ("HLA-A*02:01", 143, "T", "G-ALPHA2", "55", "pocket D/E"),
    ("HLA-A*02:01", 146, "K", "G-ALPHA2", "58", "pocket D, K146 mutant"),
    ("HLA-A*02:01", 147, "W", "G-ALPHA2", "59", "pocket D"),
    ("HLA-A*02:01", 150, "A", "G-ALPHA2", "61A", "helix insertion, pocket E"),
    ("HLA-A*02:01", 152, "V", "G-ALPHA2", "63", "TCR contact"),
    ("HLA-A*02:01", 155, "Q", "G-ALPHA2", "66", "TCR contact"),
    ("HLA-A*02:01", 156, "L", "G-ALPHA2", "67", "TCR contact"),
    ("HLA-A*02:01", 159, "Y", "G-ALPHA2", "70", "pocket C/E"),
    ("HLA-A*02:01", 163, "T", "G-ALPHA2", "73", "pocket, T163A mutant"),
    ("HLA-A*02:01", 167, "W", "G-ALPHA2", "77", "TCR/peptide, W167A mutant"),
    ("HLA-A*02:01", 170, "R", "G-ALPHA2", "80", "TCR contact"),
    ("HLA-A*02:01", 171, "Y", "G-ALPHA2", "81", "TCR contact"),
    # --- H-2-Kb groove residues ---
    ("H-2-Kb", 66, "K", "G-ALPHA1", "66", "pocket A/B, K66A mutant"),
    ("H-2-Kb", 77, "D", "G-ALPHA1", "77", "pocket F, D77S mutant"),
    ("H-2-Kb", 80, "T", "G-ALPHA1", "80", "pocket F"),
    ("H-2-Kb", 84, "Y", "G-ALPHA1", "84", "pocket F floor, Y84A mutant"),
    ("H-2-Kb", 152, "E", "G-ALPHA2", "63", "TCR contact"),
    ("H-2-Kb", 167, "W", "G-ALPHA2", "77", "TCR/peptide contact"),
]


_SEQUENCES = {
    "HLA-A*02:01": HLA_A0201_MATURE,
    "H-2-Kb": H2KB_MATURE,
}


@pytest.mark.parametrize(
    "allele, mature_pos, expected_aa, expected_domain, expected_imgt, desc",
    KNOWN_CLASS_I_MUTANTS,
    ids=[f"{row[0]}_{row[4]}_{row[2]}{row[1]}" for row in KNOWN_CLASS_I_MUTANTS],
)
def test_known_mutant_reference_residue(allele, mature_pos, expected_aa, expected_domain, expected_imgt, desc):
    """Verify that mhcseqs's reference sequence has the correct residue at each
    known mutant position, and that the IMGT mapping is consistent."""
    # Check the IMGT mapping
    domain, imgt_pos = mature_to_imgt_class_i(mature_pos)
    assert domain == expected_domain, f"{desc}: expected domain {expected_domain}, got {domain}"
    assert imgt_pos == expected_imgt, f"{desc}: expected IMGT {expected_imgt}, got {imgt_pos}"

    # Check the reference residue in the actual sequence
    seq = _SEQUENCES[allele]
    actual_aa = seq[mature_pos - 1]  # 0-indexed
    assert actual_aa == expected_aa, f"{allele} position {mature_pos} ({desc}): expected {expected_aa}, got {actual_aa}"


# ---------------------------------------------------------------------------
# Additional IMGT landmark checks
# ---------------------------------------------------------------------------


def test_conserved_cys_positions():
    """IMGT positions 11 and 74 are the conserved disulfide Cys pair."""
    assert CONSERVED_CYS_POSITIONS == (11, 74)
    # In G-ALPHA2, these map to mature 101 and 164
    assert imgt_to_mature("G-ALPHA2", "11") == 101
    assert imgt_to_mature("G-ALPHA2", "74") == 164
    # Both are Cys in HLA-A*02:01
    assert HLA_A0201_MATURE[100] == "C"
    assert HLA_A0201_MATURE[163] == "C"
    # Both are Cys in H-2-Kb
    assert H2KB_MATURE[100] == "C"
    assert H2KB_MATURE[163] == "C"


def test_helix_insertions_defined():
    assert "61A" in HELIX_INSERTIONS
    assert "72A" in HELIX_INSERTIONS
    assert "61B" in HELIX_INSERTIONS


def test_structural_element_general_dot_insertions():
    assert structural_element("51.1") == "HELIX"
    assert structural_element("51.3") == "HELIX"
    assert structural_element("20.1") == "B-STRAND"
    assert structural_element("61A.1") == "HELIX"


# ---------------------------------------------------------------------------
# Reference-anchored Class I IMGT G-DOMAIN alignment (#120)
# ---------------------------------------------------------------------------

_HLA_A01103_89AA_G1 = "SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSED"
_HLA_A01103_G2 = "GSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRCVDGLRRYLENGKETLQRTD"

_PTAL_N0101_G1 = "GFHSLRYFYTAWSRPGSGEPRFVAVGYVDDTQFVRFDSDNASPRAEPRAPWMDLVEQQDPQYWDRNTRNARDAAQTYRVGLDNVRGYYNQSEA"
_PTAL_N0101_G2 = "GSHTIQRMYGCDVGPHGRLLRGYDQLAYDGADYIALNEDLRSWTAADLAAQNTRRKWEEAGYAERDRAYLEGECVEWLLKHLENGRETLLRAD"

_TRVU_UB0101_G1 = "GSHSLRYFYTAVSGPELREPRFLSVGYVDEQQFVRFDSASESPREEPRAKWIERVGEEDPEYWERQTGILRRNTQVFRVGLETLRGYFNQSAG"
_TRVU_UB0101_G2 = "GVHTLQTMYGCELTPELTFTRGFDQSAYDGRDYISLDTDTYTWTATAPQAVNTKRKWEADRSIAEGWKAYLEETCVLWLKKYLEMGKDTLGRTD"

_MAFA_F0201_G1 = "SHSLRYFSTAVSRPGRREPRYRYIAVEYVDDTQFLRFDSDAAIPRMEPRAPWVEQEGPQYWERTTGYAKANARTDRVALRKLLLRYNQSEA"
_MAFA_F0201_G2 = "GSHTLQGMNGCDMGPDGRLLRGYHQHAYDGKDYISLNEDPRSWTAADTVARITQRFYEAEEYAEEFRTYLEGECLELLRRYLENGKETLQRA"

_ANPL_UAA01_G1 = "EPHSLRYFYTAVSDPSPGVPQFVTVGSVDGEVFVRYDSETRKMEPRVDWIVANVDQQYWDRETETSRGNEQIFRVNLDTARERYNQSR"
_ANPL_UAA01_G2 = "GSHTWQCMHGCDLEDGSIRGFQQCGYDGKDFIALDKDTLTYTAADAAAQITKRKWEQEGTVAEQWKNYLENTCIEWLRKYVSYGKDVLERRE"

_GAGA_BF1_G1 = "ELHTLRYISTAMTDPGPGQPWYVDVGYVDGELFTHYNSTARRAVPRTEWIAANTDQQYWDRETQIAQGNEQIDRENLDIRQQRHNQTG"
_GAGA_BF1_G2 = "GSHTAQWMYGCDILEDGTIRGYHQMACDGRDFIALAEDMKTFTAAVPEAVPTKRKWEEGGYAERKKQYLEETCVEGLRRYVEYGKAELGRRE"

_ONMY_UBA_G1 = "GLLHIASAVTHSLRYFYTATTGIPDFPEFVDVGVVNGKVISYYDSITKRKVPKQSWMEENFNQQYWNQGTDLLKGTEQTFKANIQVAQTRFNQTG"
_ONMY_UBA_G2 = "GVHIVQVMYGCTWDDDSGVIEGVQQHGYDGEDFLVYDMKTFTWIAPKQQAEITQRKWNNNPAGMAQRKSYITQECVEWLKKYVDFGKNTLMRTV"


def test_align_canonical_hla_a0201_groove():
    rec = decompose_class_i(HLA_A0201_MATURE, allele="HLA-A*02:01")
    aln = align_groove_to_imgt(rec)
    assert isinstance(aln, GrooveAlignment)
    assert aln.alpha1 is not None
    assert aln.alpha2 is not None
    assert not aln.has_indels
    assert not aln.is_truncated
    assert aln.summary() == "canonical"
    assert aln.alpha1.identity == 1.0
    assert aln.alpha2.identity == 1.0
    # mhcseqs groove2 includes mature pos 183 ('D'), which is C-terminal overhang past IMGT 92
    assert aln.alpha2.c_term_overhang == "D"
    assert aln.residue_at(66) == "K"
    assert aln.residue_at(101) == "C"
    assert aln.residue_at(164) == "C"
    assert aln.residue_at_imgt("G-ALPHA2", "11") == "C"
    assert aln.residue_at_imgt("G-ALPHA2", "74") == "C"
    assert aln.residue_at_imgt("G-ALPHA2", "61A") == "A"
    assert len(G_ALPHA1_POSITIONS) == 90


@pytest.mark.parametrize(
    ("seq", "domain"),
    [
        (HLA_A0201_G_ALPHA1, "G-ALPHA1"),
        (HLA_A0201_G_ALPHA2 + "D", "G-ALPHA2"),
        # SusSp-5*02:01 in v2.9.0: exact DP finds an internal del@50 + K ins@51.
        ("GPHSLSYLFTSVSRPGRGEPRFIAVGCVDDTPLVRFDSDAPKLRKEPRAWKIKQEGPEYWDEETRISKDNAQTLRGKLNTLRGYYNQSEA", "G-ALPHA1"),
        # Mafa-I*01:62 in v2.9.0: compensating indels in a 90-aa groove.
        ("GSHSLRYFGTTVSRPGRGEPRFIAVGYVDDTQFVRFDSDAESPRMEPRAPWMEQEGPEYWEEETRIAKARAQTDLGTLRTALRYYNQSEG", "G-ALPHA1"),
        # HLA-C*07:819Q (IPD-IMGT/HLA HLA25172): 93-aa groove with an internal insertion.
        ("GSHTLQRMYGCDLGPDGRLLRGYDQSAYDGKDYIALNEDLRSWTAADTAAQITQRKLEAARAAEQLRAYLEGTCVEWLRRYLEGTCVEWLRRY", "G-ALPHA2"),
    ],
)
def test_align_fast_path_matches_full_dp(seq, domain):
    expected = align_domain_to_imgt(seq, domain, use_fast_path=False)
    first = align_domain_to_imgt(seq, domain, use_fast_path=True)
    assert first == expected
    # Cached score/traceback must return fresh mutable coordinate dicts per call.
    first.imgt_to_aa.clear()
    first.fixed_slot_to_aa.clear()
    assert align_domain_to_imgt(seq, domain, use_fast_path=True) == expected


def test_align_partial_exon2_89aa_distinguishes_truncation_from_deletion():
    aln = map_groove_positions(_HLA_A01103_89AA_G1, _HLA_A01103_G2)
    assert aln.alpha1 is not None
    assert aln.alpha1.n_term_missing == 1
    assert aln.alpha1.c_term_missing == 0
    assert aln.alpha1.canonical_start_mature_pos == 2
    assert aln.alpha1.deletions == ()
    assert aln.alpha1.insertions == ()
    assert aln.is_truncated
    assert not aln.has_indels
    assert aln.residue_at(1) is None
    assert aln.residue_at(2) == "S"
    assert aln.residue_at(66) == "N"  # HLA-A*01 has N66-M67 (vs HLA-A*02 K66-V67)
    assert aln.residue_at(67) == "M"
    assert aln.residue_at(101) == "C"
    assert aln.extract_positions([1, 2, 66, 101]) == "XSNC"


def test_align_ptal_n0101_mdl_insertion():
    aln = map_groove_positions(_PTAL_N0101_G1, _PTAL_N0101_G2)
    assert aln.has_insertions
    assert not aln.has_deletions
    assert len(aln.insertions) == 1
    ins = aln.insertions[0]
    assert ins.domain == "G-ALPHA1"
    assert ins.after_mature_pos == 51
    assert ins.after_imgt_pos == "51"
    assert ins.sequence == "MDL"
    assert ins.length == 3
    assert ins.imgt_codes == ("51.1", "51.2", "51.3")
    assert ins.structural_element == "HELIX"
    # Canonical positions after the insertion still map accurately
    assert aln.residue_at(51) == "W"
    assert aln.residue_at(59) == "Y"
    assert aln.residue_at(60) == "W"
    assert aln.residue_at(84) == "Y"
    assert aln.residue_at(101) == "C"
    assert aln.residue_at(164) == "C"
    assert "MDL" in aln.summary()


def test_align_trvu_ub0101_alpha1_and_alpha2_insertions():
    aln = map_groove_positions(_TRVU_UB0101_G1, _TRVU_UB0101_G2)
    assert aln.alpha1 is not None and aln.alpha2 is not None
    assert sum(i.length for i in aln.alpha1.insertions) == 3
    assert sum(i.length for i in aln.alpha2.insertions) == 1
    assert aln.residue_at(101) == "C"
    assert aln.residue_at(164) == "C"
    assert aln.residue_at(59) == "Y"


def test_align_mafa_f0201_simultaneous_truncation_and_insertion():
    aln = map_groove_positions(_MAFA_F0201_G1, _MAFA_F0201_G2)
    assert aln.alpha1 is not None
    assert aln.alpha1.n_term_missing == 1
    assert sum(i.length for i in aln.alpha1.insertions) == 2
    assert aln.alpha1.deletions == ()
    assert aln.residue_at(1) is None
    assert aln.residue_at(2) == "S"
    assert aln.residue_at(59) == "Y"


def test_align_avian_and_fish_indels_and_overhang():
    for g1, g2 in ((_ANPL_UAA01_G1, _ANPL_UAA01_G2), (_GAGA_BF1_G1, _GAGA_BF1_G2)):
        aln = map_groove_positions(g1, g2)
        assert aln.has_deletions
        assert aln.residue_at(101) == "C"
        assert aln.residue_at(164) == "C"

    onmy = map_groove_positions(_ONMY_UBA_G1, _ONMY_UBA_G2)
    assert onmy.alpha1 is not None and onmy.alpha2 is not None
    assert len(onmy.alpha1.n_term_overhang) >= 6
    assert onmy.alpha2.has_insertions
    assert onmy.residue_at(101) == "C"
    assert onmy.residue_at(164) == "C"


def test_align_groove_to_imgt_accepts_dict_record_and_strings():
    row = {"groove1": HLA_A0201_G_ALPHA1, "groove2": HLA_A0201_G_ALPHA2}
    aln_dict = align_groove_to_imgt(row)
    aln_str = align_groove_to_imgt(HLA_A0201_G_ALPHA1, HLA_A0201_G_ALPHA2)
    assert aln_dict.mature_to_aa == aln_str.mature_to_aa

    empty_dom = align_domain_to_imgt("", "G-ALPHA1")
    assert isinstance(empty_dom, DomainAlignment)
    assert empty_dom.n_term_missing == 90

    with pytest.raises(TypeError):
        align_groove_to_imgt(12345)


# ---------------------------------------------------------------------------
# Multi-clade structural PSSM profile, fixed-index frame, a3m_seq, and residues
# ---------------------------------------------------------------------------


def test_multiclade_profile_loaded_with_loop_hotspots():
    profile = load_class_i_profile()
    assert profile["version"] == "1"
    g1 = profile["domains"]["G-ALPHA1"]
    g2 = profile["domains"]["G-ALPHA2"]
    assert g1["length"] == 90
    assert g2["length"] == 92
    # Rigid strand/helix gap penalties must exceed loop hotspot penalties
    assert g1["gap_open_ins"][50] > g1["gap_open_ins"][25]  # -7.0 > -14.0 (less negative)
    assert g1["gap_open_ins"][53] > g1["gap_open_ins"][65]
    assert g1["gap_open_del"][40] > g1["gap_open_del"][25]
    assert g2["gap_open_ins"][38] > g2["gap_open_ins"][25]
    assert g2["gap_open_ins"][59] > g2["gap_open_ins"][75]


def test_fixed_frame_schema_and_hla_a0201_gap_vs_missing():
    assert len(IMGT_G_ALPHA1_SLOTS) == 96
    assert len(IMGT_G_ALPHA2_SLOTS) == 102
    assert len(IMGT_CLASS_I_GROOVE_SLOTS) == 198

    aln = map_groove_positions(HLA_A0201_G_ALPHA1, HLA_A0201_G_ALPHA2 + "D")
    assert len(aln.fixed_seq) == 198
    assert set(aln.fixed_seq) <= set(FIXED_FRAME_ALPHABET)
    # Canonical HLA-A*02:01 has '-' (structural gap), NOT 'X', in unoccupied insertion slots
    assert aln.extract_imgt_slots([("G-ALPHA1", "54A"), ("G-ALPHA1", "54B"), ("G-ALPHA1", "54C")]) == "---"
    assert aln.extract_imgt_slots([("G-ALPHA1", "19A"), ("G-ALPHA1", "19B"), ("G-ALPHA1", "19C")]) == "---"
    assert aln.extract_imgt_slots([("G-ALPHA2", "40"), ("G-ALPHA2", "41"), ("G-ALPHA2", "61B")]) == "---"
    assert "X" not in aln.fixed_seq
    assert not aln.has_overflow_insertions


def test_fixed_frame_ptal_and_trvu_preserve_backbone_and_order():
    ptal = map_groove_positions(_PTAL_N0101_G1, _PTAL_N0101_G2)
    trvu = map_groove_positions(_TRVU_UB0101_G1, _TRVU_UB0101_G2)
    hla = map_groove_positions(HLA_A0201_G_ALPHA1, HLA_A0201_G_ALPHA2)

    assert len(ptal.fixed_seq) == len(trvu.fixed_seq) == len(hla.fixed_seq) == 198
    # Trvu-UB*01:01 alpha1 inserts VGE right after IMGT 54 -> populates 54A-54C without overflow;
    # alpha2 inserts 'S' after IMGT 62 -> stays in overflow + insertion_counts['62']=1 instead of scrambling before 62='R'.
    assert trvu.extract_imgt_slots([("G-ALPHA1", "54A"), ("G-ALPHA1", "54B"), ("G-ALPHA1", "54C")]) == "VGE"
    assert trvu.alpha1 is not None and not trvu.alpha1.has_overflow_insertions
    assert trvu.alpha2 is not None and trvu.alpha2.has_overflow_insertions
    idx_h2_62 = IMGT_PAN_CLASS_GROOVE_SLOTS.index(("HALF2", "62"))
    assert trvu.pan_class_insertion_counts[idx_h2_62] == 1
    # Ptal-N*01:01 inserts MDL after W51 (before V52/E53/Q54): stays in overflow + insertion_counts['51']=3
    # rather than teleporting MDL past VEQ into 54A-54C (issue #177 monotonic order invariant).
    assert ptal.extract_imgt_slots([("G-ALPHA1", "54A"), ("G-ALPHA1", "54B"), ("G-ALPHA1", "54C")]) == "---"
    assert ptal.has_overflow_insertions
    idx_51 = IMGT_PAN_CLASS_GROOVE_SLOTS.index(("HALF1", "51"))
    assert ptal.pan_class_insertion_counts[idx_51] == 3
    assert hla.extract_imgt_slots([("G-ALPHA1", "54A"), ("G-ALPHA1", "54B"), ("G-ALPHA1", "54C")]) == "---"

    # Downstream canonical indices in fixed_seq stay perfectly homologous across all three species
    idx_y59 = IMGT_CLASS_I_GROOVE_SLOTS.index(("G-ALPHA1", "59"))
    idx_w60 = IMGT_CLASS_I_GROOVE_SLOTS.index(("G-ALPHA1", "60"))
    idx_c101 = IMGT_CLASS_I_GROOVE_SLOTS.index(("G-ALPHA2", "11"))
    idx_c164 = IMGT_CLASS_I_GROOVE_SLOTS.index(("G-ALPHA2", "74"))
    for a in (hla, ptal, trvu):
        assert a.fixed_seq[idx_y59] == "Y"
        assert a.fixed_seq[idx_w60] == "W"
        assert a.fixed_seq[idx_c101] == "C"
        assert a.fixed_seq[idx_c164] == "C"


def test_fixed_frame_separates_gap_from_unknown_on_truncations_and_single_domain():
    # 89-aa partial exon-2: missing pos 1 is 'X', while unoccupied 54A-54C is '-'
    partial = map_groove_positions(_HLA_A01103_89AA_G1, _HLA_A01103_G2)
    assert partial.fixed_seq[0] == "X"
    assert partial.fixed_seq[1] == "S"
    assert partial.extract_imgt_slots([("G-ALPHA1", "1"), ("G-ALPHA1", "2"), ("G-ALPHA1", "54A")]) == "XS-"

    # Alpha1-only fragment: entire unsequenced G-ALPHA2 half is 'X' * 102 in fixed_seq, '-' * 92 in a3m_seq
    a1_only = map_groove_positions(HLA_A0201_G_ALPHA1, "")
    assert len(a1_only.fixed_seq) == 198
    assert a1_only.fixed_seq[96:] == "X" * 102
    assert a1_only.a3m_seq[90:] == "-" * 92

    # Avian internal deletion renders as '-' (structural gap), distinct from C-term truncation 'X'
    anpl = map_groove_positions(_ANPL_UAA01_G1, _ANPL_UAA01_G2)
    assert "-" in anpl.alpha1.fixed_seq
    assert anpl.alpha1.fixed_slot_to_aa["90"] == "X"  # C-term missing pos 90


@pytest.mark.parametrize(
    ("groove1", "groove2"),
    [
        ("", ""),
        ("", HLA_A0201_G_ALPHA2 + "D"),
        (_HLA_A01103_89AA_G1, ""),
        (_GAGA_BF1_G1, _PTAL_N0101_G2),
    ],
)
def test_a3m_gaps_preserve_match_columns_and_query_sequence(groove1, groove2):
    a3m = map_groove_positions(groove1, groove2).a3m_seq
    assert all(ch == "-" or (ch.isascii() and ch.isalpha()) for ch in a3m)
    assert sum(ch == "-" or ch.isupper() for ch in a3m) == len(G_ALPHA1_POSITIONS) + len(G_ALPHA2_POSITIONS)
    assert a3m.replace("-", "").upper() == groove1 + groove2


def test_empty_domain_a3m_preserves_all_missing_columns():
    alignment = align_domain_to_imgt("", "G-ALPHA1")
    assert alignment.a3m_seq == "-" * len(G_ALPHA1_POSITIONS)
    assert alignment.n_term_missing == len(G_ALPHA1_POSITIONS)
    assert alignment.deletions == ()


def test_a3m_seq_and_residues_are_lossless_variable_length():
    # Ptal-N*01:01 alpha1: 93 actual residues, 3 lowercase 'mdl' in a3m_seq
    ptal = map_groove_positions(_PTAL_N0101_G1, _PTAL_N0101_G2)
    assert ptal.alpha1 is not None
    assert len(ptal.alpha1.residues) == 93
    assert all(isinstance(r, AlignedResidue) for r in ptal.alpha1.residues)
    assert "".join(r.aa for r in ptal.alpha1.residues) == _PTAL_N0101_G1
    ins_tokens = [r for r in ptal.alpha1.residues if r.is_insertion]
    assert [(r.aa, r.imgt_pos, r.insertion_rank, r.structural_element) for r in ins_tokens] == [
        ("M", "51.1", 1, "HELIX"),
        ("D", "51.2", 2, "HELIX"),
        ("L", "51.3", 3, "HELIX"),
    ]
    assert "PRAPWmdlVEQQ" in ptal.alpha1.a3m_seq
    # Round-trip recovery from a3m_seq (keep letters, uppercase) == exact query groove
    recovered = "".join(ch.upper() for ch in ptal.alpha1.a3m_seq if ch.isalpha())
    assert recovered == _PTAL_N0101_G1

    # 89-aa partial exon-2: 89 actual residues in alpha1, leading '-' in a3m_seq
    partial = map_groove_positions(_HLA_A01103_89AA_G1, _HLA_A01103_G2)
    assert len(partial.alpha1.residues) == 89
    assert partial.alpha1.a3m_seq.startswith("-SHSMRYF")
    assert "".join(ch.upper() for ch in partial.alpha1.a3m_seq if ch.isalpha()) == _HLA_A01103_89AA_G1

    # Gaga-BF1: 88 actual residues in alpha1, '--' for CD-turn deletion in a3m_seq, '-' for C-term truncation
    gaga = map_groove_positions(_GAGA_BF1_G1, _GAGA_BF1_G2)
    assert len(gaga.alpha1.residues) == 88
    assert "".join(r.aa for r in gaga.alpha1.residues) == _GAGA_BF1_G1
    assert "--" in gaga.alpha1.a3m_seq
    assert gaga.alpha1.a3m_seq.endswith("-")
    assert "".join(ch.upper() for ch in gaga.alpha1.a3m_seq if ch.isalpha()) == _GAGA_BF1_G1


def test_terminal_overhangs_and_observed_x_roundtrip_losslessly():
    # C-terminal overhang (standard 93-aa HLA-A*02:01 G-ALPHA2 groove ending in D)
    seq2 = HLA_A0201_G_ALPHA2 + "D"
    a2 = align_domain_to_imgt(seq2, "G-ALPHA2")
    assert a2.c_term_overhang == "D"
    assert len(a2.residues) == len(seq2)
    assert "".join(r.aa for r in a2.residues) == seq2
    assert a2.residues[-1] == AlignedResidue(
        aa="D",
        domain="G-ALPHA2",
        query_index=92,
        mature_pos=None,
        imgt_pos="92.1",
        is_insertion=True,
        insertion_rank=1,
        structural_element="C-OVERHANG",
    )
    assert a2.a3m_seq == HLA_A0201_G_ALPHA2 + "d"
    assert "".join(ch.upper() for ch in a2.a3m_seq if ch.isalpha()) == seq2

    # N-terminal overhang (rainbow trout Onmy-UBA alpha1)
    a1_onmy = align_domain_to_imgt(_ONMY_UBA_G1, "G-ALPHA1")
    assert len(a1_onmy.n_term_overhang) >= 6
    assert "".join(r.aa for r in a1_onmy.residues) == _ONMY_UBA_G1
    assert "".join(ch.upper() for ch in a1_onmy.a3m_seq if ch.isalpha()) == _ONMY_UBA_G1
    assert a1_onmy.residues[0].imgt_pos == "0.1"
    assert a1_onmy.residues[0].mature_pos is None
    assert a1_onmy.residues[0].structural_element == "N-OVERHANG"
    assert a1_onmy.a3m_seq.startswith(a1_onmy.n_term_overhang.lower())

    # Observed unknown residue 'X' vs terminal truncation '-' must stay distinct in a3m_seq
    obs_x = align_domain_to_imgt("X" + HLA_A0201_G_ALPHA1[1:], "G-ALPHA1")
    trunc = align_domain_to_imgt(HLA_A0201_G_ALPHA1[1:], "G-ALPHA1")
    assert obs_x.a3m_seq != trunc.a3m_seq
    assert obs_x.a3m_seq.startswith("XSHSMRYF")
    assert trunc.a3m_seq.startswith("-SHSMRYF")
    assert "".join(ch.upper() for ch in obs_x.a3m_seq if ch.isalpha()) == "X" + HLA_A0201_G_ALPHA1[1:]
    assert "".join(ch.upper() for ch in trunc.a3m_seq if ch.isalpha()) == HLA_A0201_G_ALPHA1[1:]
    assert "".join(r.aa for r in obs_x.residues) == "X" + HLA_A0201_G_ALPHA1[1:]
    assert "".join(r.aa for r in trunc.residues) == HLA_A0201_G_ALPHA1[1:]


def test_domain_alignment_v29_positional_constructor_preserved():
    # Positional 18-argument constructor from v2.9.0 must bind n_term_missing..identity identically.
    dom = DomainAlignment(
        "G-ALPHA1",
        "SHS",
        "GSHS",
        "GSHS",
        "-SHS",
        {2: "S", 3: "H", 4: "S"},
        {"2": "S", "3": "H", "4": "S"},
        {2: 0, 3: 1, 4: 2},
        {0: 2, 1: 3, 2: 4},
        {0: "2", 1: "3", 2: "4"},
        (),
        (),
        1,
        2,
        "A",
        "D",
        42.5,
        0.95,
    )
    assert dom.n_term_missing == 1
    assert dom.c_term_missing == 2
    assert dom.n_term_overhang == "A"
    assert dom.c_term_overhang == "D"
    assert dom.score == 42.5
    assert dom.identity == 0.95
    assert dom.residues == ()
    assert GrooveAlignment(alpha1=dom).residues == ()


def test_multiclade_profile_verified_clades_and_training_metadata():
    import csv
    from pathlib import Path

    profile = load_class_i_profile()
    src = profile["training_source"]
    assert src["release"] == "v2.9.0"
    assert len(src["sha256"]) == 64
    rep_ids = profile["representative_ids"]
    counts = profile["clade_representatives"]
    assert set(rep_ids) == set(counts)
    for clade, ids in rep_ids.items():
        assert len(ids) == counts[clade]

    diverse_csv = Path(__file__).resolve().parent.parent / "mhcseqs" / "diverse_mhc_sequences.csv"
    with open(diverse_csv, "r", encoding="utf-8") as f:
        acc_to_group = {r["uniprot_accession"].strip(): r["source_group"].strip() for r in csv.DictReader(f)}

    # Reptile clade must contain only curated reptile sequences (no misclassified amphibians/newts)
    for acc in rep_ids["reptile"]:
        assert acc_to_group.get(acc, "").startswith("reptile_"), f"Unexpected non-reptile in reptile clade: {acc}"
    # Shark/ray clade must contain only curated shark_ray sequences (no paddlefish/bony fish)
    for acc in rep_ids["shark_ray"]:
        assert acc_to_group.get(acc) == "shark_ray", f"Unexpected non-shark_ray in shark_ray clade: {acc}"
    # Amphibian clade must contain only curated amphibian sequences
    for acc in rep_ids["amphibian"]:
        assert acc_to_group.get(acc) == "amphibian", f"Unexpected non-amphibian in amphibian clade: {acc}"

    # Explicitly check regression cases from #127: smooth newt A0A6M9U035 and paddlefish E0X5L7
    assert "A0A6M9U035" not in rep_ids["reptile"]
    assert "E0X5L7" not in rep_ids["shark_ray"]


def test_fixed_frame_tracks_overflow_insertions_losslessly():
    # Synthesize a +5 aa insertion after W51 (exceeding the 3-slot 54A-54C budget)
    long_ins_g1 = HLA_A0201_G_ALPHA1[:51] + "MDLVK" + HLA_A0201_G_ALPHA1[51:]
    aln = map_groove_positions(long_ins_g1, HLA_A0201_G_ALPHA2)
    assert aln.has_overflow_insertions
    assert len(aln.fixed_seq) == 198
    assert len(aln.alpha1.residues) == 95
    assert "".join(r.aa for r in aln.alpha1.residues) == long_ins_g1
    assert "mdlvk" in aln.alpha1.a3m_seq


def test_left_aligned_insertion_slot_filling_and_no_transposition():
    # Issue #177: Mafa-F*02:01 inserts 'RY' after IMGT 20 (P20), NOT after 19.
    # It must NOT teleport 'RY' before P20 into 19A-19C; instead 19A-19C stay '---'
    # and pan_class_insertion_counts records +2 at HALF1:20.
    mafa = map_groove_positions(_MAFA_F0201_G1, _MAFA_F0201_G2)
    ab_slots = mafa.extract_imgt_slots([("G-ALPHA1", "19A"), ("G-ALPHA1", "19B"), ("G-ALPHA1", "19C")])
    assert ab_slots == "---"
    idx_20 = IMGT_PAN_CLASS_GROOVE_SLOTS.index(("HALF1", "20"))
    assert mafa.pan_class_insertion_counts[idx_20] == 2
    assert len(mafa.alpha1.residues) == 91
    assert "".join(r.aa for r in mafa.alpha1.residues) == _MAFA_F0201_G1

    # Direct 1-aa and 2-aa insertions right after IMGT 54 fill left-aligned (54A, 54B)
    # identically for Class I G-ALPHA1 and Class II G-ALPHA.
    c1_ins1 = align_domain_to_imgt(HLA_A0201_G_ALPHA1[:54] + "K" + HLA_A0201_G_ALPHA1[54:], "G-ALPHA1")
    c2_ins1 = align_domain_to_imgt(HLA_DRA0101_G_ALPHA[:50] + "K" + HLA_DRA0101_G_ALPHA[50:], "G-ALPHA")
    assert (c1_ins1.pan_class_slot_to_aa["54A"], c1_ins1.pan_class_slot_to_aa["54B"], c1_ins1.pan_class_slot_to_aa["54C"]) == ("K", "-", "-")
    assert (c2_ins1.pan_class_slot_to_aa["54A"], c2_ins1.pan_class_slot_to_aa["54B"], c2_ins1.pan_class_slot_to_aa["54C"]) == ("K", "-", "-")


# ---------------------------------------------------------------------------
# Class II (G-ALPHA / G-BETA) and Pan-Class IMGT G-DOMAIN alignment
# ---------------------------------------------------------------------------

_HLA_DRA0101_GROOVE = "IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYTPIT"
_HLA_DQA1_0101_GROOVE = "EDIVADHVASCGVNLYQFYGPSGQYTHEFDGDEEFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNSTAAT"
_HLA_DPA1_0103_GROOVE = "IKADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGQAFSFEAQGGLANIAILNNNLNTLIQRSNHTQAT"
_HLA_DRB1_0101_GROOVE = "GDTRPRFLWQLKFECHFFNGTERVRLLERCIYNQEESVRFDSDVGEYRAVTELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVGESFTVQR"
_HLA_DQB1_0201_GROOVE = "RDSPEDFVYQFKGMCYFTNGTERVRLVSRSIYNREEIVRFDSDVGEFRAVTLLGLPAAEYWNSQKDILERKRAAVDRVCRHNYQLELRTTLQR"
_HLA_DPB1_0101_GROOVE = "RATPENYVYQGRQECYAFNGTQRFLERYIYNREEYARFDSDVGEFRAVTELGRPAAEYWNSQKDILEEKRAVPDRVCRHNYELDEAVTLQR"


def test_class_ii_profile_loaded_with_verified_clades():
    import csv
    from pathlib import Path

    profile = load_class_ii_profile()
    assert profile["version"] == "1"
    assert profile["training_source"]["release"] == "v2.9.0"
    assert len(profile["training_source"]["sha256"]) == 64
    assert profile["domains"]["G-ALPHA"]["length"] == len(G_ALPHA_POSITIONS) == 82
    assert profile["domains"]["G-BETA"]["length"] == len(G_BETA_POSITIONS) == 89

    diverse_csv = Path(__file__).resolve().parent.parent / "mhcseqs" / "diverse_mhc_sequences.csv"
    with open(diverse_csv, "r", encoding="utf-8") as f:
        acc_to_group = {r["uniprot_accession"].strip(): r["source_group"].strip() for r in csv.DictReader(f)}

    for dom_key in ("G-ALPHA", "G-BETA"):
        rep_ids = profile["representative_ids"][dom_key]
        counts = profile["clade_representatives"][dom_key]
        assert set(rep_ids) == set(counts)
        for acc in rep_ids["reptile"]:
            assert acc_to_group.get(acc, "").startswith("reptile_")
        for acc in rep_ids["shark_ray"]:
            assert acc_to_group.get(acc) == "shark_ray"
        for acc in rep_ids["amphibian"]:
            assert acc_to_group.get(acc) == "amphibian"


def test_class_ii_galpha_and_gbeta_alignments_and_lossless_roundtrip():
    assert len(IMGT_G_ALPHA_SLOTS) == 87
    assert len(IMGT_G_BETA_SLOTS) == 98
    assert len(IMGT_CLASS_II_GROOVE_SLOTS) == 185
    assert len(IMGT_PAN_CLASS_HALF1_SLOTS) == 100
    assert len(IMGT_PAN_CLASS_HALF2_SLOTS) == 102
    assert len(IMGT_PAN_CLASS_GROOVE_SLOTS) == 202

    # Canonical reference strings align with identity == 1.0 and zero indels/overhangs
    ref_a = align_domain_to_imgt(HLA_DRA0101_G_ALPHA, "G-ALPHA")
    ref_b = align_domain_to_imgt(HLA_DRB1_0101_G_BETA, "G-BETA")
    assert ref_a.identity == 1.0 and not ref_a.has_indels and not ref_a.c_term_overhang
    assert ref_b.identity == 1.0 and not ref_b.has_indels and not ref_b.n_term_overhang

    dr = map_class_ii_groove_positions(_HLA_DRA0101_GROOVE, _HLA_DRB1_0101_GROOVE)
    assert dr.resolved_mhc_class == "II"
    assert dr.g_alpha is not None and dr.g_beta is not None
    assert len(dr.fixed_seq) == 185
    assert len(dr.pan_class_fixed_seq) == 202
    assert set(dr.fixed_seq) <= set(FIXED_FRAME_ALPHABET)
    assert set(dr.pan_class_fixed_seq) <= set(FIXED_FRAME_ALPHABET)
    assert "X" not in dr.fixed_seq
    assert "X" not in dr.pan_class_fixed_seq
    # Canonical DRA / DRB1 conserved landmarks
    assert dr.residue_at_imgt("G-ALPHA", "47") == "W"
    assert dr.residue_at_imgt("G-ALPHA", "86") == "N"
    assert dr.residue_at_imgt("G-BETA", "11") == "C"
    assert dr.residue_at_imgt("G-BETA", "59") == "W"
    assert dr.residue_at_imgt("G-BETA", "61A") == "Q"
    assert dr.residue_at_imgt("G-BETA", "61B") == "K"
    assert dr.residue_at_imgt("G-BETA", "74") == "C"
    # Terminal overhangs (T83 on DRA, GDTR on DRB1) round-trip losslessly in residues and a3m_seq
    assert dr.g_alpha.c_term_overhang == "T"
    assert dr.g_beta.n_term_overhang == "GDTR"
    assert "".join(r.aa for r in dr.g_alpha.residues) == _HLA_DRA0101_GROOVE
    assert "".join(r.aa for r in dr.g_beta.residues) == _HLA_DRB1_0101_GROOVE
    assert "".join(c.upper() for c in dr.g_alpha.a3m_seq if c.isalpha()) == _HLA_DRA0101_GROOVE
    assert "".join(c.upper() for c in dr.g_beta.a3m_seq if c.isalpha()) == _HLA_DRB1_0101_GROOVE

    # DQA1*01:01 has canonical A-strand 7A insertion and preserves W47 / N86
    dq = map_class_ii_groove_positions(_HLA_DQA1_0101_GROOVE, _HLA_DQB1_0201_GROOVE)
    assert dq.residue_at_imgt("G-ALPHA", "47") == "W"
    assert dq.residue_at_imgt("G-ALPHA", "86") == "N"
    assert dq.residue_at_imgt("G-BETA", "11") == "C"
    assert dq.residue_at_imgt("G-BETA", "59") == "W"
    assert dq.residue_at_imgt("G-BETA", "74") == "C"
    assert "".join(r.aa for r in dq.g_alpha.residues) == _HLA_DQA1_0101_GROOVE
    assert "".join(r.aa for r in dq.g_beta.residues) == _HLA_DQB1_0201_GROOVE

    # DPB1*01:01 has canonical 2-aa AB-turn deletion while keeping C11, W59, 61A, 61B, C74 locked
    dp = map_class_ii_groove_positions(_HLA_DPA1_0103_GROOVE, _HLA_DPB1_0101_GROOVE)
    assert dp.g_beta is not None
    assert len(dp.g_beta.deletions) == 2
    assert dp.residue_at_imgt("G-ALPHA", "47") == "W"
    assert dp.residue_at_imgt("G-ALPHA", "86") == "N"
    assert dp.residue_at_imgt("G-BETA", "11") == "C"
    assert dp.residue_at_imgt("G-BETA", "59") == "W"
    assert dp.residue_at_imgt("G-BETA", "74") == "C"
    assert "".join(r.aa for r in dp.g_beta.residues) == _HLA_DPB1_0101_GROOVE


def test_pan_class_fixed_frame_homologous_geometry_across_class_i_and_ii():
    c1 = map_groove_positions(HLA_A0201_G_ALPHA1, HLA_A0201_G_ALPHA2 + "D")
    dr = map_class_ii_groove_positions(_HLA_DRA0101_GROOVE, _HLA_DRB1_0101_GROOVE)
    dq = map_class_ii_groove_positions(_HLA_DQA1_0101_GROOVE, _HLA_DQB1_0201_GROOVE)
    dp = map_class_ii_groove_positions(_HLA_DPA1_0103_GROOVE, _HLA_DPB1_0101_GROOVE)

    # All four molecules live in the exact same 202-slot 3D-homologous coordinate system
    for aln in (c1, dr, dq, dp):
        assert len(aln.pan_class_fixed_seq) == 202
        assert aln.pan_class_fixed_slots == IMGT_PAN_CLASS_GROOVE_SLOTS
        # Conserved N-glycosylation at HALF1:86, Cys11/W59/Cys74 at HALF2:11/59/74
        assert aln.extract_pan_class_slots([("HALF1", "86"), ("HALF2", "11"), ("HALF2", "59"), ("HALF2", "74")]) == "NCWC"

    # Class-defining structural differences appear cleanly as '-' (structural gap), NEVER 'X':
    # Class II open N-helix wall (HALF1 57-60) is '----' while Class I has 'PEYW'
    assert c1.extract_pan_class_slots([("HALF1", "57"), ("HALF1", "58"), ("HALF1", "59"), ("HALF1", "60")]) == "PEYW"
    assert dr.extract_pan_class_slots([("HALF1", "57"), ("HALF1", "58"), ("HALF1", "59"), ("HALF1", "60")]) == "----"
    # Class II G-BETA 61B occupies the same slot reserved for Class I non-mammal helix-kink insertions
    assert c1.extract_pan_class_slots([("HALF2", "61A"), ("HALF2", "61B")]) == "A-"
    assert dr.extract_pan_class_slots([("HALF2", "61A"), ("HALF2", "61B")]) == "QK"
    # Class II G-BETA shorter C-end (HALF2 89-92) is '----' while Class I has 'LQRT'
    assert c1.extract_pan_class_slots([("HALF2", "89"), ("HALF2", "90"), ("HALF2", "91"), ("HALF2", "92")]) == "LQRT"
    assert dr.extract_pan_class_slots([("HALF2", "89"), ("HALF2", "90"), ("HALF2", "91"), ("HALF2", "92")]) == "----"

    # Single-chain Class II record via align_groove_to_imgt auto-detects mhc_class="II" and marks missing partner as 'X'
    alpha_only = align_groove_to_imgt({"mhc_class": "II", "groove1": _HLA_DRA0101_GROOVE, "groove2": ""})
    assert alpha_only.resolved_mhc_class == "II"
    assert len(alpha_only.fixed_seq) == 185
    assert alpha_only.fixed_seq[87:] == "X" * 98
    assert alpha_only.pan_class_fixed_seq[100:] == "X" * 102
    assert alpha_only.a3m_seq[-89:] == "-" * 89

    beta_only = align_groove_to_imgt({"mhc_class": "II", "groove1": "", "groove2": _HLA_DRB1_0101_GROOVE})
    assert beta_only.resolved_mhc_class == "II"
    assert len(beta_only.fixed_seq) == 185
    assert beta_only.fixed_seq[:87] == "X" * 87
    assert beta_only.pan_class_fixed_seq[:100] == "X" * 100
    assert beta_only.a3m_seq[:82] == "-" * 82


def test_class_ii_cache_matches_uncached_dp():
    for seq, dom in ((_HLA_DRA0101_GROOVE, "G-ALPHA"), (_HLA_DQA1_0101_GROOVE, "G-ALPHA"), (_HLA_DPB1_0101_GROOVE, "G-BETA")):
        expected = align_domain_to_imgt(seq, dom, use_fast_path=False)
        cached = align_domain_to_imgt(seq, dom, use_fast_path=True)
        assert cached == expected
        cached.pan_class_slot_to_aa.clear()
        assert align_domain_to_imgt(seq, dom, use_fast_path=True) == expected


def test_pan_class_pocket_slot_union_and_netmhciipan_exact_pseudosequence():
    """Issue #155: exact 22/34 overlap between NetMHCpan-34 and NetMHCIIpan-34 -> 46-slot union."""
    assert len(IMGT_PAN_CLASS_NETMHCPAN_SLOTS) == 34
    assert len(IMGT_PAN_CLASS_NETMHCIIPAN_SLOTS) == 34
    overlap = set(IMGT_PAN_CLASS_NETMHCPAN_SLOTS) & set(IMGT_PAN_CLASS_NETMHCIIPAN_SLOTS)
    assert len(overlap) == 22
    assert len(IMGT_PAN_CLASS_POCKET_SLOTS) == 46
    assert set(IMGT_PAN_CLASS_POCKET_SLOTS) == set(IMGT_PAN_CLASS_NETMHCPAN_SLOTS) | set(IMGT_PAN_CLASS_NETMHCIIPAN_SLOTS)
    assert all(slot in IMGT_PAN_CLASS_GROOVE_SLOTS for slot in IMGT_PAN_CLASS_POCKET_SLOTS)

    # Exact NetMHCpan-4.1 34-aa pseudosequence for HLA-A*02:01
    c1 = map_groove_positions(HLA_A0201_G_ALPHA1, HLA_A0201_G_ALPHA2)
    assert c1.extract_pan_class_slots(IMGT_PAN_CLASS_NETMHCPAN_SLOTS) == "YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY"

    # Exact NetMHCIIpan-3.0/4.0 34-aa pseudosequence for HLA-DRA*01:01 / HLA-DRB1*01:01
    dr1 = map_class_ii_groove_positions(_HLA_DRA0101_GROOVE, _HLA_DRB1_0101_GROOVE)
    assert dr1.extract_pan_class_slots(IMGT_PAN_CLASS_NETMHCIIPAN_SLOTS) == "QEFFIASGAAVDAIMWLFLECYDLQRATYHVGFT"

    # Both extract cleanly onto the 46-slot pan-class pocket union
    assert len(c1.extract_pan_class_pocket_seq()) == 46
    assert len(dr1.extract_pan_class_pocket_seq()) == 46
    # Closed Class I A-cap (HALF1:59 = Y) vs open Class II P1 end (HALF1:59 = '-')
    idx_59 = IMGT_PAN_CLASS_POCKET_SLOTS.index(("HALF1", "59"))
    assert c1.extract_pan_class_pocket_seq()[idx_59] == "Y"
    assert dr1.extract_pan_class_pocket_seq()[idx_59] == "-"


def test_pan_class_frame_buckets_galpha_54c_and_tracks_insertion_counts():
    """Issue #155: 3-aa G-ALPHA N-helix hinge insertions fill pan-class 54A-54C without overflow; overflows record per-slot counts."""
    # Synthesize a 3-aa teleost-style insertion ('WSP') after DRA IMGT 54 (mature 50 'R')
    galpha_54c = HLA_DRA0101_G_ALPHA[:50] + "WSP" + HLA_DRA0101_G_ALPHA[50:]
    aln_a = align_domain_to_imgt(galpha_54c, "G-ALPHA")
    assert aln_a.has_overflow_insertions  # 87-slot G-ALPHA class frame only has 54A-54B
    assert not aln_a.has_pan_class_overflow_insertions  # 100-slot pan HALF1 captures all 3 in 54A-54C
    assert aln_a.pan_class_slot_to_aa["54A"] == "W"
    assert aln_a.pan_class_slot_to_aa["54B"] == "S"
    assert aln_a.pan_class_slot_to_aa["54C"] == "P"
    assert sum(aln_a.pan_class_insertion_counts) == 0

    # Synthesize a +5 aa insertion after Class I IMGT 54 (3 fill 54A-54C, 2 overflow -> count 2 at 54C)
    long_ins_g1 = HLA_A0201_G_ALPHA1[:54] + "MDLVK" + HLA_A0201_G_ALPHA1[54:]
    aln_c1 = map_groove_positions(long_ins_g1, HLA_A0201_G_ALPHA2 + "D")
    feats = aln_c1.to_pan_class_features()
    assert feats["mhc_class"] == "I"
    assert len(feats["seq"]) == 202
    assert len(feats["pocket_seq"]) == 46
    assert len(feats["insertion_counts"]) == 202
    assert aln_c1.has_pan_class_overflow_insertions
    idx_54c = IMGT_PAN_CLASS_GROOVE_SLOTS.index(("HALF1", "54C"))
    assert feats["insertion_counts"][idx_54c] == 2
    assert sum(feats["insertion_counts"]) == 2
    assert feats["c_term_overhang_half2"] == "D"


def test_pan_class_monotonic_subsequence_galpha_40_41_and_92a_regression():
    """Issue #177: zero residue transposition, no G-ALPHA 40/41 double-counting, no spurious 92A 'X'."""

    def _is_subseq(sub: str, full: str) -> bool:
        it = iter(full)
        return all(ch in it for ch in sub)

    # 1. Chicken BF2*21:01 (55.1='T', 89-aa alpha1): NT must stay in order and HALF1:92A must be '-' (not 'X')
    bf2_g1 = "ELHTLRYIRTAMTDPGPGLPWFVDVGYVDGELFMHYNSTARRAVPRTEWIAANTDQQYWDRETQIVQGSEQINRENLDILRRRYNQTG"
    bf2_g2 = "GSHTVQWMSGCDILEDGTIRGYHQAAYDGRDFVAFDKGTMTLTAAVPEAVPTKRKWEEGGYAEGLKQYLEETCVEWLRRYVEYGKAELGRRE"
    bf2 = map_groove_positions(bf2_g1, bf2_g2)
    assert bf2.alpha1 is not None and bf2.alpha2 is not None
    assert _is_subseq(bf2.alpha1.fixed_seq.replace("-", "").replace("X", ""), bf2_g1)
    assert _is_subseq(bf2.alpha1.pan_class_fixed_seq.replace("-", "").replace("X", ""), bf2_g1)
    assert _is_subseq(bf2.alpha2.pan_class_fixed_seq.replace("-", "").replace("X", ""), bf2_g2)
    assert bf2.alpha1.pan_class_slot_to_aa["90"] == "X"
    assert bf2.alpha1.pan_class_slot_to_aa["92A"] == "-"
    # Off-boundary 55.1='T' stays in insertion_counts['55'] rather than jumping ahead of 55='N'
    idx_55 = IMGT_PAN_CLASS_GROOVE_SLOTS.index(("HALF1", "55"))
    assert bf2.pan_class_insertion_counts[idx_55] == 1

    # 2. G-ALPHA 39.1 insertion populates pan-class slot 40 with ZERO double-counted overflow
    a0 = align_domain_to_imgt(HLA_DRA0101_G_ALPHA, "G-ALPHA")
    idx39 = [r.query_index for r in a0.residues if r.imgt_pos == "39"][0]
    galpha_ins39 = HLA_DRA0101_G_ALPHA[: idx39 + 1] + "Q" + HLA_DRA0101_G_ALPHA[idx39 + 1 :]
    a39 = align_domain_to_imgt(galpha_ins39, "G-ALPHA")
    assert _is_subseq(a39.pan_class_fixed_seq.replace("-", "").replace("X", ""), galpha_ins39)
    if any(ins.after_imgt_pos == "39" for ins in a39.insertions):
        assert a39.pan_class_slot_to_aa["40"] == "Q"
        assert not a39.has_pan_class_overflow_insertions
        assert sum(a39.pan_class_insertion_counts) == 0

    # 3. Direct unit check on _build_pan_class_frame for G-ALPHA 39.1 -> slot 40 without overflow
    from mhcseqs.imgt import GrooveInsertion, _build_fixed_frame, _build_pan_class_frame

    ins39 = (
        GrooveInsertion(
            domain="G-ALPHA",
            after_mature_pos=35,
            after_imgt_pos="39",
            sequence="Q",
            query_start=35,
            imgt_codes=("39.1",),
            structural_element="CD-TURN",
        ),
    )
    _, _, fixed_map, fixed_ov = _build_fixed_frame("G-ALPHA", a0.imgt_to_aa, ins39)
    _, _, pan_map, pan_ov, pan_counts = _build_pan_class_frame("G-ALPHA", fixed_map, a0.imgt_to_aa, ins39)
    assert fixed_map["40"] == "Q" and not fixed_ov
    assert pan_map["40"] == "Q" and not pan_ov and sum(pan_counts) == 0

    # 4. Monotonic subsequence check across Ptal-N, Trvu-UB, Mafa-F, and human Class II
    for g1, g2, cls in [
        (_PTAL_N0101_G1, _PTAL_N0101_G2, "I"),
        (_TRVU_UB0101_G1, _TRVU_UB0101_G2, "I"),
        (_MAFA_F0201_G1, _MAFA_F0201_G2, "I"),
        (_HLA_DQA1_0101_GROOVE, _HLA_DQB1_0201_GROOVE, "II"),
        (_HLA_DPA1_0103_GROOVE, _HLA_DPB1_0101_GROOVE, "II"),
    ]:
        aln = align_groove_to_imgt(g1, g2, mhc_class=cls)
        assert _is_subseq(aln.alpha1.fixed_seq.replace("-", "").replace("X", ""), g1)
        assert _is_subseq(aln.alpha1.pan_class_fixed_seq.replace("-", "").replace("X", ""), g1)
        assert _is_subseq(aln.alpha2.fixed_seq.replace("-", "").replace("X", ""), g2)
        assert _is_subseq(aln.alpha2.pan_class_fixed_seq.replace("-", "").replace("X", ""), g2)
