"""Tests for IMGT G-DOMAIN numbering and mutant MHC validation."""

import pytest

from mhcseqs.domain_parsing import decompose_class_i
from mhcseqs.imgt import (
    CONSERVED_CYS_POSITIONS,
    GALPHA1_POSITIONS,
    GALPHA2_POSITIONS,
    HELIX_INSERTIONS,
    HLA_A0201_GALPHA1,
    HLA_A0201_GALPHA2,
    DomainAlignment,
    GrooveAlignment,
    align_domain_to_imgt,
    align_groove_to_imgt,
    imgt_to_mature,
    imgt_to_mature_class_i,
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
    for imgt_pos in GALPHA2_POSITIONS:
        mature = imgt_to_mature_class_i("G-ALPHA2", imgt_pos)
        assert mature is not None
        domain, back = mature_to_imgt_class_i(mature)
        assert domain == "G-ALPHA2"
        assert back == imgt_pos


def test_galpha2_total_positions():
    """G-ALPHA2 should have 92 occupied positions (for HLA-A*02:01)."""
    assert len(GALPHA2_POSITIONS) == 92


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
    assert len(GALPHA1_POSITIONS) == 90


def test_align_fast_path_matches_full_dp():
    fast1 = align_domain_to_imgt(HLA_A0201_GALPHA1, "G-ALPHA1", use_fast_path=True)
    full1 = align_domain_to_imgt(HLA_A0201_GALPHA1, "G-ALPHA1", use_fast_path=False)
    assert fast1.aligned_ref == full1.aligned_ref
    assert fast1.aligned_query == full1.aligned_query
    assert fast1.mature_to_aa == full1.mature_to_aa

    fast2 = align_domain_to_imgt(HLA_A0201_GALPHA2 + "D", "G-ALPHA2", use_fast_path=True)
    full2 = align_domain_to_imgt(HLA_A0201_GALPHA2 + "D", "G-ALPHA2", use_fast_path=False)
    assert fast2.aligned_ref == full2.aligned_ref
    assert fast2.aligned_query == full2.aligned_query
    assert fast2.c_term_overhang == "D"


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
    row = {"groove1": HLA_A0201_GALPHA1, "groove2": HLA_A0201_GALPHA2}
    aln_dict = align_groove_to_imgt(row)
    aln_str = align_groove_to_imgt(HLA_A0201_GALPHA1, HLA_A0201_GALPHA2)
    assert aln_dict.mature_to_aa == aln_str.mature_to_aa

    empty_dom = align_domain_to_imgt("", "G-ALPHA1")
    assert isinstance(empty_dom, DomainAlignment)
    assert empty_dom.n_term_missing == 90

    with pytest.raises(TypeError):
        align_groove_to_imgt(12345)
