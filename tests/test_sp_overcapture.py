"""Tests for signal-peptide over-capture detection.

Over-capture happens when a sequence is translated from an in-frame ATG upstream
of the true start codon (a weak-Kozak AUG the ribosome would leaky-scan past),
prepending translated 5'-UTR to the real signal peptide.  At the amino-acid
level this looks like a long SP whose charged/polar N-terminal portion has no
h-region, followed by an internal Met beginning a clean self-contained SP that
cleaves at the same mature start.

This is a QC signal only — see ``detect_sp_internal_met_overcapture`` docstring.
"""

from mhcseqs import parse_class_ii_alpha
from mhcseqs.domain_parsing import detect_sp_internal_met_overcapture

# Mamu-DPA1*02:01:01:01 (IPD-MHC).  Full record translates from MDINYRPHNVCPEDR
# (15 aa, charged/polar, no h-region) before the real DPA1 leader
# MFQTRAIVLRTLSLAFLLSLRGAGA -> mature IKADH...  The full sequence has no
# detectable SP; the internal Met at 15 recovers a textbook one.
MAMU_DPA1 = (
    "MDINYRPHNVCPEDRMFQTRAIVLRTLSLAFLLSLRGAGAIKADHVSTYAAFVQTHRPTGEFMFEFDEDEQ"
    "FYVDLDKKETVWHLEEFGRAFSFEAQGGLANIAILNNNLNITIQRSNYTQAANDPPEVTVFPKEPVALGQP"
    "NTLICHIDKFFPPVLNVTWLCNGEPVTEGVAESLFLPRTDYSFHKFHYLTFVPSAEDYYDCRVEHWGLDQP"
    "LLKHWEAQEPIQIPETTETVLCALGLVLGLVGIIVGTVLIIKSLRSGRDPRAQGPL"
)

# A0A8C8UCV8 (UniProt/TrEMBL, class I).  Charged prefix MRRSRSEVRTMESGTR (16 aa)
# precedes the poly-Leu leader MLLLLLAAVLSPTRTRA -> mature GSHSLRYF...
TREMBL_CLASS_I = (
    "MRRSRSEVRTMESGTRMLLLLLAAVLSPTRTRAGSHSLRYFYTIVSRPGLGEPRFIIVGYVDDTQFVHYDS"
    "DAETPRMEPLAPWVEREGPEYWERETQKAKITQQTFRRNLRALSSYYNHSHEGSHTIQVISGCDVDSNGRL"
    "LRGYEQFAYDGRDYISLNQDLRTWTAVDRVAQIIRLELEQAGEAERYRAYLEGECIERLCRYLENGTTNLH"
    "ARAWGLGQLFSLPLGWGSRDEEEQTQVMELVEIRPAGDGSFQKWAAVVVPFGEELKYTCHVQHEGLPEPLT"
    "LRWDSRQSQSFPIMTMVGAVLGAVVILGFIIGAVMMWKRKNTGMEGQGLSYFSFFIQHLQKRYRLNGKTLT"
    "HTHTCP"
)

# HLA-C*03:691 (IMGT).  The canonical class I leader MRVMAPRTL... has a genuine
# 3-aa charged n-region (MRV) before the h-region.  Must NOT be flagged: a real
# n-region is short, and the internal Met sits below the prefix-length cutoff.
HLA_C_REAL = (
    "MRVMAPRTLILLLSGALALTETWVGSHSMRYFYTAVSRPGRGEPHFIAVGYVDDTQFVRFDSDAASPRGEP"
    "RAPWVEQEGPEYWDRETQKYKRQAQTDRVSLRNLRGYYNQSEARSHIIQRMYGCDVGPDGRLLRGYDQYAY"
    "DGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQLRAYLEGLCVEWLRRYLKNGKETLQRAEHPKTHV"
    "THHPVSDHEATLRCWALGFYPAEITLTWQWDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCH"
    "VQHEGLPEPLTLRWEPSSQPTIPIVGIVAGLAVLAVLAVLGAVVAVVMCRRKSSGGKGGSCSQAASSNSAQ"
    "GSDESLIACKA"
)

# HLA-DQB1*05:377 (IMGT).  A genuine 32-aa signal peptide.  An internal Met at 20
# exists but yields a *worse* SP than the full leader, so it must NOT be flagged.
HLA_DQB1_REAL = (
    "MSWKKSLRIPGDLRVATVTLMLAILSSSLTEGRDSPEDFVYQFKGLCYFTNGTERVRGVTRHIYNREEYVR"
    "FDSDVGVYRAVTPQGRPVAEYWNSQKEVLEGARASVDRVCRHNYEVAYRGILQRRVEPTVTISPSRTEALN"
    "HHNLLICSVTDFYPSQIKVRWFRNDQEETAGVVSTPLIRNGDWTFQILVMLEMTPQRGDVYTCHVEHPSLQ"
    "SPITVEWRAQSESAQSKMLSGVGGFVLGLIFLGLGLIIRQRSRKGLLH"
)


def test_overcapture_detected_mamu_dpa1():
    # Real leader starts at the internal Met (position 15).
    assert detect_sp_internal_met_overcapture(MAMU_DPA1, 40) == 15


def test_overcapture_detected_trembl_class_i():
    assert detect_sp_internal_met_overcapture(TREMBL_CLASS_I, 33) == 16


def test_real_class_i_nregion_not_flagged():
    # MRV is a genuine n-region, not over-capture.
    assert detect_sp_internal_met_overcapture(HLA_C_REAL, 26) is None


def test_genuine_long_sp_not_flagged():
    # DQB1's 32-aa SP is real; the internal Met gives a worse SP.
    assert detect_sp_internal_met_overcapture(HLA_DQB1_REAL, 32) is None


def test_short_sp_never_flagged():
    # SPs at or below 25 aa are not candidates regardless of internal Mets.
    assert detect_sp_internal_met_overcapture(MAMU_DPA1, 20) is None


def test_no_internal_met_returns_none():
    # No internal Met in the SP region -> nothing to rescue.
    seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGSHSLRYF" + "A" * 200
    assert detect_sp_internal_met_overcapture(seq, 31) is None


def test_flag_surfaces_through_parser():
    # Guards the wiring: the QC flag must reach AlleleRecord.flags via the parse
    # path, not just the standalone detector.
    record = parse_class_ii_alpha(MAMU_DPA1, allele="Mamu-DPA1*02:01", gene="DPA1")
    assert "sp_internal_met_overcapture(15)" in record.flags
