"""Binding groove extraction from MHC protein sequences.

Extracts the peptide-binding groove from MHC chains.  The groove has two
structural halves (groove1, groove2):

  Class I  – a single alpha chain contributes both halves:
               groove1  (α1 domain, ~90 aa)
               groove2  (α2 domain, ~93 aa)

  Class II – two separate chains each contribute one half:
               alpha chain →  groove1  (α1 domain, ~83 aa)
               beta  chain →  groove2  (β1 domain, ~93 aa)

Note: "α2" means something different in each class.
  Class I  α2  =  C-terminal groove half (peptide-binding)
  Class II α2  =  Ig-fold support domain (NOT part of the groove)


Domain architecture and Cys-pair anchors
=========================================

CLASS I alpha chain:

  raw pos:  0       SP   ms        ms+90  C1  C1+10     C2+20
            |--------|    |          |     |    |          |
            [ signal ]    [-- groove1 --]--[--- groove2 ---]--[-- α3 Ig --]
            [ peptide]    [     α1      ]  [      α2       ]  [  support  ]
                          |<-- 90 aa -->|  | C1.......C2   |  | C1....C2  |
                                           |<--- 93 aa --->|
                                        c1-10           c2+20

  The α2 Cys pair (C100–C163 in mature protein) is INSIDE the groove.
  groove1 = seq[mature_start : c1 - 10]      (α1 domain)
  groove2 = seq[c1 - 10     : c2 + 20]       (α2 domain, wraps around the Cys pair)
  mature_start is inferred as:  raw_c1_pos - 100

CLASS II alpha chain:

  raw pos:  0       SP   ms           ms+83  C1         C2
            |--------|    |              |    |          |
            [ signal ]    [-- groove1 ---]----[-- α2 Ig -----]
            [ peptide]    [     α1       ]    [   support    ]
                          |<-- 83 aa --->|    | C1......C2   |
                                       c1-23

  The Cys pair (C106–C162 in mature) is in the Ig SUPPORT domain.
  groove1 = seq[mature_start : c1 - 23]      (α1 groove domain)
  mature_start is inferred as:  raw_c1_pos - 106

CLASS II beta chain:

  raw pos:  0       SP   ms  C1b       C2b      C1i         C2i
            |--------|    |   |          |        |           |
            [ signal ]    [------- groove2 ------]--[- β2 Ig -----]
            [ peptide]    [         β1           ]  [  support    ]
                          |    C1b......C2b      |  | C1i....C2i  |
                          |<------- 93 aa ------>|
                                               c1i-23

  β2 Ig Cys pair (C116–C172 in mature) anchors the groove boundary.
  β1 also has its own Cys pair (C14–C78) INSIDE the groove.
  groove2 = seq[mature_start : β2_c1 - 23]   (β1 groove domain)
  mature_start is inferred as:  raw_β2_c1_pos - 116
  Fallback: if no β2 pair, use β1 pair:  groove2 = seq[ms : β1_c2 + 15]


Partial / SP-stripped sequences
================================

Sequences may be deposited without the signal peptide or with N-terminal
truncations.  The parser handles this gracefully because the Cys-pair
anchor positions are within the mature protein:

  - Full with SP:     mature_start = raw_c1 - 100 → correct SP boundary
  - SP stripped:      mature_start = raw_c1 - 100 = 0 → groove starts at pos 0
  - N-term truncated: mature_start = 0, groove1 is shorter by the missing residues

groove2 is always correct (anchored by the Cys pair which is present in the
sequence).  groove1 may lose N-terminal residues proportional to the
truncation, but the rest of the groove is correctly parsed.


Mature-position constants are calibrated against UniProt signal-peptide
annotations for canonical human alleles (see individual constant comments).
In practice these positions are universally conserved across all species in
IMGT/HLA and IPD-MHC — zero exceptions found among 25,832 verifiable entries
for class I, and 100% match for class II Ig-domain Cys positions.


Structural parsing approach (March 2026)
=========================================

In addition to Cys-pair position constants, the parser uses conserved
structural motifs to classify domain boundaries:

Signal peptide (SP): 14-42 aa in vertebrate MHC (human: typically 21-32 aa)
  - 0 (absent) for SP-stripped database deposits
  - Detection: von Heijne -3/-1 rule (A/G/S at -1, small aliphatic at -3),
    hydrophobic core upstream, KD hydrophobicity transition at boundary,
    H at mature position +3 (conserved structural His in α1 S1 strand)
  - MAX_PLAUSIBLE_SP = 50 (hard cap; no known MHC SP exceeds ~42 aa)

Groove domains: 75-95 aa per half (standard), 70-125 aa (extended range)
  - Class I α1 (groove1): no Cys pair, 90 aa in human, 88 aa in birds,
    up to 109 aa in owls/passerines with α1 insertions (Moriyama 2021)
  - Class I α2 (groove2): Cys pair at ~C101-C164 (sep ~63), flanked by
    Q@-5 (90%), M@-3 (74%), G@-1 (85%), W@+3 (72%), L@+4 (86%)
  - Class II β1 (groove2): Cys pair at ~C14-C78 (sep ~64), flanked by
    E@-1 (66%), K@-3 (52%), N@+3 (81%), Y@+4 (81%) — different from α2!
  - Class II α1 (groove1): no Cys pair, 83 aa in human (DRA), up to 97 (DMA)

Ig support domains: 75-100 aa
  - Universal Ig: Y at Cys2-2 (64-97%), V at Cys2+2 (74-93%)
  - Class I α3: LR-C-WA at Cys1 (highly distinctive), sep ~56
  - Class II β2: diverse Cys1, Y@c2-2, V@c2+2, H@c2+4, W@c2+15 (72%)
  - W at Cys2+15: 51% in Ig overall, 14% in groove — useful but not definitive

α1/α2 junction motif (class I groove1/groove2 boundary):
  - G at pos 0 (77%), S at +1 (46%), H at +2 (86%), T at +3 (54%), Q at +5 (90%)
  - Used to discover α1 domain length without per-species constants

Groove→Ig boundary motifs (class-specific):
  - Class I:   L at -5 (87%), P at +1 (82%)
  - Class II α: N at 0 (64%), PP at +2/+3 (76%/78%)
  - Class II β: R at -1 (83%), P at +3 (91%), V at +5 (89%)

Ported from presto/data/groove.py.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field, replace
from typing import Optional, Sequence

from .alleles import infer_gene, normalize_mhc_class

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Ig-fold Cys-Cys pair separation range
IG_SEP_MIN = 48
IG_SEP_MAX = 72

# Dominant Cys-Cys separations per domain (used for pair selection)
CLASS_I_ALPHA2_DOMINANT_SEP = 63
CLASS_II_ALPHA_IG_DOMINANT_SEP = 56
CLASS_II_BETA1_DOMINANT_SEP = 64
CLASS_II_BETA2_DOMINANT_SEP = 56

# Minimum sequence length to attempt groove parsing
MIN_GROOVE_SOURCE_LEN = 70

# Class I constants – validated against UniProt (P01892 HLA-A*02:01)
CLASS_I_ALPHA2_CYS1_MATURE_POS = 100  # α2 Ig-fold Cys1 at mature pos 100
CLASS_I_ALPHA3_CYS1_MATURE_POS = 202  # α3 Ig-fold Cys1 at mature pos 202
CLASS_I_ALPHA2_CYS1_OFFSET = 10
CLASS_I_ALPHA2_END_AFTER_CYS2 = 20
CLASS_I_ALPHA2_CYS1_RAW_MIN = 35  # lowered from 60 to catch N-terminal truncations
CLASS_I_ALPHA2_CYS1_RAW_MAX = 180
CLASS_I_ALPHA3_CYS1_RAW_MIN = 180

# Class II alpha constants – validated against UniProt (P01903 HLA-DRA*01:01)
CLASS_II_ALPHA_IG_CYS1_MATURE_POS = 106  # α2 Ig-fold Cys1 at mature pos 106 (default)
CLASS_II_ALPHA_GROOVE_END_BEFORE_IG_CYS = 23
CLASS_II_ALPHA_CYS1_RAW_PRIMARY_MIN = 100
CLASS_II_ALPHA_CYS1_RAW_MIN = 40  # lowered from 80 to catch DMA (Cys1 ~49)
CLASS_II_ALPHA_CYS1_RAW_MAX = 160

# Gene-specific Cys1 mature positions for class II alpha.
# The α1 groove domain varies in length across gene families, shifting the
# Ig-fold Cys1 position in the mature protein.
# Verified against UniProt signal peptide annotations:
#   DRA:  106 (P01903, SP=25)  – default
#   DQA:  109 (P01909, SP=23)  – α1 domain 3 residues longer
#   DMA:  120 (P28067, SP=26)  – α1 domain 14 residues longer; has α1 intra-domain Cys pair
#   DPA:  106 (P20036, SP=31)  – matches default
#   DOA:  106 (P06340, SP=26)  – matches default
_CLASS_II_ALPHA_CYS1_MATURE_POS_BY_GENE: dict[str, int] = {
    "DQA": 109,
    "DMA": 120,
}

# Class II beta constants – validated against UniProt (P01911 HLA-DRB1*01:01)
CLASS_II_BETA1_CYS1_MATURE_POS = 14  # β1 Ig-fold Cys1 at mature pos 14
CLASS_II_BETA2_CYS1_MATURE_POS = 116  # β2 Ig-fold Cys1 at mature pos 116 (default)

# Gene-specific Cys1 mature positions for class II beta.
# Verified against UniProt signal peptide annotations:
#   DRB:  116 (P01911, SP=29)  – default
#   DQB:  116 (P01920, SP=32)  – matches default
#   DPB:  114 (P04440, SP=29)  – β1 domain 2 residues shorter
#   DMB:  116 (P28068, SP=18)  – matches default
#   DOB:  116 (P13765, SP=26)  – matches default
_CLASS_II_BETA2_CYS1_MATURE_POS_BY_GENE: dict[str, int] = {
    "DPB": 114,
}
CLASS_II_BETA1_CYS1_RAW_MIN = 2  # lowered from 20 to catch SP-stripped entries
CLASS_II_BETA1_CYS1_RAW_MAX = 95
CLASS_II_BETA2_CYS1_RAW_MIN = 100
CLASS_II_BETA2_CYS1_RAW_MAX = 180
CLASS_II_BETA_GROOVE_END_BEFORE_BETA2_CYS = 23
CLASS_II_BETA1_ONLY_END_AFTER_CYS2 = 15

# Default groove lengths (used only in fallback paths)
DEFAULT_CLASS_I_GROOVE1_LEN = 90  # α1 groove half
DEFAULT_CLASS_I_GROOVE2_LEN = 93  # α2 groove half
DEFAULT_CLASS_II_GROOVE1_LEN = 83  # α1 groove half (class II alpha chain)
DEFAULT_CLASS_II_GROOVE2_LEN = 93  # β1 groove half (class II beta chain)

# Ig domain extraction: residues after Cys2 to include in Ig domain
IG_DOMAIN_END_AFTER_CYS2 = 20

# Fragment fallback thresholds
CLASS_I_FRAGMENT_MAX_LEN = 200  # class I groove1+groove2 = ~183 aa
CLASS_II_ALPHA_FRAGMENT_MAX_LEN = 110
CLASS_II_BETA_FRAGMENT_MAX_LEN = 120

# Maximum plausible signal peptide length.  No known MHC species has an SP
# longer than ~42 aa (DMA).  Anything beyond 50 indicates the Cys-pair
# anchor is wrong — typically because a point mutation destroyed the
# conserved Ig-fold cysteine.
MAX_PLAUSIBLE_SP = 50

# Genes whose MHC-like fold does not form a peptide-binding groove.
# These are excluded from groove extraction in the pipeline.
NON_GROOVE_GENES = frozenset({"MICA", "MICB", "MIC1", "MIC2", "HFE", "B2M", "MR1"})

# Non-classical MHC class I lineages in teleost fish.  These genes have
# MHC-like folds but do NOT form a classical peptide-binding groove.
# See Grimholt et al. 2015 (BMC Evol Biol) and Malmstrøm et al. 2019.
#   L lineage: highly variable, non-classical, some processed/intronless
#   S lineage: lacks classical groove residues
#   P lineage: extra Cys in α1, altered groove shape
#   H lineage: α3 lost, α1/α2 often deteriorated or entirely absent
NON_CLASSICAL_CLASS_I_GENE_PATTERNS = (
    "MHC1L",
    "MHC1S",
    "MHC1P",  # zebrafish-style: mhc1laa, mhc1saa, etc.
    "LLA",
    "LCA",
    "LDA",
    "LFA",
    "LGA",
    "LIA",
    "LJA",  # L lineage locus names
    "MFSD",  # NOT MHC at all — lipid transporter contaminant
)

# Known non-MHC proteins that have leaked into curated datasets via automated
# genome annotation.  Keyed by UniProt accession.
NON_MHC_ACCESSIONS = frozenset(
    {
        "Q1LUQ4",  # Dare-mfsd6a, zebrafish lipid transporter
        "B0UYT5",  # Dare-mfsd6b, zebrafish lipid transporter
    }
)

# Minimum groove half length considered potentially functional for peptide
# binding.  Below this, the α-helix + β-sheet architecture cannot form a
# complete groove wall.  Set conservatively: 70 aa allows for genuine fish
# variation while flagging molecules that almost certainly lack a functional
# groove.
MIN_FUNCTIONAL_GROOVE_HALF_LEN = 70

# Gene prefix patterns for class II chain inference.
# Covers mammalian D-series, chicken B-locus, fish (DA/DB/DC/DD/DE groups),
# and ruminant-specific DY genes.
CLASS_II_ALPHA_GENE_PREFIXES = (
    "DRA",
    "DQA",
    "DPA",
    "DMA",
    "DOA",  # mammalian standard
    "DYA",
    "DNA",  # ruminant-specific
    "DAA",
    "DBA",
    "DCA",
    "DDA",
    "DEA",  # fish class II alpha
    "BLA",  # chicken class II alpha
)
CLASS_II_BETA_GENE_PREFIXES = (
    "DRB",
    "DQB",
    "DPB",
    "DMB",
    "DOB",  # mammalian standard
    "DYB",
    "DIB",  # ruminant-specific
    "DAB",
    "DBB",
    "DCB",
    "DDB",
    "DEB",  # fish class II beta
    "BLB",  # chicken class II beta
)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RawAllele:
    """Unparsed allele record: identity, species, and sequence data.

    This is what you have before groove extraction — the raw sequence
    plus metadata from the FASTA header or CSV.
    """

    allele: str = ""
    full_allele: str = ""
    gene: str = ""
    mhc_class: str = ""
    chain: str = ""  # "alpha" or "beta"
    species: str = ""
    species_category: str = ""
    species_prefix: str = ""
    source: str = ""
    source_id: str = ""  # database accession (HLA00001, NHP00001, P01901, …)
    sequence: str = ""  # full protein (may include signal peptide)
    seq_len: int = 0
    mature_start: int = 0

    @property
    def mature_sequence(self) -> str:
        """The mature protein sequence (signal peptide removed)."""
        if self.sequence:
            return self.sequence[self.mature_start :]
        return ""


@dataclass(frozen=True)
class AlleleRecord(RawAllele):
    """Fully parsed allele record with groove domains, Ig domain, and tail.

    The peptide-binding groove has two structural halves:
      groove1  (N-terminal):  α1 in class I, α1 in class II alpha
      groove2  (C-terminal):  α2 in class I, β1 in class II beta

    Class I alpha provides both halves; class II alpha and beta each provide one.
    """

    groove_seq: str = ""  # concatenation of groove halves
    groove1: str = ""  # α1 (class I & II alpha), empty (class II beta)
    groove2: str = ""  # α2 (class I), β1 (class II beta), empty (class II alpha)
    groove1_len: int = 0
    groove2_len: int = 0
    ig_domain: str = ""  # Ig support domain: α3 (class I), α2 (class II α), β2 (class II β)
    ig_domain_len: int = 0
    tail: str = ""  # everything after Ig domain (TM + cytoplasmic)
    tail_len: int = 0
    status: str = "ok"
    anchor_type: str = ""
    anchor_cys1: Optional[int] = None
    anchor_cys2: Optional[int] = None
    secondary_cys1: Optional[int] = None
    secondary_cys2: Optional[int] = None
    mutations: tuple[str, ...] = field(default_factory=tuple)
    flags: tuple[str, ...] = field(default_factory=tuple)

    @property
    def ok(self) -> bool:
        return self.status in {
            "ok",
            "inferred_from_alpha3",
            "alpha1_only",
            "alpha2_only",
            "beta1_only_fallback",
            "fragment_fallback",
        }

    @property
    def mature_sequence(self) -> str:
        """The mature protein sequence, materialized from sequence or domain parts."""
        if self.sequence:
            return self.sequence[self.mature_start :]
        # Reconstruct from parsed domains
        if self.mhc_class == "I":
            return self.groove1 + self.groove2 + self.ig_domain + self.tail
        if self.chain == "alpha":
            return self.groove1 + self.ig_domain + self.tail
        if self.chain == "beta":
            return self.groove2 + self.ig_domain + self.tail
        return ""


# ---------------------------------------------------------------------------
# Structural template scoring — Cys pair classification & motif detection
# ---------------------------------------------------------------------------
#
# These functions classify Cys-Cys pairs as "groove domain" or "Ig domain"
# based on conserved flanking residues, and detect domain boundaries using
# structural motifs.  They replace per-gene/per-species position constants
# with universal biochemical signals.


@dataclass(frozen=True)
class CysPairAnnotation:
    """A Cys-Cys pair with structural classification scores."""

    c1: int
    c2: int
    separation: int
    groove_score: float  # higher = more likely groove domain (α2, β1)
    ig_score: float  # higher = more likely Ig domain (α3, β2, class II α2)
    domain_type: str  # "groove" | "ig" | "ambiguous"
    evidence: tuple[str, ...] = ()  # human-readable scoring details


@dataclass
class ParseTrace:
    """Full debug trace of parsing decisions.

    Attached to an AlleleRecord when trace=True is passed to a parser.
    Contains all intermediate scoring decisions for debugging ground
    truth disagreements.
    """

    sequence_len: int = 0
    mhc_class: str = ""
    cys_pairs: tuple[CysPairAnnotation, ...] = ()
    junction_pos: int = 0
    junction_score: float = 0.0
    sp_candidates: tuple[tuple[int, float], ...] = ()  # (pos, score) top candidates
    selected_sp: int = 0
    boundary_scores: dict[str, float] = field(default_factory=dict)
    log: list[str] = field(default_factory=list)

    def summary(self) -> str:
        """Human-readable summary of parsing decisions."""
        lines = list(self.log)
        return "\n".join(lines)


def _gaussian_score(value: int, mean: int, sigma: int) -> float:
    """Gaussian-like score: 1.0 at mean, decaying with distance."""
    import math

    return math.exp(-0.5 * ((value - mean) / max(sigma, 1)) ** 2)


# ---- Groove-domain Cys flanking scoring tables ----
#
# Class I α2 (n=794) and class II β1 (n=370) have almost entirely
# different flanking residues at the same structural positions because
# the two groove halves evolved independently.  The classifier scores
# against BOTH sets and takes the maximum, which also serves as a
# class-inference heuristic when class is unknown.
#
# Class I α2 groove: G-type Cys1, aromatic/Trp Cys2
# Class II β1 groove: charged Cys1, polar/aromatic Cys2

# Class I α2 Cys1-1: dominated by small residues (G 85%)
# Property: small/neutral at the Cys1 junction — structural role in
# the β-strand turn entering the Ig-fold.  P and charged are absent
# because they would disrupt the turn.
_GROOVE_I_CYS1_M1: dict[str, float] = {
    "G": 2.0,  # 85% — dominant, Gly allows tight turn
    "A": 1.0,  # 3% — small, same structural role
    "S": 1.0,  # 2% — small/polar
    "T": 0.5,  # <1% — small but bulkier
    "V": 0.3,  # 2% — occasional hydrophobic variant
    # D/E absent: charged residues disrupt the β-turn
    # P absent: Pro breaks the β-strand
}
# Class I α2 Cys1-5: polar/amide dominated (Q 90%)
# Property: exposed position on loop before Cys, hydrogen-bonding role.
_GROOVE_I_CYS1_M5: dict[str, float] = {
    "Q": 1.5,  # 90% — amide, H-bond donor/acceptor
    "N": 0.5,  # 2% — shorter amide, same property
    "E": 0.3,  # 1% — charged variant of Q
    # K/R absent: positive charge unfavorable here
}
# Class I α2 Cys1-3: hydrophobic core (M 74%)
# Property: buried position, hydrophobic packing with Ig β-sheet.
_GROOVE_I_CYS1_M3: dict[str, float] = {
    "M": 1.0,  # 74% — flexible hydrophobic, S-containing
    "L": 0.5,  # 8% — branched hydrophobic
    "I": 0.5,  # 5% — branched hydrophobic
    "F": 0.3,  # 2% — aromatic hydrophobic
    "V": 0.3,  # 3% — smaller hydrophobic
    # A/G absent: too small for this buried position
    # Charged absent: buried position requires hydrophobic packing
}
# Class I α2 Cys2+3: aromatic, often Trp (W 72%)
# Property: the conserved Ig-fold Trp (not the canonical C2+15 Trp,
# but a groove-domain-specific aromatic at the C-terminal Cys flank).
_GROOVE_I_CYS2_P3: dict[str, float] = {
    "W": 1.5,  # 72% — aromatic, stacks against β-sheet
    "F": 0.5,  # 5% — aromatic alternative
    "Y": 0.3,  # 3% — aromatic/polar
    # L/I absent: aliphatic can't stack the same way
}
# Class I α2 Cys2+4: hydrophobic (L 86%)
# Property: buried hydrophobic in the Ig-fold core.
_GROOVE_I_CYS2_P4: dict[str, float] = {
    "L": 1.5,  # 86% — dominant hydrophobic
    "I": 0.5,  # 3% — similar hydrophobic
    "V": 0.5,  # 3% — smaller hydrophobic
    "F": 0.3,  # 1% — aromatic alternative
    # Charged/polar absent: buried position
}

# Class II β1 Cys1-1: charged, dominated by Glu (E 66%)
# Property: OPPOSITE of class I G — β1 enters the groove from a
# different structural context (solvent-exposed linker vs buried turn).
_GROOVE_II_CYS1_M1: dict[str, float] = {
    "E": 2.0,  # 66% — negative charge, exposed
    "D": 1.0,  # 8% — shorter negative charge
    "Q": 0.5,  # 5% — uncharged amide variant
    # G/A absent: too small for this exposed position in β1
}
# Class II β1 Cys1-5: weakly polar (Q 34%)
_GROOVE_II_CYS1_M5: dict[str, float] = {
    "Q": 1.0,  # 34% — modest conservation
    "K": 0.5,  # 12%
    "R": 0.3,  # 5%
}
# Class II β1 Cys1-3: positive charge (K 52%)
# Property: opposite of class I M — exposed linker position.
_GROOVE_II_CYS1_M3: dict[str, float] = {
    "K": 1.5,  # 52% — positive charge
    "R": 0.5,  # 8% — alternative positive
    "Q": 0.3,  # 5% — uncharged amide
    # M absent: β1 doesn't have the buried Met of α2
}
# Class II β1 Cys2+3: polar (N 81%)
# Property: opposite of class I W — this is the exposed face of β1.
_GROOVE_II_CYS2_P3: dict[str, float] = {
    "N": 2.0,  # 81% — small polar amide
    "D": 0.5,  # 5% — charged variant
    "S": 0.3,  # 3% — small polar
    # W/F absent: β1 Cys2+3 is exposed, not aromatic-stacking
}
# Class II β1 Cys2+4: aromatic (Y 81%)
_GROOVE_II_CYS2_P4: dict[str, float] = {
    "Y": 2.0,  # 81% — hydroxyl-aromatic
    "F": 0.5,  # 5% — aromatic alternative
    "H": 0.3,  # 3% — aromatic/charged
    # L absent: β1 uses aromatic, not aliphatic, at this position
}

# ---- Ig-domain Cys flanking scoring tables ----
#
# Universal Ig signals shared across class I α3, class II α2, and β2:
#   Cys2-2: Y (64-97% across all Ig types) — structural Tyr
#   Cys2+2: V (74-93%) — buried hydrophobic
#   W@Cys2+15: present in α3 (70%) and β2 (72%), absent in α Ig (28%)
#
# Class-I-specific α3 (n=253): distinctive LR-C-WA motif at Cys1
# Class-II β2 (n=501): Cys1 diverse, Cys2 conserved (Y..C..V..H)
# Class-II α Ig (n=757): Cys1-1 G/I, W@+15 mostly absent (L 60%)
#
# Scoring: universal Ig positions apply to all, class-specific
# positions add bonus when matched.

# Universal Ig: Cys2-2 Tyr (structural: packs against β-sheet in all Ig folds)
_IG_CYS2_M2: dict[str, float] = {
    "Y": 0.8,  # 64-97% across all Ig types
    "F": 0.3,  # aromatic alternative (without hydroxyl)
    # Non-aromatic absent: this position requires an aromatic ring
}
# Universal Ig: Cys2+2 Val (buried hydrophobic in all Ig folds)
_IG_CYS2_P2: dict[str, float] = {
    "V": 0.5,  # 74-93%
    "H": 0.3,  # 27% in α3 — histidine variant
    "I": 0.2,  # occasional aliphatic
}

# Class I α3 Cys1: LR-C-WA motif (highly distinctive)
# Property: the β-strand leading into α3 Cys1 has Leu-Arg before
# and Trp-Ala after — characteristic Ig C1-set domain signature.
_IG_I_CYS1_M2: dict[str, float] = {
    "L": 0.8,  # 92% — hydrophobic before Cys
    "I": 0.3,  # 3% — aliphatic alternative
    "V": 0.3,  # 2% — smaller aliphatic
}
_IG_I_CYS1_M1: dict[str, float] = {
    "R": 0.8,  # 78% — Arg before Cys, salt bridge
    "K": 0.3,  # 5% — positive alternative
    # Neutral/negative absent: this position requires a positive charge
}
_IG_I_CYS1_P1: dict[str, float] = {
    "W": 1.0,  # 78% — conserved Ig Trp right after Cys
    "F": 0.3,  # 5% — aromatic alternative
    # Aliphatic absent: must be aromatic at this Ig position
}
_IG_I_CYS1_P2: dict[str, float] = {
    "A": 0.8,  # 96% — nearly invariant small residue
    "S": 0.3,  # 2% — small/polar variant
    "G": 0.3,  # 1% — smallest
    # Bulky residues absent: tight packing requires small side chain
}

# Class II β2 Ig Cys1: more diverse than α3
_IG_II_CYS1_M1: dict[str, float] = {
    "V": 0.5,  # 29% — hydrophobic
    "G": 0.3,  # 23% — small
    "A": 0.3,  # 19% — small
    # Much more diverse than class I — no single dominant residue
}
# Class II β2 Cys2+4: His (73%)
_IG_II_CYS2_P4: dict[str, float] = {
    "H": 0.5,  # 73% — aromatic/charged
    "N": 0.3,  # 8% — polar alternative
}


def _score_groove_flanking(seq: str, c1: int, c2: int) -> tuple[float, list[str]]:
    """Score how well flanking residues match a groove-domain Cys pair.

    Scores against BOTH class I α2 and class II β1 groove patterns,
    taking the max at each position.  This serves double duty:
    1. Classifying the pair as groove vs Ig (combined with _score_ig_flanking)
    2. Inferring class I vs class II when the MHC class is unknown

    Class I α2 groove: G@-1, Q@-5, M@-3, W@+3, L@+4
    Class II β1 groove: E@-1, K@-3, N@+3, Y@+4
    """
    score = 0.0
    evidence: list[str] = []
    n = len(seq)

    # Separation: groove pairs cluster at ~63 (both class I and II)
    sep = c2 - c1
    sep_s = _gaussian_score(sep, 63, 4) * 2.0
    score += sep_s
    evidence.append(f"sep={sep}({sep_s:+.1f})")

    # Cys1 flanking — take max of class I and class II scores at each position
    for offset, tables_I, tables_II, label in [
        (-1, _GROOVE_I_CYS1_M1, _GROOVE_II_CYS1_M1, "c1-1"),
        (-5, _GROOVE_I_CYS1_M5, _GROOVE_II_CYS1_M5, "c1-5"),
        (-3, _GROOVE_I_CYS1_M3, _GROOVE_II_CYS1_M3, "c1-3"),
    ]:
        pos = c1 + offset
        if 0 <= pos < n:
            aa = seq[pos]
            s = max(tables_I.get(aa, 0.0), tables_II.get(aa, 0.0))
            if s:
                score += s
                evidence.append(f"{aa}@{label}({s:+.1f})")

    # Cys2 flanking — same dual-class scoring
    for offset, tables_I, tables_II, label in [
        (3, _GROOVE_I_CYS2_P3, _GROOVE_II_CYS2_P3, "c2+3"),
        (4, _GROOVE_I_CYS2_P4, _GROOVE_II_CYS2_P4, "c2+4"),
    ]:
        pos = c2 + offset
        if 0 <= pos < n:
            aa = seq[pos]
            s = max(tables_I.get(aa, 0.0), tables_II.get(aa, 0.0))
            if s:
                score += s
                evidence.append(f"{aa}@{label}({s:+.1f})")

    # Groove domains almost never have W at Cys2+15 (14% groove vs 51% Ig)
    if c2 + 15 < n and seq[c2 + 15] != "W":
        score += 1.0
        evidence.append("no_W@c2+15(+1.0)")

    return score, evidence


def _score_ig_flanking(seq: str, c1: int, c2: int) -> tuple[float, list[str]]:
    """Score how well flanking residues match an Ig-domain Cys pair.

    Scores universal Ig signals (shared across all Ig domain types) plus
    class-I-specific α3 signals.  Class II Ig domains (α2, β2) match on
    the universal positions; α3 additionally matches on its distinctive
    LR-C-WA motif.

    Universal Ig signals:
      Separation ~56, Cys2-2=Y (64-97%), Cys2+2=V (74-93%)
    Class I α3 bonus: LR-C-WA at Cys1 (L@-2 92%, R@-1 78%, W@+1 78%, A@+2 96%)
    Class II β2 bonus: Cys2+4=H (73%), Cys2+15=W (72%)
    W@Cys2+15: 51% overall Ig, 14% groove — useful but not definitive
    """
    score = 0.0
    evidence: list[str] = []
    n = len(seq)

    # Separation: Ig pairs cluster at ~56
    sep = c2 - c1
    sep_s = _gaussian_score(sep, 56, 3) * 2.0
    score += sep_s
    evidence.append(f"sep={sep}({sep_s:+.1f})")

    # W at Cys2+15 — moderate Ig signal (51% Ig vs 14% groove)
    if c2 + 15 < n:
        if seq[c2 + 15] == "W":
            score += 2.0
            evidence.append("W@c2+15(+2.0)")
        elif seq[c2 + 15] in "LF":
            score += 0.3
            evidence.append(f"{seq[c2+15]}@c2+15(+0.3)")

    # Universal Ig Cys2 flanking (shared across all Ig types)
    if c2 >= 2:
        s = _IG_CYS2_M2.get(seq[c2 - 2], 0.0)
        if s:
            score += s
            evidence.append(f"{seq[c2-2]}@c2-2({s:+.1f})")
    if c2 + 2 < n:
        s = _IG_CYS2_P2.get(seq[c2 + 2], 0.0)
        if s:
            score += s
            evidence.append(f"{seq[c2+2]}@c2+2({s:+.1f})")

    # Class I α3 Cys1 flanking: LR-C-WA (take max of class I and class II)
    if c1 >= 2:
        for tables_I, tables_II, offset, label in [
            (_IG_I_CYS1_M2, {}, -2, "c1-2"),
            (_IG_I_CYS1_M1, _IG_II_CYS1_M1, -1, "c1-1"),
        ]:
            pos = c1 + offset
            if 0 <= pos < n:
                aa = seq[pos]
                s = max(tables_I.get(aa, 0.0), tables_II.get(aa, 0.0) if tables_II else 0.0)
                if s:
                    score += s
                    evidence.append(f"{aa}@{label}({s:+.1f})")
    if c1 + 2 < n:
        for table, offset, label in [
            (_IG_I_CYS1_P1, 1, "c1+1"),
            (_IG_I_CYS1_P2, 2, "c1+2"),
        ]:
            s = table.get(seq[c1 + offset], 0.0)
            if s:
                score += s
                evidence.append(f"{seq[c1+offset]}@{label}({s:+.1f})")

    # Class II β2 Cys2+4 bonus
    if c2 + 4 < n:
        s = _IG_II_CYS2_P4.get(seq[c2 + 4], 0.0)
        if s:
            score += s
            evidence.append(f"{seq[c2+4]}@c2+4({s:+.1f})")

    return score, evidence


def classify_cys_pair(seq: str, c1: int, c2: int) -> CysPairAnnotation:
    """Classify a Cys-Cys pair as groove-domain or Ig-domain using flanking motifs.

    Returns a CysPairAnnotation with groove_score, ig_score, and domain_type.
    Does not use raw position ranges — classification is purely from local
    sequence context around the Cys pair.
    """
    seq = seq.upper()
    groove_score, groove_ev = _score_groove_flanking(seq, c1, c2)
    ig_score, ig_ev = _score_ig_flanking(seq, c1, c2)

    if groove_score > ig_score + 1.5:
        domain_type = "groove"
    elif ig_score > groove_score + 1.5:
        domain_type = "ig"
    else:
        domain_type = "ambiguous"

    evidence = tuple(groove_ev + ig_ev)
    return CysPairAnnotation(
        c1=c1,
        c2=c2,
        separation=c2 - c1,
        groove_score=groove_score,
        ig_score=ig_score,
        domain_type=domain_type,
        evidence=evidence,
    )


# Tiered scoring tables for the α1/α2 junction motif.
# The groove1/groove2 boundary in class I has strong conservation at
# positions 0 (G 77%), +2 (H 86%), +5 (Q 90%) relative to the first
# α2 residue.  Derived from n=794 class I exact-match entries.
_JUNCTION_POS0: dict[str, float] = {
    "G": 3.0, "A": 1.5, "S": 1.5, "T": 1.0, "D": 0.5, "V": 0.5,
}
_JUNCTION_POS1: dict[str, float] = {"S": 2.0, "T": 1.0, "V": 1.0, "H": 0.5}
_JUNCTION_POS2: dict[str, float] = {"H": 3.0, "K": 1.0, "R": 1.0, "L": 0.5, "Y": 0.3}
_JUNCTION_POS3: dict[str, float] = {"T": 1.5, "V": 1.5, "I": 1.0, "L": 1.0}
_JUNCTION_POS5: dict[str, float] = {"Q": 2.5, "N": 1.0, "E": 0.5, "P": 0.3}


def _score_junction(seq: str, pos: int) -> float:
    """Score a candidate α1/α2 junction position (class I groove1/groove2 boundary).

    pos is the 0-indexed position of the first α2 residue.  The junction
    motif is G..H..Q at positions 0, +2, +5 (strongest signals), with
    supporting residues at +1 and +3.
    """
    if pos < 0 or pos + 5 >= len(seq):
        return -999.0
    score = _JUNCTION_POS0.get(seq[pos], -0.5)
    score += _JUNCTION_POS2.get(seq[pos + 2], -1.0)
    score += _JUNCTION_POS5.get(seq[pos + 5], -0.5)
    score += _JUNCTION_POS1.get(seq[pos + 1], 0.0)
    score += _JUNCTION_POS3.get(seq[pos + 3], 0.0)
    return score


def _find_junction(seq: str, groove_c1: int) -> tuple[int, float]:
    """Find the α1/α2 junction position by scanning near the groove Cys1.

    The junction is structurally ~10 residues before the groove-domain
    Cys1.  Scans ±5 around that expected position and picks the
    highest-scoring candidate.

    Returns (junction_pos, junction_score).
    """
    expected = groove_c1 - 10
    best_pos = expected
    best_score = -999.0
    for pos in range(max(0, expected - 5), min(len(seq), expected + 6)):
        s = _score_junction(seq, pos)
        if s > best_score:
            best_score = s
            best_pos = pos
    return best_pos, best_score


# Groove2/Ig boundary scoring tables.
# These are class-specific because the structural linker differs:
#   Class I:   L at -5 (87%), P at +1 (82%)
#   Class II α: N at 0 (64%), PP at +2/+3 (76%/78%)
#   Class II β: R at -1 (83%), P at +3 (91%), V at +5 (89%)

_BOUNDARY_CLASS1_M5: dict[str, float] = {"L": 1.5, "I": 0.5, "V": 0.5}
_BOUNDARY_CLASS1_P1: dict[str, float] = {"P": 2.0}
_BOUNDARY_CLASS1_P3: dict[str, float] = {"V": 0.8, "I": 0.5, "L": 0.5}

_BOUNDARY_CLASS2A_0: dict[str, float] = {"N": 1.5, "D": 0.5, "S": 0.5}
_BOUNDARY_CLASS2A_P2: dict[str, float] = {"P": 1.5}
_BOUNDARY_CLASS2A_P3: dict[str, float] = {"P": 1.5}

_BOUNDARY_CLASS2B_M1: dict[str, float] = {"R": 2.0, "K": 1.0}
_BOUNDARY_CLASS2B_P3: dict[str, float] = {"P": 2.0}
_BOUNDARY_CLASS2B_P5: dict[str, float] = {"V": 1.5, "I": 0.5, "L": 0.5}


def _score_groove_ig_boundary(seq: str, pos: int, mhc_class: str, chain: str = "") -> float:
    """Score a candidate groove→Ig domain boundary.

    pos is the first residue of the Ig domain.  Class-specific because
    the structural linker entering the Ig support domain differs:
      Class I:   hydrophobic...Pro linker
      Class II α: Asn...Pro-Pro linker
      Class II β: Arg...Pro linker
    """
    n = len(seq)
    if pos < 5 or pos + 5 >= n:
        return 0.0
    score = 0.0
    if mhc_class == "I":
        score += _BOUNDARY_CLASS1_M5.get(seq[pos - 5], 0.0)
        score += _BOUNDARY_CLASS1_P1.get(seq[pos + 1], -0.5)
        score += _BOUNDARY_CLASS1_P3.get(seq[pos + 3], 0.0)
    elif chain == "alpha":
        score += _BOUNDARY_CLASS2A_0.get(seq[pos], 0.0)
        score += _BOUNDARY_CLASS2A_P2.get(seq[pos + 2], 0.0)
        score += _BOUNDARY_CLASS2A_P3.get(seq[pos + 3], 0.0)
    else:  # beta
        score += _BOUNDARY_CLASS2B_M1.get(seq[pos - 1], 0.0)
        score += _BOUNDARY_CLASS2B_P3.get(seq[pos + 3], -0.5)
        score += _BOUNDARY_CLASS2B_P5.get(seq[pos + 5], 0.0)
    return score


def _score_sp_candidate(seq: str, pos: int, junction_pos: int) -> float:
    """Score a candidate signal peptide cleavage site with junction confirmation.

    Combines von Heijne cleavage rules, hydrophobicity transition, mature
    start properties (H at +3), and structural consistency with the known
    α1/α2 junction position.

    Parameters
    ----------
    seq : str
        Full protein sequence (uppercase).
    pos : int
        Candidate mature protein start position.
    junction_pos : int
        Known or estimated α1/α2 junction position (from _find_junction).
        Used to verify that the implied α1 domain length is plausible.
    """
    if pos < 3 or pos >= len(seq) - 3:
        return -999.0

    score = 0.0
    n = len(seq)

    # Signal 1: von Heijne -1 cleavage residue
    m1 = seq[pos - 1]
    if m1 in _SP_CLEAVAGE_STRONG:  # A, G
        score += 3.0
    elif m1 == "S":
        score += 2.0
    elif m1 == "C":
        score += 1.5
    elif m1 == "T":
        score += 1.0
    elif m1 in _SP_CHARGED or m1 == "P":
        score -= 3.0
    else:
        score -= 0.5

    # Signal 2: von Heijne -3 residue
    if pos >= 3:
        if seq[pos - 3] in _SP_SMALL_ALIPHATIC:
            score += 1.0
        elif seq[pos - 3] in _SP_CHARGED:
            score -= 1.0

    # Signal 3: +1 not proline
    if seq[pos] == "P":
        score -= 2.0

    # Signal 4: upstream hydrophobic density (h-region)
    if pos >= 12:
        upstream = seq[pos - 12 : pos - 3]
        hydro_frac = sum(1 for c in upstream if c in _SP_HYDROPHOBIC) / max(len(upstream), 1)
        score += (hydro_frac - 0.4) * 2.0

    # Signal 5: c-region polarity
    if pos >= 5:
        c_hydro = sum(1 for c in seq[pos - 3 : pos] if c in _SP_HYDROPHOBIC) / 3
        if c_hydro > 0.66:
            score -= 1.5

    # Signal 6: H at +3 — conserved structural His in MHC α1 S1 strand
    if pos + 2 < n and seq[pos + 2] == "H":
        score += 2.0

    # Signal 7: first mature residue properties
    if pos < n:
        p1 = seq[pos]
        if p1 in _SP_CHARGED:
            score += 1.0
        elif p1 in "GAS":
            score += 0.5

    # Signal 8: KD hydrophobicity transition
    if pos >= 4 and pos + 3 <= n:
        pre_kd = sum(_KD_SCALE.get(seq[i], 0.0) for i in range(pos - 3, pos)) / 3
        post_kd = sum(_KD_SCALE.get(seq[i], 0.0) for i in range(pos, pos + 3)) / 3
        kd_drop = pre_kd - post_kd
        if kd_drop > 0.5:
            score += min(kd_drop * 0.5, 1.5)

    # Signal 9: α1 domain length plausibility (from junction position)
    if junction_pos > 0:
        alpha1_len = junction_pos - pos
        if 80 <= alpha1_len <= 95:
            score += 2.0  # standard range
        elif 75 <= alpha1_len <= 125:
            score += 0.5  # extended range (birds, fish, DMA)
        else:
            score -= 3.0  # very unlikely

    return score


def _find_mature_start(
    seq: str,
    groove_c1: int,
    *,
    species_category: str = "",
    mhc_class: str = "I",
) -> int:
    """Find the mature protein start using structural scoring.

    Replaces the constant-subtraction approach (mature_start = c1 - CONSTANT)
    with a scored search that combines SP cleavage signals with α1/α2 junction
    motif confirmation.

    For mammals: uses the traditional ±2 scan from the Cys-based prediction
    (safe, preserves existing accuracy).
    For non-mammals: wide scored search using junction discovery.

    Parameters
    ----------
    seq : str
        Full protein sequence (uppercase).
    groove_c1 : int
        Raw position of Cys1 of the groove-domain Cys pair.
    species_category : str
        Species category (used to decide mammal vs non-mammal path).
    mhc_class : str
        "I", "II" — affects junction motif search.
    """
    # Traditional estimate as fallback
    if mhc_class == "I":
        traditional = max(0, groove_c1 - 100)
    else:
        traditional = max(0, groove_c1 - 106)

    if not seq or groove_c1 <= 0:
        return traditional

    # For mammals: use the traditional estimate + ±2 refinement (safe)
    if species_category in _MAMMAL_CATEGORIES or not species_category:
        return traditional

    # For non-mammals: junction-based scored search
    # Step 1: Find the junction (α1/α2 boundary) near groove_c1 - 10
    if mhc_class == "I":
        junction_pos, junction_score = _find_junction(seq, groove_c1)
    else:
        # Class II doesn't have a groove1/groove2 junction on one chain
        # Use the Cys-based estimate for the groove-end position
        junction_pos = max(0, groove_c1 - 23)  # groove ends 23 before Ig Cys1

    # Step 2: Search for SP boundary
    # α1 domain is 75-125 aa across all vertebrates
    search_min = max(5, junction_pos - 125)
    search_max = min(junction_pos - 75, MAX_PLAUSIBLE_SP)

    if search_max <= search_min:
        # Can't find a valid range — probably SP-stripped
        return max(0, traditional)

    best_pos = traditional
    best_score = _score_sp_candidate(seq, traditional, junction_pos)

    for pos in range(search_min, search_max + 1):
        s = _score_sp_candidate(seq, pos, junction_pos)
        if s > best_score:
            best_score = s
            best_pos = pos

    return best_pos


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _clean_seq(sequence: Optional[str]) -> str:
    return "".join(ch for ch in str(sequence or "").strip().upper() if not ch.isspace())


def _infer_mature_start(cys1_raw: int, mature_pos: int) -> int:
    return max(0, int(cys1_raw) - int(mature_pos))


# Signal peptide cleavage site residues.
#
# The -1 position (last residue of the SP) is A/G/S/C across all jawed
# vertebrates: sharks 96%, mammals 98%, birds 100%, fish 88%, reptiles 84%,
# amphibians 95% (including C at 25%).
#
# Validated against UniProt Signal features (SignalP predictions) for
# thousands of MHC entries per clade (March 2026).
_SP_CLEAVAGE_RESIDUES = frozenset("AGSC")
_SP_CLEAVAGE_STRONG = frozenset("AG")
_SP_HYDROPHOBIC = frozenset("AILMFVW")
_SP_SMALL_ALIPHATIC = frozenset("AVTSIL")
_SP_CHARGED = frozenset("DEKR")

# Kyte-Doolittle hydrophobicity scale (used for cleavage-site transition
# detection).  Positive = hydrophobic, negative = hydrophilic.
_KD_SCALE: dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
    "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3,
    "P": -1.6, "H": -3.2, "D": -3.5, "E": -3.5, "N": -3.5, "Q": -3.5,
    "K": -3.9, "R": -4.5,
}

# Species categories considered mammalian (tighter SP length distribution,
# use simple ±2 scan which is 99.7–100% accurate).
_MAMMAL_CATEGORIES = frozenset({"human", "nhp", "murine", "ungulate", "carnivore", "other_mammal"})


def _score_nonmammal_candidate(seq: str, pos: int, cys_pred: int) -> float:
    """Score a candidate SP cleavage site for a non-mammalian sequence.

    Combines universal biochemical signals rather than species-specific
    motif tables:

    1. -1 cleavage residue — von Heijne rule (A/G strongest, S/T weaker)
    2. -3 residue — von Heijne rule (small/aliphatic preferred)
    3. +1 not proline
    4. Upstream hydrophobic density (h-region present?)
    5. c-region polarity (transition from hydrophobic to polar)
    6. Distance from Cys-pair prediction (gentle anchor)
    7. Mature protein start properties:
       a. H at +3 — conserved structural His in MHC alpha1 beta-strand S1,
          present in 33-52% of sequences from sharks to humans
       b. Charged residue (E/D/K/R) at +1 — dominant in birds (E 52%)
          and common in fish (K/E)
       c. Small neutral (G/A/S) at +1 — dominant in mammals (G 50%)
    8. Hydrophobicity transition — Kyte-Doolittle drop across the
       cleavage boundary, universal across all vertebrate MHC SPs
    """
    if pos < 3 or pos >= len(seq) - 3:
        return -999.0

    score = 0.0

    # Signal 1: -1 cleavage residue
    m1 = seq[pos - 1]
    if m1 in _SP_CLEAVAGE_STRONG:  # A, G
        score += 3.0
    elif m1 == "S":
        score += 2.0
    elif m1 == "C":
        score += 1.5
    elif m1 == "T":
        score += 1.0
    elif m1 in _SP_CHARGED or m1 == "P":
        score -= 3.0
    else:
        score -= 0.5

    # Signal 2: -3 residue
    if seq[pos - 3] in _SP_SMALL_ALIPHATIC:
        score += 1.0
    elif seq[pos - 3] in _SP_CHARGED:
        score -= 1.0

    # Signal 3: +1 not proline
    if seq[pos] == "P":
        score -= 2.0

    # Signal 4: upstream hydrophobic density (h-region)
    if pos >= 12:
        upstream = seq[pos - 12 : pos - 3]
        hydro_frac = sum(1 for c in upstream if c in _SP_HYDROPHOBIC) / 9
        score += (hydro_frac - 0.4) * 2.0

    # Signal 5: c-region polarity
    if pos >= 5:
        c_hydro = sum(1 for c in seq[pos - 3 : pos] if c in _SP_HYDROPHOBIC) / 3
        if c_hydro > 0.66:
            score -= 1.5  # still hydrophobic → probably inside core

    # Signal 6: distance from Cys prediction (gentle)
    score -= abs(pos - cys_pred) * 0.3

    # Signal 7: mature protein start properties
    if pos + 2 < len(seq):
        # 7a: H at +3 — conserved structural histidine in MHC alpha1
        # first beta-strand.  Present across all vertebrate clades
        # (mammals: GSH/GPH/SPH; birds: ELH/EPH; fish: KHS/VTH).
        if seq[pos + 2] == "H":
            score += 2.0
        # 7b: first mature residue properties
        p1 = seq[pos]
        if p1 in _SP_CHARGED:  # E/D/K/R — dominant in birds and fish
            score += 1.0
        elif p1 in "GAS":  # small neutral — dominant in mammals
            score += 0.5

    # Signal 8: hydrophobicity transition at cleavage boundary
    # SP ends with a polar c-region; mature protein typically starts
    # with charged/polar residues.  A sharp Kyte-Doolittle drop across
    # the boundary is a universal cleavage signal.
    if pos >= 4 and pos + 3 <= len(seq):
        pre_kd = sum(_KD_SCALE.get(seq[i], 0.0) for i in range(pos - 3, pos)) / 3
        post_kd = sum(_KD_SCALE.get(seq[i], 0.0) for i in range(pos, pos + 3)) / 3
        kd_drop = pre_kd - post_kd
        if kd_drop > 0.5:
            score += min(kd_drop * 0.5, 1.5)

    return score


def refine_signal_peptide(
    sequence: str,
    mature_start: int,
    species_category: str = "",
) -> int:
    """Refine signal peptide cleavage site using sequence features.

    For mammals: simple ±2 scan for nearest A/G/S/C at -1.  The Cys-pair
    prediction is precise enough that this achieves 99.7–100% valid
    cleavage sites.

    For non-mammals (birds, fish, reptiles, amphibians): property-based
    scoring over a ±8 window combining von Heijne cleavage rules,
    conserved MHC structural features (H at +3 in alpha1 S1 strand),
    hydrophobicity transition detection, and distance from Cys prediction.

    Parameters
    ----------
    sequence : str
        Full protein sequence (including signal peptide).
    mature_start : int
        Predicted mature protein start from Cys-pair heuristic.
    species_category : str
        One of the mhcseqs species categories.  Determines whether to use
        simple (mammal) or scored (non-mammal) refinement.
    """
    if mature_start <= 0 or not sequence:
        return mature_start

    seq = sequence.upper()
    if mature_start >= len(seq):
        return mature_start

    if species_category in _MAMMAL_CATEGORIES or not species_category:
        # Mammals: if already at a valid cleavage residue, keep it.
        # Otherwise simple ±2 scan (tight, safe).
        if seq[mature_start - 1] in _SP_CLEAVAGE_RESIDUES:
            return mature_start
        for delta in (1, -1, 2, -2):
            candidate = mature_start + delta
            if 5 <= candidate < len(seq) and seq[candidate - 1] in _SP_CLEAVAGE_RESIDUES:
                return candidate
        return mature_start

    # Non-mammals: multi-signal scoring over ±8 window, with junction
    # motif confirmation as an additional structural signal.
    #
    # The junction motif (G..H..Q at the α1/α2 boundary) helps disambiguate
    # competing SP candidates by verifying that the implied α1 domain length
    # places the junction at a structurally plausible position.
    #
    # First, try to find the groove Cys pair and junction position for use
    # as an additional signal in scoring.
    junction_pos = 0
    pairs = find_cys_pairs(seq)
    if pairs:
        annotated = [classify_cys_pair(seq, c1, c2) for c1, c2, _sep in pairs]
        groove_pairs = [a for a in annotated if a.domain_type in ("groove", "ambiguous")]
        if groove_pairs:
            groove = min(groove_pairs, key=lambda a: abs(a.c1 - (mature_start + 100)))
            junction_pos, _ = _find_junction(seq, groove.c1)

    # Two-pass strategy:
    # Pass 1: Standard ±8 property-based search (preserves existing accuracy)
    # Pass 2: If strong junction motif found, extend to ±15 but only accept
    #         candidates outside ±8 if junction-confirmed and better scoring.
    best_pos = mature_start
    best_score = -999.0
    for pos in range(max(5, mature_start - 8), min(len(seq) - 3, mature_start + 9)):
        s = _score_nonmammal_candidate(seq, pos, mature_start)
        if s > best_score:
            best_score = s
            best_pos = pos

    # Pass 2: try wider range if junction is strong
    junction_score = 0.0
    if junction_pos > 0:
        junction_score = _score_junction(seq, junction_pos)
    if junction_score >= 6.0:
        for pos in range(max(5, mature_start - 15), min(len(seq) - 3, mature_start + 16)):
            if mature_start - 8 <= pos <= mature_start + 8:
                continue  # already scored in pass 1
            s = _score_nonmammal_candidate(seq, pos, mature_start)
            # Junction confirmation required for extended range
            alpha1_len = junction_pos - pos
            if 80 <= alpha1_len <= 95:
                s += 2.5
            elif 75 <= alpha1_len <= 125:
                s += 1.0
            else:
                continue  # not junction-confirmed, skip
            if s > best_score:
                best_score = s
                best_pos = pos

    return best_pos


def _gene_prefix(gene: str) -> str:
    """Extract the 2–3 letter gene family prefix (e.g. 'DQA' from 'DQA1')."""
    token = str(gene or "").strip().upper()
    # Strip trailing digits: DQA1 → DQA, DRB3 → DRB
    while token and token[-1].isdigit():
        token = token[:-1]
    # Strip species prefix if present: HLA-DQA → DQA, BoLA-DQA → DQA
    if "-" in token:
        token = token.rsplit("-", 1)[-1]
    return token


def _is_non_classical_class_i(gene: str, allele: str) -> bool:
    """Detect non-classical MHC class I lineages (teleost L/S/P/H, contaminants)."""
    token = str(gene or allele or "").strip().upper()
    # Strip species prefix
    if "-" in token:
        token = token.rsplit("-", 1)[-1]
    return any(token.startswith(pat) for pat in NON_CLASSICAL_CLASS_I_GENE_PATTERNS)


def _refine_status(result: AlleleRecord) -> AlleleRecord:
    """Post-parse refinement: detect non-classical lineage and short grooves.

    Applied after a successful parse to override status when the gene is
    non-classical or the groove is too short to be functional.
    """
    if not result.ok:
        return result

    flags = list(result.flags)

    # Non-classical class I detection
    if result.mhc_class == "I" and _is_non_classical_class_i(result.gene, result.allele):
        flags.append("non_classical_lineage")
        return replace(result, status="non_classical", flags=_flags_to_tuple(flags))

    # Short groove detection: groove half too small for functional peptide binding
    g1 = result.groove1_len
    g2 = result.groove2_len
    if result.mhc_class == "I" and g1 > 0 and g1 < MIN_FUNCTIONAL_GROOVE_HALF_LEN:
        flags.append(f"groove1_short({g1})")
        return replace(result, status="short", flags=_flags_to_tuple(flags))
    if result.mhc_class == "II" and result.chain == "alpha" and g1 > 0 and g1 < MIN_FUNCTIONAL_GROOVE_HALF_LEN:
        flags.append(f"groove1_short({g1})")
        return replace(result, status="short", flags=_flags_to_tuple(flags))
    if result.mhc_class == "II" and result.chain == "beta" and g2 > 0 and g2 < MIN_FUNCTIONAL_GROOVE_HALF_LEN:
        flags.append(f"groove2_short({g2})")
        return replace(result, status="short", flags=_flags_to_tuple(flags))

    return result


def _class_ii_alpha_cys1_mature_pos(gene: str) -> int:
    """Return the gene-specific Ig Cys1 mature position for a class II alpha chain."""
    prefix = _gene_prefix(gene)
    return _CLASS_II_ALPHA_CYS1_MATURE_POS_BY_GENE.get(prefix, CLASS_II_ALPHA_IG_CYS1_MATURE_POS)


def _class_ii_beta2_cys1_mature_pos(gene: str) -> int:
    """Return the gene-specific β2 Ig Cys1 mature position for a class II beta chain."""
    prefix = _gene_prefix(gene)
    return _CLASS_II_BETA2_CYS1_MATURE_POS_BY_GENE.get(prefix, CLASS_II_BETA2_CYS1_MATURE_POS)


def _flags_to_tuple(flags: Sequence[str]) -> tuple[str, ...]:
    return tuple(str(flag) for flag in flags if flag)


def _slice_or_empty(seq: str, start: int, end: int) -> str:
    lo = max(0, int(start))
    hi = max(lo, min(len(seq), int(end)))
    return seq[lo:hi]


# ---------------------------------------------------------------------------
# Mutation parsing and application
# ---------------------------------------------------------------------------

_MUTATION_RE = re.compile(r"^([A-Y])(\d+)([A-Y])$", re.IGNORECASE)


def _parse_mutation(mut: object) -> tuple[int, str, str]:
    """Parse a mutation spec into (mature_pos, original_aa, mutant_aa).

    Accepts:
      - str: "K66A" (original + position + mutant, IEDB/mhcgnomes format)
      - tuple/list: (66, "K", "A") or (66, "A") — position, [original], mutant
      - mhcgnomes Mutation object: has .pos, .aa_original, .aa_mutant

    Positions are 1-indexed mature protein numbering.
    """
    if isinstance(mut, str):
        m = _MUTATION_RE.match(mut.strip())
        if not m:
            raise ValueError(f"Cannot parse mutation string: {mut!r} (expected format like 'K66A')")
        return (int(m.group(2)), m.group(1).upper(), m.group(3).upper())

    # mhcgnomes Mutation object or similar duck-typed object
    if hasattr(mut, "pos") and hasattr(mut, "aa_mutant"):
        return (int(mut.pos), str(getattr(mut, "aa_original", "") or "").upper(), str(mut.aa_mutant).upper())  # type: ignore[union-attr]

    if isinstance(mut, (tuple, list)):
        if len(mut) == 3:
            pos, orig, mutant = mut
            return (int(pos), str(orig).upper(), str(mutant).upper())
        if len(mut) == 2:
            pos, mutant = mut
            return (int(pos), "", str(mutant).upper())
        raise ValueError(f"Mutation tuple must have 2-3 elements, got {len(mut)}: {mut!r}")

    raise TypeError(f"Unsupported mutation type: {type(mut).__name__}")


def apply_mutations(result: AlleleRecord, mutations: Sequence[object]) -> AlleleRecord:
    """Apply mutations to a AlleleRecord, returning a new result with mutated sequences.

    Mutations are specified in mature protein numbering (1-indexed).
    Accepts any sequence of mutation specs: strings ("K66A"), tuples, or
    mhcgnomes Mutation objects.

    The original groove extraction must have succeeded (result.ok must be True).
    Mutations are applied to the mature sequence reconstructed from domain parts,
    then re-sliced into the same domain boundaries.

    >>> from mhcseqs.groove import parse_class_i, apply_mutations
    >>> wt = decompose_class_i("G" * 400, allele="test")  # doctest: +SKIP
    """
    if not result.ok:
        raise ValueError(f"Cannot apply mutations to failed result (status={result.status!r})")

    parsed = [_parse_mutation(m) for m in mutations]
    if not parsed:
        return result

    # Reconstruct mature sequence from domain parts
    if result.mhc_class == "I":
        mature = result.groove1 + result.groove2 + result.ig_domain + result.tail
    elif result.chain == "alpha":
        mature = result.groove1 + result.ig_domain + result.tail
    else:  # beta
        mature = result.groove2 + result.ig_domain + result.tail

    # Apply mutations
    mature_list = list(mature)
    mutation_strs: list[str] = []
    for pos, orig, mutant in parsed:
        idx = pos - 1  # convert to 0-indexed
        if idx < 0 or idx >= len(mature_list):
            raise ValueError(f"Mutation position {pos} out of range for mature sequence (length {len(mature_list)})")
        actual = mature_list[idx]
        if orig and actual != orig:
            raise ValueError(f"Mutation {orig}{pos}{mutant}: expected {orig} at mature position {pos}, found {actual}")
        mutation_strs.append(f"{actual}{pos}{mutant}")
        mature_list[idx] = mutant

    mutated = "".join(mature_list)

    # Re-slice into domain parts using original lengths
    g1_len = result.groove1_len
    g2_len = result.groove2_len
    ig_len = result.ig_domain_len

    if result.mhc_class == "I":
        g1 = mutated[:g1_len]
        g2 = mutated[g1_len : g1_len + g2_len]
        ig = mutated[g1_len + g2_len : g1_len + g2_len + ig_len]
        t = mutated[g1_len + g2_len + ig_len :]
    elif result.chain == "alpha":
        g1 = mutated[:g1_len]
        g2 = ""
        ig = mutated[g1_len : g1_len + ig_len]
        t = mutated[g1_len + ig_len :]
    else:  # beta
        g1 = ""
        g2 = mutated[:g2_len]
        ig = mutated[g2_len : g2_len + ig_len]
        t = mutated[g2_len + ig_len :]

    return replace(
        result,
        groove_seq=g1 + g2,
        groove1=g1,
        groove2=g2,
        ig_domain=ig,
        tail=t,
        tail_len=len(t),
        mutations=tuple(mutation_strs),
    )


# ---------------------------------------------------------------------------
# Cysteine pair detection
# ---------------------------------------------------------------------------


def find_cys_pairs(
    seq: str,
    min_sep: int = IG_SEP_MIN,
    max_sep: int = IG_SEP_MAX,
) -> list[tuple[int, int, int]]:
    """Find all Cys-Cys pairs with plausible Ig-fold separation.

    Returns list of (cys1_pos, cys2_pos, separation) tuples.
    """
    cleaned = _clean_seq(seq)
    cys_positions = [idx for idx, aa in enumerate(cleaned) if aa == "C"]
    pairs: list[tuple[int, int, int]] = []
    for i, c1 in enumerate(cys_positions):
        for c2 in cys_positions[i + 1 :]:
            sep = c2 - c1
            if sep < min_sep:
                continue
            if sep > max_sep:
                break
            pairs.append((c1, c2, sep))
    return pairs


# ---------------------------------------------------------------------------
# Class II chain inference
# ---------------------------------------------------------------------------


def _class_ii_chain_from_name(
    *,
    gene: str,
    allele: str,
) -> Optional[str]:
    gene_token = str(gene or "").strip().upper()
    if not gene_token and allele:
        try:
            gene_token = infer_gene(allele)
        except Exception:
            gene_token = ""
        gene_token = str(gene_token or "").strip().upper()

    if not gene_token:
        return None
    if gene_token.startswith(CLASS_II_ALPHA_GENE_PREFIXES):
        return "alpha"
    if gene_token.startswith(CLASS_II_BETA_GENE_PREFIXES):
        return "beta"
    if gene_token.endswith("A"):
        return "alpha"
    if gene_token.endswith("B"):
        return "beta"
    return None


def _class_i_fragment_result(
    *,
    seq: str,
    allele: str,
    gene: str,
) -> AlleleRecord:
    """Return an alpha1_only or alpha2_only record for class I single-exon fragments.

    Distinguishes which groove half the fragment represents by checking for
    the α2 Ig-fold Cys pair (separation ~55–72):
      - Present → exon 3 / α2 domain → groove2
      - Absent  → exon 2 / α1 domain → groove1
    """
    cleaned = _clean_seq(seq)
    pairs = find_cys_pairs(cleaned)
    has_alpha2_pair = any(55 <= sep <= 72 for _, _, sep in pairs)

    if has_alpha2_pair:
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="I",
            chain="alpha",
            seq_len=len(cleaned),
            mature_start=0,
            groove_seq=cleaned,
            groove1="",
            groove2=cleaned,
            groove1_len=0,
            groove2_len=len(cleaned),
            status="alpha2_only",
            anchor_type="raw_fragment",
            flags=("single_exon_fragment",),
        )
    return AlleleRecord(
        allele=allele,
        gene=gene,
        mhc_class="I",
        chain="alpha",
        seq_len=len(cleaned),
        mature_start=0,
        groove_seq=cleaned,
        groove1=cleaned,
        groove2="",
        groove1_len=len(cleaned),
        groove2_len=0,
        status="alpha1_only",
        anchor_type="raw_fragment",
        flags=("single_exon_fragment",),
    )


def _class_ii_fragment_result(
    *,
    seq: str,
    allele: str,
    gene: str,
    chain: str,
) -> AlleleRecord:
    cleaned = _clean_seq(seq)
    if chain == "alpha":
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="alpha",
            seq_len=len(cleaned),
            mature_start=0,
            groove_seq=cleaned,
            groove1=cleaned,
            groove2="",
            groove1_len=len(cleaned),
            groove2_len=0,
            status="fragment_fallback",
            anchor_type="raw_fragment",
            flags=("fragment_fallback",),
        )
    return AlleleRecord(
        allele=allele,
        gene=gene,
        mhc_class="II",
        chain="beta",
        seq_len=len(cleaned),
        mature_start=0,
        groove_seq=cleaned,
        groove1="",
        groove2=cleaned,
        groove1_len=0,
        groove2_len=len(cleaned),
        status="fragment_fallback",
        anchor_type="raw_fragment",
        flags=("fragment_fallback",),
    )


# ---------------------------------------------------------------------------
# Class I groove parser
# ---------------------------------------------------------------------------


def decompose_class_i(
    seq: str,
    *,
    allele: str = "",
    gene: str = "",
) -> AlleleRecord:
    """Parse a class-I alpha chain into alpha1/alpha2 groove halves.

    Primary strategy: locate the alpha2-domain Cys pair (raw pos ~60-180),
    infer signal peptide length from it, then slice alpha1 and alpha2 domains.

    Fallback: if no alpha2 pair, use the alpha3 Cys pair (raw pos ≥180) with
    fixed-width groove boundaries.
    """
    cleaned = _clean_seq(seq)
    flags: list[str] = []
    if len(cleaned) < MIN_GROOVE_SOURCE_LEN:
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="I",
            chain="alpha",
            seq_len=len(cleaned),
            status="too_short",
        )

    pairs = find_cys_pairs(cleaned)
    if not pairs:
        if len(cleaned) <= CLASS_I_FRAGMENT_MAX_LEN:
            return _class_i_fragment_result(seq=cleaned, allele=allele, gene=gene)
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="I",
            chain="alpha",
            seq_len=len(cleaned),
            status="no_cys_pairs",
        )

    # Classify all Cys pairs by flanking motifs (groove vs Ig) instead of
    # using hardcoded position ranges.  This handles non-mammalian sequences
    # where domain lengths differ from the human calibration.
    annotated = [classify_cys_pair(cleaned, c1, c2) for c1, c2, _sep in pairs]
    groove_classified = [a for a in annotated if a.domain_type in ("groove", "ambiguous")]
    ig_classified = [a for a in annotated if a.domain_type in ("ig", "ambiguous")]

    # Pick the best groove pair (alpha2): prefer highest groove_score
    alpha2_pair = None
    if groove_classified:
        best_groove = max(groove_classified, key=lambda a: a.groove_score)
        alpha2_pair = (best_groove.c1, best_groove.c2, best_groove.separation)

    if alpha2_pair is None:
        # No groove pair found — try Ig pair as alpha3 fallback
        alpha3_ann = max(ig_classified, key=lambda a: a.ig_score) if ig_classified else None
        alpha3_pair = (alpha3_ann.c1, alpha3_ann.c2, alpha3_ann.separation) if alpha3_ann else None
        if alpha3_pair is None:
            # Fragment fallback for short class I sequences
            if len(cleaned) <= CLASS_I_FRAGMENT_MAX_LEN:
                return _class_i_fragment_result(seq=cleaned, allele=allele, gene=gene)
            return AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class="I",
                chain="alpha",
                seq_len=len(cleaned),
                status="no_alpha2_pair",
            )

        c1, c2, _ = alpha3_pair
        mature_start = _infer_mature_start(c1, CLASS_I_ALPHA3_CYS1_MATURE_POS)
        if mature_start > MAX_PLAUSIBLE_SP:
            flags.append(f"suspect_mature_start({mature_start})")
            return AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class="I",
                chain="alpha",
                seq_len=len(cleaned),
                mature_start=mature_start,
                status="suspect_anchor",
                anchor_type="alpha3_cys",
                anchor_cys1=c1,
                anchor_cys2=c2,
                flags=_flags_to_tuple(flags),
            )
        alpha1_end = mature_start + DEFAULT_CLASS_I_GROOVE1_LEN
        alpha2_end = min(len(cleaned), c1 - 20)
        half_1 = _slice_or_empty(cleaned, mature_start, alpha1_end)
        half_2 = _slice_or_empty(cleaned, alpha1_end, alpha2_end)
        if not half_1 or not half_2:
            return AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class="I",
                chain="alpha",
                seq_len=len(cleaned),
                mature_start=mature_start,
                status="inferred_from_alpha3_bad_boundaries",
                anchor_type="alpha3_cys",
                anchor_cys1=c1,
                anchor_cys2=c2,
            )
        flags.append("inferred_from_alpha3")
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="I",
            chain="alpha",
            seq_len=len(cleaned),
            mature_start=mature_start,
            groove_seq=half_1 + half_2,
            groove1=half_1,
            groove2=half_2,
            groove1_len=len(half_1),
            groove2_len=len(half_2),
            status="inferred_from_alpha3",
            anchor_type="alpha3_cys",
            anchor_cys1=c1,
            anchor_cys2=c2,
            flags=_flags_to_tuple(flags),
        )

    # Primary alpha2 strategy
    c1, c2, _ = alpha2_pair
    mature_start = _infer_mature_start(c1, CLASS_I_ALPHA2_CYS1_MATURE_POS)

    # Fragment detection: if the groove pair is so close to the start that
    # there's no room for groove1 (α1 domain), treat as a fragment.
    # This handles SP-stripped single-exon fragments that have strong groove
    # flanking motifs but are too short for a full domain parse.
    if mature_start <= 0 and len(cleaned) <= CLASS_I_FRAGMENT_MAX_LEN:
        return _class_i_fragment_result(seq=cleaned, allele=allele, gene=gene)

    if mature_start > MAX_PLAUSIBLE_SP:
        # Cys mutation likely destroyed the real anchor; this pair is wrong
        flags.append(f"suspect_mature_start({mature_start})")
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="I",
            chain="alpha",
            seq_len=len(cleaned),
            mature_start=mature_start,
            status="suspect_anchor",
            anchor_type="alpha2_cys",
            anchor_cys1=c1,
            anchor_cys2=c2,
            flags=_flags_to_tuple(flags),
        )
    if mature_start > 40:
        flags.append(f"long_sp({mature_start})")

    alpha2_start = max(mature_start, c1 - CLASS_I_ALPHA2_CYS1_OFFSET)
    alpha2_end = min(len(cleaned), c2 + CLASS_I_ALPHA2_END_AFTER_CYS2)
    half_1 = _slice_or_empty(cleaned, mature_start, alpha2_start)
    half_2 = _slice_or_empty(cleaned, alpha2_start, alpha2_end)
    # Find secondary (alpha3 Ig) pair: prefer Ig-classified pairs downstream
    ig_downstream = [a for a in ig_classified if a.c1 > c2 + 10]
    if ig_downstream:
        best_ig = max(ig_downstream, key=lambda a: a.ig_score)
        secondary = (best_ig.c1, best_ig.c2, best_ig.separation)
    else:
        secondary = next(
            (pair for pair in pairs if pair[0] > c2 + 10),
            None,
        )

    # Extract Ig domain (α3) and tail if α3 Cys pair is available
    if secondary:
        s_c1, s_c2, _ = secondary
        ig_end = min(len(cleaned), s_c2 + IG_DOMAIN_END_AFTER_CYS2)
        ig = _slice_or_empty(cleaned, alpha2_end, ig_end)
        t = cleaned[ig_end:]
    else:
        ig = ""
        t = cleaned[alpha2_end:]

    if len(half_1) < 50:
        flags.append(f"alpha1_short({len(half_1)})")
    if len(half_2) < 60:
        flags.append(f"alpha2_short({len(half_2)})")
    if not half_1 or not half_2:
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="I",
            chain="alpha",
            seq_len=len(cleaned),
            mature_start=mature_start,
            status="invalid_boundaries",
            anchor_type="alpha2_cys",
            anchor_cys1=c1,
            anchor_cys2=c2,
        )
    return AlleleRecord(
        allele=allele,
        gene=gene,
        mhc_class="I",
        chain="alpha",
        seq_len=len(cleaned),
        mature_start=mature_start,
        groove_seq=half_1 + half_2,
        groove1=half_1,
        groove2=half_2,
        groove1_len=len(half_1),
        groove2_len=len(half_2),
        ig_domain=ig,
        ig_domain_len=len(ig),
        tail=t,
        tail_len=len(t),
        status="ok",
        anchor_type="alpha2_cys",
        anchor_cys1=c1,
        anchor_cys2=c2,
        secondary_cys1=(secondary[0] if secondary else None),
        secondary_cys2=(secondary[1] if secondary else None),
        flags=_flags_to_tuple(flags),
    )


def trace_parse_class_i(
    seq: str,
    *,
    allele: str = "",
    gene: str = "",
) -> tuple[AlleleRecord, ParseTrace]:
    """Parse class I with full debug trace of scoring decisions.

    Returns (AlleleRecord, ParseTrace) — the trace contains all Cys pair
    classifications, junction scores, SP candidate scores, and boundary
    motif scores for debugging ground truth disagreements.
    """
    cleaned = _clean_seq(seq)
    trace = ParseTrace(sequence_len=len(cleaned), mhc_class="I")

    # Step 1: find and classify Cys pairs
    pairs = find_cys_pairs(cleaned)
    if pairs:
        annotated = tuple(classify_cys_pair(cleaned, c1, c2) for c1, c2, _sep in pairs)
        trace.cys_pairs = annotated
        for i, ann in enumerate(annotated):
            trace.log.append(
                f"Cys #{i+1} C{ann.c1}-C{ann.c2} sep={ann.separation} "
                f"groove={ann.groove_score:.1f} ig={ann.ig_score:.1f} → {ann.domain_type}"
            )
            trace.log.append(f"  evidence: {' '.join(ann.evidence)}")

    # Step 2: find groove pair and junction
    groove_classified = [a for a in trace.cys_pairs if a.domain_type in ("groove", "ambiguous")]
    if groove_classified:
        best_groove = max(groove_classified, key=lambda a: a.groove_score)
        jpos, jscore = _find_junction(cleaned, best_groove.c1)
        trace.junction_pos = jpos
        trace.junction_score = jscore
        trace.log.append(
            f"Junction search near C{best_groove.c1}-10={best_groove.c1-10}: "
            f"best at {jpos} score={jscore:.1f}"
        )
        # Show top 3 junction candidates
        j_candidates = []
        for p in range(max(0, best_groove.c1 - 15), min(len(cleaned), best_groove.c1 - 5)):
            j_candidates.append((p, _score_junction(cleaned, p)))
        j_candidates.sort(key=lambda x: -x[1])
        for p, s in j_candidates[:3]:
            motif = cleaned[p:p+6] if p + 6 <= len(cleaned) else cleaned[p:]
            trace.log.append(f"  junction candidate pos={p} motif={motif} score={s:.1f}")

    # Step 3: run the actual parser
    result = decompose_class_i(cleaned, allele=allele, gene=gene)
    trace.selected_sp = result.mature_start

    # Step 4: log domain boundaries
    trace.log.append(f"Parse result: status={result.status} mature_start={result.mature_start}")
    if result.ok:
        trace.log.append(
            f"  [SP 0..{result.mature_start}] "
            f"[groove1 {result.mature_start}..{result.mature_start + result.groove1_len}] "
            f"[groove2 ..+{result.groove2_len}] "
            f"[Ig ..+{result.ig_domain_len}] "
            f"[tail ..+{result.tail_len}]"
        )

    # Step 5: score boundary motifs for the trace
    if result.ok and result.anchor_cys2 is not None:
        boundary_pos = result.anchor_cys2 + CLASS_I_ALPHA2_END_AFTER_CYS2
        bscore = _score_groove_ig_boundary(cleaned, boundary_pos, "I")
        trace.boundary_scores["groove2_ig"] = bscore
        trace.log.append(f"Boundary groove2|Ig @{boundary_pos}: score={bscore:.1f}")

    return result, trace


# ---------------------------------------------------------------------------
# Class II alpha groove parser
# ---------------------------------------------------------------------------


def decompose_class_ii_alpha(
    seq: str,
    *,
    allele: str = "",
    gene: str = "",
) -> AlleleRecord:
    """Parse a class-II alpha chain into the alpha1 groove half.

    Locates the alpha2-domain Ig-fold Cys pair to infer the boundary between
    the alpha1 groove domain and the alpha2 Ig domain.
    """
    cleaned = _clean_seq(seq)
    flags: list[str] = []
    if len(cleaned) < MIN_GROOVE_SOURCE_LEN:
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="alpha",
            seq_len=len(cleaned),
            status="too_short",
        )

    pairs = find_cys_pairs(cleaned)

    primary = [pair for pair in pairs if CLASS_II_ALPHA_CYS1_RAW_PRIMARY_MIN <= pair[0] <= CLASS_II_ALPHA_CYS1_RAW_MAX]
    candidates = primary or [
        pair for pair in pairs if CLASS_II_ALPHA_CYS1_RAW_MIN <= pair[0] <= CLASS_II_ALPHA_CYS1_RAW_MAX
    ]
    if not candidates:
        if len(cleaned) <= CLASS_II_ALPHA_FRAGMENT_MAX_LEN:
            return _class_ii_fragment_result(
                seq=cleaned,
                allele=allele,
                gene=gene,
                chain="alpha",
            )
        status = "no_cys_pairs" if not pairs else "no_anchor_pair"
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="alpha",
            seq_len=len(cleaned),
            status=status,
        )

    c1, c2, _ = min(candidates, key=lambda item: (abs(item[2] - 56), -item[0]))
    cys1_mature_pos = _class_ii_alpha_cys1_mature_pos(gene)
    mature_start = _infer_mature_start(c1, cys1_mature_pos)
    if mature_start > MAX_PLAUSIBLE_SP:
        flags.append(f"suspect_mature_start({mature_start})")
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="alpha",
            seq_len=len(cleaned),
            mature_start=mature_start,
            status="suspect_anchor",
            anchor_type="alpha2_cys",
            anchor_cys1=c1,
            anchor_cys2=c2,
            flags=_flags_to_tuple(flags),
        )
    groove_end = max(mature_start, c1 - CLASS_II_ALPHA_GROOVE_END_BEFORE_IG_CYS)
    half_1 = _slice_or_empty(cleaned, mature_start, groove_end)

    # Extract Ig domain (α2) and tail
    ig_end = min(len(cleaned), c2 + IG_DOMAIN_END_AFTER_CYS2)
    ig = _slice_or_empty(cleaned, groove_end, ig_end)
    t = cleaned[ig_end:]

    if len(half_1) < 60:
        flags.append(f"alpha1_short({len(half_1)})")
    if not half_1:
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="alpha",
            seq_len=len(cleaned),
            mature_start=mature_start,
            status="invalid_boundaries",
            anchor_type="alpha2_cys",
            anchor_cys1=c1,
            anchor_cys2=c2,
        )
    return AlleleRecord(
        allele=allele,
        gene=gene,
        mhc_class="II",
        chain="alpha",
        seq_len=len(cleaned),
        mature_start=mature_start,
        groove_seq=half_1,
        groove1=half_1,
        groove2="",
        groove1_len=len(half_1),
        groove2_len=0,
        ig_domain=ig,
        ig_domain_len=len(ig),
        tail=t,
        tail_len=len(t),
        status="ok",
        anchor_type="alpha2_cys",
        anchor_cys1=c1,
        anchor_cys2=c2,
        flags=_flags_to_tuple(flags),
    )


# ---------------------------------------------------------------------------
# Class II beta groove parser
# ---------------------------------------------------------------------------


def decompose_class_ii_beta(
    seq: str,
    *,
    allele: str = "",
    gene: str = "",
) -> AlleleRecord:
    """Parse a class-II beta chain into the beta1 groove half.

    Primary strategy: locate the beta2 Ig-fold Cys pair to infer where the
    beta1 groove domain ends. Fallback: use the beta1 Cys pair directly.
    """
    cleaned = _clean_seq(seq)
    flags: list[str] = []
    if len(cleaned) < MIN_GROOVE_SOURCE_LEN:
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="beta",
            seq_len=len(cleaned),
            status="too_short",
        )

    pairs = find_cys_pairs(cleaned)
    if not pairs:
        if len(cleaned) <= CLASS_II_BETA_FRAGMENT_MAX_LEN:
            return _class_ii_fragment_result(
                seq=cleaned,
                allele=allele,
                gene=gene,
                chain="beta",
            )
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="beta",
            seq_len=len(cleaned),
            status="no_cys_pairs",
        )

    # Class II beta: position-range pair selection (class-II-specific flanking
    # motifs differ substantially from class I — E@-1, K@-3 for β1 vs G@-1,
    # M@-3 for α2 — so motif-based classification requires class-II-calibrated
    # tables which are documented above but not yet used for pair selection to
    # avoid regressions).
    beta1_candidates = [pair for pair in pairs if CLASS_II_BETA1_CYS1_RAW_MIN <= pair[0] <= CLASS_II_BETA1_CYS1_RAW_MAX]
    beta2_candidates = [pair for pair in pairs if CLASS_II_BETA2_CYS1_RAW_MIN <= pair[0] <= CLASS_II_BETA2_CYS1_RAW_MAX]

    beta1_pair = min(beta1_candidates, key=lambda item: (abs(item[2] - 64), item[0])) if beta1_candidates else None
    downstream_beta2 = [pair for pair in beta2_candidates if beta1_pair is None or pair[0] > beta1_pair[1] + 10]
    beta2_pair = min(downstream_beta2, key=lambda item: (abs(item[2] - 56), item[0])) if downstream_beta2 else None

    if beta1_pair is None and beta2_pair is None:
        if len(cleaned) <= CLASS_II_BETA_FRAGMENT_MAX_LEN:
            return _class_ii_fragment_result(
                seq=cleaned,
                allele=allele,
                gene=gene,
                chain="beta",
            )
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="beta",
            seq_len=len(cleaned),
            status="no_anchor_pair",
        )

    if beta2_pair is not None:
        c1, c2, _ = beta2_pair
        cys1_mature_pos = _class_ii_beta2_cys1_mature_pos(gene)
        mature_start = _infer_mature_start(c1, cys1_mature_pos)
        if mature_start > MAX_PLAUSIBLE_SP:
            flags.append(f"suspect_mature_start({mature_start})")
            return AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class="II",
                chain="beta",
                seq_len=len(cleaned),
                mature_start=mature_start,
                status="suspect_anchor",
                anchor_type="beta2_cys",
                anchor_cys1=c1,
                anchor_cys2=c2,
                flags=_flags_to_tuple(flags),
            )
        groove_end = max(mature_start, c1 - CLASS_II_BETA_GROOVE_END_BEFORE_BETA2_CYS)
        status = "ok"
        anchor_type = "beta2_cys"
        anchor_pair = beta2_pair
        # Extract Ig domain (β2) and tail
        ig_end = min(len(cleaned), c2 + IG_DOMAIN_END_AFTER_CYS2)
        ig = _slice_or_empty(cleaned, groove_end, ig_end)
        t = cleaned[ig_end:]
    else:
        c1, c2, _ = beta1_pair  # type: ignore[misc]
        mature_start = _infer_mature_start(c1, CLASS_II_BETA1_CYS1_MATURE_POS)
        if mature_start > MAX_PLAUSIBLE_SP:
            flags.append(f"suspect_mature_start({mature_start})")
            return AlleleRecord(
                allele=allele,
                gene=gene,
                mhc_class="II",
                chain="beta",
                seq_len=len(cleaned),
                mature_start=mature_start,
                status="suspect_anchor",
                anchor_type="beta1_cys",
                anchor_cys1=c1,
                anchor_cys2=c2,
                flags=_flags_to_tuple(flags),
            )
        groove_end = min(len(cleaned), c2 + CLASS_II_BETA1_ONLY_END_AFTER_CYS2)
        status = "beta1_only_fallback"
        anchor_type = "beta1_cys"
        anchor_pair = beta1_pair
        flags.append("beta1_only_fallback")
        # No β2 pair → cannot define Ig domain boundaries
        ig = ""
        t = cleaned[groove_end:]

    half_2 = _slice_or_empty(cleaned, mature_start, groove_end)
    if len(half_2) < 70:
        flags.append(f"beta1_short({len(half_2)})")
    if not half_2:
        return AlleleRecord(
            allele=allele,
            gene=gene,
            mhc_class="II",
            chain="beta",
            seq_len=len(cleaned),
            mature_start=mature_start,
            status="invalid_boundaries",
            anchor_type=anchor_type,
            anchor_cys1=(anchor_pair[0] if anchor_pair else None),
            anchor_cys2=(anchor_pair[1] if anchor_pair else None),
        )
    return AlleleRecord(
        allele=allele,
        gene=gene,
        mhc_class="II",
        chain="beta",
        seq_len=len(cleaned),
        mature_start=mature_start,
        groove_seq=half_2,
        groove1="",
        groove2=half_2,
        groove1_len=0,
        groove2_len=len(half_2),
        ig_domain=ig,
        ig_domain_len=len(ig),
        tail=t,
        tail_len=len(t),
        status=status,
        anchor_type=anchor_type,
        anchor_cys1=(anchor_pair[0] if anchor_pair else None),
        anchor_cys2=(anchor_pair[1] if anchor_pair else None),
        secondary_cys1=(beta1_pair[0] if beta1_pair and beta2_pair is not None else None),
        secondary_cys2=(beta1_pair[1] if beta1_pair and beta2_pair is not None else None),
        flags=_flags_to_tuple(flags),
    )


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


def decompose_domains(
    seq: str,
    *,
    mhc_class: str,
    chain: Optional[str] = None,
    allele: str = "",
    gene: str = "",
    mutations: Sequence[object] = (),
) -> AlleleRecord:
    """Dispatch groove parsing by class and chain.

    If ``mutations`` is provided, the groove is first extracted from the
    wild-type sequence (preserving Cys-pair detection), then the mutations
    are applied to the result.  Mutations use mature protein numbering
    (1-indexed) and can be strings ("K66A"), tuples, or mhcgnomes Mutation
    objects.
    """
    nc = normalize_mhc_class(mhc_class)
    if nc == "I":
        result = decompose_class_i(seq, allele=allele, gene=gene)
        if mutations and result.ok:
            result = apply_mutations(result, mutations)
        return _refine_status(result)
    if nc != "II":
        raise ValueError(f"Unsupported MHC class: {mhc_class!r}")

    chain_token = str(chain or "").strip().lower()
    result: Optional[AlleleRecord] = None
    if chain_token in {"a", "alpha", "mhc_a"}:
        result = decompose_class_ii_alpha(seq, allele=allele, gene=gene)
    elif chain_token in {"b", "beta", "mhc_b"}:
        result = decompose_class_ii_beta(seq, allele=allele, gene=gene)
    elif chain_token:
        raise ValueError(f"Unsupported class-II chain token: {chain!r}")
    else:
        # Infer chain from gene name
        name_chain = _class_ii_chain_from_name(gene=gene, allele=allele)
        if name_chain == "alpha":
            result = decompose_class_ii_alpha(seq, allele=allele, gene=gene)
        elif name_chain == "beta":
            result = decompose_class_ii_beta(seq, allele=allele, gene=gene)
        else:
            # Last resort: try both parsers
            alpha_r = decompose_class_ii_alpha(seq, allele=allele, gene=gene)
            beta_r = decompose_class_ii_beta(seq, allele=allele, gene=gene)
            if alpha_r.ok and not beta_r.ok:
                result = alpha_r
            elif beta_r.ok and not alpha_r.ok:
                result = beta_r
            elif alpha_r.ok and beta_r.ok:
                result = AlleleRecord(
                    allele=allele,
                    gene=gene,
                    mhc_class="II",
                    chain="",
                    seq_len=len(_clean_seq(seq)),
                    status="ambiguous_chain",
                    anchor_type="chain_inference",
                    flags=_flags_to_tuple(
                        (
                            f"alpha_status={alpha_r.status}",
                            f"beta_status={beta_r.status}",
                        )
                    ),
                )
            else:
                result = AlleleRecord(
                    allele=allele,
                    gene=gene,
                    mhc_class="II",
                    chain="",
                    seq_len=len(_clean_seq(seq)),
                    status="chain_inference_failed",
                    anchor_type="chain_inference",
                    flags=_flags_to_tuple(
                        (
                            f"alpha_status={alpha_r.status}",
                            f"beta_status={beta_r.status}",
                        )
                    ),
                )

    if mutations and result.ok:
        result = apply_mutations(result, mutations)
    return _refine_status(result)


def is_class_ii_alpha_gene(gene: str) -> bool:
    """Whether a gene name corresponds to a class II alpha chain.

    NOTE: This function uses a trailing-letter heuristic as a fallback
    (gene ending in "A" → alpha). This can misclassify fish class I
    genes (UBA, UCA, etc.) as class II alpha. Callers should gate on
    mhc_class == "II" before using this function.
    """
    token = str(gene or "").strip().upper()
    return token.startswith(CLASS_II_ALPHA_GENE_PREFIXES) or token.endswith("A")


# Backward compatibility aliases
extract_groove = decompose_domains
parse_class_i = decompose_class_i
parse_class_ii_alpha = decompose_class_ii_alpha
parse_class_ii_beta = decompose_class_ii_beta
