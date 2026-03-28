"""Data tables and grammar specifications for MHC domain parsing.

This module is intentionally data-heavy. It holds the canonical parser
grammars, motif tables, priors, and residue classes that drive
``mhcseqs.domain_parsing`` so the executable parsing logic can stay focused
on feature extraction, scoring, and parse assembly.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GrammarWeightSpec:
    """Weights for the fixed structural evidence of one chain grammar."""

    groove_anchor_score: float = 0.0
    groove_anchor_g_topology: float = 0.0
    groove_anchor_between: float = 0.0
    groove_secondary_g_topology: float = 0.0
    groove_secondary_between: float = 0.0
    groove_secondary_end_offset: float = 0.0
    groove_junction: float = 0.0
    groove_boundary: float = 1.0
    ig_anchor_score: float = 0.0
    ig_anchor_c_topology: float = 0.0
    ig_anchor_between: float = 0.0
    ig_secondary_c_topology: float = 0.0
    ig_secondary_between: float = 0.0
    ig_length: float = 0.35
    beta1_start_weight: float = 0.0
    tail: float = 0.35


@dataclass(frozen=True)
class GrammarSpec:
    """Parameterized domain grammar for one MHC class/chain architecture."""

    mhc_class: str
    chain: str
    description: str
    anchor_type: str
    fragment_max_len: int
    groove_output: str
    groove_segments: tuple[tuple[str, str], ...]
    support_role: str
    length_mode: str
    groove_prior_key: tuple[str, str]
    groove_boundary_from_c1: int
    groove_end_from_c2: int = 0
    uses_primary_anchor_for_ig: bool = True
    allow_self_secondary: bool = True
    short_flags: tuple[tuple[str, int], ...] = ()
    missing_no_pairs_flags: tuple[str, ...] = ("no_cys_pairs",)
    missing_no_anchor_flags: tuple[str, ...] = ()
    weights: GrammarWeightSpec = GrammarWeightSpec()


AA_PROPERTY: dict[str, str] = {}
for _aa in "AG":
    AA_PROPERTY[_aa] = "small"
for _aa in "VILM":
    AA_PROPERTY[_aa] = "aliphatic"
for _aa in "FWY":
    AA_PROPERTY[_aa] = "aromatic"
for _aa in "ST":
    AA_PROPERTY[_aa] = "hydroxyl"
for _aa in "NQ":
    AA_PROPERTY[_aa] = "amide"
for _aa in "DE":
    AA_PROPERTY[_aa] = "acidic"
for _aa in "KRH":
    AA_PROPERTY[_aa] = "basic"
for _aa in "C":
    AA_PROPERTY[_aa] = "cysteine"
for _aa in "P":
    AA_PROPERTY[_aa] = "proline"


IG_SEP_MIN = 48
IG_SEP_MAX = 72

CLASS_I_ALPHA2_DOMINANT_SEP = 63
CLASS_II_ALPHA_IG_DOMINANT_SEP = 56
CLASS_II_BETA1_DOMINANT_SEP = 64
CLASS_II_BETA2_DOMINANT_SEP = 56

MIN_GROOVE_SOURCE_LEN = 70

CLASS_I_ALPHA2_CYS1_OFFSET = 10
CLASS_I_ALPHA2_END_AFTER_CYS2 = 20
CLASS_II_ALPHA_GROOVE_END_BEFORE_IG_CYS = 23
CLASS_II_BETA_GROOVE_END_BEFORE_BETA2_CYS = 23
IG_DOMAIN_END_AFTER_CYS2 = 20

CLASS_I_FRAGMENT_MAX_LEN = 200
CLASS_II_ALPHA_FRAGMENT_MAX_LEN = 110
CLASS_II_BETA_FRAGMENT_MAX_LEN = 120

MAX_PLAUSIBLE_SP = 50
PRIMARY_PARSE_CANDIDATE_KEEP = 2

NON_GROOVE_GENES = frozenset({"MICA", "MICB", "MIC1", "MIC2", "HFE", "B2M", "MR1"})
NON_CLASSICAL_CLASS_I_GENE_PATTERNS = (
    "MHC1L",
    "MHC1S",
    "MHC1P",
    "LLA",
    "LCA",
    "LDA",
    "LFA",
    "LGA",
    "LIA",
    "LJA",
    "MFSD",
)
NON_MHC_ACCESSIONS = frozenset(
    {
        "Q1LUQ4",
        "B0UYT5",
    }
)
MIN_FUNCTIONAL_GROOVE_HALF_LEN = 70

MHC_SP_PRESENT_PREFIX_TRIADS = frozenset(
    {
        "MGP",
        "MGS",
        "MAL",
        "MAV",
        "MRV",
        "MGG",
        "MGL",
        "MAP",
        "MAA",
        "MVC",
        "MKM",
        "MRP",
        "MLR",
        "MET",
        "MSW",
        "MGA",
        "MMV",
        "MGT",
        "MER",
        "MTS",
    }
)
MHC_MATURE_PREFIX_TRIADS = frozenset(
    {
        "GSH",
        "ELH",
        "EPH",
        "IKE",
        "VTH",
        "VIH",
        "GDT",
        "CSH",
        "SSS",
        "EEL",
        "VKV",
        "RDS",
        "VKH",
        "GKE",
        "EYT",
        "RAT",
        "ETV",
        "ASS",
        "IQF",
    }
)
MHC_BOUNDARY_MATURE_TRIADS = MHC_MATURE_PREFIX_TRIADS | frozenset({"TRP", "AEL", "EET"})

PARSE_CANDIDATE_COMPAT_WEIGHTS: dict[str, float] = {
    "leaderless_conflict": -4.0,
    "leaderless_supported": 2.0,
    "sp_supported": 1.5,
    "sp_weak": -2.5,
    "sp_estimate_close": 0.8,
    "sp_estimate_far": -1.5,
    "full_missing_support": -2.5,
    "support_missing_in_long_sequence": -2.0,
    "full_support_consistent": 0.8,
    "full_tm_consistent": 0.6,
    "fragment_in_long_sequence": -2.0,
    "fragment_in_short_sequence": 0.6,
    "non_classical_gene_context": -0.5,
}

PRIMARY_CANDIDATE_CONTRADICTION_WEIGHTS: dict[str, float] = {
    "primary_topology_scale": 0.35,
    "secondary_topology_scale": 0.30,
    "support_missing_scale": 0.035,
    "sp_zone_conflict_scale": 0.10,
    "structural_mass_divisor": 40.0,
}

PRIMARY_PARSE_EXPANSION_MARGIN = 8.0
PRIMARY_PARSE_LOW_CONFIDENCE = 0.45
SP_ESTIMATE_CANDIDATE_KEEP = 3

CLASS_II_ALPHA_GENE_PREFIXES = (
    "DRA",
    "DQA",
    "DPA",
    "DMA",
    "DOA",
    "DYA",
    "DNA",
    "DAA",
    "DBA",
    "DCA",
    "DDA",
    "DEA",
    "BLA",
)
CLASS_II_BETA_GENE_PREFIXES = (
    "DRB",
    "DQB",
    "DPB",
    "DMB",
    "DOB",
    "DYB",
    "DIB",
    "DAB",
    "DBB",
    "DCB",
    "DDB",
    "DEB",
    "BLB",
)

GROOVE_I_CYS1_M1: dict[str, float] = {"G": 2.0, "A": 1.0, "S": 1.0, "T": 0.5, "V": 0.3}
GROOVE_I_CYS1_M5: dict[str, float] = {"Q": 1.5, "N": 0.5, "E": 0.3}
GROOVE_I_CYS1_M3: dict[str, float] = {"M": 1.0, "L": 0.5, "I": 0.5, "F": 0.3, "V": 0.3}
GROOVE_I_CYS2_P3: dict[str, float] = {"W": 1.5, "F": 0.5, "Y": 0.3}
GROOVE_I_CYS2_P4: dict[str, float] = {"L": 1.5, "I": 0.5, "V": 0.5, "F": 0.3}

GROOVE_II_CYS1_M1: dict[str, float] = {"E": 2.0, "D": 1.0, "Q": 0.5}
GROOVE_II_CYS1_M5: dict[str, float] = {"Q": 1.0, "K": 0.5, "R": 0.3}
GROOVE_II_CYS1_M3: dict[str, float] = {"K": 1.5, "R": 0.5, "Q": 0.3}
GROOVE_II_CYS2_P3: dict[str, float] = {"N": 2.0, "D": 0.5, "S": 0.3}
GROOVE_II_CYS2_P4: dict[str, float] = {"Y": 2.0, "F": 0.5, "H": 0.3}

IG_CYS2_M2: dict[str, float] = {"Y": 0.8, "F": 0.3}
IG_CYS2_P2: dict[str, float] = {"V": 0.5, "H": 0.3, "I": 0.2}
IG_I_CYS1_M2: dict[str, float] = {"L": 0.8, "I": 0.3, "V": 0.3}
IG_I_CYS1_M1: dict[str, float] = {"R": 0.8, "K": 0.3}
IG_I_CYS1_P1: dict[str, float] = {"W": 1.0, "F": 0.3}
IG_I_CYS1_P2: dict[str, float] = {"A": 0.8, "S": 0.3, "G": 0.3}
IG_II_CYS1_M1: dict[str, float] = {"V": 0.5, "G": 0.3, "A": 0.3}
IG_II_CYS2_P4: dict[str, float] = {"H": 0.5, "N": 0.3}

SP_EXCLUDED_M3_M1: frozenset[tuple[str, str]] = frozenset(
    {
        ("aromatic", "aliphatic"),
        ("aromatic", "aromatic"),
        ("aromatic", "hydroxyl"),
        ("aromatic", "amide"),
        ("aromatic", "acidic"),
        ("aromatic", "basic"),
        ("aromatic", "cysteine"),
        ("aromatic", "proline"),
        ("acidic", "aliphatic"),
        ("acidic", "aromatic"),
        ("acidic", "hydroxyl"),
        ("acidic", "amide"),
        ("acidic", "cysteine"),
        ("acidic", "proline"),
        ("basic", "aliphatic"),
        ("basic", "aromatic"),
        ("basic", "amide"),
        ("basic", "acidic"),
        ("basic", "basic"),
        ("basic", "cysteine"),
        ("basic", "proline"),
        ("proline", "aliphatic"),
        ("proline", "aromatic"),
        ("proline", "hydroxyl"),
        ("proline", "amide"),
        ("proline", "acidic"),
        ("proline", "basic"),
        ("proline", "cysteine"),
        ("proline", "proline"),
        ("amide", "aromatic"),
        ("amide", "hydroxyl"),
        ("amide", "amide"),
        ("amide", "acidic"),
        ("amide", "cysteine"),
        ("amide", "proline"),
        ("cysteine", "aliphatic"),
        ("cysteine", "aromatic"),
        ("cysteine", "amide"),
        ("cysteine", "acidic"),
        ("cysteine", "basic"),
        ("cysteine", "proline"),
        ("hydroxyl", "aliphatic"),
        ("hydroxyl", "aromatic"),
        ("hydroxyl", "acidic"),
    }
)

CYS2_AFTER_CLIKE: dict[tuple[str, ...], float] = {
    ("basic", "aliphatic", "acidic"): 2.0,
    ("amide", "aliphatic", "acidic"): 1.5,
    ("basic", "aliphatic", "amide"): 1.0,
}
CYS2_AFTER_GDOMAIN: dict[tuple[str, ...], float] = {
    ("aliphatic", "acidic", "aromatic"): 2.0,
    ("basic", "basic", "amide"): 1.5,
}
CYS1_BEFORE_CLIKE: dict[tuple[str, ...], float] = {
    ("aliphatic", "aliphatic", "aliphatic"): 1.0,
    ("hydroxyl", "aliphatic", "aliphatic"): 0.8,
}
CYS1_BEFORE_GDOMAIN: dict[tuple[str, ...], float] = {
    ("aliphatic", "aromatic", "small"): 1.0,
}

SP_H_REGION_WEIGHT = 0.6
SP_LATE_DRIFT_SCALE = 1.2
SP_EARLIER_TIE_MARGIN = 0.25

JUNCTION_POS0: dict[str, float] = {"G": 3.0, "A": 1.5, "S": 1.5, "T": 1.0, "D": 0.5, "V": 0.5}
JUNCTION_POS1: dict[str, float] = {"S": 2.0, "T": 1.0, "V": 1.0, "H": 0.5}
JUNCTION_POS2: dict[str, float] = {"H": 3.0, "K": 1.0, "R": 1.0, "L": 0.5, "Y": 0.3}
JUNCTION_POS3: dict[str, float] = {"T": 1.5, "V": 1.5, "I": 1.0, "L": 1.0}
JUNCTION_POS5: dict[str, float] = {"Q": 2.5, "N": 1.0, "E": 0.5, "P": 0.3}

BOUNDARY_CLASS1_M5: dict[str, float] = {"L": 1.5, "I": 0.5, "V": 0.5}
BOUNDARY_CLASS1_P1: dict[str, float] = {"P": 2.0}
BOUNDARY_CLASS1_P3: dict[str, float] = {"V": 0.8, "I": 0.5, "L": 0.5}
BOUNDARY_CLASS2A_0: dict[str, float] = {"N": 1.5, "D": 0.5, "S": 0.5}
BOUNDARY_CLASS2A_P2: dict[str, float] = {"P": 1.5}
BOUNDARY_CLASS2A_P3: dict[str, float] = {"P": 1.5}
BOUNDARY_CLASS2B_M1: dict[str, float] = {"R": 2.0, "K": 1.0}
BOUNDARY_CLASS2B_P3: dict[str, float] = {"P": 2.0}
BOUNDARY_CLASS2B_P5: dict[str, float] = {"V": 1.5, "I": 0.5, "L": 0.5}

GROOVE_PRIOR: dict[tuple[str, str], tuple[int, tuple[int, int], tuple[int, int]]] = {
    ("I", "alpha"): (90, (80, 106), (70, 120)),
    ("II", "alpha"): (80, (70, 100), (60, 115)),
    ("II", "beta"): (90, (70, 100), (60, 115)),
}
GROOVE_SEGMENT_PRIOR: dict[tuple[str, str], tuple[int, tuple[int, int], tuple[int, int]]] = {
    ("I", "groove1"): (90, (80, 100), (70, 112)),
    ("I", "groove2"): (93, (84, 104), (75, 114)),
    ("I", "total"): (183, (170, 198), (160, 214)),
    ("II", "alpha"): (80, (70, 100), (60, 115)),
    ("II", "beta"): (90, (70, 100), (60, 115)),
}
IG_PRIOR = (76, (65, 90), (55, 105))
SP_PRIOR = (23, (14, 42), (5, 50))
RELATIVE_OFFSET_PRIOR: dict[tuple[str, str], tuple[int, tuple[int, int], tuple[int, int]]] = {
    ("I", "alpha3_start_to_c1"): (18, (15, 22), (12, 26)),
    ("II", "beta1_start_to_c1"): (13, (10, 16), (6, 22)),
    ("II", "beta1_c2_to_end"): (15, (12, 18), (8, 24)),
}

GRAMMAR_SPECS: dict[tuple[str, str], GrammarSpec] = {
    ("I", "alpha"): GrammarSpec(
        mhc_class="I",
        chain="alpha",
        description="SP? -> G-alpha1 -> G-alpha2 -> C-like alpha3 -> TM? -> tail?",
        anchor_type="alpha2_cys",
        fragment_max_len=CLASS_I_FRAGMENT_MAX_LEN,
        groove_output="split_both",
        groove_segments=(("groove1", "g_alpha1"), ("groove2", "g_alpha2")),
        support_role="c1_alpha3",
        length_mode="class_i_split",
        groove_prior_key=("I", "alpha"),
        groove_boundary_from_c1=CLASS_I_ALPHA2_CYS1_OFFSET,
        groove_end_from_c2=CLASS_I_ALPHA2_END_AFTER_CYS2,
        uses_primary_anchor_for_ig=False,
        short_flags=(("alpha1_short", 50), ("alpha2_short", 60)),
        missing_no_anchor_flags=("no_plausible_groove_anchor", "no_plausible_support_anchor"),
        weights=GrammarWeightSpec(
            groove_anchor_score=0.45,
            groove_anchor_g_topology=0.85,
            groove_anchor_between=0.60,
            groove_junction=0.70,
            groove_boundary=1.0,
            ig_anchor_score=0.20,
            ig_secondary_c_topology=0.65,
            ig_secondary_between=0.20,
            ig_length=0.35,
            tail=0.35,
        ),
    ),
    ("II", "alpha"): GrammarSpec(
        mhc_class="II",
        chain="alpha",
        description="SP? -> G-alpha1 -> C-like alpha2 -> TM? -> tail?",
        anchor_type="alpha2_cys",
        fragment_max_len=CLASS_II_ALPHA_FRAGMENT_MAX_LEN,
        groove_output="groove1",
        groove_segments=(("groove1", "g_alpha1"),),
        support_role="c1_alpha2",
        length_mode="single_segment",
        groove_prior_key=("II", "alpha"),
        groove_boundary_from_c1=CLASS_II_ALPHA_GROOVE_END_BEFORE_IG_CYS,
        groove_end_from_c2=0,
        uses_primary_anchor_for_ig=True,
        short_flags=(("alpha1_short", 60),),
        missing_no_anchor_flags=("no_plausible_support_anchor",),
        weights=GrammarWeightSpec(
            groove_boundary=1.0,
            ig_anchor_score=0.35,
            ig_anchor_c_topology=0.80,
            ig_anchor_between=0.35,
            ig_length=0.35,
            tail=0.35,
        ),
    ),
    ("II", "beta"): GrammarSpec(
        mhc_class="II",
        chain="beta",
        description="SP? -> G-beta1 -> C-like beta2 -> TM? -> tail?",
        anchor_type="beta2_cys",
        fragment_max_len=CLASS_II_BETA_FRAGMENT_MAX_LEN,
        groove_output="groove2",
        groove_segments=(("groove2", "g_beta1"),),
        support_role="c1_beta2",
        length_mode="single_segment",
        groove_prior_key=("II", "beta"),
        groove_boundary_from_c1=CLASS_II_BETA_GROOVE_END_BEFORE_BETA2_CYS,
        groove_end_from_c2=0,
        uses_primary_anchor_for_ig=True,
        allow_self_secondary=False,
        short_flags=(("beta1_short", 70),),
        missing_no_anchor_flags=("no_plausible_anchor_pair",),
        weights=GrammarWeightSpec(
            groove_secondary_g_topology=0.80,
            groove_secondary_between=0.25,
            groove_secondary_end_offset=0.65,
            groove_boundary=1.0,
            ig_anchor_score=0.35,
            ig_anchor_c_topology=0.80,
            ig_anchor_between=0.35,
            ig_length=0.35,
            beta1_start_weight=0.45,
            tail=0.35,
        ),
    ),
}

MATURE_START_CLASS_I: dict[int, dict[str, float]] = {
    0: {"G": 1.0, "S": 0.8, "A": 0.6, "C": 0.4, "V": 0.3, "L": 0.3, "E": 0.2},
    1: {"S": 0.8, "P": 0.5, "L": 0.5, "T": 0.3, "H": 0.3},
    2: {"H": 2.5},
    3: {"S": 0.8, "T": 0.5, "L": 0.5, "M": 0.3},
    4: {"L": 0.5, "R": 0.5, "M": 0.3},
}
MATURE_START_CLASS_II: dict[int, dict[str, float]] = {
    0: {"R": 0.8, "D": 0.5, "E": 0.5, "K": 0.5, "I": 0.3, "F": 0.3},
    1: {"F": 0.5, "L": 0.3, "A": 0.3, "P": 0.3},
    2: {"L": 0.3, "Q": 0.3, "Y": 0.3},
}

REFINEMENT_MOTIFS: dict[str, dict[str, float]] = {
    "mammal": {"GSH": 1.5, "CSH": 1.0, "GPH": 1.0, "SPH": 0.8},
    "bird": {
        "AAE": 1.5,
        "ASE": 1.0,
        "ASG": 0.8,
        "ELH": 2.5,
        "EPH": 2.5,
        "EEL": 1.5,
        "EET": 1.5,
        "AEL": 1.5,
        "TRP": 1.5,
    },
    "fish": {"AVT": 1.0, "KHS": 1.0, "VKV": 1.5, "VDI": 1.0, "QPE": 1.5},
    "other_vertebrate": {"EEA": 1.0, "KHN": 1.0},
}
REFINEMENT_WINDOW_BY_GROUP: dict[str, int] = {
    "mammal": 5,
    "bird": 8,
    "fish": 8,
    "other_vertebrate": 10,
}
REFINEMENT_DISTANCE_PENALTY: dict[str, float] = {
    "mammal": 0.4,
    "bird": 0.6,
    "fish": 0.4,
    "other_vertebrate": 0.4,
}
SP_BOUNDARY_MODEL_PATH = Path(__file__).resolve().parent.parent / "data" / "sp_boundary_model.json"
SP_SEQUENCE_CUE_MODEL_PATH = Path(__file__).resolve().parent.parent / "data" / "sp_sequence_cue_model.json"
SP_BOUNDARY_MODEL_WEIGHT_BY_GROUP: dict[str, float] = {
    "mammal": 0.10,
    "bird": 0.45,
    "fish": 0.45,
    "other_vertebrate": 0.45,
}
SP_BOUNDARY_CLASS_WEIGHT = 0.30
SP_BOUNDARY_AA_CLASS: dict[str, str] = {
    "A": "small_aliphatic",
    "V": "small_aliphatic",
    "I": "aliphatic",
    "L": "aliphatic",
    "M": "aliphatic",
    "F": "aromatic",
    "W": "aromatic",
    "Y": "aromatic",
    "S": "polar",
    "T": "polar",
    "N": "polar",
    "Q": "polar",
    "D": "acidic",
    "E": "acidic",
    "K": "basic",
    "R": "basic",
    "H": "basic",
    "C": "cysteine",
    "G": "glycine",
    "P": "proline",
}

TM_HYDROPHOBIC = frozenset("AILMFVWCGY")
TM_CHARGED = frozenset("DEKRH")

SP_CLEAVAGE_RESIDUES = frozenset("AGSC")
SP_CLEAVAGE_STRONG = frozenset("AG")
SP_HYDROPHOBIC = frozenset("AILMFVW")
SP_SMALL_ALIPHATIC = frozenset("AVTSIL")
SP_CHARGED = frozenset("DEKR")
KD_SCALE: dict[str, float] = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "D": -3.5,
    "E": -3.5,
    "N": -3.5,
    "Q": -3.5,
    "K": -3.9,
    "R": -4.5,
}
MAMMAL_CATEGORIES = frozenset({"human", "nhp", "murine", "ungulate", "carnivore", "other_mammal"})
