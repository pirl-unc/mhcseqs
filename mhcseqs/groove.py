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
    "MHC1L", "MHC1S", "MHC1P",  # zebrafish-style: mhc1laa, mhc1saa, etc.
    "LLA", "LCA", "LDA", "LFA", "LGA", "LIA", "LJA",  # L lineage locus names
    "MFSD",  # NOT MHC at all — lipid transporter contaminant
)

# Known non-MHC proteins that have leaked into curated datasets via automated
# genome annotation.  Keyed by UniProt accession.
NON_MHC_ACCESSIONS = frozenset({
    "Q1LUQ4",   # Dare-mfsd6a, zebrafish lipid transporter
    "B0UYT5",   # Dare-mfsd6b, zebrafish lipid transporter
})

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
            "alpha3_fallback",
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
# Helpers
# ---------------------------------------------------------------------------


def _clean_seq(sequence: Optional[str]) -> str:
    return "".join(ch for ch in str(sequence or "").strip().upper() if not ch.isspace())


def _infer_mature_start(cys1_raw: int, mature_pos: int) -> int:
    return max(0, int(cys1_raw) - int(mature_pos))


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
    >>> wt = parse_class_i("G" * 400, allele="test")  # doctest: +SKIP
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
    """Return a fragment_fallback record for short class I sequences."""
    cleaned = _clean_seq(seq)
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
        status="fragment_fallback",
        anchor_type="raw_fragment",
        flags=("fragment_fallback",),
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


def parse_class_i(
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

    # Try alpha2 Cys pair first.
    # Prefer the pair whose separation is closest to the dominant value (63).
    # This avoids picking polymorphic Cys mutations (e.g. Y→C at pos N-2
    # creating a CxC motif) and intra-domain disulfides in non-classical
    # genes (e.g. HLA-G's α1 C42–C100 pair at sep=59).
    alpha2_candidates = [
        pair for pair in pairs if CLASS_I_ALPHA2_CYS1_RAW_MIN <= pair[0] <= CLASS_I_ALPHA2_CYS1_RAW_MAX
    ]
    alpha2_pair = (
        min(alpha2_candidates, key=lambda p: (abs(p[2] - CLASS_I_ALPHA2_DOMINANT_SEP), -p[0]))
        if alpha2_candidates
        else None
    )

    if alpha2_pair is None:
        # Fallback to alpha3 Cys pair
        alpha3_pair = next(
            (pair for pair in pairs if pair[0] >= CLASS_I_ALPHA3_CYS1_RAW_MIN),
            None,
        )
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
                status="alpha3_fallback_bad_boundaries",
                anchor_type="alpha3_cys",
                anchor_cys1=c1,
                anchor_cys2=c2,
            )
        flags.append("alpha3_fallback")
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
            status="alpha3_fallback",
            anchor_type="alpha3_cys",
            anchor_cys1=c1,
            anchor_cys2=c2,
            flags=_flags_to_tuple(flags),
        )

    # Primary alpha2 strategy
    c1, c2, _ = alpha2_pair
    mature_start = _infer_mature_start(c1, CLASS_I_ALPHA2_CYS1_MATURE_POS)
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


# ---------------------------------------------------------------------------
# Class II alpha groove parser
# ---------------------------------------------------------------------------


def parse_class_ii_alpha(
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
        # Fragment fallback: short sequences are likely just the groove domain
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


def parse_class_ii_beta(
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


def extract_groove(
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
        result = parse_class_i(seq, allele=allele, gene=gene)
        if mutations and result.ok:
            result = apply_mutations(result, mutations)
        return _refine_status(result)
    if nc != "II":
        raise ValueError(f"Unsupported MHC class: {mhc_class!r}")

    chain_token = str(chain or "").strip().lower()
    result: Optional[AlleleRecord] = None
    if chain_token in {"a", "alpha", "mhc_a"}:
        result = parse_class_ii_alpha(seq, allele=allele, gene=gene)
    elif chain_token in {"b", "beta", "mhc_b"}:
        result = parse_class_ii_beta(seq, allele=allele, gene=gene)
    elif chain_token:
        raise ValueError(f"Unsupported class-II chain token: {chain!r}")
    else:
        # Infer chain from gene name
        name_chain = _class_ii_chain_from_name(gene=gene, allele=allele)
        if name_chain == "alpha":
            result = parse_class_ii_alpha(seq, allele=allele, gene=gene)
        elif name_chain == "beta":
            result = parse_class_ii_beta(seq, allele=allele, gene=gene)
        else:
            # Last resort: try both parsers
            alpha_r = parse_class_ii_alpha(seq, allele=allele, gene=gene)
            beta_r = parse_class_ii_beta(seq, allele=allele, gene=gene)
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
