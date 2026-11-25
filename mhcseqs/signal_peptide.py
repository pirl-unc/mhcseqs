"""Signal peptide detection, sequence-cue/boundary models, refinement, and sequence analysis."""

from __future__ import annotations

import json
from functools import lru_cache
from typing import Optional, Sequence

from .domain_features import (
    _find_tm_span,
    _score_between_cys,
    _score_pair_topology_support,
    classify_cys_pair,
    find_cys_pairs,
    find_junction,
    score_junction,
)
from .domain_grammar import (
    AA_PROPERTY as _AA_PROPERTY,
)
from .domain_grammar import (
    KD_SCALE as _KD_SCALE,
)
from .domain_grammar import (
    MAMMAL_CATEGORIES as _MAMMAL_CATEGORIES,
)
from .domain_grammar import (
    MATURE_START_CLASS_I as _MATURE_START_CLASS_I,
)
from .domain_grammar import (
    MATURE_START_CLASS_II as _MATURE_START_CLASS_II,
)
from .domain_grammar import (
    MAX_PLAUSIBLE_SP,
    PRIMARY_PARSE_LOW_CONFIDENCE,
    SP_BOUNDARY_LOG_ODDS,
    SP_ESTIMATE_CANDIDATE_KEEP,
)
from .domain_grammar import (
    MHC_BOUNDARY_MATURE_TRIADS as _MHC_BOUNDARY_MATURE_TRIADS,
)
from .domain_grammar import (
    MHC_MATURE_PREFIX_TRIADS as _MHC_MATURE_PREFIX_TRIADS,
)
from .domain_grammar import (
    MHC_SP_PRESENT_PREFIX_TRIADS as _MHC_SP_PRESENT_PREFIX_TRIADS,
)
from .domain_grammar import (
    REFINEMENT_DISTANCE_PENALTY as _REFINEMENT_DISTANCE_PENALTY,
)
from .domain_grammar import (
    REFINEMENT_MOTIFS as _REFINEMENT_MOTIFS,
)
from .domain_grammar import (
    REFINEMENT_WINDOW_BY_GROUP as _REFINEMENT_WINDOW_BY_GROUP,
)
from .domain_grammar import (
    SP_BOUNDARY_AA_CLASS as _SP_BOUNDARY_AA_CLASS,
)
from .domain_grammar import (
    SP_BOUNDARY_CLASS_WEIGHT as _SP_BOUNDARY_CLASS_WEIGHT,
)
from .domain_grammar import (
    SP_BOUNDARY_MODEL_PATH as _SP_BOUNDARY_MODEL_PATH,
)
from .domain_grammar import (
    SP_BOUNDARY_MODEL_WEIGHT_BY_GROUP as _SP_BOUNDARY_MODEL_WEIGHT_BY_GROUP,
)
from .domain_grammar import (
    SP_CHARGED as _SP_CHARGED,
)
from .domain_grammar import (
    SP_CLEAVAGE_STRONG as _SP_CLEAVAGE_STRONG,
)
from .domain_grammar import (
    SP_EARLIER_TIE_MARGIN as _SP_EARLIER_TIE_MARGIN,
)
from .domain_grammar import (
    SP_EXCLUDED_M3_M1 as _SP_EXCLUDED_M3_M1,
)
from .domain_grammar import (
    SP_H_REGION_WEIGHT as _SP_H_REGION_WEIGHT,
)
from .domain_grammar import (
    SP_HYDROPHOBIC as _SP_HYDROPHOBIC,
)
from .domain_grammar import (
    SP_LATE_DRIFT_SCALE as _SP_LATE_DRIFT_SCALE,
)
from .domain_grammar import (
    SP_PRIOR as _SP_PRIOR,
)
from .domain_grammar import (
    SP_SEQUENCE_CUE_MODEL_PATH as _SP_SEQUENCE_CUE_MODEL_PATH,
)
from .domain_grammar import (
    SP_SMALL_ALIPHATIC as _SP_SMALL_ALIPHATIC,
)
from .domain_types import (
    SequenceFeatures,
    _charged_fraction,
    _clean_seq,
    _hydrophobic_fraction,
    _kd_mean,
    _soft_prior_score,
)

# ---------------------------------------------------------------------------
# Universal boundary scoring — species-independent structural signals
# ---------------------------------------------------------------------------
#
# These functions use coarsened amino acid properties (_AA_PROPERTY) to
# detect domain boundaries from biophysical constraints alone, without
# clade-specific motif tables.


def sp_boundary_excluded(seq: str, pos: int) -> bool:
    """Check if a candidate SP cleavage site has an impossible -3/-1 property combo.

    Returns True if the property classes at positions -3 and -1 form a
    combination that is absent from all 2406 ground-truth MHC SP cleavages.
    This is a species-independent biophysical exclusion filter layered on top
    of the classic von Heijne cleavage rule; see [S1] and [S2] in the module
    docstring.

    Args:
        seq: full protein sequence (uppercase)
        pos: candidate cleavage position (index of first mature residue)
    """
    if pos < 3 or pos >= len(seq):
        return False
    p3 = _AA_PROPERTY.get(seq[pos - 3], "")
    p1 = _AA_PROPERTY.get(seq[pos - 1], "")
    if not p3 or not p1:
        return False
    return (p3, p1) in _SP_EXCLUDED_M3_M1


def detect_h_region(seq: str, *, limit: int = 50) -> tuple[int, int]:
    """Find the signal peptide hydrophobic core (h-region) in the N-terminal region.

    Uses a windowed hydrophobicity scan (similar to TM detection) to find the
    best hydrophobic window in the first *limit* residues.  The h-region
    tolerates occasional polar interruptions (e.g., Ser/Thr in the middle of
    an otherwise hydrophobic stretch), which is common in real MHC signal
    peptides.

    Returns (h_start, h_end) as 0-indexed positions, or (0, 0) if no
    plausible h-region is found. This follows the h-region -> c-region view of
    SPase substrates rather than using only a local -3/-1 motif; see [S1] and
    [S4] in the module docstring.
    """
    region = seq[: min(limit, len(seq))].upper()
    if len(region) < 8:
        return 0, 0

    best_score = -999.0
    best_span = (0, 0)

    for win_len in range(8, 16):
        if len(region) < win_len:
            continue
        for start in range(0, len(region) - win_len + 1):
            window = region[start : start + win_len]
            hydro = sum(1 for aa in window if aa in _SP_HYDROPHOBIC)
            charged = sum(1 for aa in window if aa in _SP_CHARGED)
            frac = hydro / win_len
            charged_frac = charged / win_len
            # Score: hydrophobic fraction drives the signal, charged residues
            # penalize (normalized by window length so one charged residue
            # in an otherwise hydrophobic 8-aa window doesn't outweigh the
            # hydrophobic signal).  Prefer windows that start after position
            # ~2 (skip initial Met + n-region) but also prefer early windows
            # overall, since the SP h-region is always in the first ~25 aa.
            score = frac * 4.0 - charged_frac * 10.0
            if start < 2:
                score -= 0.5
            elif start > 20:
                score -= 0.5 * (start - 20) / 10.0
            if score > best_score:
                best_score = score
                best_span = (start, start + win_len)

    start, end = best_span
    if end <= start:
        return 0, 0

    window = region[start:end]
    hydro = sum(1 for aa in window if aa in _SP_HYDROPHOBIC)
    if hydro / len(window) < 0.50:
        return 0, 0

    # Trim trailing weakly-hydrophobic residues from the h-region end using
    # Kyte-Doolittle hydrophobicity.  The h-region core should end where the
    # hydrophobicity drops — GRAVY-like scoring is more nuanced than a binary
    # hydrophobic set.  Only trim if it leaves at least 7 residues.
    trimmed_end = end
    while trimmed_end > start + 7:
        kd = _KD_SCALE.get(region[trimmed_end - 1], 0.0)
        if kd >= 1.8:  # strongly hydrophobic (Ala=1.8 and above)
            break
        trimmed_end -= 1
    if trimmed_end > start + 6:
        end = trimmed_end

    return start, end


_SP_CORE_HYDROPHOBIC = _SP_HYDROPHOBIC | frozenset("C")


def has_plausible_sp_core(seq: str, sp_len: int) -> bool:
    """Return True if ``seq[:sp_len]`` contains a biophysically plausible SP h-region.

    A genuine eukaryotic signal peptide must contain a hydrophobic core *within*
    the cleaved prefix ``seq[:sp_len]``, not merely a hydrophobic beta-strand
    downstream in the mature domain.  Allows either a charged-free (``DEKR``-free)
    hydrophobic window or a denser core tolerating at most one charged residue
    (matching real MHC SPs like ``LLLLTRLLLL`` or ``LLLQLSEALALT``).  Leaderless
    prefixes lacking an N-terminal Methionine (``'M' not in sp[:8]``) require a
    denser hydrophobic core corresponding to an N-truncated SP tail rather than
    a hydrophilic mature N-terminus.
    """
    if sp_len < 6 or sp_len > len(seq):
        return False
    sp = seq[:sp_len].upper()
    hydro_total = sum(1 for aa in sp if aa in _SP_CORE_HYDROPHOBIC)
    if hydro_total / len(sp) < 0.36:
        return False

    has_nterm_met = "M" in sp[:8]
    max_win = min(16, len(sp))
    if has_nterm_met:
        for win_len in range(6, max_win + 1):
            for start in range(0, len(sp) - win_len + 1):
                window = sp[start : start + win_len]
                if window.count("P") > 1:
                    continue
                charged = sum(1 for aa in window if aa in _SP_CHARGED)
                if charged > 1:
                    continue
                hydro = sum(1 for aa in window if aa in _SP_CORE_HYDROPHOBIC)
                if charged == 0 and hydro >= 5 and hydro / win_len >= 0.58:
                    return True
                if charged == 1 and win_len >= 9 and hydro >= 7 and hydro / win_len >= 0.66:
                    return True
        return False

    for win_len in range(7, max_win + 1):
        for start in range(0, len(sp) - win_len + 1):
            window = sp[start : start + win_len]
            if window.count("P") > 1:
                continue
            charged = sum(1 for aa in window if aa in _SP_CHARGED)
            if charged > 1:
                continue
            hydro = sum(1 for aa in window if aa in _SP_CORE_HYDROPHOBIC)
            if charged == 0 and hydro >= 6 and hydro / win_len >= 0.72:
                return True
            if charged == 1 and win_len >= 10 and hydro >= 8 and hydro / win_len >= 0.80:
                return True
    return False


def estimate_sp_from_h_region(seq: str) -> int:
    """Estimate SP cleavage position from the hydrophobic core alone.

    Uses the biophysical rule that SP cleavage occurs ~5 residues after
    the h-region ends (in the polar c-region), at a small/neutral residue.
    Returns 0 if no plausible h-region is found.

    This provides a species-independent first-pass SP estimate that can
    anchor more detailed scoring. The intent is the same as the structural
    SP grammar described in [S1]-[S4] in the module docstring.
    """
    h_start, h_end = detect_h_region(seq)
    if h_end == 0:
        return 0

    seq_upper = seq.upper()
    # Search for a valid cleavage site in [h_end+1, h_end+12]
    # Prefer positions where -1 is A/G/S/C (von Heijne) and not excluded
    best_pos = 0
    best_score = -999.0
    for offset in range(1, 13):
        pos = h_end + offset
        if pos >= len(seq_upper) or pos < 3:
            continue
        if sp_boundary_excluded(seq_upper, pos):
            continue
        m1 = seq_upper[pos - 1]
        score = 0.0
        if m1 in "AG":
            score += 3.0
        elif m1 in "SC":
            score += 1.5
        elif m1 in "T":
            score += 0.5
        else:
            score -= 1.0
        # Prefer ~5 residues after h-end
        score -= 0.3 * abs(offset - 5)
        if score > best_score:
            best_score = score
            best_pos = pos

    return best_pos


def _score_h_region_alignment(
    seq: str,
    pos: int,
    *,
    h_region: tuple[int, int] | None = None,
) -> float:
    """Score a candidate SP cleavage site against the structural grammar:

        [n-region] → [h-region: 8-14 aa] → [c-region: 3-10 aa] → [cleavage]

    Signals scored (all species-independent):
      1. c-region length: soft prior peaked at 5-7, range [3, 10]
      2. c-region polarity: should be less hydrophobic than h-region
      3. Helix-breaking residue at c-region start (Pro/Gly = 26% in MHC SPs)
      4. Cleavage inside h-region: strong penalty
      5. Pro at -2: rare (0.6%) penalty

    Derived from 2373 ground truth MHC SPs: h-region median 8 aa,
    c-region median 6 aa (peak at 5, p5-p95 = [3, 10]). The scoring logic is
    the MHC-specific implementation of the general SP grammar summarized in
    [S1]-[S4] in the module docstring.
    """
    if pos <= 0 or pos >= len(seq):
        return 0.0

    if h_region is None:
        h_start, h_end = detect_h_region(seq)
    else:
        h_start, h_end = h_region
    if h_end <= h_start:
        return 0.0

    score = 0.0
    c_len = pos - h_end  # c-region length: h-region end to cleavage site

    # --- 1. c-region length prior ---
    # Peak at 5 (30.9%), 81% in [4, 9].  When multiple candidates remain
    # plausible, prefer the earlier coherent cut rather than drifting
    # farther downstream into a longer c-region.
    if c_len <= 0:
        # Cleavage inside or before h-region — biophysically implausible
        score -= 3.0
    elif 4 <= c_len <= 8:
        # Sweet spot: bonus peaked at 5-6, with a mild earlier preference.
        score += 1.6 - 0.2 * abs(c_len - 5)
    elif c_len == 9:
        score += 0.4
    elif 2 <= c_len <= 3:
        score += 0.3  # short but plausible (3.5% of cases)
    elif 10 <= c_len <= 12:
        score -= 0.4 + 0.25 * (c_len - 9)  # long c-regions drift late
    else:
        score -= 1.5  # very long or negative — suspicious

    # --- 2. c-region polarity ---
    # Mean hydrophobic fraction 0.38 in ground truth c-regions.
    # The c-region should be LESS hydrophobic than the h-region.
    if c_len > 0:
        c_region = seq[h_end:pos]
        c_hydro = sum(1 for aa in c_region if aa in _SP_HYDROPHOBIC) / len(c_region)
        if c_hydro <= 0.30:
            score += 0.5  # nicely polar
        elif c_hydro >= 0.65:
            score -= 1.0  # too hydrophobic for a c-region

    # --- 3. Helix-breaking at c-region start ---
    # Pro+Gly at the first c-region residue = 25.9% of MHC SPs.
    # Hydroxyl (S/T) = 29.9%.  Together these cover 55%.
    if h_end < len(seq):
        c_first = seq[h_end]
        if c_first in "PG":
            score += 0.6  # classic helix breaker
        elif c_first in "ST":
            score += 0.3  # polar, common at c-region start

    # --- 4. Pro at -2 ---
    # Only 0.6% of MHC SPs have Pro at -2 (barely above zero).
    if pos >= 2 and seq[pos - 2] == "P":
        score -= 0.8

    return score


_SP_SEQUENCE_CUE_MODEL_CACHE: Optional[dict[str, object]] = None


def load_sp_sequence_cue_model() -> dict[str, object]:
    """Load the learned MHC-specific sequence-cue model if present."""
    global _SP_SEQUENCE_CUE_MODEL_CACHE
    if _SP_SEQUENCE_CUE_MODEL_CACHE is not None:
        return _SP_SEQUENCE_CUE_MODEL_CACHE
    try:
        with _SP_SEQUENCE_CUE_MODEL_PATH.open("r", encoding="utf-8") as f:
            _SP_SEQUENCE_CUE_MODEL_CACHE = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        _SP_SEQUENCE_CUE_MODEL_CACHE = {}
    return _SP_SEQUENCE_CUE_MODEL_CACHE


def _cue_config_value(model: dict[str, object], key: str, default: float) -> float:
    """Return one numeric cue-model config value."""
    config = model.get("config", {})
    try:
        return float(config.get(key, default))
    except (TypeError, ValueError):
        return default


def _score_sp_nterm_state(
    seq: str,
    model: Optional[dict[str, object]] = None,
) -> float:
    """Score the raw N-terminus as SP-like vs mature-like from learned stats."""
    if not seq:
        return 0.0
    model = model or load_sp_sequence_cue_model()
    if not model:
        return 0.0
    nterm = model.get("nterm", {})
    score = 0.0
    start_m = nterm.get("start_m_log_odds", {})
    key = "present" if seq[:1] == "M" else "absent"
    try:
        score += float(start_m.get(key, 0.0))
    except (TypeError, ValueError):
        pass
    prefix3_scores = nterm.get("prefix3_log_odds", {})
    prefix3 = seq[:3] if len(seq) >= 3 else seq
    try:
        score += float(prefix3_scores.get(prefix3, 0.0))
    except (TypeError, ValueError):
        pass
    return score


def _exact_sp_shortcut(
    seq: str,
    model: dict[str, object],
) -> tuple[int, float, str, str]:
    """Return an exact-match early shortcut from learned repeated prefixes."""
    exact = model.get("exact_shortcuts", {})
    if not exact or not seq:
        return 0, 0.0, "", ""

    try:
        mature10_len = int(exact.get("mature10_len", 10))
    except (TypeError, ValueError):
        mature10_len = 10
    if len(seq) >= mature10_len:
        prefix10 = seq[:mature10_len]
        mature10 = exact.get("mature10_prefixes", {})
        if prefix10 in mature10:
            return 0, 0.995, "leaderless", "exact_mature10"

    if seq[:1] != "M":
        return 0, 0.0, "", ""

    try:
        prefix30_len = int(exact.get("sp_prefix30_len", 30))
    except (TypeError, ValueError):
        prefix30_len = 30
    if len(seq) < prefix30_len:
        return 0, 0.0, "", ""

    prefix30 = seq[:prefix30_len]
    payload = exact.get("sp_prefix30", {}).get(prefix30)
    if not isinstance(payload, dict):
        return 0, 0.0, "", ""
    try:
        split = int(payload.get("split", 0))
    except (TypeError, ValueError):
        split = 0
    if 3 <= split < len(seq):
        return split, 0.995, "sp_present", "exact_sp_prefix30"
    return 0, 0.0, "", ""


def _score_sp_boundary_words(
    seq: str,
    mature_start: int,
    model: Optional[dict[str, object]] = None,
) -> float:
    """Score a candidate cleavage boundary from exact and backed-off words."""
    if mature_start < 3 or mature_start + 3 > len(seq):
        return 0.0
    model = model or load_sp_sequence_cue_model()
    if not model:
        return 0.0
    boundary = model.get("boundary", {})
    weights = model.get("weights", {})
    left3 = seq[mature_start - 3 : mature_start]
    right3 = seq[mature_start : mature_start + 3]
    exact6 = left3 + right3
    score = 0.0
    try:
        score += float(boundary.get("boundary6_log_odds", {}).get(exact6, 0.0)) * float(weights.get("boundary6", 1.0))
    except (TypeError, ValueError):
        pass
    try:
        score += float(boundary.get("sp_end3_log_odds", {}).get(left3, 0.0)) * float(weights.get("sp_end3", 0.4))
    except (TypeError, ValueError):
        pass
    try:
        score += float(boundary.get("mature3_log_odds", {}).get(right3, 0.0)) * float(weights.get("mature3", 0.4))
    except (TypeError, ValueError):
        pass
    return score


def _rank_sp_word_candidates(
    seq: str,
    *,
    model: Optional[dict[str, object]] = None,
    max_candidates: int = 5,
) -> tuple[tuple[int, float], ...]:
    """Return top lexical boundary candidates in the plausible SP window."""
    model = model or load_sp_sequence_cue_model()
    if not model or len(seq) < 12:
        return ()
    config = model.get("config", {})
    try:
        search_min = int(config.get("shortcut_search_min", 8))
        search_max = int(config.get("shortcut_search_max", 45))
    except (TypeError, ValueError):
        search_min, search_max = 8, 45
    candidates: list[tuple[int, float]] = []
    for pos in range(max(3, search_min), min(search_max, len(seq) - 3) + 1):
        score = _score_sp_boundary_words(seq, pos, model)
        if score > 0.0:
            candidates.append((pos, score))
    candidates.sort(key=lambda item: item[1], reverse=True)
    return tuple(candidates[:max_candidates])


def infer_signal_peptide_candidates(
    seq: str,
    *,
    h_region: tuple[int, int] | None = None,
    max_candidates: int = SP_ESTIMATE_CANDIDATE_KEEP,
) -> tuple[tuple[int, float], ...]:
    """Return a ranked sequence-only SP cleavage zone.

    The current parser only needs a short, high-quality beam of plausible
    cleavage positions.  This keeps SP evidence flexible on divergent
    sequences without turning every parse into a wide search.
    """
    seq_upper = seq.upper() if seq else ""
    if len(seq_upper) < 15:
        return ()

    h_start, h_end = h_region if h_region is not None else detect_h_region(seq_upper)
    if h_end <= h_start or h_start > 25:
        return ()

    scored: list[tuple[int, float]] = []
    for pos in range(max(3, h_end + 3), min(len(seq_upper), h_end + 14)):
        if sp_boundary_excluded(seq_upper, pos):
            continue

        score = 0.0
        c_len = pos - h_end
        if 4 <= c_len <= 8:
            score += 2.0 - 0.3 * abs(c_len - 5)
        elif c_len == 9:
            score += 0.4
        elif c_len == 3:
            score += 0.8
        elif c_len <= 2:
            score -= 3.0
        else:
            score -= 0.8 * (c_len - 8)

        m1 = seq_upper[pos - 1]
        if m1 in _SP_CLEAVAGE_STRONG:
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

        if pos >= 3 and seq_upper[pos - 3] in _SP_SMALL_ALIPHATIC:
            score += 1.0
        elif pos >= 3 and seq_upper[pos - 3] in _SP_CHARGED:
            score -= 1.0

        if pos < len(seq_upper) and seq_upper[pos] == "P":
            score -= 2.0
        if pos >= 2 and seq_upper[pos - 2] == "P":
            score -= 0.8

        if h_end < len(seq_upper):
            c_first = seq_upper[h_end]
            if c_first in "PG":
                score += 0.4
            elif c_first in "ST":
                score += 0.2

        if c_len > 0:
            c_region = seq_upper[h_end:pos]
            c_hydro = sum(1 for aa in c_region if aa in _SP_HYDROPHOBIC)
            if c_hydro / len(c_region) <= 0.30:
                score += 0.3
            elif c_hydro / len(c_region) >= 0.65:
                score -= 0.8

        if pos >= 3 and pos + 3 <= len(seq_upper):
            pre_kd = sum(_KD_SCALE.get(seq_upper[i], 0.0) for i in range(pos - 3, pos)) / 3
            post_kd = sum(_KD_SCALE.get(seq_upper[i], 0.0) for i in range(pos, min(pos + 3, len(seq_upper)))) / 3
            kd_drop = pre_kd - post_kd
            if kd_drop > 0.5:
                score += min(kd_drop * 0.5, 1.5)

        if score > 0.0:
            scored.append((pos, score))

    scored.sort(key=lambda item: item[1], reverse=True)
    return tuple(scored[:max_candidates])


def fast_sp_triage(
    seq: str,
    *,
    h_region: tuple[int, int] | None = None,
) -> tuple[int, float, str, str, tuple[tuple[int, float], ...]]:
    """Return an optional early lexical SP shortcut, separate from general parsing.

    This path intentionally uses a stricter bar than the general parser:
    it only emits a shortcut when the N-terminus looks decisively SP-like
    or decisively mature-only. Otherwise callers should ignore the shortcut
    and let the full parser reason over the structural evidence.
    """
    model = load_sp_sequence_cue_model()
    if not model or not seq:
        return 0, 0.0, "", "", ()

    exact_est, exact_conf, exact_state, exact_kind = _exact_sp_shortcut(seq, model)
    if exact_state:
        return exact_est, exact_conf, exact_state, exact_kind, ()

    nterm_score = _score_sp_nterm_state(seq, model)
    h_start, h_end = h_region if h_region is not None else detect_h_region(seq)
    h_frac = 0.0
    if h_end > h_start:
        h_frac = _hydrophobic_fraction(seq, h_start, h_end)

    leaderless_threshold = _cue_config_value(model, "leaderless_shortcut_threshold", -3.0)
    leaderless_max_hfrac = _cue_config_value(model, "leaderless_shortcut_max_hfrac", 0.45)
    if nterm_score <= leaderless_threshold and h_frac <= leaderless_max_hfrac:
        confidence = min(1.0, max(0.75, 0.55 + 0.08 * abs(nterm_score)))
        return 0, confidence, "leaderless", "lexical_leaderless", ()

    if seq[:1] != "M":
        return 0, 0.0, "", "", ()

    strong_nterm = _cue_config_value(model, "sp_shortcut_nterm_threshold", 1.0)
    min_hfrac = _cue_config_value(model, "sp_shortcut_hfrac_min", 0.65)
    min_score = _cue_config_value(model, "sp_shortcut_score_threshold", 5.0)
    min_margin = _cue_config_value(model, "sp_shortcut_margin_threshold", 1.0)
    max_c_region = int(_cue_config_value(model, "sp_shortcut_max_c_region", 12.0))
    min_cleavage_score = _cue_config_value(model, "sp_shortcut_min_cleavage_score", 2.5)

    if nterm_score < strong_nterm or h_frac < min_hfrac:
        return 0, 0.0, "", "", ()

    candidates = _rank_sp_word_candidates(seq, model=model)
    if not candidates:
        return 0, 0.0, "", "", ()

    best_pos, best_score = candidates[0]
    second_score = candidates[1][1] if len(candidates) > 1 else -999.0
    if best_score < min_score or (best_score - second_score) < min_margin:
        return 0, 0.0, "", "", candidates
    if h_end > h_start:
        c_len = best_pos - h_end
        if c_len < 1 or c_len > max_c_region:
            return 0, 0.0, "", "", candidates
    cleavage_score = _score_sp_cleavage(seq, best_pos, h_region=(h_start, h_end))
    if cleavage_score < min_cleavage_score:
        return 0, 0.0, "", "", candidates
    confidence = min(1.0, 0.55 + 0.05 * max(0.0, best_score - second_score) + 0.10 * max(0.0, h_frac - min_hfrac))
    return best_pos, confidence, "sp_present", "lexical_boundary", candidates


def infer_signal_peptide(
    seq: str,
    *,
    h_region: tuple[int, int] | None = None,
) -> tuple[int, float]:
    """Infer signal peptide cleavage position from sequence alone.

    This is the unified, class-independent SP detector.  It uses only
    biophysical signals — no Cys pair positions, no MHC class knowledge,
    no species-specific motifs.

    The inference follows the structural grammar:
        [n-region] → [h-region: hydrophobic core] → [c-region: polar] → [cleavage]

    Returns (cleavage_pos, confidence) where cleavage_pos is the index of
    the first mature residue, or (0, 0.0) if no plausible SP is found.

    Called once per sequence by all three class-specific parsers.  The result
    constrains the mature_start search range: candidate parses whose
    mature_start is far from this estimate get penalized.
    """
    candidates = infer_signal_peptide_candidates(seq, h_region=h_region, max_candidates=1)
    if not candidates:
        return 0, 0.0
    best_pos, best_score = candidates[0]

    # Confidence from score magnitude and h-region quality
    confidence = min(1.0, max(0.0, best_score / 8.0))
    return best_pos, confidence


def _selected_sp_estimate(
    features: SequenceFeatures,
    *,
    use_early_shortcuts: bool,
) -> int:
    """Select which SP estimate to feed into structural parse enumeration.

    Exact-match shortcuts are allowed to override the grammar-based estimate
    when enabled because they are effectively dictionary lookups on highly
    repeated leaders or mature starts.  Softer lexical shortcuts still remain
    advisory only.
    """
    if use_early_shortcuts:
        if features.sp_shortcut_kind == "exact_sp_prefix30" and features.sp_shortcut_state == "sp_present" and features.sp_shortcut_estimate > 0:
            return features.sp_shortcut_estimate
        if features.sp_shortcut_kind == "exact_mature10" and features.sp_shortcut_state == "leaderless":
            return 0
    return features.sp_estimate


def _analyze_sequence_uncached(cleaned: str) -> SequenceFeatures:
    """Compute all structural features for one cleaned sequence."""
    pairs = find_cys_pairs(cleaned)
    annotations = tuple(classify_cys_pair(cleaned, c1, c2) for c1, c2, _ in pairs)

    topo: dict[tuple[int, int, str], float] = {}
    between: dict[tuple[int, int], float] = {}
    for c1, c2, _ in pairs:
        for expected in ("g_domain", "c_like"):
            topo[(c1, c2, expected)] = _score_pair_topology_support(cleaned, c1, c2, expected)
        between[(c1, c2)] = _score_between_cys(cleaned, c1, c2)

    h_start, h_end = detect_h_region(cleaned)
    sp_candidates = infer_signal_peptide_candidates(cleaned, h_region=(h_start, h_end))
    sp_est, sp_conf = infer_signal_peptide(cleaned, h_region=(h_start, h_end))
    shortcut_est, shortcut_conf, shortcut_state, shortcut_kind, shortcut_candidates = fast_sp_triage(
        cleaned,
        h_region=(h_start, h_end),
    )
    nterm_state_score = _score_sp_nterm_state(cleaned)

    tail_start = max(0, len(cleaned) - 80)
    tm_start, tm_end = _find_tm_span(cleaned[tail_start:], raw_start=tail_start)

    hydrophobic_prefix = [0]
    charged_prefix = [0]
    kd_prefix = [0.0]
    for aa in cleaned:
        hydrophobic_prefix.append(hydrophobic_prefix[-1] + int(aa in _SP_HYDROPHOBIC))
        charged_prefix.append(charged_prefix[-1] + int(aa in _SP_CHARGED))
        kd_prefix.append(kd_prefix[-1] + _KD_SCALE.get(aa, 0.0))

    return SequenceFeatures(
        seq=cleaned,
        seq_len=len(cleaned),
        cys_pairs=tuple(pairs),
        pair_annotations=annotations,
        topology_scores=topo,
        h_region=(h_start, h_end),
        sp_estimate=sp_est,
        sp_estimate_confidence=sp_conf,
        sp_estimate_candidates=sp_candidates,
        sp_shortcut_estimate=shortcut_est,
        sp_shortcut_confidence=shortcut_conf,
        sp_shortcut_state=shortcut_state,
        sp_shortcut_kind=shortcut_kind,
        sp_nterm_state_score=nterm_state_score,
        sp_word_candidates=shortcut_candidates,
        hydrophobic_prefix=tuple(hydrophobic_prefix),
        charged_prefix=tuple(charged_prefix),
        kd_prefix=tuple(kd_prefix),
        tm_span=(tm_start, tm_end),
        between_cys_scores=between,
    )


@lru_cache(maxsize=4096)
def _analyze_sequence_cached(cleaned: str) -> SequenceFeatures:
    """Cache exact-sequence feature extraction for repeated proteins."""
    return _analyze_sequence_uncached(cleaned)


def analyze_sequence(seq: str) -> SequenceFeatures:
    """Compute all structural features for *seq* in a single pass.

    The returned ``SequenceFeatures`` should be passed to
    ``decompose_class_i``, ``decompose_class_ii_alpha``,
    ``decompose_class_ii_beta``, and ``refine_signal_peptide`` via
    their ``features=`` keyword argument to avoid redundant work.
    Exact cleaned-sequence repeats reuse a cached feature bundle.
    """
    cleaned = _clean_seq(seq)
    return _analyze_sequence_cached(cleaned)


def _score_sp_candidate(
    seq: str,
    pos: int,
    junction_pos: int,
    *,
    h_region: tuple[int, int] = (0, 0),
) -> float:
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
        Known or estimated α1/α2 junction position (from find_junction).
        Used to verify that the implied α1 domain length is plausible.
    """
    if pos < 3 or pos >= len(seq) - 3:
        return -999.0

    # Fast exclusion: impossible -3/-1 property combos
    if sp_boundary_excluded(seq, pos):
        return -10.0

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

    # Signal 9: agreement with the signal-peptide hydrophobic core
    score += _score_h_region_alignment(seq, pos, h_region=h_region) * _SP_H_REGION_WEIGHT

    # Signal 10: α1 domain length plausibility (from junction position)
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

    Uses a scored search that combines SP cleavage signals with α1/α2
    junction confirmation and groove-length plausibility.

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
        junction_pos, junction_score = find_junction(seq, groove_c1)
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

    _h_region = detect_h_region(seq)

    best_pos = traditional
    best_score = _score_sp_candidate(seq, traditional, junction_pos, h_region=_h_region)

    for pos in range(search_min, search_max + 1):
        s = _score_sp_candidate(seq, pos, junction_pos, h_region=_h_region)
        if s > best_score:
            best_score = s
            best_pos = pos

    return best_pos


# ---------------------------------------------------------------------------
# Candidate-parse enumeration
# ---------------------------------------------------------------------------
#
# Enumerate candidate mature_start positions for each plausible Cys pair and
# score the full [SP]-[groove]-[Ig]-[tail] decomposition. Pick the highest-
# scoring whole parse rather than anchoring on one absolute position.
#
# Domain size expectations (aa):
#   Class I:   groove1 (α1) 80-100, groove2 (α2) ~93, Ig (α3) ~76
#   Class II α: groove1 (α1) 75-95, Ig (α2) ~76
#   Class II β: groove2 (β1) 85-100, Ig (β2) ~76
#   SP: 14-45 (0 = stripped)


def _score_boundary_family(seq: str, pos: int) -> float:
    """Score a candidate SP cleavage site by residue-class boundary family.

    Uses log-odds ratios of property classes at positions -3 to +3 relative
    to the cleavage site, comparing true GT boundaries vs nearby non-boundary
    positions.  This is a species-independent signal that captures non-mammalian
    cleavage patterns the von Heijne rules miss.
    """
    if pos < 3 or pos + 4 > len(seq):
        return 0.0
    score = 0.0
    for off in range(-3, 4):
        prop = _AA_PROPERTY.get(seq[pos + off], "")
        lo = SP_BOUNDARY_LOG_ODDS.get((off, prop), 0.0)
        score += lo
    return score * 0.15  # scale to contribute ~1-2 points for typical boundaries


def _score_sp_cleavage(
    seq: str,
    pos: int,
    *,
    h_region: tuple[int, int] | None = None,
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Score cleavage site quality at position *pos* (first mature residue).

    Universal biochemical signals: von Heijne -1/-3 rules, upstream
    hydrophobic density, c-region polarity, hydrophobicity transition,
    and conserved MHC mature-start properties (H at +3, charged +1).
    """
    if pos < 3 or pos >= len(seq) - 3:
        return -999.0

    # Fast exclusion: impossible -3/-1 property combos
    if sp_boundary_excluded(seq, pos):
        return -10.0

    score = 0.0
    n = len(seq)

    # -1 cleavage residue (von Heijne)
    # Among the strong residues (A/G), G is slightly preferred: SPase I
    # encounters the first valid site and cleaves there.  In 951 cases
    # where both GT and GT+1 have valid -1 residues, GT picks G at -1
    # 32% of the time but GT+1 has G at -1 50% of the time — the enzyme
    # prefers the earlier (G) cleavage.
    m1 = seq[pos - 1]
    if m1 == "G":
        score += 3.3  # slight G preference for "cleave earlier" tiebreaker
    elif m1 == "A":
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

    # -3 residue (von Heijne)
    if seq[pos - 3] in _SP_SMALL_ALIPHATIC:
        score += 1.0
    elif seq[pos - 3] in _SP_CHARGED:
        score -= 1.0

    # +1 not proline
    if seq[pos] == "P":
        score -= 2.0

    # Upstream hydrophobic density (h-region)
    if pos >= 12:
        hydro_frac = _hydrophobic_fraction(seq, pos - 12, pos - 3, features=features)
        score += (hydro_frac - 0.4) * 2.0

    # c-region polarity
    if pos >= 5:
        c_hydro = _hydrophobic_fraction(seq, pos - 3, pos, features=features)
        if c_hydro > 0.66:
            score -= 1.5

    # H at +3 — conserved structural His in MHC α1 S1 strand
    if pos + 2 < n and seq[pos + 2] == "H":
        score += 2.0

    # First mature residue properties
    if pos < n:
        p1 = seq[pos]
        if p1 in "AGCVIL":
            score += 0.6
        elif p1 in "ST":
            score += 0.2
        elif p1 == "D":
            score -= 0.8
        elif p1 == "E":
            score -= 0.2

    # KD hydrophobicity transition
    if pos >= 4 and pos + 3 <= n:
        pre_kd = _kd_mean(seq, pos - 3, pos, features=features)
        post_kd = _kd_mean(seq, pos, pos + 3, features=features)
        kd_drop = pre_kd - post_kd
        if kd_drop > 0.5:
            score += min(kd_drop * 0.5, 1.5)

    score += _score_h_region_alignment(seq, pos, h_region=h_region) * _SP_H_REGION_WEIGHT

    # Boundary-family scoring: residue-class log-odds at -3..+3
    score += _score_boundary_family(seq, pos)

    return score


def _score_nterm_lexical_state(
    seq: str,
    mature_start: int,
) -> float:
    """Return an MHC-specific lexical cue score from learned sequence stats.

    This is the general lexical path, not the early shortcut.  It contributes
    moderate evidence to all parses and remains active even when early lexical
    shortcuts are disabled for benchmarking.
    """
    if not seq:
        return 0.0
    model = load_sp_sequence_cue_model()
    if model:
        nterm_score = _score_sp_nterm_state(seq, model)
        if mature_start > 0:
            return nterm_score + _score_sp_boundary_words(seq, mature_start, model)
        return -nterm_score

    # Fallback if the learned cue model is unavailable.
    prefix3 = seq[:3] if len(seq) >= 3 else seq
    score = 0.0
    if mature_start > 0:
        if prefix3 in _MHC_SP_PRESENT_PREFIX_TRIADS:
            score += 2.0
        if prefix3 in _MHC_MATURE_PREFIX_TRIADS:
            score -= 3.0
        if mature_start + 3 <= len(seq):
            mature3 = seq[mature_start : mature_start + 3]
            if mature3 in _MHC_BOUNDARY_MATURE_TRIADS:
                score += 1.25
    else:
        if prefix3 in _MHC_MATURE_PREFIX_TRIADS:
            score += 2.5
        if prefix3 in _MHC_SP_PRESENT_PREFIX_TRIADS:
            score -= 3.0
    return score


def _score_sp_length(mature_start: int) -> float:
    """Score signal-peptide length with extra pressure against long leaders."""
    score = _soft_prior_score(mature_start, *_SP_PRIOR)
    if mature_start > 31:
        score -= (mature_start - 31) * 0.15
    if mature_start > 35:
        score -= (mature_start - 35) * 0.35
    return score


def _score_sp_estimate_consistency(
    pos: int,
    sp_estimate: int,
    confidence: float,
) -> float:
    """Score how compatible a boundary is with the sequence-only SP estimate.

    This is intentionally smooth and confidence-driven: strong sequence-only SP
    evidence should keep later structural scoring from drifting too far,
    especially downstream, while low-confidence estimates should stay weak.
    """
    if sp_estimate <= 0 or confidence <= 0.0:
        return 0.0
    delta = pos - sp_estimate
    span = max(4.0, 10.0 - 5.0 * min(confidence, 1.0))
    distance = abs(delta)
    if delta > 0:
        distance *= _SP_LATE_DRIFT_SCALE
    return max(-4.0, (1.0 + 2.5 * confidence) - (distance / span) * (2.0 + confidence))


def _score_sp_estimate_zone_consistency(
    pos: int,
    candidates: Sequence[tuple[int, float]],
    confidence: float,
) -> float:
    """Score compatibility with a short ranked zone of SP cleavage proposals.

    The sequence-only SP detector is most useful as a small beam of plausible
    cuts rather than a single hard point estimate.  This function rewards
    agreement with any strong candidate and turns into a contradiction only
    when *all* plausible cuts are far away.
    """
    if not candidates or confidence <= 0.0:
        return 0.0

    best_raw = max(float(score) for _pos, score in candidates)
    if best_raw <= 0.0:
        return 0.0

    span = max(3.0, 9.0 - 4.0 * min(confidence, 1.0))
    best_score = -999.0
    for cand_pos, cand_score in candidates:
        rel = max(0.25, min(1.0, float(cand_score) / best_raw))
        delta = pos - cand_pos
        distance = abs(delta)
        if delta > 0:
            distance *= _SP_LATE_DRIFT_SCALE
        candidate_score = rel * ((1.0 + 2.5 * confidence) - (distance / span) * (2.0 + confidence))
        if candidate_score > best_score:
            best_score = candidate_score
    return max(-4.0, best_score)


def _sp_zone_is_ambiguous(features: Optional[SequenceFeatures]) -> bool:
    """Whether the sequence-only SP zone has multiple close competitors."""
    if features is None:
        return False
    candidates = features.sp_estimate_candidates
    if len(candidates) < 2:
        return False
    (best_pos, best_score), (next_pos, next_score) = candidates[:2]
    if features.sp_estimate_confidence < PRIMARY_PARSE_LOW_CONFIDENCE:
        return False
    score_ratio = (float(next_score) / float(best_score)) if best_score > 0 else 0.0
    if abs(best_pos - next_pos) <= 4 and score_ratio >= 0.75 and (best_score - next_score) <= 1.5:
        return True
    return False


# ---------------------------------------------------------------------------
# N-terminal mature protein motifs (first ~5 residues after cleavage)
# ---------------------------------------------------------------------------
# Class I α1 N-terminus:
#   Mammals: G/S at +1, P/S at +2, H at +3 (GSH, GPH, SPH patterns)
#   Birds:   E at +1, L/P at +2, H at +3 (ELH, EPH patterns)
#   Fish:    K/D/E at +1, various at +2, H at +3
#   Universal: H at +3 is 33-52% across all vertebrates.
#
# Class II α1 N-terminus:
#   Mammals: variable first residue, often K/E/D/I
#   Less conserved than class I.
#
# Class II β1 N-terminus:
#   Mammals: often R/D/E at +1, F at +2
#   Less conserved than class I.


def refinement_group(species_category: str) -> str:
    if species_category in _MAMMAL_CATEGORIES or not species_category:
        return "mammal"
    if species_category == "bird":
        return "bird"
    if species_category == "fish":
        return "fish"
    return "other_vertebrate"


_SP_BOUNDARY_MODEL_CACHE: Optional[dict[str, object]] = None


def load_sp_boundary_model() -> dict[str, object]:
    """Load the learned SP motif-composition model if present."""
    global _SP_BOUNDARY_MODEL_CACHE
    if _SP_BOUNDARY_MODEL_CACHE is not None:
        return _SP_BOUNDARY_MODEL_CACHE
    try:
        with _SP_BOUNDARY_MODEL_PATH.open("r", encoding="utf-8") as f:
            _SP_BOUNDARY_MODEL_CACHE = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        _SP_BOUNDARY_MODEL_CACHE = {}
    return _SP_BOUNDARY_MODEL_CACHE


def _score_sp_boundary_composition(
    seq: str,
    mature_start: int,
    species_category: str = "",
) -> float:
    """Score a cleavage site using learned motif-composition log-odds.

    The model is trained from positives vs nearby decoys and uses both
    exact residues and coarse amino-acid classes around the boundary.
    """
    model = load_sp_boundary_model()
    if not model:
        return 0.0
    groups = model.get("groups", {})
    group = refinement_group(species_category)
    payload = groups.get(group) or groups.get("global")
    if not payload:
        return 0.0

    residue_scores = payload.get("residue_log_odds", {})
    class_scores = payload.get("class_log_odds", {})
    score = 0.0
    for offset_text, residue_table in residue_scores.items():
        try:
            offset = int(offset_text)
        except ValueError:
            continue
        pos = mature_start + offset
        if not (0 <= pos < len(seq)):
            continue
        aa = seq[pos]
        score += float(residue_table.get(aa, 0.0))
        aa_class = _SP_BOUNDARY_AA_CLASS.get(aa)
        if aa_class:
            score += float(class_scores.get(offset_text, {}).get(aa_class, 0.0)) * _SP_BOUNDARY_CLASS_WEIGHT
    return score * _SP_BOUNDARY_MODEL_WEIGHT_BY_GROUP.get(group, 0.0)


def _score_mature_nterm(seq: str, mature_start: int, mhc_class: str) -> float:
    """Score the first ~5 residues of the mature protein."""
    if mature_start <= 0:
        return 0.0
    n = len(seq)
    tables = _MATURE_START_CLASS_I if mhc_class == "I" else _MATURE_START_CLASS_II
    score = 0.0
    for offset, residue_scores in tables.items():
        pos = mature_start + offset
        if pos < n:
            score += residue_scores.get(seq[pos], 0.0)
    return score


def _score_refinement_motif(
    seq: str,
    mature_start: int,
    species_category: str,
    mhc_class: str = "",
) -> float:
    """Score mature-start motif composition around a candidate boundary."""
    if mature_start < 0 or mature_start + 3 > len(seq):
        return 0.0
    group = refinement_group(species_category)
    motif = seq[mature_start : mature_start + 3]
    score = _score_sp_boundary_composition(seq, mature_start, species_category)
    score += _REFINEMENT_MOTIFS.get(group, {}).get(motif, 0.0)
    if mhc_class in {"I", "II"}:
        score += _score_mature_nterm(seq, mature_start, mhc_class) * 0.5
    return score


# Signal peptide cleavage site residues.
#
# The -1 position (last residue of the SP) is A/G/S/C across all jawed
# vertebrates: sharks 96%, mammals 98%, birds 100%, fish 88%, reptiles 84%,
# amphibians 95% (including C at 25%).
#
def _score_refinement_candidate(
    seq: str,
    pos: int,
    initial_start: int,
    species_category: str = "",
    mhc_class: str = "",
    *,
    h_region: tuple[int, int] = (0, 0),
    features: Optional[SequenceFeatures] = None,
) -> float:
    """Score a candidate SP cleavage site around an anchor prediction."""
    if pos < 3 or pos >= len(seq) - 3:
        return -999.0

    # Fast exclusion: impossible -3/-1 property combos
    if sp_boundary_excluded(seq, pos):
        return -10.0

    score = 0.0
    group = refinement_group(species_category)
    sp_estimate = features.sp_estimate if features is not None else 0
    sp_confidence = features.sp_estimate_confidence if features is not None else 0.0
    sp_candidates = features.sp_estimate_candidates if features is not None else ()

    # Signal 1: -1 cleavage residue
    m1 = seq[pos - 1]
    if m1 in _SP_CLEAVAGE_STRONG:
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

    # Signal 2: -3 residue and canonical AxA cleavage motif.
    if seq[pos - 3] in _SP_SMALL_ALIPHATIC:
        score += 1.0
    elif seq[pos - 3] in _SP_CHARGED:
        score -= 1.0
    if m1 in "AGS" and seq[pos - 3] in _SP_SMALL_ALIPHATIC:
        score += 1.0

    # Signal 3: +1 not proline
    if seq[pos] == "P":
        score -= 2.0

    # Signal 4: upstream hydrophobic density (h-region)
    if pos >= 12:
        hydro_frac = _hydrophobic_fraction(seq, pos - 12, pos - 3, features=features)
        h_weight = 2.0 if group == "mammal" else 1.0
        score += (hydro_frac - 0.4) * h_weight
        if hydro_frac < 0.3:
            score -= 2.0

    # Signal 5: c-region polarity
    if pos >= 5:
        c_hydro = _hydrophobic_fraction(seq, pos - 3, pos, features=features)
        if c_hydro > 0.66:
            score -= 1.5

    # Signal 6: distance from the anchor prediction.
    delta_from_initial = pos - initial_start
    distance_from_initial = abs(delta_from_initial)
    if delta_from_initial > 0:
        distance_from_initial *= _SP_LATE_DRIFT_SCALE
    score -= distance_from_initial * _REFINEMENT_DISTANCE_PENALTY[group]
    if sp_candidates:
        score += _score_sp_estimate_zone_consistency(pos, sp_candidates, sp_confidence)
    else:
        score += _score_sp_estimate_consistency(pos, sp_estimate, sp_confidence)

    # Signal 7: mature protein start properties and category-specific motifs.
    if pos + 2 < len(seq):
        if seq[pos + 2] == "H":
            score += 2.0
        p1 = seq[pos]
        if p1 in "AGCVIL":
            score += 0.8
        elif p1 in "ST":
            score += 0.3
        elif p1 == "D":
            score -= 0.8
        elif p1 == "E":
            if refinement_group(species_category) == "bird":
                score += 0.3
            else:
                score -= 0.2
    score += _score_refinement_motif(seq, pos, species_category, mhc_class)

    # Signal 8: charged mature-start density.
    if pos + 5 <= len(seq):
        charged_frac = _charged_fraction(seq, pos, pos + 5, features=features)
        if charged_frac > 0.3:
            score += 1.0
        elif charged_frac < 0.1:
            score -= 1.0

    # Signal 9: hydrophobicity transition at cleavage boundary.
    if pos >= 4 and pos + 3 <= len(seq):
        pre_kd = _kd_mean(seq, pos - 3, pos, features=features)
        post_kd = _kd_mean(seq, pos, pos + 3, features=features)
        kd_drop = pre_kd - post_kd
        if kd_drop > 0.5:
            score += min(kd_drop * 0.5, 1.5)

    # Signal 10: alignment with the signal-peptide hydrophobic core.
    score += _score_h_region_alignment(seq, pos, h_region=h_region) * _SP_H_REGION_WEIGHT

    # Signal 10b: MHC-specific statistical boundary words.
    score += 0.20 * _score_sp_boundary_words(seq, pos)

    # Signal 11: don't let the cleavage drift into the SP n-region.
    if seq[:1] == "M" and pos < 12:
        score -= 2.0

    return score


def refine_signal_peptide(
    sequence: str,
    mature_start: int,
    species_category: str = "",
    mhc_class: str = "",
    *,
    features: Optional["SequenceFeatures"] = None,
    groove_anchor: Optional[tuple[int, int]] = None,
) -> int:
    """Refine signal peptide cleavage site using sequence features.

    Parameters
    ----------
    sequence : str
        Full protein sequence (including signal peptide).
    mature_start : int
        Predicted mature protein start from Cys-pair heuristic.
    species_category : str
        One of the mhcseqs species categories.  Used for motif and
        distance-prior calibration.
    mhc_class : str
        Optional MHC class hint ("I" or "II") so class-specific mature
        N-terminal signals can participate in the score.
    features : SequenceFeatures, optional
        Pre-computed features from ``analyze_sequence()``.  Avoids
        redundant Cys-pair scanning and classification.
    groove_anchor : (c1, c2), optional
        Pre-identified groove anchor pair from the parser, used for
        junction detection without re-scanning.
    """
    if mature_start <= 0 or not sequence:
        return mature_start

    seq = features.seq if features else sequence.upper()
    if mature_start >= len(seq):
        return mature_start
    if (
        features is not None
        and features.sp_shortcut_kind == "exact_sp_prefix30"
        and features.sp_shortcut_state == "sp_present"
        and features.sp_shortcut_estimate > 0
        and abs(mature_start - features.sp_shortcut_estimate) <= 1
    ):
        return features.sp_shortcut_estimate

    # Leaderless gate: every real-SP entry in the benchmark has a detectable
    # h-region (≥6 aa).  When the h-region is absent and the sequence-level
    # SP predictor also says 0, the Cys-pair-derived mature_start is almost
    # certainly a false positive on a mature-only input — override to 0.
    if (
        features is not None
        and (features.h_region[1] - features.h_region[0]) < 6
        and features.sp_estimate == 0
        and features.sp_shortcut_state != "sp_present"
    ):
        return 0

    group = refinement_group(species_category)
    window = _REFINEMENT_WINDOW_BY_GROUP[group]

    # Try to find the groove Cys pair and junction position.
    junction_pos = 0
    if groove_anchor:
        junction_pos, _ = find_junction(seq, groove_anchor[0])
    else:
        pairs = list(features.cys_pairs) if features else find_cys_pairs(seq)
        if pairs:
            annotated = list(features.pair_annotations) if features else [classify_cys_pair(seq, c1, c2) for c1, c2, _sep in pairs]
            groove_pairs = [a for a in annotated if a.domain_type in ("groove", "ambiguous")]
            if groove_pairs:
                groove = min(groove_pairs, key=lambda a: abs(a.c1 - (mature_start + 100)))
                junction_pos, _ = find_junction(seq, groove.c1)

    # Compute h_region once for all refinement candidates.
    _h_region = features.h_region if features is not None else detect_h_region(seq)

    best_pos = mature_start
    initial_score = _score_refinement_candidate(
        seq,
        mature_start,
        mature_start,
        species_category,
        mhc_class,
        h_region=_h_region,
        features=features,
    )
    best_score = initial_score
    for pos in range(max(5, mature_start - window), min(len(seq) - 3, mature_start + window + 1)):
        s = _score_refinement_candidate(
            seq,
            pos,
            mature_start,
            species_category,
            mhc_class,
            h_region=_h_region,
            features=features,
        )
        if s > best_score + _SP_EARLIER_TIE_MARGIN:
            best_score = s
            best_pos = pos
        elif pos < best_pos and s >= best_score - _SP_EARLIER_TIE_MARGIN:
            best_score = s
            best_pos = pos

    junction_score = 0.0
    if junction_pos > 0:
        junction_score = score_junction(seq, junction_pos)
    if junction_score >= 6.0:
        for pos in range(max(5, mature_start - 15), min(len(seq) - 3, mature_start + 16)):
            if mature_start - window <= pos <= mature_start + window:
                continue
            s = _score_refinement_candidate(
                seq,
                pos,
                mature_start,
                species_category,
                mhc_class,
                h_region=_h_region,
                features=features,
            )
            alpha1_len = junction_pos - pos
            if 80 <= alpha1_len <= 95:
                s += 2.5
            elif 75 <= alpha1_len <= 125:
                s += 1.0
            else:
                continue
            if s > best_score + _SP_EARLIER_TIE_MARGIN:
                best_score = s
                best_pos = pos
            elif pos < best_pos and s >= best_score - _SP_EARLIER_TIE_MARGIN:
                best_score = s
                best_pos = pos

    if best_pos != mature_start:
        move = abs(best_pos - mature_start)
        if move > 5:
            required_margin = 0.30 + 0.06 * move
            if best_pos > mature_start:
                required_margin += 0.10
            if best_score < initial_score + required_margin:
                best_pos = mature_start

    # Biophysical SP core gate: the cleaved prefix seq[:best_pos] must itself
    # contain a charged-free hydrophobic core (not just a hydrophobic beta-strand
    # downstream of mature_start that fell inside detect_h_region(limit=50)).
    if best_pos > 0 and (features is None or features.sp_shortcut_state != "sp_present") and not has_plausible_sp_core(seq, best_pos):
        return 0

    return best_pos


def _append_long_sp_flag(flags: list[str], mature_start: int) -> None:
    """Record an SP-length warning when the mature start looks unusually deep."""
    if mature_start > 40:
        flags.append(f"long_sp({mature_start})")


_OVERCAPTURE_CHARGED = frozenset("RKDE")


# Minimum SP length to even consider; over-capture inflates SP length, real SPs
# in the benchmark rarely exceed this.
_OVERCAPTURE_MIN_SP = 25


# A real SP n-region is short (positive-inside rule, <= ~8 aa); a discarded prefix
# longer than this cannot be a genuine n-region, so it is the over-capture cutoff.
_OVERCAPTURE_MIN_PREFIX = 10


# The internal-Met SP must itself be long enough to be a real signal peptide.
_OVERCAPTURE_MIN_SUB_SP = 12


# The internal-Met SP must score at least this confidently to count as real.
_OVERCAPTURE_MIN_SUB_CONF = 0.5


# h-region scan within the prefix; >= this length means the prefix is itself a
# plausible leader (so not over-capture).
_OVERCAPTURE_MAX_PREFIX_HREGION = 6


# Spurious UTR-derived prefixes are charged/polar; require this charged fraction.
_OVERCAPTURE_MIN_PREFIX_CHARGED = 0.2


# Window of residues after a candidate internal Met passed to the SP detector.
# The SP + cleavage site fall well within this; verified result-identical to
# scanning the full tail across the dataset, but avoids copying the whole protein.
_OVERCAPTURE_SUB_WINDOW = 60


def detect_sp_internal_met_overcapture(seq: str, mature_start: int) -> Optional[int]:
    """Detect signal-peptide over-capture from an upstream weak-Kozak start codon.

    Some records — overwhelmingly unreviewed UniProt/TrEMBL entries and a few
    submitter-deposited GenBank/IPD sequences — are translated from an in-frame
    ATG *upstream* of the true start codon.  Ribosomes leaky-scan past such a
    weak-Kozak AUG to a stronger downstream one (the -3 purine / +4 G context
    dominates initiation), but a naive longest-ORF translation grabs the upstream
    AUG and prepends a stretch of translated 5'-UTR to the real signal peptide.

    At the amino-acid level the artifact looks like an unusually long (>25 aa)
    "signal peptide" whose N-terminal portion is a charged/polar segment with no
    hydrophobic h-region, followed by an internal Met that begins a clean,
    self-contained signal peptide cleaving at the *same* mature start.

    A genuine SP n-region is also short and positively charged (the
    "positive-inside" rule), so composition alone cannot separate the two — e.g.
    the real HLA class I leader ``MRV|MAPRTL...`` and HLA-DPA1 ``MRPEDR|MFHIRA...``
    both have a charged head before the h-region.  The discriminator is *length*:
    a real n-region is <= ~8 aa, so this requires the discarded prefix to be at
    least 10 aa — longer than any plausible n-region.

    This is a *QC signal, not a correction*.  TrEMBL has no manually curated SP
    ground truth (where a Signal feature exists at all it is SignalP-predicted on
    the same over-captured N-terminus, so it shares the error), and the only
    independent confirmation is the source CDS / Kozak context, which is not
    available for every entry.  Callers must not rewrite the sequence on this
    signal alone.

    Parameters
    ----------
    seq : str
        Full protein sequence (including the putative signal peptide).
    mature_start : int
        Predicted mature start (== SP length) to evaluate for over-capture.

    Returns
    -------
    Optional[int]
        The internal Met position (the inferred true start codon) when the
        over-capture signature is present, else ``None``.
    """
    if mature_start <= _OVERCAPTURE_MIN_SP or mature_start >= len(seq):
        return None
    _, full_conf = infer_signal_peptide(seq)
    best: Optional[tuple[int, float]] = None
    for m in range(_OVERCAPTURE_MIN_PREFIX, mature_start):  # prefix > plausible n-region
        if seq[m] != "M":
            continue
        if mature_start - m < _OVERCAPTURE_MIN_SUB_SP:  # remaining SP too short to be real
            continue
        sub_pos, sub_conf = infer_signal_peptide(seq[m : m + _OVERCAPTURE_SUB_WINDOW])
        if sub_pos <= 0:
            continue
        if abs((m + sub_pos) - mature_start) > 2:  # must cleave at the same mature start
            continue
        if sub_conf < _OVERCAPTURE_MIN_SUB_CONF or sub_conf < full_conf:
            continue  # internal SP must be a real SP and at least as good as the full one
        if best is None or sub_conf > best[1]:
            best = (m, sub_conf)
    if best is None:
        return None
    m = best[0]
    prefix = seq[:m]
    h_start, h_end = detect_h_region(prefix)
    if h_end - h_start >= _OVERCAPTURE_MAX_PREFIX_HREGION:
        return None  # prefix itself looks like a leader -> not over-capture
    charged = sum(c in _OVERCAPTURE_CHARGED for c in prefix) / len(prefix)
    if charged < _OVERCAPTURE_MIN_PREFIX_CHARGED:
        return None  # spurious UTR-derived prefixes are charged/polar
    return m


def _append_sp_overcapture_flag(flags: list[str], seq: str, mature_start: int) -> None:
    """Record a QC flag when the SP looks over-captured from an upstream AUG."""
    rescued = detect_sp_internal_met_overcapture(seq, mature_start)
    if rescued is not None:
        flags.append(f"sp_internal_met_overcapture({rescued})")
