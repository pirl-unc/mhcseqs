#!/usr/bin/env python3
"""Train a lightweight statistical sequence-cue model for MHC signal peptides.

The model captures two kinds of MHC-specific lexical evidence:

1. N-terminal state cues:
   - whether the raw sequence starts with Met
   - exact first-3aa prefixes enriched in full-length SP-positive proteins
     versus mature-only / SP-stripped controls

2. Boundary word cues:
   - exact 6aa boundary words (last 3 aa of SP + first 3 aa of mature)
   - backed off 3aa words on each side of the boundary

The output is written to ``data/sp_sequence_cue_model.json`` and consumed by
``mhcseqs.domain_parsing``.  The parser uses these scores in two layers:

- general parsing: modest lexical evidence inside the holistic scorer
- early shortcuts: optional, high-confidence direct SP cues that can be
  disabled independently for benchmarking
"""

from __future__ import annotations

import csv
import json
import math
import statistics
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.evaluate_sp_ground_truth import GT_CSV, NEGATIVE_CONTROL_CSV

OUT_JSON = ROOT / "data" / "sp_sequence_cue_model.json"
ALPHA = 1.0
CLIP = 5.0
SEARCH_MIN = 8
SEARCH_MAX = 45
EXACT_PREFIX30_LEN = 30
EXACT_MATURE10_LEN = 10
EXACT_SP_PREFIX30_MIN_COUNT = 5
EXACT_MATURE10_MIN_COUNT = 10


def _clip(value: float, lo: float = -CLIP, hi: float = CLIP) -> float:
    return max(lo, min(hi, value))


def _log_odds(
    pos_count: int,
    pos_total: int,
    neg_count: int,
    neg_total: int,
    *,
    alpha: float = ALPHA,
) -> float:
    pos_freq = (pos_count + alpha) / max(pos_total + 2.0 * alpha, 1.0)
    neg_freq = (neg_count + alpha) / max(neg_total + 2.0 * alpha, 1.0)
    return _clip(math.log(pos_freq / neg_freq))


def _table_log_odds(
    positive: Counter[str],
    negative: Counter[str],
) -> dict[str, float]:
    pos_total = sum(positive.values())
    neg_total = sum(negative.values())
    tokens = set(positive) | set(negative)
    out: dict[str, float] = {}
    for token in sorted(tokens):
        score = _log_odds(positive[token], pos_total, negative[token], neg_total)
        if abs(score) >= 0.15:
            out[token] = round(score, 4)
    return out


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _boundary_windows(seq: str) -> list[int]:
    hi = min(SEARCH_MAX, len(seq) - 3)
    return [pos for pos in range(max(3, SEARCH_MIN), hi + 1) if pos >= 3 and pos + 3 <= len(seq)]


def _score_nterm(seq: str, model: dict[str, object]) -> float:
    nterm = model["nterm"]
    score = 0.0
    score += float(nterm["start_m_log_odds"]["present" if seq[:1] == "M" else "absent"])
    prefix3 = seq[:3] if len(seq) >= 3 else seq
    score += float(nterm["prefix3_log_odds"].get(prefix3, 0.0))
    return score


def _score_boundary(seq: str, pos: int, model: dict[str, object]) -> float:
    if pos < 3 or pos + 3 > len(seq):
        return 0.0
    boundary = model["boundary"]
    weights = model["weights"]
    left = seq[pos - 3 : pos]
    right = seq[pos : pos + 3]
    exact = left + right
    score = 0.0
    score += float(boundary["boundary6_log_odds"].get(exact, 0.0)) * float(weights["boundary6"])
    score += float(boundary["sp_end3_log_odds"].get(left, 0.0)) * float(weights["sp_end3"])
    score += float(boundary["mature3_log_odds"].get(right, 0.0)) * float(weights["mature3"])
    return score


def _best_window_scores(
    rows_pos: list[dict[str, str]],
    rows_neg: list[dict[str, str]],
    model: dict[str, object],
) -> tuple[list[float], list[float], list[float]]:
    pos_true_scores: list[float] = []
    pos_true_margins: list[float] = []
    neg_best_scores: list[float] = []

    for row in rows_pos:
        seq = row["sequence"].strip().upper()
        sp = int(row["sp_length"])
        if seq[:1] != "M" or sp < 3 or sp + 3 > len(seq):
            continue
        scored = sorted(
            ((_score_boundary(seq, pos, model), pos) for pos in _boundary_windows(seq)),
            reverse=True,
        )
        if not scored:
            continue
        best_score, best_pos = scored[0]
        second_score = scored[1][0] if len(scored) > 1 else -999.0
        true_score = _score_boundary(seq, sp, model)
        pos_true_scores.append(true_score)
        if best_pos == sp:
            pos_true_margins.append(true_score - second_score)

    for row in rows_neg:
        seq = row["sequence"].strip().upper()
        if seq[:1] != "M":
            continue
        scores = [_score_boundary(seq, pos, model) for pos in _boundary_windows(seq)]
        if scores:
            neg_best_scores.append(max(scores))

    return pos_true_scores, pos_true_margins, neg_best_scores


def train() -> dict[str, object]:
    pos_rows = _read_csv(GT_CSV)
    neg_rows = _read_csv(NEGATIVE_CONTROL_CSV)

    pos_nterm_prefix3: Counter[str] = Counter()
    neg_nterm_prefix3: Counter[str] = Counter()
    pos_start_m = 0
    neg_start_m = 0

    pos_sp_end3: Counter[str] = Counter()
    neg_sp_end3: Counter[str] = Counter()
    pos_mature3: Counter[str] = Counter()
    neg_mature3: Counter[str] = Counter()
    pos_boundary6: Counter[str] = Counter()
    neg_boundary6: Counter[str] = Counter()
    pos_prefix30_by_sp: dict[str, Counter[int]] = {}
    pos_raw_prefix10: Counter[str] = Counter()
    mature10_counts: Counter[str] = Counter()

    for row in pos_rows:
        seq = row["sequence"].strip().upper()
        sp = int(row["sp_length"])
        pos_start_m += int(seq[:1] == "M")
        if len(seq) >= 3:
            pos_nterm_prefix3[seq[:3]] += 1
        if len(seq) >= EXACT_PREFIX30_LEN:
            prefix30 = seq[:EXACT_PREFIX30_LEN]
            pos_prefix30_by_sp.setdefault(prefix30, Counter())[sp] += 1
        if len(seq) >= EXACT_MATURE10_LEN:
            pos_raw_prefix10[seq[:EXACT_MATURE10_LEN]] += 1
        if sp >= 3 and sp + 3 <= len(seq):
            left = seq[sp - 3 : sp]
            right = seq[sp : sp + 3]
            pos_sp_end3[left] += 1
            pos_mature3[right] += 1
            pos_boundary6[left + right] += 1
        if sp >= 0 and sp + EXACT_MATURE10_LEN <= len(seq):
            mature10_counts[seq[sp : sp + EXACT_MATURE10_LEN]] += 1
        for pos in _boundary_windows(seq):
            if pos == sp or pos < 3 or pos + 3 > len(seq):
                continue
            left = seq[pos - 3 : pos]
            right = seq[pos : pos + 3]
            neg_sp_end3[left] += 1
            neg_mature3[right] += 1
            neg_boundary6[left + right] += 1

    for row in neg_rows:
        seq = row["sequence"].strip().upper()
        neg_start_m += int(seq[:1] == "M")
        if len(seq) >= 3:
            neg_nterm_prefix3[seq[:3]] += 1
        for pos in _boundary_windows(seq):
            left = seq[pos - 3 : pos]
            right = seq[pos : pos + 3]
            neg_sp_end3[left] += 1
            neg_mature3[right] += 1
            neg_boundary6[left + right] += 1

    exact_sp_prefix30: dict[str, dict[str, int]] = {}
    for prefix30, split_counts in pos_prefix30_by_sp.items():
        total = sum(split_counts.values())
        if total < EXACT_SP_PREFIX30_MIN_COUNT or len(split_counts) != 1:
            continue
        split, count = split_counts.most_common(1)[0]
        exact_sp_prefix30[prefix30] = {"split": int(split), "count": int(count)}

    exact_mature10: dict[str, int] = {}
    for prefix10, count in mature10_counts.items():
        if count < EXACT_MATURE10_MIN_COUNT:
            continue
        if pos_raw_prefix10[prefix10] > 0:
            continue
        exact_mature10[prefix10] = int(count)

    model: dict[str, object] = {
        "version": 1,
        "source_gt": GT_CSV.name,
        "source_controls": NEGATIVE_CONTROL_CSV.name,
        "config": {
            "shortcut_search_min": SEARCH_MIN,
            "shortcut_search_max": SEARCH_MAX,
            "sp_shortcut_hfrac_min": 0.65,
            "sp_shortcut_max_c_region": 12,
            "sp_shortcut_min_cleavage_score": 2.5,
            "sp_shortcut_nterm_threshold": 1.0,
            "leaderless_shortcut_max_hfrac": 0.45,
        },
        "weights": {
            "boundary6": 1.0,
            "sp_end3": 0.35,
            "mature3": 0.35,
        },
        "nterm": {
            "start_m_log_odds": {
                "present": round(_log_odds(pos_start_m, len(pos_rows), neg_start_m, len(neg_rows)), 4),
                "absent": round(
                    _log_odds(
                        len(pos_rows) - pos_start_m,
                        len(pos_rows),
                        len(neg_rows) - neg_start_m,
                        len(neg_rows),
                    ),
                    4,
                ),
            },
            "prefix3_log_odds": _table_log_odds(pos_nterm_prefix3, neg_nterm_prefix3),
        },
        "boundary": {
            "sp_end3_log_odds": _table_log_odds(pos_sp_end3, neg_sp_end3),
            "mature3_log_odds": _table_log_odds(pos_mature3, neg_mature3),
            "boundary6_log_odds": _table_log_odds(pos_boundary6, neg_boundary6),
        },
        "exact_shortcuts": {
            "sp_prefix30_len": EXACT_PREFIX30_LEN,
            "sp_prefix30_min_count": EXACT_SP_PREFIX30_MIN_COUNT,
            "sp_prefix30": exact_sp_prefix30,
            "mature10_len": EXACT_MATURE10_LEN,
            "mature10_min_count": EXACT_MATURE10_MIN_COUNT,
            "mature10_prefixes": exact_mature10,
        },
    }

    pos_true_scores, pos_true_margins, neg_best_scores = _best_window_scores(pos_rows, neg_rows, model)
    neg_best_max = max(neg_best_scores) if neg_best_scores else 0.0
    pos_margin_med = statistics.median(pos_true_margins) if pos_true_margins else 1.0

    # Strong shortcut thresholds: keep controls out, then insist on a modest
    # uniqueness margin so the shortcut only fires on clean lexical calls.
    model["config"]["sp_shortcut_score_threshold"] = round(max(5.0, neg_best_max + 0.25), 4)
    model["config"]["sp_shortcut_margin_threshold"] = round(max(1.0, min(2.0, pos_margin_med * 0.5)), 4)

    # Strong leaderless threshold: sit just below the weakest positive
    # N-terminal score so the direct leaderless shortcut remains conservative.
    pos_nterm_scores = sorted(_score_nterm(row["sequence"].strip().upper(), model) for row in pos_rows if row["sequence"])
    leaderless_threshold = (pos_nterm_scores[0] - 0.25) if pos_nterm_scores else -3.0
    model["config"]["leaderless_shortcut_threshold"] = round(leaderless_threshold, 4)

    return model


def main() -> None:
    model = train()
    with OUT_JSON.open("w", encoding="utf-8") as handle:
        json.dump(model, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Wrote {OUT_JSON.relative_to(ROOT)}")
    print(f"  prefix3 terms:   {len(model['nterm']['prefix3_log_odds'])}")
    print(f"  sp_end3 terms:   {len(model['boundary']['sp_end3_log_odds'])}")
    print(f"  mature3 terms:   {len(model['boundary']['mature3_log_odds'])}")
    print(f"  boundary6 terms: {len(model['boundary']['boundary6_log_odds'])}")
    print(f"  exact sp-prefix30 shortcuts: {len(model['exact_shortcuts']['sp_prefix30'])}")
    print(f"  exact mature10 shortcuts:    {len(model['exact_shortcuts']['mature10_prefixes'])}")
    print(
        "  shortcut thresholds:"
        f" score>={model['config']['sp_shortcut_score_threshold']},"
        f" margin>={model['config']['sp_shortcut_margin_threshold']},"
        f" leaderless<={model['config']['leaderless_shortcut_threshold']}"
    )


if __name__ == "__main__":
    main()
