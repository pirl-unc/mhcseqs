#!/usr/bin/env python3
"""Train a motif-composition model for signal-peptide cleavage boundaries.

The model is deliberately simple and strongly regularized:
- exact residue log-odds around the boundary
- coarse amino-acid-class log-odds around the boundary
- positives = annotated SP boundaries from data/sp_ground_truth.csv
- negatives = nearby same-sequence decoys within +/-10 aa

The output is written to mhcseqs/sp_boundary_model.json and can be consumed by
mhcseqs.domain_parsing._score_sp_boundary_composition().
"""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mhcseqs.domain_parsing import _refinement_group
from scripts.evaluate_sp_ground_truth import GT_CSV, _row_benchmark_stratum, _row_species_category

OUT_JSON = ROOT / "mhcseqs" / "sp_boundary_model.json"
AA_ORDER = tuple("ACDEFGHIKLMNPQRSTVWY")
OFFSETS = tuple(range(-5, 6))
ALPHA = 1.0
NEGATIVE_RADIUS = 10

AA_CLASS_BY_AA = {
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
    "G": "glycine",
    "P": "proline",
    "C": "cysteine",
}
AA_CLASS_ORDER = tuple(sorted(set(AA_CLASS_BY_AA.values())))


def _new_group_counts() -> dict[str, object]:
    return {
        "n_positive": 0,
        "n_negative": 0,
        "positive_residue": {offset: Counter() for offset in OFFSETS},
        "negative_residue": {offset: Counter() for offset in OFFSETS},
        "positive_class": {offset: Counter() for offset in OFFSETS},
        "negative_class": {offset: Counter() for offset in OFFSETS},
    }


def _observe_boundary(
    seq: str,
    boundary: int,
    payload: dict[str, object],
    *,
    positive: bool,
) -> None:
    if boundary < 0 or boundary >= len(seq):
        return
    residue_key = "positive_residue" if positive else "negative_residue"
    class_key = "positive_class" if positive else "negative_class"
    count_key = "n_positive" if positive else "n_negative"
    payload[count_key] += 1
    residue_counts = payload[residue_key]
    class_counts = payload[class_key]
    for offset in OFFSETS:
        pos = boundary + offset
        if not (0 <= pos < len(seq)):
            continue
        aa = seq[pos]
        if aa not in AA_CLASS_BY_AA:
            continue
        residue_counts[offset][aa] += 1
        class_counts[offset][AA_CLASS_BY_AA[aa]] += 1


def _log_odds_table(
    positive: dict[int, Counter],
    negative: dict[int, Counter],
    alphabet: tuple[str, ...],
) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for offset in OFFSETS:
        pos_counter = positive[offset]
        neg_counter = negative[offset]
        pos_total = sum(pos_counter.values())
        neg_total = sum(neg_counter.values())
        pos_den = pos_total + ALPHA * len(alphabet)
        neg_den = neg_total + ALPHA * len(alphabet)
        out[str(offset)] = {}
        for token in alphabet:
            pos_freq = (pos_counter.get(token, 0) + ALPHA) / max(pos_den, 1.0)
            neg_freq = (neg_counter.get(token, 0) + ALPHA) / max(neg_den, 1.0)
            score = math.log(pos_freq / neg_freq)
            score = max(-2.5, min(2.5, score))
            out[str(offset)][token] = round(score, 4)
    return out


def train(path: Path = GT_CSV) -> dict[str, object]:
    """Train only from source-backed curated/gold MHC alpha/beta rows."""
    with path.open("r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    strata = Counter(_row_benchmark_stratum(row) for row in rows)
    eligible_rows = [row for row in rows if _row_benchmark_stratum(row) == "curated/gold"]

    counts = defaultdict(_new_group_counts)
    counts["global"] = _new_group_counts()

    used_rows = 0
    for row in eligible_rows:
        seq = str(row["sequence"]).strip().upper()
        boundary = int(row["sp_length"])
        if boundary < 5 or boundary >= len(seq) - 5:
            continue
        used_rows += 1

        species_group = _refinement_group(_row_species_category(row))
        for key in ("global", species_group):
            _observe_boundary(seq, boundary, counts[key], positive=True)

        neg_lo = max(5, boundary - NEGATIVE_RADIUS)
        neg_hi = min(len(seq) - 5, boundary + NEGATIVE_RADIUS + 1)
        for neg_boundary in range(neg_lo, neg_hi):
            if neg_boundary == boundary:
                continue
            for key in ("global", species_group):
                _observe_boundary(seq, neg_boundary, counts[key], positive=False)

    payload: dict[str, object] = {
        "version": 2,
        "window_offsets": list(OFFSETS),
        "negative_radius": NEGATIVE_RADIUS,
        "training": {
            "source": path.name,
            "policy": "label_status in {gold,curated} and MHC class/chain in {I alpha,II alpha,II beta}",
            "source_rows": len(rows),
            "eligible_rows": len(eligible_rows),
            "used_rows": used_rows,
            "excluded_rows": strata["excluded"],
            "unresolved_rows": strata["unresolved/inferred"],
        },
        "groups": {},
    }
    for group, group_counts in counts.items():
        payload["groups"][group] = {
            "n_positive": group_counts["n_positive"],
            "n_negative": group_counts["n_negative"],
            "residue_log_odds": _log_odds_table(
                group_counts["positive_residue"],
                group_counts["negative_residue"],
                AA_ORDER,
            ),
            "class_log_odds": _log_odds_table(
                group_counts["positive_class"],
                group_counts["negative_class"],
                AA_CLASS_ORDER,
            ),
        }
    return payload


def main() -> None:
    model = train()
    with OUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(model, f, indent=2, sort_keys=True)
        f.write("\n")
    print(f"Wrote {OUT_JSON.relative_to(ROOT)}")
    training = model["training"]
    print(
        "Training rows:"
        f" source={training['source_rows']} eligible={training['eligible_rows']} used={training['used_rows']}"
        f" excluded={training['excluded_rows']} unresolved={training['unresolved_rows']}"
    )
    for group, payload in model["groups"].items():
        print(f"{group:>16}  positives={payload['n_positive']:4d}  negatives={payload['n_negative']:5d}")


if __name__ == "__main__":
    main()
