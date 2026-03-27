#!/usr/bin/env python3
"""Analyze SP ground-truth performance by gold MHC class and failure mode.

This complements ``scripts/evaluate_sp_ground_truth.py``.

When ``data/sp_ground_truth_enriched.csv`` is present, rows are dispatched by
their gold MHC class / chain labels.  Any remaining unlabeled rows fall back to
the older classless parser-selection heuristic, but the analysis reports which
mode was used.
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from evaluate_sp_ground_truth import (
    GT_CSV,
    _row_dispatch_metadata,
    _row_species_category,
    load_ground_truth_rows,
    predict_sp_for_row,
)


def _run_all_parsers(seq: str, species_category: str) -> dict[str, dict]:
    from mhcseqs.domain_parsing import (
        decompose_class_i,
        decompose_class_ii_alpha,
        decompose_class_ii_beta,
        refine_signal_peptide,
    )

    parser_specs = (
        ("class_I", "I", decompose_class_i),
        ("class_II_beta", "II", decompose_class_ii_beta),
        ("class_II_alpha", "II", decompose_class_ii_alpha),
    )

    runs: dict[str, dict] = {}
    for name, mhc_class, parser in parser_specs:
        try:
            result = parser(seq)
            refined = 0
            if result.ok and result.mature_start > 0:
                refined = refine_signal_peptide(seq, result.mature_start, species_category, mhc_class)
            runs[name] = {
                "ok": result.ok,
                "status": result.status,
                "mature_start": result.mature_start,
                "refined": refined,
                "parse_score": getattr(result, "parse_score", 0.0),
                "anchor_type": result.anchor_type,
                "flags": tuple(result.flags),
            }
        except Exception as exc:  # pragma: no cover - defensive analysis path
            runs[name] = {
                "ok": False,
                "status": f"exception:{type(exc).__name__}",
                "mature_start": 0,
                "refined": 0,
                "parse_score": 0.0,
                "anchor_type": "",
                "flags": (),
            }
    return runs


def _chosen_class(parser_name: str) -> str:
    return "I" if parser_name == "class_I" else "II"


def _select_parser_from_runs(runs: dict[str, dict]) -> tuple[int, str]:
    """Replicate the evaluator's parser selection without reparsing."""
    typical_sp = 23
    candidates = [
        (int(info["mature_start"]), name)
        for name, info in runs.items()
        if info["ok"] and int(info["mature_start"] or 0) > 0
    ]
    if not candidates:
        return 0, ""

    in_range = [(ms, name) for ms, name in candidates if 10 <= ms <= 50]
    if in_range:
        return min(in_range, key=lambda item: abs(item[0] - typical_sp))
    return min(candidates, key=lambda item: item[0])


def _failure_mode(
    *,
    gt_sp: int,
    chosen_parser: str,
    chosen: dict,
    runs: dict[str, dict],
) -> str:
    chosen_refined = int(chosen["refined"] or 0)
    chosen_raw = int(chosen["mature_start"] or 0)
    chosen_abs = abs(chosen_refined - gt_sp)

    best_alt_name = ""
    best_alt_abs = chosen_abs
    for name, info in runs.items():
        if name == chosen_parser or not info["ok"] or not info["refined"]:
            continue
        alt_abs = abs(int(info["refined"]) - gt_sp)
        if alt_abs < best_alt_abs:
            best_alt_abs = alt_abs
            best_alt_name = name

    if best_alt_name and best_alt_abs + 2 <= chosen_abs:
        return f"better_alt_parser->{best_alt_name}"
    if abs(chosen_raw - gt_sp) <= 3 < chosen_abs:
        return "refinement_regression"
    if chosen["status"] in {
        "inferred_from_alpha3",
        "beta1_only_fallback",
        "fragment_fallback",
        "alpha1_only",
        "alpha2_only",
    }:
        return f"{chosen['status']}_miss"
    if chosen_refined - gt_sp > 3:
        return "overcall"
    if chosen_refined - gt_sp < -3:
        return "undercall"
    return "other_large_miss"


def main() -> None:
    if not GT_CSV.exists():
        raise SystemExit(f"Missing ground truth CSV: {GT_CSV}")

    gt_path, rows = load_ground_truth_rows(prefer_enriched=True)

    by_class_cat: dict[tuple[str, str], Counter] = defaultdict(Counter)
    by_parser_cat: dict[tuple[str, str], Counter] = defaultdict(Counter)
    by_dispatch = Counter()
    failure_modes = Counter()
    failure_examples: dict[str, list[dict]] = defaultdict(list)
    unparsed_signatures = Counter()
    unparsed_examples: dict[str, list[dict]] = defaultdict(list)

    parsed = 0
    exact = 0

    for row in rows:
        seq = row["sequence"]
        gt_sp = int(row["sp_length"])
        cat = _row_species_category(row)
        prediction = predict_sp_for_row(row)
        runs: dict[str, dict] = {}

        if not prediction["ok"]:
            mhc_class, chain, _gene = _row_dispatch_metadata(row)
            if mhc_class:
                signature = f"gold:{mhc_class}:{chain or 'na'}:{prediction['status']}"
            else:
                runs = _run_all_parsers(seq, cat)
                signature = " | ".join(
                    f"{name}:{runs[name]['status']}" for name in ("class_I", "class_II_alpha", "class_II_beta")
                )
            unparsed_signatures[signature] += 1
            if len(unparsed_examples[signature]) < 3:
                unparsed_examples[signature].append(
                    {
                        "accession": row["accession"],
                        "organism": row["organism"],
                        "category": cat,
                        "gt_sp": gt_sp,
                    }
                )
            continue

        chosen_parser = str(prediction["parser"])
        chosen_start = int(prediction["mature_start"])
        pred = int(prediction["predicted_sp"])
        delta = pred - gt_sp
        parsed += 1
        by_dispatch[str(prediction["dispatch_mode"])] += 1
        if delta == 0:
            exact += 1

        class_key = row.get("mhc_class", "") or str(prediction["mhc_class"]) or _chosen_class(chosen_parser)
        chosen = runs.get(chosen_parser, {
            "status": str(prediction["status"]),
            "mature_start": chosen_start,
            "refined": pred,
            "ok": bool(prediction["ok"]),
        })
        chosen["mature_start"] = chosen_start
        chosen["refined"] = pred
        chosen["status"] = str(prediction["status"])
        by_class_cat[(class_key, cat)]["total"] += 1
        by_parser_cat[(chosen_parser, cat)]["total"] += 1
        if delta == 0:
            by_class_cat[(class_key, cat)]["exact"] += 1
            by_parser_cat[(chosen_parser, cat)]["exact"] += 1
        if abs(delta) <= 1:
            by_class_cat[(class_key, cat)]["within_1"] += 1
            by_parser_cat[(chosen_parser, cat)]["within_1"] += 1
        if abs(delta) <= 2:
            by_class_cat[(class_key, cat)]["within_2"] += 1
            by_parser_cat[(chosen_parser, cat)]["within_2"] += 1
        if abs(delta) <= 3:
            by_class_cat[(class_key, cat)]["within_3"] += 1
            by_parser_cat[(chosen_parser, cat)]["within_3"] += 1

        if abs(delta) > 3:
            runs = _run_all_parsers(seq, cat)
            mode = _failure_mode(gt_sp=gt_sp, chosen_parser=chosen_parser, chosen=chosen, runs=runs)
            failure_modes[mode] += 1
            if len(failure_examples[mode]) < 5:
                failure_examples[mode].append(
                    {
                        "accession": row["accession"],
                        "organism": row["organism"],
                        "category": cat,
                        "gt_sp": gt_sp,
                        "pred": pred,
                        "delta": delta,
                        "parser": chosen_parser,
                        "status": chosen["status"],
                        "dispatch_mode": prediction["dispatch_mode"],
                    }
                )

    print(f"Loaded {len(rows)} rows from {gt_path.name}")
    print(f"Parsed: {parsed}/{len(rows)} ({100 * parsed / len(rows):.1f}%)")
    print(f"Exact:  {exact}/{parsed} ({100 * exact / parsed:.1f}%)" if parsed else "Exact: n/a")
    if by_dispatch:
        print("Dispatch mode: " + ", ".join(f"{k}={v}" for k, v in sorted(by_dispatch.items())))

    print("\nBy class and species category")
    print(f"{'Class':<6} {'Category':<18} {'Total':>6} {'Exact':>12} {'<=1':>12} {'<=2':>12} {'<=3':>12}")
    print("-" * 82)
    for class_key, category in sorted(by_class_cat.keys()):
        c = by_class_cat[(class_key, category)]
        total = c["total"]
        print(
            f"{class_key:<6} {category:<18} {total:>6} "
            f"{c['exact']:>4} ({100*c['exact']/total:5.1f}%) "
            f"{c['within_1']:>4} ({100*c['within_1']/total:5.1f}%) "
            f"{c['within_2']:>4} ({100*c['within_2']/total:5.1f}%) "
            f"{c['within_3']:>4} ({100*c['within_3']/total:5.1f}%)"
        )

    print("\nBy selected parser and species category")
    print(f"{'Parser':<16} {'Category':<18} {'Total':>6} {'Exact':>12}")
    print("-" * 58)
    for parser_name, category in sorted(by_parser_cat.keys()):
        c = by_parser_cat[(parser_name, category)]
        total = c["total"]
        print(
            f"{parser_name:<16} {category:<18} {total:>6} "
            f"{c['exact']:>4} ({100*c['exact']/total:5.1f}%)"
        )

    print("\nLarge-miss failure modes (>3 aa)")
    if not failure_modes:
        print("None")
    else:
        for mode, count in failure_modes.most_common():
            print(f"- {mode}: {count}")
            for example in failure_examples[mode]:
                print(
                    f"  {example['accession']} {example['category']} {example['parser']} "
                    f"status={example['status']} gt={example['gt_sp']} pred={example['pred']} "
                    f"delta={example['delta']:+d} dispatch={example['dispatch_mode']} "
                    f"{example['organism'][:40]}"
                )

    print("\nUnparsed status signatures")
    if not unparsed_signatures:
        print("None")
    else:
        for signature, count in unparsed_signatures.most_common(10):
            print(f"- {count}: {signature}")
            for example in unparsed_examples[signature]:
                print(
                    f"  {example['accession']} {example['category']} gt={example['gt_sp']} "
                    f"{example['organism'][:40]}"
                )


if __name__ == "__main__":
    main()
