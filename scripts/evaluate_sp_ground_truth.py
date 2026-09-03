#!/usr/bin/env python3
"""Evaluate signal peptide heuristic accuracy against UniProt ground truth.

Prefers ``data/sp_ground_truth_enriched.csv`` when available.  That enriched
benchmark carries gold MHC class / chain metadata so the evaluator can dispatch
the parser directly instead of guessing class from the sequence alone. Any
classless rows use the same whole-parse competition as the production library.
Curated/gold and unresolved/inferred results are reported separately.
Species-category strata come from the checked-in UniProt taxonomy audit rather
than open-ended organism-name hints.

Usage:
    python scripts/evaluate_sp_ground_truth.py
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.sp_ground_truth_taxonomy import species_category as _species_category

GT_RAW_CSV = ROOT / "data" / "sp_ground_truth.csv"
GT_ENRICHED_CSV = ROOT / "data" / "sp_ground_truth_enriched.csv"
NEGATIVE_CONTROL_CSV = ROOT / "data" / "sp_negative_controls.csv"
GT_CSV = GT_ENRICHED_CSV if GT_ENRICHED_CSV.exists() else GT_RAW_CSV


def _parse_cli_args(argv: list[str]) -> dict[str, bool]:
    """Return simple evaluator options from argv."""
    return {
        "use_early_shortcuts": "--no-early-shortcuts" not in argv,
    }


def load_ground_truth_rows(prefer_enriched: bool = True) -> tuple[Path, list[dict[str, str]]]:
    """Load raw or enriched GT rows, preferring the enriched file when present."""
    path = GT_ENRICHED_CSV if prefer_enriched and GT_ENRICHED_CSV.exists() else GT_RAW_CSV
    with open(path, "r", encoding="utf-8") as handle:
        return path, list(csv.DictReader(handle))


def _row_species_category(row: dict[str, str]) -> str:
    return row.get("species_category", "") or _species_category(
        row["organism"],
        row.get("taxon_id", ""),
        row.get("source_clade", ""),
    )


def _row_dispatch_metadata(row: dict[str, str]) -> tuple[str, str, str]:
    mhc_class = str(row.get("mhc_class", "") or "").strip().upper()
    chain = str(row.get("chain", "") or "").strip().lower()
    gene = str(row.get("gene", "") or "").strip()
    return mhc_class, chain, gene


def _row_benchmark_stratum(row: dict[str, str]) -> str:
    """Return the source-confidence stratum used in benchmark denominators."""
    label_status = str(row.get("label_status", "") or "").strip().lower()
    mhc_class, chain, _gene = _row_dispatch_metadata(row)
    if label_status.startswith("excluded_"):
        return "excluded"
    if label_status in {"gold", "curated"} and mhc_class in {"I", "II"} and chain in {"alpha", "beta"}:
        return "curated/gold"
    return "unresolved/inferred"


def _parser_name_for_dispatch(mhc_class: str, chain: str) -> str:
    if mhc_class == "I":
        return "class_I"
    if mhc_class == "II" and chain == "alpha":
        return "class_II_alpha"
    if mhc_class == "II" and chain == "beta":
        return "class_II_beta"
    return ""


def _try_parse(seq: str, *, features=None, use_early_shortcuts: bool = True) -> tuple[int, str]:
    """Run the production classless whole-parse competition.

    Returns ``(mature_start, parser_used)`` or ``(0, "")`` on failure.
    """
    from mhcseqs.domain_parsing import analyze_sequence, decompose_domains

    if features is None:
        features = analyze_sequence(seq)
    result = decompose_domains(
        seq,
        mhc_class="",
        features=features,
        use_early_shortcuts=use_early_shortcuts,
    )
    if not result.ok or result.mature_start <= 0:
        return 0, ""
    return int(result.mature_start), _parser_name_for_dispatch(result.mhc_class, result.chain)


def predict_sp_for_row(
    row: dict[str, str],
    *,
    use_early_shortcuts: bool = True,
) -> dict[str, str | int | bool]:
    """Predict SP length for a GT row, using gold dispatch when available."""
    from mhcseqs.domain_parsing import analyze_sequence, decompose_domains, refine_signal_peptide

    if _row_benchmark_stratum(row) == "excluded":
        return {
            "ok": False,
            "dispatch_mode": "excluded",
            "parser": "",
            "mhc_class": "",
            "chain": "",
            "gene": str(row.get("gene", "") or ""),
            "status": str(row.get("label_status", "") or "excluded"),
            "mature_start": 0,
            "predicted_sp": 0,
        }

    seq = row["sequence"]
    category = _row_species_category(row)
    mhc_class, chain, gene = _row_dispatch_metadata(row)
    features = analyze_sequence(seq)

    has_complete_dispatch = mhc_class == "I" or (mhc_class == "II" and chain in {"alpha", "beta"})
    if has_complete_dispatch:
        parser_name = _parser_name_for_dispatch(mhc_class, chain)
        try:
            result = decompose_domains(
                seq,
                mhc_class=mhc_class,
                chain=chain or None,
                gene=gene,
                species=row.get("organism", ""),
                features=features,
                use_early_shortcuts=use_early_shortcuts,
            )
        except Exception as exc:  # pragma: no cover - defensive evaluator path
            return {
                "ok": False,
                "dispatch_mode": "gold",
                "parser": parser_name,
                "mhc_class": mhc_class,
                "chain": chain,
                "gene": gene,
                "status": f"exception:{type(exc).__name__}",
                "mature_start": 0,
                "predicted_sp": 0,
            }

        resolved_chain = str(result.chain or chain or "")
        resolved_class = str(result.mhc_class or mhc_class or "")
        if not parser_name:
            parser_name = _parser_name_for_dispatch(resolved_class, resolved_chain)

        if not result.ok or result.mature_start <= 0:
            return {
                "ok": False,
                "dispatch_mode": "gold",
                "parser": parser_name,
                "mhc_class": resolved_class,
                "chain": resolved_chain,
                "gene": gene,
                "status": result.status,
                "mature_start": int(result.mature_start or 0),
                "predicted_sp": 0,
            }

        groove_anchor = (
            (int(result.anchor_cys1), int(result.anchor_cys2)) if result.anchor_cys1 is not None and result.anchor_cys2 is not None else None
        )
        refined = refine_signal_peptide(
            seq,
            result.mature_start,
            category,
            mhc_class,
            features=features,
            groove_anchor=groove_anchor,
        )
        return {
            "ok": True,
            "dispatch_mode": "gold",
            "parser": parser_name,
            "mhc_class": resolved_class,
            "chain": resolved_chain,
            "gene": gene,
            "status": result.status,
            "mature_start": int(result.mature_start),
            "predicted_sp": int(refined),
        }
    try:
        inferred_class_hint = mhc_class if mhc_class in {"I", "II"} else ""
        inferred_chain_hint = chain if inferred_class_hint == "II" and chain else None
        result = decompose_domains(
            seq,
            mhc_class=inferred_class_hint,
            chain=inferred_chain_hint,
            gene=gene,
            species=row.get("organism", ""),
            features=features,
            use_early_shortcuts=use_early_shortcuts,
        )
    except Exception as exc:  # pragma: no cover - defensive evaluator path
        return {
            "ok": False,
            "dispatch_mode": "inferred",
            "parser": "",
            "mhc_class": "",
            "chain": "",
            "gene": gene,
            "status": f"exception:{type(exc).__name__}",
            "mature_start": 0,
            "predicted_sp": 0,
        }
    if not result.ok or result.mature_start <= 0:
        return {
            "ok": False,
            "dispatch_mode": "inferred",
            "parser": _parser_name_for_dispatch(result.mhc_class, result.chain),
            "mhc_class": str(result.mhc_class or ""),
            "chain": str(result.chain or ""),
            "gene": gene,
            "status": result.status,
            "mature_start": int(result.mature_start or 0),
            "predicted_sp": 0,
        }

    inferred_class = str(result.mhc_class or inferred_class_hint or "")
    inferred_chain = str(result.chain or "")
    parser_name = _parser_name_for_dispatch(inferred_class, inferred_chain)
    groove_anchor = (
        (int(result.anchor_cys1), int(result.anchor_cys2))
        if result.anchor_cys1 is not None and result.anchor_cys2 is not None
        else None
    )
    refined = refine_signal_peptide(
        seq,
        result.mature_start,
        category,
        inferred_class,
        features=features,
        groove_anchor=groove_anchor,
    )
    return {
        "ok": True,
        "dispatch_mode": "inferred",
        "parser": parser_name,
        "mhc_class": inferred_class,
        "chain": inferred_chain,
        "gene": gene,
        "status": result.status,
        "mature_start": int(result.mature_start),
        "predicted_sp": int(refined),
    }


def evaluate(*, use_early_shortcuts: bool = True):
    if not GT_CSV.exists():
        print(f"Ground truth not found: {GT_CSV}", file=sys.stderr)
        print("Run: python scripts/fetch_sp_ground_truth.py", file=sys.stderr)
        sys.exit(1)

    gt_path, rows = load_ground_truth_rows(prefer_enriched=True)

    print(f"Loaded {len(rows)} ground truth entries from {gt_path.name}")
    print(f"Early shortcuts:  {'enabled' if use_early_shortcuts else 'disabled'}")
    print()

    # Counters
    total = 0
    excluded = 0
    parsed = 0
    unparsed = 0
    exact = 0
    within_1 = 0
    within_2 = 0
    within_3 = 0
    deltas = []
    by_category: dict[str, Counter] = {}
    by_class: dict[tuple[str, str], Counter] = {}
    by_reviewed: dict[str, Counter] = {}
    by_dispatch: Counter = Counter()
    by_stratum: dict[str, Counter] = {
        "curated/gold": Counter(),
        "unresolved/inferred": Counter(),
    }
    mismatches_gt3: list[dict] = []

    for row in rows:
        stratum = _row_benchmark_stratum(row)
        if stratum == "excluded":
            excluded += 1
            continue
        gt_sp_len = int(row["sp_length"])
        organism = row["organism"]
        reviewed = row["reviewed"]
        accession = row["accession"]
        total += 1
        by_stratum[stratum]["total"] += 1

        prediction = predict_sp_for_row(row, use_early_shortcuts=use_early_shortcuts)
        by_dispatch[str(prediction["dispatch_mode"])] += 1
        if not prediction["ok"]:
            unparsed += 1
            by_stratum[stratum]["unparsed"] += 1
            continue

        parsed += 1
        by_stratum[stratum]["parsed"] += 1

        cat = _row_species_category(row)
        refined = int(prediction["predicted_sp"])

        delta = refined - gt_sp_len
        deltas.append(delta)

        # Accuracy bins
        if delta == 0:
            exact += 1
            by_stratum[stratum]["exact"] += 1
        if abs(delta) <= 1:
            within_1 += 1
            by_stratum[stratum]["within_1"] += 1
        if abs(delta) <= 2:
            within_2 += 1
            by_stratum[stratum]["within_2"] += 1
        if abs(delta) <= 3:
            within_3 += 1
            by_stratum[stratum]["within_3"] += 1

        if abs(delta) > 3:
            mismatches_gt3.append(
                {
                    "accession": accession,
                    "organism": organism,
                    "category": cat or "unknown",
                    "gt_sp": gt_sp_len,
                    "cys_pred": int(prediction["mature_start"]),
                    "refined": refined,
                    "delta": delta,
                    "reviewed": reviewed,
                    "parser": str(prediction["parser"]),
                    "mhc_class": str(prediction["mhc_class"]),
                    "chain": str(prediction["chain"]),
                    "dispatch_mode": str(prediction["dispatch_mode"]),
                }
            )

        # Per-category stats
        cat_key = cat or "unknown"
        if cat_key not in by_category:
            by_category[cat_key] = Counter()
        by_category[cat_key]["total"] += 1
        if delta == 0:
            by_category[cat_key]["exact"] += 1
        if abs(delta) <= 1:
            by_category[cat_key]["within_1"] += 1
        if abs(delta) <= 2:
            by_category[cat_key]["within_2"] += 1

        class_key = str(row.get("mhc_class", "") or prediction["mhc_class"] or "unknown")
        if (class_key, cat_key) not in by_class:
            by_class[(class_key, cat_key)] = Counter()
        by_class[(class_key, cat_key)]["total"] += 1
        if delta == 0:
            by_class[(class_key, cat_key)]["exact"] += 1
        if abs(delta) <= 1:
            by_class[(class_key, cat_key)]["within_1"] += 1
        if abs(delta) <= 2:
            by_class[(class_key, cat_key)]["within_2"] += 1
        if abs(delta) <= 3:
            by_class[(class_key, cat_key)]["within_3"] += 1

        # Per-reviewed stats
        rev_key = "reviewed" if reviewed == "Y" else "unreviewed"
        if rev_key not in by_reviewed:
            by_reviewed[rev_key] = Counter()
        by_reviewed[rev_key]["total"] += 1
        if delta == 0:
            by_reviewed[rev_key]["exact"] += 1
        if abs(delta) <= 1:
            by_reviewed[rev_key]["within_1"] += 1

    # Report
    print("=" * 70)
    print("OVERALL RESULTS")
    print("=" * 70)
    print(f"Total entries:    {total}")
    print(f"Excluded non-MHC: {excluded}")
    print(f"Parsed:           {parsed} ({100 * parsed / total:.1f}%)")
    print(f"Unparsed:         {unparsed} ({100 * unparsed / total:.1f}%)")
    if by_dispatch:
        dispatch_summary = ", ".join(f"{k}={v}" for k, v in sorted(by_dispatch.items()))
        print(f"Dispatch mode:    {dispatch_summary}")
    print()
    if parsed > 0:
        print(f"Exact match:      {exact}/{parsed} ({100 * exact / parsed:.1f}%)")
        print(f"Within +/-1 aa:   {within_1}/{parsed} ({100 * within_1 / parsed:.1f}%)")
        print(f"Within +/-2 aa:   {within_2}/{parsed} ({100 * within_2 / parsed:.1f}%)")
        print(f"Within +/-3 aa:   {within_3}/{parsed} ({100 * within_3 / parsed:.1f}%)")

    if deltas:
        import statistics

        print(f"\nMean delta:       {statistics.mean(deltas):+.2f} aa")
        print(f"Median delta:     {statistics.median(deltas):+.1f} aa")
        print(f"Std dev:          {statistics.stdev(deltas):.2f} aa")

    print("\n" + "=" * 70)
    print("BY CURATION STRATUM")
    print("=" * 70)
    print(f"{'Stratum':<22} {'Total':>6} {'Parsed':>12} {'Exact':>12} {'<=1':>12} {'<=2':>12}")
    print("-" * 80)
    for stratum in ("curated/gold", "unresolved/inferred"):
        c = by_stratum[stratum]
        t = c["total"]
        p = c["parsed"]
        parsed_text = f"{p:>4} ({100 * p / t:5.1f}%)" if t else "   0 (  n/a)"
        exact_text = f"{c['exact']:>4} ({100 * c['exact'] / p:5.1f}%)" if p else "   0 (  n/a)"
        one_text = f"{c['within_1']:>4} ({100 * c['within_1'] / p:5.1f}%)" if p else "   0 (  n/a)"
        two_text = f"{c['within_2']:>4} ({100 * c['within_2'] / p:5.1f}%)" if p else "   0 (  n/a)"
        print(f"{stratum:<22} {t:>6} {parsed_text} {exact_text} {one_text} {two_text}")

    # Delta distribution
    delta_counts = Counter(deltas)
    print("\nDelta distribution (predicted - ground truth):")
    for d in sorted(delta_counts.keys()):
        bar = "#" * min(delta_counts[d], 60)
        print(f"  {d:+3d}: {delta_counts[d]:>4d} {bar}")

    # Per-category breakdown
    print("\n" + "=" * 70)
    print("BY SPECIES CATEGORY")
    print("=" * 70)
    print(f"{'Category':<20} {'Total':>6} {'Exact':>8} {'<=1':>8} {'<=2':>8}")
    print("-" * 56)
    for cat in sorted(by_category.keys()):
        c = by_category[cat]
        t = c["total"]
        print(
            f"{cat:<20} {t:>6} "
            f"{c['exact']:>4} ({100 * c['exact'] / t:5.1f}%) "
            f"{c['within_1']:>4} ({100 * c['within_1'] / t:5.1f}%) "
            f"{c['within_2']:>4} ({100 * c['within_2'] / t:5.1f}%)"
        )

    if by_class:
        print("\n" + "=" * 70)
        print("BY MHC CLASS AND SPECIES CATEGORY")
        print("=" * 70)
        print(f"{'Class':<8} {'Category':<18} {'Total':>6} {'Exact':>12} {'<=1':>12} {'<=2':>12} {'<=3':>12}")
        print("-" * 82)
        for mhc_class, category in sorted(by_class.keys()):
            c = by_class[(mhc_class, category)]
            t = c["total"]
            print(
                f"{mhc_class:<8} {category:<18} {t:>6} "
                f"{c['exact']:>4} ({100 * c['exact'] / t:5.1f}%) "
                f"{c['within_1']:>4} ({100 * c['within_1'] / t:5.1f}%) "
                f"{c['within_2']:>4} ({100 * c['within_2'] / t:5.1f}%) "
                f"{c['within_3']:>4} ({100 * c['within_3'] / t:5.1f}%)"
            )

    # Per-reviewed breakdown
    print(f"\n{'Review status':<20} {'Total':>6} {'Exact':>8} {'<=1':>8}")
    print("-" * 48)
    for rev in sorted(by_reviewed.keys()):
        c = by_reviewed[rev]
        t = c["total"]
        print(f"{rev:<20} {t:>6} {c['exact']:>4} ({100 * c['exact'] / t:5.1f}%) {c['within_1']:>4} ({100 * c['within_1'] / t:5.1f}%)")

    # Worst mismatches
    if mismatches_gt3:
        print(f"\n{'=' * 70}")
        print(f"MISMATCHES > 3 aa ({len(mismatches_gt3)} entries)")
        print(f"{'=' * 70}")
        print(f"{'Accession':<12} {'Category':<18} {'GT':>4} {'Cys':>5} {'Ref':>5} {'Delta':>6} {'Rev'} {'Organism'}")
        print("-" * 90)
        for m in sorted(mismatches_gt3, key=lambda x: abs(x["delta"]), reverse=True)[:30]:
            org_short = m["organism"][:30]
            print(
                f"{m['accession']:<12} {m['category']:<18} {m['gt_sp']:>4} "
                f"{m['cys_pred']:>5} {m['refined']:>5} {m['delta']:>+6} "
                f"{m['reviewed']}   {org_short}"
            )
        if len(mismatches_gt3) > 30:
            print(f"  ... and {len(mismatches_gt3) - 30} more")


if __name__ == "__main__":
    options = _parse_cli_args(sys.argv[1:])
    evaluate(use_early_shortcuts=options["use_early_shortcuts"])
