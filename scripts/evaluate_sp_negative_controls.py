#!/usr/bin/env python3
"""Evaluate SP false positives on mature-only and fragment control sequences."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
import sys

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from evaluate_sp_ground_truth import NEGATIVE_CONTROL_CSV, _parse_cli_args, _row_species_category


def main(*, use_early_shortcuts: bool = True) -> None:
    from mhcseqs.domain_parsing import analyze_sequence, decompose_domains, refine_signal_peptide

    if not NEGATIVE_CONTROL_CSV.exists():
        raise SystemExit(
            f"Missing control CSV: {NEGATIVE_CONTROL_CSV}\n"
            "Run: python scripts/enrich_sp_ground_truth.py"
        )

    with open(NEGATIVE_CONTROL_CSV, "r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    by_kind_cat: dict[tuple[str, str], Counter] = defaultdict(Counter)
    totals = Counter()

    for row in rows:
        mhc_class = row.get("mhc_class", "")
        chain = row.get("chain", "")
        gene = row.get("gene", "")
        category = _row_species_category(row)
        kind = row.get("control_type", "unknown")
        totals["total"] += 1
        by_kind_cat[(kind, category)]["total"] += 1
        features = analyze_sequence(row["sequence"])

        try:
            result = decompose_domains(
                row["sequence"],
                mhc_class=mhc_class,
                chain=chain or None,
                gene=gene,
                features=features,
                use_early_shortcuts=use_early_shortcuts,
            )
        except Exception:
            totals["abstain"] += 1
            by_kind_cat[(kind, category)]["abstain"] += 1
            continue

        if result.ok and result.mature_start > 0:
            groove_anchor = (
                (int(result.anchor_cys1), int(result.anchor_cys2))
                if result.anchor_cys1 is not None and result.anchor_cys2 is not None
                else None
            )
            predicted = refine_signal_peptide(
                row["sequence"],
                result.mature_start,
                category,
                mhc_class,
                features=features,
                groove_anchor=groove_anchor,
            )
        else:
            predicted = int(result.mature_start or 0)

        if predicted == 0:
            bucket = "correct_zero_sp" if result.ok else "abstain_zero_sp"
        else:
            bucket = "false_positive_sp"

        totals[bucket] += 1
        by_kind_cat[(kind, category)][bucket] += 1

    print(f"Loaded {len(rows)} controls from {NEGATIVE_CONTROL_CSV.name}")
    print(f"Early shortcuts:          {'enabled' if use_early_shortcuts else 'disabled'}")
    print(f"Correct zero-SP parses: {totals['correct_zero_sp']}")
    print(f"Zero-SP abstentions:    {totals['abstain_zero_sp']}")
    print(f"False-positive SPs:     {totals['false_positive_sp']}")

    print("\nBy control type and species category")
    print(f"{'Control':<22} {'Category':<18} {'Total':>6} {'Zero-SP':>10} {'Abstain':>10} {'FP SP':>10}")
    print("-" * 86)
    for key in sorted(by_kind_cat):
        counts = by_kind_cat[key]
        total = counts["total"]
        print(
            f"{key[0]:<22} {key[1]:<18} {total:>6} "
            f"{counts['correct_zero_sp']:>10} "
            f"{counts['abstain_zero_sp']:>10} "
            f"{counts['false_positive_sp']:>10}"
        )


if __name__ == "__main__":
    options = _parse_cli_args(sys.argv[1:])
    main(use_early_shortcuts=options["use_early_shortcuts"])
