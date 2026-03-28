#!/usr/bin/env python3
"""Validate the Trp41 domain-fold classifier across all ground truth sequences.

Runs classify_domain_fold on every Cys pair found in the ground truth dataset
and checks whether:
  - Ig/C-like domain pairs consistently have Trp in [C1+10, C1+22]
  - Groove/G-domain pairs consistently lack Trp in this window
  - The two signals are separable (low false-positive / false-negative)

This is the gating validation for Phase 1 of the grammar-guided decomposition
strategy.  The Trp41 marker must be >90% sensitive for C-like and >95% specific
(absent in G-domains) to proceed with grammar enumeration.

Usage:
    python scripts/validate_trp41_signal.py
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GT_CSV = ROOT / "data" / "sp_ground_truth.csv"


def main():
    from mhcseqs.domain_parsing import (
        classify_cys_pair,
        classify_domain_fold,
        find_cys_pairs,
    )

    if not GT_CSV.exists():
        print(f"Ground truth not found: {GT_CSV}", file=sys.stderr)
        sys.exit(1)

    with open(GT_CSV, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    print(f"Loaded {len(rows)} sequences from {GT_CSV.name}")
    print()

    # --- Collect all Cys pair classifications --------------------------------
    total_pairs = 0
    total_seqs = 0
    seqs_with_pairs = 0

    # Cross-tabulate: (existing_domain_type, new_fold_type) -> count
    cross_tab: Counter[tuple[str, str]] = Counter()
    # Trp41 presence by existing classification
    trp_by_existing: Counter[tuple[str, bool]] = Counter()
    # Confidence distribution
    confidence_by_fold: dict[str, list[float]] = {
        "c_like": [], "g_domain": [], "ambiguous": []
    }
    # Separation distribution
    sep_by_fold: dict[str, list[int]] = {
        "c_like": [], "g_domain": [], "ambiguous": []
    }
    # Track Trp offset distribution for c_like pairs
    trp_offsets: list[int] = []
    # Track cases where classifiers disagree
    disagreements: list[dict] = []

    for row in rows:
        seq = row["sequence"]
        total_seqs += 1
        pairs = find_cys_pairs(seq)
        if not pairs:
            continue
        seqs_with_pairs += 1

        for c1, c2, sep in pairs:
            total_pairs += 1

            # Existing classifier
            old_ann = classify_cys_pair(seq, c1, c2)
            old_type = old_ann.domain_type  # "groove" | "ig" | "ambiguous"

            # New Trp41 classifier
            new_ann = classify_domain_fold(seq, c1, c2)
            new_type = new_ann.fold_type  # "c_like" | "g_domain" | "ambiguous"

            cross_tab[(old_type, new_type)] += 1
            has_trp = new_ann.trp_position >= 0
            trp_by_existing[(old_type, has_trp)] += 1
            confidence_by_fold[new_type].append(new_ann.confidence)
            sep_by_fold[new_type].append(sep)

            if has_trp:
                trp_offsets.append(new_ann.trp_position - c1)

            # Track disagreements (excluding ambiguous)
            old_binary = {"groove": "g_domain", "ig": "c_like"}.get(old_type)
            if old_binary and old_binary != new_type and new_type != "ambiguous":
                disagreements.append({
                    "accession": row.get("accession", "?"),
                    "organism": row.get("organism", "?")[:30],
                    "c1": c1, "c2": c2, "sep": sep,
                    "old": old_type, "new": new_type,
                    "trp_pos": new_ann.trp_position,
                    "old_groove": f"{old_ann.groove_score:.1f}",
                    "old_ig": f"{old_ann.ig_score:.1f}",
                    "new_conf": f"{new_ann.confidence:.2f}",
                })

    # --- Report --------------------------------------------------------------
    print(f"Sequences: {total_seqs} total, {seqs_with_pairs} with Cys pairs")
    print(f"Total Cys pairs: {total_pairs}")
    print()

    # Cross-tabulation
    print("=== Cross-tabulation: existing vs Trp41 classifier ===")
    print(f"{'Existing':<12} {'c_like':>8} {'g_domain':>10} {'ambiguous':>10} {'total':>8}")
    print("-" * 52)
    for old_type in ("groove", "ig", "ambiguous"):
        row_counts = {
            nt: cross_tab[(old_type, nt)]
            for nt in ("c_like", "g_domain", "ambiguous")
        }
        row_total = sum(row_counts.values())
        print(
            f"{old_type:<12} {row_counts['c_like']:>8} "
            f"{row_counts['g_domain']:>10} {row_counts['ambiguous']:>10} "
            f"{row_total:>8}"
        )
    col_totals = {
        nt: sum(cross_tab[(ot, nt)] for ot in ("groove", "ig", "ambiguous"))
        for nt in ("c_like", "g_domain", "ambiguous")
    }
    total_all = sum(col_totals.values())
    print("-" * 52)
    print(
        f"{'total':<12} {col_totals['c_like']:>8} "
        f"{col_totals['g_domain']:>10} {col_totals['ambiguous']:>10} "
        f"{total_all:>8}"
    )
    print()

    # Trp presence by existing classification
    print("=== Trp41 presence by existing classification ===")
    for old_type in ("groove", "ig", "ambiguous"):
        has = trp_by_existing[(old_type, True)]
        no = trp_by_existing[(old_type, False)]
        total = has + no
        pct = 100 * has / total if total else 0
        print(f"  {old_type:<12}: {has}/{total} have Trp ({pct:.1f}%)")
    print()

    # Key metrics
    ig_with_trp = trp_by_existing[("ig", True)]
    ig_total = trp_by_existing[("ig", True)] + trp_by_existing[("ig", False)]
    groove_with_trp = trp_by_existing[("groove", True)]
    groove_total = trp_by_existing[("groove", True)] + trp_by_existing[("groove", False)]

    sensitivity = 100 * ig_with_trp / ig_total if ig_total else 0
    specificity = 100 * (groove_total - groove_with_trp) / groove_total if groove_total else 0

    print("=== Key metrics ===")
    print(f"  Trp41 sensitivity for Ig/C-like: {ig_with_trp}/{ig_total} = {sensitivity:.1f}%")
    no_trp = groove_total - groove_with_trp
    print(f"  Trp41 specificity (absent in groove): {no_trp}/{groove_total} = {specificity:.1f}%")
    print()

    # Trp offset distribution
    if trp_offsets:
        from collections import Counter as C
        offset_counts = C(trp_offsets)
        print("=== Trp offset distribution (offset from C1) ===")
        for off in sorted(offset_counts):
            bar = "#" * min(offset_counts[off], 80)
            print(f"  c1+{off:>2}: {offset_counts[off]:>5} {bar}")
        print()

    # Separation by fold type
    print("=== Separation distribution by fold type ===")
    for ft in ("c_like", "g_domain", "ambiguous"):
        seps = sep_by_fold[ft]
        if seps:
            avg = sum(seps) / len(seps)
            mn, mx = min(seps), max(seps)
            print(f"  {ft:<12}: n={len(seps):>5}, mean={avg:.1f}, range=[{mn}, {mx}]")
    print()

    # Disagreements
    if disagreements:
        print(f"=== Classifier disagreements ({len(disagreements)}) ===")
        hdr = (f"{'Accession':<14} {'Organism':<32} {'c1':>4} {'c2':>4} "
               f"{'sep':>4} {'old':<10} {'new':<10} {'trp':>4}")
        print(hdr)
        for d in disagreements[:30]:
            print(
                f"{d['accession']:<14} {d['organism']:<32} {d['c1']:>4} {d['c2']:>4} "
                f"{d['sep']:>4} {d['old']:<10} {d['new']:<10} {d['trp_pos']:>4} "
                f"{d['old_groove']:>5} {d['old_ig']:>5}"
            )
        if len(disagreements) > 30:
            print(f"  ... and {len(disagreements) - 30} more")
    else:
        print("No disagreements between classifiers!")
    print()

    # Gate check
    print("=== Gate check ===")
    gate_pass = sensitivity >= 90 and specificity >= 95
    status = "PASS" if gate_pass else "FAIL"
    print(f"  Sensitivity >= 90%: {sensitivity:.1f}% {'OK' if sensitivity >= 90 else 'FAIL'}")
    print(f"  Specificity >= 95%: {specificity:.1f}% {'OK' if specificity >= 95 else 'FAIL'}")
    print(f"  Overall: {status}")


if __name__ == "__main__":
    main()
