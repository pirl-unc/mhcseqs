#!/usr/bin/env python3
"""Generate a data funnel plot showing how sequences flow through the pipeline.

Shows: raw sources → dedup → unique proteins → groove parse outcomes.

Usage:
    python scripts/plot_funnel.py
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

RAW_CSV = Path.home() / ".cache" / "mhcseqs" / "mhc-seqs-raw.csv"
FULL_CSV = Path.home() / ".cache" / "mhcseqs" / "mhc-full-seqs.csv"
PLOT_PATH = Path(__file__).resolve().parent.parent / "data" / "pipeline_funnel.svg"


def main():
    if not RAW_CSV.exists() or not FULL_CSV.exists():
        print("No built CSVs — run `mhcseqs build` first.")
        return

    # Stage 1: Raw sources
    source_counts = Counter()
    raw_total = 0
    with open(RAW_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            src = row.get("source", "")
            # Unify UniProt sources for display
            if src.startswith("uniprot"):
                src = "UniProt"
            elif src == "imgt":
                src = "IMGT/HLA"
            elif src == "ipd_mhc":
                src = "IPD-MHC"
            source_counts[src] += 1
            raw_total += 1

    # Stage 2: Full seqs (after merge/dedup)
    full_total = 0
    groove_counts = Counter()
    functional = 0
    with open(FULL_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            full_total += 1
            status = row.get("groove_status", "")
            groove_counts[status] += 1
            if row.get("is_functional", "") == "True":
                functional += 1

    # Compute stages
    excluded = raw_total - full_total
    groove_ok = sum(
        n
        for s, n in groove_counts.items()
        if s in ("ok", "inferred_from_alpha3", "alpha1_only", "alpha2_only", "beta1_only_fallback", "fragment_fallback")
    )
    groove_fail = full_total - groove_ok - groove_counts.get("not_applicable", 0)
    not_applicable = groove_counts.get("not_applicable", 0)

    print("Pipeline Funnel")
    print("=" * 50)
    print(f"\nStage 1: Raw input ({raw_total:,})")
    for src, n in source_counts.most_common():
        print(f"  {src:15s} {n:6,}")

    print(f"\nStage 2: After merge/dedup ({full_total:,})")
    print(f"  Excluded (null/short/conflict): {excluded:,}")

    print("\nStage 3: Groove extraction")
    print(f"  Groove OK:          {groove_ok:6,} ({groove_ok / full_total * 100:.1f}%)")
    print(f"  Not applicable:     {not_applicable:6,} ({not_applicable / full_total * 100:.1f}%)")
    print(f"  Groove failed:      {groove_fail:6,} ({groove_fail / full_total * 100:.1f}%)")
    print(f"  Functional:         {functional:6,} ({functional / full_total * 100:.1f}%)")

    print("\nGroove status detail:")
    for status, n in groove_counts.most_common():
        print(f"  {status or '(empty)':25s} {n:6,}")

    # Plot
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not installed — skipping plot.")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    # Horizontal stacked bar for each stage
    stages = ["Raw input", "After merge", "Groove result"]
    y_pos = [2, 1, 0]

    # Stage 1: sources
    left = 0
    colors_src = {"IMGT/HLA": "#3498db", "IPD-MHC": "#2ecc71", "UniProt": "#e74c3c"}
    for src, n in source_counts.most_common():
        color = colors_src.get(src, "#95a5a6")
        ax.barh(2, n, left=left, color=color, edgecolor="white", label=src)
        if n > raw_total * 0.05:
            ax.text(left + n / 2, 2, f"{src}\n{n:,}", ha="center", va="center", fontsize=8, fontweight="bold")
        left += n

    # Stage 2: kept vs excluded
    ax.barh(1, full_total, color="#2ecc71", edgecolor="white")
    ax.barh(1, excluded, left=full_total, color="#e0e0e0", edgecolor="white")
    ax.text(full_total / 2, 1, f"Representatives\n{full_total:,}", ha="center", va="center", fontsize=8, fontweight="bold")
    ax.text(full_total + excluded / 2, 1, f"Excluded\n{excluded:,}", ha="center", va="center", fontsize=8, color="#888")

    # Stage 3: groove outcomes
    left = 0
    groove_colors = {
        "ok": "#27ae60",
        "alpha2_only": "#2ecc71",
        "alpha1_only": "#82e0aa",
        "beta1_only_fallback": "#58d68d",
        "fragment_fallback": "#abebc6",
        "inferred_from_alpha3": "#a9dfbf",
        "not_applicable": "#d5d8dc",
    }
    ok_statuses = [
        "ok",
        "alpha2_only",
        "beta1_only_fallback",
        "alpha1_only",
        "fragment_fallback",
        "inferred_from_alpha3",
    ]
    for status in ok_statuses:
        n = groove_counts.get(status, 0)
        if n == 0:
            continue
        color = groove_colors.get(status, "#82e0aa")
        ax.barh(0, n, left=left, color=color, edgecolor="white")
        if n > full_total * 0.03:
            ax.text(left + n / 2, 0, f"{status}\n{n:,}", ha="center", va="center", fontsize=6)
        left += n

    n_na = groove_counts.get("not_applicable", 0)
    if n_na:
        ax.barh(0, n_na, left=left, color="#d5d8dc", edgecolor="white")
        left += n_na

    fail_total = full_total - left
    if fail_total > 0:
        ax.barh(0, fail_total, left=left, color="#e74c3c", edgecolor="white")
        ax.text(left + fail_total / 2, 0, f"Failed\n{fail_total:,}", ha="center", va="center", fontsize=7, color="white")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(stages)
    ax.set_xlabel("Sequences")
    ax.set_title("mhcseqs Data Pipeline Funnel")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(PLOT_PATH, dpi=150, bbox_inches="tight")
    print(f"\nPlot saved to {PLOT_PATH}")


if __name__ == "__main__":
    main()
