#!/usr/bin/env python3
"""Benchmark mhcgnomes parsing against the mhcseqs diverse dataset.

Tracks parse rates across mhcgnomes versions and generates a comparison
report. Results are appended to data/mhcgnomes_benchmark.csv for
historical tracking.

Usage:
    python scripts/benchmark_mhcgnomes.py           # run and append
    python scripts/benchmark_mhcgnomes.py --plot     # run, append, and plot
    python scripts/benchmark_mhcgnomes.py --plot-only # just plot from existing data
"""

from __future__ import annotations

import argparse
import csv
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DIVERSE_CSV = ROOT / "mhcseqs" / "diverse_mhc_sequences.csv"
BENCHMARK_CSV = ROOT / "data" / "mhcgnomes_benchmark.csv"
PLOT_PATH = ROOT / "data" / "mhcgnomes_parse_rate.svg"

BENCHMARK_FIELDS = [
    "date",
    "mhcgnomes_version",
    "mhcseqs_version",
    "total_gene_organism_pairs",
    "parsed_as_is",
    "parsed_with_default_species",
    "wrong_species",
    "failed_known_species",
    "failed_unknown_species",
    "failed_doubled_prefix",
    "failed_other",
]


def extract_latin_binomial(organism: str) -> str:
    binomial = organism.split("(")[0].strip()
    parts = binomial.split()
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return binomial


def run_benchmark() -> dict:
    import mhcgnomes

    from mhcseqs.version import __version__ as mhcseqs_version

    mhcgnomes_version = mhcgnomes.__version__

    seen: set[tuple[str, str]] = set()
    parsed_as_is = 0
    parsed_default_sp = 0
    wrong_species = 0
    failed_known_sp = 0
    failed_unknown_sp = 0
    failed_doubled = 0
    failed_other = 0

    with open(DIVERSE_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            gene = row.get("gene", "")
            organism = row.get("organism", "")
            if not gene or (gene, organism) in seen:
                continue
            seen.add((gene, organism))

            # Strategy 1: parse as-is
            try:
                r = mhcgnomes.parse(gene)
                tp = type(r).__name__
                if tp in ("Gene", "Allele", "AlleleWithoutGene"):
                    sp = getattr(getattr(r, "species", None), "name", "")
                    if sp and sp.lower() in organism.lower():
                        parsed_as_is += 1
                        continue
                    elif sp:
                        wrong_species += 1
                        continue
                    else:
                        parsed_as_is += 1
                        continue
            except Exception:
                pass

            # Strategy 2: strip prefix, use default_species
            if "-" in gene:
                prefix = gene.split("-")[0]
                bare = gene.split("-", 1)[1]
                if bare.lower().startswith(prefix.lower()):
                    bare = bare[len(prefix) :].lstrip("-_")

                latin = extract_latin_binomial(organism)
                if latin:
                    try:
                        r2 = mhcgnomes.parse(bare, default_species=latin)
                        tp2 = type(r2).__name__
                        if tp2 in ("Gene", "Allele", "AlleleWithoutGene"):
                            sp2 = getattr(getattr(r2, "species", None), "name", "")
                            if sp2 and sp2.lower() in organism.lower():
                                parsed_default_sp += 1
                                continue
                    except Exception:
                        pass

                # Diagnose failure
                if bare.lower().startswith(prefix.lower()):
                    failed_doubled += 1
                    continue

                # Check if species known
                lp_parts = organism.split("(")[0].strip().split()
                lp = f"{lp_parts[0][:5]}{lp_parts[1][:5].capitalize()}" if len(lp_parts) >= 2 else None
                sp_known = False
                for p in [lp, prefix]:
                    if p:
                        try:
                            sp_obj = mhcgnomes.Species.get(p)
                            if sp_obj:
                                sp_known = True
                                break
                        except Exception:
                            pass

                if sp_known:
                    failed_known_sp += 1
                else:
                    failed_unknown_sp += 1
            else:
                failed_other += 1

    result = {
        "date": str(date.today()),
        "mhcgnomes_version": mhcgnomes_version,
        "mhcseqs_version": mhcseqs_version,
        "total_gene_organism_pairs": len(seen),
        "parsed_as_is": parsed_as_is,
        "parsed_with_default_species": parsed_default_sp,
        "wrong_species": wrong_species,
        "failed_known_species": failed_known_sp,
        "failed_unknown_species": failed_unknown_sp,
        "failed_doubled_prefix": failed_doubled,
        "failed_other": failed_other,
    }

    total = len(seen)
    ok = parsed_as_is + parsed_default_sp
    print(f"mhcgnomes {mhcgnomes_version} | mhcseqs {mhcseqs_version}")
    print(f"Total pairs:              {total:,}")
    print(f"Parsed (as-is):           {parsed_as_is:,} ({parsed_as_is / total * 100:.1f}%)")
    print(f"Parsed (default_species): {parsed_default_sp:,} ({parsed_default_sp / total * 100:.1f}%)")
    print(f"Total parseable:          {ok:,} ({ok / total * 100:.1f}%)")
    print(f"Wrong species:            {wrong_species:,}")
    print(f"Failed (known species):   {failed_known_sp:,}")
    print(f"Failed (unknown species): {failed_unknown_sp:,}")
    print(f"Failed (doubled prefix):  {failed_doubled:,}")
    print(f"Failed (other):           {failed_other:,}")

    return result


def append_result(result: dict):
    exists = BENCHMARK_CSV.exists()
    # Check if this version already has an entry
    if exists:
        with open(BENCHMARK_CSV, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row.get("mhcgnomes_version") == result["mhcgnomes_version"]:
                    print(f"\nVersion {result['mhcgnomes_version']} already in benchmark — updating.")
                    # Rewrite without the old entry
                    rows = []
                    f.seek(0)
                    reader = csv.DictReader(f)
                    for r in reader:
                        if r.get("mhcgnomes_version") != result["mhcgnomes_version"]:
                            rows.append(r)
                    rows.append(result)
                    with open(BENCHMARK_CSV, "w", encoding="utf-8", newline="") as wf:
                        writer = csv.DictWriter(wf, fieldnames=BENCHMARK_FIELDS)
                        writer.writeheader()
                        for r in rows:
                            writer.writerow(r)
                    return

    with open(BENCHMARK_CSV, "a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=BENCHMARK_FIELDS)
        if not exists:
            writer.writeheader()
        writer.writerow(result)
    print(f"\nAppended to {BENCHMARK_CSV}")


def plot_benchmark():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed — skipping plot. Install with: pip install matplotlib")
        return

    if not BENCHMARK_CSV.exists():
        print("No benchmark data to plot.")
        return

    rows = []
    with open(BENCHMARK_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)

    if len(rows) < 1:
        print("Need at least 1 data point to plot.")
        return

    versions = [r["mhcgnomes_version"] for r in rows]
    totals = [int(r["total_gene_organism_pairs"]) for r in rows]
    as_is = [int(r["parsed_as_is"]) for r in rows]
    with_ds = [int(r["parsed_with_default_species"]) for r in rows]
    wrong = [int(r["wrong_species"]) for r in rows]
    failed_known = [int(r["failed_known_species"]) for r in rows]
    failed_unknown = [int(r["failed_unknown_species"]) for r in rows]
    failed_doubled = [int(r["failed_doubled_prefix"]) for r in rows]
    failed_other = [int(r["failed_other"]) for r in rows]

    x = range(len(versions))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={"height_ratios": [2, 1]})

    # Top: stacked bar of parse outcomes
    ax1.bar(x, as_is, label="Parsed (as-is)", color="#2ecc71")
    ax1.bar(x, with_ds, bottom=as_is, label="Parsed (default_species)", color="#27ae60")
    ax1.bar(
        x,
        failed_known,
        bottom=[a + d for a, d in zip(as_is, with_ds)],
        label="Known species, unknown gene",
        color="#f39c12",
    )
    ax1.bar(
        x,
        failed_unknown,
        bottom=[a + d + k for a, d, k in zip(as_is, with_ds, failed_known)],
        label="Unknown species",
        color="#e74c3c",
    )
    ax1.bar(
        x,
        [d + o for d, o in zip(failed_doubled, failed_other)],
        bottom=[a + d + k + u for a, d, k, u in zip(as_is, with_ds, failed_known, failed_unknown)],
        label="Data issues (ours)",
        color="#95a5a6",
    )
    ax1.bar(
        x,
        wrong,
        bottom=[
            a + d + k + u + dd + o
            for a, d, k, u, dd, o in zip(as_is, with_ds, failed_known, failed_unknown, failed_doubled, failed_other)
        ],
        label="Wrong species",
        color="#8e44ad",
    )

    ax1.set_ylabel("Gene+organism pairs")
    ax1.set_title("mhcgnomes Parse Coverage of mhcseqs Diverse Dataset")
    ax1.set_xticks(list(x))
    ax1.set_xticklabels([f"v{v}" for v in versions], rotation=45, ha="right")
    ax1.legend(loc="upper left", fontsize=8)

    # Bottom: parse rate line
    rates = [(a + d) / t * 100 for a, d, t in zip(as_is, with_ds, totals)]
    wrong_rates = [w / t * 100 for w, t in zip(wrong, totals)]

    ax2.plot(list(x), rates, "o-", color="#2ecc71", linewidth=2, markersize=8, label="Parse rate %")
    ax2.plot(list(x), wrong_rates, "s--", color="#8e44ad", linewidth=1, markersize=5, label="Wrong species %")
    for i, (rate, v) in enumerate(zip(rates, versions)):
        ax2.annotate(f"{rate:.1f}%", (i, rate), textcoords="offset points", xytext=(0, 10), ha="center", fontsize=9)

    ax2.set_ylabel("Rate (%)")
    ax2.set_xlabel("mhcgnomes version")
    ax2.set_xticks(list(x))
    ax2.set_xticklabels([f"v{v}" for v in versions], rotation=45, ha="right")
    ax2.set_ylim(0, max(rates) * 1.3)
    ax2.legend(loc="upper left", fontsize=8)

    plt.tight_layout()
    plt.savefig(PLOT_PATH, dpi=150, bbox_inches="tight")
    print(f"Plot saved to {PLOT_PATH}")


def main():
    parser = argparse.ArgumentParser(description="Benchmark mhcgnomes parsing")
    parser.add_argument("--plot", action="store_true", help="Generate plot after benchmarking")
    parser.add_argument("--plot-only", action="store_true", help="Only generate plot from existing data")
    args = parser.parse_args()

    if args.plot_only:
        plot_benchmark()
        return

    result = run_benchmark()
    append_result(result)

    if args.plot:
        plot_benchmark()


if __name__ == "__main__":
    main()
