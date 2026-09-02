#!/usr/bin/env python3
"""Update the data summary table in README.md from the actual curated data.

Reads the completed build, which already contains IMGT/HLA, IPD-MHC and all
curated supplemental sequences. The supplemental CSV is read only to report
its species coverage; its rows must not be added to the merged build again.

Usage:
    python scripts/update_readme_counts.py
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
README = ROOT / "README.md"
DIVERSE_CSV = ROOT / "mhcseqs" / "diverse_mhc_sequences.csv"

CATEGORIES = [
    "human",
    "nhp",
    "murine",
    "ungulate",
    "carnivore",
    "other_mammal",
    "bird",
    "fish",
    "other_vertebrate",
]

BUILT_CSV = Path.home() / ".cache" / "mhcseqs" / "mhc-full-seqs.csv"


def load_built_counts() -> Counter:
    """Load counts from built CSV. Raises if no build exists."""
    if not BUILT_CSV.exists():
        raise FileNotFoundError(f"No built CSV at {BUILT_CSV} — run `mhcseqs build` first.")
    counts: Counter = Counter()
    with open(BUILT_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            cat = row.get("species_category", "")
            mc = row.get("mhc_class", "")
            if cat and mc in ("I", "II"):
                counts[(cat, mc)] += 1
    if not counts:
        raise ValueError(f"Built CSV at {BUILT_CSV} has no valid entries.")
    return counts


# Groove statuses that mean no groove could be extracted at all. Every other
# status (ok, *_fallback, *_only, inferred_from_alpha3) yielded at least a
# partial groove and counts as a successful parse.
_GROOVE_FAILURE_STATUSES = {"", "missing_groove"}


def load_groove_success_rate() -> tuple[float, int, int]:
    """Compute groove parse success rate on IMGT/IPD-MHC entries.

    Returns (percent, n_success, n_total). Success = a (possibly partial) groove
    was extracted, i.e. groove_status is not in :data:`_GROOVE_FAILURE_STATUSES`.
    """
    total = 0
    success = 0
    with open(BUILT_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("source", "") not in ("imgt", "ipd_mhc"):
                continue
            total += 1
            if row.get("groove_status", "") not in _GROOVE_FAILURE_STATUSES:
                success += 1
    pct = 100.0 * success / total if total else 0.0
    return pct, success, total


def load_diverse_species_count() -> int:
    """Return the number of species represented by the supplemental CSV."""
    species = set()
    with open(DIVERSE_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            organism = row.get("organism", "")
            if organism:
                species.add(organism)
    return len(species)


def build_table(counts: Counter) -> str:
    """Build the markdown table string."""
    lines = []
    lines.append("| Category | Class I | Class II | Total |")
    lines.append("|---|---:|---:|---:|")
    total_i = total_ii = 0
    unexpected = sorted({category for category, _ in counts if category not in CATEGORIES})
    for cat in [*CATEGORIES, *unexpected]:
        ci = counts.get((cat, "I"), 0)
        cii = counts.get((cat, "II"), 0)
        total_i += ci
        total_ii += cii
        lines.append(f"| {cat} | {ci:,} | {cii:,} | {ci + cii:,} |")
    lines.append(f"| **total** | **{total_i:,}** | **{total_ii:,}** | **{total_i + total_ii:,}** |")
    return "\n".join(lines)


def main():
    built = load_built_counts()
    num_species = load_diverse_species_count()
    total = sum(built.values())

    table = build_table(built)

    # Read README
    text = README.read_text(encoding="utf-8")

    # Replace the summary section
    header_line = (
        "All sources (IMGT/HLA, IPD-MHC, and curated UniProt/GenBank references)\n"
        "are merged into a single dataset. Categorized class I/II representatives:"
    )
    groove_pct, _, _ = load_groove_success_rate()
    summary_line = f"Covering {num_species}+ species. Groove parse success rate on IMGT/IPD-MHC\nentries: {groove_pct:.1f}%."

    # Replace between "## Current data summary" and "## Structural decomposition"
    pattern = re.compile(
        r"(## Current data summary\n\n)"
        r".*?"
        r"(\n\n## Structural decomposition)",
        re.DOTALL,
    )
    replacement = f"\\1{header_line}\n\n{table}\n\n{summary_line}\\2"
    new_text = pattern.sub(replacement, text)

    if new_text == text:
        print("No changes needed.")
        return

    README.write_text(new_text, encoding="utf-8")
    print(f"Updated README.md with {total:,} categorized entries across {num_species}+ species")


if __name__ == "__main__":
    main()
