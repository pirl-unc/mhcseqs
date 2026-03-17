#!/usr/bin/env python3
"""Update the data summary table in README.md from the actual curated data.

Reads mhcseqs/diverse_mhc_sequences.csv and merges counts with the built
IMGT/IPD-MHC counts (from a build, or hardcoded baseline if no build exists).

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

# Source group → species_category mapping (mirrors pipeline.py)
_GROUP_TO_CATEGORY = {
    "reptile_lepidosauria": "other_vertebrate",
    "reptile_crocodylia": "other_vertebrate",
    "reptile_testudines": "other_vertebrate",
    "amphibian": "other_vertebrate",
    "bird_non_chicken": "bird",
    "chicken": "bird",
    "shark_ray": "fish",
    "bony_fish": "fish",
    "marsupial": "other_mammal",
    "monotreme": "other_mammal",
    "bat": "other_mammal",
}

# Baseline IMGT/IPD-MHC counts (from last full build)
_BASELINE = {
    ("human", "I"): 17462, ("human", "II"): 7878,
    ("nhp", "I"): 4639, ("nhp", "II"): 2486,
    ("murine", "I"): 59, ("murine", "II"): 29,
    ("ungulate", "I"): 638, ("ungulate", "II"): 1128,
    ("carnivore", "I"): 166, ("carnivore", "II"): 318,
    ("other_mammal", "I"): 3, ("other_mammal", "II"): 98,
    ("bird", "I"): 28, ("bird", "II"): 0,
    ("fish", "I"): 90, ("fish", "II"): 85,
    ("other_vertebrate", "I"): 0, ("other_vertebrate", "II"): 0,
}

CATEGORIES = [
    "human", "nhp", "murine", "ungulate", "carnivore",
    "other_mammal", "bird", "fish", "other_vertebrate",
]

# Try to read counts from a real build first
BUILT_CSV = Path.home() / ".cache" / "mhcseqs" / "mhc-full-seqs.csv"


def load_built_counts() -> Counter:
    """Load counts from built CSV if available, else use baseline."""
    if BUILT_CSV.exists():
        counts: Counter = Counter()
        with open(BUILT_CSV, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                cat = row.get("species_category", "")
                mc = row.get("mhc_class", "")
                src = row.get("source", "")
                # Exclude diverse entries (counted separately from the shipped CSV)
                if src == "uniprot_diverse":
                    continue
                if cat and mc in ("I", "II"):
                    counts[(cat, mc)] += 1
        if counts:
            return counts
    return Counter(_BASELINE)


def load_diverse_counts() -> tuple[Counter, int]:
    """Load counts from diverse MHC CSV. Returns (counts, num_prefixes)."""
    counts: Counter = Counter()
    prefixes = set()
    with open(DIVERSE_CSV, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            cat = _GROUP_TO_CATEGORY.get(row.get("source_group", ""), "")
            mc = row.get("mhc_class", "")
            gene = row.get("gene", "")
            if cat and mc in ("I", "II"):
                counts[(cat, mc)] += 1
            if "-" in gene:
                prefixes.add(gene.split("-")[0])
    return counts, len(prefixes)


def build_table(merged: Counter) -> str:
    """Build the markdown table string."""
    lines = []
    lines.append("| Category | Class I | Class II | Total |")
    lines.append("|---|---:|---:|---:|")
    total_i = total_ii = 0
    for cat in CATEGORIES:
        ci = merged.get((cat, "I"), 0)
        cii = merged.get((cat, "II"), 0)
        total_i += ci
        total_ii += cii
        lines.append(f"| {cat} | {ci:,} | {cii:,} | {ci + cii:,} |")
    lines.append(f"| **total** | **{total_i:,}** | **{total_ii:,}** | **{total_i + total_ii:,}** |")
    return "\n".join(lines)


def main():
    built = load_built_counts()
    diverse, num_prefixes = load_diverse_counts()
    merged = built + diverse
    total = sum(merged.values())
    diverse_total = sum(diverse.values())

    table = build_table(merged)

    # Read README
    text = README.read_text(encoding="utf-8")

    # Replace the summary section
    header_line = (
        f"All sources (IMGT/HLA, IPD-MHC, UniProt curated references, and {diverse_total:,}\n"
        f"diverse MHC sequences from UniProt) are merged into a single dataset:"
    )
    summary_line = (
        f"Covering {num_prefixes}+ species prefixes. Groove parse success rate on IMGT/IPD-MHC\n"
        f"entries: 99.6%."
    )

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
    print(f"Updated README.md with {total:,} total entries ({diverse_total:,} diverse, {num_prefixes} prefixes)")


if __name__ == "__main__":
    main()
