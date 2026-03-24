#!/usr/bin/env python3
"""Evaluate signal peptide heuristic accuracy against UniProt ground truth.

Loads data/sp_ground_truth.csv (UniProt Signal annotations), runs our
Cys-pair + refinement heuristic on each sequence, and compares the
predicted SP length against the annotated SP length.

Usage:
    python scripts/evaluate_sp_ground_truth.py
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GT_CSV = ROOT / "data" / "sp_ground_truth.csv"

# Taxon IDs for broad clade mapping (used when normalize_mhc_species misses)
_MAMMAL_TAXIDS = {9606, 10090, 10116}  # human, mouse, rat


def _species_category(organism: str, taxon_id: str = "") -> str:
    """Map organism name to mhcseqs species category.

    Uses the species module first, then falls back to taxonomy-based
    heuristics for organisms not in the mhcseqs species table.
    """
    from mhcseqs.species import normalize_mhc_species

    cat = normalize_mhc_species(organism)
    if cat:
        return cat

    org_lower = organism.lower()

    # Broad taxonomic heuristics for vertebrate clades in the ground truth

    # Aves (birds) — common family/genus indicators
    _BIRD_HINTS = (
        "gallus", "chicken", "meleagris", "turkey", "anas ", "duck", "goose",
        "anser", "columba", "pigeon", "strix", "owl", "falco", "eagle",
        "hawk", "pipra", "parus", "corvus", "fringill", "passer", "taeniopygia",
        "zebra finch", "ficedula", "cyanistes", "sturnus", "acrocephalus",
        "phylloscopus", "lanius", "hirundo", "apus", "coturnix", "quail",
        "numida", "struthio", "dromaius", "apteryx", "rhea", "accipiter",
        "buteo", "aquila", "haliaeetus", "catharus", "turdus", "serinus",
        "lonchura", "melopsittacus", "psittac", "cacatua", "ara ",
        "amazona", "nymphicus", "pelecanus", "phalacrocorax", "ardea",
        "ciconia", "phoenicopterus", "spheniscus", "pygoscelis", "aptenodytes",
        "larus", "sterna", "alcedo", "merops", "upupa", "bucorvus",
        "picus", "dendrocopos", "indicator", "galbula", "pterocles",
    )
    if any(h in org_lower for h in _BIRD_HINTS):
        return "bird"

    # Actinopterygii + Chondrichthyes (fish)
    _FISH_HINTS = (
        "danio", "salmo", "oncorhynch", "ictalur", "oreochrom", "cyprinus",
        "carassius", "brachydanio", "poecili", "oryzias", "xiphophorus",
        "gasterosteus", "takifugu", "tetraodon", "gadus", "dicentrarchus",
        "sparus", "pagrus", "lates", "seriola", "thunnus", "rachycentron",
        "hippoglossus", "paralichthys", "scophthalmus", "solea", "mugil",
        "channa", "clarias", "pangasius", "labeo", "catla", "tor ",
        "epinephelus", "lutjanus", "acanthopagrus", "salarias", "syngnathus",
        "hippocampus", "acipenser", "polyodon", "lepisosteus", "amia ",
        "polypterus", "latimeria", "protopterus", "lepidosiren",
        "squalus", "mustelus", "triakis", "heterodontus", "carcharhinus",
        "negaprion", "ginglymostoma", "chiloscyllium", "scyliorhinus",
        "raja ", "leucoraja", "dasyatis", "rhinobatos", "callorhinchus",
        "liparis", "larimichthys", "siniperca",
    )
    if any(h in org_lower for h in _FISH_HINTS):
        return "fish"

    # Reptilia
    _REPTILE_HINTS = (
        "python", "boa ", "elaphe", "pantherophis", "naja", "bungarus",
        "crotalus", "vipera", "bothrops", "agkistrodon", "notechis",
        "pseudonaja", "ophiophagus", "micrurus", "laticauda", "hydrophis",
        "anolis", "iguana", "pogona", "eublepharis", "gekko", "lacerta",
        "podarcis", "zootoca", "varanus", "heloderma", "tiliqua",
        "thamnophis", "lampropeltis", "coluber", "natrix", "salvator",
        "tupinambis", "chamaeleo", "furcifer", "brookesia", "sphenodon",
        "crocodylus", "alligator", "caiman", "gavialis", "tomistoma",
        "chelonia", "caretta", "dermochelys", "chrysemys", "trachemys",
        "testudo", "terrapene", "emys", "mauremys", "pelodiscus",
        "amblyrhynch", "cyclura", "conolophus", "ctenosaura",
    )
    if any(h in org_lower for h in _REPTILE_HINTS):
        return "other_vertebrate"

    # Amphibia
    _AMPHIBIAN_HINTS = (
        "xenopus", "rana ", "lithobates", "bufo", "rhinella", "bombina",
        "hyla ", "pseudacris", "eleutherodactylus", "dendrobates",
        "phyllobates", "allobates", "epipedobates", "ranitomeya",
        "ambystoma", "notophthalmus", "plethodon", "cynops", "triturus",
        "ichthyophis", "typhlonectes", "caecilia", "microcaecilia",
        "andrias", "cryptobranchus", "salamandra", "lissotriton",
        "nanorana", "quasipaa", "pipa", "hymenochirus",
        "leptobrachium", "megophrys", "scaphiopus", "pelophylax",
        "pyxicephalus", "conraua", "ceratobatrachus", "atelopus",
        "mantella", "boophis",
    )
    if any(h in org_lower for h in _AMPHIBIAN_HINTS):
        return "other_vertebrate"

    # Ground truth is exclusively vertebrate (fetched from Mammalia, Aves,
    # Actinopterygii, Reptilia, Amphibia, Chondrichthyes).  Anything still
    # unrecognized is a non-mammalian vertebrate we missed above.
    return "other_vertebrate"


def _try_parse(seq: str) -> tuple[int, str]:
    """Try all domain parsers and pick the most plausible result.

    When multiple parsers succeed, pick the one whose mature_start is
    in the plausible SP range [10, 50] and closest to a typical MHC
    signal peptide length (~23 aa).  This avoids the class-I parser
    falsely matching Cys pairs deep inside class-II proteins.

    Returns (mature_start, parser_used) or (0, "") on failure.
    """
    from mhcseqs.groove import (
        decompose_class_i,
        decompose_class_ii_alpha,
        decompose_class_ii_beta,
    )

    TYPICAL_SP = 23  # median MHC SP length across vertebrates
    candidates = []
    for parser, name in [
        (decompose_class_i, "class_I"),
        (decompose_class_ii_beta, "class_II_beta"),
        (decompose_class_ii_alpha, "class_II_alpha"),
    ]:
        try:
            result = parser(seq)
            if result.ok and result.mature_start > 0:
                candidates.append((result.mature_start, name))
        except Exception:
            pass

    if not candidates:
        return 0, ""

    # Prefer candidates in plausible SP range, closest to typical length
    in_range = [(ms, n) for ms, n in candidates if 10 <= ms <= 50]
    if in_range:
        return min(in_range, key=lambda x: abs(x[0] - TYPICAL_SP))

    # Fallback: pick the smallest positive mature_start
    return min(candidates, key=lambda x: x[0])


def evaluate():
    from mhcseqs.groove import refine_signal_peptide

    if not GT_CSV.exists():
        print(f"Ground truth not found: {GT_CSV}", file=sys.stderr)
        print("Run: python scripts/fetch_sp_ground_truth.py", file=sys.stderr)
        sys.exit(1)

    with open(GT_CSV, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    print(f"Loaded {len(rows)} ground truth entries from {GT_CSV.name}")
    print()

    # Counters
    total = 0
    parsed = 0
    unparsed = 0
    exact = 0
    within_1 = 0
    within_2 = 0
    within_3 = 0
    deltas = []
    by_category: dict[str, Counter] = {}
    by_reviewed: dict[str, Counter] = {}
    mismatches_gt3: list[dict] = []

    for row in rows:
        seq = row["sequence"]
        gt_sp_len = int(row["sp_length"])
        organism = row["organism"]
        reviewed = row["reviewed"]
        accession = row["accession"]
        total += 1

        # Parse sequence
        cys_start, parser = _try_parse(seq)
        if cys_start == 0:
            unparsed += 1
            continue

        parsed += 1

        # Determine species category and refine
        cat = _species_category(organism, row.get("taxon_id", ""))
        refined = refine_signal_peptide(seq, cys_start, cat)

        delta = refined - gt_sp_len
        deltas.append(delta)

        # Accuracy bins
        if delta == 0:
            exact += 1
        if abs(delta) <= 1:
            within_1 += 1
        if abs(delta) <= 2:
            within_2 += 1
        if abs(delta) <= 3:
            within_3 += 1

        if abs(delta) > 3:
            mismatches_gt3.append({
                "accession": accession,
                "organism": organism,
                "category": cat or "unknown",
                "gt_sp": gt_sp_len,
                "cys_pred": cys_start,
                "refined": refined,
                "delta": delta,
                "reviewed": reviewed,
                "parser": parser,
            })

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
    print(f"Parsed:           {parsed} ({100*parsed/total:.1f}%)")
    print(f"Unparsed:         {unparsed} ({100*unparsed/total:.1f}%)")
    print()
    if parsed > 0:
        print(f"Exact match:      {exact}/{parsed} ({100*exact/parsed:.1f}%)")
        print(f"Within +/-1 aa:   {within_1}/{parsed} ({100*within_1/parsed:.1f}%)")
        print(f"Within +/-2 aa:   {within_2}/{parsed} ({100*within_2/parsed:.1f}%)")
        print(f"Within +/-3 aa:   {within_3}/{parsed} ({100*within_3/parsed:.1f}%)")

    if deltas:
        import statistics
        print(f"\nMean delta:       {statistics.mean(deltas):+.2f} aa")
        print(f"Median delta:     {statistics.median(deltas):+.1f} aa")
        print(f"Std dev:          {statistics.stdev(deltas):.2f} aa")

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
            f"{c['exact']:>4} ({100*c['exact']/t:5.1f}%) "
            f"{c['within_1']:>4} ({100*c['within_1']/t:5.1f}%) "
            f"{c['within_2']:>4} ({100*c['within_2']/t:5.1f}%)"
        )

    # Per-reviewed breakdown
    print(f"\n{'Review status':<20} {'Total':>6} {'Exact':>8} {'<=1':>8}")
    print("-" * 48)
    for rev in sorted(by_reviewed.keys()):
        c = by_reviewed[rev]
        t = c["total"]
        print(
            f"{rev:<20} {t:>6} "
            f"{c['exact']:>4} ({100*c['exact']/t:5.1f}%) "
            f"{c['within_1']:>4} ({100*c['within_1']/t:5.1f}%)"
        )

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
    evaluate()
