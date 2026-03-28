#!/usr/bin/env python3
"""Evaluate signal peptide heuristic accuracy against UniProt ground truth.

Prefers ``data/sp_ground_truth_enriched.csv`` when available.  That enriched
benchmark carries gold MHC class / chain metadata so the evaluator can dispatch
the parser directly instead of guessing class from the sequence alone.  When
the enriched file is absent, the evaluator falls back to the original raw GT
CSV and the older classless parser-selection heuristic.

Usage:
    python scripts/evaluate_sp_ground_truth.py
"""
from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GT_RAW_CSV = ROOT / "data" / "sp_ground_truth.csv"
GT_ENRICHED_CSV = ROOT / "data" / "sp_ground_truth_enriched.csv"
NEGATIVE_CONTROL_CSV = ROOT / "data" / "sp_negative_controls.csv"
GT_CSV = GT_ENRICHED_CSV if GT_ENRICHED_CSV.exists() else GT_RAW_CSV

# Exact taxon IDs we can categorize without ambiguity.
_TAXID_TO_CATEGORY = {
    9606: "human",
    10090: "murine",
    10116: "murine",
}


_HINT_REGEX_CACHE: dict[str, re.Pattern[str]] = {}


def _parse_cli_args(argv: list[str]) -> dict[str, bool]:
    """Return simple evaluator options from argv."""
    return {
        "use_early_shortcuts": "--no-early-shortcuts" not in argv,
    }


def _hint_match(text: str, keyword: str) -> bool:
    """Word-boundary-aware fallback matching for genus/common-name hints."""
    if keyword.endswith((" ", "-")):
        return keyword in text
    pat = _HINT_REGEX_CACHE.get(keyword)
    if pat is None:
        pat = re.compile(r"\b" + re.escape(keyword) + r"\b")
        _HINT_REGEX_CACHE[keyword] = pat
    return bool(pat.search(text))


def _species_category(organism: str, taxon_id: str = "") -> str:
    """Map organism name to mhcseqs species category.

    Uses the species module first, then falls back to taxonomy-based
    heuristics for organisms not in the mhcseqs species table.
    """
    from mhcseqs.species import extract_latin_binomial, normalize_mhc_species

    try:
        taxid_value = int(str(taxon_id or "").strip())
    except ValueError:
        taxid_value = 0
    exact_taxid_cat = _TAXID_TO_CATEGORY.get(taxid_value)
    if exact_taxid_cat:
        return exact_taxid_cat

    for candidate in (organism, extract_latin_binomial(organism)):
        cat = normalize_mhc_species(candidate)
        if cat:
            return cat

    org_lower = extract_latin_binomial(organism).lower()

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
    if any(_hint_match(org_lower, h) for h in _FISH_HINTS):
        return "fish"
    if any(_hint_match(org_lower, h) for h in _BIRD_HINTS):
        return "bird"

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
    if any(_hint_match(org_lower, h) for h in _REPTILE_HINTS):
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
    if any(_hint_match(org_lower, h) for h in _AMPHIBIAN_HINTS):
        return "other_vertebrate"

    # Ground truth is exclusively vertebrate (fetched from Mammalia, Aves,
    # Actinopterygii, Reptilia, Amphibia, Chondrichthyes).  Anything still
    # unrecognized is a non-mammalian vertebrate we missed above.
    return "other_vertebrate"


def load_ground_truth_rows(prefer_enriched: bool = True) -> tuple[Path, list[dict[str, str]]]:
    """Load raw or enriched GT rows, preferring the enriched file when present."""
    path = GT_ENRICHED_CSV if prefer_enriched and GT_ENRICHED_CSV.exists() else GT_RAW_CSV
    with open(path, "r", encoding="utf-8") as handle:
        return path, list(csv.DictReader(handle))


def _row_species_category(row: dict[str, str]) -> str:
    return row.get("species_category", "") or _species_category(row["organism"], row.get("taxon_id", ""))


def _row_dispatch_metadata(row: dict[str, str]) -> tuple[str, str, str]:
    mhc_class = str(row.get("mhc_class", "") or "").strip().upper()
    chain = str(row.get("chain", "") or "").strip().lower()
    gene = str(row.get("gene", "") or "").strip()
    return mhc_class, chain, gene


def _parser_name_for_dispatch(mhc_class: str, chain: str) -> str:
    if mhc_class == "I":
        return "class_I"
    if mhc_class == "II" and chain == "alpha":
        return "class_II_alpha"
    if mhc_class == "II" and chain == "beta":
        return "class_II_beta"
    return ""


def _try_parse(seq: str, *, features=None, use_early_shortcuts: bool = True) -> tuple[int, str]:
    """Try all domain parsers and pick the most plausible result.

    When multiple parsers succeed, pick the one whose mature_start is
    in the plausible SP range [10, 50] and closest to a typical MHC
    signal peptide length (~23 aa).  This avoids the class-I parser
    falsely matching Cys pairs deep inside class-II proteins.

    Returns (mature_start, parser_used) or (0, "") on failure.
    """
    from mhcseqs.domain_parsing import (
        analyze_sequence,
        decompose_class_i,
        decompose_class_ii_alpha,
        decompose_class_ii_beta,
    )

    if features is None:
        features = analyze_sequence(seq)

    TYPICAL_SP = 23  # median MHC SP length across vertebrates
    candidates = []
    for parser, name in [
        (decompose_class_i, "class_I"),
        (decompose_class_ii_beta, "class_II_beta"),
        (decompose_class_ii_alpha, "class_II_alpha"),
    ]:
        try:
            result = parser(seq, features=features, use_early_shortcuts=use_early_shortcuts)
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


def predict_sp_for_row(
    row: dict[str, str],
    *,
    use_early_shortcuts: bool = True,
) -> dict[str, str | int | bool]:
    """Predict SP length for a GT row, using gold dispatch when available."""
    from mhcseqs.domain_parsing import analyze_sequence, decompose_domains, refine_signal_peptide

    seq = row["sequence"]
    category = _row_species_category(row)
    mhc_class, chain, gene = _row_dispatch_metadata(row)
    features = analyze_sequence(seq)

    if mhc_class in {"I", "II"}:
        parser_name = _parser_name_for_dispatch(mhc_class, chain)
        try:
            result = decompose_domains(
                seq,
                mhc_class=mhc_class,
                chain=chain or None,
                gene=gene,
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
            (int(result.anchor_cys1), int(result.anchor_cys2))
            if result.anchor_cys1 is not None and result.anchor_cys2 is not None
            else None
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
    cys_start, parser_name = _try_parse(seq, features=features, use_early_shortcuts=use_early_shortcuts)
    if cys_start == 0:
        return {
            "ok": False,
            "dispatch_mode": "inferred",
            "parser": "",
            "mhc_class": "",
            "chain": "",
            "gene": "",
            "status": "unparsed",
            "mature_start": 0,
            "predicted_sp": 0,
        }

    inferred_class = "I" if parser_name == "class_I" else "II"
    refined = refine_signal_peptide(seq, cys_start, category, inferred_class, features=features)
    return {
        "ok": True,
        "dispatch_mode": "inferred",
        "parser": parser_name,
        "mhc_class": inferred_class,
        "chain": (
            "alpha" if parser_name == "class_II_alpha"
            else ("beta" if parser_name == "class_II_beta" else "alpha")
        ),
        "gene": "",
        "status": "ok",
        "mature_start": int(cys_start),
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
    mismatches_gt3: list[dict] = []

    for row in rows:
        gt_sp_len = int(row["sp_length"])
        organism = row["organism"]
        reviewed = row["reviewed"]
        accession = row["accession"]
        total += 1

        prediction = predict_sp_for_row(row, use_early_shortcuts=use_early_shortcuts)
        if not prediction["ok"]:
            unparsed += 1
            continue

        parsed += 1
        by_dispatch[str(prediction["dispatch_mode"])] += 1

        cat = _row_species_category(row)
        refined = int(prediction["predicted_sp"])

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
                "cys_pred": int(prediction["mature_start"]),
                "refined": refined,
                "delta": delta,
                "reviewed": reviewed,
                "parser": str(prediction["parser"]),
                "mhc_class": str(prediction["mhc_class"]),
                "chain": str(prediction["chain"]),
                "dispatch_mode": str(prediction["dispatch_mode"]),
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
    print(f"Parsed:           {parsed} ({100*parsed/total:.1f}%)")
    print(f"Unparsed:         {unparsed} ({100*unparsed/total:.1f}%)")
    if by_dispatch:
        dispatch_summary = ", ".join(f"{k}={v}" for k, v in sorted(by_dispatch.items()))
        print(f"Dispatch mode:    {dispatch_summary}")
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
                f"{c['exact']:>4} ({100*c['exact']/t:5.1f}%) "
                f"{c['within_1']:>4} ({100*c['within_1']/t:5.1f}%) "
                f"{c['within_2']:>4} ({100*c['within_2']/t:5.1f}%) "
                f"{c['within_3']:>4} ({100*c['within_3']/t:5.1f}%)"
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
    options = _parse_cli_args(sys.argv[1:])
    evaluate(use_early_shortcuts=options["use_early_shortcuts"])
