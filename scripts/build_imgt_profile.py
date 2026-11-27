#!/usr/bin/env python3
"""Build the multi-clade Class I IMGT G-DOMAIN structural PSSM + loop gap profile.

Samples clade-balanced Class I alpha-chain groove sequences spanning mammals
(human, NHP, murine, ungulate, carnivore, bat, marsupial, monotreme, other
mammals), birds, reptiles, amphibians, bony fish, and sharks/rays, and trains
per-domain log-odds PSSMs (L x 20 in half-bit units with BLOSUM62-mixture
pseudocounts) paired with position-specific gap-open penalties lowered at the
6 conserved structural loop/hinge indel hotspots of the IMGT G-DOMAIN fold.

Output: mhcseqs/imgt_class_i_profile.json (shipped as package data).

Usage:
    python scripts/build_imgt_profile.py [--full-csv PATH] [--output PATH]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

from mhcseqs.datafiles import FULL_SEQS_CSV, find_csv
from mhcseqs.imgt import (
    _BLOSUM62,
    G_ALPHA1_POSITIONS,
    G_ALPHA2_POSITIONS,
    G_ALPHA_MATURE_POSITIONS,
    G_ALPHA_POSITIONS,
    G_BETA_MATURE_POSITIONS,
    G_BETA_POSITIONS,
    HLA_A0201_G_ALPHA1,
    HLA_A0201_G_ALPHA2,
    HLA_DRA0101_G_ALPHA,
    HLA_DRB1_0101_G_BETA,
    g_alpha2_cys_anchor_bonus,
    g_beta_cys_anchor_bonus,
    imgt_to_mature_class_i,
)
from mhcseqs.species import get_species_clade

AA20 = "ARNDCQEGHILKMFPSTWYV"

# Robinson-Robinson amino acid background frequencies (normalized)
_RAW_BG = {
    "A": 0.0780,
    "R": 0.0513,
    "N": 0.0449,
    "D": 0.0536,
    "C": 0.0192,
    "Q": 0.0426,
    "E": 0.0630,
    "G": 0.0738,
    "H": 0.0220,
    "I": 0.0575,
    "L": 0.0902,
    "K": 0.0574,
    "M": 0.0224,
    "F": 0.0386,
    "P": 0.0520,
    "S": 0.0712,
    "T": 0.0584,
    "W": 0.0133,
    "Y": 0.0322,
    "V": 0.0644,
}
_BG_TOTAL = sum(_RAW_BG.values())
BG_FREQ = {aa: val / _BG_TOTAL for aa, val in _RAW_BG.items()}

# BLOSUM62 conditional substitution probabilities P_B62(b | a) for Dirichlet-style pseudocounts
B62_COND: dict[str, dict[str, float]] = {}
for _a in AA20:
    _raw = {b: BG_FREQ[b] * (2.0 ** (_BLOSUM62[_a][b] / 2.0)) for b in AA20}
    _tot = sum(_raw.values())
    B62_COND[_a] = {b: v / _tot for b, v in _raw.items()}

GAP_EXT = -1.0
END_GAP = -0.5


def make_position_gap_penalties(domain: str) -> tuple[list[float], list[float]]:
    """Return (gap_open_del, gap_open_ins) vectors encoding the 6 IMGT G-DOMAIN loop hotspots.

    Rigid beta-strands and core alpha-helices carry strong gap-open penalties (-14.0).
    The 6 structural loop/hinge hotspots carry lower penalties so indels land cleanly
    at biologically homologous sites without wobbling into adjacent strands/helices.
    """
    if domain == "G-ALPHA1":
        n_cols = len(HLA_A0201_G_ALPHA1)  # 90
        go_del = [-14.0] * n_cols
        go_ins = [-14.0] * n_cols

        # Hotspot 1: AB-turn (IMGT 14-19, 0-based 13..18)
        for idx in range(13, 19):
            go_del[idx] = -8.5
            go_ins[idx] = -8.0

        # Hotspot 2: CD-turn (IMGT 39-42, 0-based 38..41) — major avian/fish/reptile -2aa deletion site
        for idx in range(38, 42):
            go_del[idx] = -7.5
            go_ins[idx] = -8.5

        # Hotspot 3: N-helix hinge (IMGT 51-57, 0-based 50..56) — insertions flank
        # the 3_10-helix (after W51 for bat Ptal MDL, or after IMGT 54 for marsupial
        # Trvu VGE) rather than splitting the 52-54 turn internally.
        for idx in range(50, 57):
            go_ins[idx] = -11.5
            go_del[idx] = -10.5
        go_ins[50] = -6.5  # D-strand / 3_10-helix junction (after W51, e.g. Ptal-N*01:01 MDL)
        go_ins[53] = -6.5  # 3_10-helix / alpha1-helix hinge (after IMGT 54, e.g. Trvu-UB*01:01 VGE)
        go_ins[54] = -8.0
        go_ins[55] = -8.0

        # Mild flexibility at the alpha1/alpha2 linker boundary (IMGT 89-90)
        for idx in (88, 89):
            go_del[idx] = -9.5
            go_ins[idx] = -10.5
        return go_del, go_ins

    if domain == "G-ALPHA2":
        # 92-col IMGT non-gap backbone:
        #   0..38  -> IMGT 1..39 (C101 = idx 10)
        #   39..58 -> IMGT 42..61 (IMGT CD-turn gaps 40-41 sit right after idx 38)
        #   59     -> IMGT 61A
        #   60..70 -> IMGT 62..72
        #   71     -> IMGT 72A
        #   72..91 -> IMGT 73..92 (C164 = idx 73)
        n_cols = len(HLA_A0201_G_ALPHA2)  # 92
        go_del = [-14.0] * n_cols
        go_ins = [-14.0] * n_cols

        # Hotspot 4: AB-turn (IMGT 13-19, 0-based 12..18, right after Cys101)
        # Moderately lowered so a single -1/-2aa deletion is preferred over split deletions.
        for idx in range(12, 19):
            go_del[idx] = -9.5
            go_ins[idx] = -8.0

        # Hotspot 5: CD-turn (IMGT 38-44, 0-based 37..41; IMGT 40-41 gap sits after idx 38)
        for idx in range(37, 42):
            go_del[idx] = -8.5
            go_ins[idx] = -7.5
        go_ins[38] = -6.0  # Exact IMGT 40-41 CD-turn gap site

        # Hotspot 6: alpha2 helix kink (IMGT 58-65, 0-based 55..63, around 61A at idx 59)
        for idx in range(55, 64):
            go_ins[idx] = -8.0
            go_del[idx] = -10.0
        go_ins[59] = -7.0  # Exact IMGT 61A/61B helix kink site
        return go_del, go_ins

    if domain == "G-ALPHA":
        # 82-col IMGT G-ALPHA backbone (HLA-DRA*01:01 mature 1-82):
        #   0..1   -> 1.2, 1.1
        #   2..16  -> 1..15 (IMGT 7 = idx 8)
        #   17..36 -> 20..39 (IMGT 16-19 gap before idx 17; IMGT 40-41 gap after idx 36)
        #   37..51 -> 42..56 (W47 = idx 42; IMGT 54 = idx 49; IMGT 57-60 gap after idx 51)
        #   52..81 -> 61..90 (N86 = idx 77)
        n_cols = len(HLA_DRA0101_G_ALPHA)  # 82
        go_del = [-14.0] * n_cols
        go_ins = [-14.0] * n_cols

        # A-strand 7A insertion site (DQA1)
        go_ins[8] = -7.5
        go_ins[7] = -9.0
        go_ins[9] = -9.0

        # AB-turn (IMGT 14-21, idx 15..18)
        for idx in range(15, 19):
            go_del[idx] = -8.5
            go_ins[idx] = -8.0

        # CD-turn (IMGT 38-42, idx 35..37; IMGT 40-41 gap after idx 36)
        for idx in range(35, 38):
            go_del[idx] = -8.5
            go_ins[idx] = -8.0
        go_ins[36] = -6.5

        # N-helix hinge (IMGT 51-56, idx 46..51; IMGT 54A after idx 49)
        for idx in range(46, 52):
            go_del[idx] = -9.5
            go_ins[idx] = -8.5
        go_ins[49] = -6.5

        # Mild C-terminal linker flexibility (IMGT 89-90)
        for idx in (80, 81):
            go_del[idx] = -9.5
            go_ins[idx] = -10.0
        return go_del, go_ins

    if domain == "G-BETA":
        # 89-col IMGT G-BETA backbone (HLA-DRB1*01:01 mature 5-93):
        #   0..38  -> 1..39 (C11 = idx 10, N15 = idx 14)
        #   39..58 -> 42..61 (IMGT 40-41 gap after idx 38; W59 = idx 56)
        #   59..60 -> 61A, 61B
        #   61..71 -> 62..72
        #   72     -> 72A
        #   73..88 -> 73..88 (C74 = idx 74)
        n_cols = len(HLA_DRB1_0101_G_BETA)  # 89
        go_del = [-14.0] * n_cols
        go_ins = [-14.0] * n_cols

        # AB-turn (IMGT 13-21, idx 12..20 — captures DPB1's -2aa deletion cleanly)
        for idx in range(12, 21):
            go_del[idx] = -8.5
            go_ins[idx] = -8.0

        # CD-turn (IMGT 38-43, idx 37..40; IMGT 40-41 gap after idx 38)
        for idx in range(37, 41):
            go_del[idx] = -8.5
            go_ins[idx] = -7.5
        go_ins[38] = -6.0

        # Helix kink (IMGT 58-65, idx 55..64, around 61A/61B at idx 59..60)
        for idx in range(55, 65):
            go_ins[idx] = -8.5
            go_del[idx] = -10.0
        for idx in (59, 60, 61):
            go_del[idx] = -8.5
        go_ins[60] = -7.0
        return go_del, go_ins

    raise ValueError(f"Unsupported domain: {domain!r}")


def _profile_gotoh(
    pssm: list[dict[str, float]],
    go_del: list[float],
    go_ins: list[float],
    qry: str,
    *,
    anchor_bonus: dict[tuple[int, int], float] | None = None,
) -> tuple[float, str, str]:
    m = len(pssm)
    n = len(qry)
    neg = -1e9
    mat = [[neg] * (n + 1) for _ in range(m + 1)]
    ix = [[neg] * (n + 1) for _ in range(m + 1)]
    iy = [[neg] * (n + 1) for _ in range(m + 1)]
    tr_m = [[0] * (n + 1) for _ in range(m + 1)]
    tr_x = [[0] * (n + 1) for _ in range(m + 1)]
    tr_y = [[0] * (n + 1) for _ in range(m + 1)]

    mat[0][0] = 0.0
    for i in range(1, m + 1):
        ix[i][0] = END_GAP * i
        tr_x[i][0] = 1 if i > 1 else 0
    for j in range(1, n + 1):
        iy[0][j] = END_GAP * j
        tr_y[0][j] = 2 if j > 1 else 0

    for i in range(1, m + 1):
        col_scores = pssm[i - 1]
        del_open_i = go_del[i - 1]
        ins_open_i = go_ins[i - 1]
        m_i, m_p = mat[i], mat[i - 1]
        ix_i, ix_p = ix[i], ix[i - 1]
        iy_i, iy_p = iy[i], iy[i - 1]
        trm_i, trx_i, try_i = tr_m[i], tr_x[i], tr_y[i]
        is_last_ref = i == m

        for j in range(1, n + 1):
            qa = qry[j - 1]
            sub = col_scores.get(qa, -1.0)
            if anchor_bonus is not None:
                b = anchor_bonus.get((i - 1, j - 1))
                if b is not None:
                    sub += b

            v0, v1, v2 = m_p[j - 1], ix_p[j - 1], iy_p[j - 1]
            if v0 >= v1 and v0 >= v2:
                m_i[j] = v0 + sub
                trm_i[j] = 0
            elif v1 >= v2:
                m_i[j] = v1 + sub
                trm_i[j] = 1
            else:
                m_i[j] = v2 + sub
                trm_i[j] = 2

            go_x = END_GAP if j == n else del_open_i
            ge_x = END_GAP if j == n else GAP_EXT
            open_x = m_p[j] + go_x
            ext_x = ix_p[j] + ge_x
            if ext_x >= open_x:
                ix_i[j] = ext_x
                trx_i[j] = 1
            else:
                ix_i[j] = open_x
                trx_i[j] = 0

            go_y = END_GAP if is_last_ref else ins_open_i
            ge_y = END_GAP if is_last_ref else GAP_EXT
            open_y = m_i[j - 1] + go_y
            ext_y = iy_i[j - 1] + ge_y
            if ext_y >= open_y:
                iy_i[j] = ext_y
                try_i[j] = 2
            else:
                iy_i[j] = open_y
                try_i[j] = 0

    score, state = max(((mat[m][n], 0), (ix[m][n], 1), (iy[m][n], 2)))
    ar: list[str] = []
    aq: list[str] = []
    i, j = m, n
    while i > 0 or j > 0:
        if state == 0:
            prev = tr_m[i][j]
            ar.append("M")
            aq.append(qry[j - 1])
            i -= 1
            j -= 1
        elif state == 1:
            prev = tr_x[i][j]
            ar.append("M")
            aq.append("-")
            i -= 1
        else:
            prev = tr_y[i][j]
            ar.append("-")
            aq.append(qry[j - 1])
            j -= 1
        state = prev
    ar.reverse()
    aq.reverse()
    return score, "".join(ar), "".join(aq)


_DIVERSE_MHC_CSV = Path(__file__).resolve().parent.parent / "mhcseqs" / "diverse_mhc_sequences.csv"

_DIVERSE_SOURCE_GROUP_TO_CLADE: dict[str, str] = {
    "amphibian": "amphibian",
    "bat": "bat",
    "bird_non_chicken": "bird",
    "bony_fish": "bony_fish",
    "cat": "carnivore",
    "chicken": "bird",
    "marsupial": "marsupial",
    "monotreme": "monotreme",
    "murine": "murine",
    "other_mammal": "other_mammal",
    "reptile_crocodylia": "reptile",
    "reptile_lepidosauria": "reptile",
    "reptile_testudines": "reptile",
    "shark_ray": "shark_ray",
}

# Map mhcgnomes clades into the profile's vertebrate clade buckets
_MHCGNOMES_CLADE_TO_PROFILE_CLADE: dict[str, str] = {
    "human": "human",
    "nhp": "nhp",
    "murine": "murine",
    "ungulate": "ungulate",
    "carnivore": "carnivore",
    "bat": "bat",
    "marsupial": "marsupial",
    "monotreme": "monotreme",
    "cetacean": "other_mammal",
    "other_mammal": "other_mammal",
    "bird": "bird",
    "reptile": "reptile",
    "amphibian": "amphibian",
    "bony_fish": "bony_fish",
    "lobe_finned_fish": "bony_fish",
    "shark_ray": "shark_ray",
}


def load_curated_source_groups(diverse_csv: Path = _DIVERSE_MHC_CSV) -> dict[str, str]:
    """Load verified UniProt accession -> curated source_group mappings."""
    groups: dict[str, str] = {}
    if not diverse_csv.exists():
        return groups
    with open(diverse_csv, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            acc = (row.get("uniprot_accession") or "").strip()
            sg = (row.get("source_group") or "").strip()
            if acc and sg:
                groups[acc] = sg
    return groups


def classify_clade(row: dict, curated_source_groups: dict[str, str] | None = None) -> str | None:
    """Assign a verified vertebrate clade label for clade-balanced PSSM weighting.

    Uses authoritative ``mhcgnomes.Species.clade`` first (verified against NCBI
    Taxonomy), falling back to packaged curated ``source_group`` annotations by
    accession for organisms not yet in mhcgnomes.  Returns ``None`` for
    unverified organisms rather than guessing a default clade.
    """
    species = (row.get("species") or "").strip()
    mhcgnomes_clade = get_species_clade(species)
    if mhcgnomes_clade is not None:
        return _MHCGNOMES_CLADE_TO_PROFILE_CLADE.get(mhcgnomes_clade)

    source_id = (row.get("source_id") or "").strip()
    if curated_source_groups and source_id in curated_source_groups:
        sg = curated_source_groups[source_id]
        return _DIVERSE_SOURCE_GROUP_TO_CLADE.get(sg)

    cat = (row.get("species_category") or "").strip()
    if cat in ("human", "nhp", "murine", "ungulate", "carnivore", "bird", "other_mammal"):
        return cat
    # Uncurated fish/other_vertebrate rows are excluded unless verified via mhcgnomes or source_group.
    return None


def train_domain_profile(
    domain: str,
    clade_seqs: dict[str, list[str]],
    *,
    beta: float = 5.0,
    iterations: int = 2,
) -> dict:
    if domain == "G-ALPHA1":
        ref_seq = HLA_A0201_G_ALPHA1
        imgt_positions = list(G_ALPHA1_POSITIONS)
        mature_positions = [imgt_to_mature_class_i(domain, p) for p in imgt_positions]
    elif domain == "G-ALPHA2":
        ref_seq = HLA_A0201_G_ALPHA2
        imgt_positions = list(G_ALPHA2_POSITIONS)
        mature_positions = [imgt_to_mature_class_i(domain, p) for p in imgt_positions]
    elif domain == "G-ALPHA":
        ref_seq = HLA_DRA0101_G_ALPHA
        imgt_positions = list(G_ALPHA_POSITIONS)
        mature_positions = list(G_ALPHA_MATURE_POSITIONS)
    elif domain == "G-BETA":
        ref_seq = HLA_DRB1_0101_G_BETA
        imgt_positions = list(G_BETA_POSITIONS)
        mature_positions = list(G_BETA_MATURE_POSITIONS)
    else:
        raise ValueError(f"Unsupported domain: {domain!r}")

    go_del, go_ins = make_position_gap_penalties(domain)
    n_cols = len(ref_seq)

    # Seed pass 0 with BLOSUM62 around canonical human reference
    pssm = [{b: float(_BLOSUM62[ref_seq[i]][b]) for b in AA20} for i in range(n_cols)]

    for _ in range(iterations):
        col_counts: list[dict[str, float]] = [defaultdict(float) for _ in range(n_cols)]
        col_totals = [0.0] * n_cols
        for _clade, seqs in sorted(clade_seqs.items()):
            if not seqs:
                continue
            weight = 10.0 / len(seqs)  # Equal total weight per vertebrate clade
            for seq in seqs:
                if domain == "G-ALPHA2":
                    bonus = g_alpha2_cys_anchor_bonus(seq)
                elif domain == "G-BETA":
                    bonus = g_beta_cys_anchor_bonus(seq)
                else:
                    bonus = None
                _, ar, aq = _profile_gotoh(pssm, go_del, go_ins, seq, anchor_bonus=bonus)
                ref_i = 0
                for rch, qch in zip(ar, aq):
                    if rch != "-":
                        if qch in BG_FREQ:
                            col_counts[ref_i][qch] += weight
                            col_totals[ref_i] += weight
                        ref_i += 1

        new_pssm: list[dict[str, float]] = []
        for i in range(n_cols):
            tot = col_totals[i]
            obs = {a: (col_counts[i].get(a, 0.0) / tot if tot > 0 else BG_FREQ[a]) for a in AA20}
            pseudo = {b: sum(obs[a] * B62_COND[a][b] for a in AA20) for b in AA20}
            eff = min(tot, 40.0)
            col: dict[str, float] = {}
            for a in AA20:
                prob = max((eff * obs[a] + beta * pseudo[a]) / (eff + beta), 1e-6)
                col[a] = round(2.0 * math.log2(prob / BG_FREQ[a]), 2)
            new_pssm.append(col)
        pssm = new_pssm

    consensus = "".join(max(col, key=lambda aa: col[aa]) for col in pssm)
    pssm_rows = [[col[a] for a in AA20] for col in pssm]
    return {
        "length": n_cols,
        "imgt_positions": imgt_positions,
        "mature_positions": mature_positions,
        "consensus": consensus,
        "gap_open_del": go_del,
        "gap_open_ins": go_ins,
        "pssm": pssm_rows,
    }


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _sample_clades(
    rows: list[dict],
    curated_source_groups: dict[str, str],
    per_clade: int,
) -> tuple[dict[str, list[dict]], dict[str, int], dict[str, list[str]], int]:
    by_clade: dict[str, list[dict]] = defaultdict(list)
    skipped = 0
    for r in rows:
        clade = classify_clade(r, curated_source_groups=curated_source_groups)
        if clade is None:
            skipped += 1
            continue
        by_clade[clade].append(r)

    picked_by_clade: dict[str, list[dict]] = {}
    counts: dict[str, int] = {}
    rep_ids: dict[str, list[str]] = {}
    for clade, crows in sorted(by_clade.items()):
        crows.sort(key=lambda r: (r["species"], r["gene"], r["two_field_allele"], r.get("source_id", "")))
        step = max(1, len(crows) // per_clade)
        picked = crows[::step][:per_clade]
        picked_by_clade[clade] = picked
        counts[clade] = len(picked)
        rep_ids[clade] = [(r.get("source_id") or r["two_field_allele"]).strip() for r in picked]
        print(f"  clade {clade:16s}: {len(crows):6d} eligible -> {len(picked):2d} representatives")
    if skipped:
        print(f"  skipped (unverified clade): {skipped}")
    return picked_by_clade, counts, rep_ids, skipped


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mhc-class", choices=("I", "II"), default="I", help="MHC class profile to build (I or II)")
    parser.add_argument("--full-csv", type=Path, default=None, help="Path to built mhc-full-seqs.csv")
    parser.add_argument("--release", type=str, default="", help="Release tag of the input CSV (e.g. v2.9.0)")
    parser.add_argument("--output", type=Path, default=None, help="Output profile JSON path")
    parser.add_argument("--per-clade", type=int, default=20, help="Max representative sequences per clade")
    args = parser.parse_args()

    default_out = (
        Path(__file__).resolve().parent.parent / "mhcseqs" / ("imgt_class_i_profile.json" if args.mhc_class == "I" else "imgt_class_ii_profile.json")
    )
    output_path = args.output or default_out
    full_csv = args.full_csv or find_csv(FULL_SEQS_CSV)
    source_sha256 = _sha256_file(full_csv)
    curated_source_groups = load_curated_source_groups()

    with open(full_csv, "r", encoding="utf-8") as f:
        all_rows = list(csv.DictReader(f))

    if args.mhc_class == "I":
        rows = [
            r
            for r in all_rows
            if r["mhc_class"] == "I"
            and r["chain"] == "alpha"
            and r["groove_status"] == "ok"
            and 85 <= len(r["groove1"]) <= 95
            and 88 <= len(r["groove2"]) <= 96
            and "X" not in r["groove1"]
            and "X" not in r["groove2"]
        ]
        picked_by_clade, clade_counts, representative_ids, _ = _sample_clades(rows, curated_source_groups, args.per_clade)
        clade_g1 = {c: [r["groove1"] for r in rs] for c, rs in picked_by_clade.items()}
        clade_g2 = {c: [r["groove2"] for r in rs] for c, rs in picked_by_clade.items()}
        profile = {
            "version": "1",
            "alphabet": AA20,
            "gap_ext": GAP_EXT,
            "end_gap": END_GAP,
            "training_source": {
                "filename": full_csv.name,
                "release": args.release,
                "sha256": source_sha256,
                "per_clade": args.per_clade,
            },
            "clade_representatives": clade_counts,
            "representative_ids": representative_ids,
            "domains": {
                "G-ALPHA1": train_domain_profile("G-ALPHA1", clade_g1),
                "G-ALPHA2": train_domain_profile("G-ALPHA2", clade_g2),
            },
        }
    else:
        alpha_rows = [
            r
            for r in all_rows
            if r["mhc_class"] == "II"
            and r["chain"] == "alpha"
            and r["groove_status"] == "ok"
            and 76 <= len(r["groove1"]) <= 92
            and "X" not in r["groove1"]
        ]
        beta_rows = [
            r
            for r in all_rows
            if r["mhc_class"] == "II"
            and r["chain"] == "beta"
            and r["groove_status"] == "ok"
            and 84 <= len(r["groove2"]) <= 98
            and "X" not in r["groove2"]
        ]
        print("Sampling Class II G-ALPHA representatives:")
        alpha_picked, alpha_counts, alpha_rep_ids, _ = _sample_clades(alpha_rows, curated_source_groups, args.per_clade)
        print("Sampling Class II G-BETA representatives:")
        beta_picked, beta_counts, beta_rep_ids, _ = _sample_clades(beta_rows, curated_source_groups, args.per_clade)

        clade_ga = {c: [r["groove1"] for r in rs] for c, rs in alpha_picked.items()}
        clade_gb = {c: [r["groove2"] for r in rs] for c, rs in beta_picked.items()}
        profile = {
            "version": "1",
            "alphabet": AA20,
            "gap_ext": GAP_EXT,
            "end_gap": END_GAP,
            "training_source": {
                "filename": full_csv.name,
                "release": args.release,
                "sha256": source_sha256,
                "per_clade": args.per_clade,
            },
            "clade_representatives": {
                "G-ALPHA": alpha_counts,
                "G-BETA": beta_counts,
            },
            "representative_ids": {
                "G-ALPHA": alpha_rep_ids,
                "G-BETA": beta_rep_ids,
            },
            "domains": {
                "G-ALPHA": train_domain_profile("G-ALPHA", clade_ga),
                "G-BETA": train_domain_profile("G-BETA", clade_gb),
            },
        }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(profile, f, separators=(",", ":"))
        f.write("\n")

    size_kb = output_path.stat().st_size / 1024
    print(f"\nWrote {output_path} ({size_kb:.1f} KB)")
    for dname, dval in profile["domains"].items():
        print(f"  {dname:9s} consensus: {dval['consensus']}")


if __name__ == "__main__":
    main()
