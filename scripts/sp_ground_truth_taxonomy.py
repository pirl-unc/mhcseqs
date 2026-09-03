"""Taxonomy rules and audited metadata for the SP ground-truth corpus."""

from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TAXONOMY_AUDIT_CSV = ROOT / "data" / "sp_ground_truth_taxonomy.csv"

# These are the six lineage-root queries used to build the benchmark. Each name
# is the UniProt scientific name of its root taxon, verified by
# scripts/audit_sp_ground_truth_taxonomy.py, and the roots are mutually
# exclusive so every taxon resolves to exactly one of them.
#
# Note that 8504 is Lepidosauria (squamates and the tuatara), not Reptilia.
# Reptilia is only a UniProt alias for it, and true Reptilia (Sauropsida,
# 8457) would overlap Aves. Testudines and Crocodylia are therefore absent
# from the corpus; extending coverage to them is tracked in issue #71.
SOURCE_CLADES = (
    ("Mammalia", 40674),
    ("Aves", 8782),
    ("Actinopterygii", 7898),
    ("Lepidosauria", 8504),
    ("Amphibia", 8292),
    ("Chondrichthyes", 7777),
)
SOURCE_CLADE_TO_CATEGORY = {
    "Aves": "bird",
    "Actinopterygii": "fish",
    "Chondrichthyes": "fish",
    "Lepidosauria": "other_vertebrate",
    "Amphibia": "other_vertebrate",
}

# The finer mammalian reporting categories intentionally preserve the package's
# established semantics: "ungulate" covers bovids, suids, and equids, while
# other hoofed mammals and cetaceans remain "other_mammal"; "carnivore" covers
# canids and felids, as the package has historically done.
HUMAN_TAXON_ID = 9606
CATEGORY_LINEAGE_RULES = (
    ("nhp", (9443,)),  # Primates (human is handled first)
    ("murine", (39107,)),  # Murinae
    ("ungulate", (9895, 9821, 9788)),  # Bovidae, Suidae, Equidae
    ("carnivore", (9608, 9681)),  # Canidae, Felidae
    ("other_mammal", (40674,)),  # Mammalia
    ("bird", (8782,)),  # Aves
    ("fish", (7898, 7777)),  # Actinopterygii, Chondrichthyes
    ("other_vertebrate", (8504, 8292)),  # Lepidosauria, Amphibia
)


@lru_cache(maxsize=None)
def load_taxonomy_audit(path: Path = TAXONOMY_AUDIT_CSV) -> dict[str, dict[str, str]]:
    """Load the accession-independent UniProt taxon audit."""
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    result = {row["taxon_id"]: row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate taxon IDs in {path}")
    return result


def category_from_lineage(taxon_id: int, lineage_taxon_ids: set[int]) -> tuple[str, int]:
    """Return reporting category and supporting lineage root for one taxon."""
    lineage = set(lineage_taxon_ids)
    lineage.add(taxon_id)
    if HUMAN_TAXON_ID in lineage:
        return "human", HUMAN_TAXON_ID
    for category, roots in CATEGORY_LINEAGE_RULES:
        for root in roots:
            if root in lineage:
                return category, root
    raise ValueError(f"Taxon {taxon_id} is outside the benchmark vertebrate lineages")


def source_clade_from_lineage(taxon_id: int, lineage_taxon_ids: set[int]) -> tuple[str, int]:
    """Return the unique benchmark fetch clade containing one taxon."""
    lineage = set(lineage_taxon_ids)
    lineage.add(taxon_id)
    matches = [(name, root) for name, root in SOURCE_CLADES if root in lineage]
    if len(matches) != 1:
        raise ValueError(f"Expected one source clade for taxon {taxon_id}, found {matches!r}")
    return matches[0]


def species_category(
    organism: str,
    taxon_id: str = "",
    source_clade: str = "",
    *,
    audit: dict[str, dict[str, str]] | None = None,
) -> str:
    """Resolve an SP benchmark reporting category from audited taxonomy.

    The checked-in taxon audit is authoritative. A persisted source clade is
    sufficient for new non-mammalian rows. Name normalization is retained only
    as an offline compatibility fallback for rows outside this benchmark.
    """
    taxon_key = str(taxon_id or "").strip()
    source_clade = str(source_clade or "").strip()
    taxonomy = load_taxonomy_audit() if audit is None else audit
    decision = taxonomy.get(taxon_key)
    if decision:
        audited_clade = decision["source_clade"]
        if source_clade and source_clade != audited_clade:
            raise ValueError(f"Source clade mismatch for taxon {taxon_key}: row={source_clade!r}, audit={audited_clade!r}")
        return decision["species_category"]

    valid_clades = {name for name, _taxon_id in SOURCE_CLADES}
    if source_clade and source_clade not in valid_clades:
        raise ValueError(f"Unknown SP source clade: {source_clade!r}")
    if source_clade in SOURCE_CLADE_TO_CATEGORY:
        return SOURCE_CLADE_TO_CATEGORY[source_clade]

    from mhcseqs.species import extract_latin_binomial, normalize_mhc_species

    for candidate in (organism, extract_latin_binomial(organism)):
        category = normalize_mhc_species(candidate)
        if category and (source_clade != "Mammalia" or category not in {"bird", "fish", "other_vertebrate"}):
            return category
    return "other_mammal" if source_clade == "Mammalia" else "other_vertebrate"
