"""Taxonomy rules and audited metadata for the SP ground-truth corpus."""

from __future__ import annotations

import csv
import re
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TAXONOMY_AUDIT_CSV = ROOT / "data" / "sp_ground_truth_taxonomy.csv"

# These mutually exclusive roots partition the Vertebrata candidate artifact.
# Each label is the UniProt scientific name of its root taxon.  Broad parent
# roots such as Sauropsida and Sarcopterygii cannot be used because they would
# overlap narrower benchmark roots (Aves and the tetrapod clades).
SOURCE_CLADES = (
    ("Mammalia", 40674),
    ("Aves", 8782),
    ("Actinopterygii", 7898),
    ("Lepidosauria", 8504),
    ("Amphibia", 8292),
    ("Chondrichthyes", 7777),
    ("Testudines", 8459),
    ("Crocodylia", 1294634),
    ("Coelacanthimorpha", 118072),
    ("Dipnomorpha", 7878),
    ("Myxini", 117565),
)
MAMMAL_CATEGORIES = frozenset({"human", "nhp", "murine", "ungulate", "carnivore", "other_mammal"})
SOURCE_CLADE_TO_CATEGORY = {
    "Aves": "bird",
    "Actinopterygii": "fish",
    "Chondrichthyes": "fish",
    "Lepidosauria": "other_vertebrate",
    "Amphibia": "other_vertebrate",
    "Testudines": "other_vertebrate",
    "Crocodylia": "other_vertebrate",
    "Coelacanthimorpha": "fish",
    "Dipnomorpha": "fish",
    "Myxini": "fish",
}

# The finer mammalian reporting categories intentionally preserve the package's
# established semantics: "ungulate" covers bovids, suids, and equids, while
# other hoofed mammals and cetaceans remain "other_mammal"; "carnivore" covers
# canids and felids, as the package has historically done.
#
# ORDER IS LOAD-BEARING: category_from_lineage returns the first rule whose root
# appears in the lineage, so roots must run specific -> general. "other_mammal"
# (Mammalia) matches every mammal, so any finer mammal rule placed after it
# would be dead code. test_category_rule_order_is_specific_to_general pins the
# order so inserting a rule is a deliberate, reviewed change.
CATEGORY_LINEAGE_RULES = (
    ("human", (9606,)),  # Homo sapiens, including subspecies
    ("nhp", (9443,)),  # Primates
    ("murine", (39107,)),  # Murinae
    ("ungulate", (9895, 9821, 9788)),  # Bovidae, Suidae, Equidae
    ("carnivore", (9608, 9681)),  # Canidae, Felidae
    ("other_mammal", (40674,)),  # Mammalia
    ("bird", (8782,)),  # Aves
    ("fish", (7898, 7777, 118072, 7878, 117565)),
    ("other_vertebrate", (8504, 8292, 8459, 1294634)),
)


@lru_cache(maxsize=None)
def load_taxonomy_audit(path: Path = TAXONOMY_AUDIT_CSV) -> dict[str, dict[str, str]]:
    """Load the accession-independent UniProt taxon audit.

    Raises if the audit is absent. It is the authoritative source for every
    benchmark category, and returning an empty mapping here would silently
    demote audited rows to the name-matching fallback this module exists to
    replace. Callers that genuinely want the fallback pass ``audit={}``.
    """
    if not path.exists():
        raise FileNotFoundError(f"SP taxonomy audit is missing: {path}. Run scripts/audit_sp_ground_truth_taxonomy.py")
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


def classification_from_taxonomy_row(row: dict[str, str]) -> tuple[str, int, str, int]:
    """Return source clade and reporting category from one snapshot row."""
    taxon_id = int(row["Taxon ID"])
    lineage_ids = {int(value) for value in re.findall(r"\d+", row.get("Taxonomic lineage (Ids)", ""))}
    source_clade, source_root = source_clade_from_lineage(taxon_id, lineage_ids)
    category, category_root = category_from_lineage(taxon_id, lineage_ids)
    return source_clade, source_root, category, category_root


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

    # Only two cases reach here: source_clade == "Mammalia", or no clade at all.
    from mhcseqs.species import extract_latin_binomial, normalize_mhc_species

    candidates = (organism, extract_latin_binomial(organism))
    if source_clade == "Mammalia":
        # A name match that claims a non-mammalian category is wrong by
        # construction, so only mammalian answers are accepted.
        for candidate in candidates:
            category = normalize_mhc_species(candidate)
            if category in MAMMAL_CATEGORIES:
                return category
        return "other_mammal"
    for candidate in candidates:
        category = normalize_mhc_species(candidate)
        if category:
            return category
    return "other_vertebrate"
