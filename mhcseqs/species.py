"""Species normalization for MHC alleles.

Provides a 29-class fine-grained species taxonomy with roll-up mappings
to coarser 7-class MHC species categories (human, nhp, murine,
other_mammal, bird, fish, other_vertebrate).

Ported from presto/data/vocab.py.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

# 7-class MHC species categories
MHC_SPECIES_CATEGORIES = [
    "human",
    "nhp",
    "murine",
    "other_mammal",
    "bird",
    "fish",
    "other_vertebrate",
]

# 29-class fine-grained species
FINE_SPECIES = [
    # Primates
    "human",
    "macaque",
    "chimpanzee",
    "gorilla",
    "orangutan",
    "baboon",
    "other_nhp",
    # Rodents
    "mouse",
    "rat",
    # Mammals
    "cattle",
    "pig",
    "horse",
    "sheep",
    "goat",
    "dog",
    "cat",
    "rabbit",
    "other_mammal",
    # Birds
    "chicken",
    "other_bird",
    # Fish
    "salmon",
    "zebrafish",
    "other_fish",
    # Other
    "other_vertebrate",
    "invertebrate",
    "viruses",
    "bacteria",
    "fungi",
    "archaea",
]

# Keyword patterns checked in order; first match wins.
_SPECIES_PATTERNS: list[Tuple[Tuple[str, ...], str]] = [
    # Human
    (("homo sapiens", "human"), "human"),
    # NHP
    (("chimpanzee", "pan troglodytes", "pan paniscus", "patr-"), "chimpanzee"),
    (("gorilla", "gogo-"), "gorilla"),
    (("orangutan", "pongo", "popy-"), "orangutan"),
    (("baboon", "papio", "paan-"), "baboon"),
    (("macaque", "macaca", "rhesus", "mamu-"), "macaque"),
    (
        (
            "nhp",
            "aotus",
            "night monkey",
            "aona-",
            "cercopithecus",
            "saguinus",
            "callithrix",
            "saimiri",
            "ateles",
            "pithecia",
            "leontopithecus",
            "hylobates",
            "chlorocebus",
            "cercocebus",
            "mandrillus",
            "callicebus",
            "cebus",
            "semnopithecus",
            "primate",
        ),
        "other_nhp",
    ),
    # Rodents
    (("mus musculus", "mouse", "c57bl", "balb/c"), "mouse"),
    (("rattus", "rat "), "rat"),
    (("murine", "h2-", "h-2"), "mouse"),
    # Mammals
    (
        ("bos taurus", "bos ", "bovine", "cow", "cattle", "bola-", "bos grunniens", "bubalus"),
        "cattle",
    ),
    (("sus scrofa", "sus ", "porcine", "pig", "swine", "sla-"), "pig"),
    (("equus", "equine", "horse", "ela-"), "horse"),
    (("ovis", "ovine", "sheep", "ola-"), "sheep"),
    (("capra", "caprine", "goat"), "goat"),
    (("canis", "canine", "dog", "dla-"), "dog"),
    (("felis", "feline", "cat "), "cat"),
    (("rabbit", "oryctolagus"), "rabbit"),
    (
        (
            "whale",
            "dolphin",
            "porpoise",
            "cetacea",
            "cetacean",
            "balaenoptera",
            "eubalaena",
            "megaptera",
            "kogia",
            "tursiops",
            "delphinus",
            "stenella",
            "steno ",
            "grampus",
            "globicephala",
            "lagenorhynchus",
            "cephalorhynchus",
            "mesoplodon",
            "ziphius",
        ),
        "other_mammal",
    ),
    (("mammal",), "other_mammal"),
    # Birds
    (("gallus", "chicken", "gaga-"), "chicken"),
    (("duck", "turkey", "quail", "bird", "avian", "aves"), "other_bird"),
    # Pathogens
    (
        (
            "virus",
            "viral",
            "influenza",
            "sars",
            "cov",
            "hiv",
            "hcv",
            "hbv",
            "ebv",
            "cmv",
            "hsv",
            "vzv",
            "htlv",
            "dengue",
            "zika",
            "ebola",
            "measles",
            "hepatitis",
            "retrovirus",
            "coronavirus",
            "adenovirus",
            "papillomavirus",
            "herpes",
            "vaccinia",
            "poxvirus",
            "flavivirus",
            "paramyxovirus",
            "orthomyxovirus",
            "phage",
            "bacteriophage",
        ),
        "viruses",
    ),
    (
        (
            "mycobacterium",
            "tuberculosis",
            "escherichia",
            "e. coli",
            "staphylococcus",
            "streptococcus",
            "salmonella",
            "clostridium",
            "listeria",
            "helicobacter",
            "chlamydia",
            "borrelia",
            "treponema",
            "pseudomonas",
            "bacillus",
            "legionella",
            "neisseria",
            "rickettsia",
            "bartonella",
            "bacterium",
            "bacteria",
            "bacterial",
        ),
        "bacteria",
    ),
    (
        (
            "candida",
            "aspergillus",
            "cryptococcus",
            "coccidioides",
            "histoplasma",
            "blastomyces",
            "saccharomyces",
            "yeast",
            "fungus",
            "fungi",
            "fungal",
            "pneumocystis",
            "trichophyton",
        ),
        "fungi",
    ),
    (
        ("archaea", "archaeal", "methanobacterium", "halobacterium", "sulfolobus", "thermococcus"),
        "archaea",
    ),
    # Fish (after bacteria so "salmonella" doesn't match "salmon")
    (("salmo salar", "salmo ", "salmon", "trout", "oncorhynchus"), "salmon"),
    (("danio", "zebrafish"), "zebrafish"),
    (("fish", "pisces"), "other_fish"),
    # Other
    (
        (
            "reptile",
            "reptilia",
            "amphibian",
            "amphibia",
            "frog",
            "xenopus",
            "turtle",
            "lizard",
            "snake",
            "alligator",
            "crocodile",
            "salamander",
        ),
        "other_vertebrate",
    ),
    (
        (
            "drosophila",
            "insect",
            "arthropod",
            "arachnid",
            "mosquito",
            "tick",
            "worm",
            "nematode",
            "mollusk",
            "caenorhabditis",
            "c. elegans",
            "invertebrate",
            "schistosoma",
            "plasmodium",
            "toxoplasma",
            "leishmania",
            "trypanosoma",
            "parasite",
        ),
        "invertebrate",
    ),
]

FINE_SPECIES_SET = set(FINE_SPECIES)

# Fine → 7-class MHC species category
FINE_TO_MHC_SPECIES: Dict[str, Optional[str]] = {
    "human": "human",
    "macaque": "nhp",
    "chimpanzee": "nhp",
    "gorilla": "nhp",
    "orangutan": "nhp",
    "baboon": "nhp",
    "other_nhp": "nhp",
    "mouse": "murine",
    "rat": "murine",
    "cattle": "other_mammal",
    "pig": "other_mammal",
    "horse": "other_mammal",
    "sheep": "other_mammal",
    "goat": "other_mammal",
    "dog": "other_mammal",
    "cat": "other_mammal",
    "rabbit": "other_mammal",
    "other_mammal": "other_mammal",
    "chicken": "bird",
    "other_bird": "bird",
    "salmon": "fish",
    "zebrafish": "fish",
    "other_fish": "fish",
    "other_vertebrate": "other_vertebrate",
}
# Non-animal categories → None
for _fs in FINE_SPECIES:
    if _fs not in FINE_TO_MHC_SPECIES:
        FINE_TO_MHC_SPECIES[_fs] = None

# Canonical Latin names for each fine-grained species
LATIN_NAMES: Dict[str, str] = {
    "human": "Homo sapiens",
    "macaque": "Macaca mulatta",
    "chimpanzee": "Pan troglodytes",
    "gorilla": "Gorilla gorilla",
    "orangutan": "Pongo pygmaeus",
    "baboon": "Papio anubis",
    "mouse": "Mus musculus",
    "rat": "Rattus norvegicus",
    "cattle": "Bos taurus",
    "pig": "Sus scrofa",
    "horse": "Equus caballus",
    "sheep": "Ovis aries",
    "goat": "Capra hircus",
    "dog": "Canis lupus familiaris",
    "cat": "Felis catus",
    "rabbit": "Oryctolagus cuniculus",
    "chicken": "Gallus gallus",
    "salmon": "Salmo salar",
    "zebrafish": "Danio rerio",
}

# Canonical MHC naming prefixes per fine-grained species
CANONICAL_MHC_PREFIXES: Dict[str, str] = {
    "human": "HLA",
    "macaque": "Mamu",
    "chimpanzee": "Patr",
    "gorilla": "Gogo",
    "orangutan": "Popy",
    "baboon": "Paan",
    "mouse": "H2",
    "rat": "RT1",
    "cattle": "BoLA",
    "pig": "SLA",
    "horse": "ELA",
    "sheep": "OLA",
    "goat": "CLA",
    "dog": "DLA",
    "cat": "FLA",
    "rabbit": "OrCu",
    "chicken": "Gaga",
    "salmon": "Sasa",
    "zebrafish": "Dare",
}


def normalize_species(raw: Optional[str]) -> Optional[str]:
    """Normalize a species string to one of 29 fine-grained categories.

    Returns None if unrecognizable.
    """
    if raw is None:
        return None
    s = str(raw).strip().lower()
    if not s:
        return None
    if s in FINE_SPECIES_SET:
        return s
    for keywords, label in _SPECIES_PATTERNS:
        if any(kw in s for kw in keywords):
            return label
    return None


def normalize_mhc_species(raw: Optional[str]) -> Optional[str]:
    """Normalize to 7-class MHC species category.

    Returns one of: human, nhp, murine, other_mammal, bird, fish, other_vertebrate.
    Returns None for non-animal or unrecognizable inputs.
    """
    fine = normalize_species(raw)
    if fine is None:
        return None
    return FINE_TO_MHC_SPECIES.get(fine)


def get_latin_name(raw: Optional[str]) -> str:
    """Return the canonical Latin name for a species, or the raw string if unknown."""
    fine = normalize_species(raw)
    if fine and fine in LATIN_NAMES:
        return LATIN_NAMES[fine]
    # If the raw string looks like a Latin name already (two capitalized words), keep it
    if raw and " " in str(raw).strip():
        return str(raw).strip()
    return str(raw or "")


def get_canonical_prefix(raw: Optional[str]) -> str:
    """Return the canonical MHC naming prefix for a species, or empty string."""
    fine = normalize_species(raw)
    if fine and fine in CANONICAL_MHC_PREFIXES:
        return CANONICAL_MHC_PREFIXES[fine]
    return ""
