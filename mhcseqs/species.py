"""Species normalization for MHC alleles.

Provides a fine-grained species taxonomy with roll-up mappings
to MHC species categories (human, nhp, murine, ungulate, carnivore,
cetacean, other_mammal, bird, fish, other_vertebrate).

Ported from presto/data/vocab.py.
"""

from __future__ import annotations

import re
from typing import Dict, Optional, Tuple

# MHC species categories
MHC_SPECIES_CATEGORIES = [
    "human",
    "nhp",
    "murine",
    "ungulate",
    "carnivore",
    "cetacean",
    "other_mammal",
    "bird",
    "fish",
    "other_vertebrate",
]

# Fine-grained species
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
    "bat",
    "marsupial",
    "monotreme",
    "cetacean",
    "other_mammal",
    # Birds
    "chicken",
    "other_bird",
    # Fish
    "salmon",
    "zebrafish",
    "shark",
    "other_fish",
    # Reptiles & amphibians
    "reptile",
    "amphibian",
    # Other
    "other_vertebrate",
    "invertebrate",
    "viruses",
    "bacteria",
    "fungi",
    "archaea",
]

# Keyword patterns checked in order; first match wins.
# IMPORTANT: Order matters! More specific patterns must come before broader ones
# to avoid false positives (e.g., "horseshoe bat" matching "horse").
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
            "alouatta",
            "colobus",
            "lophocebus",
            "rhinopithecus",
            "sapajus",
            "theropithecus",
            "primate",
        ),
        "other_nhp",
    ),
    # Rodents
    (("mus musculus", "mouse", "c57bl", "balb/c"), "mouse"),
    (("rattus", "rat"), "rat"),
    (("murine", "h2-", "h-2"), "mouse"),
    # Bats — BEFORE horse/dog to avoid "horseshoe bat" → horse, "dogfish" → dog
    (
        (
            "chiroptera",
            "rhinolophus",
            "pteropus",
            "myotis",
            "eptesicus",
            "pipistrellus",
            "miniopterus",
            "rousettus",
            "desmodus",
            "artibeus",
            "carollia",
            "phyllostomus",
            "hipposideros",
            "noctilio",
            "tadarida",
            "molossus",
        ),
        "bat",
    ),
    # Monotremes — BEFORE duck/bird to avoid "duck-billed platypus" → bird
    (("ornithorhynchus", "platypus", "tachyglossus", "echidna", "monotreme"), "monotreme"),
    # Marsupials
    (
        (
            "marsupial",
            "monodelphis",
            "didelphis",
            "opossum",
            "macropus",
            "kangaroo",
            "wallaby",
            "phascolarctos",
            "koala",
            "sarcophilus",
            "devil",
            "vombatus",
            "wombat",
            "petaurus",
            "dasyurus",
            "antechinus",
            "potoroidae",
            "sminthopsis",
            "notamacropus",
            "wallabia",
        ),
        "marsupial",
    ),
    # Sharks and rays — BEFORE dog/fish to avoid "dogfish" → dog
    (
        (
            "chondrichthyes",
            "squalus",
            "carcharodon",
            "rhincodon",
            "triakis",
            "ginglymostoma",
            "heterodontus",
            "scyliorhinus",
            "chiloscyllium",
            "leucoraja",
            "raja",
            "pristis",
            "callorhinchus",
            "hydrolagus",
            "dogfish",
        ),
        "shark",
    ),
    # Mammals (ungulates, carnivores, cetaceans)
    (
        ("bos taurus", "bos ", "bovine", "cow", "cattle", "bola-", "bos grunniens", "bubalus"),
        "cattle",
    ),
    (("sus scrofa", "sus ", "porcine", "pig", "swine", "sla-", "phacochoerus", "warthog"), "pig"),
    (("equus", "equine", "horse", "ela-", "donkey"), "horse"),
    (("ovis", "ovine", "sheep", "ola-"), "sheep"),
    (("capra", "caprine", "goat"), "goat"),
    (("cervus", "odocoileus", "deer", "elk", "moose", "alces"), "other_mammal"),
    (("camelus", "camelid", "llama", "lama ", "vicugna", "alpaca"), "other_mammal"),
    (("bison",), "cattle"),
    (("canis lupus", "canis ", "canine", "dog", "dla-", "cuon", "lycaon", "vulpes", "fox", "dhole"), "dog"),
    (("felis", "feline", "cat"), "cat"),
    (("rabbit", "oryctolagus"), "rabbit"),
    (
        (
            "phoca",
            "seal",
            "sea lion",
            "zalophus",
            "halichoerus",
            "neomonachus",
            "pinniped",
            "mustela",
            "ferret",
            "mink",
            "neovison",
            "neogale",
            "ursus",
            "bear",
            "ailuropoda",
            "panda",
        ),
        "other_mammal",
    ),
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
            "steno",
            "grampus",
            "globicephala",
            "lagenorhynchus",
            "cephalorhynchus",
            "mesoplodon",
            "ziphius",
            "orcinus",
            "orca",
            "delphinapterus",
            "beluga",
            "eschrichtius",
            "hyperoodon",
            "inia",
            "lipotes",
            "monodon",
            "narwhal",
            "neophocaena",
            "phocoena",
            "pontoporia",
        ),
        "cetacean",
    ),
    (("mammal",), "other_mammal"),
    # Birds
    (("gallus", "chicken", "gaga-"), "chicken"),
    (
        (
            "duck",
            "turkey",
            "quail",
            "bird",
            "avian",
            "aves",
            "falcon",
            "eagle",
            "hawk",
            "owl",
            "penguin",
            "parrot",
            "pigeon",
            "sparrow",
            "finch",
            "warbler",
            "passeriformes",
            "anseriformes",
            "anas ",
            "taeniopygia",
            "ficedula",
            "parus",
            "corvus",
            "sturnus",
            "turdus",
            "passer",
            "columba",
            "struthio",
            "meleagris",
            "coturnix",
            "spheniscus",
        ),
        "other_bird",
    ),
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
    (
        (
            "*fish",
            "pisces",
            "teleost",
            "actinopterygii",
            "oryzias",
            "medaka",
            "gasterosteus",
            "stickleback",
            "gadus",
            "cod",
            "cyprinus",
            "carp",
            "oreochromis",
            "tilapia",
            "ictalurus",
            "catfish",
            "takifugu",
            "fugu",
            "tetraodon",
            "esox",
            "pike",
        ),
        "other_fish",
    ),
    # Reptiles
    (
        (
            "reptile",
            "reptilia",
            "squamata",
            "lacertilia",
            "serpentes",
            "lizard",
            "snake",
            "gecko",
            "iguana",
            "anolis",
            "pogona",
            "varanus",
            "sphenodon",
            "tuatara",
            "alligator",
            "crocodile",
            "crocodylus",
            "gavialis",
            "caiman",
            "crocodylia",
            "turtle",
            "tortoise",
            "testudines",
            "chelonia",
            "chrysemys",
            "trachemys",
            "gopherus",
            "terrapene",
        ),
        "reptile",
    ),
    # Amphibians
    (
        (
            "amphibian",
            "amphibia",
            "frog",
            "toad",
            "xenopus",
            "rana",
            "bufo",
            "salamander",
            "newt",
            "axolotl",
            "ambystoma",
            "caecilian",
            "anura",
            "caudata",
            "gymnophiona",
            "lithobates",
            "nanorana",
            "pelophylax",
        ),
        "amphibian",
    ),
    # Catch-all vertebrate
    (("vertebrate",), "other_vertebrate"),
    # Invertebrates
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

# Fine → MHC species category
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
    "cattle": "ungulate",
    "pig": "ungulate",
    "horse": "ungulate",
    "sheep": "ungulate",
    "goat": "ungulate",
    "dog": "carnivore",
    "cat": "carnivore",
    "rabbit": "other_mammal",
    "bat": "other_mammal",
    "marsupial": "other_mammal",
    "monotreme": "other_mammal",
    "cetacean": "cetacean",
    "other_mammal": "other_mammal",
    "chicken": "bird",
    "other_bird": "bird",
    "salmon": "fish",
    "zebrafish": "fish",
    "shark": "fish",
    "other_fish": "fish",
    "reptile": "other_vertebrate",
    "amphibian": "other_vertebrate",
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
    "rabbit": "RLA",
    "chicken": "Gaga",
    "salmon": "Sasa",
    "zebrafish": "Dare",
}


def _word_match(keyword: str, text: str) -> bool:
    """Match a keyword in text using word-boundary-aware logic.

    - Keywords ending with a space or hyphen are matched as prefixes (substring).
    - Keywords starting with ``*`` are matched as word suffixes
      (e.g., ``*fish`` matches "rockfish", "swordfish", "fish").
    - All other keywords are matched as whole words using regex word boundaries,
      preventing false positives like "horse" matching "horseshoe".
    """
    if keyword.startswith("*"):
        # Suffix match: the word must end with this
        suffix = keyword[1:]
        return bool(re.search(re.escape(suffix) + r"\b", text))
    if keyword.endswith((" ", "-")):
        # Explicit prefix/substring match (the trailing char is intentional)
        return keyword in text
    # Use word boundary matching
    return bool(re.search(r"\b" + re.escape(keyword) + r"\b", text))


def normalize_species(raw: Optional[str]) -> Optional[str]:
    """Normalize a species string to a fine-grained category.

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
        if any(_word_match(kw, s) for kw in keywords):
            return label
    return None


def normalize_mhc_species(raw: Optional[str]) -> Optional[str]:
    """Normalize to 10-class MHC species category.

    Returns one of: human, nhp, murine, ungulate, carnivore, cetacean,
    other_mammal, bird, fish, other_vertebrate.
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
