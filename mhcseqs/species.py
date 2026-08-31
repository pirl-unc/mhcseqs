"""Species normalization for MHC alleles.

Provides a fine-grained species taxonomy with roll-up mappings
to MHC species categories (human, nhp, murine, ungulate, carnivore,
other_mammal, bird, fish, other_vertebrate).

Ported from presto/data/vocab.py.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Optional, Tuple

# MHC species categories
MHC_SPECIES_CATEGORIES = [
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
    "cetacean",  # fine-grained, maps to other_mammal
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
        "other_mammal",
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
    "cetacean": "other_mammal",
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

# Externally established MHC naming prefixes for explicit species or genus
# scopes. Exact species entries take precedence over genus-wide systems.
#
# These are input/source aliases, not mhcseqs-generated species identities.
# HLA is maintained by IPD-IMGT/HLA; the non-human systems are represented by
# IPD-MHC or established model-organism nomenclature.  Do not add a mechanical
# 2+2/4+4 abbreviation here: new canonical output uses the self-describing
# scientific-name alias returned by ``get_canonical_prefix``.
#
ESTABLISHED_MHC_PREFIXES: Dict[str, str] = {
    "Homo sapiens": "HLA",
    "Macaca mulatta": "Mamu",
    "Macaca fascicularis": "Mafa",
    "Pan troglodytes": "Patr",
    "Gorilla gorilla": "Gogo",
    "Pongo pygmaeus": "Popy",
    "Papio anubis": "Paan",
    "Mus musculus": "H2",
    "Rattus": "RT1",
    "Bos": "BoLA",
    "Sus": "SLA",
    "Equus": "ELA",
    "Ovis": "OLA",
    "Capra": "CLA",
    "Canis": "DLA",
    "Felis": "FLA",
    "Oryctolagus": "RLA",
    "Gallus gallus": "Gaga",
    "Salmo salar": "Sasa",
    "Danio rerio": "Dare",
}

# Every short code in this compatibility table has an external database or
# literature citation. The complete source-alias registry is loaded below;
# additions to either require evidence, since mechanical 2+2/4+4/5+5
# abbreviations are not evidence. Several prefixes share a source registry URL.
ESTABLISHED_MHC_PREFIX_SOURCES: Dict[str, str] = {
    "HLA": "https://hla.alleles.org/genes/index.html",
    "Mamu": "https://www.ebi.ac.uk/ipd/mhc/taxonomy/",
    "Mafa": "https://www.ebi.ac.uk/ipd/mhc/taxonomy/",
    "Patr": "https://www.ebi.ac.uk/ipd/mhc/group/NHP/species/",
    "Gogo": "https://www.ebi.ac.uk/ipd/mhc/taxonomy/",
    "Popy": "https://www.ebi.ac.uk/ipd/mhc/taxonomy/",
    "Paan": "https://www.ebi.ac.uk/ipd/mhc/taxonomy/",
    "H2": "https://www.informatics.jax.org/downloads/datasets/misc/H2Haplotypes/H2_haplotypes.html",
    "RT1": "https://www.ebi.ac.uk/ipd/mhc/",
    "BoLA": "https://www.ebi.ac.uk/ipd/mhc/group/BoLA/nomenclature/",
    "SLA": "https://www.ebi.ac.uk/ipd/mhc/",
    "ELA": "https://www.ebi.ac.uk/ipd/mhc/",
    "OLA": "https://www.ebi.ac.uk/ipd/mhc/",
    "CLA": "https://pubmed.ncbi.nlm.nih.gov/8134633/",
    "DLA": "https://www.ebi.ac.uk/ipd/mhc/",
    "FLA": "https://pubmed.ncbi.nlm.nih.gov/2492667/",
    "RLA": "https://pubmed.ncbi.nlm.nih.gov/32522857/",
    "Gaga": "https://www.ebi.ac.uk/ipd/mhc/group/CHICKEN/species/",
    "Sasa": "https://www.ebi.ac.uk/ipd/mhc/group/FISH/",
    "Dare": "https://pubmed.ncbi.nlm.nih.gov/18439319/",
}


@dataclass(frozen=True)
class MhcPrefixAlias:
    """One externally attested MHC prefix and its taxonomic scope."""

    species: str
    prefix: str
    scope: str
    status: str
    evidence: str
    notes: str = ""


_MHC_PREFIX_REGISTRY_PATH = Path(__file__).resolve().parent / "mhc_prefix_aliases.csv"
_PREFIX_STATUS_PRIORITY = {
    "ipd_current": 0,
    "formal_system": 1,
    "literature_historical": 2,
    "external_database": 3,
}


@lru_cache(maxsize=1)
def get_mhc_prefix_registry() -> Tuple[MhcPrefixAlias, ...]:
    """Return every source-attested short prefix shipped by mhcseqs.

    The registry contains current IPD-MHC designations, historical systems
    from nomenclature papers, and spellings present in external source
    databases. It intentionally excludes mechanically generated aliases.
    """
    with _MHC_PREFIX_REGISTRY_PATH.open(encoding="utf-8") as handle:
        return tuple(MhcPrefixAlias(**row) for row in csv.DictReader(handle))


@lru_cache(maxsize=1)
def _mhc_prefixes_by_spelling() -> Dict[str, Tuple[MhcPrefixAlias, ...]]:
    """Index the source registry once for case-insensitive prefix matching."""
    grouped: Dict[str, list[MhcPrefixAlias]] = {}
    for entry in get_mhc_prefix_registry():
        grouped.setdefault(entry.prefix.casefold(), []).append(entry)
    return {prefix: tuple(entries) for prefix, entries in grouped.items()}


def _normalized_taxon_name(raw: Optional[str]) -> str:
    """Normalize whitespace/case while retaining an explicit subspecies."""
    latin = get_latin_name(raw).split("(", 1)[0].strip()
    return re.sub(r"\s+", " ", latin).casefold()


def _prefix_alias_applies(entry: MhcPrefixAlias, raw_species: Optional[str]) -> bool:
    taxon = _normalized_taxon_name(raw_species)
    if not taxon:
        return False
    if entry.scope == "genus":
        registered_genus = entry.species.strip().casefold()
        return taxon.split(" ", 1)[0] == registered_genus
    registered = _normalized_taxon_name(entry.species)
    return taxon == registered


def get_established_mhc_prefixes(raw: Optional[str]) -> Tuple[str, ...]:
    """Return all current and historical source aliases for a species.

    Current IPD designations and formal system names sort before historical
    and raw-database spellings. Prefixes are compared case-insensitively, but
    the evidence-backed spelling is retained.
    """
    matches = [entry for entry in get_mhc_prefix_registry() if _prefix_alias_applies(entry, raw)]
    matches.sort(key=lambda entry: (_PREFIX_STATUS_PRIORITY[entry.status], entry.prefix.casefold()))
    result = []
    seen = set()
    for entry in matches:
        normalized = entry.prefix.casefold()
        if normalized not in seen:
            result.append(entry.prefix)
            seen.add(normalized)
    return tuple(result)


def find_mhc_prefix_aliases(
    value: Optional[str],
    *,
    species: Optional[str] = None,
) -> Tuple[str, str, Tuple[MhcPrefixAlias, ...]]:
    """Split an attested prefix from a molecule name.

    Returns ``(prefix, body, registry_entries)``. With species context, only
    aliases explicitly attested for that taxon are returned. Prefixes may be
    followed by ``-``/``_`` or directly by an uppercase gene token, matching
    conventions present in UniProt and older papers.
    """
    token = str(value or "").strip()
    if not token:
        return "", "", ()

    by_spelling = _mhc_prefixes_by_spelling()

    for normalized_prefix in sorted(by_spelling, key=len, reverse=True):
        prefix = by_spelling[normalized_prefix][0].prefix
        if token[: len(prefix)].casefold() != normalized_prefix:
            continue
        remainder = token[len(prefix) :]
        if remainder.startswith(("-", "_")):
            body = remainder[1:]
        elif len(prefix) >= 2 and remainder and remainder[0].isupper():
            body = remainder
        else:
            continue
        if not body:
            continue
        entries = by_spelling[normalized_prefix]
        if species:
            entries = tuple(entry for entry in entries if _prefix_alias_applies(entry, species))
        return prefix, body, entries
    return "", "", ()


# Backward-compatible category mapping. Despite its historical name, this is
# an input-alias table; canonical output comes from ``get_canonical_prefix``.
CANONICAL_MHC_PREFIXES: Dict[str, str] = {
    fine: ESTABLISHED_MHC_PREFIXES[latin] for fine, latin in LATIN_NAMES.items() if latin in ESTABLISHED_MHC_PREFIXES
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
        # Explicit prefix match (the trailing char is intentional); anchor at a
        # word boundary so e.g. "sus " only matches Sus and not the "-sus "
        # tail of "cynoglossus semilaevis".
        return bool(re.search(r"\b" + re.escape(keyword), text))
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
    """Normalize to 9-class MHC species category.

    Returns one of: human, nhp, murine, ungulate, carnivore,
    other_mammal, bird, fish, other_vertebrate.
    Returns None for non-animal or unrecognizable inputs.
    """
    fine = normalize_species(raw)
    if fine is None:
        return None
    return FINE_TO_MHC_SPECIES.get(fine)


def get_latin_name(raw: Optional[str]) -> str:
    """Return the canonical Latin name for a species, or the raw string if unknown."""
    # Preserve an explicit scientific name before applying broad category
    # matching. Otherwise, for example, Macaca fascicularis would collapse to
    # the category default Macaca mulatta.
    if raw:
        explicit = str(raw).split("(", 1)[0].strip()
        if re.match(r"^[A-Z][A-Za-z-]+\s+[a-z][A-Za-z.-]+(?:\s|$)", explicit):
            return explicit
    fine = normalize_species(raw)
    if fine and fine in LATIN_NAMES:
        return LATIN_NAMES[fine]
    # If the raw string looks like a Latin name already (two capitalized words), keep it
    if raw and " " in str(raw).strip():
        return str(raw).strip()
    return str(raw or "")


def extract_latin_binomial(organism: Optional[str]) -> str:
    """Extract the latin binomial from a UniProt organism string.

    Handles formats like:
      "Crocodylus porosus (Saltwater crocodile) (Estuarine crocodile)"
      "Homo sapiens (Human)"
      "Mus musculus"

    Returns the genus + species epithet (e.g. "Crocodylus porosus"),
    or the original string stripped if no parenthetical is found.
    """
    if not organism:
        return ""
    # Strip everything from first '(' onward, then clean up
    binomial = str(organism).split("(")[0].strip()
    # Take only genus + species (first two words) to drop subspecies
    parts = binomial.split()
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return binomial


def full_species_name_alias(raw: Optional[str]) -> str:
    """Return a self-describing concatenated scientific-name alias.

    Binomials use both words (``Homo sapiens`` -> ``HomoSapiens``). Generic
    taxon labels retain ``Sp`` (``Coregonus sp.`` -> ``CoregonusSp``), and a
    genus-only label is preserved rather than shortened. This function never
    invents a 2+2, 4+4, or 5+5 code.
    """
    latin = extract_latin_binomial(get_latin_name(raw))
    if not latin:
        return ""
    words = re.findall(r"[A-Za-z]+", latin)
    if not words:
        return ""
    return "".join(word[:1].upper() + word[1:] for word in words[:2])


def get_established_mhc_prefix(raw: Optional[str]) -> str:
    """Return the preferred source-attested short prefix, or an empty string."""
    prefixes = get_established_mhc_prefixes(raw)
    return prefixes[0] if prefixes else ""


def get_canonical_prefix(raw: Optional[str]) -> str:
    """Return the canonical self-describing species alias.

    Short MHC codes are accepted as external/source aliases, but canonical
    mhcseqs output uses the full scientific name so generated abbreviations do
    not become new species identifiers.
    """
    return full_species_name_alias(raw)
