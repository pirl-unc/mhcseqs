"""Allele name parsing and normalization via mhcgnomes.

Provides functions to parse MHC allele names from FASTA headers, normalize them
to canonical two-field resolution, and infer gene/class/species metadata.

Ported from presto/data/allele_resolver.py.
"""

from __future__ import annotations

import re
from typing import Any, Optional

import mhcgnomes
from mhcgnomes.result import Result as MhcgnomesResult

from .species import full_species_name_alias, normalize_mhc_species

# ---------------------------------------------------------------------------
# MHC class normalization
# ---------------------------------------------------------------------------

_MHC_CLASS_I_ALIASES = {
    "I",
    "IA",  # mhcgnomes subclass "Ia" (class I alpha), NOT mouse I-A
    "IB",  # mhcgnomes subclass "Ib" (class I beta / non-classical)
    "IC",
    "CLASSI",
    "CLASS-I",
    "MHCI",
    "MHC-I",
}
_MHC_CLASS_II_ALIASES = {
    "II",
    "IIA",
    "IIB",
    "CLASSII",
    "CLASS-II",
    "MHCII",
    "MHC-II",
}

MOLECULE_RESULT_TYPES = (
    mhcgnomes.Allele,
    mhcgnomes.Gene,
    mhcgnomes.AlleleWithoutGene,
    mhcgnomes.Pair,
)
MOLECULE_AND_SEROTYPE_RESULT_TYPES = (*MOLECULE_RESULT_TYPES, mhcgnomes.Serotype)


def normalize_mhc_class(value: Optional[str], default: Optional[str] = None) -> Optional[str]:
    """Normalize MHC class labels to canonical "I" / "II"."""
    if value is None:
        return default
    normalized = str(value).strip().upper().replace("_", "").replace(" ", "")
    normalized = normalized.replace("/", "").replace("*", "")
    if normalized in _MHC_CLASS_I_ALIASES or (normalized.startswith("I") and not normalized.startswith("II")):
        return "I"
    if normalized in _MHC_CLASS_II_ALIASES or normalized.startswith("II"):
        return "II"
    return default


# ---------------------------------------------------------------------------
# mhcgnomes delegation
# ---------------------------------------------------------------------------


def is_generated_short_prefix(candidate: str, parsed: Any) -> bool:
    """Return whether *candidate* uses a mechanically truncated mhcgnomes prefix.

    Full scientific-name aliases (``HomoSapiens``, ``CyanistesCaeruleus``) and
    source-attested MHC prefixes (``HLA``, ``Mamu``, ``Acsi``) return ``False``.
    Mechanically truncated derived prefixes such as 4+4 shorthands
    (``CyanCaer``, ``HomoSapi``) return ``True``.
    """
    species = getattr(parsed, "species", None)
    species_name = getattr(species, "name", None)
    if not species_name or "-" not in candidate:
        return False

    prefix = candidate.split("-", 1)[0]
    normalized_prefix = re.sub(r"[^A-Za-z0-9]+", "", prefix).casefold()
    if normalized_prefix == full_species_name_alias(species_name).casefold():
        return False

    return bool(mhcgnomes.species.prefix_is_derived_for(prefix, species_name))


def parse_allele_name(
    allele: Optional[str],
    *,
    species: Optional[str] = None,
    require_explicit_species: bool = False,
    allow_serotypes: bool = False,
) -> Optional[Any]:
    """Parse an allele string with mhcgnomes.

    Parameters
    ----------
    allele : str
        Allele name to parse (e.g. "HLA-A*02:01", "Crpo-UA", "H2-Kb").
    species : str, optional
        Latin binomial of the species (e.g. "Homo sapiens", "Crocodylus porosus").
        When provided, passed to mhcgnomes as the ``species`` parameter so the
        parser can validate and disambiguate the allele for that organism.
    require_explicit_species : bool, optional
        Reject names whose species was inferred or supplied only by a parser
        default. This is useful for validating untrusted curation inputs.
    allow_serotypes : bool, optional
        When True, also accept ``mhcgnomes.Serotype`` results (for example
        ``"HLA-A2"`` or ``"B27"``) explicitly as serotypes instead of
        rejecting them.
    """
    if not allele:
        return None
    raw = str(allele).strip()
    if not raw:
        return None

    result_types = MOLECULE_AND_SEROTYPE_RESULT_TYPES if allow_serotypes else MOLECULE_RESULT_TYPES
    kwargs: dict[str, Any] = {
        "required_result_types": result_types,
        "raise_on_error": False,
        "require_explicit_species": require_explicit_species,
    }
    if species:
        kwargs["species"] = species
    elif require_explicit_species:
        # Keep strict provenance available without changing the ordinary
        # mhcgnomes behavior where a bare allele such as A*02:01 means HLA.
        kwargs["default_species"] = None

    parsed = mhcgnomes.parse(raw, **kwargs)
    if parsed is not None and is_generated_short_prefix(raw, parsed):
        return None
    return parsed


def canonical_allele_name(parsed: Any, allele_fields: Optional[int] = None) -> str:
    """Render an mhcgnomes parse result as an mhcseqs canonical allele name.

    The species is written as its full binomial alias (for example
    ``HomoSapiens`` rather than ``HLA``) so that names from different
    species do not collide on a shared short prefix.

    Parameters
    ----------
    parsed : mhcgnomes result
        Object returned by :func:`parse_allele_name` (an mhcgnomes ``Allele``,
        ``Gene``, ``AlleleWithoutGene``, ``Pair``, or ``Serotype``).
    allele_fields : int, optional
        When given, restrict the allele to this many fields (minimum 1)
        using mhcgnomes ``restrict_allele_fields`` before rendering.

    Returns
    -------
    str
        ``"<FullSpeciesAlias>-<name without species>"`` when the species
        has a full binomial alias
        (:func:`mhcseqs.species.full_species_name_alias` is non-empty) and the
        name without species is non-empty; otherwise mhcgnomes' own
        ``to_string()``.
    """
    if allele_fields is not None:
        target_fields = max(1, int(allele_fields))
        parsed = parsed.restrict_allele_fields(target_fields)

    species_name = parsed.species.name if parsed.species is not None else None
    species_alias = full_species_name_alias(species_name)
    if species_alias:
        body = str(parsed.to_string(include_species=False))
        if body:
            return f"{species_alias}-{body}"
    return str(parsed.to_string())


def normalize_allele_name(name: str) -> str:
    """Normalize an allele name to canonical two-field protein resolution.

    Examples:
        "HLA-A*02:01" -> "HomoSapiens-A*02:01"
        "A*02:01" -> "HomoSapiens-A*02:01"
        "A0201" -> "HomoSapiens-A*02:01"
        "HLA-A*02:01:01:02L" -> "HomoSapiens-A*02:01L"
    """
    parsed = parse_allele_name(name)
    if parsed is None:
        raise ValueError(f"mhcgnomes failed to parse allele: {name!r}")
    return canonical_allele_name(parsed, allele_fields=2)


def infer_gene(allele: Optional[str]) -> Optional[str]:
    """Extract gene name from allele via mhcgnomes (e.g. "A", "DRB1")."""
    if not allele:
        return None
    parsed = parse_allele_name(allele)
    gene = getattr(parsed, "gene", None)
    if gene is not None and gene.name:
        return str(gene.name).upper()
    return None


def infer_mhc_class(allele: Optional[str], *, species: Optional[str] = None) -> Optional[str]:
    """Infer MHC class ("I" or "II") from allele name via mhcgnomes.

    Uses ``parse_gene_class()`` first (species-aware classification), then
    falls back to full allele parsing.
    """
    if not allele:
        return None
    result = parse_gene_class(allele, species=species)
    if result is not None:
        cls = result.get("mhc_class")
        if cls in ("I", "II"):
            return cls
    parsed = parse_allele_name(allele, species=species)
    if parsed is None:
        return None
    return normalize_mhc_class(getattr(parsed, "mhc_class", None), default=None)


def parse_gene_class(gene: Optional[str], *, species: Optional[str] = None) -> Optional[dict]:
    """Classify a gene name by MHC class/chain using mhcgnomes.

    Returns a dict with keys ``mhc_class``, ``chain``, ``non_mhc``, ``source``,
    or None when the gene cannot be resolved safely.

    This is the lenient classification path that recognizes IPD-MHC
    suffixes like F10, BLB, Q9, E-S, DRA, DAB, etc. without requiring
    a full allele parse.
    """
    if not gene:
        return None
    token = str(gene).strip()
    if not token:
        return None
    kwargs: dict[str, Any] = {"species": species} if species else {"default_species": None}
    result = mhcgnomes.parse_gene_class(token, raise_on_error=False, **kwargs)
    if result is None or is_generated_short_prefix(token, result):
        return None
    mhc_class = str(result.mhc_class or "")
    chain_val = result.chain
    non_mhc = result.non_mhc
    # Normalize: mhcgnomes returns "I", "Ia", "Ib", "IIa", "IIb", "other"
    # We need "I" or "II" for dispatch.
    if mhc_class.startswith("II"):
        mhc_class = "II"
    elif mhc_class.startswith("I"):
        mhc_class = "I"
    elif mhc_class == "other":
        pass  # keep "other" for non-MHC detection
    else:
        mhc_class = None
    return {
        "mhc_class": mhc_class,
        "chain": str(chain_val) if chain_val else None,
        "non_mhc": bool(non_mhc),
        "source": str(result.source or ""),
    }


def is_non_mhc_gene(gene: Optional[str], *, species: Optional[str] = None) -> bool:
    """Check if a gene name is a known non-MHC gene in the MHC region.

    Returns True for genes like TAP1, TAP2, CIITA, HM13, PRR3 that
    appear in MHC-region datasets but don't encode groove proteins.
    Uses species-aware mhcgnomes classification, with the local set for
    unscoped helper-gene labels that do not carry a species identity.
    """
    if not gene:
        return False
    from .domain_grammar import NON_MHC_GENE_NAMES

    token = gene.strip()
    candidates = {token}
    if not token.startswith("~") and "-" in token:
        candidates.add(token.split("-", 1)[1])
    known_non_mhc = {g.upper() for g in NON_MHC_GENE_NAMES}
    if any(candidate.upper() in known_non_mhc for candidate in candidates):
        return True
    result = parse_gene_class(gene, species=species)
    if result is not None:
        return bool(result.get("non_mhc", False))
    return False


def infer_species_identity(allele: Optional[str]) -> Optional[str]:
    """Infer fine-grained species identity from an allele via mhcgnomes."""
    if not allele:
        return None
    parsed = parse_allele_name(allele)
    if parsed is None or parsed.species is None:
        return None
    return str(parsed.species.name).strip() or None


def infer_species(allele: str) -> Optional[str]:
    """Infer 9-class MHC species bucket from allele name."""
    species_identity = infer_species_identity(allele)
    if species_identity is None:
        return None
    return normalize_mhc_species(species_identity)


def allele_suffix_flags(allele: Any) -> dict[str, bool]:
    """Detect null, questionable, and pseudogene flags via mhcgnomes.

    Accepts either an already-parsed mhcgnomes result object or an allele
    name string. Delegates both expression suffixes (``N``, ``Q``, ``Ps``)
    and gene-level pseudogene status (for example ``HLA-H``, ``HLA-J``,
    ``HLA-DRB9``) to mhcgnomes. Haplotype designations such as ``Rano-A1*n``
    and ``H2-D*b`` are not treated as expression suffixes.
    """
    if allele is None:
        return {"is_null": False, "is_questionable": False, "is_pseudogene": False}
    parsed = allele if isinstance(allele, MhcgnomesResult) else parse_allele_name(str(allele))
    if parsed is None:
        return {"is_null": False, "is_questionable": False, "is_pseudogene": False}
    is_pseudogene = bool(getattr(parsed, "is_pseudogene", False) or getattr(parsed, "annotation_pseudogene", False))
    return {
        "is_null": bool(parsed.annotation_null),
        "is_questionable": bool(parsed.annotation_questionable),
        "is_pseudogene": is_pseudogene,
    }
