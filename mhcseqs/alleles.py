"""Allele name parsing and normalization via mhcgnomes.

Provides functions to parse MHC allele names from FASTA headers, normalize them
to canonical two-field resolution, and infer gene/class/species metadata.

Ported from presto/data/allele_resolver.py.
"""

from __future__ import annotations

import importlib
import re
from typing import Any, Optional

from .species import (
    extract_latin_binomial,
    find_mhc_prefix_aliases,
    full_species_name_alias,
    normalize_mhc_species,
)

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
# mhcgnomes wrapper
# ---------------------------------------------------------------------------


def _require_mhcgnomes() -> Any:
    try:
        mhcgnomes = importlib.import_module("mhcgnomes")
    except ImportError as exc:
        raise RuntimeError("mhcgnomes is required: pip install mhcgnomes") from exc
    if callable(getattr(mhcgnomes, "parse", None)):
        return mhcgnomes
    try:
        function_api = importlib.import_module("mhcgnomes.function_api")
    except ImportError as exc:
        raise RuntimeError("mhcgnomes is installed but its parse API could not be imported.") from exc
    parse_fn = getattr(function_api, "parse", None)
    if not callable(parse_fn):
        raise RuntimeError("mhcgnomes is installed but does not expose a callable parse API.")
    setattr(mhcgnomes, "parse", parse_fn)
    return mhcgnomes


def _coerce_allele_name(allele: Optional[str]) -> str:
    token = str(allele or "").strip().strip(",;")
    if not token:
        return ""
    token = token.replace("_", "-")
    upper = token.upper()
    if upper.startswith("H-2-"):
        token = "H2-" + token[4:]
        upper = token.upper()
    if upper.startswith("H-2"):
        token = "H2-" + token[3:]
        upper = token.upper()
    _HLA_SHORT_GENES = {"A", "B", "C", "E", "F", "G"}
    short_match = re.match(r"^(?:HLA-)?([A-Z]+)(\d+)$", upper)
    if short_match:
        gene, digits = short_match.groups()
        if gene in _HLA_SHORT_GENES:
            allele_digits = digits.zfill(2) if len(digits) == 1 else digits
            return f"HLA-{gene}*{allele_digits}"
    return token


def _molecule_result_types(mhcgnomes: Any) -> tuple[type, ...]:
    """Result types that denote an MHC molecule rather than a loose label."""
    return tuple(
        result_type for name in ("Allele", "Gene", "AlleleWithoutGene", "Pair") if isinstance((result_type := getattr(mhcgnomes, name, None)), type)
    )


def _source_alias_candidates(
    value: str,
    *,
    species: Optional[str],
    mhcgnomes: Any,
) -> tuple[list[str], bool]:
    """Rewrite a source-attested prefix to an unambiguous full taxon name.

    The boolean return value marks a genuinely ambiguous prefix. Such a name
    must not be handed to mhcgnomes without species context because its global
    alias table may otherwise choose one valid claimant silently.
    """
    _, body, entries = find_mhc_prefix_aliases(value, species=species)
    if not body:
        return [value], False

    if species:
        if not entries:
            return [value], False
        target = mhcgnomes.Species.get(species)
        target_name = getattr(target, "name", None) or species
        return [f"{full_species_name_alias(target_name)}-{body}"], False

    # A formal genus/group system such as DLA or BoLA is intentionally an
    # umbrella rather than an ambiguous species claim. Let mhcgnomes resolve
    # those ontology-backed group aliases directly.
    if any(entry.scope == "genus" and entry.status == "formal_system" for entry in entries):
        return [value], False

    binomials = {extract_latin_binomial(entry.species).casefold() for entry in entries}
    if len(binomials) != 1:
        return [], True

    target = mhcgnomes.Species.get(entries[0].species)
    if target is None:
        return [value], False
    return [f"{full_species_name_alias(target.name)}-{body}"], False


def _is_unattested_generated_prefix(candidate: str, parsed: Any) -> bool:
    """Return whether an explicit short prefix is only mechanically derived."""
    species = getattr(parsed, "species", None)
    species_name = getattr(species, "name", None)
    if not species_name or "-" not in candidate:
        return False

    prefix = candidate.split("-", 1)[0]
    normalized_prefix = re.sub(r"[^A-Za-z0-9]+", "", prefix).casefold()
    if normalized_prefix == full_species_name_alias(species_name).casefold():
        return False

    _, body, entries = find_mhc_prefix_aliases(candidate, species=species_name)
    if body and entries:
        return False

    words = re.findall(r"[A-Za-z]+", extract_latin_binomial(species_name))
    if len(words) < 2:
        return False
    generated = {(words[0][:width] + words[1][:width]).casefold() for width in (2, 4, 5)}
    return normalized_prefix in generated


def parse_allele_name(
    allele: Optional[str],
    *,
    species: Optional[str] = None,
    require_explicit_species: bool = False,
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
    """
    if not allele:
        return None
    mhcgnomes = _require_mhcgnomes()
    parse_fn = mhcgnomes.parse
    raw = str(allele).strip()
    if not raw:
        return None

    kwargs: dict[str, Any] = {
        "required_result_types": _molecule_result_types(mhcgnomes),
        "raise_on_error": False,
        "require_explicit_species": require_explicit_species,
    }
    if species:
        kwargs["species"] = species
    elif require_explicit_species:
        # Keep strict provenance available without changing the ordinary
        # mhcgnomes behavior where a bare allele such as A*02:01 means HLA.
        kwargs["default_species"] = None

    coerced = _coerce_allele_name(raw)
    initial_candidates = []
    if coerced:
        initial_candidates.append(coerced)
    if raw not in initial_candidates:
        initial_candidates.append(raw)

    candidates = []
    saw_ambiguous_alias = False
    for initial in initial_candidates:
        source_candidates, ambiguous = _source_alias_candidates(
            initial,
            species=species,
            mhcgnomes=mhcgnomes,
        )
        if ambiguous:
            # Coercion can turn an unambiguous historical spelling into a
            # colliding modern one (H-2 -> H2). Keep trying the original
            # source spelling; reject only when no unambiguous candidate is
            # available.
            saw_ambiguous_alias = True
            continue
        for candidate in source_candidates:
            if candidate not in candidates:
                candidates.append(candidate)

    if not candidates and saw_ambiguous_alias:
        return None

    for candidate in candidates:
        try:
            parsed = parse_fn(candidate, **kwargs)
        except Exception:
            parsed = None
        if parsed is not None:
            if _is_unattested_generated_prefix(candidate, parsed):
                continue
            return parsed

    return None


def _canonicalize_parsed_name(parsed: Any, allele_fields: Optional[int] = None) -> str:
    if allele_fields is not None:
        target_fields = max(1, int(allele_fields))
        restrict_fn = getattr(parsed, "restrict_allele_fields", None)
        if callable(restrict_fn):
            parsed = restrict_fn(target_fields)

    to_string = getattr(parsed, "to_string", None)
    species = getattr(parsed, "species", None)
    species_name = getattr(species, "name", None)
    species_alias = full_species_name_alias(species_name)
    if callable(to_string) and species_alias:
        body = str(to_string(include_species=False))
        if body:
            return f"{species_alias}-{body}"
    if callable(to_string):
        return str(to_string())
    return str(parsed)


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
    return _canonicalize_parsed_name(parsed, allele_fields=2)


def infer_gene(allele: str) -> str:
    """Extract gene name from allele (e.g. "A", "DRB1")."""
    try:
        parsed = parse_allele_name(allele)
    except Exception:
        parsed = None
    if parsed is not None and getattr(parsed, "gene", None) is not None:
        gene_name = getattr(parsed.gene, "name", None)
        if gene_name:
            return str(gene_name).upper()

    token = str(allele).strip()
    if "-" in token:
        _, remainder = token.split("-", 1)
        if remainder:
            token = remainder
    if "*" in token:
        return token.split("*")[0]
    return token


def infer_mhc_class(allele: Optional[str], *, species: Optional[str] = None) -> Optional[str]:
    """Infer MHC class ("I" or "II") from allele name via mhcgnomes.

    Uses ``parse_gene_class()`` first (species-aware classification in
    mhcgnomes >= 3.40), then falls back to full allele parsing.
    """
    if not allele:
        return None
    # Try species-aware gene-class inference first.
    result = parse_gene_class(allele, species=species)
    if result is not None:
        cls = result.get("mhc_class")
        if cls in ("I", "II"):
            return cls
    # Fall back to full allele parsing
    try:
        parsed = parse_allele_name(allele, species=species)
    except Exception:
        return None
    if parsed is None:
        return None
    return normalize_mhc_class(getattr(parsed, "mhc_class", None), default=None)


def parse_gene_class(gene: Optional[str], *, species: Optional[str] = None) -> Optional[dict]:
    """Classify a gene name by MHC class/chain using mhcgnomes (>= 3.40).

    Returns a dict with keys ``mhc_class``, ``chain``, ``non_mhc``,
    or None when the gene cannot be resolved safely.

    This is the lenient classification path that recognizes IPD-MHC
    suffixes like F10, BLB, Q9, E-S, DRA, DAB, etc. without requiring
    a full allele parse.
    """
    if not gene:
        return None
    try:
        mhcgnomes = _require_mhcgnomes()
        fn = getattr(mhcgnomes, "parse_gene_class", None)
        if fn is None:
            return None
        kwargs = {"species": species} if species else {"default_species": None}
        candidates, ambiguous = _source_alias_candidates(
            str(gene).strip(),
            species=species,
            mhcgnomes=mhcgnomes,
        )
        if ambiguous:
            return None
        result = None
        for candidate in candidates:
            result = fn(candidate, raise_on_error=False, **kwargs)
            if result is not None and not _is_unattested_generated_prefix(candidate, result):
                break
            result = None
        if result is None:
            return None
        # Result is a GeneClassInfo dataclass — extract fields
        mhc_class = str(getattr(result, "mhc_class", "") or "")
        chain_val = getattr(result, "chain", None)
        non_mhc = getattr(result, "non_mhc", False)
        # Normalize: mhcgnomes returns "I", "Ib", "IIa", "IIb", "other"
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
            "source": str(getattr(result, "source", "") or ""),
        }
    except Exception:
        pass
    return None


def is_non_mhc_gene(gene: Optional[str], *, species: Optional[str] = None) -> bool:
    """Check if a gene name is a known non-MHC gene in the MHC region.

    Returns True for genes like TAP1, TAP2, CIITA, HM13, PRR3 that
    appear in MHC-region datasets but don't encode groove proteins.
    Uses species-aware mhcgnomes 3.40 classification, with the local set only
    for unscoped helper-gene labels that do not carry a species identity.
    """
    if not gene:
        return False
    # Local curation is an override, not merely a compatibility fallback.
    # mhcgnomes 3.40 can otherwise split names such as Kdm5d and Daxx into
    # syntactically valid mouse K/D alleles (tracked upstream in #133).
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
    try:
        parsed = parse_allele_name(allele)
    except Exception:
        return None
    if parsed is None or getattr(parsed, "species", None) is None:
        return None
    species_name = getattr(parsed.species, "name", None)
    return str(species_name).strip() or None


def infer_species(allele: str) -> Optional[str]:
    """Infer 7-class MHC species bucket from allele name."""
    species_identity = infer_species_identity(allele)
    if species_identity is None:
        return None
    return normalize_mhc_species(species_identity)


def allele_suffix_flags(allele: str) -> dict[str, bool]:
    """Detect null/questionable/pseudogene suffix markers.

    HLA-style suffixes (N, Q, L, S) appear after numeric allele fields,
    e.g. ``HLA-A*02:01N``.  Haplotype-based systems like Rano (``Rano-A1*n``)
    and H-2 (``H2-D*q``) use single letters as the allele designation
    itself — these must NOT be treated as null/questionable markers.
    """
    token = str(allele or "").strip()
    # If there is a '*', check whether the designation after it contains
    # digits.  Pure-letter designations are haplotype names, not suffixes.
    if "*" in token:
        after_star = token.rsplit("*", 1)[1]
        if not re.search(r"\d", after_star):
            return {"is_null": False, "is_questionable": False, "is_pseudogene": False}
    suffix_match = re.search(r"([A-Za-z]+)$", token)
    suffix = suffix_match.group(1).upper() if suffix_match else ""
    return {
        "is_null": suffix == "N",
        "is_questionable": suffix == "Q",
        "is_pseudogene": suffix == "PS",
    }
