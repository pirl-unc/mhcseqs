"""Allele name parsing and normalization via mhcgnomes.

Provides functions to parse MHC allele names from FASTA headers, normalize them
to canonical two-field resolution, and infer gene/class/species metadata.

Ported from presto/data/allele_resolver.py.
"""

from __future__ import annotations

import importlib
import re
from typing import Any, Optional

from .species import normalize_mhc_species

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
    short_match = re.match(r"^(?:HLA-)?([A-Z]+)(\d)$", upper)
    if short_match:
        gene, field = short_match.groups()
        if gene in _HLA_SHORT_GENES:
            return f"HLA-{gene}*0{field}"
    return token


def parse_allele_name(allele: Optional[str], *, species: Optional[str] = None) -> Optional[Any]:
    """Parse an allele string with mhcgnomes.

    Parameters
    ----------
    allele : str
        Allele name to parse (e.g. "HLA-A*02:01", "Crpo-UA", "H2-Kb").
    species : str, optional
        Latin binomial of the species (e.g. "Homo sapiens", "Crocodylus porosus").
        When provided, passed to mhcgnomes as the ``species`` parameter so the
        parser can validate and disambiguate the allele for that organism.
    """
    if not allele:
        return None
    mhcgnomes = _require_mhcgnomes()
    parse_fn = mhcgnomes.parse
    raw = str(allele).strip()
    if not raw:
        return None

    # Build kwargs: use species= if mhcgnomes supports it, else default_species
    import inspect

    parse_params = inspect.signature(parse_fn).parameters
    kwargs: dict = {}
    if species:
        if "species" in parse_params:
            kwargs["species"] = species
        elif "default_species" in parse_params:
            kwargs["default_species"] = species

    coerced = _coerce_allele_name(raw)
    candidates = []
    if coerced:
        candidates.append(coerced)
    if raw not in candidates:
        candidates.append(raw)

    for candidate in candidates:
        try:
            parsed = parse_fn(candidate, **kwargs)
        except Exception:
            parsed = None
        if parsed is not None:
            return parsed

    # Fallback: try without species constraint (for IMGT/IPD alleles
    # that already embed the species prefix in the name)
    if kwargs:
        for candidate in candidates:
            try:
                parsed = parse_fn(candidate)
            except Exception:
                parsed = None
            if parsed is not None:
                return parsed

    return None


def _canonicalize_parsed_allele(parsed: Any, allele_fields: int = 2) -> str:
    target_fields = max(1, int(allele_fields))
    restrict_fn = getattr(parsed, "restrict_allele_fields", None)
    if callable(restrict_fn):
        parsed = restrict_fn(target_fields)
    to_string = getattr(parsed, "to_string", None)
    if callable(to_string):
        return str(to_string())
    return str(parsed)


def normalize_allele_name(name: str) -> str:
    """Normalize an allele name to canonical two-field protein resolution.

    Examples:
        "HLA-A*02:01" -> "HLA-A*02:01"
        "A*02:01" -> "HLA-A*02:01"
        "A0201" -> "HLA-A*02:01"
        "HLA-A*02:01:01:02L" -> "HLA-A*02:01L"
    """
    parsed = parse_allele_name(name)
    if parsed is None:
        raise ValueError(f"mhcgnomes failed to parse allele: {name!r}")
    return _canonicalize_parsed_allele(parsed, allele_fields=2)


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


def infer_mhc_class(allele: Optional[str]) -> Optional[str]:
    """Infer MHC class ("I" or "II") from allele name via mhcgnomes.

    Uses ``parse_gene_class()`` first (lenient suffix-based classification
    available in mhcgnomes >= 3.18), then falls back to full allele parsing.
    """
    if not allele:
        return None
    # Try lenient gene-class inference first (mhcgnomes >= 3.18)
    result = parse_gene_class(allele)
    if result is not None:
        cls = result.get("mhc_class")
        if cls in ("I", "II"):
            return cls
    # Fall back to full allele parsing
    try:
        parsed = parse_allele_name(allele)
    except Exception:
        return None
    if parsed is None:
        return None
    return normalize_mhc_class(getattr(parsed, "mhc_class", None), default=None)


def parse_gene_class(gene: Optional[str]) -> Optional[dict]:
    """Classify a gene name by MHC class/chain using mhcgnomes (>= 3.18).

    Returns a dict with keys ``mhc_class``, ``chain``, ``non_mhc``,
    or None if the function is not available (mhcgnomes < 3.18).

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
        result = fn(str(gene).strip())
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
        }
    except Exception:
        pass
    return None


def is_non_mhc_gene(gene: Optional[str]) -> bool:
    """Check if a gene name is a known non-MHC gene in the MHC region.

    Returns True for genes like TAP1, TAP2, CIITA, HM13, PRR3 that
    appear in MHC-region datasets but don't encode groove proteins.
    Uses mhcgnomes >= 3.18 when available, falls back to a local set.
    """
    if not gene:
        return False
    result = parse_gene_class(gene)
    if result is not None:
        return bool(result.get("non_mhc", False))
    # Fallback for mhcgnomes < 3.18
    from .domain_grammar import NON_MHC_GENE_NAMES

    return gene.strip() in NON_MHC_GENE_NAMES or gene.strip().upper() in {g.upper() for g in NON_MHC_GENE_NAMES}


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
