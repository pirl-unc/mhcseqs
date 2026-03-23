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
    "IA",
    "IB",
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
    short_match = re.match(r"^(?:HLA-)?([A-Z]+)(\d)$", upper)
    if short_match:
        gene, field = short_match.groups()
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
    """Infer MHC class ("I" or "II") from allele name via mhcgnomes."""
    if not allele:
        return None
    try:
        parsed = parse_allele_name(allele)
    except Exception:
        return None
    if parsed is None:
        return None
    return normalize_mhc_class(getattr(parsed, "mhc_class", None), default=None)


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
    """Detect null/questionable/pseudogene suffix markers."""
    token = str(allele or "").strip()
    suffix_match = re.search(r"([A-Za-z]+)$", token)
    suffix = suffix_match.group(1).upper() if suffix_match else ""
    return {
        "is_null": suffix == "N",
        "is_questionable": suffix == "Q",
        "is_pseudogene": suffix == "PS",
    }
