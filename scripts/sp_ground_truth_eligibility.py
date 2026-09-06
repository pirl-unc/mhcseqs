"""Shared MHC-chain eligibility policy for the SP benchmark corpus."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from mhcseqs.alleles import is_non_mhc_gene, parse_allele_name, parse_gene_class
from scripts.curate_diverse_mhc import _infer_class_ii_chain, classify_mhc

ROOT = Path(__file__).resolve().parent.parent
GT_LABEL_CURATION_CSV = ROOT / "data" / "sp_ground_truth_label_curation.csv"
ELIGIBLE_MHC_CHAINS = frozenset({("I", "alpha"), ("II", "alpha"), ("II", "beta")})


@dataclass(frozen=True)
class MhcLabel:
    """One source-backed class/chain decision for benchmark eligibility."""

    mhc_class: str
    chain: str
    label_status: str
    disposition: str

    @property
    def eligible(self) -> bool:
        """Return whether this is a usable MHC alpha/beta chain."""
        return (self.mhc_class, self.chain) in ELIGIBLE_MHC_CHAINS and self.disposition == "include"


def load_label_curation(path: Path = GT_LABEL_CURATION_CSV) -> dict[str, dict[str, str]]:
    """Load accession-level class/chain decisions and their provenance."""
    if not path.exists():
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    result = {row["accession"]: row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate accessions in {path}")
    return result


def resolve_mhc_label(
    row: dict[str, str],
    curation: dict[str, dict[str, str]],
) -> MhcLabel:
    """Resolve artifact metadata through the shared benchmark label policy."""
    accession = row.get("Entry", "").strip()
    decision = curation.get(accession)
    if decision is not None:
        disposition = decision.get("disposition", "").strip()
        mhc_class = decision.get("mhc_class", "").strip().upper()
        chain = decision.get("chain", "").strip().lower()
        label_status = decision.get("label_status", "").strip()
        if disposition == "include":
            if (mhc_class, chain) not in ELIGIBLE_MHC_CHAINS:
                raise ValueError(f"Invalid curated MHC label for {accession}: {(mhc_class, chain)!r}")
            return MhcLabel(mhc_class, chain, label_status or "curated", disposition)
        if disposition == "exclude_non_mhc":
            return MhcLabel("", "", label_status or "excluded_non_mhc", disposition)
        if disposition == "retain_unresolved":
            return MhcLabel(mhc_class, chain, label_status or "unresolved", disposition)
        raise ValueError(f"Invalid curation disposition for {accession}: {disposition!r}")

    protein_name = row.get("Protein names", "")
    gene_names = row.get("Gene Names", "")
    # Helper genes such as CIITA mention MHC in their names but do not encode
    # an MHC chain. Consult the same species-aware identity policy as parsing
    # before the permissive name heuristics, keeping explicit curation first.
    species = row.get("Organism", "").strip() or None
    gene_tokens = [token for token in re.split(r"[\s,;]+", gene_names) if token]
    non_mhc_tokens = {token for token in gene_tokens if is_non_mhc_gene(token, species=species)}
    if non_mhc_tokens:
        # Some source rows conflate MHC genes and nearby helper genes in one
        # synonym list. Conflicting gene evidence warrants abstention, not an
        # exclusion (or a guessed choice of the first token).
        for token in gene_tokens:
            if token in non_mhc_tokens:
                continue
            gene_class = parse_gene_class(token, species=species)
            if gene_class and not gene_class["non_mhc"] and (gene_class["mhc_class"], gene_class["chain"]) in ELIGIBLE_MHC_CHAINS:
                return MhcLabel("", "", "unresolved", "retain_unresolved")
        return MhcLabel("", "", "excluded_non_mhc", "exclude_non_mhc")
    classified = _classify_from_parsed_names(protein_name, gene_names)
    if classified is None:
        classified = classify_mhc(protein_name, gene_names)
    if classified in ELIGIBLE_MHC_CHAINS:
        mhc_class, chain = classified
        return MhcLabel(mhc_class, chain, "gold", "include")
    if classified is None:
        # Missing or vague metadata is not affirmative non-MHC evidence.
        return MhcLabel("", "", "unresolved", "retain_unresolved")
    mhc_class, chain = classified
    return MhcLabel(mhc_class, chain, "unresolved", "retain_unresolved")


def _classify_from_parsed_names(protein_name: str, gene_names: str) -> tuple[str, str] | None:
    """Use mhcgnomes-parsed source labels before free-text heuristics."""
    gene_tokens = [token for token in gene_names.replace(";", " ").replace(",", " ").split() if "*" in token or "-" in token]
    protein_tokens = [token for token in re.findall(r"[A-Za-z0-9*:.+-]+", protein_name) if "*" in token or "-" in token]
    labels: set[tuple[str, str]] = set()
    for token in (*gene_tokens, *protein_tokens):
        parsed = parse_allele_name(token, require_explicit_species=True)
        if parsed is None:
            continue
        parsed_class = str(getattr(parsed, "mhc_class", "") or "")
        if parsed_class in {"I", "Ia", "Ib"}:
            labels.add(("I", "alpha"))
        if parsed_class.startswith("II"):
            chain = _infer_class_ii_chain(token, "")
            if chain in {"alpha", "beta"}:
                labels.add(("II", chain))
    if len(labels) == 1:
        return labels.pop()
    return None
