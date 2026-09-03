#!/usr/bin/env python3
"""Recover source labels for the historical prefix audit gap.

The full-binomial migration in commit 78051fb identified 21 legacy
``(species, prefix)`` pairs whose source rows were absent from the retained
raw UniProt download.  This script reconstructs the affected accession list
from the pinned pre-migration corpus, fetches current records from UniProtKB,
falls back to a pinned UniSave entry version for deleted records, and writes
a versioned provenance snapshot.

By default this is read-only with respect to the packaged corpus. Pass
``--update-corpus`` to populate ``raw_gene_label`` and update the pair-level
audit from the recovered snapshot.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import subprocess
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parent.parent
LEGACY_CORPUS_REF = "78051fb^:mhcseqs/diverse_mhc_sequences.csv"
DEFAULT_OUTPUT = ROOT / "data" / "historical_prefix_provenance_20260831.csv"
DEFAULT_CORPUS = ROOT / "mhcseqs" / "diverse_mhc_sequences.csv"
DEFAULT_AUDIT = ROOT / "data" / "species_prefix_audit_pre_full_binomial.csv"
USER_AGENT = "mhcseqs issue-49 provenance recovery (https://github.com/pirl-unc/mhcseqs/issues/49)"

SOURCE_PREFIX = "source_attested_prefix"
GENERIC_LABEL = "generic_label_not_prefix"
FOREIGN_NOMENCLATURE = "foreign_nomenclature_not_prefix"

SNAPSHOT_FIELDS = [
    "uniprot_accession",
    "species",
    "legacy_prefix",
    "legacy_gene",
    "source_gene_names",
    "source_protein_name",
    "label_classification",
    "primary_reference",
    "source_kind",
    "source_url",
    "uniprot_entry_version",
    "uniprot_last_modified",
    "uniprot_release",
    "uniprot_release_date",
    "retrieved_on",
]


@dataclass(frozen=True)
class PairDecision:
    classification: str
    primary_reference: str = ""


# These decisions classify the exact source token, not the canonical species
# identity. Generic source labels and a mouse-system synonym on a rat record
# are preserved as provenance but deliberately excluded from the alias table.
PAIR_DECISIONS = {
    ("Arvicola amphibius", "Arte"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/17956550/"),
    ("Cryptomys hottentotus", "BLA"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/15058438/"),
    ("Ctenomys boliviensis", "Ctbol"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/18049818/"),
    ("Ctenomys bonettoi", "Ctbon"): PairDecision(SOURCE_PREFIX),
    ("Ctenomys talarum", "Cta"): PairDecision(SOURCE_PREFIX),
    ("Fukomys damarensis", "BLA"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/15058438/"),
    ("Heliophobius argenteocinereus", "BLA"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/15058438/"),
    ("Heterocephalus glaber", "BLA"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/15058438/"),
    ("Mus caroli", "H2"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/9178014/"),
    ("Mus cervicolor", "H2"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/9178014/"),
    ("Mus cookii", "H2"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/9178014/"),
    ("Mus macedonicus", "H2"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/9178014/"),
    ("Mus musculus", "I"): PairDecision(GENERIC_LABEL),
    ("Mus musculus", "MHC II H2"): PairDecision(GENERIC_LABEL, "https://pubmed.ncbi.nlm.nih.gov/3018929/"),
    ("Mus platythrix", "H2"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/9178014/"),
    ("Mus spicilegus", "H2"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/9178014/"),
    ("Mus spretus", "H2"): PairDecision(SOURCE_PREFIX, "https://pubmed.ncbi.nlm.nih.gov/9178014/"),
    ("Myodes glareolus", "Clgl"): PairDecision(SOURCE_PREFIX),
    ("Peromyscus californicus", "MHC"): PairDecision(GENERIC_LABEL, "https://pubmed.ncbi.nlm.nih.gov/22649541/"),
    ("Peromyscus maniculatus", "MHC"): PairDecision(GENERIC_LABEL, "https://pubmed.ncbi.nlm.nih.gov/22649541/"),
    ("Rattus norvegicus", "H2"): PairDecision(FOREIGN_NOMENCLATURE),
}


def latin_binomial(organism: str) -> str:
    """Return the first two scientific-name words from a source label."""
    return " ".join(organism.split(" (", 1)[0].split()[:2])


def read_legacy_corpus(ref: str = LEGACY_CORPUS_REF) -> list[dict[str, str]]:
    """Read the pinned pre-migration corpus from git history."""
    result = subprocess.run(
        ["git", "show", ref],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return list(csv.DictReader(io.StringIO(result.stdout)))


def select_affected_rows(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    """Select the 239 records represented by the 21 pair-level gaps."""
    selected = []
    for row in rows:
        gene = row["gene"]
        for species, prefix in PAIR_DECISIONS:
            if latin_binomial(row["organism"]) != species:
                continue
            if gene == prefix or gene.startswith(prefix + "-"):
                selected.append(row)
                break
    if len(selected) != 239:
        raise ValueError(f"Expected 239 historical records, found {len(selected)}")
    if len({row["uniprot_accession"] for row in selected}) != len(selected):
        raise ValueError("Historical provenance records must have unique accessions")
    return selected


def _open_url(url: str):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(request, timeout=180)


def _chunks(values: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def fetch_current_records(accessions: list[str]) -> tuple[dict[str, dict[str, str]], str, str]:
    """Fetch current UniProtKB records and current release metadata."""
    records: dict[str, dict[str, str]] = {}
    releases = set()
    release_dates = set()
    fields = "accession,protein_name,gene_names,organism_name,organism_id,date_modified,version"
    for chunk in _chunks(accessions, 100):
        query = urllib.parse.urlencode(
            {
                "accessions": ",".join(chunk),
                "format": "tsv",
                "fields": fields,
            }
        )
        with _open_url(f"https://rest.uniprot.org/uniprotkb/accessions?{query}") as response:
            releases.add(response.headers["X-UniProt-Release"])
            release_dates.add(response.headers["X-UniProt-Release-Date"])
            text = response.read().decode("utf-8")
        records.update({row["Entry"]: row for row in csv.DictReader(io.StringIO(text), delimiter="\t")})
    if len(releases) != 1 or len(release_dates) != 1:
        raise ValueError(f"UniProt release changed during retrieval: {releases}, {release_dates}")
    return records, releases.pop(), release_dates.pop()


def _flat_file_values(text: str, code: str) -> str:
    return " ".join(line[5:].strip() for line in text.splitlines() if line.startswith(code + "   "))


def _strip_evidence(value: str) -> str:
    return re.sub(r"\s*\{ECO:[^}]+\}", "", value).strip()


def parse_unisave_record(text: str) -> dict[str, str]:
    """Parse fields needed from one pinned UniSave flat-file version."""
    accession_match = re.search(r"^AC\s+([^;]+);", text, re.MULTILINE)
    entry_match = re.search(r"^DT\s+(\d{2}-[A-Z]{3}-\d{4}), entry version (\d+)\.", text, re.MULTILINE)
    if accession_match is None or entry_match is None:
        raise ValueError("UniSave record is missing accession or entry-version metadata")

    gene_text = _flat_file_values(text, "GN")
    names = []
    for value in re.findall(r"(?:Name|Synonyms|OrderedLocusNames|ORFNames)=([^;]+)", gene_text):
        names.extend(_strip_evidence(value).replace(",", " ").split())

    protein_text = _flat_file_values(text, "DE")
    protein_match = re.search(r"Full=([^;{]+)", protein_text)
    return {
        "Entry": accession_match.group(1),
        "Gene Names": " ".join(names),
        "Protein names": protein_match.group(1).strip() if protein_match else "",
        "Date of last modification": entry_match.group(1),
        "Entry version": entry_match.group(2),
    }


def fetch_archived_record(accession: str) -> dict[str, str]:
    """Fetch the latest archived version of a deleted UniProt accession."""
    metadata_url = f"https://rest.uniprot.org/unisave/{accession}?format=json"
    with _open_url(metadata_url) as response:
        versions = json.loads(response.read().decode("utf-8"))["results"]
    if not versions:
        raise ValueError(f"UniSave contains no versions for {accession}")
    latest = max(versions, key=lambda row: int(row["entryVersion"]))
    version = str(latest["entryVersion"])
    query = urllib.parse.urlencode({"format": "txt", "versions": version})
    source_url = f"https://rest.uniprot.org/unisave/{accession}?{query}"
    with _open_url(source_url) as response:
        parsed = parse_unisave_record(response.read().decode("utf-8"))
    parsed.update(
        {
            "source_url": source_url,
            "uniprot_release": latest["lastRelease"].split("/", 1)[0],
            "uniprot_release_date": latest["lastReleaseDate"],
        }
    )
    return parsed


def fetch_source_records(accessions: list[str]) -> dict[str, dict[str, str]]:
    """Fetch current records, then recover deleted accessions from UniSave."""
    current, release, release_date = fetch_current_records(accessions)
    result = {}
    for accession, row in current.items():
        result[accession] = {
            **row,
            "source_url": f"https://rest.uniprot.org/uniprotkb/{accession}",
            "uniprot_release": release,
            "uniprot_release_date": release_date,
            "source_kind": "uniprotkb_current",
        }

    missing = sorted(set(accessions) - set(current))
    with ThreadPoolExecutor(max_workers=8) as pool:
        archived = dict(zip(missing, pool.map(fetch_archived_record, missing)))
    for accession, row in archived.items():
        result[accession] = {**row, "source_kind": "unisave_archived"}
    return result


def build_snapshot_rows(
    legacy_rows: Iterable[dict[str, str]],
    source_records: dict[str, dict[str, str]],
    retrieved_on: str,
) -> list[dict[str, str]]:
    """Join historical curation rows to pinned source records."""
    result = []
    for legacy in legacy_rows:
        accession = legacy["uniprot_accession"]
        source = source_records.get(accession)
        if source is None:
            raise ValueError(f"No current or archived UniProt record for {accession}")
        species = latin_binomial(legacy["organism"])
        prefix = legacy["gene"].split("-", 1)[0]
        decision = PAIR_DECISIONS[(species, prefix)]
        gene_names = source["Gene Names"]
        if not gene_names:
            raise ValueError(f"Recovered source record {accession} has no gene label")
        result.append(
            {
                "uniprot_accession": accession,
                "species": species,
                "legacy_prefix": prefix,
                "legacy_gene": legacy["gene"],
                "source_gene_names": gene_names,
                "source_protein_name": source["Protein names"],
                "label_classification": decision.classification,
                "primary_reference": decision.primary_reference,
                "source_kind": source["source_kind"],
                "source_url": source["source_url"],
                "uniprot_entry_version": source["Entry version"],
                "uniprot_last_modified": source["Date of last modification"],
                "uniprot_release": source["uniprot_release"],
                "uniprot_release_date": source["uniprot_release_date"],
                "retrieved_on": retrieved_on,
            }
        )
    return sorted(result, key=lambda row: (row["species"], row["legacy_prefix"], row["uniprot_accession"]))


def write_rows(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def update_corpus(corpus_path: Path, snapshot_rows: list[dict[str, str]]) -> None:
    """Populate raw source labels for every recovered packaged accession."""
    by_accession = {row["uniprot_accession"]: row["source_gene_names"] for row in snapshot_rows}
    with corpus_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = list(reader)
    if fields is None or "raw_gene_label" not in fields:
        raise ValueError(f"{corpus_path} has no raw_gene_label column")
    matched = 0
    for row in rows:
        label = by_accession.get(row["uniprot_accession"])
        if label is not None:
            row["raw_gene_label"] = label
            matched += 1
    if matched != len(snapshot_rows):
        raise ValueError(f"Expected to update {len(snapshot_rows)} corpus rows, updated {matched}")
    write_rows(corpus_path, rows, fields)


def update_audit(audit_path: Path, snapshot_rows: list[dict[str, str]]) -> None:
    """Resolve the 21 pair-level audit gaps from the accession snapshot."""
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in snapshot_rows:
        grouped.setdefault((row["species"], row["legacy_prefix"]), []).append(row)

    with audit_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = list(reader)
    if fields is None:
        raise ValueError(f"{audit_path} has no CSV header")

    resolved = 0
    for row in rows:
        pair = (row["species"], row["legacy_prefix"])
        records = grouped.get(pair)
        if records is None:
            continue
        decision = PAIR_DECISIONS[pair]
        first = records[0]
        if decision.classification == SOURCE_PREFIX:
            row["status"] = "external_uniprot_label"
            row["source_attested_row_count"] = row["row_count"]
            row["evidence"] = decision.primary_reference or first["source_url"]
        elif decision.classification == GENERIC_LABEL:
            row["status"] = "source_label_artifact"
            row["source_attested_row_count"] = "0"
            row["evidence"] = first["source_url"] + "; generic gene/product token, not a species prefix"
        else:
            row["status"] = "source_foreign_nomenclature_alias"
            row["source_attested_row_count"] = "0"
            row["evidence"] = first["source_url"] + "; H2 is a mouse-system synonym of the primary rat RT1 gene"
        row["example_accessions"] = ";".join(record["uniprot_accession"] for record in records[:3])
        labels = []
        for record in records:
            label = " ".join((record["source_gene_names"], record["source_protein_name"])).strip()
            if label not in labels:
                labels.append(label)
        row["example_source_labels"] = ";".join(labels[:3])
        resolved += 1
    if resolved != len(PAIR_DECISIONS):
        raise ValueError(f"Expected to resolve {len(PAIR_DECISIONS)} audit pairs, resolved {resolved}")
    write_rows(audit_path, rows, fields)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-ref", default=LEGACY_CORPUS_REF)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--retrieved-on", default=date.today().isoformat())
    parser.add_argument("--update-corpus", action="store_true")
    args = parser.parse_args()

    legacy_rows = select_affected_rows(read_legacy_corpus(args.legacy_ref))
    accessions = sorted(row["uniprot_accession"] for row in legacy_rows)
    print(f"Recovering {len(accessions)} records for {len(PAIR_DECISIONS)} historical prefix pairs")
    source_records = fetch_source_records(accessions)
    snapshot_rows = build_snapshot_rows(legacy_rows, source_records, args.retrieved_on)
    write_rows(args.output, snapshot_rows, SNAPSHOT_FIELDS)
    print(f"Wrote {len(snapshot_rows)} source records to {args.output}")
    if args.update_corpus:
        update_corpus(args.corpus, snapshot_rows)
        update_audit(args.audit, snapshot_rows)
        print(f"Updated {args.corpus} and {args.audit}")


if __name__ == "__main__":
    main()
