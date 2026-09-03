"""Downloaded UniProt artifacts backing the SP benchmark corpus.

The benchmark pipeline used to issue live, parameterized UniProt REST queries
at run time: one capped search per hand-picked clade, taxonomy lookups, and
per-accession enrichment requests. Those queries were unreproducible (the
answer changed with the release), rate-limited, and silently truncated large
clades.

``scripts/fetch_sp_ground_truth.py --refresh`` is now the only networked path.
It saves one broad Vertebrata candidate artifact plus complete UniSave records
for benchmark accessions absent from the current stream. Ordinary runs derive
the benchmark only from those local, release-stamped files.

The current artifact is deliberately broad: no signal-peptide requirement,
length window, or hand-picked query clades. Precision and benchmark eligibility
are local derivation concerns, not download concerns.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
import sys
import tempfile
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

VERTEBRATA_TAXON_ID = 7742

# These mutually exclusive roots partition every candidate in the Vertebrata
# artifact. The last three cover lineages missed by the historical list.
CORPUS_CLADES = (
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

# This preserves the historical candidate definition. It is intentionally
# broad; issue #73 tracks the local distinction between genuine MHC chains and
# name-query contaminants.
MHC_QUERY_CLAUSE = '(keyword:KW-0491 OR protein_name:histocompatibility OR protein_name:"class I" OR protein_name:"class II")'

ARTIFACT_FIELDS = (
    "accession",
    "id",
    "organism_name",
    "organism_id",
    "lineage_ids",
    "lineage",
    "protein_name",
    "gene_names",
    "reviewed",
    "protein_existence",
    "fragment",
    "length",
    "sequence",
    "ft_signal",
    "ft_chain",
    "keyword",
    "date_modified",
    "version",
)

# Column labels emitted by the UniProt TSV API for ARTIFACT_FIELDS. Deleted
# entries use the same schema so downstream derivation does not need a second
# representation.
ARTIFACT_COLUMNS = (
    "Entry",
    "Entry Name",
    "Organism",
    "Organism (ID)",
    "Taxonomic lineage (Ids)",
    "Taxonomic lineage",
    "Protein names",
    "Gene Names",
    "Reviewed",
    "Protein existence",
    "Fragment",
    "Length",
    "Sequence",
    "Signal peptide",
    "Chain",
    "Keywords",
    "Date of last modification",
    "Entry version",
)
DELETED_COLUMNS = ARTIFACT_COLUMNS + (
    "Source clade",
    "Archive release",
    "Archive release date",
    "Archive source URL",
)

STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
UNISAVE_URL = "https://rest.uniprot.org/unisave"
USER_AGENT = "mhcseqs-sp-corpus/1.0"

CORPUS_ARTIFACT = "uniprot_mhc_vertebrata.tsv.gz"
DELETED_ARTIFACT = "uniprot_deleted_entries.tsv.gz"
RELEASE_MANIFEST = "uniprot_artifacts.json"


def artifact_dir(data_dir: Path | None = None) -> Path:
    """Return the directory holding downloaded UniProt artifacts."""
    if data_dir is None:
        env = os.environ.get("MHCSEQS_DATA")
        data_dir = Path(env) if env else Path.home() / ".cache" / "mhcseqs"
    return data_dir / "uniprot"


@dataclass(frozen=True)
class Artifact:
    """One downloaded artifact and the UniProt release used to retrieve it."""

    name: str
    path: Path
    release: str
    rows: int


def corpus_query() -> str:
    """Return the broad query saved as the current-candidate artifact."""
    return f"(taxonomy_id:{VERTEBRATA_TAXON_ID}) AND {MHC_QUERY_CLAUSE}"


def _request(url: str, *, timeout: int):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(request, timeout=timeout)


def _download(url: str, dest: Path) -> str:
    """Stream *url* to *dest* atomically and return its UniProt release."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=dest.parent, suffix=".tmp")
    try:
        os.close(fd)
        with _request(url, timeout=600) as response:
            release = response.headers.get("X-UniProt-Release", "")
            with open(tmp_path, "wb") as handle:
                while chunk := response.read(1 << 20):
                    handle.write(chunk)
        os.replace(tmp_path, dest)
    except BaseException:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise
    return release


def download_corpus_artifact(data_dir: Path | None = None) -> Artifact:
    """Download every broad MHC-named Vertebrata candidate as one file."""
    dest = artifact_dir(data_dir) / CORPUS_ARTIFACT
    url = f"{STREAM_URL}?" + urllib.parse.urlencode(
        {
            "query": corpus_query(),
            "format": "tsv",
            "compressed": "true",
            "fields": ",".join(ARTIFACT_FIELDS),
        }
    )
    release = _download(url, dest)
    rows = sum(1 for _ in read_artifact(dest))
    return Artifact(CORPUS_ARTIFACT, dest, release, rows)


def _flat_file_values(text: str, code: str) -> str:
    return " ".join(line[5:].strip() for line in text.splitlines() if line.startswith(code + "   "))


def _strip_evidence(value: str) -> str:
    return re.sub(r"\s*\{ECO:[^}]+\}", "", value).strip()


def _flat_file_feature(text: str, feature: str) -> str:
    """Render all flat-file instances of *feature* like a TSV FT field."""
    rendered: list[str] = []
    current: list[str] | None = None
    for line in text.splitlines():
        if not line.startswith("FT"):
            if current:
                rendered.append("; ".join(current))
                current = None
            continue
        key = line[5:21].strip()
        value = line[21:].strip()
        if key:
            if current:
                rendered.append("; ".join(current))
            current = [f"{feature} {value}"] if key == feature else None
        elif current and value:
            current.append(value)
    if current:
        rendered.append("; ".join(current))
    return "; ".join(rendered)


def _iso_date(value: str) -> str:
    if not value:
        return ""
    return datetime.strptime(value, "%d-%b-%Y").date().isoformat()


def parse_unisave_record(
    text: str,
    *,
    metadata: dict[str, object],
    source_url: str,
    source_clade: str = "",
) -> dict[str, str]:
    """Parse a pinned UniSave flat file into the current artifact schema."""
    accession_match = re.search(r"^AC\s+([^;]+);", text, re.MULTILINE)
    id_match = re.search(r"^ID\s+(\S+)\s+(Reviewed|Unreviewed);", text, re.MULTILINE)
    taxon_match = re.search(r"^OX\s+NCBI_TaxID=(\d+)", text, re.MULTILINE)
    if accession_match is None or id_match is None or taxon_match is None:
        raise ValueError("UniSave record is missing accession, entry, or taxon metadata")

    accession = accession_match.group(1)
    expected_accession = str(metadata.get("accession", ""))
    if expected_accession and accession != expected_accession:
        raise ValueError(f"UniSave returned {accession} while retrieving {expected_accession}")

    sequence_lines = text.partition("\nSQ   ")[2].partition("\n//")[0].splitlines()[1:]
    sequence = "".join(re.sub(r"[^A-Za-z]", "", line) for line in sequence_lines).upper()
    if not sequence:
        raise ValueError(f"UniSave record {accession} has no sequence")

    gene_text = _flat_file_values(text, "GN")
    gene_names: list[str] = []
    for value in re.findall(r"(?:Name|Synonyms|OrderedLocusNames|ORFNames)=([^;]+)", gene_text):
        gene_names.extend(_strip_evidence(value).replace(",", " ").split())

    protein_text = _flat_file_values(text, "DE")
    protein_names = [_strip_evidence(value) for value in re.findall(r"Full=([^;]+)", protein_text)]
    organism = _flat_file_values(text, "OS").removesuffix(".")
    lineage = _flat_file_values(text, "OC").removesuffix(".").replace(";", ",")
    keywords = _strip_evidence(_flat_file_values(text, "KW")).removesuffix(";")
    protein_existence = _flat_file_values(text, "PE").removesuffix(";")
    protein_existence = re.sub(r"^\d+:\s*", "", protein_existence)
    last_release = str(metadata.get("lastRelease", ""))
    last_release_date = str(metadata.get("lastReleaseDate", ""))

    return {
        "Entry": accession,
        "Entry Name": id_match.group(1),
        "Organism": organism,
        "Organism (ID)": taxon_match.group(1),
        "Taxonomic lineage (Ids)": "",
        "Taxonomic lineage": lineage,
        "Protein names": " (".join(protein_names) + ")" * max(len(protein_names) - 1, 0),
        "Gene Names": " ".join(gene_names),
        "Reviewed": id_match.group(2).lower(),
        "Protein existence": protein_existence,
        "Fragment": "fragment" if re.search(r"^ID\s+.*;\s+Fragment;", text, re.MULTILINE) else "",
        "Length": str(len(sequence)),
        "Sequence": sequence,
        "Signal peptide": _flat_file_feature(text, "SIGNAL"),
        "Chain": _flat_file_feature(text, "CHAIN"),
        "Keywords": keywords,
        "Date of last modification": _iso_date(last_release_date),
        "Entry version": str(metadata.get("entryVersion", "")),
        "Source clade": source_clade,
        "Archive release": last_release.split("/", 1)[0],
        "Archive release date": last_release_date,
        "Archive source URL": source_url,
    }


def _fetch_deleted_entry(args: tuple[str, str]) -> tuple[dict[str, str], str]:
    accession, source_clade = args
    metadata_url = f"{UNISAVE_URL}/{accession}?format=json"
    with _request(metadata_url, timeout=120) as response:
        release = response.headers.get("X-UniProt-Release", "")
        payload = json.loads(response.read().decode("utf-8"))
    versions = payload.get("results", [])
    if not versions:
        raise ValueError(f"UniSave contains no versions for {accession}")

    latest = max(versions, key=lambda row: int(row["entryVersion"]))
    version = str(latest["entryVersion"])
    query = urllib.parse.urlencode({"format": "txt", "versions": version})
    source_url = f"{UNISAVE_URL}/{accession}?{query}"
    with _request(source_url, timeout=120) as response:
        text = response.read().decode("utf-8")
    return (
        parse_unisave_record(
            text,
            metadata=latest,
            source_url=source_url,
            source_clade=source_clade,
        ),
        release,
    )


def _write_gzip_tsv(path: Path, rows: list[dict[str, str]], fieldnames: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
    try:
        os.close(fd)
        with gzip.open(tmp_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_path, path)
    except BaseException:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise


def download_deleted_artifact(
    accessions: list[str],
    data_dir: Path | None = None,
    *,
    source_clades: dict[str, str] | None = None,
) -> Artifact:
    """Save complete latest UniSave records for deleted corpus accessions.

    Every requested accession must have a history and a parseable archived
    record. The artifact is replaced only after all requested rows succeed.
    """
    dest = artifact_dir(data_dir) / DELETED_ARTIFACT
    requested = sorted(set(accessions))
    source_clades = source_clades or {}
    jobs = [(accession, source_clades.get(accession, "")) for accession in requested]
    with ThreadPoolExecutor(max_workers=8) as pool:
        fetched = list(pool.map(_fetch_deleted_entry, jobs))

    rows = [row for row, _release in fetched]
    if {row["Entry"] for row in rows} != set(requested):
        raise ValueError("UniSave artifact does not contain every requested accession")
    releases = {release for _row, release in fetched if release}
    if len(releases) > 1:
        raise ValueError(f"UniProt release changed during UniSave retrieval: {sorted(releases)}")

    _write_gzip_tsv(dest, rows, DELETED_COLUMNS)
    release = releases.pop() if releases else ""
    return Artifact(DELETED_ARTIFACT, dest, release, len(rows))


def read_artifact(path: Path):
    """Yield rows from a gzipped TSV artifact without loading it all at once."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _load_artifact(name: str, data_dir: Path | None) -> list[dict[str, str]]:
    path = artifact_dir(data_dir) / name
    if not path.exists():
        raise FileNotFoundError(f"Missing UniProt artifact: {path}. Run scripts/fetch_sp_ground_truth.py --refresh")
    return list(read_artifact(path))


def load_corpus_artifact(data_dir: Path | None = None) -> list[dict[str, str]]:
    """Read the saved current-candidate artifact."""
    return _load_artifact(CORPUS_ARTIFACT, data_dir)


def load_deleted_artifact(data_dir: Path | None = None) -> list[dict[str, str]]:
    """Read the saved complete UniSave records."""
    return _load_artifact(DELETED_ARTIFACT, data_dir)


def source_clade_from_artifact_row(row: dict[str, str]) -> str:
    """Resolve exactly one exhaustive corpus clade for an artifact row."""
    stored = row.get("Source clade", "").strip()
    if stored:
        if stored not in {name for name, _taxon_id in CORPUS_CLADES}:
            raise ValueError(f"Unknown stored source clade {stored!r} for {row.get('Entry', '')}")
        return stored

    lineage_ids = {int(value) for value in re.findall(r"\d+", row.get("Taxonomic lineage (Ids)", ""))}
    organism_id = row.get("Organism (ID)", "").strip()
    if organism_id:
        lineage_ids.add(int(organism_id))
    matches = [name for name, taxon_id in CORPUS_CLADES if taxon_id in lineage_ids]
    if len(matches) != 1:
        raise ValueError(f"Expected one source clade for {row.get('Entry', '')}, found {matches!r}")
    return matches[0]


def write_manifest(artifacts: list[Artifact], data_dir: Path | None = None) -> Path:
    """Record the release, row count, and size of each saved artifact."""
    path = artifact_dir(data_dir) / RELEASE_MANIFEST
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {a.name: {"release": a.release, "rows": a.rows, "bytes": a.path.stat().st_size} for a in artifacts}
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return path


def read_manifest(data_dir: Path | None = None) -> dict[str, dict[str, object]]:
    """Read artifact release metadata, or return an empty mapping."""
    path = artifact_dir(data_dir) / RELEASE_MANIFEST
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)
