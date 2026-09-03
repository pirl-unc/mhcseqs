"""Downloaded UniProt artifacts backing the SP benchmark corpus.

The benchmark pipeline used to issue live, parameterized UniProt REST queries at
run time: eight paginated ``uniprotkb/search`` calls whose 500-row page cap
silently truncated three clades, plus 253 ``taxonomy/{id}`` calls, plus a
per-accession lookup in the enrichment step. Those queries were unreproducible
(the answer changed with the release), rate-limited, and their filters were
baked into the URL, so changing one meant re-querying.

Everything is now derived from saved artifacts instead. Each artifact is
downloaded once by ``--refresh``, cached under the mhcseqs data directory, and
stamped with the UniProt release it came from. Downstream code reads only local
files, so a filter change is a local reprocessing step.

The artifacts are deliberately *exhaustive*: no signal-peptide requirement, no
length window, no hand-picked clade list. Every such narrowing happens locally
in ``fetch_sp_ground_truth`` against the saved file.

Artifacts:
  uniprot_mhc_vertebrata.tsv.gz  every MHC-named Vertebrata entry (~55k rows)
  uniprot_deleted_entries.tsv.gz entries the corpus references that UniProt has
                                 since deleted, recovered from the UniSave
                                 version archive (see issue #68)
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import os
import sys
import tempfile
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Vertebrata. A single root, rather than a hand-picked clade list: the eight
# hand-picked roots silently missed coelacanth and lungfish, and nothing checked
# that the list was complete. Clade assignment happens locally from the lineage
# column, and completeness is asserted rather than assumed.
VERTEBRATA_TAXON_ID = 7742

# Matches how the corpus has always identified MHC entries. Kept deliberately
# broad; precision is a local concern.
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

STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
UNISAVE_URL = "https://rest.uniprot.org/unisave"
USER_AGENT = "mhcseqs-sp-corpus/1.0"

CORPUS_ARTIFACT = "uniprot_mhc_vertebrata.tsv.gz"
DELETED_ARTIFACT = "uniprot_deleted_entries.tsv.gz"
RELEASE_MANIFEST = "uniprot_artifacts.json"


def artifact_dir(data_dir: Path | None = None) -> Path:
    """Directory holding downloaded UniProt artifacts."""
    if data_dir is None:
        env = os.environ.get("MHCSEQS_DATA")
        data_dir = Path(env) if env else Path.home() / ".cache" / "mhcseqs"
    return data_dir / "uniprot"


@dataclass
class Artifact:
    """One downloaded artifact and the UniProt release it came from."""

    name: str
    path: Path
    release: str
    rows: int


def corpus_query() -> str:
    """The exhaustive query whose result is saved as the corpus artifact."""
    return f"(taxonomy_id:{VERTEBRATA_TAXON_ID}) AND {MHC_QUERY_CLAUSE}"


def _download(url: str, dest: Path) -> str:
    """Stream *url* to *dest* atomically. Returns the UniProt release header."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    fd, tmp_path = tempfile.mkstemp(dir=dest.parent, suffix=".tmp")
    try:
        os.close(fd)
        with urllib.request.urlopen(request, timeout=600) as response:
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
    """Download every MHC-named Vertebrata entry as one saved file."""
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


def download_deleted_artifact(accessions: list[str], data_dir: Path | None = None) -> Artifact:
    """Recover entries UniProt has deleted from the UniSave version archive.

    Deleted accessions vanish from the knowledgebase but remain in the archive.
    They are recovered once, here, and saved; the pipeline never queries UniSave.
    """
    dest = artifact_dir(data_dir) / DELETED_ARTIFACT
    dest.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    release = ""
    for accession in sorted(set(accessions)):
        url = f"{UNISAVE_URL}/{accession}?format=json&download=false"
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=120) as response:
            release = release or response.headers.get("X-UniProt-Release", "")
            payload = json.loads(response.read())
        versions = payload.get("results", [])
        if not versions:
            continue
        latest = versions[0]
        rows.append(
            {
                "accession": accession,
                "entry_version": str(latest.get("entryVersion", "")),
                "released": str(latest.get("firstReleaseDate", "")),
                "last_release": str(latest.get("lastReleaseDate", "")),
                "provenance": "deleted_upstream_recovered_from_unisave",
            }
        )
    with gzip.open(dest, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["accession"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return Artifact(DELETED_ARTIFACT, dest, release, len(rows))


def read_artifact(path: Path):
    """Yield rows from a gzipped TSV artifact."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(io.StringIO(handle.read()), delimiter="\t")
        yield from reader


def load_corpus_artifact(data_dir: Path | None = None) -> list[dict[str, str]]:
    """Read the saved corpus artifact, failing loudly if it was never downloaded."""
    path = artifact_dir(data_dir) / CORPUS_ARTIFACT
    if not path.exists():
        raise FileNotFoundError(f"Missing UniProt corpus artifact: {path}. Run scripts/fetch_sp_ground_truth.py --refresh")
    return list(read_artifact(path))


def write_manifest(artifacts: list[Artifact], data_dir: Path | None = None) -> Path:
    """Record which UniProt release each saved artifact came from."""
    path = artifact_dir(data_dir) / RELEASE_MANIFEST
    payload = {a.name: {"release": a.release, "rows": a.rows, "bytes": a.path.stat().st_size} for a in artifacts}
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return path


def read_manifest(data_dir: Path | None = None) -> dict[str, dict[str, object]]:
    path = artifact_dir(data_dir) / RELEASE_MANIFEST
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)
