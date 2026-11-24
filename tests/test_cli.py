import csv
import json
import re
import sys
from types import SimpleNamespace
from unittest.mock import patch

import pytest

import mhcseqs
from mhcseqs import default_data_dir
from mhcseqs.__main__ import main


def test_default_data_dir():
    d = default_data_dir()
    assert isinstance(d, str)
    assert "mhcseqs" in d


def test_version_flag(capsys):
    with patch.object(sys, "argv", ["mhcseqs", "--version"]):
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "mhcseqs" in captured.out
    assert "." in captured.out  # version string contains a dot


def test_no_args_prints_help(capsys):
    with patch.object(sys, "argv", ["mhcseqs"]):
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1


def test_data_list_reports_versioned_dataset(capsys, tmp_path, monkeypatch):
    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    with patch.object(sys, "argv", ["mhcseqs", "data", "list"]):
        main()
    captured = capsys.readouterr()
    assert re.search(r"mhc-proteins\s+uniprot-2026_03-r1\s+not installed", captured.out)


def test_data_install_dispatches_version_and_sources(capsys, tmp_path, monkeypatch):
    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    records = tmp_path / "records.csv.gz"
    manifest = tmp_path / "manifest.json"
    source_root = tmp_path / "source-bundle"
    curation = source_root / "sp_ground_truth_label_curation.csv"
    seen = []

    def install(version, *, data_dir, force):
        seen.append(("records", version, data_dir, force))
        return SimpleNamespace(version=version, records=records, manifest=manifest)

    def install_sources(version, *, data_dir, force):
        seen.append(("sources", version, data_dir, force))
        return SimpleNamespace(root=source_root, label_curation=curation)

    monkeypatch.setattr("mhcseqs.__main__.install_mhc_protein_dataset", install)
    monkeypatch.setattr("mhcseqs.__main__.install_mhc_protein_source_bundle", install_sources)
    with patch.object(
        sys,
        "argv",
        ["mhcseqs", "data", "install", "mhc-proteins", "--version", "uniprot-2026_03-r1", "--with-sources"],
    ):
        main()

    assert seen == [
        ("records", "uniprot-2026_03-r1", tmp_path, False),
        ("sources", "uniprot-2026_03-r1", tmp_path, False),
    ]
    assert str(records) in capsys.readouterr().out


def test_data_path_prints_only_validated_cache(capsys, tmp_path, monkeypatch):
    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    records = tmp_path / "records.csv.gz"
    seen = []

    def validate(version, *, data_dir):
        seen.append((version, data_dir))
        return SimpleNamespace(version="uniprot-2026_03-r1", records=records)

    monkeypatch.setattr("mhcseqs.__main__.validate_mhc_protein_dataset", validate)
    with patch.object(sys, "argv", ["mhcseqs", "data", "path", "mhc-proteins"]):
        main()

    assert seen == [(None, tmp_path)]
    assert capsys.readouterr().out.strip() == str(records)


def test_data_path_rejects_corrupt_cache(capsys, tmp_path, monkeypatch):
    from mhcseqs.mhc_protein_dataset import mhc_protein_dataset_paths

    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    paths = mhc_protein_dataset_paths(data_dir=tmp_path)
    paths.records.parent.mkdir(parents=True)
    paths.records.write_bytes(b"truncated")

    with patch.object(sys, "argv", ["mhcseqs", "data", "path", "mhc-proteins"]):
        with pytest.raises(SystemExit) as exc_info:
            main()

    assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "missing or corrupt" in captured.err
    assert "--force" in captured.err


def test_data_list_marks_corrupt_version_invalid(capsys, tmp_path, monkeypatch):
    from mhcseqs.mhc_protein_dataset import mhc_protein_dataset_paths

    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    paths = mhc_protein_dataset_paths(data_dir=tmp_path)
    paths.records.parent.mkdir(parents=True)
    paths.records.write_bytes(b"truncated")

    with patch.object(sys, "argv", ["mhcseqs", "data", "list"]):
        main()

    output = capsys.readouterr().out
    assert re.search(rf"mhc-proteins\s+{paths.version}\s+invalid; reinstall with --force", output)


def test_lookup_no_csvs(capsys, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr("mhcseqs.datafiles.__file__", str(tmp_path / "pkg" / "datafiles.py"))
    with patch.object(sys, "argv", ["mhcseqs", "lookup", "HLA-A*02:01", "--output-dir", str(tmp_path)]):
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert "No built CSVs found" in captured.err


# ---------------------------------------------------------------------------
# lookup() API tests
# ---------------------------------------------------------------------------

# HLA-A*02:01 full sequence (with signal peptide)
_A0201_SEQ = (
    "MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLG"
    "TLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRA"
    "YLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGD"
    "GTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEPSSQPTIPIVGIIAGLVLFGAVITGAVVAAVMWRRKSSDRKGGSYSQAASSDSAQGSDVSLTACKV"
)


def _write_full_seqs_csv(path, rows):
    """Write a minimal full-seqs CSV for testing."""
    from mhcseqs.pipeline import FULL_FIELDS

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FULL_FIELDS)
        writer.writeheader()
        for row in rows:
            full = {k: "" for k in FULL_FIELDS}
            full.update(row)
            writer.writerow(full)


def test_lookup_finds_allele(tmp_path, monkeypatch):
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [
            {
                "two_field_allele": "HLA-A*02:01",
                "representative_allele": "HLA-A*02:01:01:01",
                "gene": "A",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Homo sapiens",
                "species_category": "human",
                "sequence": _A0201_SEQ,
                "seq_len": str(len(_A0201_SEQ)),
            },
        ],
    )
    monkeypatch.chdir(tmp_path)
    rec = mhcseqs.lookup("HLA-A*02:01")
    assert isinstance(rec, mhcseqs.AlleleRecord)
    assert rec.ok
    assert rec.allele == "HLA-A*02:01"
    assert rec.full_allele == "HLA-A*02:01:01:01"
    assert rec.species_category == "human"
    assert len(rec.groove1) == 90
    assert len(rec.groove2) == 93
    assert rec.anchor_cys1 is not None  # populated by extract_groove


def test_lookup_case_insensitive(tmp_path, monkeypatch):
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [
            {
                "two_field_allele": "HLA-A*02:01",
                "mhc_class": "I",
                "chain": "alpha",
                "sequence": _A0201_SEQ,
            },
        ],
    )
    monkeypatch.chdir(tmp_path)
    rec = mhcseqs.lookup("hla-a*02:01")
    assert rec.ok
    assert len(rec.groove1) == 90


def test_lookup_with_mutations(tmp_path, monkeypatch):
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [
            {
                "two_field_allele": "HLA-A*02:01",
                "mhc_class": "I",
                "chain": "alpha",
                "sequence": _A0201_SEQ,
            },
        ],
    )
    monkeypatch.chdir(tmp_path)
    rec = mhcseqs.lookup("HLA-A*02:01", mutations=["K66A"])
    assert rec.groove1[65] == "A"
    assert rec.mature_sequence[65] == "A"
    assert rec.mutations == ("K66A",)

    # Single string mutation works without character iteration (#121)
    rec_str = mhcseqs.lookup("HLA-A*02:01", mutations="K66A")
    assert rec_str.groove1[65] == "A"
    assert rec_str.mutations == ("K66A",)


def test_lookup_resolves_curated_non_human_and_genus_aliases(tmp_path, monkeypatch):
    """Issues #117 and #119: lookup resolves FLA-*, Trvu-*, Ptal-*, Anpl-*, SahaI*, H2-Kb, and UniProt IDs."""
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [
            {
                "two_field_allele": "FelisCatus-E*003:01",
                "representative_allele": "FelisCatus-E*003:01",
                "gene": "E",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Felis catus",
                "species_category": "carnivore",
                "species_prefix": "FelisCatus",
                "source": "uniprot",
                "source_id": "Q95482",
                "sequence": _A0201_SEQ,
            },
            {
                "two_field_allele": "TrichosurusVulpecula-UB*01:01",
                "representative_allele": "TrichosurusVulpecula-UB*01:01",
                "gene": "UB",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Trichosurus vulpecula",
                "species_category": "other_mammal",
                "species_prefix": "TrichosurusVulpecula",
                "source": "uniprot",
                "source_id": "Q95IT0",
                "sequence": _A0201_SEQ,
            },
            {
                "two_field_allele": "SarcophilusHarrisii-I*27",
                "representative_allele": "SarcophilusHarrisii-I*27",
                "gene": "I",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Sarcophilus harrisii",
                "species_category": "other_mammal",
                "species_prefix": "SarcophilusHarrisii",
                "source": "genbank",
                "source_id": "GQ411435.1",
                "sequence": _A0201_SEQ,
            },
            {
                "two_field_allele": "H2-K*b",
                "representative_allele": "H2-K*b",
                "gene": "K",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Mus musculus",
                "species_category": "murine",
                "species_prefix": "H2",
                "source": "uniprot",
                "source_id": "P01901",
                "sequence": _A0201_SEQ,
            },
            {
                "two_field_allele": "BoLA-1*007:01",
                "representative_allele": "BoLA-1*007:01:01:01",
                "gene": "1",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Bos sp.",
                "species_category": "ungulate",
                "species_prefix": "BoLA",
                "source": "ipd_mhc",
                "source_id": "BoLA00001",
                "sequence": _A0201_SEQ,
            },
            {
                "two_field_allele": "SLA-1*01:01",
                "representative_allele": "SLA-1*01:01",
                "gene": "1",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Sus sp.",
                "species_category": "ungulate",
                "species_prefix": "SLA",
                "source": "ipd_mhc",
                "source_id": "SLA00001",
                "sequence": _A0201_SEQ,
            },
        ],
    )
    monkeypatch.chdir(tmp_path)
    for query, expected in (
        ("FLA-E*003:01", "FelisCatus-E*003:01"),
        ("FLA-E*00301", "FelisCatus-E*003:01"),
        ("Feca-E*003:01", "FelisCatus-E*003:01"),
        ("FelisCatus-E*003:01", "FelisCatus-E*003:01"),
        ("Q95482", "FelisCatus-E*003:01"),
        ("Trvu-UB*01:01", "TrichosurusVulpecula-UB*01:01"),
        ("Q95IT0", "TrichosurusVulpecula-UB*01:01"),
        ("SahaI*27", "SarcophilusHarrisii-I*27"),
        ("Saha-I*27", "SarcophilusHarrisii-I*27"),
        ("H2-Kb", "H2-K*b"),
        ("BoLA-1*007:01", "BoLA-1*007:01"),
        ("Bota-1*007:01", "BoLA-1*007:01"),
        ("SLA-1*01:01", "SLA-1*01:01"),
        ("Susc-1*01:01", "SLA-1*01:01"),
    ):
        rec = mhcseqs.lookup(query)
        assert rec.allele == expected, f"query {query!r} resolved to {rec.allele!r}, expected {expected!r}"


def test_lookup_not_found(tmp_path, monkeypatch):
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [
            {"two_field_allele": "HLA-A*02:01", "mhc_class": "I", "sequence": _A0201_SEQ},
        ],
    )
    monkeypatch.chdir(tmp_path)
    with pytest.raises(KeyError, match="HLA-B"):
        mhcseqs.lookup("HLA-B*07:02")


def test_lookup_no_csv(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        mhcseqs,
        "find_csv",
        lambda name, **kw: (_ for _ in ()).throw(FileNotFoundError(f"{name} not found. Run mhcseqs.build() or 'mhcseqs build' first.")),
    )
    with pytest.raises(FileNotFoundError, match="build"):
        mhcseqs.lookup("HLA-A*02:01")


def _write_raw_seqs_csv(path, rows):
    """Write a minimal raw-seqs CSV for testing."""
    from mhcseqs.pipeline import RAW_FIELDS

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=RAW_FIELDS)
        writer.writeheader()
        for row in rows:
            full = {k: "" for k in RAW_FIELDS}
            full.update(row)
            writer.writerow(full)


def test_cli_lookup_text_prints_full_sequences(capsys, tmp_path, monkeypatch):
    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    monkeypatch.chdir(tmp_path)
    mature = _A0201_SEQ[24:]
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [
            {
                "two_field_allele": "HomoSapiens-A*02:01",
                "representative_allele": "HomoSapiens-A*02:01:01:01",
                "gene": "A",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Homo sapiens",
                "sequence": _A0201_SEQ,
                "mature_sequence": mature,
            },
        ],
    )
    main(["lookup", "HLA-A*02:01", "--output-dir", str(tmp_path)])
    captured = capsys.readouterr()
    assert f"  sequence: {_A0201_SEQ}" in captured.out
    assert f"  mature_sequence: {mature}" in captured.out
    assert "... (" not in captured.out
    assert captured.err == ""


def test_cli_lookup_formats_and_all_and_output(capsys, tmp_path, monkeypatch):
    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    monkeypatch.chdir(tmp_path)
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [
            {
                "two_field_allele": "HomoSapiens-A*02:01",
                "representative_allele": "HomoSapiens-A*02:01:01:01",
                "gene": "A",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Homo sapiens",
                "sequence": _A0201_SEQ,
                "mature_sequence": _A0201_SEQ[24:],
            },
        ],
    )
    _write_raw_seqs_csv(
        tmp_path / "mhc-seqs-raw.csv",
        [
            {
                "allele_raw": "A*02:01:01:01",
                "allele_normalized": "HomoSapiens-A*02:01:01:01",
                "two_field_allele": "HomoSapiens-A*02:01",
                "gene": "A",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Homo sapiens",
                "sequence": _A0201_SEQ,
            },
            {
                "allele_raw": "A*02:01:01:02",
                "allele_normalized": "HomoSapiens-A*02:01:01:02",
                "two_field_allele": "HomoSapiens-A*02:01",
                "gene": "A",
                "mhc_class": "I",
                "chain": "alpha",
                "species": "Homo sapiens",
                "sequence": _A0201_SEQ,
            },
        ],
    )

    # Default (without --all): 1 match per CSV -> 2 rows in JSON
    main(["lookup", "HLA-A*02:01", "--output-dir", str(tmp_path), "--format", "json"])
    payload = json.loads(capsys.readouterr().out)
    assert len(payload) == 2
    assert [row["csv"] for row in payload] == ["mhc-full-seqs.csv", "mhc-seqs-raw.csv"]
    assert payload[0]["sequence"] == _A0201_SEQ

    # With --all: 1 full + 2 raw -> 3 rows in JSON
    main(["lookup", "HLA-A*02:01", "--output-dir", str(tmp_path), "--format", "json", "--all"])
    all_payload = json.loads(capsys.readouterr().out)
    assert len(all_payload) == 3

    # CSV and TSV
    for fmt, delim in (("csv", ","), ("tsv", "\t")):
        main(["lookup", "HLA-A*02:01", "--output-dir", str(tmp_path), "--format", fmt, "--all"])
        out = capsys.readouterr().out
        rows = list(csv.DictReader(out.splitlines(), delimiter=delim))
        assert len(rows) == 3
        assert rows[0]["sequence"] == _A0201_SEQ

    # FASTA written to --output
    out_fasta = tmp_path / "out" / "a0201.fasta"
    main(["lookup", "HLA-A*02:01", "--output-dir", str(tmp_path), "--format", "fasta", "--output", str(out_fasta)])
    assert capsys.readouterr().out == ""
    fasta_lines = out_fasta.read_text(encoding="utf-8").strip().splitlines()
    assert fasta_lines[0].startswith(">HomoSapiens-A*02:01:01:01 csv=mhc-full-seqs.csv")
    assert fasta_lines[1] == _A0201_SEQ


def test_cli_lookup_no_match_goes_to_stderr(capsys, tmp_path, monkeypatch):
    monkeypatch.setattr("mhcseqs.__main__.default_data_dir", lambda: str(tmp_path))
    monkeypatch.setattr("mhcseqs.datafiles.__file__", str(tmp_path / "pkg" / "datafiles.py"))
    monkeypatch.chdir(tmp_path)
    _write_full_seqs_csv(
        tmp_path / "mhc-full-seqs.csv",
        [{"two_field_allele": "HomoSapiens-A*02:01", "sequence": _A0201_SEQ}],
    )
    with pytest.raises(SystemExit) as exc_info:
        main(["lookup", "HLA-B*07:02", "--output-dir", str(tmp_path)])
    assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "No match found for 'HLA-B*07:02'" in captured.err


def test_cli_and_api_build_use_same_fasta_dir_when_output_dir_is_set(tmp_path, monkeypatch):
    data_dir = tmp_path / "cache"
    out_cli = tmp_path / "out-cli"
    out_api = tmp_path / "out-api"
    monkeypatch.setenv("MHCSEQS_DATA", str(data_dir))
    seen_dirs = []

    def fake_download(dest_dir, force=False):
        seen_dirs.append(dest_dir)
        return {"imgt_hla": dest_dir / "hla_prot.fasta", "ipd_mhc": dest_dir / "ipd_mhc_prot.fasta"}

    monkeypatch.setattr("mhcseqs.__main__.download_all", fake_download)
    monkeypatch.setattr("mhcseqs.__main__.build_raw_index", lambda inputs, raw_csv: {"parsed": 0})
    monkeypatch.setattr("mhcseqs.__main__.build_full_seqs", lambda raw_csv, full_csv, report_path=None: {"total_groups": 0})
    monkeypatch.setattr("mhcseqs.__main__.validate_build", lambda raw_csv, full_csv: ([], {}))
    monkeypatch.setattr("mhcseqs.__main__.format_validation_report", lambda w, s: "ok")
    out_cli.mkdir()
    main(["build", "--output-dir", str(out_cli)])

    monkeypatch.setattr("mhcseqs.download_all", fake_download)
    monkeypatch.setattr("mhcseqs.build_raw_index", lambda inputs, raw_csv: {"parsed": 0})
    monkeypatch.setattr("mhcseqs.build_full_seqs", lambda raw_csv, full_csv, report_path=None: {"total_groups": 0})
    monkeypatch.setattr("mhcseqs.validate.validate_build", lambda raw_csv, full_csv: ([], {}))
    monkeypatch.setattr("mhcseqs.validate.format_validation_report", lambda w, s: "ok")
    out_api.mkdir()
    mhcseqs.build(output_dir=str(out_api))

    assert seen_dirs == [data_dir / "fasta", data_dir / "fasta"]


def test_cli_and_api_lookup_share_csv_search_order(capsys, tmp_path, monkeypatch):
    cwd_dir = tmp_path / "cwd"
    data_dir = tmp_path / "cache"
    cwd_dir.mkdir()
    data_dir.mkdir()
    monkeypatch.setenv("MHCSEQS_DATA", str(data_dir))
    monkeypatch.chdir(cwd_dir)

    _write_full_seqs_csv(
        cwd_dir / "mhc-full-seqs.csv",
        [
            {
                "two_field_allele": "HomoSapiens-A*02:01",
                "representative_allele": "HomoSapiens-A*02:01:01:01",
                "source_id": "FROM_CWD",
                "mhc_class": "I",
                "chain": "alpha",
                "sequence": _A0201_SEQ,
            },
        ],
    )
    _write_full_seqs_csv(
        data_dir / "mhc-full-seqs.csv",
        [
            {
                "two_field_allele": "HomoSapiens-A*02:01",
                "representative_allele": "HomoSapiens-A*02:01:01:01",
                "source_id": "FROM_CACHE",
                "mhc_class": "I",
                "chain": "alpha",
                "sequence": _A0201_SEQ,
            },
        ],
    )

    rec = mhcseqs.lookup("HLA-A*02:01")
    main(["lookup", "HLA-A*02:01", "--format", "json"])
    cli_rows = json.loads(capsys.readouterr().out)
    assert rec.source_id == "FROM_CWD"
    assert cli_rows[0]["source_id"] == "FROM_CWD"


def test_validate_signal_peptides_uses_resolved_data_dir(capsys, tmp_path, monkeypatch):
    from scripts import validate_signal_peptides

    monkeypatch.setenv("MHCSEQS_DATA", str(tmp_path))
    with pytest.raises(SystemExit) as exc_info:
        validate_signal_peptides.main([])
    assert exc_info.value.code == 1
    err = capsys.readouterr().err
    assert str(tmp_path / "mhc-seqs-raw.csv") in err
    assert "mhcseqs build" in err

    custom = tmp_path / "custom-raw.csv"
    with pytest.raises(SystemExit) as exc_info:
        validate_signal_peptides.main(["--raw-csv", str(custom)])
    assert exc_info.value.code == 1
    assert str(custom) in capsys.readouterr().err


def test_main_accepts_explicit_argv(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(["--version"])
    assert exc_info.value.code == 0
    assert mhcseqs.__version__ in capsys.readouterr().out


def test_package_all_lists_every_public_export():
    assert len(mhcseqs.__all__) == len(set(mhcseqs.__all__))
    assert all(hasattr(mhcseqs, name) for name in mhcseqs.__all__)
    for name in ("SequenceFeatures", "StructuralDomain", "analyze_sequence", "infer_structural_domains", "lookup", "build"):
        assert name in mhcseqs.__all__
