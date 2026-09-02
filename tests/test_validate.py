import csv

from mhcseqs.validate import (
    MIN_RELEASE_ARTIFACT_COUNTS,
    VALID_AA,
    _check_b2m,
    _check_groove_row,
    _check_mature_sequence,
    _check_signal_peptide,
    _check_valid_aa,
    format_validation_report,
    validate_release_artifact_counts,
)


def test_valid_aa_set():
    # All 20 standard amino acids + X
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        assert aa in VALID_AA
    assert "X" in VALID_AA
    # These should not be in the set
    assert "J" not in VALID_AA
    assert "O" not in VALID_AA
    assert "U" not in VALID_AA


def test_check_valid_aa_clean():
    warnings = _check_valid_aa("ACDEFGHIKLMNPQRSTVWY", "test")
    assert warnings == []


def test_check_valid_aa_invalid():
    warnings = _check_valid_aa("ACDEFJ", "test")
    assert len(warnings) == 1
    assert "invalid characters" in warnings[0]


def test_check_signal_peptide_no_sp():
    row = {"has_signal_peptide": "False", "signal_peptide_len": "0"}
    assert _check_signal_peptide(row) == []


def test_check_signal_peptide_valid():
    row = {
        "has_signal_peptide": "True",
        "signal_peptide_len": "24",
        "signal_peptide_seq": "MAVMAPRTLVLLLSGALALTQTWA",
        "allele_normalized": "HLA-A*02:01",
        "mhc_class": "I",
        "chain": "alpha",
    }
    warnings = _check_signal_peptide(row)
    assert warnings == []


def test_check_signal_peptide_no_met():
    row = {
        "has_signal_peptide": "True",
        "signal_peptide_len": "20",
        "signal_peptide_seq": "XAVMAPRTLVLLLSGALALT",
        "allele_normalized": "test",
        "mhc_class": "I",
        "chain": "alpha",
    }
    warnings = _check_signal_peptide(row)
    assert any("does not start with M" in w for w in warnings)


def test_check_signal_peptide_unusual_length():
    row = {
        "has_signal_peptide": "True",
        "signal_peptide_len": "50",
        "signal_peptide_seq": "M" + "A" * 49,
        "allele_normalized": "test",
        "mhc_class": "I",
        "chain": "alpha",
    }
    warnings = _check_signal_peptide(row)
    assert any("unusual SP length" in w for w in warnings)


def test_check_b2m_valid():
    row = {
        "sequence": "M" + "A" * 20 + "CC" + "A" * 77,  # 100 aa, 2 Cys
        "allele_normalized": "B2M_human",
    }
    assert _check_b2m(row) == []


def test_check_b2m_empty():
    row = {"sequence": "", "allele_normalized": "B2M_test"}
    warnings = _check_b2m(row)
    assert any("empty sequence" in w for w in warnings)


def test_check_b2m_unusual_length():
    row = {"sequence": "M" + "CC" + "A" * 7, "allele_normalized": "B2M_test"}
    warnings = _check_b2m(row)
    assert any("unusual length" in w for w in warnings)


def test_check_groove_row_not_applicable():
    row = {"groove_status": "not_applicable", "two_field_allele": "MICA"}
    assert _check_groove_row(row) == []


def test_check_groove_row_ok_class_i():
    row = {
        "two_field_allele": "HLA-A*02:01",
        "groove_status": "ok",
        "mhc_class": "I",
        "chain": "alpha",
        "groove1": "A" * 90,
        "groove2": "A" * 93,
    }
    assert _check_groove_row(row) == []


def test_check_groove_row_unusual_length():
    row = {
        "two_field_allele": "test",
        "groove_status": "ok",
        "mhc_class": "I",
        "chain": "alpha",
        "groove1": "A" * 50,  # too short
        "groove2": "A" * 93,
    }
    warnings = _check_groove_row(row)
    assert any("alpha1 half length" in w for w in warnings)


def test_check_mature_sequence_valid():
    row = {
        "two_field_allele": "test",
        "mature_sequence": "A" * 300,
        "mhc_class": "I",
        "chain": "alpha",
    }
    assert _check_mature_sequence(row) == []


def test_check_mature_sequence_too_many_x():
    row = {
        "two_field_allele": "test",
        "mature_sequence": "X" * 50 + "A" * 50,
        "mhc_class": "I",
        "chain": "alpha",
    }
    warnings = _check_mature_sequence(row)
    assert any("unknown residues" in w for w in warnings)


def test_format_validation_report_no_warnings():
    report = format_validation_report([], {"raw_total": 100, "full_total": 90, "groove_total": 80, "groove_ok": 75})
    assert "Validation Summary" in report
    assert "No warnings!" in report
    assert "100" in report


def test_format_validation_report_with_warnings():
    warnings = ["SP(test): unusual SP length 50 (expected 15-35)"]
    report = format_validation_report(warnings, {"raw_total": 10, "full_total": 10, "groove_total": 10, "groove_ok": 9})
    assert "Warnings: 1 total" in report
    assert "unusual SP length" in report


def test_format_validation_report_shows_all_by_default():
    """All warnings should be shown when max_warnings_per_type is None."""
    warnings = [f"SP(allele_{i}): unusual SP length {40 + i}" for i in range(10)]
    report = format_validation_report(warnings, {"raw_total": 10, "full_total": 10})
    # All 10 warnings should appear
    for i in range(10):
        assert f"allele_{i}" in report
    assert "... and" not in report


def test_format_validation_report_truncates_when_requested():
    """max_warnings_per_type should limit output per category."""
    warnings = [f"SP(allele_{i}): unusual SP length {40 + i}" for i in range(10)]
    report = format_validation_report(warnings, {"raw_total": 10, "full_total": 10}, max_warnings_per_type=3)
    assert "allele_0" in report
    assert "allele_2" in report
    assert "allele_3" not in report
    assert "... and 7 more" in report


def _write_count_fixture(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_validate_release_artifact_counts_accepts_healthy_counts(tmp_path):
    raw_csv = tmp_path / "raw.csv"
    full_csv = tmp_path / "full.csv"
    _write_count_fixture(
        raw_csv,
        ["source", "species"],
        [
            {"source": "imgt", "species": "Homo sapiens"},
            {"source": "imgt", "species": "Homo sapiens"},
            {"source": "ipd_mhc", "species": "Macaca mulatta"},
        ],
    )
    _write_count_fixture(
        full_csv,
        ["allele", "species", "mhc_class", "is_functional"],
        [
            {
                "allele": "a",
                "species": "Homo sapiens",
                "mhc_class": "I",
                "is_functional": "True",
            },
            {
                "allele": "b",
                "species": "Macaca mulatta",
                "mhc_class": "II",
                "is_functional": "True",
            },
        ],
    )

    errors, counts = validate_release_artifact_counts(
        raw_csv,
        full_csv,
        minimums={"raw_total": 3, "full_total": 2, "raw_imgt": 2, "raw_human": 2},
    )

    assert errors == []
    assert counts == {
        "raw_total": 3,
        "full_total": 2,
        "raw_imgt": 2,
        "raw_human": 2,
        "functional_total": 2,
        "functional_human": 1,
        "functional_class_i": 1,
        "functional_class_ii": 1,
    }


def test_validate_release_artifact_counts_reports_each_regression(tmp_path):
    raw_csv = tmp_path / "raw.csv"
    full_csv = tmp_path / "full.csv"
    _write_count_fixture(raw_csv, ["source", "species"], [{"source": "ipd_mhc", "species": "Macaca mulatta"}])
    _write_count_fixture(full_csv, ["allele"], [])

    errors, _ = validate_release_artifact_counts(
        raw_csv,
        full_csv,
        minimums={"raw_total": 2, "full_total": 1, "raw_imgt": 1, "raw_human": 1},
    )

    assert len(errors) == 4
    assert any("raw_imgt=0" in error for error in errors)
    assert any("raw_human=0" in error for error in errors)


def test_release_artifact_minimums_cover_hla_source_loss():
    assert MIN_RELEASE_ARTIFACT_COUNTS["raw_imgt"] >= 20_000
    assert MIN_RELEASE_ARTIFACT_COUNTS["raw_human"] >= 20_000
    assert MIN_RELEASE_ARTIFACT_COUNTS["functional_human"] >= 20_000
    assert MIN_RELEASE_ARTIFACT_COUNTS["functional_class_i"] >= 30_000
    assert MIN_RELEASE_ARTIFACT_COUNTS["functional_class_ii"] >= 18_000
