from collections import Counter

from scripts import update_readme_counts


def test_load_built_counts_keeps_helper_chains_out_of_class_i(monkeypatch, tmp_path):
    built_csv = tmp_path / "mhc-full-seqs.csv"
    built_csv.write_text(
        "species_category,mhc_class,chain\nfish,I,alpha\nfish,I,B2M\nfish,unknown,unknown\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(update_readme_counts, "BUILT_CSV", built_csv)

    assert update_readme_counts.load_built_counts() == Counter({("fish", "I"): 1, ("fish", "other"): 2})


def test_main_counts_only_completed_build(monkeypatch, tmp_path):
    readme = tmp_path / "README.md"
    readme.write_text(
        "# Test\n\n## Current data summary\n\nstale\n\n## Structural decomposition\n",
        encoding="utf-8",
    )
    built = Counter({("human", "I"): 1, ("human", "other"): 3, ("unclassified", "II"): 2})
    monkeypatch.setattr(update_readme_counts, "README", readme)
    monkeypatch.setattr(update_readme_counts, "load_built_counts", lambda: built)
    monkeypatch.setattr(
        update_readme_counts,
        "load_diverse_species_count",
        lambda: 7,
    )
    monkeypatch.setattr(
        update_readme_counts,
        "load_groove_success_rate",
        lambda: (99.0, 99, 100),
    )

    update_readme_counts.main()

    text = readme.read_text(encoding="utf-8")
    assert "| human | 1 | 0 | 3 | 4 |" in text
    assert "| bird | 0 | 0 | 0 | 0 |" in text
    assert "| unclassified | 0 | 2 | 0 | 2 |" in text
    assert "| **total** | **1** | **2** | **3** | **6** |" in text
