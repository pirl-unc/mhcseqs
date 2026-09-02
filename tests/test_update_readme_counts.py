from collections import Counter

from scripts import update_readme_counts


def test_main_counts_only_completed_build(monkeypatch, tmp_path):
    readme = tmp_path / "README.md"
    readme.write_text(
        "# Test\n\n## Current data summary\n\nstale\n\n## Structural decomposition\n",
        encoding="utf-8",
    )
    built = Counter({("human", "I"): 1, ("unclassified", "II"): 2})
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
    assert "| human | 1 | 0 | 1 |" in text
    assert "| bird | 0 | 0 | 0 |" in text
    assert "| unclassified | 0 | 2 | 2 |" in text
    assert "| **total** | **1** | **2** | **3** |" in text
