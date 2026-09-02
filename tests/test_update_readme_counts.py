from collections import Counter

from scripts import update_readme_counts


def test_main_does_not_add_supplemental_counts_twice(monkeypatch, tmp_path):
    readme = tmp_path / "README.md"
    readme.write_text(
        "# Test\n\n## Current data summary\n\nstale\n\n## Structural decomposition\n",
        encoding="utf-8",
    )
    built = Counter({("human", "I"): 1})
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
    assert "| **total** | **1** | **0** | **1** |" in text
