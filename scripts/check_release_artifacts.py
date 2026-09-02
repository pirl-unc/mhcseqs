#!/usr/bin/env python3
"""Fail when built release CSVs have catastrophic source-count regressions."""

from __future__ import annotations

import argparse

from mhcseqs.validate import MIN_RELEASE_ARTIFACT_COUNTS, validate_release_artifact_counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw_csv")
    parser.add_argument("full_csv")
    args = parser.parse_args()

    errors, counts = validate_release_artifact_counts(args.raw_csv, args.full_csv)
    print("Release artifact counts:")
    for name, minimum in MIN_RELEASE_ARTIFACT_COUNTS.items():
        print(f"  {name:12s} {counts[name]:8,} (minimum {minimum:,})")

    if errors:
        print("Release artifact validation failed:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("Release artifact count gates passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
