#!/bin/bash
set -euo pipefail

echo "Running lint checks..."
bash lint.sh

echo "Running tests..."
bash test.sh

echo "Building distribution..."
uv pip install build twine
rm -rf dist/
uv run python -m build

echo "Uploading to PyPI..."
uv run python -m twine upload dist/*

echo "Done!"
