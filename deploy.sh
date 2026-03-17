#!/bin/bash
set -euo pipefail

echo "Running lint checks..."
bash lint.sh

echo "Running tests..."
bash test.sh

echo "Building distribution..."
uv pip install build twine
rm -rf dist/
# Use pyproject-build instead of `python -m build` to avoid conflict
# with mhcseqs.__main__ shadowing the build module.
pyproject-build

echo "Uploading to PyPI..."
twine upload dist/*

echo "Done!"
