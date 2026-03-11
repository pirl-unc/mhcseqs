#!/bin/bash
set -euo pipefail

echo "Running lint checks..."
bash lint.sh

echo "Running tests..."
bash test.sh

echo "Building distribution..."
pip install build twine
rm -rf dist/
python -m build

echo "Uploading to PyPI..."
python -m twine upload dist/*

echo "Done!"
