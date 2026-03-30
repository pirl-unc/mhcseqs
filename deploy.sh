#!/bin/bash
set -euo pipefail

VERSION=$(python -c "from mhcseqs.version import __version__; print(__version__)")
echo "=== Deploying mhcseqs $VERSION ==="

# Verify clean state
if [[ -n "$(git status --porcelain)" ]]; then
    echo "ERROR: working tree is not clean"
    git status --short
    exit 1
fi

# Verify tag exists
if ! git tag -l "v$VERSION" | grep -q .; then
    echo "ERROR: tag v$VERSION does not exist. Run: git tag v$VERSION && git push origin v$VERSION"
    exit 1
fi

echo "Running lint checks..."
bash lint.sh

echo "Running tests..."
bash test.sh

echo "Building distribution..."
rm -rf dist/
# Use pyproject-build instead of `python -m build` to avoid conflict
# with mhcseqs.__main__ shadowing the build module.
pyproject-build

echo "Checking distribution..."
python -m twine check dist/*

echo ""
echo "Ready to upload mhcseqs $VERSION to PyPI."
echo "Distributions:"
ls -la dist/
echo ""
read -p "Upload to PyPI? [y/N] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    python -m twine upload dist/*
    echo "=== mhcseqs $VERSION deployed to PyPI ==="
else
    echo "Aborted. Distributions are in dist/ if you want to upload manually."
fi
