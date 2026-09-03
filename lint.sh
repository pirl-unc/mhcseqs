#!/bin/bash
set -euo pipefail
ruff check mhcseqs/ tests/ scripts/
ruff format --check mhcseqs/ tests/ scripts/
