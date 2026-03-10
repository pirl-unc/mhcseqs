#!/bin/bash
set -euo pipefail
ruff check mhcseqs/ tests/
ruff format --check mhcseqs/ tests/
