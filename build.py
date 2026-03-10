#!/usr/bin/env python3
"""Convenience shim — delegates to `python -m mhcseqs build`.

For the full CLI (build, lookup), use:
    python -m mhcseqs build
    python -m mhcseqs lookup "HLA-A*02:01"
"""
import sys
sys.argv = [sys.argv[0], "build"] + sys.argv[1:]

from mhcseqs.__main__ import main
main()
