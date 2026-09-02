# Issue #54: restore bare HLA parsing and prevent collapsed releases

GitHub issue: #54

## Baseline and corrected contract

- [x] Reproduce the scope of the regression. `mhcgnomes.parse("A*01:01")`
  still returns `HLA-A*01:01`, but mhcseqs 2.6.2's wrapper returns `None`
  because commit `78051fb` unconditionally supplies `default_species=None`.
- [x] Confirm this is not merely an IMGT/HLA-loader omission. The public
  `parse_allele_name()` helper stopped accepting ordinary bare HLA shorthand,
  and `_resolve_header_allele()` is only one affected caller.
- [x] Check the complete cached official sources before changing behavior.
  Removing only the mhcseqs override resolves all 45,762 IMGT/HLA headers as
  human. All 12,395 IPD-MHC names already carry a prefix; 12,019 resolve and
  the 376 existing failures are contested-prefix cases such as `Gogo`, not
  names that fall through to the human default.

## Design

- [x] Restore mhcgnomes' documented default behavior in
  `parse_allele_name()`, so bare HLA forms such as `A*01:01` work generally.
  Do not hard-code `Homo sapiens` in the IMGT/HLA ingestion path.
- [x] Preserve strict validation through the existing
  `require_explicit_species=True` option. Preserve collision rejection,
  explicit-species constraints, generated-prefix rejection, and non-molecule
  rejection.
- [x] Add regressions for the public wrapper, normalization, strict mode, and
  the real official header `HLA:HLA00001 A*01:01:01:01 365 bp`.
- [x] Gate official-source builds using minimum input and functional-output
  counts. This fails on a truncated download or parser collapse before release
  artifacts can be emitted.
- [x] Re-download the current official FASTAs, rebuild both CSV artifacts, and
  verify human/source/functional-groove counts against the last complete
  release rather than the broken 2.6.x artifacts.
- [x] Refresh README counts from the successful rebuild, bump the patch
  version, run `./lint.sh` and `./test.sh`, and inspect the final diff.
- [x] Fix and test README double-counting discovered during the refresh
  (issue #57); the completed build already contains the supplemental rows.
- [x] Fix the deploy runner's split Python/pytest environments (issue #58) and
  cover the actual shell invocation with a regression test.
- [x] Merge PR #56, publish v2.6.3 to PyPI, and independently verify its wheel,
  sdist, bare-HLA behavior, and downloadable release artifacts.
- [ ] Merge follow-up PR #60 after CI, deploy v2.6.4, and verify PyPI.

## Review

v2.6.3 is live. Its published artifacts contain 79,029 raw and 56,440 full
records. Functional coverage is 54,852 total, including 26,205 human, 33,938
class I, and 20,914 class II records. An isolated PyPI install parses bare
`A*01:01` as `HLA-A*01:01` with `Homo sapiens` from the default provenance.

The README refresh exposed issue #57: the completed build already included the
supplemental records, but the updater added their source CSV a second time.
PR #60 makes the completed build the sole count source and uses the supplemental
CSV only for its species lower bound. Deployment also exposed issue #58;
`test.sh` now detects xdist and runs pytest through the same `python` executable.
