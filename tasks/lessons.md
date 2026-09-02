# Lessons

- When a source ingest fails because a shared parser wrapper changed, test the public wrapper
  directly before adding source-specific context. A source-specific workaround can hide a
  library-wide regression: `A*01:01` must retain mhcgnomes' documented human default in normal
  mhcseqs parsing, while `require_explicit_species=True` is the explicit strict mode.
- Before combining counts from a completed build with source CSV counts, establish whether the
  build already embeds those source rows. Count each final record once and use source files only
  for metadata that is not present in the build.
