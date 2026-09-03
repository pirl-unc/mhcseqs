# Signal-peptide ground-truth label audit (2026-09-02)

## Scope

Issue #52 identified 241 rows in `sp_ground_truth_enriched.csv` whose MHC
class and chain were blank. Those rows were therefore sent through a separate
evaluation-only parser selector, which made signal-peptide results depend on a
23-residue prior instead of the production whole-parse score.

The audit also found 23 rows hidden from that count because their class or
chain was the literal sentinel `unknown`. They were audited under the same
rule, for 264 accession-level decisions in total.

This audit treats source annotations and parser output as independent evidence:
the parser is never used to assign a benchmark label.

## Evidence rule

Each of the 264 accessions was checked against its UniProtKB record on
2026-09-02. For entries that UniProt had deleted, the final active record was
retrieved from UniSave. A class/chain label was accepted only when the record
contained a chain-specific domain-family assignment from InterPro, Pfam,
PANTHER, SMART, or CDD. Generic protein names such as "Ig-like
domain-containing protein" were not considered sufficient.

The accession-level decision, domain identifiers, record version, and source
URL are preserved in `sp_ground_truth_label_curation.csv`.

## Results

| Decision | Count |
| --- | ---: |
| MHC class II alpha | 37 |
| MHC class II beta | 210 |
| Excluded as non-MHC/non-groove | 16 |
| Retained as source-indeterminate | 1 |

Within the 241 rows named in #52, the corresponding counts are 34 alpha, 204
beta, and 3 excluded.

The excluded records are:

- `Q08334` and `Q61190`: human and mouse IL-10 receptor beta. They entered the
  original query because their names contain "cytokine receptor class-II";
  UniProt and PANTHER assign them to the type II cytokine-receptor family, not
  MHC class II.
- `A0A401NJB8`: a catshark zinc-alpha-2-glycoprotein prediction. Its UniProt
  MHC-II keyword conflicts with the chain-specific PANTHER and InterPro calls,
  which identify an MHC-class-I-related fold rather than an MHC class II chain.
  It is therefore excluded from this MHC-chain benchmark.
- The additional sentinel-valued false inclusions comprise CRTAM, KIR2DL1,
  KIR2DL3, HMSD/serpin, five MHC-II gamma/CD74 invariant chains, and three
  minor-histocompatibility-antigen H13 signal-peptide peptidases.

`A0A8B9VS73` is retained as `unresolved`: UniProt calls the 118-residue
fragment an MHC class II antigen, but supplies no chain-specific domain-family
assignment. It therefore exercises inferred class-II chain competition and is
reported outside the curated/gold denominator.

`A0A8D0CD71` is now curated as MHC class II alpha from UniProtKB entry 21:
InterPro `IPR001003`, Pfam `PF00993`, and PANTHER `PTHR19944:SF86` agree on the
alpha-chain family. Its stored signal feature remains source data rather than a
manual correction; it is a SignalP annotation on an unreviewed sequence with
ambiguous residues.

## Evaluation policy

- `gold` and `curated` rows form the `curated/gold` denominator.
- Any future source-indeterminate row remains `unresolved/inferred` and uses
  `decompose_domains(..., mhc_class="")`, the production whole-parse
  competition.
- `excluded_non_mhc` rows remain visible in the corpus for provenance but are
  omitted from MHC benchmark denominators and negative-control generation.
- Both evaluation strata are always printed, even when one has zero rows.

During the audit, `Salvator merianae` was also found to match the fish hint
`Tor ` by substring. Issue #62 tracks that independent taxonomy bug; the same
PR applies a boundary-safe match and corrects the 27 stored tegu categories.
