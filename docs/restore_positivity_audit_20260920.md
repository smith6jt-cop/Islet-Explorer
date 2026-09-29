# Fresh marker positivity and historical divisor audit

Recorded 2026-09-20. This concerns marker-positive heatmap numerators and direct
marker-count endpoints. Hierarchical cell typing continues to use the corrected
continuous 31-marker measurements and reviewed cluster identities. RESTORE
positivity is not substituted for that classifier.

## Fresh computation

The candidate is `data/development/all-donor-20260920`. The independent output
`positivity_restore_v1/` contains 630 fresh donor/pair computations for all 21
donors, using 180,120 sampled QC-eligible native cells.

`scripts/research/restore_positivity.py` calls the vendored locked v8 scorer
without replacing its negative-control maximum. It requires:

- Fresh raw and REDSEA-corrected whole-cell means joined by the exact native
  detection UUID. It rejects differing raw/corrected object universes.
- Genuine native QuPath Membrane and Cell means. The read-only exporter
  `scripts/cohort_expansion/export_restore_compartments.groovy` reads existing
  hierarchies without modifying the source project or opening image pixels.
  Every native Cell mean was checked against the selected candidate's raw mean
  for all 31 markers. Source, sample, export, quality-policy and code hashes are
  saved.
- The same area and technical-region eligibility as cell typing: at least
  5 µm² and outside confirmed technical exclusion regions. Raw objects remain
  retained. No disease status, clinical outcome or historical phenotype enters
  control fitting or sampling.
- A deterministic nominal 1/20 sample using SHA256 of
  `20260920:case_id:UUID`. The older sample also selected nominally 1/20, but
  used DuckDB `hash(native UUID) % 20`. It did not use a 20,000-cell cap.
- The existing declared marker-specific clipping and raw floor methods,
  Membrane+Cell NNMF, three locked seeds and control diagnostics. Raw fitting
  floors are not production positivity thresholds. The fresh divisor is the
  **full maximum corrected target value in the inferred reference-control
  population**; positivity is strictly `corrected value > divisor`.

The fresh screen has 514 computationally supported donor/pairs with historical
biological reference acceptance, 75 reference-review-pending comparisons and
41 failed or unresolved computations. All retain a separate current calibration
review status. These counts do **not** report 514 newly image-validated gates.

Current application of calls requires an explicit donor-by-marker review,
matching evaluation hash, current image evidence and unchanged corrected input.
Numerically failed pairs cannot be promoted by an image-review switch. Missing
calls remain nullable, including in empty compartments. BActin, NaKATPase,
BCatenin and ColIV have no declared mutually exclusive RESTORE reference pair;
their continuous expression remains available while RESTORE positivity remains
unavailable. No reference pair is invented for them.

## Historical acceptance was reconciled

The first-pass CSV documents accepted ECAD←CD31, CD31←ECAD, CD45←ECAD,
PGP9.5←panCK, VIM←ECAD and CHGA←CD31 biological reference choices for the original
15 donors. The later user review recorded in
`scripts/codex_phenocycler/docs/HANDOFF.md` accepted the level-2 baseline pair web,
with named cautions about CD20 interpretation, sparse targets and problematic
donor channels; PDPN's clipped recheck remained pending. The archived level-2
CSV status still says pending, so the handoff must be read alongside the CSV.
The old immune-only counterpart set was superseded. Hormone/panCK pairs remain
comparators rather than silently replacing the approved hormone references.

Biological pair acceptance is retained as reference evidence. It is neither a
copied numerical threshold nor proof of appropriate current controls in the six
added donors. CD20 positivity is an assay measurement; it is not presented as a
B-cell identity outside the validated immune context.

## A consequential historical numerical provenance mismatch

The original15 fresh-to-archived median divisor ratios were 7.52 for CD3e and
7.86 for CD8a, compared with 1.03 for CD45 and 0.90 for CHGA. Sampling/QC alone
does not explain the largest discrepancy.

Re-evaluation of the **exact stored historical sample and historical REDSEA
values**, using the locked scorer, demonstrates that many archived level-2
`candidate_full_maximum` values are the controls' **99.5th percentile**, not
their maximum. On donor6516:

| Target | Control cells | True maximum | Control q99.5 | Archived “maximum” |
|---|---:|---:|---:|---:|
| CD3e | 3,011 | 1,249.548391 | 57.408935 | 57.408935 |
| CD8a | 3,060 | 4,793.654517 | 360.221192 | 360.221192 |

Across the original15, the exact historical-sample audit found this q99.5 match
for 14/15 CD3e, 13/15 CD8a and 2/15 CD45 computations. CHGA retained the true
maximum in all 15. The q99.5 matches occur where the preserved rejected patch's
`max > 5*q99` condition triggers. These exact numerical matches support the
inference that the archived level-2 artifacts were computed with the unversioned
adaptive patch despite their v8 labels. The patch was explicitly rejected and
preserved in `scripts/codex_phenocycler/docs/rejected/adaptive_restore_divisor.patch`.
The current implementation does not reinstate it.

Accordingly, the archived level-2 numerical divisors are withdrawn as a source
for current calibration; their biological reference-review history is a
separate record. Original files are preserved. The exact checks are saved in
`positivity_restore_v1/sensitivity/sampling_qc_sensitivity.csv`; the diagnostic
runner is `scripts/research/restore_sampling_sensitivity.py`. The first-pass
artifact generation and remaining pairs are audited separately rather than
assuming all historical generations share the same defect.

## Current image evidence and limits

The endocrine audit contains all21 donor strips and seven-donor pages for INS,
GCG, SST and CHGA in `positivity_restore_v1/endocrine_audit/`. The immune audit
similarly covers CD45, CD3e, CD8a, CD4, CD20, CD68, CD163, HLADR, CD35, CD44 and
LGALS3. Each strip shows the brightest inferred reference control, a cell just
above the current divisor and a higher-target cell, with DAPI, target,
reference and the exact source-cell boundary. Sidecars bind image, evaluation,
source and corrected-measurement hashes. The audit renderer is
`scripts/research/review_restore_panels.py`; fuller individual panels and
scatter diagnostics use `render_restore_positivity.py`.

Rendered evidence alone is not acceptance. Focal target staining or overlap
with a neighboring target-positive cell in a nominal negative control is a
concrete reason to hold a gate even when the numerical scorer reports VALID.
Several reviewed CD3e/CD8a examples show this limitation. A supported gate is
an operational background cutoff with uncertainty near the threshold, not
proof that every threshold-positive object is a biologically pure cell type.
No parameters are changed merely to produce desired counts or disease trends.
