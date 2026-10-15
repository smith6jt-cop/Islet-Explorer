# Donor expansion, phase 2 — 2026-09-23, revised 2026-09-29

> **2026-10-03 owner decision:** the accepted reconstruction rules and all six donor boundaries are approved. Use Harmony on protein PCA, retain the uncorrected comparison, and preserve original source records. Exact historical typing/encoder recovery and scVI retraining are no longer prerequisites. See the [current implementation record](protein_expansion_20261003.md); earlier directions below are historical.

The September 23 raw-only checkpoint remains the starting evidence. The owner's
September 29 direction replaces its projection-only and strict development
loader constraints. The regenerated-all-21 September 17–22 candidates remain
unaccepted experiments. The restored app from `6b7884b` and the preservation
gate from `70f8171` provide the baseline; the gate will need to distinguish
preserved source records from derived outputs allowed to change in a joint refit.

## Artifact and method contract

- Preserve the original 15 donors' identities and 5,214 islets, cell/islet IDs,
  physical geometry, measured values, missingness and phenotype labels exactly.
  Keep their baseline artifacts intact for comparison. A common 21-donor
  trajectory, UMAP and Leiden refit may change those derived coordinates,
  graph edges and cluster assignments for the original donors; report such
  changes separately from source-record parity.
- Stage donors 78 (HPAP-078), 6090, 6232, 6520, 6539 and 6558 separately.
  Select recovered 6090 and corrected 6539 from `cohort-v1.csv`. The expected
  raw universe is 943,068 objects, including 201 wholly unmeasured objects.
  Preserve source IDs, raw intensities and missingness without imputation.
- A staging bundle contains independent copies of the four original canonical
  artifacts at its root and six raw H5ADs under `additions/raw/`. It is an
  offline evidence bundle, not an expanded app release. Its original-only
  canonical files must pass the zero-tolerance parity gate. That pass makes
  no claim that the additions are scientifically comparable or app-ready.
- Raw H5AD `X` is raw mean intensity. QuPath `Classification` is source
  annotation. Neither substitutes for legacy normalized expression,
  `phenotype`, scVI, islet composition or trajectory coordinates. Missing
  methods remain explicitly unavailable in a donor availability table.
- *(Historical; superseded 2026-10-03.)* Recover and verify the exact phenotype rules and postprocessing used for the
  original donors. The original 31-marker numeric gates and input scale have
  now been recovered from the frozen H5AD, but the checked-in rules table does
  not reproduce its labels. To obtain comparable scVI coordinates, either
  locate the original model weights and preprocessing or explicitly review a
  joint retrain with the notebook method. Do not retune gates for the new
  donors or infer that the current UI rule CSV reproduces the full original
  typing method. The proposed enhancement using CHGA, CD45,
  VIM, LGALS3, CD44, PGP9.5, GAP43, CD35 and NaKATPase is deferred.
- Regenerated original-donor geometry, typing and coordinates from schema 4
  remain experimental evidence. Their 11,944 regions / 11,136 eligible islets
  are not accepted expansion counts. No combined islet count or common fitted
  axis has yet been accepted.
- The owner will visually review the new-donor annotations before scientific
  acceptance. The earlier classifier-eligibility and identity-review rules
  remain experimental audit behavior, not substitutes for the original numeric
  gates. Account for every requested donor and fail explicitly when required
  measurement or training support is insufficient.
- Once the six additions have comparable typed measurements, a common
  21-donor trajectory, UMAP and Leiden refit is allowed. Preserve the original
  donor source records and validate the resulting science; the experimental
  schema-4 diffusion axis is not an accepted refit.

## Checkpoints and acceptance

1. Reconcile the contract and stage a fresh raw-only bundle, recording source
   and copied-file hashes, per-donor census and missing capabilities.
2. Repair experimental eligibility, checkpoint dependencies, donor retention,
   display-rank validation and fresh-output protections. Verify those repairs
   on focused fixtures without running the all-donor regeneration pipeline.
3. Independently verify the staged bundle and run `audit_original_parity.py`.
   Preserve checkpoint reports outside app data; never overwrite a prior run.
4. *(Historical; superseded 2026-10-03 — boundaries approved, accepted
   reconstruction rules and Harmony protein PCA used instead.)* Have the owner visually review donor annotations. Reproduce frozen labels
   with the exact phenotype rule chain and recovered numeric gates; resolve
   original-model projection versus a reviewed joint scVI retrain. Process the
   six additions with the selected methods and audit original donor IDs,
   measurements and labels separately from derived coordinates.
5. Fit and validate a common 21-donor trajectory, UMAP and Leiden result from
   comparable inputs. Passing software tests or raw staging does not certify
   scientific acceptance.
6. Assemble a separate, simple development data root and point the local app to
   it with `ISLET_DATA_ROOT` after validation. No additional candidate-version
   manifest, source-selection gate, or optional loader flag is required for
   development. Keep the original baseline artifacts intact; production release
   gates are separate.

Historical transfer-screen evidence is
`data/development/expanded-20260911/reports/phenotyping_method_audit.json`:
the original initial-lineage rule matrix was not recovered, and neither tested
replacement met its prespecified donor-disjoint screen. No substitute method is
accepted by this phase-2 implementation.

## Completed checkpoint

The offline bundle is
`data/development/donor-additions-20260923-raw-v1/`. It has no runtime
`candidate.json`, no expanded canonical table, and no app loader. All copied
files have independent inodes; source hashes are checked against import
sidecars and again against the copies. Source geometry is archived for review;
images and regional exports remain explicit external references and are not
revalidated or promoted by this raw-cell checkpoint.

| Addition | Raw cells |
|---|---:|
| HPAP-078 | 267,993 |
| 6090 | 170,566 |
| 6232 | 92,012 |
| 6520 | 102,526 |
| 6539 (recovered) | 160,837 |
| 6558 | 149,134 |
| Total | 943,068 |

There are 942,867 fully measured cells, 201 wholly unmeasured objects and no
partially measured objects. All six have unavailable legacy phenotype,
normalized-expression, scVI, islet-composition and compatible trajectory
capabilities, with reasons in `additions/availability.json`.

The repaired experimental pipeline now:

- retains eligibility, measurement-support and training-support fields through
  identity review and records the review separately;
- fingerprints training masks, reference coordinates, parameters, software,
  deterministic seed and code before reusing validation checkpoints;
- checks the declared cohort before fitting and rejects any donor lacking
  required training support instead of filtering it away;
- verifies each display rank against raw-coordinate order and the SHA-256
  order of exact ties, including missing and duplicate identifier checks;
- requires a new trajectory output directory separate from the source,
  rejects recursive staging and broken output aliases, and copies derived
  schema-4 revision artifacts instead of hardlinking mutable files.

Verification: **91 focused Python tests pass**, including an isolated complete
trajectory run on synthetic data. The real bundle passes **852 zero-tolerance
original-preservation checks** and exact raw source/copy comparisons. The four
canonical source hashes also match the pre-restoration record. Evidence is in
`tmp/donor-additions-20260923-raw-v1-validation/` (tests, initial validation,
frozen hashes) and `tmp/donor-additions-20260923-final-validation/validation.json`
(final independent validation). No original data, app code, service or profile
was changed. Pre-existing untracked research files remain untouched.

As of 2026-09-29 the integration checkpoint was unmet: the original numeric gates were
recovered, but the complete phenotype rule chain and scVI encoder were not. The
2026-10-03 owner decision removed that prerequisite; see
[the protein expansion record](protein_expansion_20261003.md).
The searched locations, replay results and precise resume
requirements are in [the method recovery record](original_method_recovery_20260929.md).
The earlier donor-disjoint
transfer screens reported raw-feature accuracy/macro-F1 0.86668/0.8380 and
rescaled-feature accuracy/macro-F1 0.87988/0.8526 against required minima
0.97/0.90 (plus per-label recall requirements). These are historical failure
records, not new validation. Their proposed labels must not enter the bundle.
No common 21-donor refit or six-donor app integration was performed, and no
scientific acceptance is claimed. The saved 852-check parity result applies to
the raw-only baseline copies; it does not validate a future joint refit.

## Reproduction

Run from the repository root with the installed analysis interpreter. Both
output directories must be new; a partial or prior run is not overwritten.

```bash
PY=/home/smith6jt/miniconda3/envs/scvi-env/bin/python
BUNDLE=data/development/donor-additions-next
REPORT=tmp/donor-additions-next-validation

"$PY" scripts/research/stage_donor_additions.py \
  --cohort scripts/research/cohort-v1.csv --input-root . \
  --baseline-root data --output-root "$BUNDLE"
"$PY" scripts/research/validate_donor_additions.py \
  --baseline-root data --bundle-root "$BUNDLE" --report-root "$REPORT"
```

The validator performs the read-only `audit_original_parity.audit` comparison
at its default zero tolerance after checking every raw donor and source hash.
For a standalone original-only comparison, use `audit_original_parity.py` with
the same baseline and bundle roots and a separate fresh report directory.
There is deliberately no preview or promotion step for a raw-only bundle.

```bash
"$PY" -m pytest -q \
  scripts/tests/test_donor_additions.py \
  scripts/tests/test_phase2_trajectory_guards.py \
  scripts/tests/test_research_pipeline.py \
  scripts/tests/test_runtime_schema4.py \
  scripts/tests/test_composition_trajectory.py \
  scripts/tests/test_original_parity.py \
  scripts/tests/test_preview_preservation.py \
  scripts/tests/test_reprocessing_contracts.py
```
