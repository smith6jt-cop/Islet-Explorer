# Expanded cohort preparation status

Updated 2026-09-16. All six additions have completed offline segmentation,
measurement export, raw import and QC preparation. The audits validate 943,068
cell objects and independently reproduce all 2,989 islet-core counts. The
6090/6539 source mix-up is resolved from the local raw acquisitions. The original
phenotype method remains unrecovered; the revised task instead audits retained
labels and their upstream gates under the retirement policy. A separate
21-donor runtime and fixed outcome-blind trajectory have now been built, with
six provisional primary results and four unavailable endpoints. The local
preview is http://127.0.0.1:3839/islet-explorer/. It contains 20 nPOD and one
HPAP donor, retaining the frozen 15-image `5809be4` release. See the
[all-donor delivery record](all_donor_delivery.md) for complete counts,
checksums, validation, restart commands and remaining scientific work.
The preparation history below remains evidence for its named processing stages.

Processing resumed at the user's request after the 2026-09-12 pause. Completed
stage-script hashes were verified and the project was backed up before restart.
The pause and resume records remain in
`data/development/expanded-20260911/processing_pause.json` and
`processing_resume_20260916.json`.

## Corrected source paths

6090's existing stack matches every pixel of all 36 corresponding raw channel
planes in `6090AAbtest31_35p_led/`. Its correct processing source is:

```text
/data/QuPath_Projects/QuPath_project_IsletExplorer/Images/New/6090AAbtest31_35p_led.tif
```

The file originally named `Images/New/PAN-6539-ND test 34.tif` duplicates 6090.
The distinct 6539 acquisition was recovered from the 60 raw TIFF planes in
`/data/QuPath_Projects/QuPath_project_IsletExplorer/2022-11-11  PAN-6539-ND test 34/`.
The 36 selected channels were assembled without changing any pixels and checked
back against every original plane. Its corrected processing source is:

```text
/data/panc_CODEX/Islet-Explorer/data/development/expanded-20260911/recovery/sources/PAN-6539-ND test 34.tif
```

| Source | Shape, channels × height × width | SHA-256 |
|---|---|---|
| Verified 6090 | 36 × 9,072 × 12,096 | `81b5a874aebc828b5c6c5925d79736b17497bf286e0aefcdf9cebedc9275c2f8` |
| Recovered 6539 | 36 × 9,504 × 12,672 | `1866d5e8856dc096c1829e04b3622bc796cd9b81be4225bff12b864d4df1c69e` |

The corrected 6539 viewer is
`data/development/expanded-20260911/recovery/images/ND_6539.ome.tiff`.
The five other accepted viewer TIFFs remain in the candidate's main `images/`
directory. The older `images/ND_6539.ome.tiff` is rejected evidence; never stage
that directory wholesale. The corrected six TIFFs total 46,877,746,573 bytes.

The old 6539 exclusion annotation was drawn on the duplicated 6090 pixels. It
is preserved in the original project and in `recovery/historical_annotations/38/`,
and was not transferred onto the distinct recovered image. Recovered 6539 starts
with fresh annotations and needs its own exclusion/boundary review. Its channel
display ranges were derived from its own pixels for SAM/rendering only; raw
intensities and the original acquisition calibration are preserved.

## Measurements and imports

All six acquisitions completed the ten processing/export/QC stages and passed
the raw-data roundtrip, regional-join and independent spatial-count audits.
Recovered 6090/6539 results remain in the separate `recovery/` project.

| Addition | Cell objects | Fully measured | Unmeasured objects | Islet regions |
|---|---:|---:|---:|---:|
| 6232 | 92,012 | 91,878 | 134 | 128 |
| Selected PAN-6558 | 149,134 | 149,118 | 16 | 200 |
| 6520-01A | 102,526 | 102,521 | 5 | 321 |
| HPAP-078 | 267,993 | 267,978 | 15 | 1,234 |
| Verified 6090 | 170,566 | 170,556 | 10 | 868 |
| Recovered 6539 | 160,837 | 160,816 | 21 | 238 |
| Six completed additions | 943,068 | 942,867 | 201 | 2,989 |

The exports contain 2,989 composition rows, 26,871 regional marker rows and
17,988 regional target rows. These describe different tables and are not
interchangeable counts. Raw H5AD files are in `analysis/raw_cells/` for the first
four cases and `recovery/analysis/raw_cells/` for 6090 and 6539, under the candidate
root. All six have representative cell/nucleus review figures in their project's
`reports/segmentation/` directory; these support review without establishing
whole-image biological acceptance.

InstanSeg uses model 0.1.1 and the established 22 segmentation channels.
Structural measurements are scoped to Tissue annotations and preserve cell
UUIDs and parents. Imports preserve the original ordered 31-marker panel,
raw mean intensities as float32, a `log1p` layer, and float64 micrometer coordinates.
HPAP-078 remains an
HPAP source with analysis ID `78` and viewer/export directory ID `0078`.

The default QuPath exporter rounded coordinates enough to alter reconstructed
regional membership. All six measurement exports now use Java double precision
and the importer parses decimal values with round-trip fidelity. The five earlier
rounded exports/imports, their checkpoints and the previous script versions are
preserved in `precision_revision_20260916/`. Segmentation and cell identities were
retained. The independent audit exactly matches all 2,989 core counts. One cell
lies exactly on a boundary; its inclusion follows QuPath 0.7 `PolygonROI`'s
winding rule, verified directly against the saved project and installed software.

An earlier overlap between QuPath jobs overwrote 6539's shared project channel
metadata. It was restored from the intact normalization checkpoint after checking
the image source, dimensions, channels and calibration. The runner now locks each
project and processes requested cases sequentially. The repair and backup are
recorded in `recovery/channel_metadata_repair_20260916.json`.

Unmeasured objects have all 31 marker means missing. The archival import retains
them and their coordinates with `marker_measurements_complete=false`; it does
not impute values, discard cells, or assign phenotypes. Their handling must be
resolved before training. The exported `Parent` field remains the actual QuPath
parent; regional membership is checked against named islet geometry, not inferred
from an invented parent string. Islet-region counts are preparation counts,
not final retained islets after joint scientific filtering. The audit records
2,039 cells with membership in multiple cores and does not invent a single-islet
assignment. Segmentation/exclusion review, overlapping regions and the 201
unmeasured objects need a documented handling policy before joint analysis.

## Phenotype-method investigation

Searches of local acquisition/analysis directories, hidden files, notebooks,
alternate rule files, and repository history did not recover the original
21-class rule matrix. The current `data/phenotype_rules.csv` cannot produce all
21 labels in the frozen analysis. The historical
`7571b6d:scripts/phenotype_evaluator.py` describes its v2 method as a reconstruction
and records that the original ECAD-positive/Structural rules were unavailable.
Other local rule files use different panels or methods.

Two candidate marker-based transfer models were screened against frozen labels
using five donor-disjoint folds, 5,000 held-out cells per donor, and capped
training samples per donor/label. No new-cell labels were written.

| Marker feature space | Accuracy | Macro F1 | Automatic integration screen |
|---|---:|---:|---|
| Raw `log1p` marker means | 86.67% | 0.8380 | Failed |
| Frozen rescaled marker matrix | 87.99% | 0.8526 | Failed |

The predeclared screen required accuracy ≥97%, macro F1 ≥0.90, and recall ≥0.85
for labels with at least 50 held-out examples. These limited experiments reject
the tested replacements; they neither reconstruct the original method nor prove
that no suitable alternative exists. These historical benchmarks are not
required to reconstruct retired labels. The revised candidate uses raw marker
summaries for diffusion ordering and preserves uncertainty in retained subtype
assignments. No scVI/Harmony refit or disease-order tuning defines the new axis.
Biological gate review remains necessary for unavailable phenotype-dependent
analyses, independently of the completed local runtime.

## Evidence and migration implications

Restricted artifacts are under `data/development/expanded-20260911/`:

- `processing_completion_20260916.json`: final six-case preparation totals,
  evidence hashes, validation results and remaining requirements.
- `reports/image_conversion/source_recovery_6090.json`: all-plane raw-source match.
- `recovery/source_reconstruction_6539.json`: source hashes, channel selection,
  calibration, reconstruction and pixel checks.
- `recovery/recovered_project_provenance.json`: source reconciliation and preserved
  historical annotations.
- `reports/image_conversion/reconciled_inventory.json`: corrected six viewer paths,
  current hashes and distinct source/pixel signatures. It supersedes the initial
  duplicate audit for source selection; the initial report remains preserved.
- `reports/preparation/summary.json` and `recovery/reports/preparation/summary.json`:
  all six cases' measurements, raw roundtrips, exact spatial counts, regional
  joins, stage checksums, QC images and unmeasured-object inventories.
- `reports/preparation/spatial_predicate_evidence/`: installed QuPath boundary
  semantics and actual ROI-type checks for the six acquisitions.
- `reports/preparation/capacity_inputs_20260916.json`: measured viewer/frozen
  runtime sizes and the provisional storage lower bound.
- `reports/phenotype_transfer/` and `reports/phenotype_transfer_rescaled/`: benchmark
  design, fold results, confusion matrices and rejection criteria.
- `reports/implementation/`: local validation, hosted CI failure diagnosis and
  dated read-only UFRC correspondence references.

All 34 targeted metadata, import/audit and project-lock tests pass with zero
skips. Both workflow YAML files and their 30 shell blocks parse successfully.
These are the preparation-stage local checks. Subsequent candidate-specific
tests and browser evidence are recorded in the [delivery record](all_donor_delivery.md);
updated hosted CI and an immutable production release remain pending.

Processing commands and restart rules are in the
[workflow README](../scripts/cohort_expansion/README.md). Large images,
measurements and restricted correspondence remain outside Git.

The release manifest now declares its reviewed image count and retains support
for the original 15-image rollback release. Synthetic 21-image activation and
15-image rollback tests pass. The later `000e970` hosted jobs failed during
package installation because the runner's unrelated Chrome repository returned
a hash mismatch. Local workflow changes restrict package installation to Ubuntu
repositories and preserve initial evidence; a hosted rerun remains pending.

UFRC request `83810`, configuration `SYSCONFIG-8173`, and assessment `IRM-456246`
remain the existing records. On 2026-09-16 an assigned ISM requested a brief
Microsoft Teams meeting for the risk review. Scheduling/completion is not
established. Comment 8, received 2026-09-11, had reported the ISO waiting for
that local review. Comment 7 requested the
UFRC-designated technical contact be added to the assessment; completion is not
established. Assigned host/account/SSH/proxy details and production acceptance
remain pending. No ticket or email was sent during this continuation. See the
[migration handoff](../PUBAPPS_UFRC_HANDOFF.md) and
[readiness index](ufrc_readiness.md).
