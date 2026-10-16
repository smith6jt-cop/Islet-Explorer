# Protein donor expansion — 2026-10-03

> **2026-10-15:** donor 78 (HDL078 / HPAP-078) was removed from the app and all
> analysis by owner directive. The current build has five additions and 20
> donors (7 ND, 7 Aab+, 6 T1D; 7,028 islets) and is described under
> [Five-addition rebuild](#five-addition-rebuild--2026-10-15). The local app is
> offline at the owner's direction until the reprocessed root is accepted.
> Mentions of six
> additions below record the 2026-10-03 build. See the
> [directive record](owner_directive_20261015.md).

The 2026-10-14 [expanded-cohort grouping and comparability audit](data_comparability_20261014.md)
records which fields are consistently grouped and which retain method
differences between the original 15 donors and the six additions.

The owner approved the proposed phenotype rules and the six native boundary
sets, and selected Harmony on protein PCA with the uncorrected result retained.
This supersedes earlier requests to recover the complete historical phenotype
chain or train another scVI model. The original 15-donor source files remain
the frozen baseline. No scVI training or inference is required by this build.

## Current implementation

- `prepare_protein_expansion.py` reads the six raw additions, retains their
  source identities and raw intensities, applies the frozen QC cutoff, stores
  log1p intensity, and rescales each donor using the recovered marker gates.
  Unmeasured objects stay in a separate census. A gate above all observed values
  keeps the entire marker negative (6232 CD20); no extreme cell is promoted to
  positive merely because it is closest to that gate.
- `protein_phenotyping.py` applies `protein_phenotype_rules.csv`. Unresolved
  top-level conflicts and unsupported cells remain Unknown. Conflicting child
  labels retain their parent. Acinar, ECAD+ and Structural are not guessed for
  the additions. Original saved labels are not reclassified. The review reports
  disagreement with the historical labels explicitly; this is an accepted
  reconstruction, not proof of historical method equivalence.
- `export_protein_review_geometry.groovy` snapshots the reviewed native QuPath
  vertices and exact cell centroids without editing the project. Stable
  `Islet_N` aliases are assigned by native UUID, with the original names and
  UUIDs retained in the snapshot. Geometry is clipped to tissue minus exclusion;
  the 20 µm peri band excludes all cores. Native coordinate roundoff can produce
  slivers below 0.0001 µm²; these are retained unchanged, and any duplicate core
  cell membership still fails validation. Larger overlaps fail.
- `assemble_protein_expansion.py` requires the recorded owner acceptance, checks
  source hashes, independently copies original display artifacts, and appends
  the additions. Islet analysis requires a QC cell in both core and peri. Hormone
  positivity for the additions is computed in Python, not exported from QuPath:
  each donor's INS/GCG/SST threshold is the mean + 1.5 SD of measured raw
  intensity, over raw-object denominators. Phenotype summaries use QC cells. The raw composition census explicitly records unmeasured cells.
- The QC total-intensity cutoff (`QC_TOTAL_MFI` in `protein_phenotyping.py`) is
  the 1st percentile of per-cell total raw MFI over the unfiltered original
  cohort (2,676,419 cells); with the >=2-detected-markers rule it reproduces the
  frozen 2,649,654 cells. `derive_original_qc_threshold()` recomputes it. It was
  first derived by the archived joint-scVI preparation script.
- `assemble_protein_expansion.py` and `build_protein_embedding.py` default to
  `--donor-correction harmony`; with Harmony, assembly also writes the
  uncorrected fit for comparison.
- Spatial clustering standardizes the 31 mean core protein signals, takes up to
  15 PCs, applies donor Harmony (theta 2, seed 0), and fits a 30-neighbor graph.
  Leiden and UMAP use that same graph. Raw and rescaled marker values are not
  replaced with corrected values. The full uncorrected protein-PCA fit is saved
  separately, and uncorrected PCs remain in the selected H5AD's `X_pca`.
- Trajectory retains the established composition features:
  `log1p(count / max(core endocrine count, 1))`, including residual Endocrine.
  The default minimum is 10 cells. Both cell-range bounds apply only when the
  user presses **Recalculate trajectory**. No disease-order acceptance target
  is imposed.

## Locations and provenance

| Purpose | Location relative to repository |
|---|---|
| Frozen original baseline | `data/` |
| Immutable original cell partition | `single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad` |
| Raw additions | `data/development/donor-additions-20260923-raw-v1/additions/raw/` |
| Accepted prepared inputs | `data/development/protein-additions-20261015-v3/` (v2 superseded) |
| Native geometry snapshot | `data/development/protein-review-geometry-20261015-v3/` (v2 superseded) |
| Owner acceptance and exact hashes | `tmp/protein-expansion-20261015/owner_acceptance.json` (derived from the 2026-10-03 record) |
| Cohort manifest | `scripts/protein_expansion_cohort.csv` (20 donors) |
| Excluded donors | `scripts/excluded_donors.csv` |
| Expanded app root | `data/development/protein-expanded-20261015-v3/` (v2, with donor 78, superseded and refused by the app) |
| Build diagnostics and uncorrected fit | `tmp/protein-expansion-20261015/expanded-v3/` |
| Native membership validation | `tmp/protein-expansion-20261003/geometry-validated/validation.json` |

The v1 prepared inputs and geometry-v1 snapshot are superseded. The first expanded build failed during H5AD serialization and was never served; expanded-v2 is the corrected build. Geometry-v2's
manifest records the review state at export time; the later owner acceptance
JSON is authoritative. The explicit cohort manifest selects recovered 6090 and
corrected 6539. Images are linked to their exact source paths and inventoried;
the assembly does not duplicate image pixels or substitute the earlier 6539.

Original single-cell preservation is established by retaining the complete
source partition and verifying its SHA-256 before and after assembly, rather
than rebuilding its millions of rows. The canonical islet parity audit aligns
IDs at zero numeric tolerance and permits only documented derived refits.
Embedded table row counts are checked against their full appended tables;
every original table row remains independently compared. The separate release
validator checks original display files byte-for-byte, original CSV prefixes,
all added raw cell values, composition formulas and corrected/uncorrected fits.
Each of the six additions also carries a complete raw-cell census parquet with
its inclusion/QC status, including objects without measured markers. The
expanded-v2 release validator passed, and the browser acceptance run passed all
45 checks on 2026-10-05 (report: `tmp/protein-expansion-20261003/expanded-browser-acceptance.json`). Automated accessibility checks found no violations; color-contrast findings on plots remain for manual review.

## Five-addition rebuild — 2026-10-15

The owner directed that donor HDL078 (case 78) be removed completely from the
app and all analysis. `stage_original_gate_calls.DONORS` now lists 6090, 6232,
6520, 6539 and 6558. Preparation, assembly, validation and the browser check
reject any donor named in `scripts/excluded_donors.csv`. Commands, run from the
repository root with the scVI environment's Python:

```bash
PY=/home/smith6jt/miniconda3/envs/scvi-env/bin/python
"$PY" scripts/prepare_protein_expansion.py \
  --frozen single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad \
  --raw-root data/development/donor-additions-20260923-raw-v1/additions/raw \
  --output data/development/protein-additions-20261015-v3 \
  --report-root tmp/protein-expansion-20261015/preparation-v3
"$PY" scripts/assemble_protein_expansion.py --baseline data \
  --prepared data/development/protein-additions-20261015-v3 \
  --geometry data/development/protein-review-geometry-20261015-v3 \
  --acceptance tmp/protein-expansion-20261015/owner_acceptance.json \
  --cohort scripts/protein_expansion_cohort.csv \
  --clinical data/development/expanded-20260911/metadata/selected_donor_metadata.csv \
  --output data/development/protein-expanded-20261015-v3 \
  --report-root tmp/protein-expansion-20261015/expanded-v3 --donor-correction harmony
"$PY" scripts/validate_protein_release.py --baseline data \
  --candidate data/development/protein-expanded-20261015-v3 \
  --raw-root data/development/donor-additions-20260923-raw-v1/additions/raw \
  --comparison tmp/protein-expansion-20261015/expanded-v3 \
  --report tmp/protein-expansion-20261015/expanded-v3/release_validation.json
```

- The five prepared donor files are byte-identical to v2 (same SHA-256), as
  expected for per-donor preparation.
- The geometry snapshot is a five-donor copy of v2; geometry, centroid and
  native-project hashes are unchanged.
- The acceptance record carries the 2026-10-03 rule, embedding and boundary
  approvals over unchanged for the five donors and records the removal. No new
  review was performed.
- Assembly produced 7,028 paired islets (5,214 original + 443 6090, 137 6232,
  428 6520, 211 6539, 595 6558). The Harmony protein PCA, UMAP and Leiden fits
  were refit without donor 78, so derived coordinates and clusters differ from
  v2. Original-subset parity and the release validator passed: 5,277 original
  files byte-identical, 20 images, no excluded donor in any table.
- A scan of every table, H5AD (including embedded QuPath tables), workbook,
  file name, image link and acceptance/boundary record found no donor-78 entry.
- Browser check on the live app (`tmp/protein-expansion-20261015/expanded-v3-browser-acceptance.json`):
  41 of 45 checks passed, including the 20-donor census, absence of 78 and the
  five additions' tissue and image checks. The four failures do not involve
  donor 78. Three come from the 2026-10-14 app change (`b1ce947`), which removed
  the hormone fractions from Plot and made the PAGA layout build only on request;
  the harness was not updated for it. The fourth is a Trajectory color-contrast
  finding.

## Local service and rollback

*Since 2026-10-15 the service is stopped and disabled and the static image route
is closed, at the owner's direction; do not restart it before the reprocessed
root is accepted. The description below is how it ran.*

The user service `islet-explorer-local.service` runs Shiny on loopback port 3840;
the existing nginx endpoint is `http://10.15.152.7:8080/islet-explorer/`.
Its checked-in definition is `deploy/islet-explorer-local.service` and config is
`deploy/shiny-server-local.conf`. The service is enabled with user lingering.
Shiny Server forwards only a minimal environment to R. `deploy/run-local-r.sh` loads the selected root and runtime settings for each R worker; `R` in the service override selects this wrapper. The data-root override is `tmp/local-app.env`; after changing it, restart with
`systemctl --user restart islet-explorer-local.service`.

To restore the frozen 15-donor data, set
`ISLET_DATA_ROOT=/data/panc_CODEX/Islet-Explorer/data` in that override and restart.
No source data needs to be rewritten. This is the local research app, not a new
UFRC production deployment. Existing production release requirements still
apply to any later external publication.
