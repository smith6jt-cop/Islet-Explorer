# Data Provenance — Islet Explorer

**Current selection, 2026-11-02:** the owner selected frozen v19 endocrine calls,
with v19/v20 compared by component. Both roots share 3,199,606 eligible UUIDs and
the same Tissue/cell/endocrine geometry. The DG Lymphatic-only rerun unintentionally
changed shared calibration; corrected code now binds calibration independently
and verifies all markers and identities. Diagnostic replays are evidence, not
assembled app roots. No corrected mixed candidate has been accepted or deployed.
See [the component comparison](../docs/v19_v20_component_comparison_20261102.md),
[workflow guide](../docs/workflow_user_guide.md), and [handoff](../docs/handoff_20261102.md).
Older lineage sections below are retained as historical records.

> **2026-10-15 owner directive:** donor 78 (HDL078 / HPAP-078) is removed from the app and all analysis ([record](../docs/owner_directive_20261015.md); `scripts/excluded_donors.csv`). Rows and counts below that include HPAP-078 describe historical staging. The [recovery record](../docs/owner_directive_20261015.md#original-15-donor-files) lists which files behind the original 15 donors were located.

This document describes the canonical data lineage for the Islet Explorer app.

## Revised local candidate — 2026-09-17

`data/development/all-donor-20260917/` is a separate schema-3 revision of the
candidate below. Raw inputs/images are immutable hardlinks; derived files are
independent copies. Nine fitted trajectory artifacts remain byte-identical.
New count/evidence artifacts use all-donor parent-gate auditing independently
of archival labels. Unreviewed gates remain unvalidated and ratios unavailable.
The owner's QuPath review also identified suspect islet candidates: 348 HPAP
regions have no cells above exported hormone thresholds, including 143 in the
previous fit. The fit is now explicitly provisional pending identity review.
No review-project objects were modified. Exact analysis-region overlays and
raw-channel contact sheets are stored separately under `data/reviews/`.

Research results are copied without refitting into
`data/reports/all-donor-research/20260917/`, outside the runtime. The revised
manifest declares 13 artifacts and 164,147,382,976 runtime bytes; scientific
and UFRC acceptance remain false. See the
[revision and hash record](../docs/trajectory_revision_20260917.md).
The dated build history below is preserved for provenance.

## All-donor candidate — completed local build 2026-09-16

The separate `data/development/all-donor-20260916/` candidate uses
[`scripts/research/cohort-v1.csv`](../scripts/research/cohort-v1.csv) as its sole
21-donor input selection: seven ND, seven Aab-positive and seven T1D, including
HPAP-078 with HPAP provenance. It combines original raw H5AD marker means,
six corrected additions, geometry, structural masks, regional exports and
clinical metadata. The frozen `5809be4` release remains unchanged.

All 3,592,722 raw objects are retained; 201 wholly unmeasured objects remain
NaN and are explicitly excluded from marker-dependent analyses. New compartment
summaries cover 8,371 islets; 3,499 satisfy the fixed trajectory eligibility
criteria. The retired phenotype policy replaces the earlier goal of rebuilding
all 21 historical labels. Retained assignments require complete upstream-gate
validation; unresolved assignments and ratio denominators remain unavailable.

The new trajectory is fit offline from core/peri median `log1p` raw intensities,
excluding INS/GCG/SST/CHGA, clinical outcomes, phenotype features and old
embeddings. Disease labels and endocrine ratios cannot select or orient its
root. Its transformations, landmarks, root, graph, coordinates and validation
are saved under `trajectory/`. Prior disease-oriented code is archived under
`archive/trajectory-pre-20260916/`; no new scVI/Harmony fit defines this axis.

Schema-2 `manifest.json` records artifact/source SHA-256 hashes and the
21-image contract. The app uses these candidate files without fallback to
original-cohort artifacts. `research/` contains the registry, six provisional
primary estimates, four explicit unavailable tests and exploratory branches.
Raw inputs/images are locally hardlinked for space efficiency and must never
be modified. Production copies require 164,130,834,281 bytes for the declared
runtime. Scientific and UFRC acceptance remain false. Complete lineage,
validation limits and local operation are in the
[all-donor delivery record](../docs/all_donor_delivery.md).

## Expanded cohort preparation history — 2026-09-16

The new candidate targets the original 15 donor images plus six selected
additions from the separate QuPath acquisition sources. The 6090/6539 source
mix-up has been reconciled against the raw channel folders:

| Addition | Inclusion rule |
|---|---|
| 6090 | Existing selected stack verified against all 36 original raw channel planes |
| 6232 | Include selected new acquisition |
| 6539 | Use the reconstructed stack from its distinct raw acquisition; preserve the old duplicate and its annotation separately |
| 6558 | Include `PAN-6558-AAb_test35_35p-LED.tif`; exclude the alternate `CX_PAN-6558` image |
| 6520 | Include 6520-01A; exclude 6520-01B |
| HPAP-078 | Include and retain HPAP provenance; normalize spreadsheet identifier `78.0` explicitly to this donor |

The corrected six viewer TIFFs total 46,877,746,573 bytes. All original-resolution
pixels and all five reduced pyramid levels were validated against each selected
source. Source and pixel signatures are distinct across the corrected inventory.
The preparation selected 20 nPOD images plus one HPAP image; the subsequent
local runtime is documented above. Stage images from the explicit corrected
inventory, not from the whole development image directory.

The restricted, Git-ignored candidate record at
`data/development/expanded-20260911/reports/image_conversion/` contains
`conversion_summary.json`, per-image reports, and the initial
`cohort_identity_audit.json`. Those original reports preserve the rejected source
selection. The corrected `reconciled_inventory.json` supersedes that selection
and checks current source/TIFF/sidecar hashes. The verified 6090 source hash is
`81b5a874aebc828b5c6c5925d79736b17497bf286e0aefcdf9cebedc9275c2f8`;
the reconstructed 6539 source hash is
`1866d5e8856dc096c1829e04b3622bc796cd9b81be4225bff12b864d4df1c69e`.
The 6539 processing source is
`data/development/expanded-20260911/recovery/sources/PAN-6539-ND test 34.tif`;
its corrected viewer is in `recovery/images/`. Its raw planes, acquisition
calibration, channel mapping and pixel checks are recorded in
`recovery/source_reconstruction_6539.json`. The old exclusion annotation was
drawn on duplicated 6090 pixels and was retained separately, not transferred
to a different acquisition. Corrected absolute paths and current counts are in
[the expansion status record](../docs/cohort_expansion_status.md).

The `KRT8` workflow alias for 6520/HPAP-078 was checked against the existing
classifier convention for that acquisition panel. Its evidence is in
`reports/channel_mapping/krt8_workflow_alias_evidence.json`. Source and viewer
TIFF channel labels remain unchanged.

Use the established tissue/islet classification, boundary refinement,
InstanSeg cell segmentation, region assignment, measurement/export, and
pyramidal OME-TIFF/sidecar workflow. The approved InstanSeg model is **0.1.1**;
retain that version change and full-resolution segmentation review in the
candidate record. These are offline preparation steps and are not installed as
PubApps runtime services.

Cell-phenotyping provenance remains unresolved. The historical
`7571b6d:scripts/phenotype_evaluator.py` note states that the original rules
producing shipped ECAD-positive/Structural labels were unavailable and that v2
was a reconstruction requiring biological review. The current
`data/phenotype_rules.csv` cannot generate all 21 labels in the shipped analysis.
This historical gap does not require reconstruction of retired labels. Audit
retained labels and validate their upstream gates before assigning certified
new-cell phenotypes. Existing segmentation, marker measurements, or viewer
conversion do not demonstrate equivalent phenotype assignment. Two donor-disjoint
transfer benchmarks failed their predeclared integration screen: 86.67% and
87.99% accuracy. Their feature transforms, sampling, per-class results and
limitations are retained; no new-cell phenotypes were assigned.

The earlier scVI/Age/Gender/Harmony trajectory proposal is superseded by the
fixed outcome-blind compartment-level method described above. New coordinates
have arbitrary orientation and are not numerically comparable with the former
INS-oriented axis. Preserve the frozen `5809be4` source, immutable release,
original measurements and reference results independently.

Preparation counts, raw measurement roundtrips, regional joins and calibration
are recorded in the candidate audits linked from the expansion status.
All six additions completed offline preparation: 943,068 cell objects, 942,867
fully measured objects, 201 unmeasured objects and 2,989 islet regions.
Measurement TSVs now retain Java double precision; the five earlier rounded
exports/imports are archived in `precision_revision_20260916/`. The raw H5AD
stores marker means as float32 and centroids as float64. Independent spatial
checks use unrounded centroids and QuPath's exact polygon boundary convention,
retaining overlapping memberships. Raw imports preserve wholly unmeasured
objects as flagged NaN rows for review;
the current candidate excludes them from marker analyses without dropping raw
objects. Retained-gate validation and immutable production source-to-release
acceptance remain pending; local hashes and capacity are recorded above.
All counts and validation outcomes in the
historical lineage below retain their original cohort and date; they are not
results for the expanded candidate. Keep HPAP-specific source and applicable
usage/disclosure terms distinct from the nPOD metadata lineage.

For explicit-cohort clinical metadata, `scripts/build_donor_metadata.py` now
accepts `--cohort scripts/cohort_expansion/cohort.csv` and `--donor-key` pointing
to the source CSV/XLSX with `Case ID`, `Donor Status`, `Age`, and `Gender`.
Supply a separate candidate `--output`; the original canonical CSV is protected
from this mode. Optional `--source` supplements only rows whose `donor_source`
is `nPOD`. HPAP-078 remains an HPAP donor with no nPOD lookup; absent extended
clinical values and unreported autoantibody results remain null. Donor-key and
nPOD source hashes, source cohort, normalized source identity, and
`inclusion_status` accompany the derived rows. Identity-held rows are retained
as held metadata in the original selection record. The separately verified
6090/6539 project uses `scripts/cohort_expansion/cohort_recovered.csv`; original
cohort/checkpoint fingerprints are preserved. The reconciled inventory joins
those two recovered cases with the four original unambiguous additions.
The original command without `--cohort` retains its 15-donor defaults and schema.

Current source-preparation note (2026-09-08): the 5,214 H5AD observation count
below is distinct from the embedded QuPath/prepared-table universe. Canonical
preparation currently returns `core_area` 5,382 rows, `comp` 5,382 rows,
`targets_all` 48,438 rows, `markers_all` 64,584 rows, and `combined_pairs` 21
definitions. Historical pipeline/check counts below describe their named
artifacts. They do not establish parity for a new reviewer release; see
[the current evidence index](../docs/ufrc_readiness.md).

## Canonical H5AD Lineage

```
CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, CANONICAL single-cell)
  Location: single_cell_analysis/
  Contents: scVI single-cell CODEX data, 31 protein markers. NOTE: trained
            with categorical_covariate_keys=["Age","Gender"] but no batch_key,
            so X_scVI carries donor-specific technical variance. Downstream
            Harmony correction in reaggregate_islets.py (Steps 5b/8) handles
            this without retraining scVI.
  Key fields: .obsm['X_scVI'], .obs['phenotype'], .obs['imageid'], .obs['Donor Status']
  │
  ↓ scripts/reaggregate_islets.py  (executed 2026-02-23)
  │   (region='separate', min_cells=0, require_paired=True → aggregate + trajectory + Leiden)
  │
islets_core_fixed.h5ad  (5,214 islets, core only)
  Location: islet_analysis/
  Contents: Aggregated islet-level data, proteins in .X, phenotype proportions in .obsm
  Key fields: .obsm['X_scVI_mean'], .obsm['phenotype_proportions'], .obs['donor_status']
  │
  ↓ (trajectory computed within reaggregate_islets.py)
  │   (neighbors n=30 cosine on X_scVI_harmony_pt, PAGA → UMAP, DPT with INS root)
  │
data/adata_ins_root.h5ad  (5,214 islets + pseudotime + UMAP)
  Location: data/
  Contents: Trajectory-annotated islet data
  Key fields: .obs['dpt_pseudotime'], .obsm['X_umap'], .obsm['X_diffmap']
  Validation: INS r=-0.570, GCG r=0.143, donor ordering PASS, 6533 PASS (191 islets)
  NOTE: measured from the file 2026-08-14. Earlier revisions of this block said
    "n=15" and "INS r=-0.601, GCG r=0.260" — both stale (they described a raw-scVI
    variant that is not what shipped). See scripts/CLAUDE.md § Step 5.
  │
  ↓ scripts/build_h5ad_for_app.py  (executed 2026-02-23)
  │   (+ groovy TSV data + donor metadata + neighborhood + Leiden)
  │
data/islet_explorer.h5ad  (5,214 observation rows plus larger embedded QuPath tables, ~70 MB)
  Location: data/
  Contents: Complete app data — trajectory + groovy + donor + neighborhood + Leiden clustering
  Groovy data stored in .uns as column arrays (groovy_targets_*, groovy_markers_*, etc.)
  Neighborhood metrics stored in .obs (peri_prop_*, immune_*, enrich_z_*, min_dist_*)
  Leiden clustering stored in .obs (leiden_0.3, leiden_0.5, leiden_0.8, leiden_1.0, leiden_umap_1/2)
  Donor metadata derived from h5ad obs (not external donor key)
  │
  ↓ Shiny app loads via app/shiny_app/R/data_loading.R
  │   (load_master_auto: H5AD → Excel fallback)
```

### Neighborhood Metrics Branch (Phase 7, 2026-02-19)

```
CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells)
  + islets_core_fixed.h5ad  (qualified islet list)
  │
  ↓ scripts/compute_neighborhood_metrics.py
  │   (Parent column parsing, peri-islet composition, immune metrics, enrichment z-scores, distances, per-cell signed-distance KDE table)
  │
data/neighborhood_metrics.csv  (5,214 rows × ~100 columns)
data/signed_distance_kde.csv   (~349K per-cell rows: signed_dist + case_id — Card C donor-averaged KDE)
data/cell_distance.csv         (~2.65M ALL-cell rows: signed_dist + phenotype — Card C per-cell density/composition)
  │   All 5,214 islets have peri-islet data (100% coverage with require_paired=True)
  │   Apr 2026: enrich_z_* and min_dist_* extended from 7-immune-subset to all
  │   20 phenotypes (drops "Unknown"). Aggregate immune ratios retained for
  │   compatibility but no longer surfaced in any UI selector.
  │
  ↓ Merged into islet_explorer.h5ad .obs by build_h5ad_for_app.py (step 4.5)
```

### Single-Cell Drill-Down Branch (Phase 8, 2026-02-19)

```
CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells)
  + islets_core_fixed.h5ad  (qualified islet list)
  │
  ↓ scripts/extract_per_islet_cells.py
  │   (core + peri cells, 37 columns: coords, phenotype, region, morphology, 31 markers)
  │
data/cells/*.csv  (5,214 files, ~203 MB total, ~66 cells/file avg)
  │
  ↓ Loaded by app/shiny_app/R/drilldown_helpers.R → load_islet_cells()
```

### Per-Donor Tissue CSV Branch (Phase 9, 2026-02-20)

```
CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells)
  │
  ↓ scripts/extract_per_donor_tissue.py
  │   (ALL cells: core + peri + tissue, 5 columns: X/Y coords, phenotype, cell_region, islet_name)
  │
data/donors/*.csv  (15 files, ~78 MB total, ~177K cells/donor avg)
  │
  ↓ Loaded by app/shiny_app/R/spatial_helpers.R → load_donor_tissue()
```

### Leiden Clustering Branch (Phase 9, 2026-02-20)

```
islets_core_fixed.h5ad  (5,214 islets)
  │
  ↓ Leiden clustering (4 resolutions: 0.3, 0.5, 0.8, 1.0) + UMAP
  │   (computed within scripts/reaggregate_islets.py)
  │
islet_analysis/islets_core_clustered.h5ad
  │   Columns: leiden_0.3, leiden_0.5, leiden_0.8, leiden_1.0, X_umap (2D)
  │
  ↓ Merged into islet_explorer.h5ad .obs by build_h5ad_for_app.py (step 4.7)
```

## scVI Batch Correction — Covariate Selection (2026-02-23)

### Model Configuration
The canonical single-cell H5AD uses scVI (v1.4.0) with the following configuration:
- **Architecture**: n_latent=10, n_layers=2, n_hidden=128, dropout=0.1
- **Likelihood**: negative binomial (`gene_likelihood='nb'`, `dispersion='gene-batch'`)
- **Covariates**: `categorical_covariate_keys=["Age", "Gender"]` — NO `batch_key`
- **Training**: max_epochs=400, early_stopping patience=45, NVIDIA RTX 4000 Ada GPU

### Why Age/Gender Covariates?
Each of the 15 nPOD donors has a unique Age value (float64 precision, e.g., 20.59 vs 20.70 vs 20.75). Combined with Gender, Age+Gender forms a perfect bijection to donor identity. This provides **partial donor correction** via the categorical covariate mechanism (one-hot encoded to encoder/decoder) without aggressively removing disease-associated biological signal.

> **Do not read this as "scVI is batch-corrected".** The covariates absorb *some* donor
> variance — the retrain experiment below proves it — but scVI's batch decoder is activated
> only by `batch_key`, which was never passed. Enough donor structure survived to make donor
> 6356 a standalone Leiden cluster (80–92% dominance, Cramér's V ≈ 0.4), which is exactly why
> Phase 17 added downstream Harmony. The original wording here also asserted that
> `batch_key=imageid` "would" remove disease signal because imageid is 1:1 with
> donor_status — that was never tested, only asserted, and it should not be cited as a
> finding. See `scripts/CLAUDE.md` § "don't claim scVI is batch-corrected" for the full
> reconciliation.

### Validation: Covariate vs No-Covariate Comparison
A full retraining experiment (scripts/retrain_scvi_no_covariates.py) compared the original model (with Age/Gender) against an identically-configured model without any covariates:

| Metric | Original (Age/Gender) | No Covariates | Better Model |
|--------|----------------------|---------------|--------------|
| R² imageid (donor variance) | 0.397 | 0.554 | Original (lower = less donor signal) |
| R² donor_status (biology) | 0.107 | 0.149 | No covariates (higher = more signal) |
| Excess donor eta² | 0.291 | 0.405 | Original (lower = less confounding) |
| Donor silhouette (cosine) | -0.041 | +0.250 | Original (near 0 = well-mixed) |
| k-NN same-donor enrichment | 8.4x | 11.9x | Original (lower = better mixing) |
| PC0 pseudotime gradient | Smooth (r=0.71) | Buried (r=0.10) | Original (biology drives PC0) |

**Conclusion**: The original model wins 4/5 key metrics. Removing covariates causes donor identity to dominate the latent space, burying the disease-progression signal. The small increase in biological R² (0.107→0.149) is offset by a much larger increase in donor confounding.

### Supporting Files
- `scripts/retrain_scvi_no_covariates.py` — Retraining script (identical architecture, no covariates)
- `scripts/plot_scvi_comparison.py` — 8-panel supplemental figure generator
- `scripts/scvi_comparison_results/comparison_results.json` — All numeric results
- `scripts/scvi_comparison_results/supplemental_scvi_comparison.png/.pdf` — Publication figure
- `scripts/scvi_comparison_results/FigureGuide.txt` — Panel descriptions and suggested caption

### T Cell Subtype Fix & Pipeline Rebuild (Phase 17, 2026-03-24)

```
CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells)
  │
  ↓ scripts/fix_tcell_subtypes.py
  │   (Reclassifies T cell lineage: CD8a/CD4 gate=0.5, only 2 cells changed — both tissue background)
  │   (Backup: CODEX_scvi_BioCov_phenotyped_newDuctal.backup_20260324_124235.h5ad)
  │
  ↓ Full downstream pipeline rebuild:
  │   reaggregate_islets.py → compute_neighborhood_metrics.py
  │   extract_per_islet_cells.py → extract_per_donor_tissue.py
  │   build_h5ad_for_app.py
  │
data/islet_explorer.h5ad  (rebuilt, 60.1 MB)
```

Also fixed `phenotype_rules.csv` header marker names:
- `CD8` → `CD8a`, `HLA-DR` → `HLADR`, `Podoplanin` → `PDPN`
- These now match `adata.var_names` in the canonical H5AD

Trajectory validation: INS r=-0.628, GCG r=0.256, donor ordering PASS.
Leiden clustering: 5/8/12/15 clusters at resolutions 0.3/0.5/0.8/1.0.

## Data Sources

### Single-cell (canonical)
- **File**: `single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad`
- **Provenance**: scVI batch correction with Age/Gender categorical covariates (see covariate validation above), phenotyped with Rules1
- **Donors**: 15 nPOD donors (ND, Aab+, T1D)
- **Markers**: 31 protein markers

### QuPath exports (groovy)
- **Location**: `panc_CODEX/results/groovy_exports/`
- **Types**: islet_targets, islet_markers, islet_composition, LGALS3
- **Per donor**: 15 donors × 4 types = 60 TSV files
- **Master Excel**: `data/master_results.xlsx` (built by `scripts/build_master_excel.py`)

### Donor metadata
- **Source**: Derived from h5ad obs (imageid → donor_status, age, gender, AAb flags)
- **Note**: `panc_CODEX/results/CODEX_Pancreas_Donors.xlsx` has a different donor cohort (not used)

### Segmentation boundaries
- **GeoJSON**: `data/json/*.geojson` (uncompressed, per-donor)
- **Compressed**: `data/gson/*.geojson.gz` (gzipped, per-donor)
- **Lookup**: `data/islet_spatial_lookup.csv` (centroid coordinates)

### Neighborhood metrics (Phase 7, updated Phase 11 / Apr 2026 systematic)
- **File**: `data/neighborhood_metrics.csv` (5,214 rows, ~100 columns)
- **Provenance**: Computed from single-cell H5AD by `scripts/compute_neighborhood_metrics.py`
- **Contents**: Peri-islet composition (21 phenotypes), immune infiltration (5 metrics), enrichment z-scores (20 types), distance metrics (per-phenotype `min_dist_*` + `min_dist_immune_mean`)
- **Coverage**: All 5,214 islets have peri-islet data (100% with require_paired=True)

### Signed-distance KDE (Spatial Card C, Jul 2026)
- **File**: `data/signed_distance_kde.csv` (~349K rows: one per islet-associated cell)
- **Provenance**: Emitted by `scripts/compute_neighborhood_metrics.py` alongside `neighborhood_metrics.csv` (regenerated in lockstep, so it can't drift from the current phenotyping)
- **Contents**: `phenotype` (safe-name), `case_id`, `donor_status`, `signed_dist` (QuPath `Dist to Closest Islet`; neg = inside core, pos = peri), `cell_region`, `is_immune`
- **Consumed by**: Spatial Card C KDE, which averages per-donor densities (equal donor weight) instead of pooling all cells — removing the pseudoreplication a single large donor (ND 6356) caused

### All-cells cell-distance (Spatial Card C per-cell views, Jul 2026)
- **File**: `data/cell_distance.csv` (~2.65M rows: one per cell — ALL cells, not just islet-associated)
- **Provenance**: Emitted by `scripts/compute_neighborhood_metrics.py` alongside the other two outputs
- **Contents**: `case_id`, `donor_status`, `cell_region` (core/peri/tissue), `phenotype` (raw name), `signed_dist` (QuPath `Dist to Closest Islet`; −182 → +1131 µm), and `{marker}_pos` int columns for ~16 lineage markers (scimap-rescaled X>0.5 — the phenotyping gate — for the app's "Marker+" mode)
- **Consumed by**: Card C-right "Per-cell distance" view — symlog density + phenotype composition. Ported from Islet-Explorer-Senior; CODEX phenotypes only (no RESTORE marker positivity, no DuckDB)

### Per-islet cell data (Phase 8, updated Phase 11)
- **Directory**: `data/cells/` (5,214 CSV files, ~203 MB)
- **Provenance**: Extracted from single-cell H5AD by `scripts/extract_per_islet_cells.py`
- **Contents**: 37 cols per cell: X/Y centroids (µm), phenotype, core/peri region, Cell/Nucleus Area, 31 markers
- **Naming**: `{imageid}_Islet_{N}.csv` matches `combined_islet_id`

### Per-donor tissue data (Phase 9)
- **Directory**: `data/donors/` (15 CSV files, ~78 MB)
- **Provenance**: Extracted from single-cell H5AD by `scripts/extract_per_donor_tissue.py`
- **Contents**: 5 cols: X_centroid, Y_centroid, phenotype, cell_region (core/peri/tissue), islet_name
- **Naming**: `{imageid}.csv` (e.g., `6505.csv`)
- **Coverage**: All 15 donors with core/peri/tissue cells.

### Leiden clustering (Phase 9)
- **Source**: `islet_analysis/islets_core_clustered.h5ad`
- **Provenance**: Leiden clustering on scVI embeddings at 4 resolutions
- **Contents**: leiden_0.3, leiden_0.5, leiden_0.8, leiden_1.0 (categorical cluster assignments) + X_umap (2D UMAP coordinates)
- **Merged into**: `islet_explorer.h5ad` .obs via `build_h5ad_for_app.py` step 4.7

## Deprecated Files

| File | Status | Replacement |
|------|--------|-------------|
| `data/adata_ins_root.h5ad` (old) | Replaced | Rebuilt from fixed pipeline |
| `islet_analysis/islets_codex_panc.h5ad` | Legacy | Use islets_core_fixed.h5ad |
| `islet_analysis/islets_codex_panc_all_islet_scvi.h5ad` | Legacy | Use islets_core_fixed.h5ad |
| `islet_analysis/islets_codex_panc_all_islet_scvi_anno.h5ad` | Legacy | Use islets_core_fixed.h5ad |

## Known Issues

### Old adata_ins_root.h5ad (pre-Phase 3, replaced 2026-02-17)
- Mixed morphology + protein features in .X (47 vars instead of 31)
- Used n_neighbors=5 (non-standard)
- Missing X_scVI_mean (no batch correction in embeddings)
- Missing phenotype_proportions
- 5,158 islets (included peri-islet, combined regions)
- Donor 6533 was missing (fixed Mar 2026 — Parent annotations patched)
- **Replaced by**: 1,015 islets, 31 protein vars, scVI embeddings (note: scVI training was not batch_key-aware; Apr 2026 added downstream Harmony)

### Duplicate files
- `islet_analysis/` files exist in both `panc_CODEX/` and `Islet-Explorer/`
- `annotations.tsv` exists in both `data/` and `islet_analysis/`
- Cleanup deferred to after pipeline validation

## Validation Criteria

### scVI QC (`notebooks/scvi_qc_validation.ipynb`) — executed 2026-02-17
- [PASS] Silhouette batch: PCA=-0.021, scVI=-0.068 (scVI better)
- [PASS] Silhouette cell-type: PCA=-0.077, scVI=-0.013 (preserved)
- [PASS] Batch LISI: PCA=2.73, scVI=8.07 (ideal=15)
- [PASS] Donor 6533 present (205K cells)

### Trajectory QC (`notebooks/rebuild_trajectory.ipynb`) — executed 2026-02-17 (949 islets, superseded)
- [PASS] INS vs pseudotime: r=-0.741 (threshold: r < -0.3)
- [PASS] GCG vs pseudotime: r=0.379 (threshold: r > 0.2)
- [PASS] Donor ordering: ND=0.491, Aab+=0.516, T1D=0.743
- [PASS] No donor dominance: max 38% in any quintile (threshold: <80%)
- [PASS] Donor 6533 included

### App QC — validated 2026-02-17
- [PASS] H5AD loads via load_master_auto() (17s)
- [PASS] prep_data() produces identical output from H5AD and Excel
- [PASS] ANOVA computes without errors
- [PASS] Trajectory H5AD: 1,015x31, UMAP + pseudotime + combined_islet_id
- [PASS] Excel fallback works when H5AD is removed

### Neighborhood QC (`scripts/compute_neighborhood_metrics.py`) — executed 2026-02-19 (949 islets, superseded)
- [PASS] Output: 1,015 rows × 62 columns
- [PASS] 949 islets (93.5%) have peri-islet data (total_cells_peri > 0)
- [PASS] Biological signal: T1D immune_frac_peri (0.155) > Aab+ (0.106) > ND (0.069)
- [PASS] H5AD rebuild: 60 new .obs columns added, islet_explorer.h5ad = 47.6 MB
- [PASS] All R files parse cleanly, app loads without errors

### Single-Cell Drill-Down QC (`scripts/extract_per_islet_cells.py`) — executed 2026-02-19 (949 islets, superseded)
- [PASS] Output: 949 CSV files, 180,019 total cells, 110.8 MB
- [PASS] Per-cell columns: 37 (X/Y, phenotype, region, morphology, 31 markers)
- [PASS] File naming matches combined_islet_id (e.g., `6505_Islet_284.csv`)
- [PASS] drilldown_helpers.R loads and caches CSV files correctly

### Phase 11: Full Pipeline Rebuild (min_cells=0) — executed 2026-02-23

#### Reaggregation (`scripts/reaggregate_islets.py`)
- [PASS] Aggregation: 5,214 paired islets (from 5,126 core, 5,054 peri, 5,214 paired)
- [PASS] Trajectory: INS r=-0.591 (threshold: r < -0.3)
- [PASS] GCG correlation: r=0.308 (threshold: r > 0.2)
- [PASS] Donor ordering: ND=0.488, Aab+=0.506, T1D=0.571
- [PASS] Donor 6533: 191 islets (Parent annotations patched in canonical H5AD, Mar 2026)
- [PASS] Leiden clustering: 6/8/14/14 clusters at resolutions 0.3/0.5/0.8/1.0

#### Neighborhood metrics
- [PASS] Output: 5,214 rows × 62 columns
- [PASS] 100% peri-islet coverage (require_paired=True eliminates NaN)

#### Per-islet cell CSVs
- [PASS] Output: 5,214 files, 329,413 total cells, 203 MB
- [PASS] Avg 66 cells/islet (down from 190 due to inclusion of small islets)

#### App H5AD
- [PASS] islet_explorer.h5ad: 5,214 islets, 70 MB
- [PASS] All neighborhood, Leiden, and groovy data merged successfully

### Phase 17: Donor Batch Correction & Systematic Phenotype UI — executed 2026-04-25

Root-cause: `single_cell_analysis/SCVI_CODEX_v2.ipynb` cell 16 trained scVI without `batch_key="imageid"` (only Age/Gender as covariates). Donor-specific technical variance (especially donor 6356, ND, 18.8% of cells) leaked into `X_scVI_mean`, manifesting as donor-dominated Leiden clusters (Cramer's V ≈ 0.4, AMI(donor) > AMI(status)) and a noisy pseudotime within Aab+ status.

#### Reaggregation v2 (`scripts/reaggregate_islets.py` Apr 2026 rebuild)
- [PASS] Step 5 (core pseudotime): medoid-of-top-50-ND root (idx=703, donor=6516, INS=0.663). **INS r=−0.570, GCG r=0.143** on the shipped file, built on `X_scVI_harmony_pt` (Harmony θ=0.1), n_neighbors=30. The figures previously recorded here — "raw scVI", INS −0.694, GCG 0.274 — described a raw-scVI variant that did not ship; re-measured 2026-08-14.
- [PASS] Step 5b (combined pseudotime): 36-d input [core scVI(10) | peri scVI(10) | α·z(struct(16))], α=0.15, Harmony θ=0.1 on imageid. INS r=−0.51, CD8a r=+0.225 (vs +0.190 core), CD3e r=+0.178 (vs +0.125 core), CD31 r=+0.131, CD34 r=+0.204, PDPN r=+0.073.
- [PASS] Step 6 validation: ND<Aab+<T1D ordering, eta²(donor|pt) per status — ND 0.033, Aab+ 0.090, T1D 0.024.
- [PASS] Step 8 (Leiden): Harmony θ=2 on imageid. Cramer's V donor↔cluster (leiden_0.5) 0.38 → 0.135; donor 6356 worst-cluster fraction 0.83 → 0.29; AMI(donor) 0.19 → 0.028; 0 clusters with ≥60% one donor.

#### Systematic neighborhood metrics (`scripts/compute_neighborhood_metrics.py`)
- [PASS] `enrich_z_*` and `min_dist_*` now computed for all 20 phenotypes (drops "Unknown") — was 7-immune subset.
- [PASS] Output: 5,214 rows × ~100 columns. `IMMUNE_TYPES` constant retained for aggregate metrics only (immune_frac_*, immune_ratio, etc.).
- [PASS] Spatial tab Cards A/B/C selectors populate dynamically from `phenotype_options()` reactive — no more hand-curated immune subset.
- [PASS] Plot tab Composition selector: "Immune Metrics" option group removed.

## Donor clinical metadata (Jul 2026)

`data/donor_metadata.csv` — 15 rows, built by `scripts/build_donor_metadata.py`.

**Source:** an nPOD "All-Download" export (`All-Download_From_nPOD_YYYY-MM-DD_HHMMSS.xlsx`, 741
donors x 47 columns) kept in the repo root and **deliberately gitignored** (`*.xlsx`). It carries
free-text admission courses, clinical histories, medication lists and HLA types for every nPOD
donor; only 15 donors and ~12 numeric/categorical fields are relevant here, so only the derived CSV
is committed. Re-run the script against a fresh export to refresh.

**Case-ID mapping:** donor `112` is filed upstream as **`HDL112`** (89 nPOD donors carry that
prefix). All other 14 match on the bare numeric ID. The script fails loudly if any donor is unresolved.

**Fields added over the H5AD's age + sex:** `t1d_duration_years` (nPOD `diabetes_hx_years`),
`BMI`, `HbA1c_percent`, `C_peptide_ng_mL` (+ `C_peptide_censored` — one donor is reported as
`<0.02`), the four AAb results as logicals, `n_aab`, and `age_group` (<13 / >=13).

**Two caveats recorded at build time:**

1. **Disease duration is a descriptor, not a regressor.** All five T1D donors are recent-onset
   (0, 0, 0, 0.58, 0 years) and non-diabetic donors have no value. There is no variance to model.
   It establishes cohort scope — no long-duration T1D donors — and nothing more.
2. **One age discrepancy.** Donor 112 is 3.00 years in the H5AD and 3.72 in nPOD; the other 14
   match to the decimal. The H5AD `age` remains canonical for all modelling so existing results do
   not shift; `age_years_npod` is carried alongside and the Inflammation tab's donor table surfaces
   the mismatch. Both values fall under the <13 age-group cut, so binning is unaffected.

**Cross-check:** the H5AD and nPOD autoantibody flags agree for all 15 donors x 4 assays.
