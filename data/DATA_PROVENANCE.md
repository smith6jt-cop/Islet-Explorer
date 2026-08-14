# Data Provenance — Islet Explorer

This document describes the canonical data lineage for the Islet Explorer app.

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
data/islet_explorer.h5ad  (5,214 islets, all app data, ~70 MB)
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
