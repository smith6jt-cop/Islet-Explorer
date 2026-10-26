# CLAUDE.md - Pipeline Scripts

> **2026-10-15 owner directives:** donor 78 (HDL078 / HPAP-078) is removed from all analysis. `excluded_donors.csv` is the record and `excluded_donors.py` the guard: cohort loaders (`research/common.py`, `cohort_expansion/prepare_analysis.py`, `build_donor_metadata.py`) and the protein pipeline reject it, so the historical 21-donor `research/cohort-v1.csv` no longer loads. The current protein cohort is `protein_expansion_cohort.csv` (20 donors). All 20 donors are to be reprocessed through one scripted workflow; see [the reprocessing plan](../docs/uniform_reprocessing_plan_20261015.md).

> **2026-10-03 owner decision:** the accepted reconstruction rules and all six donor boundaries are approved. Use Harmony on protein PCA, retain the uncorrected comparison, and preserve original source records. Exact historical typing/encoder recovery and scVI retraining are no longer prerequisites. See the [current implementation record](../docs/protein_expansion_20261003.md); earlier directions below are historical.

Pipeline scripts that build the data files consumed by the Shiny app.

## Restoration checkpoint — 2026-09-23

The original 15-donor builders and trajectory entry points have been restored
from `000e970`. Their restored defaults are for the frozen original lineage;
do not run them in place or apply them to candidate inputs. The research
pipeline below remains isolated and untrusted pending revision and evaluation.
The owner requires preservation of original IDs, physical geometry,
measurements, missingness and phenotype labels. A validated common 21-donor
trajectory, UMAP and Leiden refit may change derived outputs. See the current
phase-2 contract before any app integration.

## Historical owner correction — 2026-09-19, revised 2026-09-29

[Typing must precede endocrine-normalized trajectory features](../docs/trajectory_typing_reference_20260919.md).
The nine-marker enhancement recorded in
`research/cell_typing_requirements-v1.json` is deferred for the current expansion.
As of 2026-09-29 the original numeric gates were recovered, but the complete
phenotype rules and scVI encoder were unavailable, and the plan was to resolve
model projection versus a reviewed joint retrain. The 2026-10-03 owner decision
superseded that plan: the six additions use the accepted reconstruction rules
and Harmony on protein PCA, and `retrain_joint_scvi.py` is archived. Historical
schema 4 ran cell typing before composition and fitted trajectory from typed
counts; its outputs are not accepted for the restored app.

## Historical offline defaults — 2026-09-22 (superseded)

Current additions follow [phase 2](../docs/donor_expansion_phase2.md).
The historical 21-donor experiment is documented in [research/README.md](research/README.md).
The schema-4 path is `research/reprocess_candidate.py` followed by correction,
typing, composition, trajectory, `build_islet_embedding.py`, runtime
materialization, validation and finalization. Every builder requires explicit
input/output roots. Research runners require `--report-root` outside app data.
The historical local candidate was
`data/development/all-donor-20260922-trajectory-rank-v2`; preserve the native QuPath
project, source images, prior candidates and frozen release.
The fixed 0–1 reference rank is a display/binning artifact: raw pseudotime
orders observations, and exact ties are spread by SHA-256 of islet ID. Outcome
trajectory inference uses `pseudotime_raw`, where the ties remain intact.
The former ND/INS-root, disease-order validation and small/large refits are
archived under `archive/trajectory-pre-20260916/`. Preserve those historical
sources and frozen outputs; do not run their defaults against the new candidate.
Retained gates replace reconstruction of retired labels; unresolved outputs
remain unavailable. See [the delivery record](../docs/all_donor_delivery.md).

## Historical original-cohort methods

The following descriptions and commands document the frozen original pipeline;
their old parameters and entry-point arguments are not current defaults.

## Donor Batch Correction & Pseudotime Modes (Phase 17, Apr 2026)

### scVI batch_key bug (root cause)

`single_cell_analysis/SCVI_CODEX_v2.ipynb` cell 16 set `categorical_covariate_keys=["Age", "Gender"]` but **never passed `batch_key="imageid"`**, despite comments claiming otherwise. Result: scVI's batch decoder was never activated and donor-specific technical variance (staining, antibody lot, imaging gain) leaked into `X_scVI_mean`. Symptoms: donor 6356 (ND, 18.8% of cells) formed a standalone Leiden cluster at every resolution (80–92% dominance); Cramer's V donor↔cluster ≈ 0.4 with AMI(donor) > AMI(status).

**Fix:** downstream Harmony correction (avoids scVI retrain on 2.6M cells). Two different θ values for two different goals.

### Step 8: Harmony for Leiden clustering (θ=2)

`scripts/reaggregate_islets.py` Step 8 applies `harmonypy.run_harmony` on `X_scVI_mean` with `batch=imageid`, `theta=2.0` before `sc.pp.neighbors`. Embedding stored as `obsm['X_scVI_harmony']`. Goal: donor-invariant islet cell-types. Post-fix metrics: Cramer's V dropped from 0.38 → 0.135, donor 6356's worst-cluster fraction 0.83 → 0.29, AMI(donor) 0.19 → 0.028. Status-orthogonality is fine/desired (clusters are cell-types, disease lives on pseudotime).

### Step 5: Trajectory pseudotime (Harmony θ=0.1 core, robust medoid root)

**The shipped `data/adata_ins_root.h5ad` was built with light Harmony (θ=0.1), not raw scVI.**
Measured from the file itself: `uns['neighbors']['params']` = `{use_rep: 'X_scVI_harmony_pt',
n_neighbors: 30, metric: 'cosine'}`. An earlier revision of this section claimed "raw
`X_scVI_mean` (no Harmony)" and `n_neighbors=15`; both were stale — never updated after θ=0.1
was adopted — and the numbers quoted below the old text belonged to the raw variant. The θ
sweep that drove the choice is recorded as a comment at `reaggregate_islets.py:124-134`; the
shipped values match its θ=0.1 row exactly.

Aab+ donor variance accepted as biology (Aab+ donors exhibit a heterogeneous range of islet
types).

**Verify rather than trust this paragraph** — read the file:

```bash
python -c "import anndata as ad; a=ad.read_h5ad('data/adata_ins_root.h5ad'); print(a.uns['neighbors']['params'])"
```

**Robust root selection** (replaces single-islet `argmax(INS)` over ND): root = medoid of top-50 ND islets by INS in scVI space. Anchored by the center-of-mass of "healthy ND, INS-high" islets and naturally draws from multiple donors. Current root: index 703, donor 6516, INS 0.663 (99.7th percentile).

### Step 5b: Combined-mode pseudotime (Harmony θ=0.1 + structures α=0.15)

Computes a SECOND pseudotime saved as `dpt_pseudotime_combined`. Input matrix: `[core_scVI(10) | peri_scVI(10) | α_struct · z(struct(16))]` = 36 dims. Structural block: (`prop`, `count`) × (Neural, Blood Vessel, Endothelial, Lymphatic) × (core, peri), z-scored, scaled by `ALPHA_STRUCT = 0.15`. Disease ordering breaks at α ≥ 0.20; α=0.15 is the largest weight preserving ND<Aab+<T1D ordering while strengthening immune/vascular marker correlations (CD3e, CD8a, CD31, CD34, PDPN). Light Harmony (θ=0.1) on the combined matrix; same medoid root as core mode.

### Step 6 validation thresholds

- INS Spearman r < −0.3
- GCG Spearman r > 0.1 (relaxed from 0.2 — see the caveat below)
- ND<Aab+<T1D mean ordering
- Donor eta² within status: ≤ 0.15 for ND/T1D, **≤ 0.25 for Aab+** (relaxed because Aab+ donors genuinely span heterogeneous islet types)
- 6533 included

**Measured on the shipped file** (not the stale figures this section used to carry — it
claimed INS −0.69 raw / −0.55 combined and GCG ~0.27):

| | core (`dpt_pseudotime`) | combined (`dpt_pseudotime_combined`) |
|---|---|---|
| INS Spearman | **−0.570** | −0.511 |
| GCG Spearman | **0.143** | **0.103** |
| SST Spearman | −0.206 | −0.225 |
| ND / Aab+ / T1D mean | 0.3040 / 0.3634 / 0.3883 | 0.2245 / 0.2619 / 0.2867 |

Core-vs-combined agreement is only Spearman ρ = **0.846** — the two modes are not
interchangeable.

⚠️ **The GCG gate is the weak one.** At r = 0.143 (core) and 0.103 (combined) it clears
only the *relaxed* 0.1 threshold and would **FAIL the original 0.2**. The relaxation is
justified solely by a comment at `reaggregate_islets.py:381-383` ("genuine shared effect is
r ≈ 0.14"); nothing independent backs that number. Combined mode sits 0.003 above the
threshold — and Step 6 never actually tests combined mode (it validates `dpt_pseudotime`
only), so that margin is unmonitored.

⚠️ **Step 6 is advisory and islet-level.** On failure it prints a warning and saves anyway,
and it aggregates islets rather than donors — so the gate is itself pseudoreplicated. See
`app/shiny_app/docs/trajectory.md`.

### Important: don't claim "scVI is batch-corrected"

Old comments in `reaggregate_islets.py` and `islet_analysis/CLAUDE.md` claimed "scVI latent has donor-level variation corrected out (Age+Gender covariates bijectively map to donors)". This is wrong — Age/Gender as `categorical_covariate_keys` don't activate batch correction; only `batch_key` does. Comments updated Apr 2026.

**Reconciling this with `data/DATA_PROVENANCE.md`.** That file reports a controlled retrain
experiment (`scripts/retrain_scvi_no_covariates.py`) showing the Age/Gender model *beats* a
no-covariate model on 4 of 5 donor-mixing metrics — donor η² 0.291 vs 0.405, PC0 pseudotime
gradient r = 0.712 vs 0.101. Read alongside this section that can look contradictory. It is
not. Both are true, and they answer different questions:

- **Do the Age/Gender covariates absorb *some* donor variance?** Yes — they are one-hot fed
  to the encoder/decoder and Age is a bijection to donor, so they demonstrably do. That is
  what the retrain experiment measured, and it is why removing them made things worse.
- **Is that the same thing as scVI batch correction?** **No.** scVI's batch decoder is
  activated only by `batch_key`, which was never passed. The residual donor structure that
  survived (donor 6356 as a standalone Leiden cluster at 80–92% dominance, Cramér's V ≈ 0.4)
  is what Phase 17's downstream Harmony exists to fix.

So: the covariates help and should stay; they are **not** batch correction; and the latent
feeding DPT is best described as *partially* donor-corrected by covariates plus explicitly
Harmony-corrected at θ=0.1 for the trajectory. Do not write "scVI is batch-corrected", and do
not write "the covariates do nothing".

## Pipeline Scripts & Notebooks

- `scripts/fix_tcell_subtypes.py` - Reclassifies CD4/CD8 T cells by marker gating (CD8a/CD4 >= 0.5); `--dry-run` flag for preview
- `scripts/reaggregate_islets.py` - **Primary pipeline**: aggregation (min_cells=0) + trajectory + Leiden in one script
- `scripts/compute_neighborhood_metrics.py` - Computes per-islet peri-islet metrics from single-cell H5AD. Emits **three** files: `data/neighborhood_metrics.csv` (per-islet); `data/signed_distance_kde.csv` (per islet-associated cell: `signed_dist` + `case_id` + `is_immune`, feeds the donor-averaged Card C KDE); `data/cell_distance.csv` (ALL ~2.65M cells: `case_id, donor_status, cell_region, phenotype, signed_dist` + `{marker}_pos` int columns for ~16 lineage markers via `build_marker_positivity()` = scimap-rescaled X>0.5, the phenotyping gate; tissue-wide, feeds the Card C per-cell density + composition views incl. their phenotype/marker toggle)
- `scripts/extract_per_islet_cells.py` - Extracts per-islet cell CSVs for drill-down viewer
- `scripts/extract_per_donor_tissue.py` - Extracts per-donor tissue CSVs (ALL cells: core+peri+tissue) for Spatial tab scatter
- `scripts/build_h5ad_for_app.py` - Builds enriched H5AD (trajectory + groovy + neighborhood + Leiden)
- `scripts/build_master_excel.py` - Builds master_results.xlsx from groovy TSV exports
- `notebooks/scvi_qc_validation.ipynb` - Validates scVI batch correction (silhouette, LISI, UMAP)
- `notebooks/rebuild_trajectory.ipynb` - Regenerates trajectory from fixed pipeline (superseded by reaggregate_islets.py)
- `islet_analysis/fixed_islet_aggregation.py` - Core aggregation module (islet-level from single-cell, default min_cells=0)

## Python Environment

Pipeline scripts use the `scvi-env` conda environment:
```bash
conda activate scvi-env  # scanpy, anndata, scvi-tools, sklearn, scib-metrics, scipy, harmonypy
```

`harmonypy==0.2.0` is a downstream batch-correction dependency (added Apr 2026). One-time install: `pip install harmonypy` inside `scvi-env`.

## Running Pipeline Scripts

```bash
conda activate scvi-env

# Fix T cell subtypes in canonical H5AD (dry-run first, then apply)
python scripts/fix_tcell_subtypes.py --dry-run
python scripts/fix_tcell_subtypes.py

# Full rebuild from single-cell H5AD (aggregate + trajectory + Leiden, ~15 min)
python scripts/reaggregate_islets.py

# Compute neighborhood metrics (reads 2.6M-cell H5AD, ~3 min)
python scripts/compute_neighborhood_metrics.py

# Extract per-islet cell CSVs (reads 2.6M-cell H5AD, ~5 min)
python scripts/extract_per_islet_cells.py

# Extract per-donor tissue CSVs for Spatial tab scatter (~2 min, independent of islet filtering)
python scripts/extract_per_donor_tissue.py

# Build enriched H5AD (merges trajectory + groovy + neighborhood + Leiden)
python scripts/build_h5ad_for_app.py
```

## compute_split_pseudotime.py (Aug 2026)

Recomputes DPT **within** each cell-count stratum (`total_cells <= 10` / `> 10`) and writes
`dpt_pseudotime_small` / `dpt_pseudotime_large` / `cell_count_stratum` back into
`data/adata_ins_root.h5ad`. Repeats `reaggregate_islets.py` Step 5 exactly, including re-deriving
the root inside each stratum. Leaves the shipped `dpt_pseudotime` bit-identical.

```bash
conda activate scvi-env
python scripts/compute_split_pseudotime.py            # --split 10 by default
```

**The two strata are independent diffusion maps and are NOT on a common scale.** See
`app/shiny_app/docs/trajectory.md` §10 for the measured results and the reasoning — read that
before quoting either column.

## build_barlow_features.py (Aug 2026)

Builds `data/barlow_islet_features.csv` — the islet composition matrix for the **single reviewer
trajectory** adapted from Barlow et al. (eLife 2024;13:RP100535): per-islet cell-type counts for core and
the 20 µm peri band, divided by the islet's endocrine-cell count, log1p-transformed.

```bash
conda activate scvi-env
python scripts/build_barlow_features.py
```

Reads `islet_analysis/islets_core_fixed.h5ad` + `islets_peri_fixed.h5ad` (joined on `islet_id` —
they are **not** in the same row order). The cell-count range is deliberately NOT baked in: the app
recomputes the diffusion map for the current lower and upper bounds. Barlow et al.'s >=10-cell
criterion motivates the default lower bound. See `app/shiny_app/docs/trajectory.md` §11.
