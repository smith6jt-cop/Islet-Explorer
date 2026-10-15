# Frozen 15-donor non-trajectory processing

> **Scope (2026-10-14):** this traces the frozen 15-donor root under `data/`. Since 2026-10-03 the local app serves the protein-expanded root, whose Spatial UMAP and Leiden use Harmony on protein PCA instead of the scVI latent; see [the protein expansion record](protein_expansion_20261003.md).

Read-only trace performed on 2026-09-29, starting from the restored app's actual
files. This describes the original 15 donors under `data/`, not the isolated
21-donor experiments. It distinguishes verified numerical lineage from notebook
recipes whose complete historical run cannot be recovered.

The app contains two measurement branches and a separate clustering branch.
QuPath hormone-positive counts, scimap cell phenotypes, and scVI-derived Leiden
clusters are different outputs. The app's older “Phenotype Proportions (scVI)”
wording and “scVI-denoised expression” descriptions do not identify the actual
processing of the phenotype counts or marker measurements.

## Trace backwards from the outputs

The local startup log identifies `data/islet_explorer.h5ad` as the loaded source.
`app/shiny_app/R/data_loading.R` reconstructs QuPath tables from `uns`, extracts
phenotype and neighborhood columns from `obs`, then joins them in `prep_data()`.
The H5AD contains 5,214 observation rows, but its embedded QuPath tables cover
5,382 islet cores. Single-cell spatial displays also read the CSVs below.

| Output | Immediate source | Processing before storage |
|---|---|---|
| Plot: INS/GCG/SST-positive counts, percentages and densities | `uns/groovy_composition_*` in `data/islet_explorer.h5ad` | Direct QuPath hormone thresholds and region counts; separate from scimap phenotyping |
| Plot: nerve, capillary and lymphatic area fractions | `uns/groovy_targets_*` | QuPath structure-polygon intersections with core/peri regions, divided by region area |
| Core/peri cell-type proportions, counts and densities | H5AD `obs/prop_*`, `peri_prop_*`, `count_*`, `peri_count_*` | Count saved single-cell phenotypes after QC, grouping by donor and QuPath `Parent` |
| Drill-down phenotype colors and 31 protein-marker colors | `data/cells/{donor}_Islet_{N}.csv` | Selected canonical cells; markers copied from scimap-rescaled `X` |
| Spatial whole-section cells | `data/donors/{donor}.csv` | All QC-retained canonical cells; copy centroids and phenotype, parse `Parent` into core/peri/tissue |
| Spatial enrichment and nearest-cell distances | H5AD neighborhood columns, merged from `data/neighborhood_metrics.csv` | Phenotype counts and cell coordinates; formulas below |
| Spatial signed distances to islet edge | `data/cell_distance.csv`, `data/signed_distance_kde.csv` | Carry forward QuPath signed distance and canonical phenotype |
| Spatial UMAP coordinates | H5AD `obs/leiden_umap_1/2`, copied from `islets_core_clustered.h5ad` | PCA and UMAP of mean scimap-rescaled protein values |
| Spatial Leiden cluster identities and cluster composition | H5AD `obs/leiden_0.3/0.5/0.8/1.0` | Mean single-cell scVI latent per core, donor Harmony correction, neighbor graph, Leiden |
| Boundary overlays and image viewer | GeoJSON under `data/gson`/`data/json`; pyramidal OME-TIFFs | QuPath annotation export and image conversion; separate from the numerical cell-analysis branches |

The regional CD3e/CD163/CD68/LGALS3 positivity tables are also embedded and loaded
as `markers_all`. The current visible Plot module's selectors use the hormone,
phenotype and structure branches above; the presence of a regional marker table
in the H5AD does not mean it defines the single-cell phenotype calls.

## 1. Direct QuPath summaries

`scripts/build_h5ad_for_app.py` reads
`/data/panc_CODEX/results/groovy_exports/`, joins donor metadata, and stores the
source columns in H5AD `uns`. It does not apply scVI, CLR, or scimap to those
measurements. The current embedded source fields match all 60 source TSVs:

| Files | Files checked | Stored rows |
|---|---:|---:|
| `islet_composition_*.tsv` | 15 | 5,382 |
| `islet_targets_*.tsv` | 15 | 35,523 |
| `islet_markers_*.tsv` plus `LGALS3_*.tsv` | 30 | 64,516 |

The original QuPath scripts remain under
`/data/QuPath_Projects/QuPath_project_IsletExplorer/scripts/`:

- `3_Islet_Counter.groovy`: for each image, use whole-cell raw MFI and a hormone
  threshold of **mean + 1.5 × sample SD across that image's cells**. The script
  has `perIsletThresholds = false`. Count INS/GCG/SST-positive cells inside each
  core by centroid containment. `Ins_any`, `Glu_any`, and `Stt_any` include
  co-positive cells; the single-positive categories, `Multi_Pos`, and
  `Triple_Neg` form a mutually exclusive partition. All 45 stored donor/hormone
  thresholds reconstruct from the original rounded cell export within
  **4.28 × 10^-7 raw MFI**. Thresholds are constant within each donor's TSV.
- `4_Islet_Cells.groovy` and `Old/4_Islet_Cells_addMac2.groovy`: for CD3e,
  CD163, CD68 and LGALS3, calculate **mean + 1.5 × sample SD separately within
  each region**, then count cells with MFI at or above that threshold. All
  49,637 finite threshold rows satisfy the exported mean/SD formula. These
  thresholds have a different scope from the per-image hormone thresholds.
- `5_Islet_Objects.groovy`: intersect Nerve/Capillary/Lymphatic annotation
  polygons with each islet region; convert intersection and region areas using
  image calibration. Export area, region area, area fraction and object counts.
  The visible Microenvironment mode uses **100 × structure area / region area**.
  Core/peri/union geometry is handled by the exporter, including its documented
  band-resolution fallbacks.

At runtime `prep_data()` joins by donor and islet name, calculates equivalent
core diameter as `2 * sqrt(core_area / pi)`, and synthesizes missing combined
rows from core and band components. Plot percentages are count/appropriate cell
total × 100; densities are count/QuPath area in µm² × 10^6. Combined phenotype
percentages pool counts and denominators rather than averaging percentages.

## 2. Single-cell import and QC

The original cell export is
`/data/panc_CODEX/results/old_data_plots/all_data.tsv`.

1. It contains **2,676,466 QuPath cell objects**. The imported
   `/data/panc_CODEX/results/panc_CODEX.h5ad` retains **2,676,419**: the omitted
   **47** have all 31 analysis-marker measurements missing. There are no
   partially missing rows in those 31 columns.
2. The H5AD selects 31 `Cell: <marker>: Mean` columns. DAPI, Ki67 (exported as a
   median), epCAM and KRT8-18 are absent from this analysis matrix. The imported
   values remain raw whole-cell MFI. Every one of **82,968,989** retained marker
   measurements matches the TSV exactly, with zero absolute numerical
   difference. Object IDs, donor IDs, parents and both centroid coordinates
   also match across every retained row.
3. `single_cell_analysis/SCVI_CODEX_v2.ipynb`, cells 8–9, calculates QC on raw
   MFI before transformation: require at least **2** positive markers and total
   MFI at or above the pooled original cohort's first percentile,
   **5234.796272**. The executable value is 2 even though its adjacent comment
   incorrectly says 15. The marker-count criterion removes no additional cells;
   the total-MFI criterion removes **26,765**.
4. The remaining **2,649,654 cells** exactly match the cell IDs and raw marker
   matrix in `CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad`. The notebook also
   requires a marker to be positive in `int(0.001 * n_cells)` cells; all 31
   markers survive that filter.

The exact original import script was not identified. The TSV-to-H5AD selection
and value preservation above are verified from all rows, rather than inferred
from a later importer.

## 3. Protein transformation and phenotyping

`single_cell_analysis/CODEX_Panc_scimap_Analysis.ipynb`, cells 13–26, rebuilds
AnnData from the QC-retained **raw MFI**, retains it in `raw.X`, computes
`log1p(raw MFI)`, applies the manual gates using `scimap.pp.rescale`, and runs
`scimap.tl.phenotype_cells` with an external `phenotype_rules.csv`.

The saved source confirms these representations:

- `raw.X`: original, QC-retained whole-cell MFI, unchanged.
- `layers['log']`: `log1p(raw.X)`, verified over all retained marker values;
  maximum numerical discrepancy is 1.78 × 10^-15.
- `X`: marker-specific, **per-image piecewise min–max rescaling of log MFI**.
  Values below the gate map to 0–0.5 and those above it to 0.5–1. The scaler
  anchors its split to the observed value nearest the supplied gate. Each
  image's observed endpoints determine its two linear scales. This is not
  scVI-decoded expression. Numerical reconstruction matched 31,000 sampled
  saved values within 8.89 × 10^-16.
- `uns['gates']`: 31 markers × 15 donors. All donor columns carry the same
  manual thresholds; provenance records manual adjustment on donor 6450.
  Examples in log1p-MFI units: INS 5.3, GCG 7.75, ECAD 9.2, CD45 6.9, CK19 8.25.
- `obs['phenotype']`: 21 saved labels, consumed directly by aggregation.
- `obsm['X_scVI']`: a separate, ten-dimensional saved latent representation.

The scVI notebook also creates CLR and scaled-CLR representations, but the
phenotyping notebook explicitly returns to `raw_mfi`. Neither CLR nor the
scVI latent defines the scimap marker values stored in `X`.

The exact historical phenotype rule matrix and complete edit chain remain
unrecovered. The checked-in rule CSV does not reproduce the saved labels.
`scripts/fix_tcell_subtypes.py` documents an additional CD8a/CD4 refinement at
the 0.5 rescaled threshold, and all saved T-cell-lineage labels are consistent
with that refinement. The existing outputs and gates are available; recovering
the exact historical label-generation recipe is the remaining lineage gap.
See [the method-recovery evidence](original_method_recovery_20260929.md).

## 4. Aggregation, spatial displays and denominators

`scripts/reaggregate_islets.py` calls
`islet_analysis/fixed_islet_aggregation.py` with separate core/peri aggregation,
`min_cells=0`, and `require_paired=True`. A retained islet needs at least one
QC-retained cell in each compartment. `Parent=Islet_N` means core;
`Parent=Islet_N_exp20um` means peri. This yields **5,214 paired islets**,
containing **138,718 core cells** and **214,011 peri cells**.

For every retained islet the code averages scimap-rescaled `X` across all core
cells, counts each saved phenotype, and divides phenotype counts by total
cells in that region. It separately averages `X_scVI`. All 5,214 served marker
means, core and peri phenotype counts, region totals and mean scVI vectors
reproduce from the canonical single-cell file, including the float32 casts
used for saved means. These means include all core cell types.

The direct QuPath branch covers 5,382 cores; **168** lack a corresponding
paired-cell H5AD row. Even among shared islets, **3,124** have different
QuPath versus QC/Parent-based core cell totals. The two pipelines use different
selection/membership and gating rules, so neither their denominators nor their
labels should be interchanged. For example, `6505_Islet_284` has:

- QuPath: 260 core cells, of which 58 are INS-positive.
- Canonical phenotyped H5AD: 245 core cells, of which 48 are labeled Beta cell.

`scripts/extract_per_islet_cells.py` copies `X`, labels, centroids and morphology
into 5,214 per-islet CSVs. Samples covering the largest islet in every donor
(13,065 cells) match the canonical marker values, labels and centroids.
`scripts/extract_per_donor_tissue.py` copies all 2,649,654 cells into 15 donor
CSVs; all five output fields were checked across every row and match.

`scripts/compute_neighborhood_metrics.py` computes peri phenotype counts and
proportions. Its enrichment score is `(observed - expected) / sqrt(expected)`,
where expected is peri cell total × phenotype fraction **pooled across all
2,649,654 cells in the 15-donor cohort**. All 20 served enrichment columns
reproduce from this pooled baseline. Nearest-cell distances use the mean core
cell centroid and the nearest peri cell of the selected phenotype. The separate
signed-distance tables inherit QuPath's distance to the islet boundary, with
negative values inside; a 1,000-row sample matches the saved source after the
exporter's float32 distance cast.

## 5. scVI's live non-trajectory dependency

In the frozen 15-donor root, the Spatial Leiden assignments use scVI:

`saved single-cell X_scVI -> mean over each core -> Harmony(imageid, theta=2)
-> 15-neighbor cosine graph -> Leiden at 0.3, 0.5, 0.8 and 1.0`.

This is Step 8 of `scripts/reaggregate_islets.py`.
`islets_core_clustered.h5ad` records `use_rep='X_scVI_harmony'`; every served
cluster assignment matches that file. The notebook recipe for the earlier
cell-level fit uses raw MFI in a `counts` layer, a negative-binomial likelihood,
ten latent dimensions, and categorical Age/Gender covariates. Its executable
setup omits a donor `batch_key`, despite comments claiming one. Saved encoder
weights and a complete executed training record were not recovered, so those
notebook settings should not be presented as independently verified historical
fit parameters.

The displayed Spatial UMAP uses a different representation: PCA with 15
components on the per-islet mean **scimap-rescaled markers**, a 30-neighbor
Euclidean graph, then UMAP (`min_dist=0.5`, `spread=2`). The script calls this
“raw marker PCA,” but its actual input is the rescaled marker mean matrix.
These coordinates are copied into the clustered file and then into the app;
the copies were verified numerically. Thus cluster membership and displayed
UMAP coordinates come from different feature representations.

## Other inputs and evidence

Donor status, age, sex and autoantibody flags enter through metadata attached to
the canonical H5AD; they are not QuPath intensity measurements. The app builder
uses that H5AD metadata, not the currently present `CODEX_Pancreas_Donors.xlsx`
(which now lists a different donor set). Additional nPOD clinical fields come
from `data/donor_metadata.csv` through `scripts/build_donor_metadata.py`.

The measurement comparisons are saved in
`tmp/nontrajectory_lineage_evidence_20260929.json`, produced by the temporary
read-only audit scripts under `tmp/`. No served data, original exports, model
weights or application code were changed during this trace, and no training
job was submitted.
