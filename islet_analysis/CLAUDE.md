# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains analysis code for pancreatic islet CODEX (CO-Detection by indEXing) multiplexed imaging data. The primary goal is to analyze islet composition, microenvironment, and progression trajectories in Type 1 Diabetes (T1D) using single-cell protein expression data from nPOD (Network for Pancreatic Organ Donors with Diabetes) samples.

## Data Architecture

### Input Data
- **Single-cell AnnData**: `CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad` in `../single_cell_analysis/`
  - 2.6M+ cells × 31 protein markers
  - Contains scVI latent space in `.obsm['X_scVI']`. **Caveat:** the model in
    `SCVI_CODEX_v2.ipynb` was trained with `categorical_covariate_keys=["Age","Gender"]`
    but *no `batch_key`*, so donor-specific technical variance was **not** removed.
    Downstream clustering in `scripts/reaggregate_islets.py` Step 8 applies Harmony
    (`batch=imageid`, θ=2) on `X_scVI_mean` to correct this before Leiden.
  - Cell phenotypes in `.obs['phenotype']` (Alpha cell, Beta cell, Acinar, Macrophage, etc.)
  - Spatial coordinates in `.obs['X_centroid']` and `.obs['Y_centroid']`

### Region Definition
Cells are categorized by the `Parent` column in `.obs`:
- **Islet core**: `Parent` starts with `Islet_` but does NOT end with `_exp20um`
- **Peri-islet**: `Parent` ends with `_exp20um` (20μm expansion around islets)
- Islet IDs are unique per image: combine `imageid` + `base_islet_id`

### Output Datasets
- `islets_core_fixed.h5ad`: Aggregated islet core data (5,214 islets, min_cells=0, require_paired=True)
- `islets_peri_fixed.h5ad`: Aggregated peri-islet microenvironment (5,214 islets)
- `islets_merged_fixed.h5ad`: Paired islet+peri data for comparative analysis (5,214 islets)

## Key Analysis Patterns

### AnnData Structure (Islet-level)
```python
adata.X                              # Mean protein expression only (n_islets × 31)
adata.obsm['X_scVI_mean']            # Mean scVI latent representation (use for neighbors)
adata.obsm['phenotype_proportions']  # Cell type proportions per islet
adata.obsm['phenotype_counts']       # Cell type counts per islet
adata.obsm['morphology']             # Mean morphological features
adata.uns['phenotype_names']         # List of phenotype category names
adata.obs['donor_status']            # ND, Aab+, or T1D
```

### Creating Islet Datasets
```python
from fixed_islet_aggregation import create_islet_dataset_fixed

# Separate core and peri-islet (recommended)
adata_islet, adata_peri = create_islet_dataset_fixed(
    adata_scvi,
    region='separate',
    min_cells=0,         # Keep all islets (default, no size filter)
    require_paired=True  # Must have both core and peri-islet
)
```

For production rebuilds, prefer `scripts/reaggregate_islets.py` which runs
aggregation + trajectory + Leiden in a single script.

### Trajectory Analysis
```python
# Neighbors from scVI latent (used for trajectory/DPT; NOT batch-corrected).
# For Leiden clustering the reaggregate pipeline applies Harmony first.
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_scVI_mean', metric='cosine')
sc.tl.paga(adata, groups='donor_status')
sc.tl.umap(adata, init_pos='paga')

# Pseudotime: root at ND islet with highest INS
ins_idx = list(adata.var_names).index('INS')
nd_mask = adata.obs['donor_status'] == 'ND'
root_idx = np.where(nd_mask)[0][np.argmax(adata.X[nd_mask, ins_idx])]
adata.uns['iroot'] = root_idx
sc.tl.diffmap(adata, n_comps=10)
sc.tl.dpt(adata)
```

## Key Protein Markers

| Marker | Cell Type/Function |
|--------|-------------------|
| INS | Insulin (Beta cells) |
| GCG | Glucagon (Alpha cells) |
| SST | Somatostatin (Delta cells) |
| CD45 | Pan-immune |
| CD3e, CD4, CD8a | T cells |
| CD68 | Macrophages |
| panCK, CK19 | Epithelial/Ductal |
| SMA | Smooth muscle/myofibroblasts |
| PGP9.5 | Neural |

## Disease Progression Insights

From pseudotime analysis:
- **Decreasing** with progression: INS, PGP9.5, Beta cell proportion
- **Increasing** with progression: GCG, CD3e, CD4, Alpha cell proportion, Ductal cells
- T1D islets show beta cell loss and immune infiltration in peri-islet regions

## Enhanced Trajectory Analysis

The notebook includes improved trajectory analysis based on published methods:

### Methods Comparison
| Method | Used For | Reference |
|--------|----------|-----------|
| Diffusion Pseudotime (DPT) | Core trajectory | Scanpy default |
| Palantir (optional) | Fate probabilities | [Nature Biotech 2019](https://github.com/dpeerlab/Palantir) |
| PAGA | Trajectory topology | Scanpy |

### Validation Metrics
The trajectory is validated by checking:
1. Pseudotime correlation with disease stage (ND → Aab+ → T1D)
2. Biological markers: INS decreases, GCG increases along trajectory
3. Diffusion map eigenvalue gap (trajectory structure quality)

### Visualization Outputs
- `pseudotime_distributions.png`: Violin plots and diffusion map by disease state
- `expression_trends_pseudotime.png`: Smooth protein expression trends
- `phenotype_trends_pseudotime.png`: Cell type proportion changes
- `expression_heatmap_pseudotime.png`: Ordered heatmap of all markers
- `phenotype_by_disease_*.png`: Box plots with statistical tests

## Running Analysis

The main analysis notebook is `Panc_codex_islet_analysis_fixed.ipynb`. The workflow:
1. Load single-cell scVI data
2. Create separate islet core and peri-islet datasets using `fixed_islet_aggregation.py`
3. Compare core vs peri-islet composition (with statistical tests)
4. Compute embeddings and trajectory on islet core
5. Validate trajectory with biological markers
6. Analyze protein-pseudotime correlations with smooth trends
7. Generate comprehensive visualizations
8. Save processed datasets

## Statistical Methods

- **Group comparisons**: Kruskal-Wallis test (non-parametric) with pairwise Mann-Whitney U
- **Correlations**: Spearman rank correlation with p-values
- **Multiple testing**: Benjamini-Hochberg FDR correction
