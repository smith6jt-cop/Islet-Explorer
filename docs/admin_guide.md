# Islet Explorer — Administration & Operations Guide

Developer- and maintainer-facing reference for rebuilding the app data, deploying the Shiny app, and troubleshooting. End-user feature documentation lives in [user_guide.md](user_guide.md).

**Requested reviewer URL (pending UFRC production-access confirmation)**:
`https://islet-explorer.rc.ufl.edu/islet-explorer/`

> **Deployment authority:** [deploy/README.md](../deploy/README.md) and
> [PUBAPPS_UFRC_HANDOFF.md](../PUBAPPS_UFRC_HANDOFF.md) define the August 2026 reviewer
> deployment. Any older mutable-code, direct-nginx, Apptainer, or AI-key instruction in
> historical notes is not a production procedure.

> **Current migration status, 2026-10-16:** The local app is offline pending
> uniform reprocessing of the 20-donor target; donor 78/HDL078 is permanently
> excluded. The PubApps instance request form was submitted with ticket `83810`
> on 2026-08-31, and UFIT canceled individual assessment `IRM-456246` on
> 2026-10-15 while routing the project to the master assessment. UFRC comment 13
> says VM provisioning is next and asks for DNS/proxy requirements; the reply is
> drafted in the [launch-preparation note](ufrc_launch_preparation_20261015.md).
> See the [reprocessing plan](uniform_reprocessing_plan_20261015.md) for the
> processing gates. Do not start or deploy a local candidate before the
> reprocessed root is accepted.

---

## Data Pipeline

The 2026-09-08 source preparation distinguishes the 5,214-observation trajectory
H5AD from prepared app tables: `core_area` and `comp` each have 5,382 rows,
`targets_all` has 48,438, `markers_all` has 64,584, and `combined_pairs` has 21.
The counts in the historical pipeline below refer to their named artifacts,
not to every runtime summary. Use [the evidence index](ufrc_readiness.md) and
the generated release manifest for current release validation.

### Rebuilding the App Data

The canonical pipeline rebuilds `data/islet_explorer.h5ad` from source:

```bash
# Activate the Python environment
conda activate scvi-env

# Step 1: Reaggregate islets + trajectory + Leiden (~15 min)
python scripts/reaggregate_islets.py
# → islet_analysis/islets_core_fixed.h5ad (5,214 islets)
# → data/adata_ins_root.h5ad (+ pseudotime + UMAP)
# → islet_analysis/islets_core_clustered.h5ad (+ Leiden at 4 resolutions)

# Step 2: Compute neighborhood metrics (from 2.6M-cell single-cell H5AD)
python scripts/compute_neighborhood_metrics.py
# → data/neighborhood_metrics.csv (5,214 rows × ~100 cols, per-phenotype enrich_z_* and min_dist_*)

# Step 3: Extract per-islet cell CSVs (for drill-down viewer)
python scripts/extract_per_islet_cells.py
# → data/cells/*.csv (5,214 files, ~203 MB)

# Step 4: Extract per-donor tissue CSVs (for Spatial tab scatter, independent of islet filtering)
python scripts/extract_per_donor_tissue.py
# → data/donors/*.csv (15 files, ~78 MB)

# Step 5: Build enriched H5AD (trajectory + groovy + neighborhood + Leiden + donor metadata)
python scripts/build_h5ad_for_app.py
# → data/islet_explorer.h5ad (~70 MB)
```

### Full Lineage

```
CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, canonical single-cell)
  ↓ scripts/reaggregate_islets.py (min_cells=0, require_paired=True)
islets_core_fixed.h5ad  (5,214 islets, proteins + scVI embeddings + trajectory + Leiden)
  ↓ scripts/build_h5ad_for_app.py (+ groovy + neighborhood + donor metadata)
data/islet_explorer.h5ad  (complete app data, incl. Leiden clustering)
```

Branch pipelines:
- `scripts/compute_neighborhood_metrics.py` — peri-islet composition, immune metrics, enrichment z-scores, distances
- `scripts/extract_per_islet_cells.py` — individual cell CSVs for drill-down
- `scripts/extract_per_donor_tissue.py` — per-donor tissue CSVs for Spatial tab scatter
- Leiden clustering from `islet_analysis/islets_core_clustered.h5ad` (4 resolutions) — merged by `build_h5ad_for_app.py`

### Upstream Data Sources

| Source | Location | Description |
|--------|----------|-------------|
| Single-cell H5AD | `single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad` | 2.6M cells, scVI-embedded and phenotyped; downstream pipeline artifacts retain audited legacy Harmony/DPT experiments |
| Groovy TSV exports | `panc_CODEX/results/groovy_exports/` | QuPath islet measurements (15 donors × 4 types) |
| GeoJSON boundaries | `data/json/*.geojson` | Islet segmentation polygons |
| Spatial lookup | `data/islet_spatial_lookup.csv` | Islet centroid coordinates |

See `data/DATA_PROVENANCE.md` for the complete lineage documentation including validation results.

---

## Deployment & Administration

### Architecture

```
reviewer HTTPS + Basic Auth
  ↓ UFRC high-availability reverse proxy
Requested PubApps backend 127.0.0.1:3838
  (final address and port pending UFRC confirmation)
  ↓ rootless Podman forwarding
immutable Shiny container 0.0.0.0:3838
  ├── manifest-backed data + images (read-only)
  └── bounded temporary files + operational logs (read/write)
```

### Updating the App

Never edit or `git pull` into an active reviewer release. Build and scan a clean,
digest-pinned OCI image; build and fully checksum-verify a new data release; then use
`deploy/activate-release.sh`. It health-checks the new Quadlet service and restores the
prior unit on failure. Use `deploy/rollback-release.sh` for an explicit rollback. See the
deployment runbook for complete commands and required evidence.

### Local Development

```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
# Then open http://localhost:7777 in browser
```

Port 7777 is for local testing only and is not accessible to end users.

### R Dependencies

```
shiny, shinyjs, plotly, ggplot2, dplyr, tidyr, readxl, sf, jsonlite,
RColorBrewer, scales, anndata, reticulate, lmerTest, lme4, emmeans
```

### Python Dependencies

```bash
conda activate scvi-env
# scanpy, anndata, scvi-tools, sklearn, scib-metrics, scipy, pandas, numpy
```

The production app environment is independently pinned by
`deploy/requirements-python.lock`; do not substitute the mutable pipeline environment in
the reviewer image.

---

## Troubleshooting

### Production app fails to load

1. Run `deploy/verify-deployment.sh` with the active release root and expected image
   reference.
2. Inspect the sanitized service journal and Shiny log; do not print environment values or
   reviewer request URLs containing query strings.
3. Verify the active data manifest in full mode. Production must fail closed on a missing
   artifact or hash mismatch; do not enable the development fallback.

### Local development app fails to load

1. Check that `data/islet_explorer.h5ad` exists and is ~70 MB
2. Check that R packages are installed: `library(anndata); library(reticulate)`
3. Check that Python is accessible via reticulate: `reticulate::py_available()`
4. For development only, the Excel fallback requires `data/master_results.xlsx` and has
   reduced functionality.

### No phenotype or demographic filters

The app is running from Excel fallback. Rebuild the H5AD:
```bash
conda activate scvi-env
python scripts/build_h5ad_for_app.py
```

### No neighborhood metrics in Plot selector

1. Verify `data/neighborhood_metrics.csv` exists (5,214 rows)
2. Rebuild the H5AD: `python scripts/build_h5ad_for_app.py`
3. Kill stale R workers and refresh the browser

### Single-cell drill-down not available

1. Verify `data/cells/` directory contains ~5,023 CSV files
2. If missing, regenerate: `python scripts/extract_per_islet_cells.py`

### Segmentation viewer shows "No GeoJSON found"

1. Check that `data/json/` or `data/gson/` has `.geojson` / `.geojson.gz` files
2. Case ID zero-padding: GeoJSON files use 4-digit padded IDs (`0112.geojson`), data uses unpadded (`112`). The app handles this automatically via `sprintf("%04d", ...)` fallback.

### Statistics tab shows no results

The Statistics tab uses data from the Plot sidebar. Ensure:
1. A feature is selected in the Plot sidebar
2. At least 2 disease groups are checked
3. The diameter range includes some islets

### Spatial tab tissue scatter is empty

1. Verify `data/donors/` directory contains 15 CSV files
2. If missing, regenerate: `python scripts/extract_per_donor_tissue.py`
3. If Leiden panel says "not available", rebuild the H5AD with Leiden data: `python scripts/build_h5ad_for_app.py` (requires `islet_analysis/islets_core_clustered.h5ad`)

### AI assistant shows error or no response (development only)

The reviewer release intentionally has no AI panel or key; this is expected. Do not add a
Navigator/OpenAI credential to production. Development AI configuration is outside the
approved reviewer data flow and requires a separate local secret plus appropriate review.

### Leiden UMAP shows a blob

The Leiden UMAP visualization uses raw marker PCA coordinates (same as the Trajectory tab). If it shows a blob, the clustered H5AD may have outdated scVI-based UMAP coordinates. Fix:
1. Copy visualization UMAP from `data/adata_ins_root.h5ad` into `islet_analysis/islets_core_clustered.h5ad`
2. Rebuild: `python scripts/build_h5ad_for_app.py`
3. Restart the local development app. For production, build and activate a new immutable
   release through `deploy/README.md`; do not replace code inside the running container.

### Drill-down does nothing or Viewer flickers on macOS (fixed Jun 2026)

Earlier builds had two macOS-only problems: clicking a plot point did nothing (macOS suppresses plotly's `plotly_click` event), and the Viewer tab flickered/reloaded continuously (the WebGL `<iframe>` was destroyed and recreated on every reactive update — heavy flicker under Firefox's WebRender compositor on macOS, Mozilla bug #1555544). Both are fixed:

- **Drill-down** clicks are read from native pointer events by the linked-plot script (`app/shiny_app/www/trajectory-links.js`), so they work even when `plotly_click` is suppressed; the Plot and Trajectory plots route through the shared `select_islet_from_key()` resolver.
- **Viewer** mounts the iframe once and updates its `src` in place via `session$sendCustomMessage`, so the WebGL canvas is never recreated.

If these symptoms reappear on macOS, confirm the active release manifest points to the
expected image digest and hard-refresh the browser. Do not kill individual production R
workers or modify container code in place; activate a rebuilt release through the runbook.
On the Plot tab, enable **Show individual points** so there are per-islet points to click.
If Viewer flicker remains specific to a Retina display, investigate Avivator's
device-pixel backing-buffer behavior in a separate tested image build.
