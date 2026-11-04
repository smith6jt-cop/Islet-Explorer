# Islet Explorer — Administration & Operations Guide

Developer- and maintainer-facing reference for rebuilding the app data, deploying the Shiny app, and troubleshooting. End-user feature documentation lives in [user_guide.md](user_guide.md).

**Current workflow, 2026-11-03:** use [the workflow guide](workflow_user_guide.md)
and [selected-component v2 build record](component_build_v2_20261102.md).
Scientific producer `67ee6c9` preserves all 110,721 selected named v19 endocrine
calls across 20 donors and 3,199,606 eligible cells. Its controlled DG revision
has 353 permitted Lymphatic-related label changes and no unrelated changes.
The final 7,398 Unclassified cells retain their evidence; 716,451 Tissue calls
remain unchanged. Zero unnamed endocrine definitions or boundary placeholders
remain. The classifier-lineage proof passes while corrected-cell and centroid
diagnostics remain visible; they are not evidence of negative original pixels.

Production `pubapps-2026-11-ui-refresh-1659c9d` is live from `1659c9d`, image
`16c4979f…`. Its full data parity passed 25/25 and authenticated proxy access
passed 15/15; the current service, hashes and checks are bound in the
[production receipt](evidence/pubapps_live_20261104_ui_refresh.json). The
separate strict-DuckDB dev release `pubapps-2026-11-dev-duckdb-v1` remains live
from `378843e`, image `e704ac2a…`; see its [receipt](evidence/dev_duckdb_deployment_20261103.json).
The November 3 display fixes and card controls are included in production. The
local preview has been stopped; its image, tests and visual observations are
preserved in the [local review record](local_duckdb_review_20261103.md). The
older `islet-explorer-local.service` remains inactive. A designated clean
off-campus acceptance run and VM reboot test remain pending.

**Live reviewer URL (authentication required):**
`https://islet-explorer.rc.ufl.edu/islet-explorer/`

> **Deployment authority:** [deploy/README.md](../deploy/README.md) and
> [PUBAPPS_UFRC_HANDOFF.md](../PUBAPPS_UFRC_HANDOFF.md) define the reviewer
> deployment contract and its dated history. Any older mutable-code, direct-nginx, Apptainer, or AI-key instruction in
> historical notes is not a production procedure.

The receipts establish actual active services, not merely provisioned endpoints.
Production retains its dated DuckDB exception through **2026-11-10 UTC**;
the separate dev fix passes the strict decoder gate without an exception.
Neither security nor browser checks replace scientific validation. The separate
outside-network acceptance and reboot tests remain pending; see [UFRC readiness](ufrc_readiness.md).
Historical failed reports and earlier v1 diagnostics remain preserved.

---

## Data Pipeline

### Current uniform workflow

Follow [workflow_user_guide.md](workflow_user_guide.md) and
[scripts/uniform/README.md](../scripts/uniform/README.md). The old H5AD rebuild below is
historical compatibility processing, not the uniform workflow.

The frozen v19/v20 roots each contain 20 donors and 3,199,606 eligible cell UUIDs. Their
assembled metadata records 4,848 analysed Islet/Cluster objects, 5,227 Singles counted
separately and 1,392 saved embedding objects. These counts belong to specific artifacts;
upstream geometry and filtered app views need not have the same row count. Use each
root's `uniform_release.json` and label groups. Donor 78 must remain absent.

A coherent candidate needs one identity per UUID and explicit component provenance.
The owner selected v19 endocrine; v20 implements later measurement and structure rules
but also contains the demonstrated Lymphatic calibration regression. Do not copy v20
labels wholesale or mix already aggregated v19/v20 summaries. The component builder has
regenerated the selected-component v2 app and validated its exact UUID/identity
contract. Unsupported rows remain Unclassified; original-classifier lineage and
Single whole-cell association are the release object proof. The older strict
corrected-cell diagnostic remains separately visible.

For a **Lymphatic-only update**:

1. Prepare a fresh root from an explicitly selected baseline, preserving measurements,
   QC and endocrine membership.
2. Bind it with `scripts/uniform/typing/structure_update.py`. Its immutable calibration
   context is separate from the new live Lymphatic geometry.
3. Regenerate permitted structure context and endpoints. Reuse calibration;
   `levels.py` and `canary.py` refuse recalibration of a bound run.
4. Run `identity.py` for every bound image into its fresh identity directory. The CLI
   automatically runs `audit_structure_update.py` in a fresh external report directory.
   Named endocrine identities, marker positivity/canary/folds, measurements/QC and
   unrelated context must stay unchanged; changed labels must involve Lymphatic cells
   within the old/new footprint.
5. Review image evidence and reconcile selected components. Passing the isolation audit
   establishes consistency, not scientific acceptance.

Use the workflow guide for commands and prerequisites. The earlier diagnostic in
[the component comparison](v19_v20_component_comparison_20261102.md) used initial-v20
calibration to test isolation. The later [component build](component_build_20261102.md)
preserves the selected v19 endocrine contract and tests the DG change against its matched
control. Do not rerun the historical all-marker v20 driver for a structure-only revision.

After reconciliation, rebuild summaries, marker scales, embeddings, spatial tables and
geometry mappings from the coherent candidate. Endocrine-object and release guards must
pass; do not hide unnamed cells by broadening an aggregate. Full hashes and source/release
parity establish packaging fidelity, not biological accuracy.

### Historical 15-donor pipeline — compatibility only

The 2026-09-08 source preparation distinguishes the 5,214-observation trajectory
H5AD from prepared app tables: `core_area` and `comp` each have 5,382 rows,
`targets_all` has 48,438, `markers_all` has 64,584, and `combined_pairs` has 21.
The counts in the historical pipeline below refer to their named artifacts,
not to every runtime summary. Use [the evidence index](ufrc_readiness.md) and
the generated release manifest for current release validation.

#### Rebuilding the historical app data

The historical pipeline rebuilt `data/islet_explorer.h5ad` from its older source labels.
The following commands are retained for reproducing that artifact, not repairing a current
uniform root:

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

#### Historical lineage

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

#### Historical upstream data sources

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
Provisioned VM private listeners: production :3838 / dev :3839
  (both active under their separately verified release receipts)
  ↓ rootless Podman forwarding
immutable Shiny container 0.0.0.0:3838
  ├── manifest-backed data + images (read-only)
  └── bounded temporary files + operational logs (read/write)
```

### Updating the App

Production and dev have completed their recorded validation and activation.
Production retains the dated November 10 DuckDB deadline; dev passes the strict
decoder gate. Future updates still need their own matching image/data evidence
and must preserve the prior release for rollback.

For the authorized release, never edit or `git pull` into an active release. Build and scan a clean,
digest-pinned OCI image; build and fully checksum-verify a new data release; then use
`deploy/activate-release.sh`. It health-checks the new Quadlet service and restores the
prior unit on failure. Use `deploy/rollback-release.sh` for an explicit rollback. See the
deployment runbook for complete commands and required evidence.

### Local Development

```bash
cd app/shiny_app
ISLET_DATA_ROOT=/absolute/path/to/explicit-review-root Rscript -e 'shiny::runApp(".", port = 7777)'
# Then open http://localhost:7777 in browser
```

Port 7777 is for local testing only. Select the data root explicitly and retain its
scientific status in review notes. These example commands do not change the older,
inactive `islet-explorer-local.service`, the separate persistent, enabled review
service on loopback port 3845, or either active PubApps service, or establish an
unaccepted root as a release.

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
the reviewer image. The scVI packages above describe legacy pipeline capability; the
frozen v19/v20 manifests explicitly record `scvi_required: false`.

---

## Troubleshooting

### Production app fails to load

1. Confirm the selected instance and its recorded release. Both PubApps instances are
   active under their separate receipts; inspect the current state before changing it.
   Run `deploy/verify-deployment.sh` with that release root and expected image reference.
2. Inspect the sanitized service journal and Shiny log; do not print environment values or
   reviewer request URLs containing query strings.
3. Verify the active data manifest in full mode. Production must fail closed on a missing
   artifact or hash mismatch; do not enable the development fallback.

### Local development app fails to load

1. Confirm the explicit `ISLET_DATA_ROOT` and its manifest/provenance.
2. Check required H5AD, summary, geometry, donor and cell artifacts against that root's
   inventory. Historical fixed file sizes are not a completeness check.
3. Check the pinned R/Python environment and the error log.
4. The Excel loader is historical development compatibility only; production must not
   fall back to a different dataset.

### No phenotype or demographic filters

Check that the loaded root and prepared summaries match. Current phenotypes follow the
frozen identity/label-group contract, not the old 19-rule CSV. Missing fields require a
source correction and a coherent rebuild; do not fill them from another version's tables.

### No neighborhood metrics in Plot selector

Check the selected root's neighborhood artifacts and prepared summary fields, then the
feature's current tab placement. Correct source mismatches through the uniform workflow
and recheck parity; the historical 5,214-row file is not the current contract.

### Single-cell drill-down not available

Check the selected object identifier against the root's cell/geometry inventory and
production spatial tables. Do not infer completeness from an old fixed CSV count or
repair a uniform root with the historical extraction script.

### Segmentation viewer shows "No GeoJSON found"

1. Check the selected root's `json/` or `gson/` and its release manifest.
2. Case ID zero-padding: GeoJSON files use 4-digit padded IDs (`0112.geojson`), data uses unpadded (`112`). The app handles this automatically via `sprintf("%04d", ...)` fallback.

### Statistics tab shows no results

Statistics is disabled in the publication profile. If explicitly enabled for local
development, its historical module uses the Plot sidebar. Ensure:

1. A feature is selected in the Plot sidebar
2. At least 2 disease groups are checked
3. The diameter range includes some islets

### Spatial tab tissue scatter is empty

Check the selected donor, current filters and root-specific donor/spatial inventory.
The current comparison cohort has 20 donors. If Leiden data are missing, inspect the
recorded embedding cohort and artifacts before rebuilding through the uniform workflow.

### AI assistant shows error or no response (development only)

The reviewer release intentionally has no AI panel or key; this is expected. Do not add a
Navigator/OpenAI credential to production. Development AI configuration is outside the
approved reviewer data flow and requires a separate local secret plus appropriate review.

### Leiden or trajectory data appear inconsistent

Check the root's recorded embedding method, eligible object cohort and identity hashes.
Frozen v19/v20 use endocrine composition; historical raw-marker/scVI layouts are different
artifacts. Do not copy coordinates from another H5AD to make a plot appear correct.
Rebuild derived outputs from the chosen identities and rerun parity. Visual separation
alone does not establish correct processing or biological progression.

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
