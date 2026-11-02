# Islet Explorer — Administration & Operations Guide

Developer- and maintainer-facing reference for rebuilding the app data, deploying the Shiny app, and troubleshooting. End-user feature documentation lives in [user_guide.md](user_guide.md).

**Current workflow, 2026-11-02:** use [the complete workflow user guide](workflow_user_guide.md)
and [current handoff](handoff_20261102.md) for processing, component selection and release state.
The Lymphatic correction freezes calibration and automatically audits all components; the old
all-marker rerun is not a structure-only procedure. Dev last served unaccepted v20; dev and
production are now stopped and persistently masked during
[security remediation](container_security_20261102.md). The staged whole-v19 reference was
never activated. V19 endocrine is selected within the
[component comparison](v19_v20_component_comparison_20261102.md); no corrected mixed root exists.

**Provisioned reviewer URL (production currently stopped)**:
`https://islet-explorer.rc.ufl.edu/islet-explorer/`

> **Deployment authority:** [deploy/README.md](../deploy/README.md) and
> [PUBAPPS_UFRC_HANDOFF.md](../PUBAPPS_UFRC_HANDOFF.md) define the reviewer
> deployment contract and its dated history. Any older mutable-code, direct-nginx, Apptainer, or AI-key instruction in
> historical notes is not a production procedure.

The provisioned endpoints do not imply an active service. The owner stopped dev on
2026-11-02; production has been stopped since 10-26. Preserve both service masks during
the current remediation. A security replacement image is not an accepted scientific
release or an instruction to restart. Earlier October provisioning notes are historical;
use the current [security record](container_security_20261102.md) and
[readiness record](ufrc_readiness.md).

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
labels wholesale or mix already aggregated v19/v20 summaries. No corrected mixed
candidate has been assembled or accepted.

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

Use the workflow guide for commands and prerequisites. The successful diagnostic in
[the component comparison](v19_v20_component_comparison_20261102.md) uses initial-v20
calibration to test isolation; it does not implement the selected v19 endocrine component.
Do not rerun the historical all-marker v20 driver for a structure-only revision.

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
  (both currently stopped and masked)
  ↓ rootless Podman forwarding
immutable Shiny container 0.0.0.0:3838
  ├── manifest-backed data + images (read-only)
  └── bounded temporary files + operational logs (read/write)
```

### Updating the App

Neither service is authorized to restart during remediation. Preserve its mask; an older
activation or rollback command does not authorize unmasking or restarting it.

For a future authorized release, never edit or `git pull` into an active release. Build and scan a clean,
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
scientific status in review notes. This does not start the stopped in-house or PubApps
services or establish an unaccepted root as a release.

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

1. Confirm whether a service is expected to be running. Both are currently stopped and
   masked; do not restart them to troubleshoot. For a future active release, run
   `deploy/verify-deployment.sh` with its release root and expected image reference.
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
