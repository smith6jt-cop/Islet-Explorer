# Owner directives — 2026-10-15

> Donor HDL078 needs to be removed completely from the app and all analysis.
> Then the files used to generate the data for the first 15 donors needs to be
> located to process the additional 5 in the same way. If the files cannot be
> located, a plan to reanalyze the first 15 in the same way as the additional 5
> needs to be generated since we cannot have data without a reproducible
> workflow/pipeline.

Outcome: donor 78 is removed and the app serves a rebuilt 20-donor root. The
files behind the original 15 donors were only partly located. The phenotype
rule file and the trained scVI model are lost, so the five additions cannot be
processed by the original method. The plan is in
[the uniform reprocessing plan](uniform_reprocessing_plan_20261015.md).

## Donor 78 removal

**Identity.** The repository's case 78 is the only donor matching HDL078. Its
acquisition file is `CX_PANT-HPAP-078-T1D_test28_34p_16LED50.tif` and its
staging records call it HPAP-078, T1D, female, age 1.21. The nPOD export
(`All-Download_From_nPOD_2026-07-27_163844.xlsx`) lists case HDL078 as No
Diabetes, female, 1.69 years. The removal does not depend on which label is
correct; the discrepancy is recorded in `scripts/excluded_donors.csv`.

**What changed**

- `scripts/excluded_donors.csv` is the single record and
  `scripts/excluded_donors.py` the shared guard. IDs compare after dropping an
  HPAP/HDL prefix and leading zeros, so 78, 0078, HPAP-078 and HDL078 all match.
  A loader that meets an excluded donor stops; it does not filter silently.
- Protein pipeline: `stage_original_gate_calls.DONORS` lists 6090, 6232, 6520,
  6539 and 6558. Preparation, assembly and the release validator reject
  excluded donors and no longer assume six additions or 21 donors. The 20-donor
  manifest is `scripts/protein_expansion_cohort.csv`.
- Other cohort loaders: `research/common.read_cohort()` rejects excluded donors,
  so the historical 21-donor `research/cohort-v1.csv` and the schema-4/phase-2
  tools built on it no longer load. `cohort_expansion/cohort.csv` marks 78
  `excluded` and its reader requires that marking; the importer already refuses
  non-selected cases. `build_donor_metadata.py` rejects it. The completed QuPath
  runners (`run_qupath.py`, `convert_images.py`, `prepare_qupath.py`,
  `resume_gpu.py`) no longer offer `0078`.
- App: `ISLET_EXCLUDED_DONORS` in `00_globals.R` makes
  `load_shared_prepared_data()` refuse a data root that contains donor 78, and
  the Viewer never lists its image. Cohort notes now read 20 donors (7 ND, 7
  Aab+, 6 T1D) and five additions.
- Browser check (`deploy/browser/check-app.mjs`): expects 20 donors and the five
  additions and fails if 78 appears in the donor selector or image inventory.
- Tests: `scripts/tests/test_excluded_donors.py` covers the registry, every
  guard and the app constant. Synthetic fixtures that used case 78 now use other
  IDs. Python: 264 passed. R app suite: 233 of 233 assertions passed.

**Rebuild and verification.** The 20-donor root
`data/development/protein-expanded-20261015-v3/` (7,028 islets) was assembled
with the existing pipeline; commands and results are in
[the protein expansion record](protein_expansion_20261003.md#five-addition-rebuild--2026-10-15).
The release validator passed, and a scan of every table, H5AD, workbook, file
name, image link and acceptance record found no donor-78 entry. The local app
served v3 (selected in `tmp/local-app.env`; previous file
`tmp/local-app.env.before-20261015`) until it was taken offline later the same
day. The browser check passed 41 of 45; the four failures do not involve donor
78 and are explained in the expansion record.

**Kept, unused.** Source and historical material for donor 78 is still on disk
and is no longer read by any current pipeline or by the app:

| Material | Location | Size |
|---|---|---:|
| Source acquisition and QuPath entry 42 | `/data/QuPath_Projects/QuPath_project_IsletExplorer/Images/New/CX_PANT-HPAP-078-T1D_test28_34p_16LED50.tif` | 7.8 GB |
| Viewer image | `data/development/expanded-20260911/images/T1D_0078.ome.tiff` | 7.1 GB |
| QuPath exports and raw cells | `data/development/expanded-20260911/exports/0078/`, `analysis/raw_cells/78.h5ad` | 2.7 GB |
| Raw-only staging copy | `data/development/donor-additions-20260923-raw-v1/additions/raw/78.h5ad` | 585 MB |
| Superseded protein inputs and roots | `protein-additions-20261003-v1`/`v2`, `protein-review-geometry-20261003-v1`/`v2`, `protein-expanded-20261003-v1`/`v2` | about 6 GB |
| Historical candidates and reviews | `data/development/all-donor-2026092*`, `islet-regeneration-20260917*`, several `data/reviews/*` folders | several GB |

The owner decided to keep this material on disk, unused. The QuPath review
projects under `data/reviews/` and the original QuPath project were not edited.

## Follow-up decisions — 2026-10-15

1. **Donor 78's files are kept, unused.** Nothing listed above is deleted; the
   current pipelines and the app refuse the donor.
2. **The local app is offline until the reprocessed root is accepted.**
   - `islet-explorer-local.service` is stopped and disabled
     (`systemctl --user disable --now`). Nothing listens on port 3840 and the
     endpoint returns 502.
   - nginx (`/etc/nginx/conf.d/islet-explorer.conf`, root-owned) also served
     `/islet-explorer/local_images/` straight from
     `data/development/all-donor-20260920/images/`, bypassing Shiny. That folder
     holds all 21 earlier viewer images, so `T1D_0078.ome.tiff` was still
     downloadable there (HTTP 206 on 2026-10-15) after the app stopped listing
     it. Other-user access to the folder was removed (`chmod o-rwx`); nginx runs
     as `www-data` and now returns 404 for every file there.
   - Root access was not available to edit the nginx config. With `sudo`, disable
     that server block or make it return 503, so the closure does not depend on
     folder permissions, and point any future static image route only at an
     accepted root.
   - To restore later, only after the reprocessed root is accepted: set
     `ISLET_DATA_ROOT` in `tmp/local-app.env`, run
     `systemctl --user enable --now islet-explorer-local.service`, and retarget
     the nginx alias (or `chmod o+rx` the old folder if it is still wanted).

## Original 15-donor files

Each stage is listed with what was found. "Located" means the file or the exact
executed code exists now.

| Stage | Original method | Status | Evidence |
|---|---|---|---|
| Source images | 15 PhenoCycler TIFF/OME-TIFF | Located | QuPath project `Images/`, 15 entries |
| Cell segmentation | QuPath v0.6, InstanSeg `fluorescence_nuclei_and_cells` 0.1.0, 512-px tiles, 32-px padding; 22- or 35-channel input lists | **Detections located and intact**; settings recorded | The project's 15 images hold exactly the 2,676,466 cells of the original export, including all 2,676,419 imported cell IDs (checked 2026-10-15), so segmentation need not be rerun. Settings are in the per-image `.qpdata` workflows in `/data/QuPath_Projects/QuPath_project_IsletExplorer/data/`. A model 0.1.0 copy exists under `/data/KINTSUGI_*/notebooks/instanseg/bioimageio_models/`. QuPath 0.6 is absent and its InstanSeg extension version is unrecorded. The additions used QuPath 0.7 with model 0.1.1. |
| Tissue, islet and structure annotations | QuPath pixel classifiers (Tissue, Islets, Nerve, Capillary, Lymphatics), refinement, 20 µm expansion | Partly located; not regenerable | Classifier JSONs exist, but `Islets.json` and `Tissue.json` were re-saved on 2025-09-03 during the export run and earlier versions were not kept. Manual edits are unrecorded. The saved boundaries exist as GeoJSON. |
| Regional QuPath summaries | `3_Islet_Counter`, `4_Islet_Cells`, `5_Islet_Objects` | Scripts present but edited later | `3_Islet_Counter.groovy` changed 2026-03-09 and `1_Panc_Objects.groovy` 2026-07-30, after the Aug–Sep 2025 exports. Stored thresholds reproduce from the cell export ([trace](current_nontrajectory_processing.md)). |
| Whole-cell export | QuPath export to `all_data.tsv` (2,676,466 objects) | Output located; export step unrecorded | `/data/panc_CODEX/results/old_data_plots/all_data.tsv` |
| Import to `panc_CODEX.h5ad` | Mean columns; drop DAPI, epCAM, KRT8-18 and unmeasured rows; merge donor sheet | Code located | HiPerGator home IPython session 89 (2025-09-16) and the executed `/blue/maigan/smith6jt/KINTSUGI_old/notebooks/Pancreas_CODEX.ipynb` |
| QC | At least 2 markers; total MFI at or above the pooled 1st percentile (5234.796272) | Located | HiPerGator session 79; reproduces the 2,649,654 retained cells |
| scVI latent (`X_scVI`) | scvi-tools 1.4.0; NB likelihood; 10 latent dimensions; Age/Gender covariates; seed 0 | Training code located; **trained model lost** | `/orange/maigan/smith6jt/Xenium/scratch/tmp/ipython/profile_default/history.sqlite`, session 79 (2025-11-13), saved the model to `/blue/clive/smith6jt/KINTSUGI/notebooks/CODEX_scvi_model_v2/`. That folder no longer exists; it appears to have been dropped in the 2026-03-23 `/blue/clive` to `/orange` move. A retrain would be a new fit, not the saved latent. |
| Gating | Interactive napari gates (set on 6450), 31 thresholds applied to every image | Located (values) | `uns['gates']` in the gated and phenotyped H5ADs; local IPython sessions of 2025-11-16 and 17 |
| Rescaling | `scimap.pp.rescale(adata, gate=gates)` per image | Call located; scimap version unrecorded | Local IPython sessions 672, 680, 681 |
| Cell phenotypes | `scimap.tl.phenotype_cells` with `/home/smith6jt/panc_CODEX/phenotype_rules.csv`, edited between runs; the final run on 2025-11-18 ("newDuctal") made the saved labels | **Rule file lost** | That path became a symlink to the app's `data/phenotype_rules.csv` on 2026-08-25. No copy exists in git, VS Code history (workstation or HiPerGator), LibreOffice backups, agent transcripts, the lab share, Google Drive or HiPerGator. Only its structure is known, from scimap's log on 2025-11-17: top level ECAD+, SMA+, Endothelial, Immune, Endocrine, Exocrine, Structural, Neural; Immune to B cell, T cell, Macrophage, APCs; T cell to CD8a Tcell, CD4 Tcell; Endocrine to Beta, Alpha, Delta; Exocrine to Ductal, Acinar; Endothelial to Lymphatic, Blood Vessel. That matches the archived `phenotype_rules2.csv` template. The checked-in Rules1 file, and the 1,394-byte copy on HiPerGator, is a different table. |
| T-cell refinement | `scripts/fix_tcell_subtypes.py` (2026-03-24) | Located | Repository; it explains the only label differences between the 2025-11-18 output and the frozen file |
| Aggregation, neighborhoods, app build | `reaggregate_islets.py`, `compute_neighborhood_metrics.py`, extraction scripts, `build_h5ad_for_app.py` | Located | Repository (restored from `000e970`); these consume the lost labels and scVI latent |

The five additions therefore cannot be processed by the original method. Cell
segmentation is the exception: the original detections are intact, and the
additions were segmented with model 0.1.1. Its trained network is identical to
0.1.0 (all 677 learned tensors match); only the built-in post-processing
changed. The [InstanSeg v0.1.1 release](https://github.com/instanseg/instanseg/releases/tag/instanseg_models_v0.1.1)
describes a merging change to avoid overflow errors and an Apple-GPU fix; the
compiled code also turns fragment cleanup on by default, merges by overlap with
the smaller object instead of IoU, and drops 0.1.0's 10-pixel minimum-size
filter (diff: `tmp/original-pipeline-recovery-20261015/instanseg_0.1.0_vs_0.1.1_forward_code.diff`). The original islet boundaries cannot be regenerated
by script. The rule file that defines the saved cell types and the model behind
the saved latent are gone.

An exact-replay reconstruction of the rules from the saved labels and rescaled
values was considered and is not recommended. It would be a reverse-engineered
table rather than the original file, it cannot be shown to behave identically on
new donors, and it would leave the other differing stages unchanged.

**Search coverage.** Workstation: name searches of `/home`, `/data` and `/tmp`;
the IPython history database (18,750 inputs, Nov 2024 to Oct 2026); VS Code local
history; LibreOffice backups; shell histories; Claude and Codex transcripts; git
history on all branches; and hard-linked copies of the frozen H5AD. Lab share
`/mnt/ahc_share/PERSONAL/smith6jt` to depth 7. Google Drive title and full-text
search. HiPerGator: home, `/blue/maigan/smith6jt` and `/orange/maigan/smith6jt`
(bounded), three IPython history databases and the VS Code server history. The
only `model.pt` files there are Xenium gene models. Not searched: other users'
folders, several `/orange/maigan` group folders, and `/orange/clive` folders
that denied access.

Evidence copies (IPython sessions, notebook snapshots, HiPerGator extracts) are
in `tmp/original-pipeline-recovery-20261015/` with `SHA256SUMS`.
