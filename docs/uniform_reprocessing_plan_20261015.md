# Uniform reprocessing plan — all 20 donors — 2026-10-15

**Why.** The owner requires every dataset to come from a reproducible
workflow. The [file search](owner_directive_20261015.md#original-15-donor-files)
found that the phenotype rule file behind the original 15 donors' saved labels
and the trained scVI model behind their latent are lost. The original
islet-boundary recipe also cannot be rerun as it was. As directed, the original
15 donors will be reanalyzed the same way as the five additions, so all 20
donors pass through one scripted workflow from their saved cell detections to
the app root. Cell segmentation is already complete for every donor and is
reused, not rerun.

**Status (updated 2026-10-15):** owner directives later the same day decided D1 and D2: every image is re-segmented with InstanSeg 0.1.1 on a 9-channel Senior-style set, and Senior's newest (v3) islet method is applied to every image, reusing the existing Tissue outlines. That work runs through [`scripts/uniform/`](../scripts/uniform/README.md). D6 is decided (the app stays offline); D3–D5 remain open.

**Canary, 2026-10-15.** Aab_6505 and ND_6232 completed every step, but the islet classifier fit missed most beta cells (8 % and 52 % of the earlier-typed beta cells inside the fitted islets). Two causes, both in the port's cell test rather than in Senior's method: PGP9.5 had been mapped to Senior's Neural lineage although beta cells carry it most strongly, and on raw CODEX counts E-cadherin outranks insulin, so insulin was almost never a cell's highest marker. The port now leaves PGP9.5 out of the identity markers and reads the highest marker with each marker divided by its own 99.9th percentile in the image, insulin's and somatostatin's never below the cohort median (glucagon keeps the image's own: alpha cells are plentiful in every image, and a glucagon floor penalized the dimmer acquisitions). This adaptation awaits the owner's review with the islet objects. The first canary run is kept under `superseded/canary-v1-20261015/`.

**Islet run, 2026-10-15.** All 20 images were segmented again; 19 have islet objects ready for review in `data/development/uniform-20donor-20261015/qupath/project.qpproj`. T1D_6563's fit was refused by Senior's minimum of 50 confirmed cells (its alpha cells carry little CHGA) and needs an owner decision. Step 1b drops endocrine cells that cross script 1's outline by more than 0.75 µm, which costs from about 1 to 43 points of beta-cell coverage depending on the image. See [the README](../scripts/uniform/README.md#islet-run-on-all-20-images-2026-10-15).

**CHGA aid, 2026-10-16.** At the owner's direction CHGA no longer gates which cells train the islet classifier: a hormone (INS, GCG or SST) as the cell's highest marker confirms it, and CHGA only adds cells beside confirmed ones. All 20 images were refitted and now have islet objects, including T1D_6563 (glucagon and somatostatin). Coverage of the earlier endocrine cells rose on every image, with somewhat wider islet margins. See [the README](../scripts/uniform/README.md#chga-is-an-aid-not-a-gate-owner-2026-10-16). All 20 await the owner's review.

**Whole-cell objects, 2026-10-16.** At the owner's direction step 1b now builds every endocrine object from whole cell
segmentations: no cell cut by an edge, no cell in two objects, no object smaller than a cell, no Hormone Fragments. All
20 images pass `check_islet_objects.py`. The owner reviewed the objects the same day: "The islet annotations look good."
See [the README](../scripts/uniform/README.md#endocrine-objects-are-made-of-whole-cells-owner-2026-10-16).

**Cell typing, 2026-10-16.** As the owner decided (D3-D5, A6), every cell is typed with Islet-Explorer-Senior's
method, ported in `scripts/uniform/typing/`: compartment REDSEA, UniFORM alignment and the cohort scale, context
from the endocrine and structure objects, Senior's structure, endocrine, immune and Exocrine calls, refinement and
immune subtypes. Nerve, blood-vessel and lymphatic cells come only from the structure objects. The first run on
all 20 images is in `$OUT/typing/`; see [the README](../scripts/uniform/README.md#cell-typing-islet-explorer-seniors-method-owner-2026-10-16).

**App root, 2026-10-16 (P4 built, awaiting review).** `scripts/uniform/assemble_uniform.py` builds the whole
cohort from the uniform run, with nothing appended to the frozen release:
`data/development/uniform-20donor-20261016-v1/`. It holds 20 donors and 10,113 endocrine objects (1,491 islets,
3,415 clusters and 5,207 singles), of which 9,437 have analysed cells in both core and 20 um band. That is 3.08
million analysed cells, with 28,247 Excluded and 320,737 without a nucleus. The app serves all three object
classes behind an app-wide selector (A5). Local checks: R tests pass. The browser harness passes 37 of 45 checks,
including the 20-donor census, no donor 78 and the Harmony method. The four Viewer checks need the Range-capable
image route that only the nginx deployment provides. The other four (PAGA render and download, Plot vocabulary,
Trajectory contrast) also failed on v3, because the harness predates the 2026-10-14 app change. The service
stays stopped until the owner accepts the root (P6); the build report is `tmp/uniform-20donor-20261016-v1/`.

## Cohort

| Status | Original 15 | Additions |
|---|---|---|
| ND | 112, 6356, 6479, 6516, 6548 | 6232, 6539 (recovered) |
| Aab+ | 6450, 6505, 6521, 6538, 6549 | 6090 (recovered), 6558 |
| T1D | 6533, 6534, 6550, 6551, 6563 | 6520 |

Donor 78 (HDL078 / HPAP-078) is excluded (`scripts/excluded_donors.csv`). The
T1D group has six donors; the others have seven. The input list is
`scripts/protein_expansion_cohort.csv`, extended with each donor's QuPath
source acquisition path and SHA-256.

## The workflow

This is the scripted workflow already used for the five additions, run the same
way for every donor:

1. **Image intake.** Verify each source acquisition's channels (31 analysis
   markers plus DAPI and the unused channels) and the 0.3774 µm calibration;
   build viewer OME-TIFFs with `scripts/cohort_expansion/convert_images.py`.
2. **Cells (re-segmented, 2026-10-15 directive).** Every image is segmented again with InstanSeg
   `fluorescence_nuclei_and_cells` 0.1.1, the newest cell-and-nucleus model, on 9 channels (DAPI, CD68,
   CD44, VIM, ECAD, LGALS3, panCK, CD3e, BActin: Senior's 13 that exist in this panel), inside each
   image's existing Tissue outline, with Senior's detection script (duplicate checks) and measurement
   script. The earlier segmentations (0.1.0 for the original 15, verified intact; 0.1.1 on 22 channels for
   the additions) stay in their projects as references.
3. **Islets (2026-10-15 directive).** Senior's v3 method on every image: a per-image islet pixel
   classifier fitted on the image's own cells (CHGA in CD99's role; PGP9.5 not an identity marker; the
   highest marker read on each marker's image top, insulin's and somatostatin's never below the cohort median), script 1 (size-adaptive
   erode/dilate, 50 µm² minimum) and step 1b, which builds every object from whole cells and labels
   it Islet, EndocrineCluster or EndocrineSingle. You review all 20 in QuPath
   before acceptance. As in Senior v3, QuPath draws no peri-islet bands; the app's 20 µm peri region is
   computed downstream. Structure objects (nerve, capillary, lymphatic) are not yet in the uniform chain.
4. **Raw import.** `prepare_analysis.py --retain-unmeasured`: raw mean intensity,
   `log1p`, centroids, and a census of every object.
5. **Preparation.** Fixed QC (at least 2 detected markers; total raw MFI ≥
   5234.796272, the original cohort's 1st percentile), `log1p`, and per-donor
   rescaling against the 31 recovered gates
   (`prepare_protein_expansion.py`).
6. **Typing.** `protein_phenotyping.py` with the accepted
   `protein_phenotype_rules.csv`.
7. **Assembly.** Regions (tissue minus exclusions; a 20 µm peri band that
   excludes all cores), cell membership, hormone positivity (mean + 1.5 SD of
   raw INS/GCG/SST per donor), regional CD3e/CD163/CD68/LGALS3 positivity
   (mean + 1.5 SD per region), composition, neighborhood and distance tables,
   per-islet and per-donor cell files, Harmony protein PCA (θ 2, seed 0) with
   UMAP and Leiden, and the trajectory composition features.
8. **Validation.** Release validator, browser check, and a rerun to confirm
   reproducibility.

## What changes for the original 15 donors

Every value currently shown for them is regenerated:

- New cell objects, IDs and measurements from the new segmentation.
- New islet boundaries and `Islet_N` names. The count will differ from 5,214.
- Cell types come from the accepted rule table, not the saved scimap labels.
  That table does not assign Acinar, ECAD+ or Structural (see D3).
- Marker values use per-donor gate rescaling instead of scimap's per-image
  scaling. Hormone and regional immune positivity are recomputed in Python
  rather than taken from the historical QuPath exports.
- Spatial UMAP and Leiden come from Harmony protein PCA; scVI is not used.

The frozen 15-donor release and the protein-expanded roots stay on disk as
comparison baselines. They are not deleted and are not served once the new root
is accepted.

## Decisions needed first

| | Decision | Recommendation |
|---|---|---|
| D1 | Cell segmentation | **Decided 2026-10-15:** InstanSeg 0.1.1 on 9 channels for all 20 images; the additions are re-run too, since they used 22 channels. |
| D2 | Islet boundaries | **Decided 2026-10-15:** Senior's v3 islet method on all 20 images, CHGA in CD99's role, existing Tissue outlines reused; owner review of all 20. |
| D3 | Cell-type rules | **Decided 2026-10-16:** use Islet-Explorer-Senior's cell typing as far as this panel allows. Earlier recommendation: freeze one table before typing. Either keep `protein_phenotype_rules.csv` as accepted on 2026-10-03, or revise it once (for example explicit Acinar, ECAD+ and Structural rules, or the deferred CHGA/CD45/VIM/LGALS3/CD44/PGP9.5/GAP43/CD35/NaKATPase enhancement) and accept it again. |
| D4 | Hormone and regional immune positivity | **Decided 2026-10-16:** as in Senior's typing. Earlier options: keep mean + 1.5 SD for every donor, which preserves the current definitions, or switch to gate-based positivity with the recovered gates so these measures match typing. |
| D5 | Marker gates | **Decided 2026-10-16:** as in Senior's typing. Earlier recommendation: keep the 31 recovered thresholds as fixed parameters. They were set interactively in November 2025; re-gating would add a new manual step. |
| D6 | What the app serves meanwhile | **Decided 2026-10-15: offline** until the reprocessed root is accepted (P6). |
| A3 | Endocrine object categories | **Decided 2026-10-16: keep** EndocrineSingle = one cell; EndocrineCluster = two or more cells, mean diameter under 35 µm; Islet = two or more cells, 35 µm or more. |
| A5 | App analysis unit | **Decided 2026-10-16:** all three categories, with a filter. |
| A6 | Structures | **Decided 2026-10-16:** drop the cell classes that stand for nerve, blood vessel or lymphatic and use structure objects instead; the five additions already have such objects. |

## Engineering before any compute

- **E1 Manifest-driven runners.** Replace the hard-coded case maps in
  `run_qupath.py`, `convert_images.py`, `prepare_qupath.py` and `resume_gpu.py`
  with the cohort manifest, so the original donors take exactly the same path.
- **E2 Whole-cohort assembly.** Have preparation and assembly build all 20
  donors from prepared inputs instead of appending additions to the frozen
  root. The parity audit against the frozen root becomes a comparison report.
- **E3 Run records.** Each stage writes to a fresh directory and records input
  SHA-256 values, the git commit, parameters, software versions and output
  SHA-256 values.
- **E4 Environment lock.** Lock the Python environment and record the QuPath
  0.7.0 build, extension versions, and InstanSeg and SAM checkpoint checksums.
- **E5 Determinism.** Fix seeds and thread counts. If D1 re-segments the five
  additions, measure rerun differences on one donor, GPU against CPU
  (`--device cpu` exists), and use CPU segmentation if exact reruns are
  required.
- **E6 Tests.** Cover the manifest-driven path and exclusion guard on synthetic
  fixtures, plus a one-crop end-to-end smoke test.

## Execution

| Phase | Work | Gate |
|---|---|---|
| P0 Freeze | Manifest with checksums, D1–D6 recorded, environment lock | Owner sign-off |
| P1 Pilot | Full workflow for one original donor (6505) and one addition (6232), plus a rerun; compare 6505 with the frozen release | Owner reviews the pilot |
| P2 Re-export and regions | Re-export saved cells and measurements for the remaining 18 donors and run the region and boundary stages, one at a time (the project lock is per run); re-segment the five additions only if D1 requires it | Each export manifest complete; original-donor values identical to `all_data.tsv` |
| P3 Boundary review | Owner reviews all 20 donors; reviewed geometry exported; acceptance record written | Owner acceptance |
| P4 Build | Import, preparation, typing and assembly into a new root such as `data/development/uniform-20donor-<date>-v1` | Validator passes |
| P5 Validation | Rebuild from the accepted P3 outputs and compare hashes; browser check; comparison report against the frozen release and v3 | All checks pass |
| P6 Acceptance | Owner scientific review; switch the app; mark the older roots historical; update docs | Owner acceptance |

## Acceptance criteria

- All 20 donors and no excluded donor appear in every app table, file and image
  list.
- Each donor's census reconciles every exported object as included, below QC or
  unmeasured.
- Rebuilding P4 from the accepted P3 outputs reproduces every app artifact's
  SHA-256, or matches within a documented tolerance where GPU arithmetic is
  involved.
- Every stage has a run record naming code commit, inputs, parameters, software
  and outputs.
- Boundaries for all 20 donors and the rule table are accepted by the owner.
- Differences from the frozen release are reported per donor (cells, islets,
  label composition) for information; matching the old values is not required.

## Effort

- Engineering E1–E6: about 3–5 working days.
- Compute: segmentation is not rerun for the original 15. The additions'
  automated QuPath stages took about 15–40 minutes per donor, of which
  segmentation was 5–11, so re-export and boundary stages for 20 donors should
  fit in about a day. Import through assembly takes under an hour.
- Owner review of 20 donors' boundaries is the longest step and depends on the
  reviewer's time.

## Risks

- Original-donor cells and raw measurements are kept, but islet boundaries
  (D2), cell types and every derived value change. Figures, statistics or
  manuscripts based on the frozen data would need regenerating.
- Without Acinar, ECAD+ and Structural rules, exocrine composition is coarser
  than in the frozen labels unless D3 adds them.
- If D1 re-segments the additions, InstanSeg on GPU may not rerun bit-for-bit;
  E5 decides between a tolerance and CPU execution.
- The SAM step depends on the separate `samapi` environment and a local VIT_B
  checkpoint; both must be locked (E4).
