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
cohort from the uniform run, with nothing appended to the frozen release. The current build is
`data/development/uniform-20donor-20261019-v7/` (v1-v6 are superseded by the fixes below). It holds 20 donors
and 10,113 endocrine objects (1,491 islets, 3,415 clusters and 5,207 singles). Islets and clusters with analysed
cells in both core and 20 um band are the 4,763 analysis rows on v5 (1,487 islets, 3,276 clusters). Singles are only
counted (v4, below). That is 3.07 million analysed cells, with 33,843 Excluded (11,095 for an unresolved hormone, v5) and
320,737 without a nucleus. The app serves Islets and Clusters behind an app-wide selector (A5), which since 2026-10-17 also
replaces every minimum-cells control. The Trajectory instead fits every object with ten or more core cells,
whatever its class (owner, 2026-10-18). Local checks on v4 (2026-10-18): 232 R tests pass, and the browser
harness passes all 45 checks against a shiny-server instance with the local service's configuration
(`deploy/shiny-server-local.conf` on another port), whose image route answers Range requests; `shiny::runApp()`
cannot, so the four Viewer checks fail there. The PAGA and Plot vocabulary checks were brought up to the
2026-10-14 app (PAGA built on request; hormone-positivity fractions withdrawn). The local env file now sets
`RETICULATE_PYTHON` to the app venv, which has `fa2-modified` (the scvi-env interpreter in the service unit does
not, so PAGA silently fell back to Fruchterman-Reingold). The build report is `tmp/uniform-20donor-20261018-v5/`; on v5 the harness again passes 45 of 45.

**Owner ruling, 2026-10-19 (v6).** The UMAPs separated disease poorly. The investigation found the original
15-donor UMAP's separation was donor batch (72 % same-donor neighbours; same disease among other-donor neighbours at
chance; donor explained 86 % of HLADR variance, disease 1 %), while v5's protein + Harmony layout lost much of the
real disease difference, T1D beta-cell loss, because Harmony on donor also removes donor-level disease shifts. On
the owner's choice ("ten-cell cutoff", composition refit), the UMAP and Leiden are now the endocrine composition of
the 1,215 objects with ten or more core cells, with no donor correction (`scripts/uniform/composition_embedding.py`;
measurements in [the README](../scripts/uniform/README.md)). v6 was rebuilt by the assembler from the same run
(report `tmp/uniform-20donor-20261019-v6/`) and is identical to v5 apart from the embedding: every other table and
boundary file has the same content, and the protein record differs only by Harmony's run-to-run rounding (at most
0.0001; the protein Leiden assignments are identical).

**v7, 2026-10-19.** The v6 composition UMAP was one dense blob with two detached islands. The islands are real
objects (31 all-beta small islets, mostly ND and Aab+; 21 all-alpha insulin-deficient islets, mostly T1D; each far
above chance from the donors' pooled make-up, spread over 10-12 donors, not fragments of one islet and not a typing
artefact) but their detachment came from identical features: tied points pull together in a 30-neighbour graph. On
the owner's direction (adjust min_dist and spread; points must not appear to overlap), the UMAP is refitted at
min_dist 3, spread 3 on the same graph, so those objects are the ends of the cloud and on-screen overlap falls from 42 %
to 7 % of points within half a marker of a neighbour. The fit now runs on one thread, because the thread count changed
the layout and Leiden (adjusted Rand 0.63); v6's Leiden cannot be reproduced exactly, and v7's clusters agree with
v6's at adjusted Rand 0.57-0.73. The app draws each UMAP as one interleaved trace, so no disease state is painted over
another. v7 was rebuilt by the assembler (report `tmp/uniform-20donor-20261019-v7/`): X, every table and every
boundary file have v6's content; only `X_umap`, `leiden_umap_1/2` and `leiden_<r>` (and the record-only protein UMAP)
differ.

**The two open items, 2026-10-17** (owner: "use common sense and available tools to resolve the scaling and
threshold"). (1) Aab_6558's Excluded cells inside objects: its INS stain is weak (2.4x its tissue p99 in the
objects against 15-86x on the other ND and Aab+ images; C-peptide 8 ng/mL), and the tissue-p99 hormone level cut
through its beta cells. Each hormone is now read over its own tissue background, and inside the objects the level
is lowered to a local FDR of 10 % (Beta 722 -> 1,769; Excluded 35 % -> 21 % of object cells; well-stained images
gain 0.2-3 %). (2) T1D_6520's lymphatics at 2.6x the earlier area: the additions' structure thresholds had been
moved by UniFORM factor ratios, the worst transfer in leave-one-out on the 15. Nerve and Capillary now use the
15's raw median unscaled; Lymphatic keeps the original random-forest classifier's objects (T1D_6520 back to its
audited 43,000 um2). Details in [the README](../scripts/uniform/README.md).

**Owner rulings, 2026-10-17 (v3).** "Blood vessel, lymphatic, and nerve should not be cell types, but B cell should
be": structure cell calls are off (the structures stay as objects and areas), and B cells are called on CD20 with
CD45 in place of CD79a (2-795 per image). "Remove the minimum cells control sliders and wire all functions to the
new classifications": the Plot, Spatial and Trajectory minimum-cells controls are gone; the Islet / Cluster /
Cluster selection decides which objects every tab uses. Follow-up 2026-10-18, "Use 10 or more cells as the
cutoff instead of a category": the Trajectory fits every object with ten or more core cells (1,265 objects:
1,242 Islets and 23 Clusters; ρ with the beta-cell ratio −0.70, against −0.20 on all objects), and the class
selector is hidden on that tab.

**Owner ruling, 2026-10-18 (v4).** "There should not be any singles in the core. Singles should not be included
in any islet-level analysis except for quantifying or showing how many there are since a single cell cannot also
include immune cells or vessels/nerves/lymphatics." Singles are no longer rows, cores, bands or signed-distance
targets; their typed cell is an ordinary tissue cell (1,264 of them lie in an islet's or cluster's 20 um band and
count there). `endocrine_objects_by_donor.csv` counts every class per donor and per mm2 of tissue (median Singles
per mm2: ND 10.0, Aab+ 19.0, T1D 32.7), shown in the Plot tab's "Endocrine objects per donor" card; the
segmentation map draws them. Of the 5,207 Singles, 4,720 cells are typed (1,081 Beta, 2,533 Alpha, 1,106 Delta),
472 are Excluded and 15 have no nucleus. The trajectory still fits 1,265 objects (rho -0.69, ND 0.30 < Aab+ 0.38
< T1D 0.74).

**Owner ruling, 2026-10-18 (v5).** "These hormones are mutually exclusive so there needs to be resolution or
exclusion based on the approved islet cell typing methods. The 'only' and 'at all' designations are not biology,
they are indicators of technical problems that need to be solved. We need to add a CHGA+INS-GCG-SST- cell type."
The second hormone is neighbour bleed (see [the README](../scripts/uniform/README.md)); a hormone now names a cell
only at 2x the runner-up (Senior's confirmation margin), else the cell is Excluded (11,095), and CHGA+
hormone-negative is a cell type (2,350; mean share of Islet endocrine cells ND 2.3 %, Aab+ 2.8 %, T1D 4.7 %). The
trajectory fits 1,215 objects (rho -0.67; ND 0.29 < Aab+ 0.36 < T1D 0.73).

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
   UMAP and Leiden (kept for the record), the served UMAP and Leiden on the
   endocrine composition of objects with ten or more core cells (2026-10-19),
   and the trajectory composition features.
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
- Spatial and Trajectory UMAP and Leiden come from endocrine composition of the
  objects with ten or more core cells (owner, 2026-10-19; Harmony protein PCA
  before that); scVI is not used.

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
