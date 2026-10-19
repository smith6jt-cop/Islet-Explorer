# CODEX uniform segmentation and islets

Owner directives of 2026-10-15: every cell and nucleus segmentation comes from the most recent InstanSeg model
with a small channel set, as in Islet-Explorer-Senior, and the newest islet method is applied to every image.
Donor 78 is excluded (`../excluded_donors.csv`). This folder applies one scripted chain to all 20 images;
it replaces both the original 15 donors' 2025 segmentation (InstanSeg 0.1.0) and the additions' 22-channel
segmentation.

## Decisions

| Item | Choice | Source |
|---|---|---|
| Model | InstanSeg `fluorescence_nuclei_and_cells` 0.1.1, the newest cell+nucleus model (v0.1.2 added only a DAPI nuclei model) | owner, 2026-10-15 |
| InstanSeg channels | DAPI, CD68, CD44, VIM, ECAD, LGALS3, panCK, CD3e, BActin: Senior's 13 with a CODEX counterpart; CD99, HLA-A, CD79a and Iba1 are not in this panel and are dropped, as IO60pLN dropped CD99 | owner, 2026-10-15 |
| Islet method | Senior v3: per-image islet pixel classifier (`islet_classifier_fit.py`), script 1 (size-adaptive erode/dilate, 50 µm² minimum), step 1b (objects made of whole cell segmentations; Islet / EndocrineCluster / EndocrineSingle) | owner, 2026-10-15; step 1b rule 2026-10-16 |
| CD99's role in the islet classifier | CHGA | owner, 2026-10-15 |
| Tissue outlines | Reused unchanged: the original project's Tissue for the original 15, the six-donor review project's for the additions | owner, 2026-10-15 |
| Peri-islet bands | None in QuPath, as in Senior v3 (per-cell signed distances instead); the app's 20 µm peri regions are computed downstream | Senior v3 |

## What differs from Senior v3 (CODEX parameters only)

Copied scripts and their source checksums are listed in `PORTED_FROM.json` (Senior scripts commit `340bf0c`).

- Pixel size 0.3774 µm (Senior 0.5078) in preflight, detection and script 1.
- The islet classifier works on pyramid level 1 of these images, 0.7548 µm (Senior 1.0156 µm). Pixel-count
  parameters keep their physical size: seed 8 px = 6.04 µm (Senior 6.09), closing 4 px = 3.02 µm (Senior 3.05);
  opening and margin stay one pixel (0.75 µm; Senior 1.02). Step 1b's fit tolerance is one classifier pixel
  (0.7548 µm).
- The fit reads its four planes from QuPath (`qupath/export_classifier_planes.groovy`) instead of a `.qptiff`,
  uses QuPath's raw `Cell: <marker>: Mean` (this pipeline has no REDSEA correction), and its identity markers are
  Senior's lineage table on this panel (`islet_cell_test.LINEAGE_MARKERS`). Structure objects are optional and
  are not used yet.
- PGP9.5 is not an identity marker (`islet_cell_test.SHARED_WITH_ENDOCRINE`): on this panel it is brightest in
  beta cells. Senior's Neural lineage is B3TUBB alone; GAP43, which beta cells do not carry, stays beside it.
- The fit reads "the cell's highest identity marker" with each marker divided by its own 99.9th percentile in the
  image (`islet_cell_test.common_scale`), as Senior's typing lane scales markers before its argmax. Raw CODEX
  counts differ by two orders of magnitude between channels (E-cadherin p99.9 about 37,000, insulin about 800 on
  6505), so on raw counts E-cadherin, not insulin, is highest in beta cells. Insulin's and somatostatin's tops are
  never below the cohort's median top (`islet_cell_test.FLOORED_HORMONES`, `cohort_scale.py`,
  `fit/cohort_scale.json`). A p99.9 lands on a hormone's cells only where they are more than 1 in 1,000 of the
  image's cells; in the earlier typing beta cells fall below that in four T1D images (0.02–0.24 per 1,000) and delta
  cells in 14 of the 20 (0.04–0.91). There the image's own top is background, and the division turns it into
  "hormone highest": on its own insulin top T1D_6534 had 8,255 insulin-highest cells more than 30 µm from any earlier
  endocrine cell, 23 with the floor (T1D_6551: 77,508 and 11). Glucagon keeps the image's own top. Alpha cells are
  5–35 per 1,000 cells in every image, T1D included, and a glucagon floor, the first rule (its all-image run was
  stopped after two images), only penalized the dimmer acquisitions: the additions' glucagon tops are 9,249–17,425
  against a median of 18,288, and the floor halved the share of earlier alpha cells with glucagon highest on ND_6232
  (53 → 25 %) and Aab_6090 (57 → 23 %). Every other marker keeps the image's own top, so channel differences between
  acquisitions cancel (a cohort-wide scale for every marker made 2,008 such cells on ND_6232, whose SMA, VIM and CD3e
  are much dimmer than in the original images). A scale built under another rule is rebuilt and every fit made on it
  reruns. CHGA is still compared with its level on raw values. These changes followed the first canary (below) and
  await the owner's review.

## CHGA is an aid, not a gate (owner, 2026-10-16)

Senior's fit trains on "confirmed" endocrine cells: a hormone as the cell's highest identity marker **and** its CD99
(here CHGA) at or above the image's tissue level. With CHGA in that role the gate discarded most hormone-highest cells
in every image (Aab_6505 kept 669 of 8,033; ND_6548 168 of 1,616) and refused T1D_6563 outright (45 of 3,007: its
alpha cells carry glucagon, median about 4,000, but little CHGA, median 132). On the original 15 images 94–99.5 % of
hormone-highest cells lie within 10 µm of a cell the earlier typing called endocrine (T1D_6533 81 %). The owner's
direction: "CHGA is supposed to be an aid to endocrine identification, not a blocker ... GCG and SST can be used."

So a hormone as the cell's highest marker (on the per-image scale, insulin and somatostatin floored) now confirms it.
CHGA still adds hormone-less, CHGA-bright cells within 30 µm of a confirmed cell, and stays a channel the classifier
may use where a hormone seeds it. An image with no usable CHGA level is fitted on its hormones alone. A fit made by
other fit code is stale (`run_uniform.py` compares the recorded script checksums), so the change reran every image's
islet steps; cells and segmentation were not touched. The 2026-10-15 islet outputs are kept in
`$OUT/superseded/chga-gate-20261015/`.

Rerun of 2026-10-16 (`islets_all_20_v3.log`): **all 20 images have islet objects.** T1D_6563 fits on glucagon and
somatostatin (8 insulin-highest cells; the earlier typing had 3 beta cells): 42 islets, 95 clusters, 226 singles,
91 % / 89 % of the earlier alpha / delta cells inside the final objects. Earlier beta / alpha / delta cells inside the
final objects now 75–99 / 74–96 / 85–100 % (were 25–98 / 51–92 / 22–100 %); the largest gains are where the gate had
hurt most (beta: T1D_6551 25 → 75 %, Aab_6558 37 → 85 %, T1D_6550 35 → 77 %, Aab_6505 63 → 98 %, ND_6548 54 → 93 %).
Islets broken into many small objects joined up again on Aab_6505 and Aab_6558. The cost is a wider margin: in the
fit's emulation of script 1 the islet area grew 5–15 % on most images and about 50 % on Aab_6505 and Aab_6558; the
non-endocrine cells taken in lie a median 8–10 µm from an earlier endocrine cell. On T1D_6563 a third of them lie
farther than 20 µm, many earlier-typed Ductal: duct-associated alpha cells or over-reach, for the owner's review.

## Endocrine objects are made of whole cells (owner, 2026-10-16)

"The endocrine objects outlines need to also pick up the cell segmentations. There should not be multiple annotations
in one cell outline; among other problems this causes a single-cell to be counted as a cluster. There should never be
an annotation for endocrine object that is smaller than a cell."

Senior's step 1b took a cell only when it reached no farther than one classifier pixel beyond the outline, cut every
other cell out, and turned an outline holding no such cell into a Hormone Fragment shaped like the hormone's own pixels.
On the 2026-10-16 refit that left, per image, 62–522 Hormone Fragments, 42–307 objects holding no whole cell, 106–905
cells cut by an object's edge and up to 45 cells inside two objects (`check_islet_objects.py`; the objects are kept in
`$OUT/superseded/fragments-20261016/`). Step 1b now:

- picks up every cell an outline covers by at least half, whole, also where it reaches beyond the outline (a cell two
  outlines cover goes to the one covering more of it);
- lets an outline that picks up no cell that way (hormone signal inside one cell or across a few) pick up the cell
  holding most of it, unless another outline has that cell -- largest overlaps first, so one cell never holds two
  objects;
- removes an outline that ends up with no cell, and makes no Hormone Fragments;
- keeps pieces that hold parts of one cell together: about 7 % of InstanSeg cell outlines come in parts.

`check_islet_objects.py` reads the exports and fails an image with any object holding no whole cell, any cell in two
objects, any cell cut by an object's edge (5–95 % inside), any single that is not one cell, or any cluster / islet with
fewer than two.

Rerun of 2026-10-16 (`islets_all_20_v4.log`, fits unchanged): **all 20 images pass every check**
(`$OUT/review/object_check.csv`). Earlier beta / alpha / delta cells inside the final objects are now 80–100 / 88–99 /
93–100 % (were 75–99 / 74–96 / 85–100 %; T1D_6563 had 3 earlier beta cells). Clusters and singles roughly double, since
the former Hormone Fragments now pick up the cells under them (e.g. T1D_6533 100 → 316 clusters, 297 → 422 singles);
islets rise slightly (e.g. ND_0112 130 → 146).

Owner review, 2026-10-16: "The islet annotations look good." The four images saved during the review (Aab_6450, 6505,
6521, 6538) hold the same objects, class by class, as `$OUT/islets/`, so the exports stand for the reviewed project.

## Cell typing: Islet-Explorer-Senior's method (owner, 2026-10-16)

Owner decisions of 2026-10-16: D3, D4 and D5, "use the cell typing from Islet-Explorer-Senior as much as possible";
A6, "drop cell classes that indicate neural, blood vessel or lymphatic and use the objects instead". The port is
`typing/` (sources and checksums in `PORTED_FROM.json` under `typing`; the Senior files are byte-identical to the
snapshots of its 2026-10-08 run). It writes to `$OUT/typing/`:

| Stage | Script | Senior source | Output |
|---|---|---|---|
| Geometry QC, compartment REDSEA | `redsea_cells.py` | `phenocycler/redsea.py`, `geometry_qc.py` | `redsea/<image>.parquet` |
| UniFORM shift, cohort scale | `align_scale.py` | `celltype/run_step2.py`, `run_step3.py` | `aligned/`, `scale.parquet`, `shifts.json` |
| Context: endocrine object, structure object, tissue | `cell_context.py` | `typing/context.py`, `islet_membership.py` | `context/`, `$OUT/structures/` |
| Structure, endocrine, immune, Exocrine calls | `classify.py` | `typing/subpanel.py` | `calls/<image>.parquet`, `.json` audit |
| Refinement, once | `refine.py` | `typing/refine.py` | `refine/` (`refined_class` beside `cell_class`) |
| Immune class and subtype | `immune_subtype.py` | `typing/immune_subtype.py` | `immune/` |
| Additions' structure objects (owner) | `structure_thresholds.py fit / make` | -- | `$OUT/structures/` |
| Background check (owner) | `background_check.py` | -- | `background_check.csv` |
| Endocrine checks (owner) | `endocrine_checks.py` | -- | `endocrine_checks.json` |

```bash
PY=/home/smith6jt/miniconda3/envs/scvi-env/bin/python
SENIOR_PY=$HOME/IO60panc2nd/Islet-Explorer-Senior-v3/.venv/bin/python   # has scikit-fda, which UniFORM needs
$PY typing/redsea_cells.py --out $OUT --images all --workers 6        # ~1.5 min and ~12 GB per image
MPLBACKEND=Agg $SENIOR_PY typing/align_scale.py --out $OUT
$PY typing/structure_thresholds.py --out $OUT fit                     # once; then `make` for the five additions
$PY typing/structure_thresholds.py --out $OUT make
$PY typing/cell_context.py --out $OUT --images all
$PY typing/background_check.py --out $OUT
$PY typing/classify.py --out $OUT --images all
$PY typing/refine.py --out $OUT
$PY typing/immune_subtype.py --out $OUT
$PY typing/endocrine_checks.py --out $OUT
$PY typing/summary.py --out $OUT
```

What Senior's method means here:

- **Nerve, capillary and lymphatic structures are objects, not cell types** (A6; owner 2026-10-17: "Blood
  vessel, lymphatic, and nerve should not be cell types"). Senior's structure call is switched off
  (`codex_typing.STRUCTURE_CELL_CALLS = False`): a cell inside a structure object keeps its own identity through
  the endocrine, immune and Exocrine questions, and its centroid membership and covered fraction are still
  recorded (`context/`, `passed_through_*`). The structures enter the app as areas (the Microenvironment layer).
  Endothelial and nerve cells therefore fall to the residual Exocrine class, and inside an endocrine object to
  Excluded: Excluded rose 0.7-6 % per image (Aab_6505 2,000 -> 2,121). The objects are each image's Nerve, Capillary and Lymphatic annotations from the original
  project's pixel classifiers, copied to `$OUT/structures/` with their source and the 2026-09-21 audit verdict
  (`structures/provenance.json`). That audit rejected Nerve on Aab_6505 and T1D_6520 and reviewed few Capillary
  and Lymphatic sets.
- **Endocrine cells are called only inside the accepted endocrine objects**, on each hormone over its own
  ordinary-tissue background against the image's hormone level (a departure, below); an EndocrineSingle's cell is named by its highest hormone. Exocrine
  inside an endocrine object is Excluded (`exocrine_inside_islet_border`). The CHGA route (Senior's CD99) needs 500
  hormone-confirmed cells per islet-size stratum and fired in no stratum, so admission is the hormone route, as on
  96 of Senior's 100 donor strata.
- **Gates (D5)**: no fixed gates. UniFORM aligns each marker per image; one cohort floor/top scale puts it on
  [0, 1]; structure levels are each image's p99 of ordinary tissue, the hormone level that p99 lowered inside the
  objects to a local FDR of 10 % (below); immune rules use scaled >= 0.10.
- **Positivity (D4)**: hormone positivity is endocrine admission; immune positivity is a marker at scaled >= 0.10.
- **Labels (D3)**: Beta, Alpha, Delta, Endocrine (unresolved / ambiguous subtype), Immune (with its definition,
  class and subtype), Exocrine (one label, the residual: Senior never calls acinar or ductal; endothelial and
  neural cells fall here too), Excluded, and after refinement "Unresolved (no lineage evidence)".

Departures, all for this panel (`typing/codex_typing.py`):

- CHGA in CD99's place; Capillary objects (drawn on CD31) are Senior's Blood vessel, read on CD31.
- Compartments for this panel's own markers by analogy (granule and cytoskeletal markers Cytoplasm, surface
  markers Membrane).
- KRT8-18 and epCAM (absent from 7 images) and Ghrelin (4) are left out of the typing on every image: Senior's
  rule test reads an unmeasured negative as unknown, so no immune rule could fire on those seven donors.
- Rule terms on markers this panel never acquired are removed and recorded
  (`typing/immune_rules_on_this_panel.json`); a rule that loses a required positive is unavailable. Available:
  Imm-T-cell, Imm-macrophage; T-CD4-helper, T-CD8, Mac-M2-like (LGALS3 is Senior's M2Gal3). Unavailable: NK,
  plasma, dendritic, neutrophil; regulatory, exhausted, cytotoxic, progenitor-memory T; M1-like macrophage.
- **B cell (owner, 2026-10-17: "B cell should be" a cell type).** Senior's Imm-B-cell needs CD20 and CD79a, and
  CD79a is not on this panel. CD20 alone at scaled 0.10 called 1,392-34,773 per image (CD20's background stays
  11-fold apart between images after UniFORM). The rule asked here (`codex_typing.SUBSTITUTED_RULES`) takes CD45
  as the second required positive and requires CD3e, CD68, E-cadherin and panCK absent: 2-795 B cells per image,
  the T1D donors highest (T1D_6534 795, T1D_6533 716, T1D_6520 393), against 0-1,182 in the earlier typing
  (1,182, 651, 357). Reading CD20 over each image's epithelial background instead left 0-10 on 11 of the 15.
- **Mutually exclusive hormones and CHGA+ hormone-negative cells (owner, 2026-10-18).** "These hormones are
  mutually exclusive so there needs to be resolution or exclusion based on the approved islet cell typing methods.
  The 'only' and 'at all' designations are not biology, they are indicators of technical problems that need to be
  solved. We need to add a CHGA+INS-GCG-SST- cell type." Measured first: at the hormone level, 18-38 % of the
  typed endocrine cells in the objects cleared it for two or more hormones (37 % on ND_6548 even at the unlowered
  tissue p99, because the level sits about 1.3x tissue background while islet hormone runs 15-86x it). The second
  hormone is a touching neighbour's: a cell's runner-up / top ratio is 4-7x larger when at least half its
  neighbours within 15 um have that runner-up as their top hormone (ND_6548 0.21 against 0.05), which is the boundary bleed
  Senior's `typing/validate.py` records for tightly packed, cytoplasm-poor islet cells. REDSEA, run with Senior's
  parameters (`codex_typing.REDSEA`, identical to Senior's `config.py`), removes part of it: on ND_6548 the share of
  typed endocrine cells with two hormones above their tissue p99 falls from 65 % (QuPath's uncorrected cytoplasm
  means) to 41 %, the median runner-up / top ratio from 0.13 to 0.07, and the share under a 2x margin from 15 % to
  11 % (Aab_6505 60 -> 35 %, 16 -> 12 %; T1D_6533 45 -> 30 %, 31 -> 23 %); UniFORM leaves all three unchanged.
  Senior's RESTORE intensity step, which tested exactly this mutual-exclusivity premise, is retired in Senior
  (2026-08-14, never-revive) and is not ported. No compartment removes it
  (raw nucleus, cytoplasm, cell and membrane means leave 13-37 % under a 2x margin; the REDSEA-corrected
  cytoplasm used here 11 %). So the approved method's resolution -- the strongest hormone names the cell -- now
  requires that hormone to be at least `HORMONE_MARGIN` = 2x the runner-up (Senior's confirmation margin,
  `islet_identity.REFERENCE_MARGIN`); under it the cell is Excluded (`hormone_unresolved_margin`), Singles
  included. That excluded 11,095 cells in the objects (4,288 Alpha, 3,546 Beta, 3,182 Delta, 79 ties; 5 % in 1-3
  cell objects to 12 % in islets of 30-100 cells) and moves the Beta share of endocrine cells by at most 2 points
  per status. Senior's Endo-hormone-negative is added as **CHGA+ hormone-negative**: inside an Islet or Cluster,
  CHGA at or above the image's tissue CHGA level and all three hormones measured and below the hormone level,
  called ahead of the immune rules as Senior's islet cascade orders it (hormone, then CD99, then the rest): 2,350
  cells, 6-433 per image (2,235 formerly Excluded, 115 formerly Immune). It counts as an endocrine cell. The
  composition table's per-hormone positivity fields (`Ins_any`, `Ins_single`, `Multi_Pos`, ...) are kept only as
  that technical indicator; the app does not show them.
- **The hormone level (2026-10-17, for the Aab_6558 open item).** Senior admits a cell inside an object when
  max(scaled INS, GCG, SST) clears the tissue p99 of that max. Here each hormone is first divided by its own
  ordinary-tissue p99 in the image (so the cohort scale no longer pits INS against GCG's background), and inside
  the objects the level is lowered, never raised, to the lowest level from which every bin up to the p99 has a
  local false-discovery rate of at most 10 %, background being the tissue distribution scaled to the object cells
  below the tissue median (`classify.object_hormone_level`, `codex_typing.HORMONE_LEVEL_DEPARTURE`). The tissue
  p99 holds false positives among exocrine cells to 1 %, which inside the objects is a global FDR of 0.1-0.2 % on
  every image, and on a weak stain it cuts through the endocrine population. Aab_6558 is that case: its INS stain
  is weak (in the objects INS reaches 2.4x its tissue p99, against 15-86x on the other ND and Aab+ images; the
  donor is GADA+ only with C-peptide 8 ng/mL), a distinct INS-high / GCG-low cluster sits below the p99, and 35 %
  of its object cells were Excluded, 1,284 of them with INS their top hormone. The per-hormone division alone
  changed little (6,876 -> 6,888 endocrine cells); the local-FDR level lowers its level from 1.71 to 0.83 and
  admits 1,474 more (Beta 722 -> 1,769; Excluded 35 % -> 21 %). Well-stained images gain 0.2-3 %; the weak-stain
  images T1D_6533 19 %, T1D_6534 9 %, T1D_6563, T1D_6520 and ND_6232 5-7 %. The subtype winner stays Senior's
  (highest scaled hormone): ranking the hormones by their folds instead would change 0.4-7 % of calls and none
  from Alpha to Beta on Aab_6558, whose alpha-rich make-up (5,899 Alpha, 1,769 Beta) is what its stains show.
- No reference donor is pinned for UniFORM (Senior pinned its own); no hormone is masked (step 1b makes no
  fragments); Senior's unshifted-CD8 refinement exception is not ported.

Owner rulings later on 2026-10-16, applied:

- **Cell quality**: only cells without a nucleus (and the technical failures that leave no measurement: empty
  raster, invalid centroid or area, duplicates) are excluded; size and solidity are recorded in `qc_flags` and
  exclude nothing. 86-94 % of cells are kept. Dropping the nucleus-less cells does not remove most small cells:
  of 302,619 cells under 20 um2, 23,602 have no nucleus, so 278,794 small cells (9 % of kept cells) stay.
- **CD45 rule**: `Imm-CD45-leukocyte` (CD45 >= 0.10, E-cadherin and panCK below) is asked after every Senior
  definition and before the supporting tier, so it never takes a cell from Imm-T-cell or Imm-macrophage; it names
  574-10,762 leukocytes per image that the panel cannot name further (most in T1D_6533 and T1D_6534).
- **Structure objects of the five additions** (`typing/structure_thresholds.py`): the original classifiers are
  linear in smoothed intensity, and a plain threshold on the Gaussian-smoothed marker (sigma 2 px) reproduces the
  original 15's objects at median Dice 0.88 (Nerve, B3TUBB), 0.89 (Capillary, CD31) and 0.73 (Lymphatic, PDPN; the
  original was a random forest). Those thresholds agree across the 15 in RAW units (p90/p10 1.4, 1.4, 2.0) and
  not after the UniFORM factors (4.6, 2.8, 14): the classifiers were one raw cut on images acquired alike.
  **Revised 2026-10-17.** Each addition takes the 15's median raw threshold unscaled (Nerve 475, Capillary
  1,951): predicting each original from the other 14, the raw median is the best transfer (|log error| median
  0.07 and 0.17) and the UniFORM ratio used on 2026-10-16 among the worst (0.32, 0.27; UniFORM's PDPN factors
  run 0.64-15.9 across images whose raw thresholds agree within 2x). The objects are drawn as
  `4_Islet_Cells.groovy` drew them (inside Tissue, split, 2 / 4 um2 minimum, no object touching a 10 um edge
  ring) and match the original classifiers' own objects on the additions (Dice 0.81-0.90).
  **Lymphatic keeps the original classifier's objects on the additions** (`KEEP_CLASSIFIER`). That classifier is
  a random forest on PDPN texture, not a threshold, and no transferred threshold reproduces it on both cohorts:
  the raw median drew lymphatics on 0.06-0.15 % of tissue in four additions, below every original (0.15-3.22 %,
  median 1.10 %), because their PDPN bright tail is dimmer; the UniFORM ratio drew T1D_6520 at 2.6x the
  classifier (112,000 against 43,000 um2, the open item); rules that fit the additions fail on the 15. The
  classifier's objects cover 0.39-0.54 % on all five, and T1D_6520's passed the 2026-09-21 audit. `make` still
  draws and records the threshold alternative. The classifier objects are in `$OUT/structures/superseded_additions/`;
  records in `structures/thresholds_fit.json` and `thresholded.json`.
- **Background vs variation between images** (`typing/background_check.py`, `background_check.csv`): the
  UniFORM bar reads each image's histogram mode, which for markers most cells carry is the positive population.
  Measured instead in cells of a lineage that does not carry the marker (epithelial cells for hormone, immune,
  vascular, stromal and neural markers; immune cells for epithelial markers), the background is aligned within
  1.2-1.5x for INS, CD8a, CD35, SMA and Ki67, and still differs 5-14x between images for CD3e, CD45, CD68, CD163,
  HLADR, PDPN, B3TUBB and GAP43. It stays below the typing cut: in those negative cells CD3e, CD8a, CD163, B3TUBB,
  PDPN, CD31 and the hormones reach scaled 0.10 in under 1.5 % of cells on every image, CD68 and CD45 in up to
  2-3 %, CD4 in up to 3.9 % (LGALS3 is excluded from this reading: ductal epithelium carries galectin-3).
- **Endocrine checks** (`typing/endocrine_checks.py`): of the cells inside the endocrine objects, 64-85 % are
  called endocrine (2026-10-17, after the hormone level: 75-86 %). The Excluded cells there (exocrine fall-through
  inside an object) split into CHGA-positive (other endocrine; 4-405 per image), a hormone statistic between the
  tissue p95 and the level, and no endocrine evidence. Before 2026-10-17 the middle band was "any hormone above its
  floor", which two-thirds of ordinary exocrine cells pass; now `hormone_below_level_expected` gives the count
  background alone puts there, and on every image the observed count is within 2x of it (Aab_6558 30 against 91):
  what stays Excluded carries no more hormone than exocrine tissue does. Islets almost always hold an endocrine
  call (2 of 1,491 do not); clusters without one are 2-22 per image (up to 75 before) and singles 4-48.

Run of 2026-10-17 after these rulings and the two open items: `$OUT/typing/summary.csv`, written by
`typing/summary.py` (per image: kept cells, Beta / Alpha / Delta,
Excluded and its three kinds, Immune by definition, structure calls, Exocrine, refined Unresolved).

## App root (P4), 2026-10-16

```bash
$PY scripts/uniform/assemble_uniform.py --run $OUT --donor-source data/development/protein-expanded-20261015-v3 \
    --output data/development/uniform-20donor-<date>-v<N> --report-root tmp/uniform-20donor-<date>-v<N>
```

Every Islet and Cluster (Islet_N, Cluster_N) is a row with `object_class`. Its core is the cells whose
centroid it contains. Its peri band is the 20 um ring inside Tissue, minus every Islet and Cluster. Singles
(Single_N) are not rows (owner, 2026-10-18: "Singles should not be included in any islet-level analysis except
for quantifying or showing how many there are"): they have no core or band, signed distances ignore them, and
their typed cell is an ordinary tissue cell, so it counts in a neighbouring object's band when it lies there.
They are counted in `census.csv` and `endocrine_objects_by_donor.csv` (per donor: tissue area, Islets, Clusters,
Singles, Singles by cell type, each class per mm2 of tissue) and drawn in `gson/` as `EndocrineSingle`. Cells are
shown under the app's labels: Beta / Alpha / Delta cell; Endocrine; CD8a Tcell and CD4 Tcell (Senior's
T-CD8 / T-CD4-helper); T cell; B cell; Macrophage; Immune (CD45-only or supporting tier); and the new
Exocrine (the residual). Since 2026-10-17 there are no Blood Vessel, Neural or Lymphatic cells. Endocrine, T cell and Immune hold the cells left without a subtype; the app
does not offer them beside their subtypes (`PHENOTYPE_RETIRED`, owner 2026-10-17), though their cells stay in every
denominator. Excluded cells and cells without a nucleus are counted in `census.csv` and
kept out of every table. Hormone positivity is the image's hormone level on each scaled hormone (D4; the level
times that hormone's own tissue p99, `levels.hormone_scaled` in the calls audit), and the
regional and per-cell `_pos` flags are scaled >= 0.10.

Two things are carried from v3 because they were never computed from cells: donor demographics and
autoantibodies, and the viewer image links. The served UMAP and Leiden clusters are the endocrine composition of the objects with ten or more core cells (the Beta, Alpha, Delta and CHGA+ hormone-negative shares of each core's endocrine cells; `scripts/uniform/composition_embedding.py`): one
30-neighbour graph for both, no donor correction, and no point or cluster for smaller objects (owner, 2026-10-19,
from v6). The Harmony protein PCA, UMAP and Leiden (`build_protein_embedding.fit_embedding`) are still fitted and
kept in the H5AD for the record (`X_protein_pca`, `X_pca`, `X_umap_protein`, `protein_leiden_<r>`). Seven objects
made only of cells whose nucleus fills the cell have no cytoplasmic mean; X keeps the NaN, and the protein
embedding is fitted with the marker's cohort mean in its place (`embedding_imputation.json`). The composition
diagnostics are `composition_embedding.json` in the report root.

Why composition (measured on v5, 2026-10-19). Disease is a donor property, so only neighbours from *other* donors
test whether an embedding places objects by disease. Among each object's 30 nearest neighbours from other donors,
same disease exceeded chance by +0.098 on composition (17 % same-donor neighbours, 6.5 % by chance), +0.081 on the
v5 protein + Harmony layout, +0.068 on protein PCA without Harmony (45 % same donor) and +0.025 on protein + Harmony
refitted on the same objects. Harmony removes donor shifts and, with them, much of the disease difference; the
original 15-donor UMAP's apparent disease separation was donor separation (72 % same-donor neighbours; other-donor
same disease at chance). Trained on the other donors, endocrine composition names a held-out donor's disease for
14 of 20 donors (T1D 6 of 6, ND 6 of 7, Aab+ 2 of 7), the 31 markers for 10 of 20. Objects under ten core cells
take only a few composition values and carried no more disease signal (+0.13 at two cells, as discrete blobs).
Composition Leiden is less seed-stable than the protein Leiden (adjusted Rand 0.53-0.78 against 0.81-0.99): its
clusters cut a continuum of beta-cell share.

From v7 the UMAP uses min_dist 3, spread 3 (owner, 2026-10-19), and the fit runs on one thread. At v6's 0.5 / 2,
80 % of points lay within one 4-pixel marker of a neighbour in the app's plot (42 % at 3 / 3), and the 31 all-beta
and 21 all-alpha objects, which have identical features, drew as islands detached from the cloud; they are real
(small beta-only islets in ND and Aab+ donors, 31 against 0.2 expected from each donor's pooled make-up; alpha-only,
insulin-deficient islets, 12 of the 21 in T1D donors, against 3.8 expected) but sit at the ends of the composition
continuum, not apart from it. The neighbour graph and Leiden do not use min_dist or spread. With the BLAS and numba
thread counts left to the machine, the same input gave a different layout and Leiden adjusted Rand 0.63 against a
one-thread run, so v6's layout cannot be reproduced exactly; v7's can.

The app reads `object_class` through an app-wide "Endocrine objects" selector (Islet / Cluster). The
Trajectory tab does not use it: it fits every object with ten or more core cells (`total_cells_core` in
`barlow_islet_features.csv`), whatever its class.

## Canary, 2026-10-15

Aab_6505 and ND_6232 went through every step. Cells matched to the earlier typing show what the first fit did with
beta cells:

| Image | Fit | INS threshold | Earlier beta cells inside the fitted islets | alpha | delta |
|---|---|---|---|---|---|
| Aab_6505 | first run (raw highest marker, PGP9.5 Neural) | not used | 8.0 % of 3,071 | 81.7 % | 2.6 % |
| Aab_6505 | corrected | 177 | 87.3 % | 83.4 % | 99.4 % |
| ND_6232 | first run | 1,411 | 52.1 % of 1,171 | 99.6 % | 76.2 % |
| ND_6232 | corrected | 129 | 98.0 % | 99.7 % | 100 % |

In the first run PGP9.5 was the highest marker of 67 % (6505) and 77 % (6232) of beta cells; as a Neural marker it
barred them from the CHGA route, and with CHGA below the image's level they counted as `other`, so the fit learned
to leave insulin out. Percentages are of the earlier labelled cells matched to a new cell within 3 µm, against the
fit's own emulation of script 1; the QuPath objects are what the owner reviews. The "corrected" rows use each
image's own tops, before the floor; the run below gives the canaries' QuPath objects.
- No pLN annotations exist in these images; exclusion annotations, if a source image has any, become `Artifact`.

## Islet run on all 20 images, 2026-10-15

Every image was segmented again (no failures), then fitted and run through the islet chain with the insulin and
somatostatin floor (cohort medians 883 and 2,005). `summarize_review.py` writes the per-image table to
`$OUT/review/summary.csv` (objects, thresholds, floors, and where the earlier beta / alpha / delta cells lie: inside
script 1's outlines and inside the final objects).

- **19 of 20 images have islet objects.** Final objects hold 75–98 % of the earlier beta cells and 66–92 % of the
  alpha cells on 12 images. Canaries: Aab_6505 63 / 63 / 86 % (beta / alpha / delta; INS threshold 232), ND_6232
  94 / 87 / 95 % (INS 129).
- **T1D_6563 has none: Senior's guard refused the fit.** Only 45 cells have a hormone highest and CHGA at its
  level (minimum 50). Its alpha cells carry clear glucagon (median about 4,000) but little CHGA (median 132), while
  acinar cells and macrophages far from any islet set the level at 335; leaving out nerve-like cells moves it to
  338. With the minimum lowered to 40 the fit uses glucagon and CHGA only and covers 51 % of the earlier alpha
  cells. Awaiting the owner's decision (hand-set thresholds for this image, or no islets for now).
- **Step 1b drops edge cells.** Script 1's outlines enclose 80–100 % of the earlier beta and alpha cells on most
  images; 1b then keeps only the cells that fit whole. The cells it drops lie mostly inside (median 84–87 % of their
  area) but reach a median 1.7–2.0 µm beyond the outline (Aab_6558, Aab_6505, T1D_6520, Aab_6450). Senior's
  physical tolerance (1.02 µm instead of one CODEX classifier pixel, 0.75 µm) would keep 15–17 % of them, so the
  loss is the rule rather than the port. It is largest where islets break into many small objects: beta cells
  inside script 1 against the final objects, Aab_6558 80 → 37 %, T1D_6520 82 → 57 %, Aab_6505 83 → 63 %,
  ND_6539 89 → 71 %; Aab_6450 100 → 97 %.
- **Weaker fits:** ND_6548 (beta 54 %, delta 52 %; CHGA threshold 1,367, somatostatin unused), T1D_6550 and
  T1D_6551 (alpha 58 % and 52 %; T1D_6551 uses neither insulin nor somatostatin).

## Run

```bash
PY=/home/smith6jt/miniconda3/envs/scvi-env/bin/python
OUT=data/development/uniform-20donor-20261015
"$PY" scripts/uniform/make_project.py --out "$OUT"                  # once: 20 fresh entries
"$PY" -u scripts/uniform/run_uniform.py --out "$OUT" --images all       # or a list of images; a fit needs every image's cells
```

`run_uniform.py` runs, per image: `tissue_export` (read-only, on the source project) → `tissue_import` →
`preflight` → `0` InstanSeg (GPU) → `2` measurements → `C_cells` (Parquet) → `C_planes` → `C_geojson`; then, once
every cohort image has its cells, the cohort scale (`cohort_scale.py`, built by the first `fit`) and per image `fit`
→ `1` → `1_outlines` → `1b` → `4` → `3` → `8` → final exports. A run works in those two phases, so `--images all`
finishes every image's cells before the first fit. Each QuPath step runs on a private copy of the project file
with Senior's footer, which saves only when the step's check passed and writes `$OUT/provenance/<image>/<step>.json`.
A step whose prerequisites have not passed is refused; `--force` reruns a passed step. Read-only source steps fail
if any file of the source entry changes.

Reruns are exact, not incremental:

- Running a step moves the records of every step that rests on it to `provenance/<image>/history/`, so those steps
  run again.
- Before its first islet run an image's QuPath data is copied to `snapshots/<image>/before_islets.qpdata`, and
  every run of step `1` starts from that copy. Running a cell step (`tissue_import`, `preflight`, `0`, `2`) retires
  the copy.
- A fit is stale when the cohort scale changed or no longer holds every image's current cells.

### Starting an image again

Cell detection keeps existing annotations, so an image whose islet steps ran cannot go back through the cell steps
in place; `run_uniform.py` refuses. Move the image's outputs aside (nothing is deleted) and run it from the start:

```bash
"$PY" scripts/uniform/reset_entry.py --out "$OUT" --label <name> --images <image> [...]
```

The first canary run is kept that way in `$OUT/superseded/canary-v1-20261015/`. Segmenting both canaries again
from fresh entries gave identical cell tables (148,003 and 93,727 cells; every centroid, area and mean intensity
equal), so InstanSeg 0.1.1 through QuPath on this GPU reruns bit for bit.
