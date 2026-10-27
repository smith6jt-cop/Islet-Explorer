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
| Tissue outlines | The original project's Tissue for the original 15, reused unchanged; for the five additions, redrawn from the image (`redraw_tissue.py`, every hole filled) after the review project's outlines were found defective (below) | owner, 2026-10-15; additions 2026-10-27 |
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

## Additions' Tissue redrawn (owner, 2026-10-27)

The additions' Tissue outlines, reused from the six-donor review project's `Tissue` pixel classifier, were never
compared with the images until the owner found ND_6232's "horribly inaccurate" (2026-10-26; standing rule in the root
CLAUDE.md). Measured against DAPI (`tissue_qc.py`): ND_6232's outline held 39 % of the nucleated tissue (the original
15: 93.9-98.9 %), half the other donors' cells per mm2, and 89 % of its endocrine objects touched it; ND_6539, Aab_6558,
Aab_6090 and T1D_6520 held 289-3,525 holes under 500 um2 each (the 15: none). Cell detection, islet objects and
structures run only inside Tissue.

- `redraw_tissue.py` draws Tissue from the image: DAPI (nuclei) or ColIV / VIM / SMA (stroma) at the original
  classifier's 3.02 um, pieces under 80 um2 removed, clipped to the image less 5 um as the original Tissue was. Without
  hole filling the same rule reproduces the original 15 outlines at Dice 0.958-0.991. The owner approved the outlines
  "if all holes are filled", so every hole is filled. The files are `tissue_outlines/<image>.geojson` (tracked despite
  the `*.geojson` ignore rule), named in `cohort.csv`'s `tissue_outline`; `run_uniform.py`'s `tissue_export` copies
  such a file instead of exporting from the source project and records both checksums. The annotation must be unnamed:
  a named Tissue changes every cell's `Parent` from "Annotation (Tissue)", which the pipeline keys tissue cells on, and
  the first import (named) left every cell ineligible; the runner now refuses a named outline and `tissue_qc.py`
  reports the proportion of cells with that Parent.
- The additions' Lymphatic objects (KEEP_CLASSIFIER) were the original Lymphatics classifier's objects inside the old
  outlines; `additions_lymphatics.py` (with `qupath/lymphatics_in_tissue.groovy`, read-only) draws them again with the
  same classifier file inside the new Tissue, before `structure_thresholds.py make`.
- Run `data/development/uniform-20donor-20261027` is a full copy of `uniform-20donor-20261015` (which v10 came from and
  is untouched), with the five images reset (`superseded/review-tissue-outlines-20261027/`, and
  `superseded/named-tissue-annotation-20261027/` for the named-outline pass) and run again from `tissue_export`. Their
  new cells changed the cohort scale's insulin and somatostatin floors (883 -> 817, 2,005 -> 1,916), so all 20 fits and
  islet chains ran again (the 15 from their pre-islet snapshots), then typing (`redsea_cells.py --force` for the five)
  and `assemble_uniform.py` into `data/development/uniform-20donor-20261027-v11`.
- After the rerun the five hold 99.6-99.8 % of their nucleated tissue, no holes, 11,345-12,607 cells per mm2 of
  nucleated tissue (the 15: 10,676-14,171) and 1.1-5.9 % of objects touching the outline. ND_6232: 93,727 -> 172,332
  cells. The 15's cell counts are unchanged; their typed counts move by 0-40 endocrine and up to about 120 immune cells
  per image (the cohort scale and UniFORM alignment include the five), and four of them changed by 1-6 objects.
  Checks and renders: `output/tissue_outline_check_20261026/`.

## Whole-cell fallback and CD8a/CD4 double positives (owner, 2026-10-27; root v12)

- **No cell is without cytoplasm** (owner: "there should be no such thing as 'no cytoplasm'. This was already
  determined."). Where the nucleus fills the cell REDSEA leaves the Cytoplasm compartment NaN (345,328
  analysis-eligible cells, about 10 %); every cytoplasm-measured marker of such a cell is now read over the whole cell
  (REDSEA's `Cell__<marker>`), in `typing/align_scale.py` `read_measured`, so alignment, scale, typing and the app's
  marker values all see it. Before, a missing panCK left Senior's T-cell rule (which forbids panCK) undecided: CD3e+
  cells fell to the CD45 fallback or to Exocrine, and islet cells had no hormone. Only the cytoplasm markers' scale
  moved (floor 0.80-1.15x, top 0.94-1.07x); GCG's alignment now ends just over the bar (mode SD 0.111 against 0.10).
- **CD8a and CD4 both positive**: Senior's T-CD8 forbids CD4 and T-CD4-helper forbids CD8, so such a T cell kept no
  subtype ("T cell"). It now takes the subtype of the higher scaled marker (`codex_typing.T_SUBTYPE_BY_HIGHER`,
  `immune_subtype.py`, `immune_subtype_resolved_by`): 2,453 cells; 11 exact ties stay "T cell".
- Typing reran from `align_scale.py` on the run (v11's typing kept in `superseded/typing-v11-20261027/`) and
  `assemble_uniform.py` built `data/development/uniform-20donor-20261027-v12` (its `git_commit` field reads c3ac421;
  these changes were in the working tree and are committed in the next commit). v11 -> v12: analyzed cells
  3,162,613 -> 3,172,248; 34,325 cells change label: Exocrine -> Macrophage 9,898, -> Immune 7,901, -> CD8a Tcell 1,093;
  formerly excluded -> Alpha 5,130, Beta 2,821, Delta 1,519; T cell -> CD8a Tcell 1,566, -> CD4 Tcell 705. Totals:
  T cell 5,685 -> 3,795, CD8a Tcell 13,944 -> 16,943, CD4 Tcell 4,260 -> 5,391, Macrophage 61,140 -> 72,041.
  Endocrine objects are unchanged (4,853 outlines); embedding objects 1,207 -> 1,278.
- **No 0.10 level below a broad marker** (owner, the same day: "Once there is positivity for a broad defining marker,
  we do not need the 0.10 positivity level"). A T cell (CD3e+) that neither T-CD8 nor T-CD4-helper names now takes
  the higher scaled of CD8a and CD4 at any level; exact ties (both at the scale floor) keep "T cell". 5,763 cells are
  named this way, 485 tie. T-CD8 16,943 -> 17,891, T-CD4-helper 5,391 -> 7,753. Only `immune_subtype.py` and
  `summary.py` reran (v12's `immune/` kept in `superseded/typing-v12-immune-20261027/`); root
  `data/development/uniform-20donor-20261027-v13`. Mac-M2-like, the one other subtype on this panel, has no
  competing subtype (Mac-M1-like is unavailable) and still requires CD163 and LGALS3 at 0.10; it is not an app label.
- **No exact ties** (owner: "I don't believe there is biological support for exact ties, for these a new review of
  the total marker profile is indicated"). The v13 ties (485) were made by the scale's 0-1 clip (474 with CD8a and CD4
  both below the floor, 11 both at the top). The CD8a/CD4 comparison now uses scaled values before the clip; a T cell
  tied after the clip stays a T cell only when its highest marker (lineage markers and CD45, before the clip) is
  CD3e, CD8a, CD4 or CD45 (`codex_typing.T_TIE_KEEP_MARKERS`), otherwise `classify.py` types it again without the
  T-cell rules (`resolved_by` 'tie_profile_*'; `immune_subtype.py` leaves Imm-T-cell out for it): 486 ties, 423 kept
  (310 -> CD8a Tcell, 112 -> CD4 Tcell, 1 kept another call), 63 retyped (39 Immune, 23 Exocrine, 1 Excluded). No
  "T cell" remains. Review table and profile heatmap: `output/tcell_tie_review_20261027/`.
- **Macrophages by marker profile** (owner: "We don't want to get into M1 and M2-like labels, but we can indicate
  macrophages and LGALS3 as defined in the senior workflow. We can also define CD4+ macrophages", then "Use
  CD68+CD163+Macrophage_LGALS3+CD4+ with each plus also able to be a minus", and a macrophage is CD68 or CD163). The
  macrophage rule on this panel is `SUBSTITUTED_RULES['Imm-macrophage']`: CD68 or CD163 >= 0.10, CD3e, ECAD and panCK
  absent (Senior's needs CD68 or Iba1, and Iba1 is not on this panel). No M1/M2 subtypes are resolved
  (`IMMUNE_LEVEL3`); each macrophage takes `macrophage_profile` (`codex_typing.macrophage_label`, each marker + at
  >= 0.10), which `assemble_uniform.py` uses as its label: 12 labels (CD68-CD163- cannot occur). In text the + and -
  are superscripts. Macrophages 72,041 -> 109,589 (25,212 from Exocrine, 12,192 from Immune, 144 formerly Excluded,
  all CD68-CD163+); `cd8_to_macro_ratio` counts all twelve. Classify onward reran (v13's calls, refine and immune in
  `superseded/typing-v13-20261027/`); root `data/development/uniform-20donor-20261027-v14`.
- **No default class** (owner: "I don't want any more misclassified cells falling through the cracks. The established
  cell typing should classify cells appropriately and there should be no default of any class such as exocrine or
  immune. If typing has 'left over' cells and a check of top markers can't place them, call them by their top marker
  e.g. CD45+"; top marker at "Any level"; "ColIV cannot be a cell type, if it is highest, it is a reason for
  exclusion"). Senior's supporting tier (Immune when CD4, CD8a, CD163, HLADR or LGALS3 >= 0.8715 and no rule fired;
  never put to the owner) and the Exocrine fallthrough are gone. Exocrine is a rule: E-cadherin or panCK >= 0.10. A
  cell no rule claims is placed by its highest marker (`codex_typing.TOP_MARKER_PROFILE`: lineage markers and CD45,
  scaled before the clip): CD3e T cell, CD20 B cell, CD68 or CD163 macrophage, E-cadherin or panCK Exocrine; ColIV
  excludes it (`top_marker_ColIV`); any other marker labels it `<marker>+` (`TOP_MARKER_LABELS`, 14 labels).
- **Zero floor for markers most cells express** (owner: "Floor at zero"). The scale's floor is the most common
  positive value; for E-cadherin, KRT8-18, beta-catenin, beta-actin and EpCAM, whose background REDSEA already
  compensates to 0 (17-28 % of cells), that is the expressing population: E-cadherin's floor 34,035 sat above 85 % of
  cells (median 9,044), so only 18-27 % of the former Exocrine cells reached 0.10 in any image. Their floor is now 0
  (`ZERO_FLOOR_MARKERS`); `scale.parquet` keeps the former floor as `mode_floor`. Owner, the same day: the change
  applies only to cells the rules leave ("I thought we were addressing cells that were not already typed"): the
  immune rules, their subtype step and the T-cell tie review read `classify.rule_scale` (the mode floor), so all
  318,869 cells v14's rules typed keep their v14 call, immune class and subtype (checked cell by cell). A first run
  with the zero floor everywhere moved ~26,000 rule-typed immune cells to Exocrine and was discarded.
- The 2,878,457 cells v14 left to its defaults (Exocrine fallthrough, supporting tier, retyped ties) became:
  Exocrine 1,876,363 (1,603,778 by the rule, 272,585 by top marker), Excluded 322,174 (319,287 ColIV), B cell by
  CD20 62,090, macrophage 21,655, T cell 5,762, and by top marker LGALS3+ 143,349, CD34+ 135,679, VIM+ 92,949, CD31+
  53,670, B3TUBB+ 30,452, CD4+ 28,057, SMA+ 23,698, PDPN+ 21,239, CD45+ 16,649, INS+ 14,917, GCG+ 11,983, CD8a+ 9,818,
  SST+ 7,915. Typing reran from `align_scale.py` (v14's typing in `superseded/typing-v14-20261027/`); root
  `data/development/uniform-20donor-20261027-v15` (2,862,676 cells analyzed, 334,650 excluded).

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

v8 (2026-10-20, test version only; the live app stays on v7) differs from v7 only in `cell_distance.csv`: its
`<marker>_pos` flags cover all 31 panel markers instead of 18 (owner: the distance profile "needs to have all markers
that are common to all donors"; every marker is measured in every image, and `build_image` stops if one is not).
The flag is unchanged: the marker at or above 0.10 on the cohort scale.

v9 (2026-10-20, test version only) adds `signed_dist_islet` to `cell_distance.csv`: each cell's signed distance to the
nearest Islet only (negative inside one), beside `signed_dist`, the distance to the nearest Islet or Cluster. Clusters
are small and scattered through the tissue, so 90 % of the cells beyond 20 um of any object were measured to a
Cluster (median tissue-cell distance 133 um to any object, 204 um to an Islet), which blurred the profile around
Islets (T1D CD3e+ at 0-5 um: 5.0 % by any object, 6.1 % by Islet). The distance-to-islet profile lets the reviewer
choose. Nothing else changes from v8.

v10 (2026-10-22, test version only) takes every donor field from the nPOD export (owner: "All data can be found in
All-Download_From_nPOD_2026-07-27_163844.xlsx"): `scripts/build_donor_metadata.py --uniform-cohort` writes
`scripts/uniform/donor_metadata.csv`, which `--donor-metadata` serves and from which each donor's age is taken. The one
change to v9's values is donor 112's age, 3.72 years as nPOD records it instead of the 3.00 carried from v3 (age group
unchanged); status, sex and autoantibodies are checked against the donor source and agree for all 20 donors. The
table adds height, weight, BMI percentile, pancreas weight, admission and peak glucose, ICU days, organ transit time,
RIN, race/ethnicity, recovery type (DCD or not), cause of death and ABO group. nPOD's 0.0 HbA1c (6090, 6232) is
written as missing, and 6534's transit time (-1034 min) is left out as impossible.

```bash
$PY scripts/build_donor_metadata.py --uniform-cohort scripts/uniform/cohort.csv \
    --source All-Download_From_nPOD_2026-07-27_163844.xlsx --output scripts/uniform/donor_metadata.csv
$PY scripts/uniform/assemble_uniform.py --run $OUT --donor-source data/development/protein-expanded-20261015-v3 \
    --donor-metadata scripts/uniform/donor_metadata.csv \
    --output data/development/uniform-20donor-20261022-v10 --report-root tmp/uniform-20donor-20261022-v10
```

The app reads `object_class` through an app-wide "Endocrine objects" selector (Islet / Cluster). The
Trajectory tab does not use it: it fits every object with ten or more core cells (`total_cells_core` in
`barlow_islet_features.csv`), whatever its class.

### Production release (data schema 5), 2026-10-22

A uniform root is packaged for PubApps by `scripts/build_release_data.R --schema 5` (`deploy/README.md` § 2):
the five prepared tables as Parquet, the spatial tables as Parquet, the geometry as GeoPackage, `cells/`, the
images with their XML and offsets sidecars resolved beside each linked TIFF, the passthrough files
(`adata_ins_root.h5ad`, `barlow_islet_features.csv`, `islet_spatial_lookup.csv`, `donor_metadata.csv`,
`uniform_release.json`, `endocrine_objects_by_donor.csv`, `protein_phenotype_rules.csv`) and the H5AD's
`composition_embedding` / `protein_embedding` records as `summary/embedding_metadata.json`. `census.csv`,
`cell_types/`, `image_sources.json` and `neighborhood_metrics.csv` are provenance of the root and are not
served. The builder refuses `--schema 1` for a uniform root; the runtime (`ISLET_DATA_SCHEMA=5`) refuses a
release that lacks any of the uniform files. The first candidate is built from v10 at commit `d55680a`
(`tmp/ufrc-release-20261022/`; evidence in `docs/ufrc_readiness.md`).

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
