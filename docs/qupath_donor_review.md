# Review the six newly processed donors in QuPath

The current [six-donor QuPath project](../data/reviews/qupath-endocrine-rollout-20260918/project.qpproj)
contains 1,431 owner-approved HPAP-078 boundaries and 2,030 cyan `Endocrine
pixel candidate` objects across 6090, 6232, 6520, recovered 6539 and 6558.
Open the project in QuPath 0.7, double-click each donor, inspect the candidate
outlines against the raw INS/GCG/SST signal, and save any review edits in that
project. The [offline gallery](../data/reviews/endocrine-rollout-20260918/index.html)
contains whole-field and size-stratified examples. All 943,068 source cells are
retained. The five donors still need the owner's own visual decision; the
HPAP-078 decision does not transfer to them. The earlier HPAP-only project and the
[original15 pixel refinement](endocrine_reference_refinement_20260917.md) and
[existing-classifier comparison](endocrine_pixel_classifier_pilot_20260917.md)
are retained as baseline evidence, along with their editable review projects.

The sections below describe earlier September 17 review projects and are kept
for provenance. Their counts and R5 masks are not the current five-donor review.

**Boundary correction in progress:** the R5 candidates linked below have a
confirmed SAM background-overreach failure and are withdrawn as a usable
segmentation result. The native files and any review edits remain intact. See
the [quality status record](../data/reviews/qupath-regenerated-islets-20260917-r5/quality_status.json).

**Withdrawn R5 candidate, retained as historical evidence:**
[`qupath-regenerated-islets-20260917-r5/project.qpproj`](../data/reviews/qupath-regenerated-islets-20260917-r5/project.qpproj)
for the 1,162 regenerated boundary candidates across these six donors. The
[regeneration record](islet_regeneration_20260917.md) explains the colors,
raw-channel comparisons and verification. That independent project retains all
943,068 original cells. The earlier project described below, including any
review edits, is preserved.

The six additions have complete native QuPath object files. A separate,
editable project is available at:

```text
/data/panc_CODEX/Islet-Explorer/data/reviews/qupath-six-donors-20260917/project.qpproj
```

Open that project with QuPath 0.7 and double-click a donor in the project list.
Its objects load with the image; no segmentation rerun or object import is
needed. The project references the existing full-resolution source images and
contains independent copies of the saved object data. Save review edits in this
project to preserve the processing evidence and app candidate.

| Donor | Status | Cell objects | Islet annotations | Object file within the review project |
|---|---|---:|---:|---|
| HPAP-078 | T1D | 267,993 | 1,234 | `data/42/data.qpdata` |
| 6090 | Aab+ | 170,566 | 868 | `data/39/data.qpdata` |
| 6232 | ND | 92,012 | 128 | `data/43/data.qpdata` |
| 6520-01A | T1D | 102,526 | 321 | `data/40/data.qpdata` |
| 6539, recovered acquisition | ND | 160,837 | 238 | `data/38/data.qpdata` |
| PAN-6558 | Aab+ | 149,134 | 200 | `data/44/data.qpdata` |
| Total | | 943,068 | 2,989 | |

These are the six donors with the newly executed segmentation pipeline. The
original 15 donors in the all-donor analysis reuse their earlier processing.

**Islet review update, 2026-09-17:** the owner's HPAP-078 review found many small
annotations without hormone expression. These are unvalidated candidates, not
accepted islet detections. Raw-channel examples, flags for all 21 donors, SAM
alignment and the required correction workflow are in the
[trajectory revision record](trajectory_revision_20260917.md). The review
project's objects and the owner's edits have not been overwritten.

## What is available to inspect

- Actual cell outlines and available nuclear outlines from InstanSeg, with
  original object UUIDs and parent hierarchy.
- Tissue, SAM-refined islets, existing IsletExpanded annotations, and
  nerve, capillary and lymphatic annotations from the processing pipeline.
- Saved cell measurements, including the 31 shared marker means, and existing
  QuPath classifications such as INS_pos, GCG_pos and Multi_Pos.

The QuPath classifications are upstream marker calls. They are not the
validated retained cell-type assignments needed for the trajectory denominator.
The 201 wholly unmeasured cells remain present with their missing measurements.
Some cells have no separate nuclear ROI; the native verification records the
number with nuclear outlines for each donor.

The existing IsletExpanded objects are upstream QuPath annotations. The
all-donor analysis additionally constructs tissue-clipped peri-islet bands
excluding neighboring cores; those derived bands are not added to this native
review copy. Their exact overlays are now exported separately under
`data/reviews/analysis-regions-20260917/`, with every compartment area checked
against the count artifact. Import only the matching donor file into a separate
review entry if needed; importing twice duplicates the overlay annotations.

The project uses the corrected 6539 acquisition from
`expanded-20260911/recovery/sources/`. The earlier file named 6539 in
`Images/New/` contains duplicated 6090 pixels and must not receive these
objects. Recovered 6539 still needs its own exclusion and boundary review;
the earlier exclusion annotation belongs to the other acquisition.

## Import into an image you have already opened

QuPath 0.7 can import objects directly from a `.qpdata` file. On a fresh project
entry containing the matching full-resolution image, run
`importObjectsFromFile(...)` in the QuPath script editor. For example, for
**6232 only**:

```groovy
importObjectsFromFile('/data/panc_CODEX/Islet-Explorer/data/reviews/qupath-six-donors-20260917/data/43/data.qpdata')
```

Choose the donor's file from the table. Use a fresh entry to avoid duplicating
existing objects. Object coordinates are full-resolution image pixels, so the
image must have the same orientation, origin and dimensions. The exact source
and corresponding app TIFF paths are recorded in the review project's
`donors.csv`. If working on another machine, make those images accessible and
update the image locations in the copied project.

The native project is the most direct way to retain the saved hierarchy,
channel setup and measurement calibration. The analysis `cells.tsv` files
contain measurements rather than complete boundaries, and
`annotations.geojson.gz` contains annotation geometry rather than all cell
segmentations; neither alone replaces the native cell objects.

## Provenance and verification

`review_manifest.json` maps each donor to the correct processing project,
source image, export manifest and copied `.qpdata` SHA-256. The project contains
only the four completed entries from the first project and the two completed
entries from the recovery project, avoiding the stale 6090/6539 entries in the
first project. Object files occupy approximately 3.54 GiB; source image pixels
were not duplicated.

Verification lives in the review project's `verification/` folder:

- `native-summary.json`: QuPath reads every saved hierarchy, checks cell/islet
  counts, counts nuclear ROIs and checks completeness of all 31 marker means.
- `analysis-alignment.json`: matches every cell UUID to the analysis inputs
  using the documented `case_id:qupath_uuid` analysis ID format, reconciles
  measured/unmeasured cell counts and checks that reading preserved file hashes.
- `CASE-native.json` and `CASE-cell-ids.txt`: per-donor evidence.

These are native file and analysis-alignment checks, not a completed visual
or biological review. Review edits should be recorded separately from the
unchanged starting manifest.

The project was prepared using
`scripts/cohort_expansion/prepare_review_project.py`, which refuses to
overwrite an existing review folder. Native reading is reproducible with:

```bash
env JAVA_TOOL_OPTIONS='-Xmx24g -Djava.awt.headless=true' \
  /home/smith6jt/QuPath/bin/QuPath script \
  --args /data/panc_CODEX/Islet-Explorer/data/reviews/qupath-six-donors-20260917 \
  scripts/cohort_expansion/verify_review_project.groovy
```
# Original15 reference-trained pixel diagnostic (2026-09-17)

A separate [HPAP-078 diagnostic project](../data/reviews/qupath-endocrine-pixel-diagnostic-20260917/project.qpproj)
contains 1,056 conservative endocrine pixel candidates and all 267,993 original
cells. It uses a new INS/GCG/SST forest trained from signal-supported reference
regions in the original 15 donors. This is not an accepted replacement delivery;
weak/background staining and boundary extent remain unresolved. See the
[method, checks and plots](endocrine_reference_refinement_20260917.md).
