# Owner-approved HPAP boundaries and five-donor rollout

On 18 September 2026 the owner confirmed the `Endocrine pixel candidate`
objects in the HPAP-078 QuPath project as visually good and instructed us to
proceed. The approval is bound to the saved native object file, SHA-256
`d342f3105a5727c3dde42bb4517e4ecbe5e0fa652d37b8f2fd1f778b94f6af75`.
The reviewed project remains unchanged.

Open the independent [six-donor QuPath project](../data/reviews/qupath-endocrine-rollout-20260918/project.qpproj)
with **File → Project → Open project**. Each image already has its objects;
no import or segmentation rerun is needed. The [offline gallery](../data/reviews/endocrine-rollout-20260918/index.html)
provides full-field overlays, deterministic small/medium/large examples, and
PNG/PDF plots. Research results remain outside the application.

| Donor | Current boundaries | Original cells retained | Boundary status |
|---|---:|---:|---|
| HPAP-078 | 1,431 | 267,993 | Owner approved; membership and measurements rebuilt |
| 6090 | 474 | 170,566 | Awaiting donor-specific visual review |
| 6232 | 139 | 92,012 | Awaiting donor-specific visual review |
| 6520 | 480 | 102,526 | Awaiting donor-specific visual review |
| 6539, recovered acquisition | 227 | 160,837 | Awaiting donor-specific visual review |
| 6558 | 710 | 149,134 | Awaiting donor-specific visual review |
| Total | 3,461 | 943,068 | HPAP approved; 2,030 new candidates for review |

In HPAP, cyan `Islet` objects retain the exact approved IDs, names and outlines.
The 9,569 cells inside them are reparented using native ROI containment. All
cell/nucleus shapes, cell classifications and raw measurements are preserved;
the two islet-distance measurements are recalculated. Tissue-clipped
`IsletExpanded` objects provide 20-µm context. Grey `Previous segmentation`
objects retain the older geometry.

The other five donors have cyan `Endocrine pixel candidate` objects. Their
existing cell parents and regional measurements still refer to the grey
previous segmentation and must not be used as measurements of the new outlines.
HPAP approval does not automatically approve these donors or any cell phenotype
thresholds. Each candidate remains editable in the new project.

## Fixed method and source preservation

The forest trained on original-15 reference material is unchanged, SHA-256
`a9fcdfc2cb38350946d9b762b61305166c06105722a55a22d9ad70c4d1135186`.
All five predictions use the same 160/160-vote operating point as the approved
HPAP objects and INS/GCG/SST only. The [IO60 size-dependent morphology](size_adaptive_islet_boundaries_20260917.md)
is unchanged: 1/2/3 µm erosion, 85% restoration, 50 µm² fragments, hole filling,
5 µm image-edge exclusion, and valid-tissue clipping at the verified
0.3774 µm/pixel calibration. No disease labels or phenotype gates enter this fit.
The incompatible CD99-dependent IO60 classifier is not substituted.

Predictions, unfiltered foreground, excluded fragments, geometry audits,
source fingerprints, logs and copied method sources are under
`data/development/endocrine-rollout-20260918/`. The corrected, distinct 6539
acquisition is retained. Original projects and previous review edits are
preserved. The original application release and all-donor preview are not
replaced by this review delivery.

## Native coordinate precision and verification

Independent checks found four very thin intersections introduced by native
serialization in 6090, 6520 and 6558. Their combined area is approximately
0.000023 µm². All disappear when eroded by the coordinate-displacement bound
derived from QuPath's 0.01-unit grids and source-coordinate float precision.
Actual overlap before serialization remains zero. Each residual is reported;
the check rejects intersections extending beyond that bound. Geometry validity,
minimum fragment size, tissue coverage and source identity are also checked.

A separate cell-membership check identified one HPAP cell exactly on the right
edge of `AdaptiveIslet_78_00247`. The native float coordinate is 7055.58984375 px;
the ordinary GeoJSON writer represents the edge as 7055.59 px, incorrectly
making this cell interior in reconstructed geometry. The measurement export
now preserves native vertices and applies each native ROI's boundary convention:
PolygonROI winding containment versus GeometryROI boundary-inclusive containment.
The 1,430 PolygonROIs and one GeometryROI then reproduce all 9,569 native
memberships exactly. No outline or cell coordinate was moved to make counts agree.
The initial rounded-export discrepancy is retained as audit evidence.

HPAP core/peri tables use all 31 raw markers, `median(log1p(raw intensity))`,
and tissue-clipped 10/20/30-µm peri bands that exclude all core interiors.
Structural densities use unioned positive mask area divided by valid area.
Overlapping 30-µm footprints are flagged for primary spatial exclusions.
All 267,993 raw cells remain available; the 15 with no marker measurements
are explicitly excluded from marker-dependent summaries. Phenotype availability
comes from the separate gate registry; boundary approval never changes it.

Final evidence is recorded in
`data/reviews/qupath-endocrine-rollout-20260918/delivery_manifest.json` and
`data/development/endocrine-rollout-20260918/regions/78/verification.json`.
It includes fresh six-image project loading, native object save/reload checks,
independent geometry checks, source hashes, regional cross-checks, measured
storage, and runtime configuration. Historical passing evidence is not reused
as acceptance of these new files.

The fresh checks passed: 32 Python regression tests, 10 native morphology cases,
four import-stage cases, all six native object round trips, and all five new
donors' independent geometry checks. HPAP's 31 core marker summaries and every
core membership agree between the native project and offline tables.

Measured delivery storage is 3.55 GiB for the native project, 0.81 GiB for
analysis/diagnostics, and 0.043 GiB for the gallery, excluding linked source
images. Fixed-model prediction took 241–261 seconds per new donor; adaptive
geometry took 7.8–35.1 seconds per donor. HPAP regional calculation took 117.2
seconds. These phase timings exclude native loading, object hashing, rendering
and development time. Peak RSS was not measured; heap limits are recorded below.

## Reproduction and rollback

Use explicit fresh roots with these scripts, in dependency order:

1. `stage_adaptive_rollout.py`: snapshot the saved approved HPAP review and stage
   independent copies of all six donors. Record the owner's exact approval text.
2. `promote_reviewed_islets.groovy`: promote only that hash-bound HPAP snapshot;
   rebuild native membership, core marker summaries, expansions and distances.
3. `refine_endocrine_pixel_classifier.py predict --minimum-votes 160`: predict
   each of the other five explicit target specifications with the fixed model.
4. `vectorize_endocrine_pixel_regions.groovy`, then
   `build_size_adaptive_islets.groovy`: trace foreground and apply IO60 rules.
   Verify each output with `verify_size_adaptive_output.py`.
5. `import_endocrine_pixel_review.groovy PROJECT CANDIDATES CASE`: import each
   selected donor. Independent native files may run concurrently; a file lock
   serializes the shared project/palette metadata updates. Do not run two
   imports for the same donor or edit the project during import.
6. `export_exact_reviewed_geometry.groovy`: export approved native coordinates
   with ROI type provenance. `rebuild_reviewed_islet_regions.py` takes that
   explicit export plus cohort, gate registry and fresh output root.
7. `render_endocrine_rollout.py`: render the offline gallery after the regional
   verification passes. Finalize hashes and capacity with
   `finalize_endocrine_rollout.py`.

QuPath uses its installed bundled JVM, headless AWT, and workspace-local Java
preferences. Predictions run on CPU; each native import was limited to a
24 GiB Java heap. Raw full-resolution images are referenced, not duplicated.
The gallery uses the existing loopback-only server on port 3840. Its log is
`data/development/endocrine-pixel-refinement-20260917/review_http/server.log`.
If that server is stopped and port 3840 is free, restart it with:

```bash
/usr/bin/python3 -u -m http.server 3840 --bind 127.0.0.1 \
  --directory /data/panc_CODEX/Islet-Explorer/data/reviews
```

Preserve services on ports 3838 and 3839. Rollback consists of reopening the
earlier independent project; no restoration or overwrite of source data is
necessary. Generated data, native objects, plots and logs remain Git-ignored;
source scripts, regression checks and documentation are the reviewable changes.
