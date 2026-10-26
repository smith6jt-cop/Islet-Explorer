# Size-adaptive endocrine boundaries

**18 September update:** the owner approved the 1,431 HPAP-078 conservative
objects. See the [approved-boundary rebuild and five-donor rollout](endocrine_boundary_rollout_20260918.md)
for the current project and measurements. The comparison below records the
17 September delivery and retains its original, then-unreviewed status.

The owner directed this correction to use the fixes in
`/home/smith6jt/IO60panc2nd/scripts/1_panc_objects.groovy`. The active boundary
workflow now applies its size-dependent geometry rules to the saved endocrine
pixel foreground. The classifier and its previous operating points are held
fixed for comparison. It does not require one uniform geometry operation to
serve every islet size.

| Foreground region section area | Erosion | Subsequent dilation |
|---|---:|---:|
| Below 2,000 µm² | 1 µm | 0.85 µm |
| 2,000–9,999 µm² | 2 µm | 1.7 µm |
| At least 10,000 µm² | 3 µm | 2.55 µm |

The source workflow also specifies a 5 µm image-edge exclusion, a 50 µm²
minimum fragment area, and filling of enclosed holes (`maxHoleSizeMicrons=-1`
in QuPath). The installed QuPath 0.7 implementation confirms that a negative
hole threshold fills all enclosed holes. These rules replace the previous
uniform 2 µm closing / 80 µm² fragment-and-hole cleanup for current boundaries.
Raw foreground, rejected fragments, parent-to-split relationships, filled area,
and the old comparator remain available in the offline artifacts.

## Adaptations needed for this cohort

- Keep the verified 0.3774 µm/pixel calibration. Do not copy the IO60 script's
  `setPixelSizeMicrons(0.5078,0.5078)`. Its buffer calls operate directly on
  pixel coordinates despite using variables named in micrometers. The adapted
  helper transforms X and Y to physical coordinates before buffering, then
  transforms back; anisotropic calibration is covered by a regression check.
- Keep the INS/GCG/SST foreground model. IO60's `Islets2.json` also needs CD99,
  which is absent from this cohort. That classifier cannot be transferred by
  silently replacing or dropping its CD99 input.
- Union nested foreground islands after hole filling to avoid duplicate
  annotations. Clip to valid tissue and the image-edge exclusion again after
  dilation, so hole filling cannot restore excluded tissue.
- Partial dilation can make split regions overlap again near narrow necks.
  Remove the shared restored area from both candidates, retain the original
  foreground, and report the excluded area. No assignment depends on object
  order, disease status, cell gates, or clinical metadata.

Hole filling constructs region outlines; it does not turn every enclosed pixel
into hormone-positive tissue. Each candidate carries its foreground fraction
and filled/background area. Structural or hormone-positive area measurements
must continue to use actual classified positive pixels.

## Current delivery

The fresh HPAP-078 project is
`data/reviews/qupath-size-adaptive-islets-20260917/project.qpproj`.
The comparison gallery is served at
[the local review URL](http://127.0.0.1:3840/endocrine-adaptive-20260917/).
The existing app services and prior editable QuPath projects are preserved.

On HPAP-078, the unchanged 160/160 vote sensitivity produced **1,431** adaptive
candidate regions from 1,632 initial regions above the fragment threshold.
There were 205 split parents; 261 parents had no surviving fragment after
refinement. Ten dilation collisions accounted for 23.55 µm² of shared restored
area, excluded before delivery. The final outlines have no positive-area
overlaps. The previous count was 1,056 using a larger 80 µm² minimum and uniform
closing; counts therefore need not decrease after correction.

At the unchanged 80/160 cutoff, the same geometry rules produced 2,040 regions.
Both versions are shown at the same 30 previously selected locations. The
QuPath copy uses the conservative 160/160 version. The six historical
hormone-negative examples have zero overlap with that version; one has a small
overlap with the majority-cutoff version. The latter still visibly includes
some weak/background patches. These examples are not a specificity estimate.

The reported failure is now `AdaptiveIslet_78_00853`. All 267,993 original cells
are retained with their nuclei, measurements, classifications, and parents.
Existing regional measurements and cell memberships refer to the previous
geometry; they were not recalculated. The new regions remain image-review
candidates, rather than an accepted replacement for downstream analysis.
The classifier's prior weak-reference validation limits still apply.

The identical geometry helper was also exercised on all **5,382** preserved
reference shapes from all **15** original donors, including small, medium and
large regions. This is a geometry check using existing reference shapes, not a
new classifier run or manual boundary validation. Its outputs are comparisons,
not replacements for source annotations. The other five new donors were not
regenerated in this correction.

## Reproduction and verification

`data/development/endocrine-adaptive-20260917/` holds the two HPAP candidate
sets, per-parent audits, reference checks, script snapshots, logs and delivery
verification. Initial attempts rejected dilation overlaps and are retained as
failed logs. Only the completed manifests declare usable review outputs.

The active implementation is
`scripts/cohort_expansion/size_adaptive_islet_geometry.groovy`. It is shared by
target construction and reference checks. The following command consumes a
hash-verified, natively traced pixel prediction and writes a **fresh** folder:

```bash
env JAVA_TOOL_OPTIONS='-Xmx12g -Djava.awt.headless=true' \
  /home/smith6jt/QuPath/bin/QuPath script \
  --args /absolute/traced_pixel_prediction \
  --args /absolute/fresh/adaptive_boundaries \
  --args /absolute/scripts/cohort_expansion/size_adaptive_islet_geometry.groovy \
  --args /home/smith6jt/IO60panc2nd/scripts/1_panc_objects.groovy \
  scripts/cohort_expansion/build_size_adaptive_islets.groovy
```

For a new image, `refine_endocrine_pixel_classifier.py predict` creates the
foreground and a legacy uniform comparator; `vectorize_endocrine_pixel_regions.groovy`
traces their source coordinates. Then run the builder above. New prediction
manifests mark adaptive refinement as required, and the importer rejects those
predictions until it is complete. Older saved manifests remain readable for
reproduction. Import only into a fresh independent project created by
`prepare_review_project.py`; never rerun the importer on a delivered project.

Run the native geometry regressions with QuPath's bundled runtime, since the
system Java version is older than the installed QuPath classes:

```bash
env JAVA_TOOL_OPTIONS='-Xmx4g -Djava.awt.headless=true' \
  /home/smith6jt/QuPath/bin/QuPath script \
  --args /absolute/scripts/cohort_expansion/size_adaptive_islet_geometry.groovy \
  scripts/tests/test_size_adaptive_islet_geometry.groovy
```

Ten native checks cover exact size thresholds, physical calibration, splitting,
dilation collisions, hole filling, nested islands, tissue exclusions, image
edges, fragment accounting, empty input, and deterministic ordering. The nine
existing endocrine pixel regressions also pass. Native review serialization
must preserve the complete cell-record hash and every original annotation
geometry; its evidence is `pixel_review_verification.json` in the new project.
Plots and CSVs are saved outside the app, with HTTP link verification recorded
alongside the delivery evidence.

Independent Shapely checks also verify the serialized GeoJSON, including
non-overlap, area thresholds and tissue/image clipping. QuPath's finite
coordinate precision uses its installed 0.01-unit geometry grid in both the
physical and source-pixel frames, plus float-coordinate rounding. The resulting
displacement tolerance is about 0.0101 µm for HPAP-078; actual subpixel fringe
areas and pre/post-serialization areas are reported rather than silently
equated. This is much smaller than the 0.7548 µm prediction pixel.

At delivery verification the new candidates, reference checks, plots and native
review copy occupied 1,143,696,021 bytes (1.065 GiB), excluding linked source
images and the previously fitted model. Java heap ceilings were 12 GiB for
target construction, 8 GiB for reference checks, and 24 GiB for native review
import. All 59 local page/asset links passed HTTP checks. This offline review
adds no app runtime or UFRC deployment requirement.
