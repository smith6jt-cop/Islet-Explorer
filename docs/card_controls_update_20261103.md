# Card controls update — 2026-11-03

This local review successor implements the owner's November 3 card-layout request.
The [user guide](user_guide.md) explains how to use it; the
[local runtime record](local_duckdb_review_20261103.md) identifies the running image.
Public dev and production are outside this update's scope.

## Behavior

- **Plot card titles:** both left and right headings identify the chosen data
  layer and plot type. Both plot titles identify the selected measure, with count
  or density wording when those metrics are selected.

- **Composition By Donor:** a compact Region / Order Donors By / Phenotype Colors
  toolbar and collapsed **Include Phenotypes** panel replace Leave Out. The grouped
  checkboxes match Spatial. Existing region defaults are preserved by inversion;
  an empty selection remains empty. The subtitle is removed, and the y-axis reads
  **Proportion of selected phenotypes**. Palette changes stay synchronized with
  Spatial and the cell viewer, including when object details are first opened or
  recreated on another tab. The CSV follows the selected rows, retaining original
  `proportion` and adding `proportion_selected` and `shown_cells`.
- **PAGA Force-Directed Layout:** two equally sized, independently colored plots
  share one fitted graph. Each offers Pseudotime, Donor Status, Leiden Cluster,
  Donor and Selected Feature. They retain linked hover/click selection and CSV
  export. The shared guide placement uses coordinates supported by the pinned
  Plotly runtime and reserves equal space below both panels.
- **Tissue Scatter Plot:** responsive legend columns sit below a map that uses
  its available card width. Raster height comes from measured title, axes and
  legend dimensions plus the current spatial aspect, including fixed-aspect zoom
  limits. The cache uses the exact output dimensions. Zoom/reset and downloads
  retain their existing controls.
- **Spatial donor:** **Donor For Tissue And Immune Gradient** is visible above
  Data Controls. It selects the tissue donor and dotted line in Cell Distribution
  Around Islets. Other cohort summaries retain their own donor-status controls.
- **Distance-To-Islet Profile:** titled groups separate distance target, plot
  type, disease status, grouping, y-axis and displayed population. Marker-grid
  styling is limited to marker mode so it cannot override phenotype subgroups.
- **Macrophages:** **Separate by:** and CD68 / CD163 / LGALS3 / CD4 boxes appear
  beneath the Macrophage checkbox. When profiles are split, that checkbox toggles
  only the displayed children; it does not contribute an aggregate phenotype value.
  Marker selections remain synchronized across cards.
- **Headings:** every card-heading word is capitalized while acronyms retain their
  spelling. Explanations remain in the relevant **?** popovers.

## Scientific boundaries

This change builds on the Unclassified display projection at `18df4f8`. It does
not retype cells or modify accepted data. The v19 endocrine identities, original
counts and denominators, marker measurements, object membership and saved
embeddings remain unchanged. Additive selectors contain disjoint subclasses;
parent aggregates and hidden Conflict cells are not reintroduced by UI controls.

Composition bars normalize to selected categories. Their CSV retains both that
display percentage and the original-denominator percentage. PAGA coloring changes
are joined by object identity and do not refit, reorder or filter the graph.

## Source and verification

The first integrated implementation is `2fba122`, comprising grouped commits:

| Commit | Change |
|---|---|
| `6e00d8a` | Composition Include controls, palette and selected CSV |
| `62ad81a` | Paired PAGA views and shared guide placement |
| `79225a5` | Shared macrophage controls, palette synchronization and headings |
| `f98a682` | Spatial width, legend, donor and distance controls |
| `2fba122` | Measured fixed-aspect zoom sizing |

Targeted R checks cover selected-composition normalization and original denominators,
PAGA keyed joins, the absence of additive macrophage-parent values, and actual raster
map width with sparse/dense legends and fixed-aspect zoom. Browser and native image
acceptance are recorded separately; source checks alone do not establish either.

Later corrections preserve the current palette when dynamic object details are
created, apply the Macrophage master selection before Shiny synchronizes child
inputs, and wrap compact profile labels between complete marker names. The
packaged-data projection fixes (`337513e` and its All-T predecessor) are included;
they apply the same display mapping to both loading routes while retaining raw
biological aggregates. See [the Unclassified contract](unclassified_display_20261103.md).

The integrated `ba0944c` R run passed **795 assertions in 121 tests**, with **36 R
files parsed** and no test-framework failures, errors, warnings or skips. Its
initial predecessor run exposed three outdated display-test fixtures; their
correction retains source-immutability and hidden-conflict assertions while
checking the new unresolved leaf projection. The later palette/title source
passed the same full suite. The R startup message about tibble being built under
R 4.5.2 is separate from test-framework warnings.

Diagnostic browser runs use a read-only source snapshot over the `9633f4f1…`
image, whose baked source remains `18df4f8`. Those runs are not native image
acceptance. The PAGA `ba0944c` diagnostic passes 15 checks and 10 inspected card
screenshots; its two occluded ancillary segmentation captures are excluded from
visual-pass claims. Spatial `ba0944c` passes 11 checks and four directly reviewed
tissue rasters; the `ea177fd` Distance follow-up passes six checks and four
inspected screenshots. Composition diagnostics exposed the master-toggle timing
and dynamic palette-reset defects described above. Failed application checks and
harness-preparation failures are preserved alongside their successors.

Native image `c7608905…`, built from `ccda497`, passed **15 PAGA checks** and
**20 Plot-title checks**, with 10 PAGA and 20 Plot screenshots directly inspected.
The title checks select two populations across proportion, count and density,
two structural measures, both data layers and swapped plot types at 1200 and
1600 pixels. They wait for the actual requested titles and rendered dimensions.
Earlier title-harness attempts used the wrong Selectize label field or sampled
stale reactive outputs; those failed attempts remain preserved.

The same native image exposed two final presentation defects: the Distance
grouping selector lacked an accessible label, and excess Region-control width
forced the Composition ordering selector onto an avoidable extra row. `d6405d6`
adds a screen-reader label; `3c81075` reallocates the Composition toolbar width.
The latter removes 40 pixels above the plot at 1200 pixels in the browser probe
and the final native run; both tested widths now have a 34-pixel toolbar.
The full source comparison proves these three UI substitutions explain the
entire difference from `ccda497`: the other 102 tracked app files are identical.
The intermediate `d6405d6` image was built but never activated.

The final native `3c81075` / `eb6f3c32…` build passes the strict security gate
with **zero High/Critical findings**, **49/49 standard browser checks**, **12/12
Composition checks** and **12/12 Spatial checks**. Composition validates the
actual selected colors through first detail creation and Plot/Spatial remounts.
Spatial validates default and dense legends, brush/reset, the shared donor
selector's exact scope and compact Distance controls at both widths. A prior
concurrent Spatial attempt timed out waiting for raster stability; its cause
was not established. The final isolated run passes the unchanged stability
criterion and records stable pixels, dimensions and marker-input copies.

Final native acceptance is recorded in the linked local runtime record and
[card-controls receipt](evidence/card_controls_validation_20261103.json). Prior
checks retain their original source/image identities; they are not relabeled as
successor executions. These checks do not authorize or report a public deployment.
