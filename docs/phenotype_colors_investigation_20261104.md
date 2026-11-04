# Composition-by-donor phenotype colors — 2026-11-04

The production PubApps container inspected on November 4 was image
`sha256:16c4979fff4d7e8d72befbaade013ce668153dc6161f938c8bc6a49a60e5ae24`,
built from `1659c9d`. The earlier color change `e2b514e` is an ancestor of
that source and its shared `phenotype_palette_for()` function is present in
the running container. The issue was not a lost commit or a stale image.

That change wired the shared resolver into Spatial and the single-cell views,
but did not change `donor_composition_figure()`. The donor chart still looked
up names directly and used `#8f8f8f` for every missing color. The deployed
phenotype inventory contains `Endothelium Blood Vessel`,
`Endothelium Lymphatic` and `Stromal`, none of which had explicit entries in
the four palettes. All three consequently appeared gray in that chart.

The prior composition acceptance checked that some bar colors changed when a
palette changed and that the selected palette stayed synchronized; it did not
require every phenotype to have a distinct color. The native predecessor report
for `3c81075` recorded a passing palette check while its before/after arrays
both contained the same three gray values. That report predates `e2b514e` and
is evidence of the test gap, not evidence of a run on the later release.

This correction assigns explicit colors to the three deployed labels in all
four palettes and routes donor bars through the shared resolver. The resolver
also spaces fallback hues for new labels around the color wheel. The browser
acceptance now checks for repeated colors and the gray fallback, and compares
the affected colors with the Spatial cluster chart. Unit coverage checks all
16 macrophage marker combinations across the four palettes and verifies that
known colors remain stable under reordered and narrowed selections.

The checkout's existing R suite initially failed on stale test fixtures for the
current Unclassified-selector policy. The same failures reproduced against the
unchanged baseline. The affected fixtures now load the current helpers and
assert the current selection behavior; source cell records and scientific
calculations were not changed.

Local verification used the production schema-5 data in a read-only container
with the candidate application source mounted as an overlay. It exercised all
three regions and all four palettes in the browser, checked a narrowed
selection and Spatial agreement, and checked all 192 region, palette and
macrophage-split combinations against the production composition summaries.
These checks are local source-overlay checks, not evidence that a new image
has been built or that production has been updated.
