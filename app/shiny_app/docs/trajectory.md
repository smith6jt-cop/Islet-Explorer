# Fixed all-donor trajectory

The 2026-09-17 revision uses runtime schema 3, with dedicated
`endocrine_counts.csv`, `phenotype_availability.csv`, `marker_panel.json` and
`count_manifest.json` artifacts. `quality_review.json` marks the preserved fit
as requiring islet-identity review; passing numerical stability checks does not
certify candidate regions. See the [current delivery record](../../../docs/trajectory_revision_20260917.md).

The application reads the selected release's `trajectory/coordinates.csv`,
`features.csv`, `counts.csv`, `edges.csv` and `validation.json`. It never fits
pseudotime inside a Shiny session. The old composition-, ND-root- and
beta-orientation-based implementation is archived under
`archive/trajectory-pre-20260916/` and is not a fallback.

The y-axis is only `log1p(selected count / core endocrine union count)`. Core,
peri-islet and combined numerators share the same core denominator. That
denominator is the unique union of validated alpha, beta and delta assignments;
the retired Endocrine residual is excluded. Unknown and zero denominators both
produce missing ratios. Downloads retain the numerator, denominator, ratio and
plotted value, including unavailable rows.
Downloads also expose validation flags and the reason a ratio is unavailable.

All 31 panel markers are selected in the heatmap by default, including
INS/GCG/SST/CHGA. Select all/Clear controls preserve the panel order. Bins show
the mean of per-islet median log1p core intensity. Marker exclusion from fitting
the axis does not exclude it from descriptive display.

The donor menu defaults to all donors and lists donor IDs with disease status.
Donor and cell-count filters affect points, descriptive curves, marker heatmaps,
graph highlights and all downloads, while preserving the full-reference
coordinates, percentile and axis limits. The full graph remains faintly visible.
Graph point clicks retain donor/islet keys and open the correct image, including
HPAP-078. Single-donor curves have no donor-population confidence intervals or
between-donor significance tests.

Current candidate subtype ratios are unavailable because upstream phenotype
gates are not certified. The graph and separate marker heatmap remain available
for provisional inspection while islet identity is reviewed. Research results
are written outside the app; no Research tab or report route is registered.
The operational QuPath hormone-positive union used in H1 and
H9 does not certify historical subtype assignments and is not silently used as
the trajectory denominator.

Unavailable ratios appear as an informational notice, without an empty plot
or a Shiny error. A donor/cell-count selection with no islets returns empty
downloads with their column headers and keeps the full reference graph visible.
Clearing all heatmap markers prompts for a marker; restoring the selection
restores the heatmap without changing any fitted coordinates.

The [offline workflow](../../../scripts/research/README.md) documents features,
landmarks, scaling, diffusion graph, root assumptions, projection support,
donor-held-out validation and uncertainty. Stability is visible in the UI and
does not establish an unbiased disease axis or a measurement of elapsed time.
