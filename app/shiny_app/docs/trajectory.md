# Fixed all-donor trajectory

**Current schema-4 implementation, 2026-09-22:** cell typing precedes
trajectory fitting, and the trajectory uses endocrine-normalized typed-cell
composition. The expanded typing requirements and nine explicitly required
marker contributors are recorded in the
[reference/decision record](../../../docs/trajectory_typing_reference_20260919.md)
and `scripts/research/cell_typing_requirements-v1.json`.

The application reads fixed offline `coordinates`, `composition_features`,
`counts`, `endocrine_counts`, `edges`, `graph_layout`, marker summaries and
validation artifacts. Parquet is preferred and cached per worker; CSV mirrors
remain auditable. Release identity is part of the cache key. Shiny never fits,
rescales or reranks pseudotime during interaction. The old ND-root and beta-orientation
implementations are archived under `archive/trajectory-pre-20260916/` and are
not fallbacks.

Eligibility has exactly two requirements: the endocrine pixel-classifier core
does not contact an image edge, and at least one canonical cell centroid lies
inside that core. There are 11,136 eligible islets across all 21 donors, with at
least 137 per donor. Downstream typing does not override this classification.
Of the eligible islets, 7,493 currently have complete typed-composition
features. The largest graph component supplies 2,626 supported primary
coordinates; 14 disconnected components remain explicit and are not forced
onto a comparable axis. Median leave-one-donor-out absolute Spearman is 0.131
and supported primary coverage is 23.6%, so the saved ordering is **unstable**
under the prespecified checks. Parameters were not changed to improve disease
ordering.

The raw coordinate has 731 unique values among the 2,626 supported islets; its
largest exact tie contains 1,896 islets. The old empirical midrank collapsed
that tie to one x-position and produced large empty gaps. The default saved
reference coordinate is now the fixed midpoint grid `(rank + 0.5) / 2,626`:
raw pseudotime determines order, and SHA-256 of the immutable islet ID spreads
exact ties without donor, disease or outcome information. This makes every
supported islet visible across 0–1, but the order inside a raw tie has no
biological meaning. Hovers and downloads retain raw pseudotime, the raw-scale
control exposes the tie directly, and inferential code uses raw pseudotime.

The y-axis is only `log1p(selected typed-cell count / canonical cell count in
the classifier-determined endocrine core)`. Core, peri-islet and combined
numerators share this fixed core denominator. Each canonical core cell is
counted once regardless of its downstream phenotype support. A classifier core
with no cells is ineligible, and an unavailable phenotype numerator remains
missing. Downloads retain the numerator, denominator, raw ratio, plotted value,
eligibility, support flags and availability reason.
For the current primary component, all 2,626 default beta/core ratios are valid
measured zeros. The app states that fact and uses a narrow nonnegative y-range
instead of visually implying missing data or positive variation.

All 31 panel markers are selected in the heatmap by default, including
INS/GCG/SST/CHGA. Select all/Clear controls preserve panel order. The default
measurement is mean log1p corrected expression within the selected typed
population, initially the endocrine union. Each bin first averages islets
within donor and then weights donors equally; colors are within-marker z-scores.
The marker-positive-count/E alternative remains explicit and gray where
reviewed positivity calls are unavailable. Heatmap measurements never enter or
change the saved trajectory.

The donor menu defaults to all donors and lists donor IDs with disease status.
Donor and cell-count filters affect points, descriptive curves, marker heatmaps,
graph highlights and all downloads, while preserving the full-reference
coordinates, percentile and axis limits. The full graph remains faintly visible.
The historical column name remains `pseudotime_percentile` for schema
compatibility, but the UI calls it a fixed reference-cohort rank.
Graph point clicks retain donor/islet keys and open the correct image, including
HPAP-078. Single-donor curves have no donor-population confidence intervals or
between-donor significance tests.

Direct islet search follows the current donor and cell-count filters, joins the
schema-4 donor cell partition to exact core/peri membership and typed identity,
provides outline and per-cell downloads, and opens the corresponding donor
image. Exact single-islet GeoJSON downloads avoid reading a whole donor's
geometry into the R session. Research results are written outside the app; no
Research tab or report route is registered. The operational QuPath
hormone-positive union used by historical research results is not substituted
for the schema-4 trajectory denominator.

Unavailable ratios appear as an informational notice, without an empty plot
or a Shiny error. A donor/cell-count selection with no islets returns empty
downloads with their column headers. A graph is shown only when the release
actually supplies supported coordinates.
Clearing all heatmap markers prompts for a marker; restoring the selection
restores the heatmap without changing any fitted coordinates.

The [offline workflow](../../../scripts/research/README.md) documents features,
landmarks, scaling, diffusion graph, root assumptions, projection support,
donor-held-out validation and uncertainty. Stability is visible in the UI and
does not establish an unbiased disease axis or a measurement of elapsed time.
