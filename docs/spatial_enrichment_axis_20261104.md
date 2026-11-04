# Spatial enrichment axis and donor means — 2026-11-04

The Phenotype Composition & Enrichment card formerly defaulted to dropping
object-level enrichment z scores outside ±5 before calculating each donor's
mean. The control has been removed. Donor means now include every finite
object value for the selected region and baseline, so summaries, heatmaps and
CSV values may differ from the earlier filtered defaults.

The bar chart now shows each donor mean as a point alongside the status
summary. Its **Axis range** selector defaults to **Fit 99% of donor values**;
**Fit every donor value** shows the full range. The fitted option uses the Plot
tab's 99th-percentile rule for each tail, while keeping summary bars and their
intervals visible. Donor values beyond that display range appear as open
triangles with their actual values in hover. The selector changes no donor
means, summaries, heatmap values or CSV data rows.

The use of donor values in the axis label is deliberate: this chart summarizes
one mean per donor and phenotype, while the Plot tab shows individual objects.
The CSV includes the selected display range in its provenance and retains the
unfiltered donor and summary values under either setting.

In the production schema-5 composition table, the CD8 phenotype has 4,848
finite object z scores, including 220 outside ±5. A browser check of the local
source overlay matched all 20 displayed CD8 donor means to independent means
computed from those unfiltered source values. It also verified identical bars,
heatmaps and CSV rows across the two axis settings for cohort and own-tissue
enrichment, core proportion and peri-islet proportion. This is local preview
evidence and does not establish a new public deployment.

**Deployment update, 2026-11-04:** the native image and data bundle were subsequently verified and activated on PubApps; see the [live receipt](evidence/pubapps_live_20261104_ui_corrections.json). The local-only observations above retain their original scope.
