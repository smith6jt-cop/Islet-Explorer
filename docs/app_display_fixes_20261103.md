# App display corrections — 2026-11-03

The owner requested consistent parent/subtype categories, no visible Conflict
classes, readable Plot dropdowns, removal of the Trajectory citation suffix,
even Trajectory card spacing, and Spatial scatter legends below the plot.
These are app changes against the accepted component-v2 scientific data.
They do not retype cells or change the selected v19 endocrine evidence.
The frozen app source is **`81b3857b45e7870fd4d43dee3751c8d956619052`**.
It includes the reviewed DuckDB security patch and a metadata-only Containerfile
change that permits dependency caching across future app releases. The patch,
dependency installation and strict runtime security gates are unchanged.

The current live dev release is documented separately in
[the strict DuckDB deployment receipt](dev_duckdb_deployment_20261103.md).
Its `378843e` source predates these app corrections. This record does not claim
that the corrections are live until a new immutable image is built and verified.

The [review follow-up receipt](evidence/app_display_followup_20261103.json)
records two corrections found by Senior's independent review of `5e6580f`:
mixed raw/profile macrophage inputs now retain all merged metrics and calculate
enrichment from the original member counts; a wider right plot margin keeps the
last tissue-map axis label inside the image at 1200 px. The accepted cohort has
profile labels only, and all 16 marker selections remain bit-identical to the
previous source. No scientific source records or frozen data changed.

On the corrected source, all 35 R files parse and the full R suite passes
**601 assertions**, with zero test-framework failures, warnings or skips.
All 48 real-data donor-composition scenarios pass again. Four static tissue-map
checks (default/all 27 labels at 1200/1600 px) show complete axes and legends,
with visible maps **401–483 px** high. Brush zoom and exact reset were tested
for the all-27, 1600 px case. Both download links remain available below the
plot in all four cases; downloaded file contents were not tested in this check.
The earlier `5e6580f` receipts below remain unchanged historical evidence.

Senior / MultiPlex5000 independently closed both findings in recipient commit
`20f869ed839fe9874b99c4834dbbe58cb35191c0`: 138 focused assertions, 2,048
synthetic helper cases and inspection of all four final tissue images. The
[recipient closure receipt](evidence/app_display_senior_closure_20261103.json)
binds that review and its limits. Senior reviewed the supplied real-cohort
results without independently recomputing them, and did not approve or perform
a deployment. The image, installed-runtime and public-site acceptance are
separate release gates.

A subsequent [execution-provenance supplement](evidence/app_donor_execution_supplement_20261103.json)
repeats all 48 real-data scenarios against exact `81b3857` source snapshots,
with 960 donor-keyed comparisons. The runner, all ten consumed application
sources and all eight data inputs are bound before, during and after execution.
All comparisons pass and all inputs remain unchanged. This closes the earlier
report's execution-provenance gap with a new run; it does not retroactively
change the old evidence or Senior's recorded review scope.

## Presentation policy

The older retirement list covered Endocrine, T cell and Immune, but omitted the
newer hematopoietic pan-marker residual labels. Some paths also filtered only
palette entries, leaving hidden cell points on the plot. Group totals were
available to additive checkbox views, where they overlapped their children.
Shared helpers now enforce the policy at those boundaries.
The immutable pair lookup also omits Conflict labels whose slashes were
sanitized in peri-islet column names. Donor-composition accounting therefore
reads hidden counts directly from the original regional columns, independently
of that display lookup.

- Broad residual types are absent from phenotype choices and legends. They stay
  in source counts and the existing **No subtype** donor-composition category.
- Parent totals such as **All T cells** are single-measure queries over their
  source members. They are absent from additive stacks and checkbox categories.
  Changing the macrophage marker split preserves a disjoint partition.
- Conflict classes are absent from app choices, phenotype plots, tables, hover
  labels and downloads. Marker-mode display subsets also omit their cells;
  all-cell denominators and cached source measurements remain intact.
- Conflict counts are excluded from **No subtype**. Donor bars sum to 100% of
  displayed cells; hover and CSV preserve proportions of all original cells.
  Those original-denominator CSV proportions can therefore sum to less than 100%.
- Tissue, Unclassified, Doublet suspected and terminal anatomical labels are
  preserved. Hiding a class is not evidence for assigning another cell identity.
- Live DPT/PAGA fitting excludes the same retired and Conflict feature columns.
  This can change a newly fitted ordering; stored feature ratios, original
  denominators and frozen UMAP/Leiden embeddings are unchanged.

The full-cohort check accounts for all **3,199,606** source cells:

| Presentation category | Cells |
|---|---:|
| Visible named or terminal labels, before No subtype | 3,091,614 |
| General residuals: T cell / hematopoietic CD45+ / CD45− | 98,779 |
| Hidden Conflict labels | 9,213 |

The residual counts are 3,397 / 50,698 / 44,684 respectively. Source parent totals
remain 110,721 endocrine cells, 36,864 T cells and 190,096 macrophages. Additive
memberships are disjoint. The source manifest and labels/counts are unchanged.
The [cohort receipt](evidence/app_phenotype_hierarchy_20261103.json) records the
membership check; the [donor-composition receipt](evidence/app_donor_composition_20261103.json)
records the separate object-region checks and exact input/helper hashes.

## Layout

Plot selections wrap long labels and reserve a gutter for the down arrow,
including the donor-composition controls. Pseudotime and UMAP cards form an
aligned row; the Segmentation map follows in a full-width row with consistent
20 px spacing. The visible PAGA title and export title omit the requested
`(Barlow et al. 2024)` suffix; scientific method citations remain in help.
Spatial tissue and islet scatter guides sit below their plots, with the UMAPs
using the full horizontal plot domain. The two UMAPs reserve the same measured
guide space, including long donor legends, and colorbar titles clear their
gradients by 6 px. Long tissue legend labels wrap; the tissue output grows
vertically with dense legends instead of shrinking the map. The cached plot
uses an explicit height policy because Shiny's plot cache otherwise overrides
the renderer with the client's initial zero/short auto height.

## Verification

Focused checks cover retained labels, exact parent membership, mutually
exclusive macrophage splits, export filtering, original donor denominators,
drilldown points, both CSV and DuckDB Spatial queries, and live trajectory
invariance to hidden feature changes. Browser checks use the immutable accepted
schema-5 data on a transient loopback app; no host service or source data is
changed by that check. The [validation receipt](evidence/app_display_validation_20261103.json)
binds the frozen source, test logs and the separate scientific/presentation
receipts. All 73 captured tracked app/test source files match the frozen commit;
two ignored test reproducers are explicitly excluded from that binding.

At `5e6580f`, all 35 R files parse; the complete R suite passes **588 assertions,
zero test-framework failures, warnings or skips**. The three Python trajectory
tests pass. Logs retain a tibble R-build-version startup warning and an
anndata/zarr deprecation warning outside those summaries. A
read-only check of the real release covers **48 donor-composition scenarios**
(all 16 macrophage marker selections × core/peri/both), with original
denominators and displayed proportions verified for every donor. Recomputed
aggregate columns match the frozen release exactly. Both runtime loader paths
and all 18 release-builder helper dependencies pass.

The [visual receipt](evidence/app_display_visual_20261103.json) records the
1200/1600 px checks: at least 7 px dropdown-arrow clearance, exact 20 px
Trajectory gaps, no hidden phenotype choices or Conflict plot text, and six
paired Spatial color/width cases with equal areas and 1,392 points per plot.
Default and full 27-phenotype tissue legends retain a **407–483 px visible map**,
complete legend labels, download links below the output, the same output
element, and working brush zoom with exact reset restoration in the all-27,
1600 px case. The subsequent review found the clipped 1200 px axis tick closed
by the follow-up above. Failed intermediate layout
receipts remain preserved and are superseded by the measured final check.

The expanded loopback smoke passed 45/49 checks. Its four Viewer/static-Range
checks need the production Shiny Server aliases, which the transient `runApp`
process does not provide. Full container and public-dev acceptance are required
before activation; this local receipt makes no deployment claim. Production
remains on its separately recorded release while the replacement dev build runs.

The real-release donor totals below count object memberships, not distinct
whole-cohort cells; a cell can belong to more than one object's neighborhood.

| Region | Original memberships | Hidden QC memberships | No subtype | Displayed memberships |
|---|---:|---:|---:|---:|
| Core | 129,249 | 1,269 | 1,190 | 127,980 |
| Peri-islet | 228,409 | 1,304 | 8,080 | 227,105 |
| Both | 357,658 | 2,573 | 9,270 | 355,085 |
