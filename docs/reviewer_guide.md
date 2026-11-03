# Islet Explorer — Journal Reviewer Guide

This guide describes the publication-profile interface. It does not describe dormant
development modules or certify a scientific data release.

**Current status, 2026-11-02:** the [selected-component v2 build](component_build_v2_20261102.md)
is live as `pubapps-2026-11-components-v1`, verified at 2026-11-03 02:42 UTC
(November 2 local) with 15/15 authenticated-access and 49/49 browser checks.
It preserves the exact named v19 endocrine cells and the isolated DG Lymphatic revision. Consult
the [deployment record](container_security_20261102.md) and the supplied review
round's release/access notice for actual availability.

The [workflow user guide](workflow_user_guide.md) explains data processing, component
selection and release checks; the [user guide](user_guide.md) describes controls in more
detail. A passing software, security or packaging check is not scientific acceptance.

## Access and data handling

- Live URL: `https://islet-explorer.rc.ufl.edu/islet-explorer/` (authentication required;
  use the access link and release notice supplied by the journal).
- Username: supplied in the journal's confidential reviewer instructions.
- Password: supplied separately through the confidential reviewer channel.
- Supported clients: current desktop Firefox, Chrome, Edge, or Safari. The interface is
  designed for a desktop viewport at least 1,200 pixels wide.

The site is classified **Sensitive** during peer review. Do not share its credential,
direct image links, screenshots containing case-level information, or downloaded data
outside the journal's review process. Report a suspected credential disclosure to the
support contact immediately.

The application accepts no uploads and cannot modify its scientific release. CSV and PNG
downloads are generated in temporary session storage and returned through the browser;
downloading a file does not alter the hosted data. Files saved to your computer remain
Sensitive reviewer material and should be handled according to the journal's policy.

## Reviewer-visible tabs

The publication profile contains **Plot, Trajectory, Spatial and Viewer**, in that order.
Counts, phenotype choices and aggregate definitions follow the loaded root's manifest and
frozen label groups. The v19/v20 comparison cohort has 20 donors (7 ND, 7 Aab+, 6 T1D);
older 15-donor counts and the legacy 19-rule CSV are not the current data contract.
The frozen comparison roots are historical references. In the selected-component
release, endocrine cell totals use exactly Beta, Alpha, Delta and CHGA+
hormone-negative. Anatomical Islet/Cluster/Single objects are distinct from cell
identities. Cells without a supported assignment remain visible as Unclassified
in overall counts, plots and downloads; they do not enter endocrine cell totals.
An object with no named endocrine core cells retains its structure and counts,
but its endocrine composition is unavailable. A positive region from the
original image classifier does not imply a positive corrected mean for every
cell inside it. The method help records the selected assignment rule.

### Plot

Explore marker expression, target measurements, hormone fractions, reviewer-facing cell
phenotypes, peri-islet proportions, and structural microenvironment measurements by islet
diameter and donor status.
Age, sex, autoantibody, donor status, and diameter filters constrain the displayed cohort.
Click an islet point to open segmentation and single-cell drill-down. The plot and its
filtered summary/distribution data can be downloaded.

### Trajectory

Explore the composition-derived islet trajectory. The release uses phenotype-composition
features, scanpy diffusion pseudotime, and a PAGA-initialized ForceAtlas2 layout. Marker
expression or a cell-type/endocrine-cell ratio can be displayed against raw or
rank-normalized pseudotime.

The current fitted cohort comprises endocrine objects with at least ten core cells,
regardless of Islet/Cluster class; the other tabs' endocrine-object selector does not
change this definition. The selected marker changes the displayed measurement, not the
pseudotime axis. Pseudotime is an ordering of compositional similarity, not
elapsed time or evidence of a causal transition. Donor identity is an important source of
variation, and the view should not be read as proof that every donor follows one uniform
ND → Aab+ → T1D path.

The PAGA coordinates, UMAP coordinates, and multi-feature heatmap rows can be downloaded
as CSV; publication-sized plots can be downloaded as PNG. Click an islet to open its
segmentation and cell view.

For frozen v19/v20, marker values are per-image scaled measurements averaged over core
cells, not scVI expression. The heatmap separately standardizes its displayed markers.
CODEX INS intensity is not an approved quantitative between-donor hypothesis endpoint
under the recorded acquisition restriction; see the
[measurement review](ins_measurement_review_20261026.md).

### Spatial

Inspect donor tissue cells, islet Leiden clusters, peri-islet phenotype enrichment,
infiltration, and signed distance from islets. Large tissue views are rendered as static
plots for stability. Filtered spatial rows and the displayed aggregate results can be
downloaded where a CSV button is present.

Structure area, cell-area overlap, distance and structure-qualified cell identity are
different measurements. In the v20 method, marker positivity alone does not assign a
generic Endothelial or Neural type: those identities require their corresponding
structure qualification. The current component-selection review must reconcile this
rule with the selected v19 endocrine cells before a new scientific release.

### Viewer

Inspect the donor OME-TIFF image with the bundled Avivator viewer and configurable marker
channels. The browser requests only needed byte ranges from the large image. First access
to an image may be slower than later navigation. Clicking an islet in Plot or Trajectory
selects the corresponding donor image when that mapping is available.

## Features intentionally absent

AI, Statistics, and Inflammation are disabled in the journal-review release. The service
makes no Navigator/OpenAI request and has no production API key. Their source modules may
remain in the repository, but they are neither mounted in the reviewer UI nor initialized
by its server.

## Reproducibility and problem reports

A reviewer release is pinned to one Git commit, OCI image digest, and checksummed data/image
release for the duration of a review round. The project will document any emergency
security update that changes reviewer-visible behavior.

If an action fails, contact `[SUPPORT NAME / EMAIL]` with:

1. Date, time, and time zone.
2. Browser and version.
3. Active tab and selected filters/feature.
4. The visible error text or a screenshot that does **not** contain the password.

Do not send the password or downloaded scientific rows in an ordinary email or issue.
