# Islet Explorer — Journal Reviewer Guide

This guide describes the frozen PubApps reviewer release. It does not describe dormant
development modules.

## Access and data handling

- Requested URL: `https://islet-explorer.rc.ufl.edu/islet-explorer/` (pending UFRC
  provisioning confirmation; use the access link supplied by the journal).
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

The cell-count range is part of the analytical definition: changing it refits the
trajectory on the selected islets rather than merely hiding points. The default lower
bound is 10 cells per islet. Pseudotime is an ordering of compositional similarity, not
elapsed time or evidence of a causal transition. Donor identity is an important source of
variation, and the view should not be read as proof that every donor follows one uniform
ND → Aab+ → T1D path.

The PAGA coordinates, UMAP coordinates, and multi-feature heatmap rows can be downloaded
as CSV; publication-sized plots can be downloaded as PNG. Click an islet to open its
segmentation and cell view.

### Viewer

Inspect the donor OME-TIFF image with the bundled Avivator viewer and configurable marker
channels. The browser requests only needed byte ranges from the large image. First access
to an image may be slower than later navigation. Clicking an islet in Plot or Trajectory
selects the corresponding donor image when that mapping is available.

### Spatial

Inspect donor tissue cells, islet Leiden clusters, peri-islet phenotype enrichment,
infiltration, and signed distance from islets. Large tissue views are rendered as static
plots for stability. Filtered spatial rows and the displayed aggregate results can be
downloaded where a CSV button is present.

## Features intentionally absent

AI, Statistics, and Inflammation are disabled in the journal-review release. The service
makes no Navigator/OpenAI request and has no production API key. Their source modules may
remain in the repository, but they are neither mounted in the reviewer UI nor initialized
by its server.

## Reproducibility and problem reports

The reviewer site is pinned to one Git commit, OCI image digest, and checksummed data/image
release for the duration of a review round. The project will document any emergency
security update that changes reviewer-visible behavior.

If an action fails, contact `[SUPPORT NAME / EMAIL]` with:

1. Date, time, and time zone.
2. Browser and version.
3. Active tab and selected filters/feature.
4. The visible error text or a screenshot that does **not** contain the password.

Do not send the password or downloaded scientific rows in an ordinary email or issue.
