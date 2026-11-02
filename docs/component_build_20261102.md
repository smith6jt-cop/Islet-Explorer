# Component candidate build — 2026-11-02

The owner requested the builds after selecting v19 endocrine cells and requiring
the Lymphatic correction to preserve unrelated components. This build is a fresh
offline review candidate. It is not a scientific acceptance or a deployment.

The input selection is explicit:

- The archived v19 endocrine UUIDs and their subtype decisions, together with
  their INS, GCG, SST and CHGA measurements, levels and coverage.
- The initial v20 calibration for the other markers, before the Lymphatic DG
  geometry changed shared reference populations.
- The instructed v20 structure-only neural/endothelial rules, whole-cell
  positivity with membrane differentiation, and nuclear Ki67 state.
- The verified DG Lymphatic geometry and unchanged common geometry, raw cell
  measurements, QC and donor inventory.

The component builder records the natural mixed-calibration classifier result
before applying the selected endocrine component. It checks the archived named
endocrine calls against their original hormone evidence. It does not copy old
non-endocrine identity fields or silently restore v19 marker-only neural and
endothelial calls.

The exact archived v19 versus isolated v20 comparison contains 110,604 identical
named endocrine calls, 11 subtype changes, 106 losses and 2,519 gains. An
unselected natural endocrine call retains its archived unresolved object class
only where that was its actual v19 state. Other such cells are recorded as
unresolved component boundary conflicts. Those records are technical conflicts,
not new biological cell types.

The v19 archive also retains 7,151 unnamed cells and 953 populated objects
without a positive eligible member under its corrected-cell thresholds. The
build preserves that evidence. No threshold adjustment, invented subtype,
whole-version rollback or changed aggregate denominator resolves it.

All app tables, normalized values, spatial summaries and embeddings must be
regenerated from the new cell records. The export now carries the original
binary `Ki67_pos` state with its per-image measurement and threshold provenance.
Ki67 remains a state; it is not added to the 31-marker embedding panel or used
to create an identity.

`component_selection.json` and `component_selection_audit.json` accompany the
review root and are hash-bound in `uniform_release.json`. The release builder
checks their hashes, donor inventory, per-image sums, preservation of the selected
source and assembly census. Nonzero unresolved assignments or boundary conflicts
fail even if someone changes the displayed group list. The owner has authorized
the current live build: an older `diagnostic_pending_review` status by itself is
not a new approval requirement. Actual data and runtime checks still apply.

The matching released-fixes image is `41e62a97…`, built from `b15f0fc`. Its package
scan has zero High/Critical findings and its dependency and installed Shiny/R
runtime checks pass. The known DuckDB finding remains visible and is nonblocking
under the owner's [dated remediation window](duckdb_remediation_window_20261102.md)
through 2026-11-10 UTC. The local DuckDB fix is kept in the separate dev build.
Both PubApps services remain stopped while the current data checks are resolved.

The subsequent [separate dev DuckDB build](dev_duckdb_fix_20261102.md) is also
complete: image `b26558a9…` from branch commit `92c63ba` passes the strict decoder
and security gate without a deferral, plus read-only BROTLI round trips on real
v19/v20 release tables. It is retained locally and has not been activated.
This does not change the selected production image or the data failures below.

The executable entry point is `scripts/hpg/build_component_candidate.py --config
CONFIG.json`, run through `scripts/hpg/jobs/site_job.sbatch`. The exact config and
stage commands are copied into the fresh report root. The [workflow guide,
sections 6–7](workflow_user_guide.md#6-reconcile-selected-components-and-inspect-them-in-qupath)
describes the full sequence and config keys.

Pre-build validation passed across all 20 donors: 3,199,606 eligible UUIDs,
110,721 selected named endocrine calls, exact archived endocrine
positivity/canary/fold fields, and both source coverage receipts for every image.
The Ki67 export checks also pass over all eligible UUIDs. The updated isolation
auditor reproduces the previous all-20 passing diagnostic result. These checks
validate preservation and implementation; they do not establish image accuracy.

The initial-v20 source package contains 173 copied files (1,486,929,911 bytes)
and pins 303 shared files (44,081,533,658 bytes). Package manifest SHA256:
`a42b5576f2096f03f2c1493a2986f27f0ebd6a29fac3e5cdc0a21a316bebfd3a`.
The copied package is verified on HPG. Shared-file verification runs on the
allocated compute node before any source links are created.

## Completed builds

The scientific build used source `b15f0fcb4566dcfaf8f7a02fcaa88c550e3f3beb`
through the sealed HPG environment. Job `49141854` completed all eight stages:
source materialization, control preparation and identity, DG preparation and
identity, Lymphatic invariance, app assembly, and export-state verification.

| Artifact | New path relative to the repository |
|---|---|
| Mixed-calibration control | `data/development/uniform-20donor-20261102-component-control-v1` |
| DG candidate processing | `data/development/uniform-20donor-20261102-component-candidate-v1` |
| App data | `data/development/uniform-20donor-20261102-component-app-v1` |
| Build and export checks | `output/component_candidate_build_20261102_v1` |
| QuPath project | `data/development/qupath-review-20261102-components-v1` |
| QuPath assignment reports | `output/component_candidate_qupath_20261102_v1` |
| Strict object audit | `output/component_candidate_object_audit_20261102_v1` |

The control and candidate retain all 3,199,606 eligible UUIDs and the exact
110,721 selected named v19 endocrine calls. The mixed-control→DG comparison
passes for all 20 donors: **351 label changes**, all involving Lymphatic identity
inside the changed Lymphatic footprint; zero unrelated label changes. All nine
calibration tables, marker positivity/canary/fold fields, Ki67, raw/REDSEA inputs,
QC and non-Lymphatic context are unchanged. This is distinct from the earlier
initial-v20 diagnostic's 352 changes and 113,134 named endocrine cells.

App export checks pass for all 20 donors, 3,199,606 cells, 4,848 object CSVs
(357,658 rows), and 78,680 Ki67-positive cells. UUID/label/state associations,
donor-coordinate exports, distance exports and measurement/level provenance
agree with the candidate. The embedding input remains 31 markers. Five geometry
lookup objects without paired analysis observations correctly have no object
CSV; the inventory follows the H5AD's analysis-object IDs.

QuPath job `49145767` completed all 20 images. The separate strict object audit,
job `49146295`, exited 4: **953 populated objects** (915 Singles, 38 Clusters,
zero Islets) lack a positive eligible cell under the selected corrected-cell
thresholds. Two Singles have no eligible member, and one has two. This does not
mean the original pixel classifier detected signal-free objects: its saved
masks replay exactly. See the [classifier discrepancy](endocrine_classifier_discrepancy_20261102.md)
for the differing measurements and polygon-centroid effects.

The assembled candidate contains **8,629 unresolved object-class cell labels**
and **2,472 component boundary conflicts**, including all 7,151 archived unresolved
v19 labels. The latter conflicts identify natural endocrine calls outside the
exact selected v19 membership. They remain explicit pending a consistent
component-boundary rule; they are not exported as invented biological subtypes.
These actual data failures prevent release packaging. Successful assembly and
state preservation do not resolve them.

The matching offline image is
`sha256:41e62a979b7ebb88ee97048d63cebae9f87d8202660a4767217861bdc4330162`,
tagged `localhost/islet-explorer:component-review-duckdb-deferred-20261102`.
It passes the installed Shiny/R/JavaScript/static-image runtime smoke and
the security policy at `7941df1` with the dated DuckDB exception. It has not been
exported, staged on PubApps or activated. Source inputs and historical outputs
are preserved. The compact [build evidence](evidence/component_build_20261102.json)
records artifact hashes, completed checks and the separate failed data checks.
