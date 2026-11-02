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
review root and are hash-bound in `uniform_release.json`. Diagnostic component
roots are rejected by the release builder and parity gate, even if someone
changes their displayed label-group definitions. There is no acceptance override
for this diagnostic protocol.

The [released-fixes security image](container_security_20261102.md) is already
built and staged offline. DuckDB remains deferred at the owner's direction and
fails the full image gate. Both PubApps services remain stopped and masked.

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

Execution results and artifact identities will be recorded here after the
builds finish. Source inputs and historical outputs are preserved; new scientific
computation runs through the sealed HPG SLURM environment.
