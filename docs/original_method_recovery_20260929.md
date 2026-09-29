# Original typing and scVI method recovery — 2026-09-29

The six additions cannot yet be processed with the exact method used for the
original 15 donors. The raw-only bundle has 943,068 objects and preserves its
copied 15-donor app artifacts, but it contains no accepted phenotype,
normalized-expression, scVI or common trajectory result. The frozen app remains
the default. This is a scientific dependency, not a runtime-loader problem.

## What is needed

- The exact original phenotype rule matrix and the numeric per-marker,
  per-image gates, with their input scale and application order. Applying the
  same labels to new donors requires checking that this combination reproduces
  the frozen original labels first.
- The trained CODEX scVI model weights and matching preprocessing: input feature
  order, count/intensity layer, transforms, categorical covariates, library and
  model settings. Saved latent coordinates alone cannot project new cells.
- The owner's visual decision on the cyan endocrine candidates for 6090, 6232,
  6520, recovered 6539 and 6558 in the
  [current QuPath project](../data/reviews/qupath-endocrine-rollout-20260918/project.qpproj).
  HPAP-078's 1,431 boundaries are already owner-approved; this does not approve
  the other five donors or their phenotype gates.

## Evidence checked

The repository has `single_cell_analysis/phenotype_rules1.csv` and the
app-facing `data/phenotype_rules.csv`, but the latter uses 15 markers and does
not reproduce all 21 phenotype labels in the shipped single-cell H5AD. The
original rule matrix has not been verified. `data/phenotype_gates_provenance.csv`
records later gate provenance; its existence does not establish that the entire
original typing chain is recovered.

`single_cell_analysis/SCVI_CODEX_v2.ipynb` contains a saved-model path pattern,
`$ISLET_DATA_DIR/CODEX_scvi_model_v2/CODEX_scvi_model_v2`, but its training cells
are unexecuted in the checked-in notebook. Its code/comments also disagree on
latent dimension and batch handling, so the notebook alone cannot substitute
for the trained model and matching run record. The frozen H5ADs carry saved
scVI representations, not the model weights needed to process the six donors.

Read-only searches under `/data/panc_CODEX`, `/blue/maigan/smith6jt` and
`/orange/maigan/smith6jt`, including plausible CODEX project/model locations,
found no verified original CODEX scVI weights or full rule matrix. A 1,394-byte
Blue rule CSV under `Islet-Explorer-Senior-main/data/` remains an unexamined
candidate; its content could not be transferred to this host. The only remote
`model.pt` hits identified were for Xenium, not CODEX. Large Orange acquisition
directories could not be exhaustively listed through the available Globus
endpoint, and creating a local transfer endpoint required interactive login.
Absence across all remote storage is therefore unproven.

## Resume condition

Once the exact artifacts are located, verify their identity and reproduce the
original donor labels and scVI representation before applying them to the six
new raw H5ADs. Keep original donor IDs, physical geometry, measurements,
missingness and phenotype labels fixed. A validated common 21-donor trajectory,
UMAP and Leiden refit may then change derived coordinates and clusters. Record
source-cell parity separately from the existing app-artifact parity check;
that check alone does not compare the multi-gigabyte original single-cell H5ADs.
The nine-marker typing enhancement remains deferred. Final 21-donor PubApps
instructions must wait for this integration and its review evidence.
