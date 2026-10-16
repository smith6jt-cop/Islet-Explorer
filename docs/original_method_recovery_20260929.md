# Original typing and scVI method recovery — 2026-09-29

> **2026-10-15 owner directive:** donor 78 (HDL078 / HPAP-078) is removed from the app and all analysis ([record](owner_directive_20261015.md); `scripts/excluded_donors.csv`). The 2026-10-15 search, which recovered the scVI training and import code and confirmed the rule file and model weights are lost, is in [the directive record](owner_directive_20261015.md#original-15-donor-files).

> **2026-10-03 owner decision:** the accepted reconstruction rules and all six donor boundaries are approved. Use Harmony on protein PCA, retain the uncorrected comparison, and preserve original source records. Exact historical typing/encoder recovery and scVI retraining are no longer prerequisites. See the [current implementation record](protein_expansion_20261003.md); earlier directions below are historical.

As of 2026-09-29, the six additions could not yet be processed with the exact method used for the
original 15 donors. The raw-only bundle has 943,068 objects and preserves its
copied 15-donor app artifacts, but it contains no accepted phenotype,
normalized-expression, scVI or common trajectory result. The frozen app remains
the default. This is a scientific dependency, not a runtime-loader problem.

## What is needed

- The exact original phenotype rule matrix and its postprocessing order.
  Applying the same labels to new donors requires first reproducing the frozen
  original labels from the saved marker values.
- The trained CODEX scVI model weights and matching preprocessing: input feature
  order, count/intensity layer, transforms, categorical covariates, library and
  model settings. Saved latent coordinates alone cannot project new cells.
- The owner's visual decision on the cyan endocrine candidates for 6090, 6232,
  6520, recovered 6539 and 6558 in the
  [current QuPath project](../data/reviews/qupath-endocrine-rollout-20260918/project.qpproj).
  HPAP-078's 1,431 boundaries are already owner-approved; this does not approve
  the other five donors or their phenotype gates.

## Evidence checked

The user's suggested
`single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad` contains the
frozen 2,649,654-cell, 21-label phenotype output and a 31-marker gate table
for all 15 original donors. After matching marker order, all 465 gate values
equal the earlier single-column gate table in
`CODEX_scvi_BioCov_gated.h5ad`: the original numeric thresholds are recovered.
For example, the `log1p(raw MFI)` thresholds are INS 5.3, GCG 7.75, ECAD 9.2,
CD45 6.9 and CK19 8.25. On 1,000 evenly spaced frozen cells, all 31,000
comparisons of saved rescaled `X >= 0.5` to `layers['log'] >= gate` agreed, and
the saved log layer matched `log1p(raw)` to floating-point precision.

`single_cell_analysis/CODEX_Panc_scimap_Analysis.ipynb` describes the order:
raw MFI, `log1p`, saved Napari gates, `scimap.pp.rescale`, then an external
`phenotype_rules.csv` through `scimap.tl.phenotype_cells`. The notebook's typing
cells are unexecuted in the checked-in copy; it is a method recipe rather than
proof of the historical run. The repository has
`single_cell_analysis/phenotype_rules1.csv` and the app-facing
`data/phenotype_rules.csv`, each with 19 rule rows. The latter uses 15 markers
and is only a partial upstream reference, not a verified complete explanation
of the 21 labels in the shipped H5AD. The exact rule matrix and subsequent
edits must reproduce the frozen labels before use on six additions. In a
read-only replay of 9,254 cells sampled across all 21 frozen labels and 15
donors, installed scimap 2.3.5 could not run Rules1 because it references
absent `CD56` (and `CD99`). Removing those two columns produces the app CSV.
That runnable table agreed with only 26.17% of frozen labels literally, or
37.98% after favorable label aliases and the documented CD8a/CD4 T-cell
correction. It called 5,851 sampled cells Ductal, including all 450 sampled
frozen Acinar cells, and cannot emit `ECAD+` or `Structural`. These are
diagnostic sample results, not population-wide accuracy estimates.

`scripts/fix_tcell_subtypes.py` describes a later 0.5 CD8a/CD4 refinement.
All 22,442 frozen T-cell-lineage labels are consistent with it. The saved
output alone does not prove when or in which run that correction was applied.

`single_cell_analysis/SCVI_CODEX_v2.ipynb` contains a saved-model path pattern,
`$ISLET_DATA_DIR/CODEX_scvi_model_v2/CODEX_scvi_model_v2`, but its training cells
are unexecuted in the checked-in notebook. Its code/comments also disagree on
latent dimension and batch handling, so the notebook alone cannot substitute
for the trained model and matching run record. The frozen H5ADs carry saved
scVI representations, not the model weights needed to project the six donors.
The canonical H5AD has ten-dimensional `X_scVI` means for the original cells;
in 100 sampled matching IDs these were bitwise equal to the pre-gating
`CODEX_panc_scvi_BioCov.h5ad`. A joint scVI retrain with the notebook method
would be a new fit, with new original-donor derived coordinates; it is not a
recovery of the trained encoder.

Read-only searches under `/data/panc_CODEX`, `/blue/maigan/smith6jt` and
`/orange/maigan/smith6jt`, including plausible CODEX project/model locations,
found no verified original CODEX scVI weights or full rule matrix. A 1,394-byte
Blue rule CSV under `Islet-Explorer-Senior-main/data/` remains an unexamined
candidate; its content could not be transferred to this host. The only remote
`model.pt` hits identified were for Xenium, not CODEX. Large Orange acquisition
directories could not be exhaustively listed through the available Globus
endpoint, and creating a local transfer endpoint required interactive login.
Absence across all remote storage is therefore unproven.

## Resume condition (historical)

*Withdrawn 2026-10-03: the owner removed exact recovery and scVI retraining as
prerequisites. This is the condition as written on September 29.*

First reproduce the original donor labels from the recovered numeric gates,
complete rule chain and postprocessing. Then either locate and verify the
original scVI encoder or explicitly select and validate a joint 21-donor
retrain using the notebook method. Keep original donor IDs, physical geometry,
measurements, missingness and phenotype labels fixed. A validated common
21-donor trajectory, UMAP and Leiden refit may then change derived coordinates
and clusters. Record source-cell parity separately from the existing
app-artifact parity check; that check alone does not compare the multi-gigabyte
original single-cell H5ADs.
The nine-marker typing enhancement remains deferred. Final 21-donor PubApps
instructions must wait for this integration and its review evidence.
