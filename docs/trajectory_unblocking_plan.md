# Trajectory unblocking plan — 2026-09-17

**Superseding measurement/typing correction, 2026-09-19:** see
[the reference and decision record](trajectory_typing_reference_20260919.md).
Typing must precede the endocrine-normalized composition inputs used to fit
trajectory, and the nine newly specified markers must contribute meaningfully
to cell phenotypes. Instructions below to preserve median-intensity fitting
inputs no longer define the target implementation.

Status updated 2026-09-22: the schema-4 runtime now uses the most recent
approved boundaries and typed-cell outputs. Eligibility is determined only by
no image-edge contact and at least one canonical cell inside the classifier
core. The saved fit records 11,136 eligible islets across all 21 donors, 7,493
feature-supported rows and 2,626 supported primary coordinates. Spatial UMAP
and Leiden use the 7,493 feature-supported typed-composition rows. The ratio
scatter, pseudotime strip, graph, display controls, 31-marker typed-expression
heatmap and schema-4 islet drilldown pass browser acceptance. Marker-positive
heatmap ratios remain separately unavailable where positivity gates are absent.
Research outputs belong in a local report folder and are absent from the app.

## Resolved trajectory blockers and remaining scientific limits

Earlier checkpoints contained 3,499 marker-median eligible rows, a 54-islet
typed gate and a 1,510-islet multi-gate set. All are superseded. Current
eligibility has only the two classifier rules above, and every donor contributes
at least 137 eligible islets. Downstream typing can make a trajectory feature or
coordinate unavailable but cannot remove eligibility. The saved trajectory is
unstable by its prespecified thresholds; its root orientation, topology, donor
association and incomplete positivity gates remain substantive limits.

The table records the original causes that motivated the implementation:

| Verified condition | Consequence |
|---|---|
| `data/phenotype_gates_provenance.csv` covers 15 donors | The six additions have no entries in the gate registry used by the builder. |
| The six additions lack archival `phenotype` values | The old `existing_label & rule_match` audit could not produce their new assignments; this dependency has now been removed. |
| All 252 donor × retained-phenotype audit rows have `validated=False` | All current validated subtype counts and core endocrine denominators are unavailable. |
| Endocrine subtype rules inherit several immune-negative gates | Valid INS/GCG/SST thresholds alone do not satisfy the current full rule chain. |
| `trajectory/features.csv` already contains all 31 core marker columns | The heatmap needs control/default/layout changes, not new marker measurements. |
| The heatmap defaults to INS, GCG, SST and CD45 | Most panel markers are initially hidden. |
| Runtime schema 2 requires a `research` artifact | Removing only the Research tab would leave an unnecessary data dependency. |

The additions are 78, 6090, 6232, 6520, 6539 and 6558. Their raw intensities,
QuPath classifications and regional exports are available. Exported hormone
thresholds can help reconstruct marker-positive flags. Some immune thresholds
are region-specific: the HPAP-078 marker export contains 3,666 distinct
thresholds for each of CD3e, CD68 and CD163. Those thresholds cannot simply be
imported as one donor-wide gate.

The 12–13 donor counts in the primary research models are a separate issue:
those models use the ND-versus-T1D contrast and endpoint-specific exclusions.
They must not restrict the trajectory donor selector or its heatmap.

## 1. Make the heatmap show the complete panel

- Derive the marker list from the explicit 31-marker panel and check it against
  the saved feature columns. Eliminate both hard-coded abbreviated lists in
  `mod_trajectory_ui.R` and `mod_trajectory_server.R`.
- Present all 31 names in a visible, compact checkbox panel, selected by
  default, with Select all and Clear controls. Preserve the panel's marker
  order and give every selected marker a readable heatmap row.
- Include INS, GCG, SST and CHGA in the heatmap. Their exclusion from fitting
  the primary axis does not exclude them from descriptive display.
- Preserve the existing core-intensity measurement and fixed pseudotime bins.
  Label bin aggregation accurately: the mean of per-islet median log1p core
  intensities. Missing bins stay missing.
- Apply the selected donor and display filters to both the heatmap and its
  measurement download. Clearing selections must produce a recoverable empty
  state. Marker selection must not refit or change the coordinate scale.

Panel: CD31, CD8a, CD3e, SMA, CD4, CD34, HLADR, PDPN, panCK, ECAD,
CD163, SST, ColIV, VIM, CD20, LGALS3, B3TUBB, GCG, KRT14, GAP43,
CD35, CHGA, PGP9.5, INS, CD44, CD45, NaKATPase, BCatenin, CD68,
BActin and CK19.

## 2. Repair the phenotype pipeline, starting with the denominator

1. Build a versioned gate-and-rule evidence registry for all 21 donors from
   the explicit cohort sources. Record marker, acquisition, intensity scale,
   threshold, spatial scope, parent-rule dependencies, source and validation
   evidence. Check raw-to-log1p conversions against the original measurements.
   Existing exported thresholds establish reproducibility; they do not alone
   establish biological validity.
2. Resolve alpha, beta and delta assignments first, including their upstream
   conditions. Produce a concrete review pack per donor: intensity and
   coexpression plots, gate overlays, representative cells on corrected images,
   and reconciliation with applicable QuPath exports. Record pass, fail or
   missing evidence for each dependency. Do not broadcast donor 6450's manual
   gates to the additions or tune gates to disease labels or desired curves.
   Use the independent native QuPath review project for the six additions;
   provide exact analysis core/peri-region overlays with the gate review pack,
   distinguishing them from the upstream IsletExpanded annotations.
3. Separate application of validated retained rules from auditing archival
   assignments. Apply the reviewed rules to measured cells in every donor,
   including cells with no historical label. Save new assignment flags and
   their rule version separately; retain original labels for comparison.
4. Document double- and triple-positive handling. Compute the endocrine
   denominator as the union of validated alpha/beta/delta cell IDs, counted
   once per core. Do not map an ambiguous `Multi_Pos` class wholesale into a
   subtype or substitute the retired `Endocrine` label.
5. Validate the remaining retained immune and structural types with the same
   evidence contract. Track availability by donor and selected type, so a
   missing unrelated gate does not suppress a supported endocrine ratio.
   Correct the negative-only ductal rule through reviewed CK19/panCK evidence;
   do not certify it merely because its negative conditions pass. PDPN remains
   unresolved until its provisional evidence is addressed.

Implementation starts in `scripts/research/phenotypes.py` and
`scripts/research/build_candidate.py`; reuse applicable logic in
`scripts/derive_gates.py` with explicit cohort inputs instead of its original
15-donor paths. No reconstruction of all historical labels is required.

If a dependency needs scientific review, the deliverable is a specific gate,
its proposed interpretation and the supporting plots/cells for that decision.
A blanket request to “validate the phenotypes” is not a completed step. An
uncertified gate must remain uncertified until the evidence supports it.

## 3. Generate count artifacts and reconnect the ratio plot

- Stage a revised candidate at `data/development/all-donor-20260917/` using
  immutable existing inputs. Write new derived files; never overwrite
  hardlinked raw cells, images or previous candidate artifacts.
- Reuse validated geometry and cell membership calculations. Generate retained
  cell-type counts in cores and tissue-clipped 20-µm peri-islet regions, with
  numerator counts for core, peri or their union. Check units, compartment
  boundaries and unique cell IDs independently against source cells.
- Write `trajectory/counts.csv`, `trajectory/endocrine_counts.csv` and
  `trajectory/phenotype_availability.csv`, plus a versioned assignment/evidence
  manifest. Keep large cell-assignment audits outside the app payload.
- Read the denominator from the dedicated count artifact. Leave the saved
  marker feature matrix, fitted model, landmarks, graph, root, coordinates,
  reference percentiles and trajectory validation artifacts unchanged.
  This preserves the old fit for inspection. The subsequently discovered
  islet-identity defect requires a new offline fit once reviewed region
  corrections are available; changed geometry also requires new measurements.
- Plot `log1p(numerator / core_endocrine_union_count)`. A validated zero
  numerator plots as zero when the denominator is positive. A zero denominator
  or unavailable required count remains missing. Retain every wholly
  unmeasured object and its exclusion reason (201 in the earlier checkpoint;
  248 in the final native-cell universe).
- Make availability messages reflect the selected donor/type and actual reason.
  Export numerator, denominator, raw ratio, plotted value and availability
  reason. Preserve all retained types in the control, including immune and
  structural types, with clear status where evidence remains unresolved.

Keep “All donors” as the default. Donor selection must affect points, curves,
heatmaps, graph highlights and downloads consistently. The full graph remains
faintly visible and the axis limits stay fixed. Single-donor curves remain
descriptive. Retain image click-through and clear empty selections.

## 4. Remove Research results from the app and package offline reports

- Remove the Research tab, server registration, result downloads, representative
  result links, related app copy and the `research_results` static resource
  route. Remove `mod_research.R` from the active app. Keep ordinary trajectory
  and marker exploration downloads.
- Export existing and future hypothesis outputs under
  `data/reports/all-donor-research/20260917/`, with an index, a combined results
  table, and folders for H1–H10 and exploratory follow-ups. Include effect
  estimates, confidence intervals, p/q values, donor/islet counts, exclusions,
  relevant plots, diagnostics and representative-image references. Use PNG
  plus a suitable vector format for generated charts. Describe unavailable
  tests explicitly; do not fabricate estimates or plots for them.
- Give the offline runner an explicit report output root. Preserve the earlier
  result artifacts and keep the report folder outside app resource mappings
  and the required runtime payload. This packaging change does not itself
  require refitting the primary hypothesis models.
- Introduce runtime schema 3 for the revised artifact contract: cohort,
  fixed trajectory, required count/availability artifacts and cell drill-downs;
  research reports are not required. Keep readers/verifiers compatible with
  existing schema 1 and 2 artifacts. Update schema-specific app branches,
  manifest generation, preview configuration, deployment verification and
  tests together. Store report hashes in a separate analysis manifest.

The heatmap and report removal can be completed while gate evidence is being
assembled. Actual ratio availability depends on steps 2 and 3.

## 5. Verify actual availability and publish the revised local preview

Acceptance requires all of the following, using new evidence for this revision:

- All 21 donors are discoverable. A donor/type availability table explains
  eligible islets, validated counts, positive/zero denominators and exclusions.
  A biological zero must never be “fixed” to make every donor show a curve.
- Source-checked alpha/beta/delta counts produce real ratio points wherever
  the denominator is valid and positive. Endocrine gate review covers all 21
  donors; unresolved donors are reported as remaining work. Passing an
  all-missing availability notice is not completion of the ratio repair.
- Tests cover cells without archival labels, complete parent-rule dependencies,
  coexpression union deduplication, missing and zero counts, retired labels,
  compartment membership and independent source-count reconciliation.
- All 31 heatmap controls and rows appear by default. All/Clear, donor changes,
  empty states and downloads agree with the selected markers and source values.
- Hash comparisons show the fitted trajectory artifacts are unchanged. Tests
  confirm changing labels, ratios, donors or heatmap controls cannot refit,
  rerank or rescale the axis.
- Browser checks exercise populated ratio plots, core/peri/both values,
  representative donors including the additions, single-donor curves, downloads
  and the correct islet images. Check the complete donor availability table
  programmatically, including legitimate zero-denominator cases.
- The app starts without the report folder. Research tabs, result downloads
  and the former HTTP resource route are absent. Offline reports can be opened
  directly from their folder.

Restart only the candidate preview, normally at
`http://127.0.0.1:3839/islet-explorer/`; retain the existing service on 3838.
Record the actual URL, logs, restart command, artifact identities and new
verification results. Update the delivery/readiness documentation to
distinguish working trajectory ratios, any remaining type-specific limitations,
offline research reports and UFRC acceptance. Preserve rollback to the previous
candidate and the frozen original release.
## Superseded 2026-09-21 checkpoint

This checkpoint predates the 2026-09-22 classifier eligibility correction
recorded at the top of this document.

The post-typing composition runner has been executed against the all-donor
candidate. It writes typed core/peri counts, the endocrine denominator, all-31
marker summaries, and explicit exclusion reasons. The fixed trajectory runner
was then invoked and stopped before fitting because donors 6516, 6520, 6549,
6550, 6558, and 6563 currently have no supported eligible compositions under
the primary rules (at least ten measured core cells, valid geometry, no primary
region overlap, and supported typed core/peri measurements). This is recorded in
`trajectory_fit/validation.json` as `status: unavailable`.

No zero denominator or unsupported peri-islet measurement was substituted. The
next unblock is biological lineage review/application for the held cells and
rebuilding composition; changing the diffusion parameters would not address the
upstream missing typing. The fitted trajectory remains absent until every donor
has supported eligible compositions.
