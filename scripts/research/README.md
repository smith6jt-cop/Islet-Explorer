# All-donor candidate and offline research

**2026-09-22 implementation status:** trajectory follows cell typing and uses
cell-type composition per canonical cell in the classifier-determined endocrine
core. The nine explicitly required
markers participate through the versioned typing requirements and refinement
paths. See the
[reference/decision record](../../docs/trajectory_typing_reference_20260919.md)
and `cell_typing_requirements-v1.json`. Whole-islet marker averages are not
trajectory or Spatial-embedding features.

The authoritative input selection is `cohort-v1.csv`: 21 donors, seven ND,
seven Aab-positive and seven T1D, including HPAP-078 (analysis ID `78`, image
ID `0078`). Paths are resolved against an explicit repository root; external
regional exports have explicit absolute paths. Corrected 6539 is selected from
the recovery tree. The rejected duplicate is never selected.

The current candidate is `data/development/all-donor-20260922-trajectory-fit-support-v1/`; the preceding
dated candidates and frozen original are preserved. Raw H5ADs,
geometry, regional exports and TIFFs are immutable hardlinks; revised derived
files are independent copies. **Never edit linked inputs.** Compatibility schema 3
requires cohort, trajectory, counts/availability, marker panel and cell
drill-downs. The current runtime is schema 4 and additionally requires typed
cell partitions plus the fixed Spatial embedding and its manifest. Research
reports live outside the app data root. Schema 1/2/3 readers
remain compatible with earlier releases. See the
[integration record](../../docs/app_boundary_integration_20260918.md) for measured
capacity, exact evidence and outstanding scientific decisions.

HPAP-078's 1,431 owner-approved corrected outlines replace the rejected set.
The app uses the most recent approved all-donor geometry and exposes no
identity-review blocker. Any future geometry edit still requires rebuilding
memberships, composition, trajectory and dependent analyses. Absence of
threshold-positive hormone cells alone does not establish non-islet identity.

`integrate_reviewed_boundaries.py` stages the verified native geometry and its
dependent measurements in a fresh explicit output root. `build_cell_lookup.py
--case 78` respects native PolygonROI/GeometryROI containment and replaces no
existing cell files. `index_candidate_geometry.R` writes missing indexed
geometry only. After full hashing, `optimize_candidate_runtime.py` adds Parquet
trajectory tables and exact single-islet outline downloads, then refreshes their
manifest entries. `validate_candidate.py` checks their CSV parity, native cell
universes, classifier eligibility, typed counts, embedding and trajectory
validation. The app only reads the saved offline fit.

`revise_trajectory_reference_rank.py` creates a fresh hardlinked schema-4
revision while replacing only the fixed display rank, its model metadata and
provenance. It preserves the raw diffusion coordinate byte-for-value and
requires a new output root. Run finalization and independent validation before
serving that revision.

## Reproduction

Use the installed analysis interpreter or an environment matching the versions
in `trajectory/specification.json`. This workflow does not require Hugging Face
services, remote training or scVI.

```bash
PY=/home/smith6jt/miniconda3/envs/scvi-env/bin/python
REPO=/data/panc_CODEX/Islet-Explorer
CANDIDATE="$REPO/data/development/all-donor-next"
REPORTS="$REPO/data/reports/all-donor-research/next"
ISLET_REVIEW="$REPO/data/reviews/islet-identity-next/islet_review_decisions.csv"
GATE_REGISTRY="$REPO/data/reviews/phenotype-gates-next/gate_registry.csv"

# Use a NEW output directory for a fresh build. --resume accepts the same
# manifest and builder identity and reuses completed donor checkpoints.
"$PY" scripts/research/build_candidate.py --cohort scripts/research/cohort-v1.csv \
  --input-root "$REPO" --output-root "$CANDIDATE" --gate-registry "$GATE_REGISTRY"
"$PY" scripts/research/prepare_preview.py --cohort scripts/research/cohort-v1.csv \
  --input-root "$REPO" --output-root "$CANDIDATE"
# Supply completed identity decisions tied to this exact source geometry.
"$PY" scripts/research/trajectory.py --input-root "$CANDIDATE" \
  --output-root "$CANDIDATE/trajectory" --islet-review "$ISLET_REVIEW"
"$PY" scripts/research/run_hypotheses.py --input-root "$CANDIDATE" \
  --registry scripts/research/hypotheses-v1.json --report-root "$REPORTS"
"$PY" scripts/research/endocrine_profiles.py --input-root "$CANDIDATE" --report-root "$REPORTS"
"$PY" scripts/research/structural_coupling.py --cohort scripts/research/cohort-v1.csv \
  --input-root "$REPO" --output-root "$CANDIDATE" --report-root "$REPORTS"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$PY" scripts/research/extended_tests.py --input-root "$CANDIDATE" --report-root "$REPORTS"
"$PY" scripts/research/check_influence.py --input-root "$CANDIDATE" --report-root "$REPORTS"
"$PY" scripts/research/build_cell_lookup.py --input-root "$CANDIDATE"
RETICULATE_PYTHON="$PY" Rscript scripts/research/build_runtime.R "$REPO" "$CANDIDATE"
"$PY" scripts/research/validate_candidate.py --input-root "$CANDIDATE" --report-root "$REPORTS"
"$PY" scripts/research/finalize_candidate.py --input-root "$CANDIDATE"
python3 scripts/research/serve_preview.py --data-root "$CANDIDATE" \
  --app-root "$REPO/app/shiny_app" --python "$PY"
```

The preview can be staged before offline analysis completes. Its initial
`candidate.json` explicitly prohibits production use. The launcher chooses the
first free loopback port at 3838 or above, writes `preview/service.json`, and
starts a separate Shiny Server with its own logs and application snapshot.
`--restart` stops only the process matching that candidate's recorded server
configuration. Preserve other services. Do not pass the unfinished candidate
to production: both the R backend and deployment verifier reject schema-2/3
data without scientific acceptance.

## Revise counts and prepare review material

For a count/UI revision that preserves the existing fit, use
`stage_revision.py --input-root PREVIOUS_CANDIDATE --output-root NEW_CANDIDATE`,
then `rebuild_phenotype_counts.py --input-root PREVIOUS_CANDIDATE
--output-root NEW_CANDIDATE --review-root NEW_REVIEW_FOLDER
--rules data/phenotype_rules.csv --legacy-gates data/phenotype_gates_provenance.csv`.
Add `--gate-review COMPLETED_DECISIONS.csv` only for recorded donor-wide gate
reviews with matching evidence hashes. The registry and assignment audit are
kept outside the runtime. This is how the current count revision was prepared;
it did not make missing scientific evidence available.

`islet_identity_review.py` takes `--input-root`, a fresh `--output-root` and
`--review-project` pointing to the independent native QuPath review project.
It writes all-islet decision templates and flagged overlays. Follow with
`audit_sam_candidates.py --input-root CANDIDATE --identity-review REVIEW_FOLDER
--review-project NATIVE_REVIEW_PROJECT --output-root NEW_SHAPE_AUDIT` for
small-candidate raw-channel figures and SAM alignment that accounts for QuPath's
collinear-vertex simplification. `export_review_regions.py --input-root CANDIDATE
--output-root NEW_REGION_REVIEW` exports the exact core/peri20 analysis regions
and checks their areas against the counts.

`package_reports.py --input-root CANDIDATE
--source-report-root PREVIOUS_CANDIDATE/research --identity-review REVIEW_FOLDER
--output-root REPORT_FOLDER`
packages the prior hypothesis results into an offline index without refitting.
Keep `REPORT_FOLDER` outside the runtime. The current report package explicitly
records the new islet-identity concern. New biological decisions belong in new
candidate/report versions, never in overwritten historical outputs.

## Trajectory contract

Cell typing precedes trajectory aggregation. The primary features are
`log1p(validated base-type count / classifier-core cell count)` in the core;
separate 20-µm peri-islet counts are retained for display and sensitivity work.
The denominator is every canonical cell centroid inside the endocrine
pixel-classifier core, counted once. Disease status, clinical outcomes, section
area, whole-islet marker averages and earlier pseudotime do not enter the fit.

Classifier eligibility has exactly two requirements: no image-edge contact and
at least one canonical cell in the classifier core. Downstream typing cannot
override it. The current cohort has 11,136 eligible islets. Complete supported
typed-composition features are available for 7,493; feature and projection
support are reported separately from eligibility and are never relabeled as an
eligibility gate.

Graph training requires valid geometry, complete typed-composition features and
at least ten measured classifier-core cells. This technical support rule does
not change classifier eligibility. Smaller eligible islets remain in the
canonical tables and are Nyström-projected after fitting when their feature
vector is available. Landmarks have equal numbers per donor, capped at 200 and
the smallest donor's graph-training set. Stable ID hashes distribute draws
across within-donor section-area quartiles. Robust landmark median/IQR scaling drops constant features; a
nonconstant feature with zero IQR uses unit scale, recorded in the model.
Cosine kNN (k=15) uses an adaptive Gaussian kernel, symmetric union edges,
alpha=1 density normalization and ten nontrivial diffusion components. The
largest component receives multiscale diffusion distance from a hash-selected
two-sweep endpoint starting at its diffusion medoid. Orientation has no
biological meaning. The method follows the diffusion ordering framework of
[Haghverdi et al.](https://www.nature.com/articles/nmeth.3971); topology/root
assumptions and method sensitivity remain material, as discussed in the
[trajectory benchmark](https://www.nature.com/articles/s41587-019-0071-9).

Nyström extension fixes the training representation. Support requires all k
neighbors to belong to one component and nearest distance no greater than twice
the nearest landmark's kth-neighbor radius. Other components retain membership
without a primary coordinate. Full-reference ranks, graph, root, transformations,
landmark IDs, eigenfunctions, feature manifest, source hashes and versions are
saved. The app never fits, reranks or rescales them when filters change.

Validation includes 21 donor-held-out reconstructions, 200 donor bootstrap
reconstructions, k=10/30, cap=100 and two alternative seeds. Shared reconstructed
coordinates are compared with absolute Spearman correlation; omitted donors are
projected using only fold-training preprocessing. Bootstrap duplicate donors
are separate sampled units; no disease-order target is used. Report all
replicates, including failures. The specified instability thresholds are median
LODO correlation below .8 or supported primary coverage below .8.

For the 2026-09-22 fit-support revision, all requested reconstructions
completed. The 1,514 graph-training islets yield one landmark component, and
7,490 of 7,493 feature-supported islets receive a primary coordinate. Median
LODO absolute Spearman correlation is 0.754 and total classifier-eligible
coverage is 67.3% (7,490 of 11,136), so the ordering remains recorded as
unstable. Parameters were not changed to obtain a preferred disease ordering.

Outcome-exclusion helpers expand the marker module and complete phenotype parent
chain. **No inferential pseudotime association tests are currently published**;
the app's curves and heatmaps are descriptive. Future tests must reconstruct
these outcome-excluded features inside donor folds and donor bootstraps, fit
donor-aware smooth associations and use a separate multiplicity family.

## Phenotype and hypothesis interpretation

Nine retired labels are removed from selectors and tests, while archival labels
remain intact. The most recent typed-cell outputs are the current approved
assignments. Unsupported phenotype numerators and marker-positive calls remain
missing rather than becoming zeros. Their availability does not affect the
classifier-core denominator or classifier eligibility, and the retired
Endocrine label is never substituted.

H1 uses a different, explicitly defined operational endpoint: the independently
exported QuPath hormone-positive union (`cells_total - Triple_Neg`) and INS-any
count. It does not certify historical subtype assignments. H9 recomputes raw
hormone positivity using the literal exported thresholds, requires exact union
count agreement, at least five measured hormone-positive cells per islet and
three islets per donor, then uses all 31 marker medians. H9 performs donor-balanced
linear size adjustment and estimates donor dispersion; its bootstrap repeats
preprocessing, and technical cell-area/count sensitivity is reported.

H2/H3 compare paired peri-minus-core structural area fractions using continuous
models with donor-balanced weights, donor-clustered inference and donor bootstrap
intervals. H1 uses a beta-binomial model. Size is decomposed into within- and
between-donor effects. H5 uses an analytic compartment-uniform spatial null:
expected nerve proximity equals the capillary 20-µm-buffer fraction of valid
compartment area. Its model adjusts structural densities separately within and
between donors. This null does not preserve nerve autocorrelation; a
shape-preserving randomization remains a scientific sensitivity requirement.

H10 uses donor-level medians of assessable islets. Baseline features are age,
sex, log section area and hormone composition. The full model adds all 31 core
and peri-islet marker measurements, structural area fractions and proximity
summaries. Every inner and outer donor-held-out fold fits its own median
imputation, scaling and L2 logistic regularization (C=.01/.1/1; strongest
regularization wins a tie). Two hundred donor-label permutations repeat the
complete nested procedure. Two hundred donor bootstrap resamples also repeat
fitting; all copies of a sampled donor stay together in both fold levels.

H4/H6/H7/H8 remain unavailable with specific missing-gate/mask reasons. Their BH
slots remain among the ten; no p-values are invented. F1 is generated after the
primary results; F2/F3 await their missing primary prerequisites. F4 selects the
prespecified fallback when there is no stable added predictive information, but
the tissue-burden test awaits a validated INS-positive tissue-area measurement.
All follow-ups on this cohort are exploratory.

The primary spatial set excludes overlapping 30-µm footprints and recovered
6539's pending annotation review. H1–H3 include 10/30-µm, flagged-region and HPAP
exclusion sensitivities; additional donor-influence and HPAP checks cover
H5/H9/H10. Complete structural-null/band and segmentation sensitivities remain
part of scientific acceptance. One acquisition per donor prevents separation of
donor and acquisition effects. Aab-positive comparisons and clinical associations
are secondary, currently untested, and must retain missingness and censoring.
