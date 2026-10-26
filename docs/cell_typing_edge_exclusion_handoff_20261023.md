# Cell typing and image edge exclusion implementation handoff

Owner directive recorded 23 October 2026. **Hypothesis testing is paused.** The next session should implement and validate the typing, state-calling, counting and exclusion requirements below before any H1–H10 analysis resumes. This file is a ready-to-use new-session prompt and implementation plan. No typing or exclusion changes have been applied by writing this handoff.

## Prompt for the new session

Implement the directives in `/data/panc_CODEX/Islet-Explorer/docs/cell_typing_edge_exclusion_handoff_20261023.md` in the Islet Explorer project. Hypothesis testing is paused; do not run or extend H1–H10. Prioritize marker-based cell typing, macrophage subtypes, calibrated binary state calls including proliferating, exclusive cell counting, and exclusion of cells and endocrine objects truncated by the image boundary.

Use the existing reproducible 20-donor pipeline and the accepted v10 data as the baseline. Inspect current project instructions and the relevant source code. Build a separately versioned development data root; preserve source images, segmentation records, current accepted data and the existing research outputs. Carry the changes through pipeline outputs, app summaries, drilldowns, exports, schema validation and tests, rather than implementing only an analysis-side relabeling. Produce a reviewable app/data result with QC figures and a reproducible run record. Keep donor 78 excluded. Do not deploy or switch the live data root as part of this task.

Each eligible physical cell must have one canonical identity/subtype leaf. State calls are attributes of that cell. Detailed counting categories must be a mutually exclusive identity-by-state partition, so a proliferating macrophage is counted once in its detailed category. Parent totals must be calculated from children, not added as extra categories. Missing or uncallable measurements must never become negatives. Exclude edge-truncated cells and whole edge-truncated objects using source geometry and pre-clipping provenance, not only final centroids or final clipped outlines.

Use all acquired markers relevant to identity, macrophage subtyping and measurable states, with their proper measurement compartments. Calibrate and validate the binary state thresholds using source images and technical/biological evidence. Implement positive/negative proliferation calls on callable cells; retaining only a continuous Ki67 rank does not satisfy the requested result. Do not invent a threshold, functional macrophage identity, or negative result when required evidence/measurements are unavailable. Complete independent implementation work and report precisely any calibration evidence still required. Deliver the files and acceptance evidence listed below. Leave hypothesis testing paused.

## Governing decisions and baseline

Read the latest owner decisions at the top of [CLAUDE.md](/data/panc_CODEX/Islet-Explorer/CLAUDE.md), [uniform pipeline README](/data/panc_CODEX/Islet-Explorer/scripts/uniform/README.md), and applicable subdirectory instructions. The current baseline is:

- Accepted app root: `/data/panc_CODEX/Islet-Explorer/data/development/uniform-20donor-20261022-v10/`.
- Uniform source run: `/data/panc_CODEX/Islet-Explorer/data/development/uniform-20donor-20261015/`.
- Cohort: 20 donors, 7 ND / 7 Aab+ / 6 T1D. The excluded-donor registry remains authoritative.
- Retain the approved whole-cell endocrine-object geometry contract and mutually exclusive hormone resolution, except where the new image-edge directive requires exclusion. Blood-vessel, lymphatic and nerve annotations remain structural objects, per the existing owner direction.
- All raw inputs and previous results remain available for provenance. Old H1/H1-extension results describe the earlier typing/eligibility rules and must not be treated as validation of the new rules or silently overwritten.
- Preserve unrelated working-tree changes. Do not start another hypothesis analysis, tune typing to disease separation or Beta-density associations, or automatically refit hypothesis models after rebuilding data.

The new owner direction explicitly requires macrophage subtypes, binary states and edge exclusions; older documentation describing the coarser typing is not a reason to defer these requirements.

## Implementation sequence

### 1. Inventory data and trace the current code path

Trace cells from persistent source IDs and segmentation geometry through REDSEA, alignment/scaling, lineage typing, refinement, immune subtyping, assembly and every app/export consumer. Record input paths, channels/compartments, versions, hashes, dimensions and excluded-donor checks.

Starting points:

| Component | Files to inspect |
|---|---|
| Marker compartments, rules and refinement features | `scripts/uniform/typing/codex_typing.py`, `classify.py`, `refine.py`, `immune_subtype.py` |
| Measurement, masks and geometry QC | `scripts/uniform/typing/redsea_cells.py`, `cell_context.py` |
| Object creation and whole-cell snapping | `scripts/uniform/qupath/1_islet_objects.groovy`, other islet steps, `scripts/uniform/check_islet_objects.py` |
| App labels and aggregate construction | `scripts/uniform/assemble_uniform.py`, `scripts/uniform/app_phenotype_rules.csv` |
| App counting, marker flags, palettes and exports | `app/shiny_app/R/`, H5AD/master/release builders and validators; search all phenotype and `_pos` consumers |

Known issues to verify and resolve:

- `immune_subtype.py` already resolves one winning subtype, but `assemble_uniform.py::app_labels()` flattens macrophage subtypes into `Macrophage`. A subtype implementation must survive assembly and app loading.
- Imported rules contain `Mac-M1-like`/`Mac-M2-like`. The current M1-like rule lacks required PD_L1/iNOS measurements; M2-like can reduce to an acquired CD163/LGALS3/CD68 pattern. Do not carry these names forward as established functional identities.
- The generic positivity export currently uses comparisons that can convert NaN to false/zero. Repair missingness semantics throughout the pipeline, files and UI.
- Ki67 is included in current refinement features. However, assembly currently maps `cell_class` into the served phenotype and stores `refined_class` alongside it; do not assume refinement already determines served identity. Explicitly choose and document the final authoritative label source. Separate state measurement from any authoritative lineage/refinement stage so a Ki67 state is not partly selected by that same signal. Apply the same independence check to other state markers, while retaining legitimate lineage-defining markers.
- Anonymous per-object cell CSVs omit persistent IDs. New cell/drilldown/distance exports must carry source cell IDs so category and membership reconciliation is possible.

### 2. Specify one identity hierarchy and measurable states

Create a versioned, machine-readable marker/rule registry before cohort rebuilding. For each identity, subtype and state specify required/forbidden markers, OR groups, compartment, threshold/calibration provenance, missingness behavior, conflict handling and parent/leaf mapping. Write the scientific rationale beside the registry. Check relevant primary literature/documentation when settling biological definitions.

Common immune-associated channels are **CD45, CD3e, CD8a, CD4, CD20, CD68, CD163, HLADR, CD35, LGALS3 and CD44**, with **VIM** as a shared state/context channel. **Ki67 is acquired in all 20 donors in the nucleus**, although omitted from the served 31-marker matrix. Use other measured channels for lineage exclusion/context where appropriate. epCAM/KRT8-18/Ghrelin are not available on every donor; unavailable measurements cannot be assumed negative. CD206, iNOS, PD_L1, FOXP3 and the other absent channels cannot be used to manufacture measured phenotypes.

Retain appropriate endocrine and residual identities and explicitly resolve contradictory lineage evidence. A single marker-positive event must not automatically assign a lineage: CD4 can occur outside T cells, and LGALS3/CD44/VIM/HLADR/CD35 are not lineage-exclusive. Use a documented evidence hierarchy and an explicit unresolved/conflicting leaf when evidence cannot distinguish identities. Alphabetical ordering must not decide a biological conflict.

Within confirmed macrophage identity, implement mutually exclusive **marker-defined macrophage subtypes**. A concrete starting design is the four CD163± × HLADR± combinations, with LGALS3 status retained or incorporated into exclusive refinements if calibration and support justify it. Validate the actual design against source measurements before finalizing it. Label these as measured patterns rather than inferring M1/M2 function. Include a macrophage-subtype-unresolved leaf when necessary; never count one macrophage in several subtype bins. Export legacy subtype matches as provenance only.

Identify which additional binary states the acquired panel can support and implement those with explicit rules. Avoid unsupported activation, exhaustion, Treg, senescence or cytotoxic-function labels. Proliferation is a required state implementation, addressed next.

### 3. Implement calibrated binary proliferation calls

Read the [source Ki67 audit](/data/panc_CODEX/Islet-Explorer/output/h1_immune_size_20261023/ki67/audit.md) and source-linked files in that folder. Ki67 has complete nuclear measurements for the previously curated immune set, but alignment/background agreement alone did not validate biological positivity. Extend the source join and state calling to **all eligible cell lineages**, rather than limiting proliferation to the old immune union.

Use nuclear Ki67 with DAPI/nucleus geometry and source image review. Assemble positive/negative evidence, representative donors/acquisition groups and nuclear-localization examples. Do not use epithelial background as an automatically nonproliferating control. Do not define proliferating as a fixed donor percentile, a donor-relative rank tail, or generic scaled Ki67 ≥0.10 without calibration.

Choose and document a reproducible threshold procedure on a calibration set, freeze it, and validate on held-out or blinded cell/image examples. The method may account for technical alignment/acquisition differences with recorded evidence; it must not force the same positive fraction per donor or optimize disease-group separation. Record threshold sensitivity, nuclear/mask artifacts, error cases and callable coverage. Preserve raw/aligned/scaled continuous Ki67 alongside the binary result.

Required fields include `proliferating` as nullable Boolean (or 0/1/NA), `proliferating_callable`, `proliferating_reason`, and a threshold/rule version. Callable cells must receive an actual positive or negative call. NA is reserved for a recorded missing, failed-QC or genuinely uncallable measurement; it is not a substitute for implementing the binary caller. A Ki67-positive operational proliferation-associated state must retain its assay definition and must not be described as direct proof of completed cell division. If available evidence cannot support a defensible threshold, implement the complete caller/validation workflow and report the specific missing calibration evidence without inventing a biological call or resuming hypotheses.

### 4. Enforce exclusive counting at the cell and app levels

Create one canonical cell table keyed by stable donor/image/source `object_id`, with at least:

- Unique identity leaf, parent lineage and macrophage subtype where applicable.
- Marker measured/callable flags, calibrated binary state fields, continuous values and provenance.
- Cell/nucleus geometry QC, image-edge flags, exclusion reason(s) and overall analysis eligibility.
- One canonical detailed category key: identity leaf × the complete implemented state vector, including explicit uncallable-state combinations where needed.
- Unique primary spatial compartment/object ownership for additive censuses, plus any separate many-to-many neighborhood relation table.

For example, a CD163+HLADR− macrophage with proliferating=true occupies one detailed leaf. A plot must not stack `Macrophage`, its CD163+HLADR− child, and `Proliferating` as three additive categories containing that same cell. Parent rollups are sums of their children. State marginal percentages use their own callable denominators and are explicitly non-additive across different states. The main categorical counts and percentages always form a partition.

Default tissue/donor censuses count distinct persistent cell IDs. The same peri cell near two objects cannot inflate a donor census. Where additive object/compartment totals are needed, use a deterministic exclusive owner with a documented spatial rule and tie break. Keep overlapping neighborhood membership, if needed for spatial inspection, in a separate relation and never add it as independent cell categories or biological replicates.

### 5. Exclude image-edge-truncated cells and objects

Obtain authoritative full-resolution image bounds and physical calibration from source image/OME/QuPath metadata, cross-checking axes, coordinate transforms and dimensions. Tissue annotations, nonzero stain bounds, current centroids and thumbnail extents are not substitutes for image acquisition bounds. Treat genuine missing-acquisition/stitching boundaries separately from internal processing tile boundaries.

**Cells:** test complete cell **and nucleus** polygons or segmentation masks against the acquisition border. Use a documented conservative border-contact rule and pixel tolerance, in source coordinates. A centroid safely inside the image does not make the cell complete. Record flags before any geometry clipping/rasterization that could hide the cut. Exclude flagged cells from analysis counts, type/state denominators and calibration/reference fitting; preserve their source records and reasons for QC. Technical spillover/segmentation context may retain them if the measurement algorithm requires it, but they must never re-enter analytical counts through that route.

**Objects:** exclude the whole endocrine object when its source/pre-clipping component is truncated by the acquisition boundary or contains an edge-truncated member cell. Do not clip away the incomplete part and present the remnant as a complete Islet, Cluster or Single. Propagate flags through merges, splits, morphology and whole-cell snapping. Do not silently shrink or reclassify an edge-damaged object to rescue it.

Critical current behavior: `1_islet_objects.groovy` intersects classifier polygons with an approximately 5 µm inset rectangle. A truncated final object may therefore lie away from the physical image edge. Trace the pre-intersection components and carry clipping provenance into the final objects. If original/pre-clip geometry cannot be recovered from retained intermediates, regenerate the relevant object step with provenance capture. A distance check on the existing final outline is insufficient. Explicitly flag objects altered by acquisition-border/inset clipping and exclude incomplete results.

Keep **image-edge truncation** distinct from **Tissue-outline contact** and ordinary endocrine-object boundaries. Do not drop every irregular cell or every object near a Tissue outline under the new rule. Preserve the existing whole-cell assignment contract for complete interior objects.

Maintain separate `core_complete` and `neighborhood_complete_<radius>` eligibility. An intact core can have a 20/100/250 µm neighborhood cut by the image. Require complete acquisition support for the corresponding neighborhood summary; report excluded support and correct its area denominator. Do not call a partial peri band a fully observed negative immune neighborhood. Core-valid objects can remain eligible for core-only summaries when their requested neighborhood is incomplete, with that distinction explicit in the app/export.

### 6. Rebuild the data contract and app consistently

Implement the rules and eligibility in the reproducible pipeline, then rebuild all 20 donors into a new versioned development root. Determine which upstream normalization, references, segmentation/object derivations and caches must be invalidated when excluded observations or feature sets change. Record deliberate reuse and input hashes; do not reuse stale cached summaries.

Update canonical cell tables, immune subtype/state files, object membership, paired drilldowns, H5AD matrices/metadata, master tables, donor/tissue summaries, spatial and distance tables, phenotype counts/proportions, object sizes, embeddings that depend on changed composition, legends/palettes, selectors, tooltips, downloads and release validation. Derived embeddings may be recomputed as data preparation with appropriate checks, while hypothesis fits remain paused. Avoid introducing subtype-dependent hard-coded assumptions that break the existing immune union or endocrine denominator.

Expose identity/subtype and state views clearly. All app views and downloaded files must use the same rule version, eligibility, persistent IDs and counting semantics. Add state coverage/uncallable counts and edge-exclusion audits. Preserve state columns through assembly; do not reduce the new results back to the old coarse label list. Any schema migration must update builders, readers and validators together. Build a local review result; keep the accepted live root available for rollback.

Audit hidden/retired label filters as well as hard-coded immune lists. `00_globals.R::PHENOTYPE_RETIRED` currently hides coarse T-cell/Immune labels that can contain unresolved cells. Replace this behavior with explicitly named, visible unresolved leaves whose counts reconcile with denominators. Inspect `spatial_islet_immune.R`, `mod_spatial_server.R`, `mod_inflammation_server.R`, palettes and `drilldown_helpers.R`; derive eligible immune parents/leaves from the canonical registry rather than maintaining inconsistent label lists.

## Required acceptance evidence

Write meaningful automated checks and run the existing applicable pipeline/app/release checks. At minimum demonstrate:

1. One canonical row per physical source cell; one identity leaf and one detailed category per eligible cell. Mutually exclusive detailed counts sum exactly to eligible unique cells, overall and within each donor/object/primary compartment. Parent sums reconcile with children, including unresolved and uncallable leaves.
2. Macrophage leaves sum exactly to the macrophage parent. Subtypes and binary states survive source-to-canonical-to-app-to-export round trips. Competing rule matches are audited and resolve reproducibly without duplicate assignment.
3. Missing markers, zero cytoplasm, failed nuclear QC and uncallable states remain unknown/excluded as specified; they never become negative because of a numeric comparison or Boolean cast. Callable proliferation values are strictly binary, with validated thresholds and versioned provenance.
4. Geometric tests cover all four edges/corners, a border-crossing cell with an interior centroid, a nucleus touching the border, a pre-clipped endocrine object whose final outline is interior, merges/splits, and a complete interior negative control. Final counted cells/objects contain no prohibited image-edge cases.
5. Incomplete neighborhoods cannot masquerade as complete negative neighborhoods. Area and cell support agree. Independent donor/tissue counts deduplicate shared peri memberships.
6. A before/after census reports exclusions and identity/subtype/state transitions by donor, object class, compartment and size. Report denominator changes explicitly; do not preserve old cohort/object totals by forcing labels or hiding exclusions. Donor 78 remains absent from every output and selector.
7. Image QA includes DAPI/Ki67 examples across donors, positive/negative/uncallable and near-threshold calls, macrophage subtype marker overlays, conflicting labels, and all edge-exclusion types. Include ordinary interior examples as well as difficult cases. Validation sampling/rules must not use H1 outcome associations to select favorable examples.
8. Canonical summaries match rendered app counts and downloaded counts. Applicable regression tests pass, input/output/rule hashes are saved, and clean rerun/resume behavior reproduces the result.

## Deliverables and stopping point

Deliver the implemented code/configuration and focused tests; a new development data root; a canonical data dictionary and exclusive-category registry; threshold/rule calibration and image-validation evidence; donor/object/type/state/edge QC tables and figures; app review/export checks; and a run record with exact commands, dependency versions, checksums and remaining limits. Include the migration/difference report against v10 and a short next-session summary.

Stop at a validated, reviewable typing/edge-exclusion implementation. **Do not resume H1–H10 or deploy the result automatically.** Hypothesis work remains paused until the typing and eligibility prerequisites have been completed and the user directs resumption. The previous research scripts/results stay archived as analyses of the earlier data definition.
