# Cell typing before trajectory: reference and decision record

Recorded 2026-09-19 from the owner's corrections and direct review of the two
supplied literature PDFs. This changes the required scientific workflow; it
does not claim that a new cell classifier or trajectory has been implemented.

## Owner requirements

Cell typing must precede trajectory calculation. The primary islet trajectory
uses cell-type abundance per core endocrine cell. Merely running the existing
median-intensity fitter after a typing job would not implement this correction.

CHGA, CD45, VIM, LGALS3, CD44, PGP9.5, GAP43, CD35 and NaKATPase must each
participate meaningfully in cell typing. Measurement-completeness checks,
unused positivity thresholds, heatmap rows and generic quality flags do not
satisfy that requirement. The machine-readable requirements are in
[cell_typing_requirements-v1.json](../scripts/research/cell_typing_requirements-v1.json).

All 31 markers remain available in the trajectory heatmap. Marker-specific
summaries must identify their cell population and denominator. No raw
fluorescence-sum-per-endocrine formula has been established by these papers or
by the historical count-based implementation.

## What the papers actually do

| Reference | Cell typing | Trajectory observation and inputs | Relevant distinction |
|---|---|---|---|
| [Barlow2025.pdf](../literature/Barlow2025.pdf), Table 2 p.6; methods pp.20-21; Figure 2 p.8 | Donor-normalized marker profiles, hierarchical clustering and annotation, checked on original images | One islet; cell-type counts separately in the core and 20-um surround, divided by core endocrine-cell count, then log-transformed; PAGA with cosine neighbors, k=15 and Leiden resolution 1 | Direct reference for our composition measurement. Table 2 deliberately assigns different marker subsets to broad and immune/endocrine typing; some markers serve downstream functional analyses. |
| [Bodenmiller2025.pdf](../literature/Bodenmiller2025.pdf), methods pp.27-29; Figure 1 p.40, legend pp.15-16 | Hierarchical annotation using PhenoGraph/Leiden and two normalizations; three-of-four agreement, otherwise ambiguous | One identified beta cell; 28 lineage/functional marker expression features, diffusion map followed by Slingshot | Supports typing before analysis of marker state within a cell population. This is not Barlow's islet composition trajectory or endocrine-count normalization. |

The Barlow file is the April 2025 version of record of *The extra-islet
pancreas supports autoimmunity in human type 1 diabetes*. Its Figure 2 also
shows expression overlays; those do not change the composition fitting inputs.
The paper used a non-T1D reference and temporarily omitted 25 islets to isolate
a path. Those choices do not override the owner's existing requirement for
outcome-independent root selection and no disease-order tuning.

The Bodenmiller file is Steenbuck et al.'s March 11, 2025 preprint, *Imaging mass
cytometry reveals early beta-cell dysfunction and changes in immune signatures
during type 1 diabetes progression in human pancreata*: 88 donors, approximately
10,000 islets and 16 million cells. It reports no supplied trajectory start,
leaf or cluster labels. That does not establish an assumption-free ordering.
Its heatmap summarizes expression within beta cells by ROI. Its separate
infiltration score uses total image-cell count, not endocrine-cell count.
Tables S10/S11 are not in the supplied PDF. The text is inconsistent about
arcsine versus arcsinh transformation; exact replication would require code.
Its beta-cell panel includes markers absent from ours, including IAPP, HLA-ABC,
ProINS, PDX1 and NKX6.1. Our HLADR is not HLA-ABC.

## Prior typing work must be reconciled, not discarded

The active all-donor rules in `data/phenotype_rules.csv` use 15 markers. None of
the nine newly specified markers is a biological dependency in that engine.

Earlier work already provides a broader foundation. The July handoff and
`scripts/codex_phenocycler/codex_lineage.py` describe CD45 immune, CHGA endocrine,
VIM mesenchymal and PGP9.5 neural interpretation, with GAP43/B3TUBB neural
refinement and CD44/LGALS3 state annotations. Actual first-pass acceptance
records exist at
`/data/codex_phenocycler/restore_validation/reviews/pair_reviews_first_pass.csv`.
The associated `screen_first_pass_valid/run_manifest.json` covers the original
15 donors only. Accepted first-pass evidence is not a finished subtype pipeline
and is not evidence for the six added donors.

The first-pass results contain 88 valid pairs out of 90; VIM remains
computationally unstable on 6505 and 6516. The 2,676,466-row saved compartment
table contains only object ID, donor ID and broad compartment, without finished
subtype assignments. These gates operate on REDSEA-corrected measurements;
they cannot be copied directly into the current raw-log1p gate registry.
Recorded first-pass method provenance and current source versions also differ.
The handoff records later user acceptance, while corresponding level-2 CSV rows
remain pending on disk. Reconcile that evidence mismatch without discarding the
recorded acceptance or claiming that the full subtype runner already exists.

The prior screen excluded NaKATPase, and its handoff treated CD44/LGALS3 states
as separate from cell typing. The owner's current requirement takes precedence:
these markers must contribute to the resulting cell phenotypes, through
lineage/subtype inference or explicit, lineage-conditioned state inference.

## Operational roles to implement and verify

These roles guide implementation; they are not newly validated thresholds.

| Marker | Required contribution to the expanded cell phenotype |
|---|---|
| CHGA | Endocrine identity evidence alongside hormone and epithelial profiles; retain discordant/hormone-low cases for explicit interpretation. |
| CD45 | Immune identity evidence before immune subtype assignment and interpretation of shared vascular/stromal markers. |
| VIM | Mesenchymal/stromal and within-lineage expression-pattern refinement, interpreted with epithelial, vascular and immune evidence. |
| LGALS3 | Marker-defined phenotypic refinement within the assigned lineage, including macrophage and epithelial contexts. |
| CD44 | Marker-defined phenotypic refinement within the assigned lineage; do not equate positivity alone with a specific activation mechanism. |
| PGP9.5 | Neural identity support jointly with B3TUBB, GAP43 and tissue/cellular context; account for endocrine staining. |
| GAP43 | Neural phenotype refinement and growth-associated marker state; lack of GAP43 alone must not reject mature neural structures. |
| CD35 | Complement-receptor phenotype interpreted with immune/B-cell and stromal context; positivity alone does not establish a dendritic lineage. |
| NaKATPase | An actual expression feature in multivariate cell-phenotype refinement, interpreted with epithelial/endocrine and other lineage markers; quality screening alone is insufficient. |

These proposed roles are informed by the supplied papers and direct human
pancreatic evidence: [PGP9.5 staining in endocrine cells](https://pmc.ncbi.nlm.nih.gov/articles/PMC4049838/),
[different PGP9.5/GAP43 neural distributions](https://pubmed.ncbi.nlm.nih.gov/7898650/),
and [Na/K-ATPase in human beta cells, Figure 5](https://www.nature.com/articles/s41467-022-34166-z).
They must be implemented and evaluated on this panel rather than assumed to
be validated merely because a marker has a published biological association.

Broadly expressed markers need not become independent lineage labels or hard
AND gates. Preserve base cell identity and explicit marker-defined states.
Do not revive retired plotted labels by using them internally as parent groups.
The owner subsequently approved an expanded denominator on 2026-09-20: the
unique validated core alpha/beta/delta and CHGA-supported hormone-unresolved
endocrine union. Only evidence-supported assignments contribute; ambiguous
objects remain explicit. This supersedes the earlier three-type definition.

## Required downstream sequence and acceptance evidence

1. Reconcile the accepted typing work with current raw measurements, marker
   names, donor/acquisition provenance and the retained-label policy. Define
   and implement meaningful use of all nine markers, including NaKATPase.
2. Export cell identity, subtype/state, ambiguity and assignment provenance.
   Verify representative cells on images. Reuse applicable accepted evidence;
   do not demand repetition of an already accepted review without a concrete
   change or unresolved issue.
3. Count typed cells in corrected cores and tissue-clipped 20-um peri-islet
   bands. Compute `log1p(cell_type_count / core_endocrine_union_count)` as
   separate core and peri fitting features. Unknown or zero denominators remain
   unavailable. This normalization does not prove that all size effects vanish.
4. Fit and validate the fixed offline ordering on these composition features.
   Keep donor balance, outcome-independent orientation, component/projection
   checks and donor-held-out uncertainty. Donor display filters never refit it.
5. Keep marker-state displays tied to identified cell populations with explicit
   measurement definitions. Preserve all 31 choices and image click-through.
   Outcome association tests must exclude their dependent features and gating
   modules from the supporting trajectory to address circularity.

Completion evidence must show each required marker entering the relevant
assignment/state algorithm, its parameters and data provenance, and measured
influence or ablation results. Inclusion in a manifest alone is insufficient;
do not manufacture label changes when a marker carries little information.

Implementation began on 2026-09-20 in the independent candidate
`data/development/all-donor-20260920`. See
[the execution record](reprocessing_execution_20260920.md). The served release
will change only after the new data and runtime pass verification. The frozen
original remains preserved.
