# Final component candidate: Unclassified review, 2026-11-03 UTC

The final v2 candidate's **7,398 Unclassified cells** (0.2312% of 3,199,606 eligible cells) pass an independent replay of the existing rules. None replaces a named non-endocrine phenotype in the immediately preceding app-v1. All **110,721 selected v19 endocrine UUID/subtype calls**, every eligible cell's Ki67 flag and all **716,451 Tissue calls** remain unchanged against app-v1.

This is rule consistency, not independent biological validation. It does not replace object-lineage, assembled-export, release or runtime checks. The producer is scientific source `67ee6c9`, HPG job `49168359`.

## Exact preceding-app comparison

Every final Unclassified UUID was already an anatomical class or technical boundary conflict in app-v1:

| Actual app-v1 phenotype | Final Unclassified cells |
|---|---:|
| Islet | 5,171 |
| EndocrineCluster | 437 |
| EndocrineSingle | 153 |
| Conflict Endocrine component boundary | 1,637 |
| Named non-endocrine phenotype | **0** |

Across all cells, 11,101 labels change from those four prior classes: 7,398 become Unclassified, 3,672 gain a supported non-endocrine label and 31 retain an explicit conflict/doublet result. No other prior label changes. The full label census and transition matrix are in the [review receipt](unclassified_component_review_20261103.json).

**Unclassified is not a count of every cell lacking a specific biological identity.** The unchanged 716,451 Tissue cells and the separately retained conflict/doublet labels are distinct. Original object class, source decision, natural current-policy decision, marker evidence and resolution reason remain recorded.

## Current evidence

The review rehashed all 20 returned identities against the completed component audit, bound that audit to its manifest, and checked the underlying REDSEA, context, coverage, calibration and registry/code hashes. It matched historical v19 decision metadata to the archived UUIDs. All current positivity, canary and float32 fold fields replay exactly for all 7,398 cells from the actual measurements; nuclear Ki67 also matches.

Recomputed compartment and structure gates find:

- 5,761 cells have no detected compartment; 1,637 win Epithelial through endocrine evidence without positive exocrine residual markers.
- None has a winning non-endocrine compartment or satisfies a current structure-qualified cell call.
- Retained reasons are 2,870 shared-marker-only, 1,412 no-marker-above-background, 1,479 structure-marker-only and 1,637 endocrine-marker-only.

Hormone evidence stays in compartment competition even when the UUID is outside the selected v19 named set. The review did not remove that evidence to force a different cell type.

Nine cells carry historical v19 immune names, but app-v1 already held seven as technical boundary conflicts and two as Islet. Their complete current gates were inspected: two no longer meet whole-cell macrophage positivity; four have endocrine canary evidence while immune defining markers are below canary; two equal-depth compartment comparisons favor Epithelial by both brightness and the existing 0.15 interior-coverage lead; one is placed by the existing plain-positivity depth rule. The historical T-cell case also fails current whole-cell CD3e positivity. Restoring old marker-only identities would introduce a rule decision unsupported by these current gates.

## The two additional cells versus the control

The completed initial-Lymphatic control has 7,396 Unclassified cells. The final DG geometry adds exactly two; none disappears from that set. Both retain positive PDPN but fail the existing requirement of more than half the cell area inside Lymphatic geometry:

| Image / cell UUID | Control overlap | Final DG overlap | Final centroid membership | Actual app-v1 phenotype |
|---|---:|---:|---|---|
| Aab_6090 / `858e1f35-875a-41ff-8ce8-8941f1a29602` | 0.592605 | 0 | outside | EndocrineCluster |
| T1D_6551 / `c8f2fcc8-eaad-4b74-99a3-9abb4958d065` | 0.581461 | 0.404962 | inside | Islet |

These are geometry-dependent losses of a Lymphatic call in the controlled comparison. Neither removes a named identity from the preceding app. The independent all-component audit reports the corresponding authorized Lymphatic transitions and no unrelated changes.

The [compact receipt](unclassified_component_review_20261103.json) contains final per-image identity hashes, calibration bindings, the full census and transitions, and hashes of the generated review artifacts. Full cell rows, all nine immune cases, representative rows, exact measurements and the reproducible read-only script remain in `output/unclassified_candidate_review_20261103/`; the control review is in `output/unclassified_control_review_20261103/`. No scientific inputs or frozen runtime code were changed by this review.
