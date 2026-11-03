# Unclassified display hierarchy — 2026-11-03

At the owner's request, the application groups Tissue, literal Unclassified,
Doublet suspected and the former No subtype residual under **Unclassified**,
with four distinct subclasses. The corrected, tested source is
`337513e4ed2a4ac5775faf71fbf4590c559f7864`. The accepted data producer remains
`5e6580fe68af12dd18bee430988bc1ed64061526`; raw phenotype labels, identity
evidence, measurements and scientific artifacts are unchanged.

| Display subclass under Unclassified | Raw source labels | Interpretation |
|---|---|---|
| Insufficient typing evidence | `Tissue` | No permitted typing rule assigned a more specific identity. This can include positive shared markers, structure markers without a qualifying structure, or other insufficient/context-inappropriate evidence. It does not mean every marker is negative. |
| Insufficient evidence after reassessment | Literal `Unclassified` | A targeted cell reconsidered during endocrine-component reconciliation still lacked a supported identity under the existing nonendocrine rules. The label does not assert that the cell is biologically endocrine. |
| Lineage or subtype unresolved | `Hematopoietic CD45+`, `Hematopoietic CD45-`, `T cell` | The three source categories represented by No subtype in the current dataset. A broader assignment is retained without the displayed lineage or subtype detail. |
| Competing lineage/compartment signals | `Doublet suspected` | The classifier detected three or more competing compartments, or three or more lineages within a compartment, at its canary thresholds, subject to its object/structure and protected endocrine rules. |

The last subclass describes competing evidence. It does **not** establish that
exactly two physical cells were segmented together, and it is not triggered
merely by any three positive markers. Separate `Conflict …` labels remain
hidden QC categories and are not included in this family.

The parent is displayed as **Unclassified (all)** and uses the internal key
`All unclassified cells`. That distinct key prevents the parent from being
confused with the literal raw `Unclassified` label, which belongs only to the
reassessment subclass.

The [read-only whole-cohort census](../tmp/unclassified-help-20261103/actual-census.json)
counts all tissue-cell records across **20 donors and 3,199,606 cells**, without
card, donor, object, region or phenotype selections:

| Classification | ND | Aab+ | T1D | Total |
|---|---:|---:|---:|---:|
| **Unclassified (all)** | **333,948** | **234,223** | **255,776** | **823,947** |
| Insufficient typing evidence | 301,875 | 197,764 | 216,812 | 716,451 |
| Insufficient evidence after reassessment | 1,701 | 3,190 | 2,507 | 7,398 |
| Lineage or subtype unresolved | 30,061 | 32,679 | 36,039 | 98,779 |
| Competing lineage/compartment signals | 311 | 590 | 418 | 1,319 |

Thus the earlier **7,398 Unclassified cells** remain the literal-source
reassessment subclass; **823,947** is the newly consolidated parent. The larger
display total does not represent newly unclassified or retyped source cells.

These global counts differ from composition charts that pool selected object
memberships in core, peri-islet or combined regions. A cell can contribute to
more than one object's membership, so those pooled values need not equal the
whole-tissue census. The four subclasses are mutually exclusive within each
display population. Additive stacks and multi-selection lists use their leaves
once; the parent is available as a total and is never added beside its children.
Hidden conflict cells remain in the scientific all-cell denominator, so shown
proportions need not sum to 100%.

Display membership does not replace the scientific immune definition. Raw
T-cell and CD45-positive queries still expand to those exact source labels.
The unchanged **All immune cells** population contains **279,516 cells** and
continues to exclude **44,684 Hematopoietic CD45-negative cells**, even though
the display's unresolved-lineage subclass includes CD45-negative cells.
Likewise, **All T cells** continues to include CD4 T cells, CD8 T cells and
the raw unresolved `T cell` population. Moving that last population into an
Unclassified display subclass does not remove it from the biological T-cell
total or add either CD45 residual to that total.

The earlier candidate applied the projection during data preparation, which
packaged production releases bypass. Source `337513e…` corrects the shared
loading step so both packaged and development data expose the new family in
their selectors and combined core/peri metrics. Biological totals are computed
from the original labels first; **All T cells** is preserved after the raw
T-cell column is projected into its display subclass. The
[loader regression test](../app/shiny_app/tests/testthat/test-runtime-data-backend.R)
covers both loading routes and the rebuilt paired metrics, while the
[taxonomy regression test](../app/shiny_app/tests/testthat/test-unclassified-taxonomy.R)
checks raw T-cell membership and preservation of its total across projection.

The [actual-data integration receipt](../tmp/unclassified-local-review-20261103/actual-data-run2/actual-data-integration.json)
records **100 passing checks on the preceding `18df4f8…` source**. It verifies the parent and leaf counts across
all 20 donors, matching in-memory and DuckDB queries for both distance
definitions, disjoint additive stacks, and object/donor composition against
raw source populations. The Barlow input matrix is unchanged from baseline
`d7e20a6b6f6261b95fc6ec93c05cc18e9c31f001`: **4,848 rows and 54 feature columns**,
matrix SHA256 `4fa1e811b49ee506c83bd980345e35c42f6e4b5b5ae1ccec114d872835959513`.
No new graph fit was performed.
Those checks predate the production-loader correction and do not by themselves
establish that the packaged app presented the corrected selectors.

The [integration invocation](../tmp/unclassified-local-review-20261103/actual-data-run2/invocation.json)
binds source `18df4f8…` and its file hashes to a separate ephemeral container
using the existing `846962f7…` image, with the checkout and accepted data mounted
read-only and networking disabled. This was a source/data integration check,
not final native-image browser acceptance. The separate
[preceding candidate source record](../tmp/unclassified-local-review-20261103/image-18df4f8/source-identity.txt)
records that clean `18df4f8…` build input. **Final native browser acceptance of
the combined taxonomy and UI candidate is pending.** This document does not
stamp a current runtime identity or a PubApps deployment.

Evidence SHA256 values:

- Whole-cohort census: `f0aeaca560b0108ed813a0191af32d935e4b41cac24821ad233920de262ca415`.
- Actual-data integration: `6f494b169c3d7563266f4372c6028bfd3cbd4364606d863ca79e452ebe222364`.
- Integration invocation: `581ff88289f73323e532a855ca557a1c58dd52ed7051c467bdedc1ff9ee0cd59`.
- Preceding candidate source record: `ea93dfd9a49a209495ec8faf5731b3953cd641096f3c342666436ec1922b8242`.
