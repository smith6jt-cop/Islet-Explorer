# Unclassified display hierarchy — 2026-11-03

At the owner's request, the application groups Tissue, literal Unclassified,
Doublet suspected and the former No subtype residual under **Unclassified**,
with four distinct subclasses. The accessibility- and toolbar-corrected local candidate source is
`3c81075307c97356cf6ded5df0ac10e810c8b9fa`, with installed local image
`sha256:eb6f3c32779a22e4c1ca7a843bd9de5011de253ac42fb3fd4e574d7b34b2c920`.
The native standard suite passed 49/49, targeted family smoke passed and final
Composition and Spatial checks each passed 12/12, with the recorded visual
reviews complete. The [review receipt](evidence/local_unclassified_ui_review_20261103.json)
records the completed checks and their provenance. The accepted data producer remains
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
Hidden conflict cells remain in scientific all-cell totals and in distance-profile
all-cell denominators. The donor-composition chart normalizes the
included categories to 100%. Its CSV preserves original `cells`, `total_cells`
and `proportion`, and adds `proportion_selected` for the selected-category
percentage. The latter should not replace the original-denominator fraction.

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

The [production-data integration receipt](../tmp/unclassified-local-review-20261103/actual-data-run3-production-loader/actual-data-integration.json)
records **131 passing checks on source
`337513e4ed2a4ac5775faf71fbf4590c559f7864`**. It exercises the actual packaged
production loader with core artifact verification, confirms all four subclasses
and their parent across 20 donors, and verifies in-memory/DuckDB agreement for
both distance definitions and disjoint additive stacks. The biological All T
core/peri metrics match baseline `d7e20a6…` byte-for-byte, SHA256
`41d8b4f0f1ba02120625d73de8dd0f1163ee623bdc4dbb59d8c4d10c179af11d`.

The Barlow input matrix remains identical to baseline
`d7e20a6b6f6261b95fc6ec93c05cc18e9c31f001`: **4,848 rows and 54 feature columns**,
SHA256 `4fa1e811b49ee506c83bd980345e35c42f6e4b5b5ae1ccec114d872835959513`.
No new graph fit was performed. The
[integration invocation](../tmp/unclassified-local-review-20261103/actual-data-run3-production-loader/invocation.json)
used a separate ephemeral container with the source checkout and accepted data
mounted read-only and networking disabled. These source/data checks retain
their tested source identity; they are not a native-browser result for the
successor `3c81075…` candidate.

The [native family acceptance](../tmp/unclassified-native-20261103/native-ccda497-acceptance-receipt.json)
passed on predecessor source `ccda497…` / image `c7608905…`, including exported
Plot and Trajectory arithmetic, family toggles, four help tables and eight
inspected images. Its separate standard browser suite passed 48 of 49 checks;
the Spatial distance-group selector lacked an accessible label. The successor
adds an off-screen “Group by” label and changes two inline flex styles to
remove an observed blank row in the donor-composition toolbar. Classification
and count code are unchanged; the [exact source bridge](../tmp/unclassified-native-20261103/ccda497-to-3c81075-source-bridge.json)
records only those three substitutions, with 102 other app files identical.
The intermediate `d6405d6…` / `aceb89dc…` build
was never activated. The predecessor evidence retains its original identity;
the [successor smoke receipt](../tmp/unclassified-native-20261103/native-3c81075-smoke-receipt.json)
separately verifies the parent/four choices, 4,842 plotted objects, exact readable
five-row help census and no page errors. Its standard native suite passed 49/49
with zero Spatial accessibility violations, and final Composition checks passed
12/12, with three images reviewed. Final Spatial checks also passed 12/12,
including default/dense layouts at 1200/1600px, brush/reset, donor scope and
Distance controls. Direct visual review passed for four complete raw maps,
two default browser views and four Distance views. Two dense browser crops
were excluded for sticky-control capture occlusion, not application clipping;
the [card-controls receipt](evidence/card_controls_validation_20261103.json)
preserves these visual limits and the original identities of predecessor checks.

The [earlier 100-check receipt](../tmp/unclassified-local-review-20261103/actual-data-run2/actual-data-integration.json)
and [preceding candidate source record](../tmp/unclassified-local-review-20261103/image-18df4f8/source-identity.txt)
remain historical evidence for `18df4f8…`. Those checks preceded the packaged
loader correction. These prior records remain separate from the passed native
successor checks; no PubApps deployment is claimed.

Evidence SHA256 values:

- Whole-cohort census: `f0aeaca560b0108ed813a0191af32d935e4b41cac24821ad233920de262ca415`.
- Production-data integration (131 checks): `7da3f3f7e421a6737b7ded509cefa99919541341089609c68c0d3fe53dc174ab`.
- Production integration invocation: `ed116d83d448a9d862b889e93d7d61978fd57cf58e322f1f8e4f9696c5a23e40`.
- Historical integration (100 checks): `6f494b169c3d7563266f4372c6028bfd3cbd4364606d863ca79e452ebe222364`.

- Predecessor native family acceptance: `d904806c3bd543297d2c014754d613b365b6568679d01349290fa0f7d3392672`.
- Native predecessor-to-successor source bridge: `5a7df6eaf293a9309a30c89836bea1e590df5cdacf93015738861ee8a7fd6f6f`.
- Targeted native successor smoke: `8f1088ebb96990839913c1548a72483cf971e3a279940a00668b583a1abc55ee`.
