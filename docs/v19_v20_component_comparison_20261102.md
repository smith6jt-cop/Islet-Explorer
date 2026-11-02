# V19 versus v20: component comparison — 2026-11-02

Owner: “I meant use v19 for the endocrine cells since they were done properly. compare v19 with v20 to see which has the better of each component.”

The task is to evaluate components separately, using v19 as the selected endocrine source. A wholesale rollback to v19 was the assistant's mistaken interpretation. The prepared whole-v19 release is retained as reference evidence only; no site was switched. Dev still serves v20, production remains stopped, and no mixed candidate has been built. The main-tree exception that allowed unresolved v19 identities into a new release has been withdrawn.

## Comparison contract

* Compare frozen `uniform-20donor-20261029-v19` with frozen `uniform-20donor-20261030-v20` by exact cell UUID, donor, geometry and source hashes.
* V19's matching typing is archived in `uniform-20donor-20261027/superseded/typing-v19-20261030`. Final v20 processing is authoritative on HPG. The local working `v19-qc-work` directory is a stale intermediate and must not stand in for final v20.
* A label-count increase, fewer conflicts, or a passing assertion does not independently prove better cell identity. Distinguish compliance with an explicit owner rule, quantitative consistency and image-supported accuracy.
* Preserve v19 endocrine as the owner's selection; compare its actual calls and measurements rather than presuming every discrepancy is fixed. No arbitrary new labels, positivity thresholds or geometric changes follow from selecting a version.
* Any eventual combination must assign each UUID one identity, resolve cross-component overlaps explicitly, and recompute summaries, normalization, embeddings and spatial outputs from the selected identities. Mixing already aggregated tables is not a reproducible correction.

Evidence is being assembled under `output/v19_v20_component_comparison_20261102/`. The completed source audit for the frozen roots remains in `output/v19_selection_20261102/source_audit_r2.json`.

## Verified starting point

Both frozen roots contain the same 3,199,606 eligible UUIDs and 20 donors, with donor 78 absent. Exactly 509,408 cell labels differ. The two versions retain the same endocrine-object geometry. V19 contains 7,151 unresolved cells in the Islet/EndocrineCluster/EndocrineSingle labels; v20 contains 7,800. These figures are cells, not counts of negative islets. V19 has 953 fully negative populated objects (915 singles and 38 clusters), and v20 has 948 (910 singles and 38 clusters); all 1,468 islets in each version contain positive member cells. These measurements do not establish which individual changed call is correct.

Detailed component findings and the evidence-based selection table will be added after the matched comparisons and targeted image checks finish.
