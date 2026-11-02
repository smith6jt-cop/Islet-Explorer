# CODEX fixes for the Senior, monorepo, HPG and MultiPlex5000 agents

This supersedes the current-status portion of the 10-27 notice. The owner requested
delivery of the recent fixes after the component comparison. The exchange follows
Senior's `docs/HANDOFF_2026-10-28_SENIOR_CODEX_COORDINATION.md`: a committed notice
and a named immutable package, with the Senior agent owning any port into
`core/multiplex5000`. Neither active Senior checkout nor its scientific data is
overwritten by this delivery.

## What is ready

The **Lymphatic isolation correction** is implemented in CODEX commit `a425c05`.
Changing Lymphatic geometry had changed the shared ordinary-tissue/structure
reference populations, recalibrating endocrine and unrelated markers. A controlled
initial-v20→DG comparison, with the same typing code and policy, found 205 changed
positivity levels across 15 markers, 265 changed canary thresholds, 715 changed
labels, and 303 changed labels outside both Lymphatic footprints.

The corrected all-20 replay passed its automatically invoked invariant audit:
3,199,606 eligible UUIDs; all marker positivity/canary/fold fields, Ki67,
REDSEA, calibration and non-Lymphatic context unchanged; all 113,134 named
endocrine decisions of that **initial-v20 diagnostic baseline** preserved.
Only 352 label changes remain, all involving Lymphatic and changed Lymphatic
geometry. The unchanged-context control reproduces its baseline Parquet exactly.
See `docs/evidence/lymphatic_isolation_20261102.json` and
`docs/v19_v20_component_comparison_20261102.md` for exact evidence paths/hashes.

Shared corrections to review for MultiPlex5000:

| CODEX source | Behavior to carry across with the recipient's profile |
|---|---|
| `typing/calibration_context.py` | Copy and hash an explicitly selected complete reference context; validate UUID, coordinate, QC and endocrine membership; never derive calibration silently from a changed structure. |
| `typing/structure_update.py` | Bind a fresh structure revision to calibration, registry, raw measurement and baseline identity hashes. Refuse rebinding or existing identity outputs. |
| `typing/levels.py`, `typing/canary.py` | Read pinned reference populations; refuse recalibration of a bound structure-only run. |
| `typing/cell_context.py` | Require the baseline before regenerating selected structure membership. Live geometric context remains distinct from calibration context. |
| `typing/identity.py` | Keep baseline named endocrine cells safe from structure overwrite, including islet members; reject added/lost/retyped endocrine calls rather than covering inconsistent evidence with copied labels. Automatically run the full component audit. |
| `audit_structure_update.py` | Require fixed measurements, QC, all markers, non-target context and endocrine decisions. Every changed label must involve Lymphatic and lie inside its old/new footprint. Unknown fields remain immutable. |

These modules use CODEX registry/schema conventions. They are source evidence for
an adapted, tested shared implementation, not drop-in replacements for Senior's
typing producers. Map marker roles through `core/panels`. Do not export CODEX
numeric thresholds, endocrine membership, cohort exclusions or labels into
PhenoCycler, Xenium, pLN or spleen as shared defaults.

## Other recent fixes included in the source history

* `34b6853`: include `Hematopoietic CD45+` in the app immune-gradient definition;
  the frozen v20 omission affected 50,625 cells. Review the equivalent label mapping
  in the receiving app.
* `6cee485`: reject unnamed endocrine definitions and upstream/UI aggregate
  disagreement at release construction. Do not count generic object labels as
  endocrine to make totals agree.
* `8035134`, `da3a153`, `62ab0af`, `5629465`: trace original classifier pixels
  through cell geometry, native/filter/REDSEA measurements and per-image levels;
  correct the diagnostic QuPath hole-fill replay and grid sampling, with complete
  input hashes and actual image panels. These are diagnostic corrections, not new
  scientific labels.
* `9127f52`: reject final endocrine objects with no positive eligible member,
  empty membership, or invalid single-cell membership after level selection.
* `f854856`: keep per-run calibration bindings/snapshots, audit outputs and transfer
  archives out of source control; preserve their receipts with each run.
* `8be9391`: render app typing descriptions from the recorded identity generation
  and measurement metadata instead of showing its stale legacy rules CSV. Real
  v19/v20 descriptions distinguish their methods; historical v10 HTML is unchanged.
  This updates explanatory text only, not scientific identities or calculations.

The temporary whole-v19 release exception at `190addc` was withdrawn by `82afb31`.
Do not port that exception or interpret its built image as the selected release.
The current source archive reflects the net result, including its removal.

## Scientific data status — do not infer acceptance from the fix

The owner selected **v19 endocrine cells**, then asked for the better component of
v19/v20. The isolation replay uses initial-v20 only to establish a controlled
counterfactual; it is not the selected combined dataset. Common Tissue, cell,
endocrine, Nerve, Capillary and Collagen IV geometries match. V20 implements the
later whole-cell/membrane separation and structure-only rules; DG is the instructed
Lymphatic method. Counts alone do not establish better cell accuracy.

No corrected mixed candidate has been assembled or accepted. V19 itself retains
7,151 unnamed endocrine-object cell labels and 953 populated objects with no
positive eligible member under its corrected-cell levels; v20 has 7,800 such
cell labels and 948 such objects. All 1,468 islets in both versions have positive
member cells. The original classifier masks replay correctly: positive pixels
and a positive corrected whole-cell mean are different tests. Do not invent
subtypes, lower thresholds to meet expected counts, or suppress those discrepancies.
Both frozen app roots also omit the upstream Ki67 state in examined cell exports.

Dev last served v20 (`c974946`); dev and production are now stopped and masked
under the owner's parallel security-remediation direction. No restart is part of this handoff. No new app root,
scientific acceptance or deployment is delivered by this notice. Fresh
authoritative processing belongs on HPG and must rebuild all dependent outputs
from one coherent identity set. The four inode-pinned measurement/refresh modules
have not been edited.

## Delivery and recipient responsibilities

The immutable source/evidence packages are placed at:

* workstation: `/home/smith6jt/IO60panc2nd/CODEX-handoff-20261102-lymphatic-isolation`
* HPG: `/blue/maigan/smith6jt/CODEX-handoff-20261102-lymphatic-isolation`

The package contains the committed source snapshot, patch history, concise evidence,
this notice and the complete workflow guide, with `SHA256SUMS` and a source-commit
record. `docs/evidence/agent_delivery_20261102.json` records the verified delivery;
file delivery is distinct from an agent reading it or completing a port.

Verified recipients/contexts at discovery:

| Workspace | State |
|---|---|
| `/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior` | `fix/portability-blockers-20261102` at `610a6a1`; actively changing. Its `core/multiplex5000/` is the shared package. |
| HPG `/blue/maigan/smith6jt/Islet-Explorer-Senior-main` | `main` at `40fde07`; dirty source and untracked handoffs. |
| HPG `/blue/maigan/smith6jt/IES-typing` | `typing/per-image-method` at `a9549cc`; clean at discovery, current typing work. |
| Local `Islet-Explorer-Senior-v3` | Older, dirty reference checkout; preserved. |

Senior should review the shared dependency fix, adapt it in its own branch, pin
its own assay/profile reference generation, and replay real data against exact
UUIDs. The invariant audit verifies change scope, not biological acceptance. CODEX
retains ownership of the current CODEX candidate and app; HPG jobs/data are not
changed by the sibling handoff package. The CODEX HPG source tree is synchronized
separately using its explicit changed-file rule and a recorded checksum receipt.

Start with [the workflow user guide](workflow_user_guide.md), then the component
comparison and the evidence receipt. Retain this notice when replacing the old
coordination status so another session does not resume from the superseded v18 wait.
