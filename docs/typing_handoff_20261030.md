# Typing handoff — 2026-10-30

Written before a context clear. Read with `docs/typing_rulings_ledger.md` (section "Decisions — 2026-10-30") and
`docs/typing_resume_20261029.md` (the Codex agent's record of the v19 rebuild).

## Where things are

- Current candidate generation: `data/development/uniform-20donor-20261029-v19-qc-work` (typing, context, levels,
  identity); assembled candidate `data/development/uniform-20donor-20261029-v19` (final reconciliation not run).
- QuPath project the owner inspects: `data/development/qupath-review-20261029-v19/project.qpproj` (3,525,977
  detections: the unchanged 3,522,872 original cells plus 3,105 reviewed additions). The 10-27 review copy is superseded.
- Baseline run (v18 era): `data/development/uniform-20donor-20261027`. Production app still serves v10.
- The Codex agent's gated drivers and records: `output/typing_resume_20261029/` (start with
  `remaining_qc_pipeline_instructions.txt` and `run_remaining_qc_pipeline.py`; freshness and review gates apply).
- Uncommitted work from the Codex agent is in the tree (`scripts/uniform/assemble_uniform.py`, `review_project.py`,
  `blind_review_project.py`, `typing/cell_context.py`, `typing/structure_metrics.py`, `structure_generation.py`,
  `scripts/tests/test_blind_review_project.py`, `test_structure_generation.py`, `docs/typing_resume_20261029.md`,
  `CLAUDE.md`, `docs/ufrc_readiness.md`). Do not commit or rewrite it without the owner.

## To implement next (owner rulings of 2026-10-30, ledger)

1. Positivity on the whole cell, membrane for differentiating.
   - `levels.py` / registry: positivity levels on the whole cell for every marker; Ki67 on the nucleus (registry row
     `Ki67` measurement and candidates). CD3e, CD68, CD163 and ECAD currently take membrane levels because their
     whole-cell options failed the epithelial check (BCatenin/ECAD positive among the immune cells it called): that
     check now describes what the membrane comparison must resolve, so whole-cell levels need methods that pass the
     remaining checks (as CD8a/CD4 needed known-negative levels on 10-28). Keep a membrane level for these four.
   - `identity.py`: P and C (canary) from whole-cell values; the fold used to rank and compare competing groups
     (`groups_matrix` score and brightness, conflicts, the interior-coverage tie rule) from the membrane reading for
     the four markers. Confirm the tie rule's interior coverage stays whole-cell pixels against the whole-cell level.
2. Neural and Endothelial only from structure objects.
   - Cells get Neural only through the Nerve structure step and Endothelium Blood Vessel only through the Capillary
     (blood vessel) step; outside those objects B3TUBB/GAP43 and CD31 must not name a cell (no marker-route
     "Neural" or "Endothelial" label). The overlap stays a feature (`crossed_frac_nerve`, `crossed_frac_vessel`).
   - Design to settle in code: outside the objects, treat those markers as non-detecting (like shared markers) so the
     cell is placed by its other evidence (the example 494269d7 in Aab_6450 is also ECAD/BCatenin positive and would
     most likely become Exocrine), or keep the start class with a status. Say which, with counts, before rerunning.
   - Check the labels file rows (Endothelial lineage, Neural lineage) and `name_lineage`; add synthetic test cases
     (CD31-positive cell outside a vessel; B3TUBB-positive cell with partial nerve overlap).
3. Then rerun on the v19 generation: levels, canary, coverage, identity, acceptance checks, display levels, assembly,
   QuPath import. Use the Codex agent's gated drivers or reproduce their freshness checks.

## Constraints

- The owner is working in QuPath: no heavy computation and no QuPath import until the owner allows it. Check QuPath
  with `pgrep -x QuPath || pgrep -f '^[^ ]*/QuPath/(bin/QuPath|lib/)'` (a `ps | grep` pipeline is broken here).
- Lymphatic work is paused (owner, to the Codex agent): previous Lymphatic objects are used.
- Collagen IV is all tissue collagen (mesh included), structural area and distance only; it names no cell.
- Verified this session: object membership is by centroid in v18 and v19 (no fault). The blood vessel (Capillary)
  objects come from the original `Capillary.json` pixel classifier on CD31 alone (ANN_MLP; the five additions:
  CD31 threshold fitted to it); Nerve from B3TUBB and INS; Lymphatic from PDPN.
- Owner rules in memory still apply: no biological interpretation, "proportion" not "share", American English,
  fix methods instead of flagging images, endocrine objects settled.
