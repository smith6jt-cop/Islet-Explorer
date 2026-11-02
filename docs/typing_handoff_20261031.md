# Typing handoff — 2026-10-31

> **Status 2026-11-02:** steps 1–6 below are done. The rest of this file is the 10-31 state.
>
> - **DG installed.** The owner chose **DG**, the despeckled detection plus the PDPN level gate. All processing
>   from the install onward ran on HiPerGator from a sealed copy; see
>   [the HPG processing record](hpg_processing_20261031.md) and memory `hpg-processing-site`.
> - **Validation.** The validator branch is `lymphatic_despeckled_level_gate/1` (`structure_generation.py`,
>   c010ad4). HPG reproduced the workstation's 10-30 outputs byte-for-byte before anything new ran.
> - **Typing.** No acceptance assertion was broken. Endothelium Lymphatic has 2,705 cells (option 1: 3,033).
> - **v20 delivered and returned to the workstation:**
>   - root `data/development/uniform-20donor-20261030-v20`;
>   - QuPath project `data/development/qupath-review-20261030-v20`;
>   - the final reconciliation passed on all 20 images.
> - **v20 on the dev site (2026-11-02 06:07 UTC):** release `pubapps-2026-11-uniform-v20` (commit c974946) runs on
>   the staging instance at `https://devislet-explorer.rc.ufl.edu/islet-explorer/`; production is unchanged and
>   stopped. Record: `docs/ufrc_readiness.md`, "Candidate c974946 with the v20 root".
> - **Step 6 answered (2026-11-02):** "They don't need despeckle"; Nerve and Capillary detection stays on raw
>   pixels (`docs/typing_rulings_ledger.md`, 2026-11-01/02).
> - **Open:**
>   - the owner's review of v20 in QuPath.

Written before a context clear. Read with `docs/typing_rulings_ledger.md` (sections "Decisions — 2026-10-30",
"Decisions — 2026-10-30/31" and "Decisions — 2026-11-01/02"), the previous handoff `docs/typing_handoff_20261030.md`, and the approved plan
`/home/smith6jt/.claude/plans/make-a-plan-to-virtual-church.md`. Task-local scripts and records live in
`output/typing_20261030/` (never committed). Run: `data/development/uniform-20donor-20261029-v19-qc-work` (`$RUN`).

## Done

1. **The 10-30 typing rulings (commit `2075d31`).**
   - **Positivity.** Every typing marker is positive on the whole cell; Ki67 is read on the nucleus.
   - **Registry.** A new column `differentiation` (= `membrane`) for CD3e, CD68, CD163 and ECAD names the reading
     that decides between groups competing for a cell.
   - **`levels.py`: pass 1.** It chooses the differentiation readings exactly as v19 did. Those readings also define
     every reference population.
   - **`levels.py`: pass 2.** A `positivity` stage chooses the whole-cell levels of those four markers:
     - CD3e: background peak among CD45⁺ cells;
     - CD68 and CD163: narrow peak;
     - E-cadherin: edge, floored at the empty-slide null.

     The isolated-hormone background check gates them. The 10-28 agent tests (CD45⁺ among CD3e⁺ ≥ 0.80;
     E-cadherin among isolated leukocytes ≤ 0.10) are reported for these levels and still gate the membrane levels.
   - **`levels.py` outputs:** `levels.csv` (positivity) and `levels_differentiation.csv`, each with a `_selection`
     file.
   - **`canary.py`.** It reads its reference cells on the differentiation readings.
   - **`identity.py`.**
     - The comparison fold (scores, brightness, conflict order, tie rule, structure top marker) uses the membrane
       reading for the four markers.
     - Structure-defined compartments (Neural; Endothelial) are never detected by their markers. Only the structure
       step names them.
     - A cell left with only such markers keeps its start class with status `structure_marker_only` and the tag
       `<compartment>_marker_outside_structure`.
   - **Labels file.** The marker-route Endothelial and Neural rows are gone.
   - **New acceptance assertions** cover both rulings.
   - **Owner choices.** Design A (markers outside the objects place nothing; other evidence types the cell). No new
     synthetic tests.
2. **Lymphatic, first resolution (owner option 1, the PDPN level gate).**
   - **Rule.** Keep the reviewed re-run objects (`$RUN/structures/Lymphatic`) whose mean despeckled PDPN reaches the
     image's own PDPN level (v19 `levels.csv`) and whose Tissue area is at least 4 µm².
   - **Evidence:** `output/typing_20261030/lymphatic_candidates/`.
     - Fragment and hole sizes alone did not separate ND_6548/T1D_6534 from the other 18 images.
     - The gate brought both into the cohort's range and kept 94 % and 88 % of their previous objects' area.
     - Every image was rendered and viewed; the renders I reviewed are byte-identical to the installed views.
   - **Selection.** `output/typing_20261030/lymphatic_pdpn_gate_v1/` (objects, per-image records, renders and reviews,
     `owner_instruction.json`, `restoration.json`, `selection/manifest.json`). It is active through
     `$RUN/structure_generation.json`.
   - **Validator.** `scripts/uniform/structure_generation.py`, branch `lymphatic_pdpn_level_gate/1` (commit `20b0285`).
     The previous binding (the exact previous Lymphatic objects) is kept at `selection/previous_binding.json`.
3. **Reruns** (driver `output/typing_20261030/run_v20.py`, the Codex agent's guards):
   - The structures phase (`runs/20261030T210113.520762Z`: cell_context, structure_metrics) and the typing phase
     (`runs/20261030T212756.342638Z`) completed, with all assertions 0.
   - Comparison with v19: `output/typing_20261030/typing_report/` (`summary.json`, `label_totals.csv`,
     `label_by_image.csv`, `transitions.csv`, `levels_changes.csv`).
   - Main counts: 508,972 labels changed. Endothelial 124,106 → 0; Neural 86,635 → 17,341, all from structures.
     Tissue → Exocrine 178,061, all newly E-cadherin-positive on the whole cell.
   - This typing still uses the option-1 gated objects, not despeckled detection.
   - **v19 typing archive:** `$RUN/superseded/typing-v19-20261030` with `SHA256SUMS`. `superseded/` is a symlink into
     `uniform-20donor-20261027/superseded`.
4. **Commits.** On 10-31 the owner asked for a logical commit series. It is `222cd5d` … `578767e`, including the Codex
   agent's v19 work. `.gitignore` now ignores `outputs/`, `.codex-artifacts/` and the testthat failure snapshots.

## In progress: "the lymphatic object detection needs the despeckle" (owner, 2026-10-31)

- **Method.** The unchanged `Lymphatics.json` classifier now sees despeckled PDPN.
  1. `despeckle_pdpn.py` writes PDPN through `redsea_cells.despeckle` (the cells' 3×3 median) as raw uint16:
     `lymphatic_despeckle/pdpn/<image>.u16` with a `.json` of hashes and check sums.
  2. `lymphatics_despeckled.groovy` is a QuPath `TransformingImageServer` that serves those values in place of the
     PDPN band. It verifies the served sums and that every other channel is unchanged.
  3. Objects are created exactly as `structures_in_tissue.groovy` does: 4 µm² minimum, holes kept, SPLIT,
     edge-touching objects kept, each clipped to the exact Tissue.
  4. The driver is `lymphatic_despeckle_run.py`.
- **Validation.** With `--mode original`, the same script reproduces the reviewed re-run objects exactly on Aab_6505
  and ND_6548 (`lymphatic_despeckle/runs/original/`).
- **Run.** `--mode despeckled`, 8 workers, started 2026-10-31 09:02 as a background task of the old session.
  - At 09:51, 3 of 20 were done: Aab_6505 1,877 objects (raw re-run 2,570), ND_6356 6,029 (7,285), ND_6516 9,791
    (11,952).
  - The other eight sat in QuPath's contour tracing at 1–2 cores each, probably bound by input/output (8 JVMs
    reading TIFFs from `/data`).
- **Check completion.** Run `ls output/typing_20261030/lymphatic_despeckle/runs/despeckled/*/result.json | wc -l`
  (20 when done). `lymphatic_despeckle/despeckled_run.log` ends with `exit=N`; a nonzero N names the failed image in
  the traceback.
- **Rerunning an image.** The Groovy script refuses to overwrite. Remove that image's folder (only one without
  `result.json`), or move it aside, then rerun `--mode despeckled --images <...> --workers 4`. Fewer workers are
  likely faster on this disk.

## Next steps

1. **Evaluate.** Run `python output/typing_20261030/lymphatic_despeckle_eval.py`, which writes `evaluation.csv` and
   `size_bins.csv`. It compares D (despeckled detection) and DG (D plus the PDPN level gate) with the previous, raw
   re-run and installed gated objects.
2. **Image check, every image** (standing rule). Render with `lymphatic_render.draw_views` and view all 20. Report
   the numbers and renders to the owner; the owner chooses D or DG. Nothing downstream runs before that.
3. **Install the choice as a new selection.**
   - The validator's `lymphatic_pdpn_level_gate/1` branch requires `source_objects` to be
     `$RUN/structures/Lymphatic` (the raw re-run). A despeckled source needs a new branch or schema that pins the
     despeckled run's `result.json`, the Groovy script, the raw PDPN hashes and the validation.
   - Adapt `activate_pdpn_gate_selection.py`. It moves the current binding aside, restores it if validation fails,
     and keeps the old binding.
   - Update `run_v20.py`, which currently requires `lymphatic_selection == 'pdpn_level_gate'`.
4. **Rerun.** `run_v20.py --phase structures`, then `--phase typing`, then `typing_report.py`. Report to the owner.
5. **Deliver, on the owner's go-ahead, with QuPath closed.** `run_v20.py --phase deliver`:
   - assembly → `data/development/uniform-20donor-20261030-v20`;
   - QuPath → `data/development/qupath-review-20261030-v20`;
   - then the final reconciliation (`output/typing_resume_20261029/final_reconciliation.py`).
6. **Ask the owner** whether the Nerve and Capillary detection also needs the despeckle. They still read raw pixels;
   only Lymphatic was named.

## Constraints

- **QuPath check before heavy steps.** Run `pgrep -x QuPath || pgrep -f '^[^ ]*/QuPath/(bin/QuPath|lib/)'`. It also
  matches our own headless QuPath runs.
- **Untouched outputs.** Production v10, the v19 root and the v19 QuPath project stay unchanged. New outputs always
  go to fresh destinations.
- **Owner rules:**
  - no interpretation;
  - "proportion", American English;
  - fix methods rather than flag images;
  - check every gating geometry against the image (whole image and six tiles, every image viewed, outliers reported
    before downstream).
