# Final typing run: handoff (2026-10-28, end of session)

## Resume prompt (paste into a new session)

> Continue the approved final typing run for the uniform 20-donor CODEX data. Read, in order:
> `docs/final_typing_handoff_20261028.md` (this file), `docs/typing_rulings_ledger.md`,
> `docs/final_typing_plan_20261028.txt`, the full plan `/home/smith6jt/.claude/plans/the-agent-sessions-have-keen-summit.md`,
> and `output/typing_final_checks_20261028/README.md`. Steps 1-4 are done; the owner decided the two Step 4 stop items
> on 2026-10-28 (below). Apply those two decisions in `scripts/uniform/typing/levels.py` and the registry, rerun
> `levels.py`, then do Steps 5-7 of the plan and stop at Step 7 with the report. Do not assemble a root, refresh QuPath
> or touch the app before the owner has seen Step 7. Ask the owner only about items no ruling in the ledger covers.

## Where things are

| What | Where |
|---|---|
| Approved plan | `/home/smith6jt/.claude/plans/the-agent-sessions-have-keen-summit.md`; owner copy `docs/final_typing_plan_20261028.txt` |
| Owner rulings (read before asking anything) | `docs/typing_rulings_ledger.md` |
| Marker registry (roles, compartments, measurements, required) | `scripts/uniform/typing/registry/codex-uniform-32__pancreas.csv`; constants `...__constants.csv`; loader `scripts/uniform/typing/registry.py`; run panel/tissue in `scripts/uniform/cohort.csv` |
| Positivity selection (Step 4) | `scripts/uniform/typing/levels.py`; outputs `data/development/uniform-20donor-20261027/typing/levels.csv`, `levels_options.csv`, `levels_selection.csv` |
| Step 2-4 report, histograms, tiles | `output/typing_final_checks_20261028/` (README.md, histograms/, tiles/); renders `scripts/research/typing_final_evidence_20261028.py` |
| Islet object check (Step 2) | `output/islet_object_check_20261028/`; `scripts/research/islet_object_check_20261028.py` |
| v17 typing (never assembled) | `data/development/uniform-20donor-20261027/superseded/typing-v17-20261028/` (README); its `positivity_levels.csv` and `shifts.json` hold the empty-slide nulls `levels.py` reads |
| Inputs of the run | `typing/redsea/` (REDSEA compartments, `Cell__`/`Membrane__`/`Nucleus__` per marker), `typing/context/`, `cells/`, `islets/`, `structures/` |
| Notice to other agents | `docs/typing_change_notice_20261027.md` (status: final run in progress; adopt nothing until root v18) |
| Commits this session | 34c7f7e, 0417264, 1cf5fe6, 7d006c0 (+ this handoff) on `spatial-review-20261020` |

## Owner decisions at the Step 4 stop (2026-10-28, verbatim: "1. use this", "2. use (a)")

1. **SMA**: positive above Barlow's known-negative 99th percentile, per image: the 99th percentile of SMA (REDSEA
   whole cell) among endocrine-object hormone-positive cells and leukocytes (CD45+ and CD3e+, CD68+ or CD163+) that have
   no SMA-positive cell (at the background-peak level) within 20 um; validated against the Capillary objects (positive
   >= 5x as often inside as in ordinary tissue). Measured 2026-10-28: SMA-positive cells 1-5 % of all cells (was 2-22 %);
   vessel enrichment 6.7-16x on every ND image, > 5x in 19 of 20 images (T1D_6520 4.4x, flag).
2. **B cells, option (a)**: B cell = CD20 above the 99.9th percentile of CD20 (REDSEA whole cell) among epithelial cells
   (E-cadherin or panCK positive, CD45-negative) more than 20 um from any leukocyte, AND CD45-positive, AND a single
   hematopoietic lineage in the canary step (no CD3e, CD68 or CD163). Report per image the number of B-cell calls, the
   proportion of CD20-positive cells that are CD45-positive, and the expected proportion of false calls (by
   construction 0.1 % of non-B cells exceed the cutoff: expected false calls = 0.001 x the image's CD45+ non-B cells).
   Measured 2026-10-28 at this cutoff: CD45+ among CD20+ 80-85 % in T1D_6520/6533/6534, 2-38 % in ND images.

## Implementation for the next session (Step 4 completion)

- `levels.py`: add a method `negative_p99` (SMA) and `negative_p999` (CD20) computing the level from the named negative
  population per image (isolation as in `isolated_rate`, ISOLATION_UM); SMA's check returns to the Capillary-object
  enrichment (non-circular: the level comes from endocrine cells and leukocytes); CD20 is chosen as decided (its check
  is the reported B-cell specificity, not a gate). Record both in the registry (`measurement_candidates`, a
  `background_method` column if useful) and the constants file with source "owner decision 2026-10-28".
- Rerun: `python scripts/uniform/typing/levels.py --out data/development/uniform-20donor-20261027 --workers 20`
  (about 20 s; OMP/OPENBLAS/MKL/NUMBA threads = 3 with 20 workers). Expected: no `stop`; CD34, CK19, PGP9.5, VIM stay
  dropped (shared or not required; they fail on an ND image).
- Commit; update `output/typing_final_checks_20261028/README.md`.

## Then Steps 5-7 (plan)

- Step 5: canary detection level (0-3 widths above each marker's level) chosen once on the 7 ND images against a null
  (each marker's positive rate in compartments that cannot carry it, cells farther than 20 um from any structure
  object); write CANARY_EXTRA_WIDTHS to the constants file.
- Step 6: canary on every cell (compartment groups from the registry; depth = canary-positive defining + shared markers
  of a detected compartment; single / resolved (top depth >= 2 and >= 2 above runner-up) / lineage_ambiguous /
  doublet_suspected; the same rule between lineages inside Epithelial and Hematopoietic); identity by the hierarchy
  (start class = parent object class; endocrine only in endocrine objects with the hormone 2x margin and CHGA+
  hormone-negative; T by CD3e with CD8a/CD4 by positivity, both -> higher, neither -> interim `T_subtype_not_shown`;
  B as decided; macrophage profile; CD45+; Exocrine; others named by top positive defining marker; ColIV never read);
  `scripts/uniform/typing/acceptance_checks.py` and `scripts/tests/test_uniform_typing.py` (new); tiles for every label
  and status, inspected by parallel reviewers.
- Step 7: stop and report (cells with no marker above background, T cells with neither subtype marker, CK19+ share of
  Exocrine, confirmation rates, transitions vs v16 and v17). No assembly before the owner has seen it.

## Determinations already made (agent, delegated by the owner 2026-10-28; in the constants file)

K = 3 widths; out-of-place rates on reference cells isolated by 20 um; leukocyte reference = CD45 + (CD3e|CD68|CD163);
options chosen on the 7 ND images, Aab+/T1D failures flagged; the most sensitive passing option kept; a width that
cannot be read or exceeds 1.5 log units, or a level > 10x from its cohort median, fails the image; pass bars:
object enrichment >= 5, CD45+ among CD3e+ >= 0.80, CD3e+ among CD68+/CD163+ <= 0.05, leukocyte markers in hormone
cells <= 0.05, epithelial markers in leukocytes <= 0.10, shared markers out of place <= 0.10.
