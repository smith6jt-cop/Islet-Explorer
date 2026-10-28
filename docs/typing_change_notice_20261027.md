# Typing change notice, 2026-10-27 (for agents using the uniform typing)

**For:** the Codex agent in this repository and the Claude agent on HiPerGator (spleen same-slide analysis).
Written by the Claude Code session doing the owner's typing review.

**Status (2026-10-28, evening): final typing run stopped at Step 7, before assembly, for the owner's decisions.**
Every typing marker now has a positivity method that passes its check on all 20 images (owner: "flag the methods you
are using so they will work with the images and markers"); the endocrine typing is kept as established (owner: "The
endocrine typing didn't need to change"). Report: `output/typing_final_checks_20261028/README.md` (last section);
rulings: `docs/typing_rulings_ledger.md`; code: `scripts/uniform/typing/` (levels.py, canary.py, identity.py,
acceptance_checks.py; registry/). The typing is in the run (`typing/identity/`) but not assembled. **Do not adopt any
root, level or label set until root v18 is assembled and this notice says so.**

## Current root and branch

- **Root:** `data/development/uniform-20donor-20261027-v16` (assembled 2026-10-28). v17 is paused, not assembled.
- **Run:** `data/development/uniform-20donor-20261027`. Earlier typing is kept under `superseded/`.
- **Branch:** `spatial-review-20261020`.
- **Served app and production:** still on v10; nothing is deployed.

## Rules the owner set today (all in `scripts/uniform/typing/codex_typing.py`)

| Topic | Rule | Constant or code |
|---|---|---|
| Cytoplasm | No cell is without cytoplasm. A cytoplasm-measured marker of a cell whose nucleus fills it is read over the whole cell (REDSEA `Cell__`). | `align_scale.read_measured` |
| T-cell subtype | A CD3e⁺ T cell takes the higher of CD8a and CD4, at any level. Ties are reviewed on the whole profile. | `T_SUBTYPE_BY_HIGHER`, `T_TIE_KEEP_MARKERS` |
| Macrophage | CD68 or CD163, with CD3e, E-cadherin and panCK absent. No M1/M2 labels. Labeled by profile: `CD68±CD163±Macrophage_LGALS3±CD4±` (12 labels). | `MACROPHAGE_PROFILE` |
| No default class | No Exocrine fallthrough and no supporting tier. Exocrine is a rule (E-cadherin or panCK). | |
| Leftover cells | A cell no rule claims is placed by its top marker at any level: CD3e T cell, CD20 B cell, CD68/CD163 macrophage, E-cadherin/panCK Exocrine. Otherwise it is labeled `<marker>+`. | `TOP_MARKER_*` |
| ColIV | A cell whose top marker is ColIV is excluded. | |
| CD45 classes | The CD45 rule's class is `CD45+` (was `Immune`). A leftover cell with CD45 as its top marker is `CD45+ (top marker)`. | in the code; in the data from v16 |
| Display names | Broad name first, with superscript ±: Endocrine Beta, T Cell CD4⁺, B Cell, Macrophage CD68⁺CD163⁻LGALS3⁺CD4⁻. Data labels keep ASCII `+`/`-`. | |

## In progress (2026-10-28): positivity on the log scale, per image

v15's level, scaled ≥ 0.10 on the linear cohort scale, is a fraction of a range set by each marker's brightest 0.1 %
of cells. It sat at 7× background for CD68 and CD20, 17–28× for CD45, CD8a and CD163, 64–170× for B3TUBB and SMA,
and 2.6–19× across images for one marker (`output/positivity_method_comparison_20261028/`). The owner's rule: a
threshold comes from each sample's own background; the published pipelines log-transform and standardize per sample.

The method now in `scripts/uniform/typing/positivity.py` (run after `align_scale.py`):

- `v = log1p(x)` per image on the aligned own-compartment REDSEA values.
- **Center:** the image's background peak (the mode over x > 0).
- **Width:** the peak's left half-width at half-maximum / 1.1774, pooled over the 20 images as the median per
  marker (REDSEA's clamp at 0 removes 25–46 % of cells in 8 images and leaves no readable peak there).
- **Level:** center + 3 widths (`codex_typing.POSITIVITY_K`); a cell is positive above it. Every marker is then
  positive at the same multiple of its own image background (CD3e 7.4×, CD8a 6.1×, CD68 8.5×, CD45 10.8×).
- **Majority markers** (E-cadherin, panCK, …): no background peak in tissue. Their background is the empty-slide
  null (cell-sized windows outside the Tissue outline, p99.9); their level comes from their own expressing
  population, floored at that null: E-cadherin's lower edge (peak − 3 right-side widths), panCK's mixture gate.
- **E-cadherin alone is the epithelial condition.** Measured on v15's populations, E-cadherin's edge keeps 97 % of
  rule-typed Exocrine cells and 12 % of rule-typed immune cells above it; the null would have rejected 31 % of the
  immune cells as "epithelial" (a first run moved 54,160 macrophages and 42,185 B cells to Exocrine on neighbor
  E-cadherin). panCK separates epithelial from immune cells poorly at any level, so it is no longer a negative in
  the immune rules nor the Exocrine rule's marker; it keeps a level for the top-marker step and the `_pos` flags.
- `typing/positivity_levels.csv` holds the levels; the rules, the macrophage profile, the top-marker ranking (in
  background widths above the level) and assembly's `<marker>_pos` fields read it.

Methods measured and rejected on these values: scimap's GMM(2) gate (splits the background), per-image peak width
(unreadable with the zeros), the mirrored lower-half spread (inflated by the subtraction tail), the empty-slide null
for immune markers (far below tissue background), another lineage's p99 as the reference (too clean or contaminated).

- **E-cadherin is measured twice (v17).** Review tiles of v16 showed acinar cells named by a top marker with
  E-cadherin at exactly 0: REDSEA subtracts a junctional marker's own signal (6.4 % of all cells are epithelial-level
  on QuPath's membrane mean yet zero after REDSEA; no other rule marker loses more than 0.3 % this way). Read
  uncompensated everywhere, the junction became neighbor signal in immune cells (17,884 rule immune cells lost). So
  `ECAD` (REDSEA) stays the immune rules' epithelial condition, and `ECAD_uncompensated` (QuPath membrane mean) joins
  it in the Exocrine rule (either measurement at its level) and stands for E-cadherin in the top-marker step
  (`codex_typing.UNCOMPENSATED_COPIES`). For the spleen lane: check junctional markers for REDSEA zeroing.
- **Top marker by the Xenium canary logic (v17; owner, 2026-10-28).** Cells no rule claims are resolved like
  notebook 02 cell 22: lineage depth (markers positive per panel), mutually exclusive groups, dominance margin 2 with
  top depth >= 2; ambiguous and suspected-doublet cells are kept and tagged (as cell 22 does), named by their top
  detected marker; a lineage names an approved class only with two of its markers including the defining one;
  otherwise, or with no lineage detected, the cell is named by its top marker (`classify.canary_step`,
  `codex_typing.CANARY_*`; the status is in `typing/calls` `resolved_by`). New labels: `CD3e+`, `CD20+`, `CD68+`,
  `CD163+`, `ECAD+`, `panCK+`. The release record carries `top_marker_method: canary`.

**Expect rule-typed cells to change in v16 and v17.** The transitions and the label-free co-expression checks are
in `output/typing_accuracy_checks_20261028/`.

## For downstream use (e.g. spleen)

- **Use one root throughout:** take labels, positivity and levels from a single root and its run, and don't mix in v10–v15.
- **Labels can change:** the label set changes again with v16 (`CD45+`, `CD45+ (top marker)`).
- **Positivity:** with v16, the marker positivity flags (`<marker>_pos`, `pos_frac`) use the per-image levels.
- **Record:** `scripts/uniform/README.md` records every rule with the owner's words.
