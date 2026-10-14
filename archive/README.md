# archive/

Files kept for reference but **not part of the running system**. The app does not
load, source, import, or serve them. Some archived research scripts can be run
explicitly to inspect historical evidence and may then import shared helpers
from `scripts/`; that does not make them current processing or acceptance steps.
Tracked moves retain their history through `git log --follow`; newly archived
untracked work acquires history when committed.

The app does not depend on these files. Preserve the historical record when a
method or delivery document still cites it.

## joint-scvi-retrain-superseded-20261003/

`retrain_joint_scvi.py` and its HiPerGator `retrain_joint_scvi_hpg.sbatch` staged
a joint 21-donor scVI retrain (September 29). Only the `prepare` step ran
(`tmp/joint-scvi-20260929/preparation.json`, status `prepared_not_trained`);
no model was trained. The 2026-10-03 owner decision replaced scVI with Harmony
on protein PCA. The script's `original_qc()` first derived the QC cutoff now
pinned as `QC_TOTAL_MFI` in `scripts/protein_phenotyping.py`, whose
`derive_original_qc_threshold()` carries that derivation forward. The script
expects to sit in `scripts/`; its repository-relative paths do not resolve from
here.

## sam-r5-withdrawn-20260917/

R5-only SAM inference, finalization, rendering, QuPath import, and mask-geometry
test. The candidate was withdrawn after repeated HPAP-078 boundary errors. See
[the archive note](sam-r5-withdrawn-20260917/README.md) before interpreting or
running these scripts.

## docs/

`HANDOFF_2026-07-25_shiny-correctness-pass.md` — the Shiny app correctness pass of
2026-07-24/25 (zero-exclusion z-scores, phenotype availability, plotly in-plot text,
Sex/Gender, CSV exports). Superseded by `docs/HANDOFF.md` (2026-07-27/28), which
reworked the same Statistics code and closed this document's one open item (the
uncommitted Phenocycler divisor files — rejected, see `63f07a8`).

Kept because it still holds the only write-up of the plotly in-plot text findings
(`hoverlabel.font`, `uniformtext`, `automargin`) and the chromote selectize gotcha.

## shiny_app/

Superseded Shiny entrypoints. The app has used a modular architecture since
Feb 2026: `app/shiny_app/app.R` is the only entrypoint, and Shiny auto-sources
`app/shiny_app/R/`. **Neither file here was ever loaded** — `app.R` contains no
`source()` calls, and Shiny does not auto-source arbitrary `app_*.R` at the app
root. They sat in the deployed directory (the shiny-server symlink target) where
a reader could easily mistake 5,400 lines of dead code for live code.

| File | Why archived |
|---|---|
| `app_original.R` | Pre-modularization monolith (5,397 lines). Last real change May 2026, and that was only the identity-scrubbing sweep. Superseded by `app.R` + `R/*.R`. |
| `app_modular.R` | Stale snapshot of the modular entrypoint, last touched Feb 2026. Had drifted 373 lines from the live `app.R`. |

## misc/

| File | Why archived |
|---|---|
| `shiny-Ex.R` | Auto-generated **Shiny package** example file (`pkgname <- "shiny"; source(examples-header.R)`) — R CMD check output, committed by accident in the Sep 2025 restructure. Nothing in this project references it. |
| `shiny-Ex.Rout` | Its companion output. Gitignored (`*.Rout`), so local-only. |
| `FileReview.txt` | One-off Phase-2 implementation summary from May 2026. Superseded by `CLAUDE.md`, `scripts/CLAUDE.md`, and the phase history in Claude's auto-memory. |

## Root

| File | Why archived |
|---|---|
| `phenotype_rules2.csv` | Superseded phenotyping rule set (Dec 2025). The live rules are `data/phenotype_rules.csv`, parsed at startup into `PHENOTYPE_RULES_CSV` (`app/shiny_app/R/00_globals.R`). |

## Deliberately NOT archived

- **`app/docs/logo.png`** — looks like a stray doc asset, but `app/shiny_app/www/logo.png`
  is a **symlink** to it, and `mod_ai_assistant_ui.R` serves it as `src = "logo.png"`.
  Moving it breaks the AI Assistant panel's logo in production.
- **`scripts/codex_phenocycler/docs/HANDOFF.md`** — active handoff for the in-flight
  CODEX Phenocycler port, with an open blocker and a pending image-review gate. Read
  it before touching that workstream's gates, floors, or pairs.
