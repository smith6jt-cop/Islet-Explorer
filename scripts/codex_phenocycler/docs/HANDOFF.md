# CODEX Phenocycler port — handoff (2026-07-23)

Branch `claude/spatial-cardc-per-cell-distance`, 8 commits ahead of `7571b6d`.
Plan: `~/.claude/plans/review-the-last-session-s-sharded-acorn.md`.

---

## Goal

CODEX Islet-Explorer replaces its scimap rules-gating phenotyping with the
`smith6jt-cop/Phenocycler_Analysis` method (REDSEA → RESTORE → ordered-residual lineage).
Labels **overwrite** the existing 21. The public PubApps launch is **gated** on it landing.

Upstream clone: `~/IO60panc2nd/Islet-Explorer-Senior/Phenocycler_Analysis`,
branch **`redsea-restore-curated-5compartment`** (NOT `main` — `main` is ~4 commits behind and
does not have the CD20 work). Consumed via `PHENOCYCLER_ROOT`; not yet a submodule.

---

## THE BLOCKER — read this first

**The agreed v0 no longer exists.** The 2026-07-22 decision was to port the parked
`restore.py` (SSC + mean+3σ) as a CODEX v0 and launch on it. Upstream `6d49ae7` **deleted**
`MARKER_PAIRS`, `FUNCTIONAL_PAIRS`, `run_restore` and `run_restore_functional`, because ten of
those pairs referenced CD20 and there is no mechanical substitute. It is enforced:

- `tests/test_restore.py:71` → `assert not hasattr(restore, "run_restore")`
- `pipeline.py:123` → the `restore` stage **raises**, so a run stops rather than "skipping
  silently to lineage and typing cells from gated parquet that no accepted method produced".

**Do not write a CODEX `run_restore` over the retained machinery.** That was recommended once
and withdrawn — it recreates for CODEX exactly the known-wrong web just deleted for Senior.

Retained and explicitly untouched: `idx_select`, `neg_stat`, `fit_clusters`,
`build_restore_input`, `apply_thresholds`, the threshold LUT, `run_restore_proliferation`.

### Open decision (needs the user)

1. **Follow upstream** — launch gate becomes Gate 4 + the human pair review. Correct, but the
   date is set by a review round that has not started (0 VALID pairs, 543 `PAIR_REVIEW_PENDING`,
   zero reviewer edits as of the last audit).
2. **Run CODEX's own pair review** — review scales with *pairs*, not donors × pairs, and CODEX
   has **5** first-pass pairs against Senior's 27. CODEX could reach an accepted divisor before
   Senior does. **Unverified prerequisite:** whether `restore_validation.py` can consume CODEX
   input without the per-compartment QuPath means it currently expects. Check this first.

---

## Done and verified

| Stage | State | Location |
|---|---|---|
| Stage 1 cells parquet | **DONE** 2,676,466 cells × 15 donors, 404 MB, 19s | `/data/codex_phenocycler/cells/donor_id=*/` |
| Per-cell GeoJSON | **DONE** 15/15, 4.3 GB, all verified | `/data/codex_phenocycler/redsea_scratch/geojson/` |
| REDSEA | **DONE** all 15, 843 MB, ~110s/donor (8-way) | `/data/codex_phenocycler/cells_redsea/` |
| Input floors | **DONE** measured post-REDSEA | `codex_lineage.OBSERVED_FLOORS` |
| Compartment gates + pairs | **DONE**, validators enforce | `codex_lineage.py` |
| RESTORE | **BLOCKED** — see above | — |
| Lineage / app integration | not started | — |

Everything above is upstream of the retired entry point, so none of it is invalidated.
`object_id` joins raw↔corrected **totally** (2,676,466 both sides, all unique).

### Modules (`scripts/codex_phenocycler/`)

- `codex_panel.py` — marker/channel resolution. 11 tests run against the **real** OME-XML
  sidecars and QuPath TSV header, not fixtures.
- `codex_lineage.py` — compartment gates, RESTORE pairs, floors, immune sub-split. All
  validators run at import.
- `build_cells_parquet.py` — QuPath TSV → hive parquet (`--dry-run`, `--min-cell-area`).
- `export_cells_geojson.{groovy,sh}` — per-cell GeoJSON from the QuPath project.
- `verify_geojson.py` — the four silent-failure checks.
- `codex_redsea.py` — REDSEA driver (`--check` resolves all 15 without reading pixels).

---

## Verified facts — do not relitigate

**Inputs (all present; an earlier session wrongly concluded REDSEA was impossible by searching
only `/home` at shallow depth — search `/data`).**
- Raw linear MFI: `/home/smith6jt/panc_CODEX/results/old_data_plots/all_data.tsv` — 1.39 GB,
  2,676,466 rows all `Object type = Cell`, 15 donors. **RESTORE needs raw linear; never feed it
  `layers['log']`** — no transform is applied and `_norm`/`_log2r` are ratio-scale.
- QuPath project: `/data/QuPath_Projects/QuPath_project_IsletExplorer/` (15 entries, 12 GB).
- Images: `app/shiny_app/www/local_images/*.ome.tiff` is the **source of truth** over
  `/data/.../Images/`, which has per-donor marker-name variants, a HLADR/PDPN order swap in its
  11 plain `.tif` entries, and labels the KRT14 slot `TRAM2.tif` on 6521/6548/6551. Both copies
  verified dimension-identical across all 15 donors.
- Disk: `/data` has 7.2 TB free. `/home` is 97% full — never target it.

**Panel.**
- Channel count is **not uniform**: 8 donors 35-ch, 7 donors 33-ch. Extras (`epCAM`, `KRT8-18`)
  are appended, so 35-ch is a strict prefix extension. **Never key a marker by index.**
- `epCAM`/`KRT8-18` cannot be used cohort-wide (8/15 donors).
- `Ki67` is the **only** marker QuPath exported as `Median`, not `Mean`.
- CD99 and MPO are **absent** — Senior's Endocrine-bright-CD99 and neutrophil branches are
  unreachable, as are DC/Plasma/NK (no CD11c/CD209/CD38/CD56/CD57).
- `redsea.qptiff_channels()` **cannot** read CODEX names (Akoya `<Biomarker>` *elements* vs
  Bio-Formats `Name` *attribute*) — it returns `ch0..chN` silently. `codex_redsea.py` patches it.

**Method.**
- The `>50` floor is **abandoned upstream**, replaced by donor-local marker-specific raw-linear
  floors (Otsu; `max(Otsu, Triangle)` for CD3e/CD20). It was inert on CODEX anyway — 99.2% of
  cells exceed it on ECAD.
- Floors must be measured **post-REDSEA**: correction lowers every one (CHGA −40%, ECAD −27%,
  panCK −21.7%). Method assignment is unchanged; **CHGA is CODEX's only `max(Otsu, Triangle)`
  marker**.
- **CD20 is never a reference** (`INVALID_REFERENCE_MARKERS`), and never a compartment anchor.
- **CD3e cannot substitute for CD4/CD8a** — they are T-cell co-receptors. Upstream reached the
  same conclusion independently. CODEX uses `CD4 ← CD31`, `CD8a ← CD68`.
- `epCAM` was specified as the PGP9.5 counterpart but is unavailable on 7/15 donors; **panCK**
  substitutes and measured *cleaner* on both metrics.

---

## Open issues, in priority order

1. **The RESTORE path decision** (above).
2. **CD8a and HLADR prevalence are implausible** — CD8a 17.66% cohort with a 0.59–59.80%
   per-donor range; HLADR 36.54% (0.81–62.83%). The 2-component Otsu fit is landing *inside the
   background* on some donors. **CD8a anchors the T-subtype split, so this blocks any CD4/CD8
   call** independent of the counterpart question.
3. **REDSEA zero-clamping is ~2× Senior's and hits two gates** — VIM 39.7% cohort median
   (54.6% on 6534), ECAD 31.1%; Senior saw CD31 25.7%, VIM 18.1%. Subtract-only cannot create
   signal, so these were neighbour bleed — but VIM is the last gate before `Other`. **Compare
   `Other%` and the Mesenchymal fraction pre/post-REDSEA before accepting the first pass**; both
   parquets exist, so this is cheap.
4. **App integration** (plan Phase 3) — `cell_type` must be the primary column;
   compartment-only deletes Beta/Alpha/Delta and with it the insulin-loss thesis. Blast radius
   is in the plan file. Note `ALPHA_STRUCT = 0.15` tuning is invalidated by any label change.
5. **Deployment (Phase 0)** — code fixes landed in `1159ebe` but the container has **never been
   built** (no Apptainer here; needs a HiPerGator node). User must also rotate the UF Navigator
   key in `app/shiny_app/.Renviron` and fill in `UFIT_PubApps_Instance_Request.docx`.

---

## Reproduce

```bash
conda activate scvi-env
cd ~/panc_CODEX/Islet-Explorer/scripts/codex_phenocycler
python -m pytest test_codex_panel.py -q          # 11 pass
python build_cells_parquet.py --dry-run          # schema + per-donor counts
python codex_redsea.py --check                   # resolves all 15, reads no pixels
python verify_geojson.py                         # 15/15
# upstream: 132 tests pass under scvi-env (no phenocycler_analysis env on this machine)
```

Gotchas: import `phenocycler` **before** numpy (BLAS pinning avoids a documented 15 h
`futex_wait` deadlock); pass top-level functions to `map_donors`, never lambdas (they don't
pickle); QuPath 0.7 rejects bare `-D` on `script` (use `--args`) and **exits 0 even when the
script errors**, so drive exit status off the output file.
