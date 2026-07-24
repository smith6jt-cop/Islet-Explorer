# CODEX Phenocycler port — handoff (2026-07-23)

Branch `claude/spatial-cardc-per-cell-distance`, 8 commits ahead of `7571b6d`.
Plans: `~/.claude/plans/continue-from-handoff-parsed-kay.md` (Phase 0, DONE below);
`~/.claude/plans/make-a-plan-to-eager-phoenix.md` (the DONE/verified table, still valid).

---

## 2026-07-24 — Phase 0 done + level-2 spec received (READ FIRST)

**Phase 0 (submodule + v8) — DONE, verified.** Upstream is now a git submodule at
`external/Phenocycler_Analysis` pinned at `8db37af` (METHOD **v8**); all 6 consumers resolve it
(`PHENOCYCLER_ROOT` override → else the submodule); pins bumped + audited (v7→v8 is additive /
contract-safe — see the Pin note near the bottom). Verified: 18 tests pass, `codex_floors.py --check`
green, `--smoke 6516 6534` load-gate passes against `external/`. Staged but NOT committed
(`.gitmodules` + the gitlink); source edits unstaged. **v8 adds `phenocycler/reference_selection.py`**
(Gate-2 `reference≥target frequency` rule + `ACCEPTED_REFERENCE_OVERRIDES` + a frozen `ACCEPTED_PAIRS`),
which is the machinery the level-2 sub-typing consumes.

**Level-2 sub-typing spec (authoritative, from the user 2026-07-24) — supersedes my earlier ad-hoc
pair guesses.** RESTORE thresholds are **global per scene** (whole-cell), one per target; the
compartment hierarchy is an *interpretation* layer applied afterward. Counterpart rules:
epithelial-family→CD31, endocrine-confounded→epCAM, hormone→other hormone, else→ECAD (LGALS3←CD31 as
the stated exception). 24 level-2 targets appended to the 5 level-1 compartment targets. Dual-use
markers (PDPN, CD34, LGALS3) get ONE threshold each, disambiguated by parent node. Full spec + node
resolutions + runtime caveats are in the session transcript; key deltas from the current code:
- **Immune pairs move to ECAD counterparts** (CD3e/CD68/CD163/CD4/CD8a/HLADR ← ECAD), matching v8's
  `ACCEPTED_PAIRS`. This **retires the current `IMMUNE_PAIRS`** (CD4←CD31, CD8a←CD68, …) and the 79/90
  immune screen built on them. CD4 read only within CD3e+; CD163 nested under CD68 — at interpretation.
- **Endocrine β/α/δ = RESTORE pairs** INS←GCG, GCG←INS, SST←INS (+ an `Endocrine-hormone-neg` class
  and a min-margin guard so all-hormone-low insulin-loss cells are never argmax'd to Beta). Acinar has
  NO positive marker → ECAD+/CHGA−/CK19− residual.
- **epCAM/KRT8-18 are 8/15-donor only** → DECIDED: substitute **panCK** for the neural counterparts
  (GAP43←panCK, B3TUBB←panCK) via `ACCEPTED_REFERENCE_OVERRIDES`, and **drop** epCAM/KRT8-18 as
  confirmatory targets. (Including epCAM as a screened marker would silently drop the 7 33-ch donors
  via the not-null filter.)
- **CD20 is re-included as a B target** (`CD20←ECAD`) — reverses the pinned both-role exclusion. This
  is a CODEX-side unwind (`RESTORE_EXCLUDED_MARKERS` + the two validators); upstream does NOT block it
  at runtime (its exclusion is import-time over its own `CANDIDATE_PAIRS`, which CODEX overrides).
- **OPEN (unresolved):** prevalence-collapse handling for {CD20, KRT14, CD35, SST, GAP43} — proposed
  per-marker split: borrow-strength-from-accepted-divisors for the ubiquitous {SST, GAP43, CD20} vs
  abstain-per-donor for the conditional {KRT14, CD35}. Resolve before running the full screen.

**`reference_selection.py` (v8) reviewed — how it feeds the level-2 design:**
- It is a Gate-2 **review artifact**, NOT a production step: reads a screen's `pair_evaluations.csv` →
  writes a per-`(target,reference)` decision matrix + 3-panel per-target figures + `frozen_reference_matrix.csv`.
  No divisor / no call. Rule = `reference_n/target_n >= min_reference_to_target_ratio` (default **1.0**,
  pure proportion, on the **median** across donors). Precedence: a frozen `overrides` entry wins outright
  (can accept a failing reference, caveat auto-recorded) → keep passing baseline → highest-ratio passing → `open`.
- **API for CODEX:** `build_reference_matrix`/`build_bundle` take `overrides=` (pass CODEX's dict — no
  monkeypatch needed for overrides). BUT the module binds `CANDIDATE_PAIRS`/`BASELINE_PAIRS` at import, so a
  CODEX adapter must **rebind `reference_selection.CANDIDATE_PAIRS`/`.BASELINE_PAIRS`** after
  `crv.patch_upstream()` (same import-namespace trap the `rpr` adapter already handles).
- **Finding A:** the existing v7 screens LACK `reference_to_target_ratio`/`reference_undersized`
  (build_reference_matrix raises without them) → the level-2 screen MUST run at v8 (already planned).
- **Finding B (non-obvious):** the endocrine hormone→hormone refs read as `undersized` in exactly the
  states of interest — `INS←GCG` is undersized in ND (alpha < beta), `GCG←INS`/`SST←INS` in T1D (beta-loss).
  This is the `E_cadherin←CD31` precedent: accept via `ACCEPTED_REFERENCE_OVERRIDES` with the caveat
  auto-recorded. Consider adding ONE exocrine comparator (panCK/ECAD) for the 3 hormones so the dossier
  has something to weigh where undersizing bites.
- **Bearing on the OPEN question:** reference_selection handles REFERENCE sparsity, NOT TARGET
  positive-mode collapse — so the prevalence-collapse decision stands unchanged. But its called-fraction
  + divisor-vs-consensus panels are the per-donor evidence for that acceptance call at review time.

Next (flow, once the OPEN decision lands): design the level-2 pair web (predeclare one reference per
target per the counterpart rules + CODEX `overrides`), register floors for the new markers, run ONE
combined screen **at v8** → the `reference_selection` dossier (pick/record references) → the
sparse-target acceptance policy → ONE image-review bundle (`codex_restore_pair_review`, extend its
shortlist) → `codex_positivity` for level-2 → `codex_lineage_run.py` (the cell_type interpretation
tree) → app integration.

---

## 2026-07-23 continuation — decisions + corrections (READ THIS OVER THE OLD BODY BELOW)

A 16-agent investigation + 3 adversarial verifiers re-measured every load-bearing claim. Net:

**Decisions taken:** (1) launch gate = **CODEX runs its own pair review**, not Senior's Gate 4;
(2) run the **QuPath per-compartment intensity pass** so the NNMF is genuine — Cell-only NNMF is
provably degenerate (n×2 matrix, rank-2 NMF = exact factorization, H = identity, labels = trivial
argmax at 99.99%); (3) **port the full `restore_pair_review.py` bundle**; (4) **do not reopen the
gate marker set** — judge composition on RESTORE output, not floors.

**Prerequisite ANSWERED (was "unverified"):** upstream v6 `evaluate_locked_cohort` runs on CODEX
today and emits a divisor on ACCEPT (proven: donor 6516 ECAD←CD31, divisor 9023.12). The human
pair review is the last gate, not a separate stage.

**Three recorded facts were WRONG — now corrected in `codex_lineage.py`:**
1. **Floors are computed on RAW, not REDSEA.** `pair_input_qc` floors on `raw__{m}`, and
   `load_validation_sample` allclose-asserts `raw__{m}` == uncorrected QuPath means. So
   `OBSERVED_FLOORS_RAW` is OPERATIVE; the post-REDSEA table is reference-only. The HANDOFF's
   "measure floors post-REDSEA" instruction (below) was wrong on this premise.
2. **`OBSERVED_FLOORS["CHGA"]` was a chimera** — plain-Otsu floor+spread while CHGA is registered
   as max(Otsu,Triangle). Corrected to `(811.9, 1.33, 4.3)`.
3. **`floor_method()`/`MARKER_INPUT_FLOOR_METHODS` have zero callers** — a doc registry, not an
   execution path. `codex_floors.py` (plan Phase 2) is the module that makes them run.

**Issue 3 (REDSEA zero-clamp) CLOSED, not a blocker:** clamp lands on background; Mesenchymal
2.95%→2.79%; prevalence-matched pre/post Other differs 0.15 pp. Moot anyway — RESTORE floors on RAW.

**Issue 2 (CD8a/HLADR) — mechanism found, NOT threshold-fixable:** `raw_otsu_floor` clips at p99.5;
markers whose positives live above the clip lose their mode and Otsu bisects background. Predictor
is dynamic range (p99.5/median). Recorded on the wrong (REDSEA) table; ~2.4× milder on operative
RAW. `max(Otsu,Triangle)` is a disguised quantile cut (CD8a+ ≈ CD3e+; HLADR misses 93–99% of
CD68+). Route to RESTORE + image review; do not "fix" with a floor statistic. CHGA is a third case.

**Compartment composition (new, not in old body):** at floors, Other ~61% / Epithelial ~28%, 67%
of app-Beta → Neural. This is a floors-as-positivity artifact, NOT the verdict — it becomes the
**acceptance test on RESTORE output** (plan Phase 5). If Other stays ~60% there, the gate set reopens.

**Upstream moved during the investigation** to `7227836` + 9 uncommitted files + a live screen run.
Load-bearing: **CD20 now excluded from RESTORE in BOTH roles** (`INVALID_REFERENCE_MARKERS` →
`RESTORE_EXCLUDED_MARKERS`, both-role import guard, new `DEFERRED_EXTRACTION_MARKERS`; B cells are
an abstention → `Immune-other`, and CODEX has no CD79a). CODEX's `IMMUNE_PAIRS` still lists
`["CD20","CD3e"]` — **must be dropped** (12 pairs, not 13). Not yet done here: it is bundled with
the upstream-pinning decision (commit upstream, or vendor a snapshot — CODEX must not build against
a working tree).

**Security:** the UF Navigator key IS in local history (`FEATURES_RESTORED_20251119.md`, 21
commits, stash `49b0f25`) — the earlier "never committed" was wrong. **Never pushed.** Rotate the
key FIRST, then purge (`git filter-repo` + stash/reflog expiry). Rotation is the user's action.

---

## Goal

CODEX Islet-Explorer replaces its scimap rules-gating phenotyping with the
`smith6jt-cop/Phenocycler_Analysis` method (REDSEA → RESTORE → ordered-residual lineage).
Labels **overwrite** the existing 21. The public PubApps launch is **gated** on it landing.

Upstream is now **VENDORED as a git submodule** at repo-root `external/Phenocycler_Analysis`
(remote `github.com/smith6jt-cop/Phenocycler_Analysis`, branch
**`redsea-restore-curated-5compartment`**, pinned at **`8db37af`** = METHOD v8). All 6 consumers
resolve it via `PHENOCYCLER_ROOT` env override → else the submodule; set `PHENOCYCLER_ROOT` to a
working-tree checkout only for upstream development. `main` is behind and lacks the CD20 work — do
not switch the submodule off this branch.

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
| QuPath per-compartment means | **DONE** 2,676,466 cells × 4 compartments, 11.5 GB, verified (15/15 donors, joins parquet totally, 10.3% no-nucleus) | `CellCompartments.tsv`; `verify_compartment_means.py` |
| Input floors (operative, RAW) | **DONE** measured + `--check` green; 13 screened markers registered | `codex_lineage.OBSERVED_FLOORS_RAW`, `codex_floors.py` |
| Compartment gates + pairs (12, no CD20) | **DONE**, validators + both-role guard enforce | `codex_lineage.py` (aligned to upstream 55b904d) |
| **RESTORE v6 screen (faithful Membrane+Cell NNMF)** | **DONE (first pass + immune)** — first-pass (0.995 clip): **88/90 PENDING**, 2 MODEL_UNSTABLE (VIM←ECAD 6505/6516). Immune (UNCLIPPED, see below): **79/90 PENDING** — CD8a/CD3e/CD68/CD163 15/15, CD4 14/15; HLADR←CD31 is the only weak pair (8 LOW_SIGNAL + 2 degenerate, genuine weak signal). | `codex_restore_validation.py`; `screen_first_pass`, `screen_immune` |
| CD8a/HLADR "anomaly" | **RESOLVED — was the p99.5 clip, NOT the markers.** Distributions are normal (bright tail above the clip on 2-3 low-rate donors). Fix = screen immune pairs unclipped (`IMMUNE_SCREEN_UPPER_QUANTILE=1.0`); validated on the RESTORE fit (6356 CD8a MODEL_UNSTABLE→VALID). `max(Otsu,Triangle)` rejected. | `codex_lineage.py` resolution note |
| **Review bundle (Phase 4)** | **DONE (first pass)** — 6 five-panel PDFs + 32-row `review_decisions.csv` + 32 QuPath overlay CSVs, visually verified (image panels populated: ECAD magenta / CD31 cyan). Immune bundle pending (needs a complete `evaluate_locked_cohort` results dir at the unclipped config). | `codex_restore_pair_review.py`; `review_first_pass/` |
| **First-pass review + gated parquet (Phase 5)** | **DONE** — reviewer ACCEPTED all 6 first-pass pairs → `--pair-reviews` re-run → **88/90 VALID** + per-donor divisors (ECAD 8582/CD31 2580/CD45 2623/PGP9.5 3411/VIM 3817/CHGA 650). `codex_positivity.py` wrote the gated compartment parquet (2.65M cells, `redsea>divisor`, ordered residual). Composition Epithelial 29.8% / Other 58.9% — **descriptive only**. | `reviews/pair_reviews_first_pass.csv`, `screen_first_pass_valid/`, `/data/codex_phenocycler/gated_restore/` |
| Immune bundle + lineage `cell_type` + app integration | not started | — |

**VALIDATION RULE (user directive 2026-07-23): do NOT use the app's scimap labels or scimap-derived
tests to judge the new phenotyping — they ARE the flawed method being replaced. Acceptance = the
image-based pair review; escalate uncertain calls to images/QuPath.** See memory
`feedback_no_app_labels_for_validation.md`.

Everything above is upstream of the retired entry point, so none of it is invalidated.
`object_id` joins raw↔corrected **totally** (2,676,466 both sides, all unique). The RESTORE blocker
is CLEARED: the adapter runs upstream `evaluate_locked_cohort` verbatim on the faithful Membrane+Cell
feature space (recorded `feature_compartments: ['Membrane','Cell']`, 4 features/pair — NOT Cell-only),
and emits reviewable rows + a divisor on ACCEPT. Remaining gate: the human pair review (Phase 4 bundle).

### Modules (`scripts/codex_phenocycler/`)

- `codex_panel.py` — marker/channel resolution. 11 tests run against the **real** OME-XML
  sidecars and QuPath TSV header, not fixtures.
- `codex_lineage.py` — compartment gates, RESTORE pairs, floors, immune sub-split. All
  validators run at import.
- `build_cells_parquet.py` — QuPath TSV → hive parquet (`--dry-run`, `--min-cell-area`).
- `export_cells_geojson.{groovy,sh}` — per-cell GeoJSON from the QuPath project.
- `verify_geojson.py` — the four silent-failure checks.
- `codex_redsea.py` — REDSEA driver (`--check` resolves all 15 without reading pixels).
- `verify_compartment_means.py` — fails closed on the compartment export (columns, 15/15 donors,
  total `object_id` join, NaN-nucleus fraction).
- `codex_floors.py` — measures the operative RAW floors + picks each marker's statistic by the
  dynamic-range × sparseness guard; `--check` regression-asserts `OBSERVED_FLOORS_RAW`.
- `codex_restore_validation.py` — the screen adapter. `patch_upstream()` injects CODEX's
  markers/floors/12-pairs + a pinned-SHA tripwire; `build_sample()`/`build_raw_view()` shape the
  inputs (donor-id normalized, all-NaN cells dropped, new-vs-old allclose preserved); `run()` calls
  upstream `evaluate_locked_cohort` on the default Membrane+Cell config. `--smoke` proves the load
  gate without scoring. `--immune` uses `run_robust()` (per-pair try/except) at the unclipped floor.
- `codex_restore_pair_review.py` — the review-bundle adapter. `patch_review()` overrides the Senior
  shortlist/rationale, rebinds `rpr.MARKER_INPUT_FLOOR_METHODS` (rpr has its own import), and patches
  `_qptiff_channel_names` (CODEX OME-XML attributes via `codex_panel`); `build_images_dir()` makes a
  no-extension symlink dir for image resolution; `build_bundle_inputs()` builds the full-compartment
  tree + matching clean raw/redsea for the strict `n == canonical_n` audit. 8 tests total.
- `codex_positivity.py` — applies accepted per-donor divisors (`redsea > divisor`) → ordered-residual
  compartment call → gated parquet + provenance. Composition is descriptive only (no app-label
  comparison). The immune sub-split + endocrine Beta/Alpha/Delta split are the not-yet-built lineage step.

Pin: `UPSTREAM_PINNED_SHA = 8db37af` (the submodule commit) + `EXPECTED_METHOD_VERSION =
restore-pair-validation-v8`. Now that upstream is a submodule the SHA pin is authoritative, but the
METHOD_VERSION tripwire is kept for the case where someone points `PHENOCYCLER_ROOT` at a working-tree
checkout. v7→v8 (2026-07-24) was audited ADDITIVE / contract-safe: adds `reference_selection.py`
(Gate-2 reference≥target frequency rule + `ACCEPTED_REFERENCE_OVERRIDES`) and a frozen `ACCEPTED_PAIRS`
tuple; no floor/divisor/clip/pair change, so the accepted first-pass divisors and 88/90 screen states
stand. Re-audit + bump both when a tripwire fires.

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
