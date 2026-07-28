# Handoff — Modernist redesign (2026-07-28)

Scope: the whole Shiny app. Branch `claude/modernist-redesign`, **12 commits**.
**Not deployed.** Production is untouched and still serving the pre-redesign app
from the original checkout.

Previous handoff archived to
`archive/docs/HANDOFF_2026-07-27_statistics-metadata-inflammation.md`.

---

## Where the work lives, and why it is a worktree not a branch

```
~/panc_CODEX/Islet-Explorer               main checkout, branch claude/spatial-cardc-…
                                          → this is what production serves
~/panc_CODEX/Islet-Explorer-modernist     git worktree, branch claude/modernist-redesign
                                          → all of this work
```

`/srv/shiny-server/islet-explorer` is a symlink **to the main checkout's
`app/shiny_app`**. Checking a branch out there would rewrite the files
shiny-server serves, and the next worker spawn would pick them up without anyone
touching `restart.txt`. A branch does not isolate production; a worktree does.

The worktree's `data/` entries are symlinks to the main checkout's (those files
are gitignored) and `.Renviron` was copied. `git status` is clean in both.

**Dev instance:** `http://10.15.152.7:7780/` — from the worktree's
`app/shiny_app`, `Rscript -e 'shiny::runApp(".", port=7780, host="0.0.0.0")'`.

---

## What prompted it

`app_redesign.R` appeared at the repo root, produced by the claude.ai **Design**
feature, with its supporting files in `Shiny app aesthetic redesign.zip`. Both
are preserved under `design/` with a README recording what was and was not
carried across. The design system is **Modernist** — Archivo, flat, zero radius,
2 px rules, flush left, one accent used sparingly.

It was a shell, not a replacement: synthetic data, no drill-down, no downloads,
no phenotype modal, no AI panel, no Inflammation tab, a 12-entry demo phenotype
palette against the app's 21, a theme at 10 pt against the app's documented 18,
and a cohort filter that z-scored every feature before any panel saw it.

---

## The port

Phases, in the order they were committed: worktree + vendored Archivo + CI
fences → vendored stylesheets (six defects fixed in the delivered pair) → dead
CSS and the dark-mode toggle that never was → Archivo reaches the plots +
palette aliases → **the flip** (Bootstrap 5, `page_fillable`, left rail,
`navset_hidden`) → inline-style strip + shared blocks → frozen control panel,
later split into a top bar and a right panel.

**Measured facts it rests on.** Shiny emits byte-identical markup under BS3 and
BS5, and bslib's BS5 bundle ships a BS3 compatibility layer, so the delivered
adapter's legacy selectors survive. But the BS5 **runtime** rewrites `li.active`
to `a.nav-link.active`, so the old layout jQuery returns `""` and the sidebar
would be permanently hidden — which is why the framework flip and that script's
deletion had to be a single commit.

`renderPlot` and `ggsave` already go through **ragg**, so
`systemfonts::register_font()` is sufficient; showtext is neither installed nor
wanted (its `dpi = 96` would have shrunk every download to a third scale).

`navset_hidden` mounts **every** panel and toggles CSS visibility, exactly like
`tabsetPanel` — so the drill-down duplicate-non-namespaced-ID guards are still
required. `IE_PIPELINE` / `IE_TAB_KEYS` in `00_globals.R` are the tab-identity
contract, and `app.R` asserts them at startup.

---

## Changes requested after the port

- **AI chat card removed.** Module files remain unreferenced under `R/`; delete
  `mod_ai_assistant_*.R` + `ai_helpers.R` to remove it entirely.
- **Donor-cohort table removed** from what is now Composition & Infiltration.
- **edgeR → compositional test.** The input was genuine integer cell counts, but
  every edgeR calibration assumes a transcriptome: `estimateDisp` borrows
  dispersion across features (empty across 21 types), TMM assumes most features
  are unchanged (false when β cells collapse), and "library size" here is tissue
  area imaged, not sequencing depth. Replaced with a centred log-ratio on **one
  value per donor**, computed on the FULL composition and only then subset.
  Scope is eight pre-specified endocrine/immune types, every pairwise stage
  contrast, and **no multiple-testing correction** — the types are pre-specified
  so there is no selection to correct for, and at 5 donors/stage BH would report
  "not significant" whatever the data said.
- **Nine phenotypes retired** — Endocrine, Acinar, Endothelial, T cell, Immune,
  Unknown, ECAD+, APCs, Structural. 21 → 12. Excluded via `PHENOTYPE_RETIRED` +
  `phenotype_keep()` / `phenotype_filter_rows()` applied at the two cell-level
  **loaders**, not deleted from the palettes: the data still carries the labels
  and every lookup is by exact name, so a palette without them would return `NA`
  and silently grey out real cells rather than dropping them.
- **Spatial split** into Spatial + Neighborhoods (same module namespace, so
  every input and output id is unchanged).
- **Pipeline reordered** to the order of work: Spatial → Neighborhoods →
  Composition & Infiltration → Trajectory → Plot → Statistics.
- **Significance is a star, never colour.** Red/green reads stop/go and
  bad/good; a p-value crossing a threshold is neither.
- Rail is logo-only, help consolidated into it, Viewer section renamed
  "The images themselves". Plots are transparent (the missing piece was
  `theme_islet()` inheriting `theme_minimal()`'s opaque white `plot.background`,
  not the device or the CSS). Frozen top bar holds what is identical on every
  tab; the right panel holds what changes.

---

## The Viewer — root cause found, not fully confirmed

**It has never worked, and the cause is not app code.** Avivator reads OME-TIFF
tile by tile with HTTP byte-range requests. Requesting bytes 100–199 of a 7.1 GB
image:

```
nginx        :8080/islet-explorer/   206 Partial Content, 100 bytes    works
shiny-server :3838                   200, Content-Length 7647738808    fails
runApp       :7780                   200, Content-Length 7647738808    fails
```

A server that ignores `Range` answers every tile fetch with the whole file. That
is the `AggregateError: Request failed` in the iframe console.

The image URL was derived from `session$clientData` — whichever origin the
browser used — so the Viewer worked or silently failed depending on how the app
was opened, with no message either way. **`IMAGE_BASE_URL` in `.Renviron` now
pins it** to a Range-capable origin. Verified the URL points at nginx and the
fetch errors stop.

**NOT confirmed: whether the image actually renders.** No canvas appeared within
25 s. That may be a 7 GB file fetching its first tiles, or a second problem
behind the first. **Check this first.** If it is still blank, the durable fix is
OME-Zarr — chunked, many small files, no Range dependency, which is what
Viv/Avivator expect.

---

## Verification — and the gap that cost real time

`parse()` and the 71 tests **both pass on an app that dies on every connect.**
So does `curl → 200`, which only proves the static page was served. A broken
`spatial_server` signature reached the user exactly that way.

Always green:

```bash
cd app/shiny_app
Rscript -e 'for (f in c("app.R", list.files("R","[.]R$",full.names=TRUE))) parse(f)'
Rscript tests/run_tests.R                    # 71 pass
../../scripts/check_ui_style_fences.sh       # fence 2 clean; fence 1 ratchets
```

Plus two browser checks, written this session but **not committed — worth
committing**:

1. **Smoke test** — open a session, click every rail entry, fail on a
   `#shiny-disconnected-overlay` or any non-validation `.shiny-output-error`.
2. **Mount audit** — for every module, every `output$x <-` must have a matching
   `ns("x")` in its UI file. This caught the infiltration card, whose plot,
   table, note and download button had been deleted by an over-broad regex,
   leaving a section heading with nothing under it. All modules currently clean.

**Two self-inflicted bugs, same cause: regex edits spanning block boundaries.**
One turned every `df[rows, , drop = FALSE]` into `df[rows, drop = FALSE]` —
which selects columns, not rows — and it parsed cleanly with all tests passing.
The other deleted the infiltration card. Both were caught and reverted, but:
**use line-anchored edits with an explicit before/after diff for structural
removals**, never pattern matching across blocks.

---

## Open, highest value first

1. **Confirm the Viewer renders.** See above.
2. **The cohort filter reaches Plot and Spatial only.** The top bar renders on
   every tab but `trajectory_server`, `statistics_server` and
   `inflammation_server` still read their own filter state — so it is *visible*
   there and does not yet *drive* them. `cohort$apply()` / `cohort$donor_ids()`
   / `cohort$provenance()` are ready; this is wiring. Statistics additionally
   needs the inherit/override control (user-approved: the bar sets the default,
   the Advanced block may override, and an override is announced in the banner
   and the CSV header).
3. **38 inline skin declarations remain**, all in multi-line `paste0()`-built
   style strings a single-line rewrite cannot reach — 12 in
   `mod_statistics_server.R`, 10 in `mod_trajectory_server.R`, 5 in
   `mod_plot_server.R`.
4. **Raw column encodings still leak into the UI** (`prop_Beta cell`,
   `Capillary|islet_core`). Needs one `feature_display_label()` covering every
   prefix plus the `className|region` form, kept strictly separate from the
   lookup key — every lookup here is by exact column name, so a display label
   must never become one.
5. **The `statistics.md` open list** — `eta_sq` is really epsilon-squared;
   sex-stratified tests use pooled-SD while Section 2 uses Welch, with reversed
   contrast labels; an over-restrictive filter returns `NULL` and is
   indistinguishable from "you never clicked Run"; the covariate card says "does
   NOT survive adjustment" when nothing was significant; `outlier_summary` is
   counted before the min-cells/diameter filters.
6. **Shapiro-Wilk still colours** green "normal" / orange "non-normal". Not
   significance, so it was left — but it is the same good/bad framing applied to
   something that is neither.
7. **`INFLAM_FDR_FAMILY`** is 4 names after the retirement. Revisit when the
   phenotypes are re-gated; it is fixed in source deliberately.

**Blocked on the user, not on us:** the phenotyping is being reworked. Nothing
in the composition card should be read as a finding until it is.

---

## Deploying, when it is time

```bash
cd ~/panc_CODEX/Islet-Explorer
git merge claude/modernist-redesign
touch app/shiny_app/restart.txt
curl -sI http://127.0.0.1:3838/islet-explorer/  | head -1
curl -sI http://10.15.152.7:8080/islet-explorer/ | head -1
```

**Set `IMAGE_BASE_URL` in the deployed `.Renviron`** to the nginx origin, or the
Viewer breaks for anyone reaching the app another way. `.Renviron` is gitignored
— it will not travel with the merge.

Docs to update at merge: `CLAUDE.md` (shell architecture, tab keys),
`app/shiny_app/CLAUDE.md` (palette syncing → rail-owned; add the `ie_segmented`
`updateRadioButtons(choices=)` rule), `docs/statistics.md` (significance is a
star now).
