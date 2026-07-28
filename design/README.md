# Delivered design package (Modernist)

Source material for the Jul 2026 UI redesign. **Nothing here is loaded by the app** — it is the
as-delivered artefact, kept so the port can be diffed against its origin.

| Path | What |
|---|---|
| `app_redesign.R` | The shell as it landed at the repo root on 2026-07-28. Byte-identical to `delivered/shiny/app.R`. |
| `Shiny app aesthetic redesign.zip` | The package as delivered. |
| `delivered/shiny/` | The same zip, unpacked, so the files are readable and greppable in-repo. |

Produced by the claude.ai **Design** feature. The design system it derives from is the
**"Modernist"** project on claude.ai/design — `delivered/shiny/www/styles.css` is that project's
stylesheet, unmodified (it still carries the system's original `#ec3013` accent).

## What the app actually adopts

The delivered shell runs on synthetic `demo_islets()` data and is not a replacement — it omits the
drill-down segmentation panel, all 21 download buttons, the phenotype-rules modal, the AI assistant,
and the whole Inflammation tab. It is the *design*, not the app.

Three things in it are deliberately **not** carried across, each for a measured reason:

- **`PAL_PHENO`** is a 12-entry demo palette (`Beta`, `CD8 T`, `Vessel`, …). The app has 21 real
  phenotypes across four selectable colourblind-safe palettes. Only 5 of the 12 keys match, and with
  `scale_*_manual(drop = FALSE)` the rest render `NA`. The port aliases it to `PHENOTYPE_COLORS`.
- **`theme_modernist(base_size = 12)`** puts axis text at 10 pt and blanks legend titles. The app's
  scale is base 18 / axis 16, raised deliberately in Jul 2026 for legibility, and four legends are
  named (`"Structures"`, `"Phenotype"`, `"Region"`, `"Cluster"`). The port keeps `theme_islet()` as
  the single source of truth and takes only the geometry.
- **`normalize_features()`** z-scores every feature before any panel sees the frame. The Statistics
  tab reports unstandardised differences in the feature's own units, and Inflammation's edgeR path
  needs integer counts. The port drops it; the cohort bar owns row selection only.

`modernist-shiny.css` overrides the design system's accent to crimson `#b71533` + gold `#c8a02e`.
That override is authoritative — it is further from the Aab+ donor colours than the original
orange-red was.

See `docs/HANDOFF.md` for the port itself.
