# Shiny App: Islet Area Distributions

Interactive Shiny app to explore islet area distributions across donor status (ND → Aab+ → T1D) with:

- Band-region Islet_Markers: % positive per islet vs. size (mean ± SE)
- Band-region Islet_Targets: density/count vs. size (mean ± SE)
- Islet composition vs. size (INS/GCG/SST fractions)
- One-way ANOVA per size bin across donor groups (BH-adjusted p-values)

Data source: `../data/master_results.xlsx` (generated from TSVs in `../data/results` with donor key `CODEX_Pancreas_Donors.xlsx`).

## Run

Prereqs: R >= 4.2 with packages: `shiny`, `readxl`, `dplyr`, `tidyr`, `stringr`, `ggplot2`, `plotly`, `broom`.

```r
shiny::runApp("shiny")
```

Or from the shell:

```bash
R -q -e 'shiny::runApp("shiny", launch.browser=TRUE)'
```

If you are missing packages, install them via:

```r
source("scripts/install_shiny_deps.R")
```

## Segmentation Viewer (Vitessce)

The Trajectory tab embeds a Vitessce-based OME-TIFF viewer which overlays one or more segmentation label rasters on the base image.

- Base images live under `app/shiny_app/www/local_images/`, named like: `<anything>_<YYYY>.ome.tiff` (YYYY is a 4-digit sample ID)
- Segmentation label images live under `app/shiny_app/www/LabelExports/`, named like: `YYYY_<label>.ome.tiff` (e.g., `6505_Islet_42.ome.tiff`)
- Multiple segmentations per sample are supported; they're auto-detected and shown as toggleable layers.

Tips:
- Ensure both base and segmentation files are OME-TIFFs. Non-OME TIFFs won’t load.
- Click a point in the trajectory scatter to load the corresponding case image; if an islet ID is present in the filename, matching segmentations are preferred.
- Channel names can be optionally provided via `app/shiny_app/Channel_names` or `Channel_names.txt` to color INS/GCG/SST/DAPI by default.

## Notes

- Islet diameter is standardized to CORE area: derived from `Islet_Targets` rows where `type == "islet_core"` (core `region_um2`). Diameter = `2 * sqrt(core_area/pi)`.
- Band-region rows are used for markers (`region_type == "islet_band"`) and targets (`type == "islet_band"`).
- Binning is by diameter (µm); default bin width = 50 µm. Adjust in the sidebar.
- Markers use `% positive` (100 × `pos_frac`). Targets use `area_density` by class.
- Composition uses `*_any / cells_total` for INS/GCG/SST.
- Stats tab runs one-way ANOVA per bin across groups (ND, Aab+, T1D) and reports BH-adjusted p-values.

