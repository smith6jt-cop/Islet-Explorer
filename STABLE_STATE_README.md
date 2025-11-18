# Islet-Explorer Stable State

**Date**: November 18, 2025
**Branch**: `stable`
**Status**: ✅ STABLE - Core functionality working

## What's Working

### Core Shiny App
- **Plot Tab**: Marker/Target/Composition analysis with interactive plotly charts
- **Statistics Tab**: One-way ANOVA per size bin across donor groups
- **Trajectory Tab**: Pseudotime analysis with correlation statistics
- **Viewer Tab**: Whole-slide AVIVATOR viewer for OME-TIFF images

### Data Infrastructure
- Master Excel data: `data/master_results.xlsx` (11MB)
- Annotations: `data/annotations.tsv` (72MB spatial data)
- Islet spatial lookup: `data/islet_spatial_lookup.csv` (661KB)
- OME-TIFF images: 15 whole-slide images in `app/shiny_app/www/local_images/` (109GB)
- AnnData files: 4 trajectory analysis files in `scripts/` (~5MB each)

### Features
- Filter by donor status (ND, Aab+, T1D)
- Autoantibody filtering for Aab+ group
- Size-based binning and analysis
- Normalization options (None, Global z-score, Robust per-donor)
- Export-ready statistics
- AVIVATOR whole-slide image viewer with channel configuration

## What Was Removed

The following features were attempted but failed and have been removed:
- **VitessceR integration** - Package has schema validation bugs
- **Cytomapper viewer** - Integration issues
- **Segmentation overlays in trajectory tab** - Broken due to vitessceR

## Running the App

### Local Development
```bash
cd /home/smith6jt/panc_CODEX/Islet-Explorer/app/shiny_app
R -e "shiny::runApp('.', launch.browser=TRUE)"
```

### With Shiny Server
The app is configured for deployment with shiny-server + nginx reverse proxy.

## Git Repository State

### Current Branch Structure
- **stable** (current) - Clean, working app
- **main** - Base state before experimental features
- **master** - Old state with GeoJSON experiments
- **backup/** branches - From failed rollback attempts (can be deleted)

### Recent Commits on Stable Branch
```
8e828fe - Viewer: use app-absolute image_url for reverse proxy
7a712c3 - Update .gitignore to exclude workspace files
d02e822 - Enhance interactive islet lookup feature (base)
```

## Next Steps

### Phase 1: Implement Simple JavaScript Overlays (Recommended)
Goal: Add segmentation boundary overlays to AVIVATOR viewer using canvas/SVG

**Approach**:
1. Use existing GeoJSON/annotation data
2. Draw polygons on canvas overlay above AVIVATOR iframe
3. Add layer controls to toggle segmentation types
4. Much simpler than vitessceR - ~200 lines of JavaScript

**Timeline**: 2-3 days

**Benefits**:
- No R package dependencies
- Direct control over rendering
- Fast and lightweight
- Works with existing AVIVATOR viewer

### Phase 2: Enhanced Features (Future)
- Click-to-zoom on islets in trajectory plot
- Export annotated images
- Real-time measurement tools
- Integration with anndata for single-cell overlays

## File Organization

```
Islet-Explorer/
├── app/shiny_app/
│   ├── app.R                    # Main Shiny app (2,336 lines)
│   ├── Channel_names            # Channel configuration
│   └── www/
│       ├── avivator/            # AVIVATOR viewer
│       └── local_images/        # 15 OME-TIFF files (109GB)
├── data/
│   ├── master_results.xlsx      # Main analysis data
│   ├── annotations.tsv          # Spatial annotations
│   └── islet_spatial_lookup.csv # Islet centroids
├── scripts/
│   ├── *.h5ad                   # AnnData trajectory files
│   ├── *.R                      # R analysis scripts
│   └── *.py                     # Python preprocessing
└── analysis_output/             # Generated plots and results
```

## Dependencies

### R Packages
- shiny, shinyjs
- readxl, dplyr, stringr, tidyr
- ggplot2, plotly, broom
- jsonlite, base64enc

### Optional
- anndata (for trajectory analysis in R)

## Troubleshooting

### App won't start
1. Check that `data/master_results.xlsx` exists
2. Verify R packages are installed: `scripts/install_shiny_deps.R`
3. Check working directory is `app/shiny_app/`

### Images don't load in viewer
1. Verify images exist in `app/shiny_app/www/local_images/`
2. Check Channel_names file matches image channel count
3. Confirm AVIVATOR index.html exists

### Data errors
1. Run NA audit: `source('../../scripts/na_audit.R')`
2. Check Excel file integrity
3. Verify islet_key generation logic

## Contact & Support

For issues with this stable branch, check:
- Git commit history for recent changes
- AI_Markdown/ directory for implementation notes (gitignored)
- analysis_output/ for generated results (gitignored)

---

**This stable state prioritizes reliability over features. All experimental code has been removed.**
