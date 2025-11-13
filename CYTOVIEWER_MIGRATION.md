# Migration from vitessce/AVIVATOR to cytoviewer

## Summary

This document describes the migration from the ineffective vitessce/AVIVATOR implementation to cytoviewer for multichannel islet visualization in the Shiny app.

## Changes Made

### 1. Removed vitessce/AVIVATOR Code

**Package Loading (app.R lines 25-66)**
- Removed vitessceR package loading and verification
- Removed AVIVATOR setup code
- Added cytoviewer package check

**Functions Removed:**
- `discover_islet_assets()` - Asset discovery for vitessce
- `choose_islet_asset()` - Asset selection for vitessce
- `build_islet_vitessce()` - Vitessce viewer configuration
- `resolve_avivator_base()` - AVIVATOR base path resolution
- `build_channel_config_b64()` - AVIVATOR channel config builder

**UI Components Removed:**
- "Viewer" tab with AVIVATOR embed (lines 1829-1834)
- vitessce fallback UI (lines 3346-3361)
- vitessce_assets reactive (lines 3364-3381)
- traj_vitessce_viewer renderUI (lines 3389-3597)
- viewer_info reactive for AVIVATOR (lines 3599-3674)
- viewer-mode observe block (lines 3676-3687)
- local_image_picker renderUI for AVIVATOR (lines 4394-4491)

**Total Lines Removed:** ~859 lines of code

### 2. Removed Files and Directories

```
app/shiny_app/www/vitessce_data/          # All JSON annotations
app/shiny_app/www/vitessce_libs/          # Vitessce JavaScript libraries
app/shiny_app/www/vitessce_*.html         # Test HTML files
app/shiny_app/www/vitessce_*.js           # Loader scripts
test_vitessce*.R                          # Test R scripts
TRAJECTORY_VITESSCE_TEST_GUIDE.md         # Documentation
scripts/Shell_scripts/*vitessce*.sh       # Shell scripts
scripts/Shell_scripts/install_avivator.sh # AVIVATOR installer
scripts/Python/*vitessce*.py              # Python conversion scripts
tmp/*vitessce*                            # Temporary test files
```

### 3. Added cytoviewer Integration

**Package Loading**
```r
# Load cytoviewer for multichannel image visualization
suppressPackageStartupMessages({
  if (!requireNamespace("cytoviewer", quietly = TRUE)) {
    message("Package 'cytoviewer' not found. Install with: BiocManager::install('cytoviewer')")
  }
})
```

**UI Update**
Changed from:
```r
uiOutput("traj_vitessce_viewer")
```

To:
```r
uiOutput("islet_viewer")
```

**Server Logic**
Added new `output$islet_viewer` renderUI that:
1. Shows placeholder when no islet is selected
2. Checks for required packages (cytomapper, EBImage)
3. Finds image files for the selected case
4. Provides setup instructions for full cytomapper integration
5. Displays selected islet information (case ID, coordinates, annotations)

### 4. Updated Configuration

**APP_VERSION**
- Changed from: `"vitessce-trajectory-v1"`
- Changed to: `"cytoviewer-v1"`

**.gitignore**
Added patterns to exclude removed vitessce/AVIVATOR files

## Implementation Status

### ✅ Completed
- All vitessce/AVIVATOR code removed
- All vitessce/AVIVATOR files deleted
- Basic cytoviewer UI placeholder implemented
- App loads without syntax errors
- Selected islet information displays correctly

### 🔄 Next Steps for Full Integration

To complete the multichannel viewer implementation, you need to:

1. **Install Required Packages**
   ```r
   BiocManager::install('cytoviewer')
   BiocManager::install('cytomapper')
   BiocManager::install('EBImage')
   ```

2. **Prepare Data in Required Format**
   - Convert OME-TIFF files to `CytoImageList` objects
   - Load segmentation masks as `CytoImageList`
   - Create `SingleCellExperiment` with cell metadata

3. **Implement Viewer Rendering**
   Use cytomapper's plotting functions:
   ```r
   # Example using plotPixels
   output$islet_viewer <- renderPlot({
     cytomapper::plotPixels(
       image = images_list,
       mask = masks_list,
       object = sce_object,
       img_id = selected_case_id,
       colour_by = c("DAPI", "INS", "GCG"),
       bcg = list(DAPI = c(0, 1, 1.5),
                  INS = c(0, 1, 2),
                  GCG = c(0, 1, 2))
     )
   })
   ```

4. **Resources**
   - [cytomapper documentation](https://bodenmillergroup.github.io/cytomapper/articles/cytomapper.html)
   - [cytoviewer GitHub](https://github.com/BodenmillerGroup/cytoviewer)

## Benefits of cytoviewer

1. **Built for Bioconductor** - Better integration with R/Bioconductor ecosystem
2. **Designed for Multiplexed Imaging** - Specifically built for imaging mass cytometry
3. **Shiny Integration** - Has proper Shiny functions (plotPixels, plotCells)
4. **Maintained** - Active development by Bodenmiller Group
5. **No JavaScript dependencies** - Pure R implementation

## Files Modified

- `app/shiny_app/app.R` - Main application file
- `.gitignore` - Added patterns for removed files

## Backup

A backup of the original app.R was created at:
`app/shiny_app/app.R.backup`

## Testing

The app was tested and successfully loads with:
```
[APP LOAD] Version=cytoviewer-v1
[CHANNELS] Loaded 35 channel names
[SETUP] Found 15 images in www/local_images
[SEGMENTATION] Loaded 576589 segmentation records
[SPATIAL] Loaded 10736 islet spatial records
```

No syntax errors were encountered during loading.
