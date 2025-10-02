# Getting Started

This guide will walk you through your first session with Islet-Explorer.

## Prerequisites

- Completed [Installation](installation.md)
- Preprocessed GeoJSON overlays (if using Viewer tab)
- Data file: `data/master_results.xlsx` (included with app)

## Launching the Application

### Local Development

```bash
cd /path/to/Islet-Explorer
R -q -e "shiny::runApp('app/shiny_app', launch.browser=FALSE)"
```

The app will start on a random port. Look for output like:

```
Listening on http://127.0.0.1:4329
```

Open this URL in your web browser.

### Using VS Code Task

If using VS Code, run the pre-configured task:

1. Press `Ctrl+Shift+P` (or `Cmd+Shift+P` on Mac)
2. Type "Run Task"
3. Select "Run Shiny app to reproduce issues"

## Understanding the Interface

### Main Layout

The app consists of four main tabs:

```
┌─────────────────────────────────────────────────────┐
│  [Plot] [Statistics] [Trajectory] [Viewer]          │
├─────────────────────────────────────────────────────┤
│  Sidebar          │  Main Panel                     │
│  ────────────     │  ─────────────                  │
│  • Filters        │  • Plots                        │
│  • Options        │  • Tables                       │
│  • Parameters     │  • Visualizations               │
└─────────────────────────────────────────────────────┘
```

### Donor Groups

Data is organized by three donor status groups:

- **ND** (Non-Diabetic) - Control group, shown in blue
- **Aab+** (Autoantibody Positive) - At-risk group, shown in orange
- **T1D** (Type 1 Diabetes) - Disease group, shown in red

These colors are consistent throughout the app.

## Basic Workflow

### 1. Explore Composition Data

**Navigate to**: Plot tab (default)

1. **Select mode**: "Composition" (default)
2. **Choose metric**: `Ins_any`, `Glu_any`, or `Stt_any`
3. **View plot**: Scatter plot shows values by islet size and donor status
4. **Examine distribution**: Click "Distribution plots" to see violin/box plots

**What you're seeing**: Percentage of cells expressing insulin (Ins), glucagon (Glu), or somatostatin (Stt) across different islet sizes.

### 2. Analyze Markers

1. **Switch mode**: Select "Markers"
2. **Choose marker**: e.g., "KI67" (proliferation marker)
3. **Select metric**: "Fraction positive" or "Counts"
4. **Choose region**: "islet_core", "islet_band", or "islet_union"

**What you're seeing**: Cell marker expression patterns across donor groups and islet sizes.

### 3. Run Statistical Tests

**Navigate to**: Statistics tab

1. **Select data**: Choose mode (Targets/Markers/Composition)
2. **Pick metric**: Choose what to test
3. **View results**:
   - Global model (donor status + islet size effects)
   - Pairwise comparisons (ND vs Aab+, etc.)
   - Adjusted p-values (BH correction)

**Interpretation**: Significant differences between donor groups after accounting for islet size.

### 4. Explore Pseudotime Trajectories

**Navigate to**: Trajectory tab

1. **View UMAP**: Two panels show cell clustering
   - Left: Colored by donor status
   - Right: Colored by selected feature
2. **Select feature**: Choose from markers or targets
3. **Examine trends**: See how features change along pseudotime axis

**What you're seeing**: Computational lineage analysis showing beta cell differentiation/dedifferentiation trajectories.

### 5. View Images with Overlays

**Navigate to**: Viewer tab

1. **Select image**: Choose from dropdown (e.g., "ND_0112.ome.tiff")
2. **Wait for load**: Large images may take 10-20 seconds
3. **Verify overlays**: Green cell boundaries should appear automatically
4. **Toggle overlays**: Use checkbox to show/hide
5. **Navigate image**:
   - **Pan**: Click and drag
   - **Zoom**: Mouse wheel or pinch gesture
   - **Channels**: Use viewer controls (top-left)

**What you're seeing**: Multiplexed fluorescence image with ~187,000 segmented cells overlaid.

## Common Tasks

### Filter by AAb Status (Aab+ Group Only)

1. In Plot or Statistics tab, expand "AAb filters"
2. Select "Any" (has any autoantibody) or "All" (has all autoantibodies)
3. Choose specific antibodies: GADA, IA2A, ZnT8A, IAA, mIAA

**Note**: AAb filters only affect the Aab+ group; ND and T1D remain unchanged.

### Compare Specific Islet Sizes

1. Adjust "Islet diameter range" slider
2. Plots update to show only islets within selected range
3. Use "Bin width" to change smoothing/resolution

### Normalize Data

1. Select "Curve normalization" option:
   - **None**: Raw values
   - **Global z-score**: Standardize across all donors
   - **Robust per-donor**: Median/MAD within each donor

2. Distribution plots respect the same normalization

### Export Data

Currently, export screenshots or plots:
1. Right-click on plot → "Save image as..."
2. Or use browser print function to save as PDF

## Example Analysis

### "Are beta cells reduced in T1D islets?"

1. **Plot tab** → Composition → `Ins_any`
2. **Observe**: Red points (T1D) lower than blue (ND)
3. **Statistics tab** → Select same metric
4. **Result**: Significant reduction (p < 0.001)
5. **Trajectory tab** → Feature: INS expression
6. **Interpretation**: T1D beta cells show lower INS even after accounting for size

## Tips and Best Practices

### Performance

- **Large images**: Initial load is slow; subsequent pans/zooms are fast
- **Many overlays**: Reduce feature limit if browser slows down
- **Big plots**: Use "Distribution plots" for cleaner comparisons

### Data Interpretation

- **Always check islet size**: Many effects correlate with islet diameter
- **Use statistics**: Visual differences may not be significant
- **Consider AAb subtypes**: Aab+ is heterogeneous
- **Look at regions**: Core vs band can show different patterns

### Troubleshooting

- **Plot doesn't update**: Try switching tabs and back
- **Overlays missing**: Check browser console (F12) for errors
- **App crashes**: May need to restart R session
- **Data NA warnings**: Check console output, usually benign

## Next Steps

- [Plot Tab Detailed Guide](plot_tab.md) - Advanced plotting features
- [Statistics Guide](statistics_tab.md) - Statistical methods explained
- [Viewer Tab Guide](viewer_tab.md) - Image viewer and overlays in depth
- [Data Format](../technical/data_format.md) - Understanding input data structure

## Getting Help

- Check [FAQ](../reference/faq.md) for common questions
- See [Troubleshooting](../reference/troubleshooting.md) for known issues
- Open an issue on GitHub for bugs or feature requests
