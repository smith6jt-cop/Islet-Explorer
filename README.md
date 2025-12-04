# Islet Explorer

Interactive web applications for exploring human pancreatic islet measurements from multiplex imaging studies. Analyze islet composition, cell type distributions, and marker expression patterns across donor groups (Non-Diabetic, Autoantibody-positive, Type 1 Diabetic).

<div align="center">

[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue.svg)](https://github.com/smith6jt-cop/Islet-Explorer)
[![R Shiny](https://img.shields.io/badge/R-Shiny-blue.svg)](https://shiny.rstudio.com/)
[![Python Panel](https://img.shields.io/badge/Python-Panel-green.svg)](https://panel.holoviz.org/)

</div>

## Overview

Islet Explorer provides two complementary interactive web applications for analyzing pancreatic islet data:

- **R Shiny App** (`app/shiny_app/`) - Production-ready application with trajectory analysis and embedded OME-TIFF viewer
- **Python Panel App** (`app/panel_app/`) - Alternative implementation with AI-powered analysis assistant

Both applications work with the same data format (`master_results.xlsx`) and provide:
- Interactive scatter plots of islet metrics by size and donor group
- Statistical analysis (ANOVA, pairwise comparisons, effect sizes)
- Islet composition analysis (insulin, glucagon, somatostatin fractions)
- Quality assurance and data diagnostics
- OME-TIFF image viewing capabilities

## Key Features

### Data Analysis
- **Donor Group Comparison**: Analyze differences between ND (Non-Diabetic), Aab+ (Autoantibody-positive), and T1D (Type 1 Diabetic) donors
- **Size-Based Analysis**: Explore how islet size correlates with cellular composition and marker expression
- **Autoantibody Filtering**: Filter Aab+ donors by specific autoantibodies (GADA, IA2A, ZnT8A, IAA, mIAA)
- **Normalization Options**: Apply global z-score or robust per-donor normalization
- **Statistical Testing**: One-way ANOVA with Benjamini-Hochberg correction, pairwise t-tests

### Visualization
- **Interactive Plots**: Mean ± SE plots with plotly/ggplot2 interactivity
- **Distribution Comparisons**: Violin plots and density distributions by donor group
- **Trajectory Analysis** (R Shiny): UMAP and pseudotime analysis with slingshot
- **Image Viewer**: Embedded Avivator/Viv viewer for local OME-TIFF files

### Data Processing
- **Region Synthesis**: Automatically compute union regions from core + band measurements
- **Quality Checks**: Built-in validation for data completeness and consistency
- **Flexible Metrics**: Support for markers, targets, and composition data

## Quick Start

### R Shiny Application (Recommended)

```bash
# Install R dependencies
Rscript scripts/install_shiny_deps.R

# Run the application
R -q -e "shiny::runApp('app/shiny_app', launch.browser=TRUE)"
```

The app will open in your default web browser at `http://127.0.0.1:XXXX`.

### Python Panel Application

```bash
# Create virtual environment
cd app/panel_app
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run the application
panel serve app.py --show --autoreload
```

## Installation

### Prerequisites

**For R Shiny App:**
- R >= 4.2.0
- CRAN packages: `shiny`, `readxl`, `dplyr`, `tidyr`, `stringr`, `ggplot2`, `plotly`, `broom`
- Optional: Bioconductor packages for trajectory analysis (`slingshot`, `SingleCellExperiment`)

**For Python Panel App:**
- Python >= 3.9
- See `app/panel_app/requirements.txt` for full dependency list

**For Data Processing:**
- Python >= 3.9 with `pandas`, `openpyxl` for building `master_results.xlsx`

### Detailed Installation

#### R Shiny App

```bash
# Automated installation of all dependencies
Rscript scripts/install_shiny_deps.R

# Optional: Install Bioconductor packages for trajectory analysis
INSTALL_VITESSCER=1 Rscript scripts/install_shiny_deps.R
```

#### Python Panel App

```bash
cd app/panel_app
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

#### Avivator Image Viewer (Optional)

```bash
# Install static Avivator bundle for OME-TIFF viewing
./scripts/install_avivator.sh
```

This installs the Avivator viewer into `app/shiny_app/www/avivator/`.

## Data Requirements

### Input Data Format

The applications require a `master_results.xlsx` file with the following sheets:

1. **Islet_Targets**: Target density and count measurements
   - Columns: `Case ID`, `Donor Status`, `region`, `class`, `region_um2`, `area_um2`, `area_density`, `count`
   - Region format: `Islet_<number>_<type>` (e.g., `Islet_64_core`, `Islet_64_band`)

2. **Islet_Markers**: Marker positivity per islet region
   - Columns: `Case ID`, `Donor Status`, `region`, `marker`, `n_cells`, `pos_count`, `pos_frac`

3. **Islet_Composition**: Cell type composition
   - Columns: `Case ID`, `Donor Status`, `islet_id`, `cells_total`, `Ins_any`, `Glu_any`, `Stt_any`

4. **LGALS3** (Optional): Additional LGALS3 marker data

### Donor Metadata

Donor information with:
- `Case ID`: Zero-padded 4-digit string (e.g., "0001", "0064")
- `Donor Status`: One of ND, Aab+, T1D
- Autoantibody flags: `AAb_GADA`, `AAb_IA2A`, `AAb_ZnT8A`, `AAb_IAA`, `AAb_mIAA`

### Building master_results.xlsx

```bash
# From TSV files in data/results/
python scripts/build_master_excel.py

# From a ZIP archive of donor data
python scripts/compile_donors.py --input donors.zip --output data/master_results.xlsx
```

## Project Structure

```
Islet-Explorer/
├── app/
│   ├── shiny_app/           # R Shiny application
│   │   ├── app.R            # Main application file
│   │   ├── README.md        # Shiny-specific documentation
│   │   └── www/             # Static assets and image viewer
│   │       └── avivator/    # Embedded OME-TIFF viewer
│   └── panel_app/           # Python Panel application
│       ├── app.py           # Main application file
│       ├── requirements.txt # Python dependencies
│       ├── components/      # UI components
│       └── utils/           # Data processing utilities
├── scripts/
│   ├── install_shiny_deps.R      # Install R dependencies
│   ├── install_avivator.sh       # Install image viewer
│   ├── build_master_excel.py     # Build Excel from TSVs
│   ├── compile_donors.py         # Compile donor data from ZIP
│   ├── compute_pseudotime.R      # Trajectory analysis
│   ├── convert_to_ometiff.sh     # Image format conversion
│   └── ...                       # Additional utility scripts
└── data/                    # Data directory (not in repo, see .gitignore)
    └── master_results.xlsx  # Main data file
```

## Usage

### Application Tabs

#### Plot Tab
Visualize islet metrics as a function of size:
1. Select metric type (Markers, Targets, or Composition)
2. Choose specific marker/target/component
3. Apply normalization (optional)
4. Filter by autoantibodies (for Aab+ donors)
5. View mean ± SE plots and distribution comparisons

#### Statistics Tab
Perform statistical analysis:
- Global model: `value ~ donor_status + islet_diam_um`
- Pairwise comparisons with size correction
- Effect size calculations
- Benjamini-Hochberg adjusted p-values

#### Trajectory Tab (R Shiny)
Analyze developmental trajectories:
- UMAP visualization of islet space
- Pseudotime computation with slingshot
- Gene expression patterns along trajectories

#### Viewer Tab
View OME-TIFF multiplex images:
- Embedded Avivator viewer
- Channel-specific visualization
- Configurable color mappings

#### QA Tab
Data quality diagnostics:
- Region count validation
- Missing value analysis
- Donor group summaries

### Keyboard Shortcuts & Tips

- **R Shiny**: Use the sidebar to control filtering and normalization
- **Python Panel**: AI assistant available for data interpretation (requires API key)
- **Plots**: Hover over points for detailed information
- **Export**: Download plots and tables using built-in export features

## Development

### Running in Development Mode

**R Shiny:**
```bash
# With automatic reloading (in RStudio)
shiny::runApp('app/shiny_app')

# From command line
R -q -e "shiny::runApp('app/shiny_app', launch.browser=FALSE)"
```

**Python Panel:**
```bash
cd app/panel_app
panel serve app.py --show --autoreload --dev
```

### Testing Data Preparation

```bash
# Test data loading and preparation
Rscript scripts/test_shiny_prep.R

# Diagnose islet data quality
Rscript scripts/diagnose_islets.R
Rscript scripts/diagnose_marker_regions.R

# Audit for missing values
Rscript scripts/na_audit.R
```

### Working with Images

```bash
# Convert images to OME-TIFF format
./scripts/convert_to_ometiff.sh input.tiff output.ome.tiff

# Convert to OME-ZARR format (alternative)
./scripts/convert_to_omezarr.sh input.tiff output.zarr

# Annotate channel names
python scripts/annotate_ome_tiff_channels.py input.ome.tiff channel_names.txt

# Auto-color channels
python scripts/auto_color_ome_tiff_channels.py input.ome.tiff
```

### Environment Variables

**R Shiny App:**
- `LOCAL_IMAGE_ROOT`: Directory containing OME-TIFF files for viewer

**Python Panel App:**
- `ISLET_DATA_PATH`: Path to `master_results.xlsx` (default: `../../data/master_results.xlsx`)
- `LOCAL_IMAGE_ROOT`: Directory containing OME-TIFF files
- `NAVIGATOR_API_URL`: AI API endpoint for analysis assistant
- `NAVIGATOR_API_KEY`: API key for AI features

## Data Conventions

### Islet Keying
All data joins use three keys:
- `Case ID`: Donor identifier (zero-padded 4-digit string)
- `Donor Status`: ND, Aab+, or T1D
- `islet_key`: Derived from region name (e.g., `Islet_200` from `Islet_200_core`)

### Region Types
Three normalized region types:
- `islet_core`: Central region of islet
- `islet_band`: Peripheral band region
- `islet_union`: Combined core + band region

### Islet Diameter
Calculated from CORE region area:
```
islet_diam_um = 2 * sqrt(core_region_um2 / pi)
```

### Region Synthesis
If union rows are missing:
- **Targets**: Synthesize `count` only as `union = core + band`
- **Markers**: Synthesize `n_cells` and `pos_count` as sums, compute `pos_frac`

## Troubleshooting

### Common Issues

**"Package not found" errors:**
```bash
# R packages
Rscript scripts/install_shiny_deps.R

# Python packages
pip install -r app/panel_app/requirements.txt
```

**"Cannot find master_results.xlsx":**
- Ensure data file exists at `data/master_results.xlsx`
- Build from TSVs: `python scripts/build_master_excel.py`
- Check that data is not ignored by .gitignore

**Image viewer not loading:**
```bash
# Install Avivator viewer
./scripts/install_avivator.sh

# Verify static files
./scripts/verify_static_images.sh
```

**Trajectory analysis not available:**
```bash
# Install Bioconductor packages
INSTALL_VITESSCER=1 Rscript scripts/install_shiny_deps.R
```

**App performance issues:**
```bash
# Clear Shiny cache
./scripts/clear_shiny_cache.sh

# Reduce data size or increase bin width in UI
```

### Getting Help

1. Check application-specific READMEs:
   - [R Shiny App Documentation](app/shiny_app/README.md)
   - [Python Panel App Documentation](app/panel_app/README.md)

2. Review diagnostic output in the console/terminal

3. Run data quality scripts to identify issues:
   ```bash
   Rscript scripts/na_audit.R
   Rscript scripts/diagnose_islets.R
   ```

## Contributing

### Development Workflow

1. **Code Style**
   - R: Follow tidyverse style guide
   - Python: Follow PEP 8 conventions

2. **Data Conventions**
   - Preserve three-key system (`Case ID`, `Donor Status`, `islet_key`)
   - Use normalized region labels (`islet_core`, `islet_band`, `islet_union`)
   - Maintain donor group ordering: ND < Aab+ < T1D

3. **Testing**
   - Test with representative data
   - Verify plots and statistics output
   - Check edge cases (missing data, single donors, etc.)

4. **Documentation**
   - Update README when adding features
   - Document new scripts in their headers
   - Add inline comments for complex logic

### Project Conventions

**Normalization Options:**
- `none`: Raw values
- `global z-score`: (value - mean) / sd across all data
- `robust per-donor`: (value - median) / (MAD * 1.4826) per donor

**Donor Colors:**
- ND: `#1f77b4` (blue)
- Aab+: `#ff7f0e` (orange)
- T1D: `#d62728` (red)

**Statistical Models:**
- Global: `value ~ donor_status + islet_diam_um`
- Pairwise: t-tests on residuals from `value ~ islet_diam_um`
- Multiple testing correction: Benjamini-Hochberg (BH)

## License

Research use only. Contact the repository maintainers for licensing and collaboration inquiries.

## Citation

If you use Islet Explorer in your research, please cite this repository:

```
Islet Explorer: Interactive web applications for exploring human pancreatic islet measurements
https://github.com/smith6jt-cop/Islet-Explorer
```

## Acknowledgments

- Developed for analysis of human pancreatic islet multiplex imaging data
- Built with R Shiny and Python Panel frameworks
- Avivator/Viv viewer for OME-TIFF visualization
- Statistical analysis with R stats, broom, and scipy/statsmodels

## Contact

For questions, issues, or collaboration opportunities, please open an issue on GitHub or contact the repository maintainers.

---

**Version**: 1.0  
**Last Updated**: December 2024
