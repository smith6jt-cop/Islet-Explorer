# Islet Explorer - Panel Application

A Python-based spatial omics visualization platform for analyzing pancreatic islet data from multiplex imaging studies. Built with Panel for interactive web applications.

## Features

- **Plot Tab**: Interactive scatter plots of islet metrics by size, with distribution comparisons and AI-powered analysis assistant
- **Trajectory Tab**: UMAP visualizations and pseudotime trajectory analysis
- **Viewer Tab**: Embedded AVIVATOR viewer for OME-TIFF multiplex images
- **Statistics Tab**: ANOVA, pairwise comparisons, and effect size calculations

## Quick Start

### Installation

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt
```

### Running the Application

```bash
# Development mode
panel serve app.py --show --autoreload

# Production mode
panel serve app.py --address 0.0.0.0 --port 8080 --allow-websocket-origin=your-domain.com
```

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `ISLET_DATA_PATH` | Path to master_results.xlsx | `../../data/master_results.xlsx` |
| `LOCAL_IMAGE_ROOT` | Directory containing OME-TIFF files | `../shiny_app/www/local_images` |
| `NAVIGATOR_API_URL` | UF Navigator AI API endpoint | https://api.rc.ufl.edu/ai/ufos/v1/chat/completions |
| `NAVIGATOR_API_KEY` | API key for Navigator AI | (required for AI features) |

## Project Structure

```
panel_app/
├── app.py                 # Main application entry point
├── requirements.txt       # Python dependencies
├── README.md
├── components/
│   ├── __init__.py
│   ├── ai_assistant.py    # AI chat interface component
│   └── image_viewer.py    # AVIVATOR viewer embedding
├── utils/
│   ├── __init__.py
│   ├── data_loader.py     # Data loading and preprocessing
│   ├── plotting.py        # Plotly/HoloViews utilities
│   └── statistics.py      # Statistical analysis functions
└── assets/
    └── (static files)
```

## Data Requirements

The application expects a master_results.xlsx file with these sheets:

- **Islet_Markers**: Marker positivity data per islet region
- **Islet_Targets**: Target density and count measurements
- **Islet_Composition**: Cell type composition fractions
- **LGALS3** (optional): Additional LGALS3 marker data

## nginx Configuration

For production deployment behind nginx:

```nginx
location /islet-explorer {
    proxy_pass http://127.0.0.1:8080;
    proxy_http_version 1.1;
    proxy_set_header Upgrade $http_upgrade;
    proxy_set_header Connection "upgrade";
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
    proxy_read_timeout 86400;
}

# Static images with Range support for OME-TIFF
location /islet-explorer/images/ {
    alias /path/to/local_images/;
    add_header Accept-Ranges bytes;
    add_header Access-Control-Allow-Origin *;
}
```

## Key Differences from R Shiny Version

| Feature | R Shiny | Panel (Python) |
|---------|---------|----------------|
| Reactive system | `reactive()` / `observe()` | `param.depends()` / `pn.bind()` |
| Plots | ggplot2 + plotly | Plotly + HoloViews |
| Tables | DT::datatable | pn.widgets.Tabulator |
| Chat | Custom HTML | pn.chat.ChatInterface |
| Data loading | readxl | pandas + openpyxl |
| Statistics | broom, stats | scipy, statsmodels |

## Development

### Adding New Features

1. Add data processing to `utils/data_loader.py`
2. Add plot functions to `utils/plotting.py`
3. Add UI components as methods in `IsletExplorerApp` class
4. Wire up reactivity with `@pn.depends()` decorators

### Testing

```bash
# Run the app in development mode
panel serve app.py --show --autoreload --dev

# Check data loading
python -c "from utils.data_loader import IsletDataLoader; dl = IsletDataLoader(); print(dl.validate_file())"
```

## License

Research use only. Contact the Islet Explorer team for licensing information.
