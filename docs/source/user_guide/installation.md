# Installation Guide

## System Requirements

### Minimum Requirements

- **Operating System**: Linux, macOS, or Windows with WSL2
- **Memory**: 8 GB RAM (16 GB recommended)
- **Storage**: 20 GB free space (for data and cache files)
- **R**: Version 4.0.0 or higher
- **Python**: Version 3.6 or higher (for GeoJSON preprocessing)

### Recommended Specifications

- **Memory**: 32 GB RAM
- **Storage**: 50 GB+ free space (SSD recommended)
- **CPU**: Multi-core processor (4+ cores)
- **Display**: 1920×1080 or higher resolution

## Installing R Dependencies

The application requires several R packages. Install them using the provided script:

```bash
cd /path/to/Islet-Explorer
Rscript scripts/install_shiny_deps.R
```

### Manual Installation

If the script fails, install packages manually in R:

```r
# Core packages
install.packages(c(
  "shiny", "shinyjs", "readxl", "dplyr", "stringr",
  "tidyr", "ggplot2", "plotly", "broom", "jsonlite",
  "base64enc"
))

# Optional: For trajectory analysis
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("anndata")
```

## Installing Python Dependencies

Python is required for GeoJSON preprocessing:

```bash
python3 --version  # Verify Python 3.6+ is installed

# No additional packages needed - uses stdlib only!
```

## Setting Up GeoJSON Overlay System

The overlay system requires one-time preprocessing:

### 1. Verify Data Structure

```bash
ls gson/*.geojson.gz
# Should show: 0112.geojson.gz, 6356.geojson.gz, etc.
```

### 2. Run Preprocessing

```bash
./scripts/preprocess_all_geojson.sh
```

**Time Required**: ~1-2 hours for 15 files  
**Disk Space**: Adds ~3 GB in `data/gson_cache/`

### 3. Verify Preprocessing

```bash
ls -lh data/gson_cache/*.pkl
# Should show 15 files, each ~200 MB
```

## Installing Avivator Viewer

The embedded image viewer requires a static build of Avivator:

```bash
./scripts/install_avivator.sh
```

**Requirements**: Node.js ≥ 18 and pnpm (optional, for building from source)

The script will:
1. Try to download pre-built bundle from public site
2. If Node.js available, build from source
3. Install to `app/shiny_app/www/avivator/`

## Server Deployment (Optional)

### Using Shiny Server

1. Install Shiny Server:

```bash
# Ubuntu/Debian
sudo apt-get install gdebi-core
wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.20.1002-amd64.deb
sudo gdebi shiny-server-1.5.20.1002-amd64.deb
```

2. Deploy app:

```bash
sudo cp -r /path/to/Islet-Explorer/app/shiny_app /srv/shiny-server/islet-explorer
sudo systemctl restart shiny-server
```

### Using Docker

Create `Dockerfile`:

```dockerfile
FROM rocker/shiny:4.3.0

# Install system dependencies
RUN apt-get update && apt-get install -y \\
    python3 \\
    python3-pip \\
    && rm -rf /var/lib/apt/lists/*

# Copy app
COPY app/shiny_app /srv/shiny-server/islet-explorer
COPY data /srv/shiny-server/data
COPY scripts /srv/shiny-server/scripts

# Install R dependencies
RUN Rscript /srv/shiny-server/scripts/install_shiny_deps.R

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
```

Build and run:

```bash
docker build -t islet-explorer .
docker run -p 3838:3838 islet-explorer
```

## Reverse Proxy Setup (nginx)

For production deployment with SSL:

```nginx
server {
    listen 80;
    server_name your-domain.com;
    return 301 https://$server_name$request_uri;
}

server {
    listen 443 ssl;
    server_name your-domain.com;
    
    ssl_certificate /path/to/cert.pem;
    ssl_certificate_key /path/to/key.pem;
    
    location / {
        proxy_pass http://localhost:3838;
        proxy_redirect http://localhost:3838/ https://$host/;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_read_timeout 20d;
        proxy_buffering off;
    }
}
```

## Verification

Test your installation:

```bash
# Start the app
R -q -e "shiny::runApp('app/shiny_app', launch.browser=FALSE)"

# In another terminal, test accessibility
curl http://localhost:XXXX
```

Expected output: HTML content from the Shiny app

## Troubleshooting Installation

### R Package Installation Fails

**Problem**: Package compilation errors

**Solution**:
```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
```

### GeoJSON Preprocessing Fails

**Problem**: "Out of memory" error

**Solution**: Process files one at a time:
```bash
for f in gson/*.geojson.gz; do
    python3 scripts/preprocess_geojson.py "$f" "data/gson_cache/$(basename ${f%.geojson.gz})_simplified.pkl"
done
```

### Avivator Not Found

**Problem**: Viewer shows "Static build not found"

**Solution**:
1. Verify `app/shiny_app/www/avivator/index.html` exists
2. Re-run `./scripts/install_avivator.sh`
3. Or manually download from https://avivator.gehlenborglab.org/

## Next Steps

After installation, proceed to:
- [Getting Started](getting_started.md) - Learn the basic workflow
- [Plot Tab Guide](plot_tab.md) - Explore data visualization
- [Viewer Tab Guide](viewer_tab.md) - Use the image viewer
