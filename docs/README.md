# Islet-Explorer Documentation

This directory contains the comprehensive Sphinx documentation for Islet-Explorer.

## Building the Documentation

### Prerequisites

```bash
# Create virtual environment (first time only)
python3 -m venv docs_env

# Activate environment
source docs_env/bin/activate

# Install dependencies (first time only)
pip install sphinx sphinx-rtd-theme myst-parser sphinx-autodoc-typehints
```

### Build HTML Documentation

```bash
# From project root
cd docs
source ../docs_env/bin/activate
make html
```

Output will be in `build/html/`.

### View Locally

```bash
# Start local server
cd build/html
python3 -m http.server 8000

# Open browser to http://localhost:8000
```

### Clean Build

```bash
make clean
make html
```

## Documentation Structure

```
docs/
├── source/
│   ├── index.rst                  # Main landing page
│   ├── conf.py                    # Sphinx configuration
│   │
│   ├── user_guide/               # User-facing documentation
│   │   ├── installation.md       # Setup instructions
│   │   ├── getting_started.md    # First-time user guide
│   │   ├── plot_tab.md           # Plot features
│   │   ├── statistics_tab.md     # Statistical analysis
│   │   ├── trajectory_tab.md     # Pseudotime analysis
│   │   └── viewer_tab.md         # Image viewer
│   │
│   ├── technical/                # Technical documentation
│   │   ├── architecture.md       # System design
│   │   ├── data_format.md        # Data specifications
│   │   ├── geojson_overlay.md    # Overlay system (COMPLETE)
│   │   ├── preprocessing.md      # Data preprocessing
│   │   └── api_reference.md      # Function reference
│   │
│   ├── development/              # Developer documentation
│   │   ├── contributing.md       # Contribution guidelines
│   │   ├── testing.md            # Testing procedures
│   │   └── deployment.md         # Deployment guide
│   │
│   └── reference/                # Reference materials
│       ├── glossary.md           # Term definitions
│       ├── troubleshooting.md    # Common issues
│       ├── faq.md                # Frequently asked questions
│       └── changelog.md          # Version history
│
├── build/                         # Generated documentation
│   └── html/                      # HTML output
│
├── Makefile                       # Build commands
└── README.md                      # This file
```

## Documentation Sections

### ✅ Complete

- **Main Index** - Landing page with overview
- **Installation Guide** - Complete setup instructions
- **Getting Started** - First-time user walkthrough
- **GeoJSON Overlay System** - Comprehensive technical documentation

### 🚧 In Progress

- **User Guide** - Tab-specific guides (stubs created)
- **Technical Docs** - Architecture and data format
- **Reference** - Glossary, FAQ, troubleshooting

## Writing Documentation

### File Format

Use Markdown (`.md`) or reStructuredText (`.rst`). Markdown is preferred for new docs.

### MyST Markdown Extensions

The docs support enhanced Markdown via MyST:

````markdown
# Heading

## Subheading

**Bold** and *italic* text

- Bullet list
- Item 2

1. Numbered list
2. Item 2

```python
# Code block with syntax highlighting
def hello():
    print("world")
```

:::{note}
This is a note admonition
:::

:::{warning}
This is a warning
:::

[Link text](https://example.com)

![Image alt](path/to/image.png)
````

### Adding New Pages

1. Create `.md` or `.rst` file in appropriate directory
2. Add to `toctree` in `source/index.rst`:

```rst
.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/installation
   user_guide/your_new_page
```

3. Rebuild: `make html`

## Deploying to Read the Docs

### Setup (One-time)

1. Create account at https://readthedocs.org/
2. Import project from GitHub
3. Configure webhook (automatic builds on push)

### Configuration

The project includes `.readthedocs.yaml` for automatic deployment:

```yaml
version: 2

build:
  os: ubuntu-22.04
  tools:
    python: "3.12"

sphinx:
  configuration: docs/source/conf.py

python:
  install:
    - requirements: docs/requirements.txt
```

### Manual Deployment

```bash
# Commit documentation changes
git add docs/
git commit -m "Update documentation"
git push origin master

# Read the Docs will automatically rebuild
```

## Best Practices

### Writing Style

- **Be concise**: Get to the point quickly
- **Use examples**: Show, don't just tell
- **Include code**: Runnable examples are best
- **Add images**: Screenshots for UI documentation
- **Link heavily**: Cross-reference related content

### Organization

- **User-first**: Most common tasks in User Guide
- **Progressive detail**: Simple → Complex
- **Searchable**: Use descriptive headings and terms
- **Maintainable**: Keep related content together

### Code Examples

```python
# Good: Complete, runnable example
import pandas as pd

df = pd.read_excel("data/master_results.xlsx", sheet_name="Islet_Markers")
filtered = df[df['donor_status'] == 'T1D']
print(f"Found {len(filtered)} T1D samples")
```

```r
# Good: R example with context
library(dplyr)

# Load and filter islet data
targets <- read_excel("data/master_results.xlsx", sheet = "Islet_Targets")
large_islets <- targets %>% filter(islet_diam_um > 100)
```

### Images

Place in `source/_static/images/`:

```markdown
![Plot tab interface](/_static/images/plot_tab_screenshot.png)
```

## Continuous Improvement

### Feedback Loop

1. User asks question → Check if docs cover it
2. If not, add to FAQ or relevant section
3. If unclear, rewrite or add examples
4. Update troubleshooting for common errors

### Version Control

- Document new features in Changelog
- Update Getting Started for workflow changes
- Archive old versions if significant changes

## Troubleshooting Build Errors

### "Module not found"

```bash
# Reinstall dependencies
pip install --upgrade -r requirements.txt
```

### "Lexing error" warnings

Ignore - these don't affect output. Or fix syntax in code blocks.

### "Document not in toctree"

Either add to `index.rst` or add to `exclude_patterns` in `conf.py`.

### Theme issues

```bash
# Rebuild clean
make clean
make html
```

## Resources

- [Sphinx Documentation](https://www.sphinx-doc.org/)
- [MyST Parser Guide](https://myst-parser.readthedocs.io/)
- [Read the Docs](https://docs.readthedocs.io/)
- [reStructuredText Primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
