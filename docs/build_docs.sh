#!/bin/bash
# Build and optionally serve Sphinx documentation

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCS_DIR="$SCRIPT_DIR"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VENV_DIR="$PROJECT_ROOT/docs_env"

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}[DOCS]${NC} Building Islet-Explorer documentation..."

# Check if virtual environment exists
if [ ! -d "$VENV_DIR" ]; then
    echo -e "${BLUE}[DOCS]${NC} Creating virtual environment..."
    python3 -m venv "$VENV_DIR"
    source "$VENV_DIR/bin/activate"
    echo -e "${BLUE}[DOCS]${NC} Installing dependencies..."
    pip install --quiet --upgrade pip
    pip install --quiet -r "$DOCS_DIR/requirements.txt"
else
    source "$VENV_DIR/bin/activate"
fi

# Clean previous build
echo -e "${BLUE}[DOCS]${NC} Cleaning previous build..."
cd "$DOCS_DIR"
make clean > /dev/null 2>&1

# Build HTML documentation
echo -e "${BLUE}[DOCS]${NC} Building HTML documentation..."
make html

if [ $? -eq 0 ]; then
    echo -e "${GREEN}[DOCS]${NC} Documentation built successfully!"
    echo -e "${GREEN}[DOCS]${NC} Output: $DOCS_DIR/build/html/index.html"
    
    # Ask to serve
    if [ "$1" != "--no-serve" ]; then
        echo ""
        read -p "Start local server? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            echo -e "${BLUE}[DOCS]${NC} Starting documentation server..."
            echo -e "${GREEN}[DOCS]${NC} Open: http://localhost:8000"
            echo -e "${BLUE}[DOCS]${NC} Press Ctrl+C to stop"
            cd build/html
            python3 -m http.server 8000
        fi
    fi
else
    echo -e "${RED}[DOCS]${NC} Build failed! Check errors above."
    exit 1
fi
