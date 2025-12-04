#!/bin/bash
# Run script for Islet Explorer Panel Application

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default configuration
PORT="${PORT:-8080}"
HOST="${HOST:-0.0.0.0}"

# Activate virtual environment if it exists
if [ -d "$SCRIPT_DIR/venv" ]; then
    source "$SCRIPT_DIR/venv/bin/activate"
elif [ -d "$SCRIPT_DIR/../../venv" ]; then
    source "$SCRIPT_DIR/../../venv/bin/activate"
fi

# Set data path if not already set
export ISLET_DATA_PATH="${ISLET_DATA_PATH:-$SCRIPT_DIR/../../data/master_results.xlsx}"

# Determine image directory
if [ -n "$LOCAL_IMAGE_ROOT" ] && [ -d "$LOCAL_IMAGE_ROOT" ]; then
    echo "Images: $LOCAL_IMAGE_ROOT"
elif [ -d "$SCRIPT_DIR/local_images" ]; then
    export LOCAL_IMAGE_ROOT="$SCRIPT_DIR/local_images"
    echo "Images: $LOCAL_IMAGE_ROOT (default)"
else
    echo "WARNING: No images directory found!"
    echo "  Set LOCAL_IMAGE_ROOT=/path/to/images"
fi

echo "============================================"
echo "  Islet Explorer Panel Application"
echo "============================================"
echo "Data: $ISLET_DATA_PATH"
echo "Server: http://$HOST:$PORT"
echo "============================================"

cd "$SCRIPT_DIR"

# Run using Python directly (enables CORS support)
python app.py
