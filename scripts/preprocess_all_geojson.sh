#!/bin/bash
# Batch preprocess all GeoJSON files for efficient viewer access
# This creates simplified, spatially-indexed pickle files

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
GSON_DIR="$PROJECT_ROOT/gson"
OUTPUT_DIR="$PROJECT_ROOT/data/gson_cache"

echo "[PREPROCESS] Starting batch GeoJSON preprocessing"
echo "[PREPROCESS] Input:  $GSON_DIR"
echo "[PREPROCESS] Output: $OUTPUT_DIR"

if [ ! -d "$GSON_DIR" ]; then
    echo "[ERROR] GeoJSON directory not found: $GSON_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Check for Python 3
if ! command -v python3 &> /dev/null; then
    echo "[ERROR] python3 is required but not found"
    exit 1
fi

# Run Python preprocessing script
python3 "$SCRIPT_DIR/preprocess_geojson.py" --batch "$GSON_DIR" "$OUTPUT_DIR"

echo "[PREPROCESS] Complete! Processed files are in: $OUTPUT_DIR"
echo "[PREPROCESS] You can now run the Shiny app with efficient overlay support"
