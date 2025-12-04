#!/bin/bash
# Run script for Islet Explorer Panel Application

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default configuration
PORT="${PORT:-8080}"
HOST="${HOST:-0.0.0.0}"
MODE="${MODE:-development}"

# Activate virtual environment if it exists
if [ -d "$SCRIPT_DIR/venv" ]; then
    source "$SCRIPT_DIR/venv/bin/activate"
elif [ -d "$SCRIPT_DIR/../../venv" ]; then
    source "$SCRIPT_DIR/../../venv/bin/activate"
fi

# Set data path if not already set
export ISLET_DATA_PATH="${ISLET_DATA_PATH:-$SCRIPT_DIR/../../data/master_results.xlsx}"
export LOCAL_IMAGE_ROOT="${LOCAL_IMAGE_ROOT:-$SCRIPT_DIR/../shiny_app/www/local_images}"

echo "============================================"
echo "  Islet Explorer Panel Application"
echo "============================================"
echo "Mode: $MODE"
echo "Data path: $ISLET_DATA_PATH"
echo "Image root: $LOCAL_IMAGE_ROOT"
echo "Server: http://$HOST:$PORT"
echo "============================================"

cd "$SCRIPT_DIR"

if [ "$MODE" = "production" ]; then
    # Production mode
    panel serve app.py \
        --address "$HOST" \
        --port "$PORT" \
        --allow-websocket-origin="*" \
        --num-procs 0 \
        --log-level info
else
    # Development mode
    panel serve app.py \
        --address "$HOST" \
        --port "$PORT" \
        --show \
        --autoreload \
        --allow-websocket-origin="*" \
        --log-level debug
fi
