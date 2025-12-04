#!/bin/bash
# Run script for Islet Explorer Panel Application

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHINY_WWW="$SCRIPT_DIR/../shiny_app/www"

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
export LOCAL_IMAGE_ROOT="${LOCAL_IMAGE_ROOT:-$SHINY_WWW/local_images}"

# Build static dirs arguments
STATIC_ARGS=""
if [ -d "$SHINY_WWW/avivator" ]; then
    STATIC_ARGS="$STATIC_ARGS --static-dirs avivator=$SHINY_WWW/avivator"
    echo "AVIVATOR: $SHINY_WWW/avivator"
fi
if [ -d "$LOCAL_IMAGE_ROOT" ]; then
    STATIC_ARGS="$STATIC_ARGS --static-dirs images=$LOCAL_IMAGE_ROOT"
    echo "Images: $LOCAL_IMAGE_ROOT"
fi

echo "============================================"
echo "  Islet Explorer Panel Application"
echo "============================================"
echo "Mode: $MODE"
echo "Data path: $ISLET_DATA_PATH"
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
        --log-level info \
        $STATIC_ARGS
else
    # Development mode (no --show to avoid browser issues)
    panel serve app.py \
        --address "$HOST" \
        --port "$PORT" \
        --autoreload \
        --allow-websocket-origin="*" \
        --log-level debug \
        $STATIC_ARGS
fi
