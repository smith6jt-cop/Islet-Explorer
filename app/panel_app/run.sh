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

# Determine image directory
if [ -n "$LOCAL_IMAGE_ROOT" ] && [ -d "$LOCAL_IMAGE_ROOT" ]; then
    IMAGE_DIR="$LOCAL_IMAGE_ROOT"
    echo "Images: $IMAGE_DIR (from LOCAL_IMAGE_ROOT)"
elif [ -d "$SCRIPT_DIR/local_images" ]; then
    IMAGE_DIR="$SCRIPT_DIR/local_images"
    export LOCAL_IMAGE_ROOT="$IMAGE_DIR"
    echo "Images: $IMAGE_DIR (default)"
else
    echo "WARNING: No images directory found!"
    echo "  Set LOCAL_IMAGE_ROOT=/path/to/images"
    echo "  Or create: $SCRIPT_DIR/local_images"
    IMAGE_DIR=""
fi

# Build static dirs argument for images only
STATIC_ARGS=""
if [ -n "$IMAGE_DIR" ]; then
    STATIC_ARGS="--static-dirs images=$IMAGE_DIR"
fi

echo "============================================"
echo "  Islet Explorer Panel Application"
echo "============================================"
echo "Mode: $MODE"
echo "Data: $ISLET_DATA_PATH"
echo "Server: http://$HOST:$PORT"
echo "Static: $STATIC_ARGS"
echo "============================================"

cd "$SCRIPT_DIR"

if [ "$MODE" = "production" ]; then
    panel serve app.py \
        --address "$HOST" \
        --port "$PORT" \
        --allow-websocket-origin="*" \
        --num-procs 0 \
        --log-level info \
        $STATIC_ARGS
else
    panel serve app.py \
        --address "$HOST" \
        --port "$PORT" \
        --autoreload \
        --allow-websocket-origin="*" \
        --log-level debug \
        $STATIC_ARGS
fi
