#!/usr/bin/env bash
set -euo pipefail

# Convert a TIFF to OME-Zarr. Prefers Docker bioformats2raw if available; otherwise uses local JAR.

usage() {
  echo "Usage: $0 input.tif output.zarr [--tile_width 1024 --tile_height 1024]" >&2
}

if [[ $# -lt 2 ]]; then
  usage; exit 1
fi
IN="$1"; shift
OUT="$1"; shift
EXTRA_ARGS=("$@")

if [[ ! -f "$IN" ]]; then
  echo "Input not found: $IN" >&2; exit 1
fi
OUT_DIR="$(dirname "$OUT")"
mkdir -p "$OUT_DIR"

# Try Docker container (recommended)
if command -v docker >/dev/null 2>&1; then
  echo "Using Docker: ghcr.io/glencoesoftware/bioformats2raw"
  docker pull ghcr.io/glencoesoftware/bioformats2raw:latest >/dev/null 2>&1 || true
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "$(readlink -f "$(dirname "$IN")")":/in \
    -v "$(readlink -f "$OUT_DIR")":/out \
    ghcr.io/glencoesoftware/bioformats2raw:latest \
    /in/"$(basename "$IN")" /out/"$(basename "$OUT")" \
    --file_type=ome-zarr "${EXTRA_ARGS[@]}"
  echo "Created: $OUT"
  exit 0
fi

# Fallback to local JAR
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR="$SCRIPT_DIR/tools/bioformats2raw-0.8.0.jar"
if [[ ! -f "$JAR" ]]; then
  echo "bioformats2raw JAR not found at $JAR" >&2
  echo "Run: $(dirname "$0")/get_bioformats2raw.sh" >&2
  exit 1
fi
if ! command -v java >/dev/null 2>&1; then
  echo "Java not found. Install with: sudo apt-get install -y openjdk-17-jre" >&2
  exit 1
fi
echo "Using local JAR: $JAR"
java -Xms2g -Xmx8g -jar "$JAR" "$IN" "$OUT" --file_type=ome-zarr "${EXTRA_ARGS[@]}"
echo "Created: $OUT"

