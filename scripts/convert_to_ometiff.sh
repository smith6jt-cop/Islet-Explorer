#!/usr/bin/env bash
set -euo pipefail

# Convert a TIFF to pyramidal OME-TIFF using bioformats2raw + raw2ometiff (Docker or local JARs).

usage() {
  echo "Usage: $0 input.tif output.ome.tif [--tile_width 1024 --tile_height 1024 --compression LZW]" >&2
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
TMPDIR="$(mktemp -d)"
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

RAW_DIR="$TMPDIR/raw"

if command -v docker >/dev/null 2>&1; then
  echo "Using Docker: ghcr.io/glencoesoftware/bioformats2raw + raw2ometiff"
  docker pull ghcr.io/glencoesoftware/bioformats2raw:latest >/dev/null 2>&1 || true
  docker pull ghcr.io/glencoesoftware/raw2ometiff:latest >/dev/null 2>&1 || true
  docker run --rm -u "$(id -u):$(id -g)" \
    -v "$(readlink -f "$(dirname "$IN")")":/in \
    -v "$TMPDIR":/work ghcr.io/glencoesoftware/bioformats2raw:latest \
    /in/"$(basename "$IN")" /work/raw
  docker run --rm -u "$(id -u):$(id -g)" \
    -v "$TMPDIR":/work \
    -v "$(readlink -f "$(dirname "$OUT")")":/out \
    ghcr.io/glencoesoftware/raw2ometiff:latest \
    /work/raw /out/"$(basename "$OUT")" "${EXTRA_ARGS[@]}"
  echo "Created: $OUT"
  exit 0
fi

# Fallback to local JARs
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BF2RAW_JAR="$SCRIPT_DIR/tools/bioformats2raw-0.8.0.jar"
R2OME_JAR="$SCRIPT_DIR/tools/raw2ometiff-0.8.0.jar"
if [[ ! -f "$BF2RAW_JAR" || ! -f "$R2OME_JAR" ]]; then
  echo "Required JARs not found in $SCRIPT_DIR/tools" >&2
  echo "Run: $(dirname "$0")/get_bioformats2raw.sh" >&2
  exit 1
fi
if ! command -v java >/dev/null 2>&1; then
  echo "Java not found. Install with: sudo apt-get install -y openjdk-17-jre" >&2
  exit 1
fi
echo "Using local JARs"
java -Xms2g -Xmx8g -jar "$BF2RAW_JAR" "$IN" "$RAW_DIR"
java -Xms2g -Xmx8g -jar "$R2OME_JAR" "$RAW_DIR" "$OUT" "${EXTRA_ARGS[@]}"
echo "Created: $OUT"

