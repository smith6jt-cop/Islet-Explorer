#!/usr/bin/env bash
set -euo pipefail

# Verify that the shiny-server app symlink is correct, static folders are present,
# permissions are world-readable/executable, and that HTTP Range requests work.
#
# Usage:
#   sudo bash scripts/verify_static_images.sh [APP_DIR] [APP_URL]
# Defaults:
#   APP_DIR=/srv/shiny-server/islet-explorer
#   APP_URL=http://localhost:3838/islet-explorer

APP_DIR=${1:-/srv/shiny-server/islet-explorer}
APP_URL=${2:-http://localhost:3838/islet-explorer}

echo "==> Checking app directory: $APP_DIR"
if [[ -L "$APP_DIR" ]]; then
  TARGET=$(readlink -f "$APP_DIR")
  echo "    Symlink target: $TARGET"
else
  echo "    WARNING: $APP_DIR is not a symlink (this is OK if intentional)"
  TARGET="$APP_DIR"
fi

echo "==> Checking static folders under app"
ls -ld "$TARGET" || true
ls -ld "$TARGET/www" || { echo "ERROR: Missing $TARGET/www"; exit 1; }
ls -ld "$TARGET/www/local_images" || { echo "ERROR: Missing $TARGET/www/local_images"; exit 1; }
ls -ld "$TARGET/www/avivator" || { echo "ERROR: Missing $TARGET/www/avivator"; exit 1; }

echo "==> Checking directory execute permissions to local_images (all parents must be +x)"
IMG_DIR=$(readlink -f "$TARGET/www/local_images")
P="$IMG_DIR"
while [[ "$P" != "/" ]]; do
  ls -ld "$P"
  P=$(dirname "$P")
done | tac

echo "==> Listing a few image files and permissions"
find "$IMG_DIR" -maxdepth 1 -type f \( -iname "*.ome.tif" -o -iname "*.ome.tiff" -o -iname "*.tif" -o -iname "*.tiff" \) -printf "%m %u:%g %p\n" | head -n 5 || true

SAMPLE=$(find "$IMG_DIR" -maxdepth 1 -type f \( -iname "*.ome.tif" -o -iname "*.ome.tiff" -o -iname "*.tif" -o -iname "*.tiff" \) | head -n 1 || true)
if [[ -z "$SAMPLE" ]]; then
  echo "ERROR: No image files found in $IMG_DIR"
  exit 1
fi
BASE=$(basename "$SAMPLE")
echo "==> Using sample file: $BASE"

URL_STATIC="$APP_URL/local_images/$BASE"
echo "==> Range test (should be HTTP/1.1 206 Partial Content)"
set +e
curl -sD - -o /dev/null -H 'Range: bytes=0-1' "$URL_STATIC" | egrep -i 'HTTP/|accept-ranges|content-range|content-type' || true
set -e

echo "==> If you see HTTP 200 and no Accept-Ranges, your static server isn't honoring Range."
echo "    Ensure images are served by a web server with Range support (e.g., nginx alias)."

