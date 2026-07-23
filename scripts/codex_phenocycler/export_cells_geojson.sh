#!/usr/bin/env bash
# Batch driver: export per-cell GeoJSON for every CODEX donor, for REDSEA.
#
# One QuPath process per image, serially — each entry's data.qpdata is 500-830 MB and holds
# 142k-220k cell objects, so running them concurrently competes for RAM with no throughput gain.
# The groovy skips an image whose output already exists, so this is safely re-runnable after a
# partial failure.
#
#   bash export_cells_geojson.sh              # all 15 donors
#   bash export_cells_geojson.sh 0112 6533    # named subset (QuPath image names, zero-padded)
#
# Output: /data/codex_phenocycler/redsea_scratch/geojson/cells__<image>.geojson
set -euo pipefail

QUPATH="${QUPATH:-/home/smith6jt/QuPath/bin/QuPath}"
PROJECT="${PROJECT:-/data/QuPath_Projects/QuPath_project_IsletExplorer/project.qpproj}"
OUT_DIR="${OUT_DIR:-/data/codex_phenocycler/redsea_scratch/geojson}"
SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/export_cells_geojson.groovy"
LOG_DIR="${LOG_DIR:-${OUT_DIR}/logs}"

# QuPath image names are ZERO-PADDED (0112), unlike the unpadded donor_id (112) used everywhere
# downstream. These must stay in QuPath's vocabulary — the groovy derives the filename tag from
# the image name, and redsea.resolve_paths reconstructs the same tag from the parquet's `image`
# column.
ALL_IMAGES=(0112 6356 6450 6479 6505 6516 6521 6533 6534 6538 6548 6549 6550 6551 6563)
IMAGES=("$@")
[ ${#IMAGES[@]} -eq 0 ] && IMAGES=("${ALL_IMAGES[@]}")

for f in "$QUPATH" "$PROJECT" "$SCRIPT"; do
  [ -e "$f" ] || { echo "ERROR: not found: $f" >&2; exit 1; }
done
mkdir -p "$OUT_DIR" "$LOG_DIR"

# /home is ~97% full; per-donor GeoJSON runs 0.9-6 GB, so refuse to start on a full volume.
avail_gb=$(df -BG --output=avail "$OUT_DIR" | tail -1 | tr -dc '0-9')
echo "Output: $OUT_DIR  (${avail_gb} GB free)"
if [ "$avail_gb" -lt 100 ]; then
  echo "ERROR: <100 GB free at $OUT_DIR; per-donor GeoJSON is 0.9-6 GB and 15 donors will not fit." >&2
  exit 1
fi

failed=()
for img in "${IMAGES[@]}"; do
  echo "=== $img ==="
  # --args is the ONLY way to pass the output dir: QuPath 0.7's `script` subcommand rejects
  # bare -D system properties ("Unknown option").
  #
  # QuPath exits 0 even when the script prints an error, so the exit code alone is not a
  # success signal — the output file is checked below and drives the exit status.
  if "$QUPATH" script --project "$PROJECT" --image "$img" \
       --args "$OUT_DIR" "$SCRIPT" > "$LOG_DIR/$img.log" 2>&1; then
    tail -1 "$LOG_DIR/$img.log"
  else
    echo "  QuPath exited non-zero (see $LOG_DIR/$img.log)"; tail -3 "$LOG_DIR/$img.log"
  fi
  [ -s "$OUT_DIR/cells__${img}.geojson" ] || failed+=("$img")
done

echo
echo "=== summary ==="
for img in "${IMAGES[@]}"; do
  f="$OUT_DIR/cells__${img}.geojson"
  if [ -s "$f" ]; then printf '  %-6s %s\n' "$img" "$(du -h "$f" | cut -f1)"
  else printf '  %-6s MISSING\n' "$img"; fi
done
if [ ${#failed[@]} -gt 0 ]; then
  echo "FAILED (no output written): ${failed[*]}" >&2
  echo "Logs: $LOG_DIR" >&2
  exit 1
fi
echo "All ${#IMAGES[@]} exports present. Verify the UUID join with verify_geojson.py before REDSEA."
