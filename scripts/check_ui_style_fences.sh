#!/usr/bin/env bash
# Two fences for the Modernist redesign. Run from the repo root.
#
#   ./scripts/check_ui_style_fences.sh
#
# Fence 1 — no inline skin styles.
#   An inline `style=` attribute beats any stylesheet without `!important`, so
#   every one that sets colour/background/border/shadow/radius is a hole in the
#   design system. Layout declarations (display, flex, gap, height, overflow)
#   are fine and deliberately NOT matched.
#   Expected to FAIL until Phase 3 strips them; it is a ratchet, not a gate,
#   until then. Baseline when the fence was written: 212 matching lines across
#   14 files (a line can carry more than one declaration).
#
# Fence 2 — chrome accent must never appear inside a plot.
#   Crimson #b71533 sits at deltaE 4.3 from #B2182B, the Aab+ colour in the
#   "Diverging" donor palette — below the threshold of reliable discrimination.
#   Gold #c8a02e collides with #CCBB44 (Colorblind Safe APCs) at deltaE 13.8.
#   Chrome accent stays lines and borders; a filled swatch reads as a legend key.
#
# Exit 0 = both clean.

set -uo pipefail
cd "$(dirname "$0")/.." || exit 2

FAIL=0
APP="app/shiny_app"
SRC=("$APP"/R/*.R "$APP"/app.R)

echo "== Fence 1: inline skin styles =="
HITS=$(grep -nE "style *=.*(background|linear-gradient|border-radius|box-shadow|(^|[^-])color: *#)" \
        "${SRC[@]}" 2>/dev/null)
if [ -n "$HITS" ]; then
  echo "$HITS" | awk -F: '{print $1}' | sort | uniq -c | sort -rn | sed 's/^/  /'
  echo "  TOTAL: $(echo "$HITS" | wc -l) inline skin declarations"
  FAIL=1
else
  echo "  clean"
fi

echo
echo "== Fence 2: chrome accent inside plot code =="
HITS2=$(grep -nEi 'b71533|c8a02e' "$APP"/R/mod_*_server.R "$APP"/R/drilldown_helpers.R \
         "$APP"/R/segmentation_helpers.R 2>/dev/null)
if [ -n "$HITS2" ]; then
  echo "$HITS2" | sed 's/^/  /'
  echo "  Chrome accent leaked into plot code -- use a data palette or --color-accent-700."
  FAIL=1
else
  echo "  clean"
fi

exit $FAIL
