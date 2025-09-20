#!/usr/bin/env python3
"""
Generate offsets.json sidecars for OME‑TIFF files so Avivator/Viv can load them efficiently.

Usage:
  python scripts/make_offsets_json.py local_images/*.ome.tif*

The script writes a JSON file next to each TIFF: <name>.offsets.json
Requires: tifffile (pip install tifffile)
"""
import json
import os
import sys
from typing import Dict, Any

try:
    import tifffile
except Exception as e:
    print("ERROR: tifffile is required. Install with: pip install tifffile", file=sys.stderr)
    sys.exit(1)


def compute_offsets(path: str) -> Dict[str, Any]:
    info = {"offsets": []}
    with tifffile.TiffFile(path) as tif:
        for page in tif.pages:
            # Record byte offsets for strips or tiles
            o = {}
            if page.dataoffsets is not None:
                o["dataoffsets"] = list(map(int, page.dataoffsets))
            if page.databytecounts is not None:
                o["databytecounts"] = list(map(int, page.databytecounts))
            # Useful metadata
            o["shape"] = tuple(int(x) for x in getattr(page, 'shape', ()) or ())
            o["dtype"] = str(getattr(page, 'dtype', ''))
            info["offsets"].append(o)
    return info


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    rc = 0
    for pattern in argv[1:]:
        for p in sorted([pattern] if os.path.isfile(pattern) else []):
            pass
        # Expand globs manually (to support shell on Windows, etc.)
        import glob
        for path in glob.glob(pattern):
            if not os.path.isfile(path):
                continue
            base = os.path.splitext(path)[0]
            out = base + ".offsets.json"
            try:
                print(f"Computing offsets for {path} -> {out}")
                data = compute_offsets(path)
                with open(out, "w") as f:
                    json.dump(data, f)
            except Exception as e:
                print(f"ERROR on {path}: {e}", file=sys.stderr)
                rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

