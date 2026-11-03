#!/usr/bin/env python3
"""Register the pinned Matplotlib fonts already present in the image."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import zlib

import matplotlib


FILES = (
    "DejaVuSans.ttf",
    "DejaVuSans-Bold.ttf",
    "DejaVuSans-Oblique.ttf",
    "DejaVuSans-BoldOblique.ttf",
    "LICENSE_DEJAVU",
)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pdf_superscripts(path):
    """Check Cairo's embedded Unicode maps, not merely a nonempty PDF."""
    payload = path.read_bytes()
    maps = []
    for match in re.finditer(rb"stream\r?\n(.*?)\r?\nendstream", payload, re.S):
        stream = match.group(1)
        try:
            stream = zlib.decompress(stream)
        except zlib.error:
            pass
        if b"begincmap" in stream:
            maps.append(stream.lower())
    joined = b"\n".join(maps)
    if not all(code in joined for code in (b"<207a>", b"<207b>")):
        raise RuntimeError(f"Cairo PDF lacks both superscript Unicode mappings: {path}")
    if b"DejaVuSans" not in payload:
        raise RuntimeError(f"Cairo PDF did not embed the expected font: {path}")
    return {"file": path.name, "sha256": sha256(path), "superscript_unicode_maps": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=Path("/usr/local/share/fonts/islet-dejavu"))
    parser.add_argument("--evidence", type=Path, default=Path("/opt/islet-explorer-build/raster-fonts"))
    parser.add_argument("--validation-script", required=True, type=Path)
    args = parser.parse_args()
    source = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
    args.destination.mkdir(parents=True, exist_ok=True)
    args.evidence.mkdir(parents=True, exist_ok=False)
    records = []
    for name in FILES:
        original, installed = source / name, args.destination / name
        shutil.copyfile(original, installed)
        installed.chmod(0o644)
        if sha256(original) != sha256(installed):
            raise RuntimeError(f"Copied font differs from installed Matplotlib: {name}")
        records.append({"name": name, "source": str(original), "installed": str(installed),
                        "sha256": sha256(installed), "bytes": installed.stat().st_size})
    cache = subprocess.run(["fc-cache", "--force", str(args.destination)], check=True,
                           text=True, capture_output=True)
    subprocess.run(["Rscript", "--vanilla", str(args.validation_script), str(args.evidence)], check=True)
    validation = json.loads((args.evidence / "validation.json").read_text())
    if validation.get("result") != "pass":
        raise RuntimeError("Raster-font validation did not pass")
    pdfs = [pdf_superscripts(args.evidence / f"cairo-{face}.pdf")
            for face in ("plain", "bold", "italic", "bolditalic")]
    receipt = {
        "schema": "islet.raster_fonts/1", "result": "pass",
        "matplotlib_version": matplotlib.__version__, "matplotlib_data_path": str(source.parent.parent),
        "files": records, "fontconfig_cache_exit_code": cache.returncode,
        "installer_sha256": sha256(Path(__file__)),
        "validation_script_sha256": sha256(args.validation_script),
        "raster_validation_sha256": sha256(args.evidence / "validation.json"),
        "cairo_pdf_validation": pdfs,
        "artifacts": [{"file": path.name, "sha256": sha256(path), "bytes": path.stat().st_size}
                      for path in sorted(args.evidence.iterdir()) if path.is_file()],
        "dependency_changes": False,
    }
    (args.evidence / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"result": "pass", "receipt": str(args.evidence / "receipt.json")}))


if __name__ == "__main__":
    main()
