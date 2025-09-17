#!/usr/bin/env python3
"""
Annotate OME-TIFF files with channel names and colors (in-place) by patching the OME-XML.

- Input: a CSV with columns: filename,channel_index,channel_name,color_hex
  * filename: basename or relative path to the OME-TIFF under a data root
  * channel_index: 0-based channel index within the OME image
  * channel_name: desired channel name (e.g., CD45)
  * color_hex: hex RGB like #FF0000 (optional); if blank, color entry is skipped

- Usage:
    python scripts/annotate_ome_tiff_channels.py \
        --data-root local_images \
        --csv channel_map.csv

Requires: tifftools (pip install tifftools), lxml
This updates only the OME-XML metadata; pixel data is untouched.
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

try:
    import tifftools as tt
except Exception as e:
    print("ERROR: tifftools is required: pip install tifftools", file=sys.stderr)
    sys.exit(1)

try:
    from lxml import etree
except Exception:
    print("ERROR: lxml is required: pip install lxml", file=sys.stderr)
    sys.exit(1)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--data-root", required=True, help="Root folder containing OME-TIFF files")
    p.add_argument("--csv", required=True, help="CSV mapping: filename,channel_index,channel_name,color_hex")
    p.add_argument("--dry-run", action="store_true", help="Print planned changes without writing")
    return p.parse_args()


def load_mapping(csv_path):
    mapping = defaultdict(dict)
    with open(csv_path, newline="") as f:
        rdr = csv.DictReader(f)
        required = {"filename", "channel_index", "channel_name"}
        if not required.issubset(rdr.fieldnames):
            raise ValueError(f"CSV must include columns: {', '.join(sorted(required))}")
        for row in rdr:
            fname = row["filename"].strip()
            if not fname:
                continue
            try:
                ci = int(row["channel_index"])
            except Exception:
                continue
            cname = row["channel_name"].strip()
            chex = (row.get("color_hex") or "").strip()
            mapping[fname][ci] = {"name": cname, "color": chex}
    return mapping


def patch_ome_xml(xml_bytes, chan_map):
    parser = etree.XMLParser(remove_blank_text=True)
    root = etree.fromstring(xml_bytes, parser=parser)
    ns = {"ome": root.nsmap.get(None) or "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    # For each Image/Channels, update names and colors
    for image in root.xpath("//ome:Image", namespaces=ns):
        pixels = image.find("{*}Pixels")
        if pixels is None:
            continue
        channels = [c for c in pixels.findall("{*}Channel")]
        for idx, ch in enumerate(channels):
            if idx in chan_map:
                entry = chan_map[idx]
                if entry.get("name"):
                    ch.set("Name", entry["name"])  # Name attribute
                color = entry.get("color")
                if color and color.startswith("#") and len(color) == 7:
                    # OME-XML Channel Color is commonly a 24-bit RGB integer (RRGGBB) in decimal.
                    try:
                        r = int(color[1:3], 16)
                        g = int(color[3:5], 16)
                        b = int(color[5:7], 16)
                        rgb_int = (r << 16) | (g << 8) | b
                        ch.set("Color", str(rgb_int))
                    except Exception:
                        pass
    return etree.tostring(root, pretty_print=True, xml_declaration=True, encoding="UTF-8")


def annotate_tiff(tif_path, chan_map, dry_run=False):
    if not os.path.exists(tif_path):
        print(f"SKIP (missing): {tif_path}")
        return
    info = tt.read_tiff(tif_path)
    # Find OME-XML in ImageDescription of first IFD
    ifds = info.get("ifds", [])
    if not ifds:
        print(f"SKIP (no IFDs): {tif_path}")
        return
    # Search ImageDescription
    desc_tag = 270
    found = None
    for ifd in ifds:
        tags = ifd.get("tags", {})
        if desc_tag in tags:
            val = tags[desc_tag].get("data")
            if isinstance(val, (bytes, bytearray)):
                txt = val
            else:
                txt = str(val).encode("utf-8")
            if b"OME" in txt and b"Pixels" in txt:
                found = (ifd, tags[desc_tag])
                break
    if not found:
        print(f"SKIP (no OME-XML found): {tif_path}")
        return
    ifd, tag = found
    try:
        new_xml = patch_ome_xml(tag["data"], chan_map)
    except Exception as e:
        print(f"ERROR parsing or patching OME-XML in {tif_path}: {e}", file=sys.stderr)
        return
    if dry_run:
        print(f"DRY-RUN would update: {tif_path}")
        return
    # Update tag and write back
    tag["data"] = new_xml
    tt.write_tiff(tif_path, info)
    print(f"UPDATED: {tif_path}")


def main():
    args = parse_args()
    mapping = load_mapping(args.csv)
    if not mapping:
        print("No mapping entries loaded; check CSV.")
        return
    for fname, chan_map in mapping.items():
        tif_path = os.path.join(args.data_root, fname)
        annotate_tiff(tif_path, chan_map, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
