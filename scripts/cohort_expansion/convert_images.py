#!/usr/bin/env python3
"""Convert the selected cohort expansion into lossless viewer OME-TIFFs.

The full-resolution pixels are retained exactly. Five SubIFD levels use repeated
2-by-2 integer-floor averaging, matching the existing viewer pyramids. Source
channel labels, explicit name/index mappings, calibration and checksums are kept
in conversion reports. KRT8 is deliberately not relabelled as KRT8-18.

Run with the scvi-env Python environment. Completed images are checkpoints;
rerunning validates their source/output hashes before skipping them. Incomplete
files are retained with a timestamp before a new conversion is attempted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

import imagecodecs
import numpy as np
import tifffile


REPO = Path(__file__).resolve().parents[2]
SOURCE_ROOT = Path("/data/QuPath_Projects/QuPath_project_IsletExplorer/Images/New")
OUTPUT_ROOT = REPO / "data/development/expanded-20260911/images"
SELECTION = {
    "6090": ("6090AAbtest31_35p_led.tif", "Aab_6090"),
    "6232": ("CX_PAN-6232-ND_test33_35p_led.tif", "ND_6232"),
    "6539": ("PAN-6539-ND test 34.tif", "ND_6539"),
    "6558": ("PAN-6558-AAb_test35_35p-LED.tif", "Aab_6558"),
    "6520": ("CX_PANT6520-01A_T1D_Test22_34p-16LED50.tif", "T1D_6520"),
    "78": ("CX_PANT-HPAP-078-T1D_test28_34p_16LED50.tif", "T1D_0078"),
}
CHANNEL_ORDER = [
    "DAPI", "CD31", "CD8a", "CD3e", "SMA", "Ki67", "CD4", "CD34", "HLADR",
    "PDPN", "panCK", "ECAD", "CD163", "SST", "ColIV", "VIM", "CD20", "LGALS3",
    "B3TUBB", "GCG", "KRT14", "GAP43", "CD35", "CHGA", "PGP9_5", "INS", "CD44",
    "CD45", "NaKATPase", "BCatenin", "CD68", "BActin", "CK19", "epCAM", "KRT8_18",
    "KRT8", "Ghrelin",
]
ALIASES = {
    "DAPI-01": "DAPI", "COL4A1": "ColIV", "ATP1A1": "NaKATPase",
    "CTNNB1": "BCatenin", "ACTB": "BActin", "KRT19": "CK19",
    "PGP9.5": "PGP9_5", "KRT8-18": "KRT8_18",
}
TILE = 512
LEVELS = 6
FORMAT_VERSION = 1
CONVERTER_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def log(message: str) -> None:
    print(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), message, flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def canonical_channel(label: str) -> str:
    name = re.sub(r"^reg\d+_cyc\d+_ch\d+_", "", label.strip())
    name = re.sub(r"\.tiff?$", "", name, flags=re.IGNORECASE)
    return ALIASES.get(name, name)


def inspect_source(path: Path) -> dict:
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        if len(tif.series) != 1 or series.axes != "CYX" or series.dtype != np.dtype("uint16"):
            raise ValueError(f"Expected one uint16 CYX image: {path}")
        metadata = tif.imagej_metadata or {}
        labels = metadata.get("Labels", [])
        if len(labels) != series.shape[0]:
            raise ValueError(f"Missing or ambiguous channel labels: {path}")
        names = [canonical_channel(label) for label in labels]
        if len(names) != len(set(names)):
            raise ValueError(f"Channel aliases produce duplicate markers: {names}")
        missing = sorted(set(CHANNEL_ORDER[:33]) - set(names))
        if missing:
            raise ValueError(f"Required core channels absent in {path.name}: {missing}")
        unknown = sorted(set(names) - set(CHANNEL_ORDER))
        if unknown:
            raise ValueError(f"Unrecognized channels require explicit review: {unknown}")
        if metadata.get("unit") != "mm":
            raise ValueError(f"Expected explicit ImageJ millimeter calibration: {path}")
        pixel_um = {}
        for axis in ("X", "Y"):
            numerator, denominator = tif.pages[0].tags[f"{axis}Resolution"].value
            if numerator <= 0 or denominator <= 0:
                raise ValueError(f"Invalid {axis} calibration: {path}")
            pixel_um[axis] = 1000.0 * denominator / numerator
        order = [names.index(name) for name in CHANNEL_ORDER if name in names]
        return {
            "source": str(path.resolve()), "source_bytes": path.stat().st_size,
            "source_mtime_ns": path.stat().st_mtime_ns,
            "shape_cyx": list(series.shape), "source_byteorder": tif.byteorder,
            "pixel_size_um": pixel_um,
            "channel_mapping": [
                {"output_index": i, "source_index": j, "source_label": labels[j],
                 "channel_name": names[j]}
                for i, j in enumerate(order)
            ],
            "warnings": (["KRT8 retained separately; equivalence to KRT8-18 is unconfirmed."]
                         if "KRT8" in names else []),
            "source_imagej_description": tif.pages[0].description,
        }


def reduce_half(plane: np.ndarray) -> np.ndarray:
    """Floor the 2x2 mean; an unpaired edge appears only at higher resolution."""
    height, width = plane.shape
    height -= height % 2
    width -= width % 2
    total = plane[:height:2, :width:2].astype(np.uint32)
    total += plane[1:height:2, :width:2]
    total += plane[:height:2, 1:width:2]
    total += plane[1:height:2, 1:width:2]
    total //= 4
    return total.astype(np.uint16)


def pyramid_plane(source: np.ndarray, source_index: int, level: int) -> np.ndarray:
    plane = source[source_index]
    for _ in range(level):
        plane = reduce_half(plane)
    return plane


def iter_tiles(source: np.ndarray, mapping: list[dict], level: int):
    for entry in mapping:
        plane = pyramid_plane(source, entry["source_index"], level)
        for y in range(0, plane.shape[0], TILE):
            for x in range(0, plane.shape[1], TILE):
                tile = np.zeros((TILE, TILE), dtype=np.uint16)
                view = plane[y:y + TILE, x:x + TILE]
                tile[:view.shape[0], :view.shape[1]] = view
                yield tile


def convert(source_path: Path, partial: Path, info: dict, threads: int) -> None:
    source = tifffile.memmap(source_path, mode="r")
    try:
        metadata = {
            "axes": "CYX", "Name": partial.name.removesuffix(".partial"),
            "PhysicalSizeX": info["pixel_size_um"]["X"], "PhysicalSizeXUnit": "µm",
            "PhysicalSizeY": info["pixel_size_um"]["Y"], "PhysicalSizeYUnit": "µm",
            "Channel": {"Name": [c["channel_name"] for c in info["channel_mapping"]]},
        }
        with tifffile.TiffWriter(partial, bigtiff=True, ome=True, byteorder="<") as writer:
            for level in range(LEVELS):
                scale = 2 ** level
                shape = (source.shape[0], source.shape[1] // scale, source.shape[2] // scale)
                if min(shape[1:]) < 1:
                    raise ValueError("Image is too small for six pyramid levels")
                log(f"{source_path.name}: writing level {level} {shape}")
                writer.write(
                    iter_tiles(source, info["channel_mapping"], level), shape=shape,
                    dtype=np.uint16, tile=(TILE, TILE), compression="lzw",
                    photometric="minisblack", subifds=LEVELS - 1 if level == 0 else None,
                    subfiletype=0 if level == 0 else 1,
                    metadata=metadata if level == 0 else None,
                    resolution=(10000 / info["pixel_size_um"]["X"] / scale,
                                10000 / info["pixel_size_um"]["Y"] / scale),
                    resolutionunit="CENTIMETER", maxworkers=threads,
                    software=f"tifffile {tifffile.__version__}; Islet-Explorer image conversion v{FORMAT_VERSION}",
                )
    finally:
        del source


def validate_pixels(source_path: Path, output: Path, info: dict, threads: int) -> dict:
    source = tifffile.memmap(source_path, mode="r")
    checked_pixels = 0
    channel_hashes = []
    try:
        with tifffile.TiffFile(output) as tif:
            if not tif.is_ome or not tif.is_bigtiff or len(tif.series) != 1:
                raise ValueError("Output must be a single-series OME BigTIFF")
            levels = tif.series[0].levels
            if len(levels) != LEVELS:
                raise ValueError(f"Expected {LEVELS} pyramid levels, got {len(levels)}")
            pixels = ET.fromstring(tif.ome_metadata).find(".//{*}Pixels")
            names = [c.attrib["Name"] for c in pixels.findall("{*}Channel")]
            if names != [c["channel_name"] for c in info["channel_mapping"]]:
                raise ValueError("OME channel names/order do not match source mapping")
            for axis in ("X", "Y"):
                if pixels.attrib[f"PhysicalSize{axis}Unit"] != "µm" or not np.isclose(
                    float(pixels.attrib[f"PhysicalSize{axis}"]), info["pixel_size_um"][axis],
                    rtol=1e-12, atol=0,
                ):
                    raise ValueError(f"OME {axis} calibration differs from source")
            for level, series in enumerate(levels):
                expected_shape = (source.shape[0], source.shape[1] // 2 ** level,
                                  source.shape[2] // 2 ** level)
                if series.shape != expected_shape or len(series.pages) != source.shape[0]:
                    raise ValueError(f"Unexpected level {level} dimensions or channel count")
                log(f"{source_path.name}: validating every pixel at level {level}")
                for entry, page in zip(info["channel_mapping"], series.pages):
                    tags = page.keyframe
                    if int(tags.compression) != 5 or (tags.tilelength, tags.tilewidth) != (TILE, TILE):
                        raise ValueError(f"Unexpected compression or tile size at level {level}")
                    decoded = page.asarray(maxworkers=threads)
                    expected = pyramid_plane(source, entry["source_index"], level)
                    if not np.array_equal(decoded, expected):
                        raise ValueError(f"Pixel mismatch: level {level}, {entry['channel_name']}")
                    checked_pixels += decoded.size
                    if level == 0:
                        channel_hashes.append({
                            "channel_name": entry["channel_name"],
                            "sha256_uint16_little_endian": hashlib.sha256(
                                decoded.astype("<u2", copy=False).tobytes()).hexdigest(),
                        })
                    del decoded, expected
            file_size = output.stat().st_size
            for series in levels:
                for page in series.pages:
                    if not all(0 < offset < file_size and 0 < size <= file_size - offset
                               for offset, size in zip(page.dataoffsets, page.databytecounts)):
                        raise ValueError("Invalid compressed-tile byte range")
            return {
                "all_pixels_equal": True, "checked_pixels_all_levels": checked_pixels,
                "levels": [list(level.shape) for level in levels],
                "base_channel_pixel_sha256": channel_hashes,
                "calibration_and_channels_valid": True, "tile_byte_ranges_valid": True,
            }
    finally:
        del source


def sidecars(output: Path) -> dict:
    stem = output.name.removesuffix(".ome.tiff")
    xml_path = output.with_name(stem + ".ome.xml")
    offsets_path = output.with_name(output.name + ".offsets.json")
    alias = output.with_name(stem + ".offsets.json")
    with tifffile.TiffFile(output) as tif:
        xml_path.write_text(tif.ome_metadata + "\n")
        offsets = [int(page.offset) for page in tif.pages]
        # Viv indexes this array by channel-plane IFD index. Tile offsets are not
        # interchangeable with IFD offsets and would silently disable indexing.
        write_json(offsets_path, offsets)
        for index, offset in enumerate(json.loads(offsets_path.read_text())):
            if offset != tif.pages[index].offset:
                raise ValueError("IFD sidecar validation failed")
    if alias.is_symlink():
        if os.readlink(alias) != offsets_path.name:
            raise ValueError(f"Refusing to replace unexpected sidecar alias: {alias}")
    elif alias.exists():
        if alias.read_bytes() != offsets_path.read_bytes():
            raise ValueError(f"Refusing to replace unexpected sidecar: {alias}")
    else:
        alias.symlink_to(offsets_path.name)
    return {p.name: sha256(p) for p in (xml_path, offsets_path, alias)}


def process(source_path: Path, output: Path, report_path: Path, threads: int) -> dict:
    started = time.time()
    info = inspect_source(source_path)
    log(f"{source_path.name}: hashing source")
    source_hash = sha256(source_path)
    if report_path.exists() and output.exists():
        previous = json.loads(report_path.read_text())
        if previous.get("validation", {}).get("all_pixels_equal") and previous.get("output_sha256"):
            if previous["source_sha256"] != source_hash:
                raise ValueError(f"Source changed since checkpoint: {source_path}")
            if (previous.get("format_version") != FORMAT_VERSION
                    or previous.get("channel_mapping") != info["channel_mapping"]
                    or previous.get("pixel_size_um") != info["pixel_size_um"]):
                raise ValueError(f"Conversion contract changed since checkpoint: {output}")
            if previous["output_sha256"] != sha256(output):
                raise ValueError(f"Output changed since checkpoint: {output}")
            if previous.get("status") == "validated":
                for name, digest in previous["sidecar_sha256"].items():
                    if sha256(output.parent / name) != digest:
                        raise ValueError(f"Sidecar changed since checkpoint: {name}")
            else:
                previous["sidecar_sha256"] = sidecars(output)
                previous["status"] = "validated"
                previous.pop("error", None)
                previous["resumed_after_pixel_validation"] = True
                write_json(report_path, previous)
            log(f"{output.name}: validated checkpoint retained")
            return previous
    if output.exists():
        raise FileExistsError(f"Output exists without a validated checkpoint: {output}")
    partial = output.with_name(output.name + ".partial")
    if partial.exists():
        retained = partial.with_name(partial.name + f".interrupted-{time.time_ns()}")
        partial.rename(retained)
        log(f"Retained interrupted output: {retained.name}")
    report = dict(info, format_version=FORMAT_VERSION, status="converting",
                  converter_sha256=CONVERTER_SHA256,
                  source_sha256=source_hash, output=str(output.resolve()),
                  created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                  software={"python": sys.version, "numpy": np.__version__,
                            "tifffile": tifffile.__version__, "imagecodecs": imagecodecs.__version__},
                  conversion={"dtype": "uint16", "bigtiff": True, "compression": "LZW",
                              "byteorder": "little", "tile_yx": [TILE, TILE], "levels": LEVELS,
                              "downsampling": "repeated 2x2 integer-floor mean; omit unpaired edge",
                              "base_pixels": "unchanged; channel order follows explicit named mapping"})
    write_json(report_path, report)
    try:
        convert(source_path, partial, info, threads)
        report["status"] = "validating"
        write_json(report_path, report)
        report["validation"] = validate_pixels(source_path, partial, info, threads)
        # Recheck the source file after conversion to reject concurrent edits.
        if sha256(source_path) != source_hash:
            raise ValueError("Source file changed during conversion")
        report["output_sha256"] = sha256(partial)
        report["output_bytes"] = partial.stat().st_size
        report["status"] = "pixels_validated"
        write_json(report_path, report)
        partial.rename(output)
        report["sidecar_sha256"] = sidecars(output)
        report["status"] = "validated"
        report["elapsed_seconds"] = round(time.time() - started, 3)
        write_json(report_path, report)
        log(f"{output.name}: complete; {report['output_bytes']:,} bytes")
        return report
    except Exception as error:
        report["status"] = "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        write_json(report_path, report)
        raise


def write_identity_audit(reports_root: Path, donors: list[str]) -> dict:
    """Detect identical source files or pixels labelled as different donors."""
    source_groups: dict[str, list[str]] = {}
    pixel_groups: dict[str, list[str]] = {}
    validated = []
    for donor in donors:
        report_path = reports_root / (SELECTION[donor][1] + ".json")
        if not report_path.exists():
            continue
        report = json.loads(report_path.read_text())
        if report.get("status") != "validated":
            continue
        validated.append(donor)
        source_groups.setdefault(report["source_sha256"], []).append(donor)
        pixels = {entry["channel_name"]: entry["sha256_uint16_little_endian"]
                  for entry in report["validation"]["base_channel_pixel_sha256"]}
        signature = hashlib.sha256(json.dumps(pixels, sort_keys=True).encode()).hexdigest()
        pixel_groups.setdefault(signature, []).append(donor)
    duplicate_sources = [{"sha256": digest, "case_ids": group}
                         for digest, group in source_groups.items() if len(group) > 1]
    duplicate_pixels = [{"pixel_signature_sha256": digest, "case_ids": group}
                        for digest, group in pixel_groups.items() if len(group) > 1]
    has_duplicates = bool(duplicate_sources or duplicate_pixels)
    complete = len(validated) == len(donors)
    audit = {
        "requested_case_ids": donors, "validated_case_ids": validated,
        "status": ("duplicate_inputs_detected" if has_duplicates else
                   "complete_unique" if complete else "incomplete"),
        "duplicate_check_passed": complete and not has_duplicates,
        "duplicate_source_groups": duplicate_sources,
        "duplicate_pixel_groups": duplicate_pixels,
        "limitation": "Unique pixels do not establish correct donor attribution or scientific acceptance.",
    }
    write_json(reports_root / "cohort_identity_audit.json", audit)
    return audit


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="islet-image-conversion-") as temp:
        root = Path(temp)
        source = root / "source.tif"
        # Distinct channels and nonzero edge pixels expose channel-order, endian,
        # padding and odd-dimension downsampling mistakes.
        rng = np.random.default_rng(51011)
        data = rng.integers(0, 65536, size=(33, 65, 97), dtype=np.uint16)
        labels = list(reversed(CHANNEL_ORDER[:33]))
        tifffile.imwrite(source, data, imagej=True, resolution=(2649.3832, 2649.3832),
                         metadata={"axes": "CYX", "unit": "mm", "Labels": labels})
        output = root / "test.ome.tiff"
        report = process(source, output, root / "report.json", 2)
        if not report["validation"]["all_pixels_equal"]:
            raise AssertionError("Synthetic pixel validation failed")
        retained = process(source, output, root / "report.json", 2)
        if retained != report:
            raise AssertionError("Resume modified a validated checkpoint")
        if canonical_channel("KRT8") != "KRT8":
            raise AssertionError("Ambiguous KRT8 channel must remain distinct")
        interrupted = dict(report, status="pixels_validated")
        interrupted.pop("sidecar_sha256")
        write_json(root / "report.json", interrupted)
        recovered = process(source, output, root / "report.json", 2)
        if not recovered.get("resumed_after_pixel_validation"):
            raise AssertionError("Sidecar-finalization interruption was not recovered")
        offsets = root / "test.ome.tiff.offsets.json"
        offsets.write_text("[]\n")
        try:
            process(source, output, root / "report.json", 2)
        except ValueError as error:
            if "Sidecar changed" not in str(error):
                raise
        else:
            raise AssertionError("Corrupted checkpoint sidecar was accepted")
        for donor in ("6090", "6539"):
            write_json(root / (SELECTION[donor][1] + ".json"), report)
        identity = write_identity_audit(root, ["6090", "6539"])
        if identity["duplicate_check_passed"] or not identity["duplicate_source_groups"]:
            raise AssertionError("Identical images with different donor labels were accepted")
        log("Synthetic conversion, channel reorder, pyramid, sidecar, resume and duplicate tests passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    parser.add_argument("--donors", nargs="+", choices=list(SELECTION), default=list(SELECTION))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--inspect-only", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.threads < 1:
        parser.error("--threads must be positive")
    if args.inspect_only:
        for donor in args.donors:
            print(json.dumps(dict(donor=donor, **inspect_source(args.source_root / SELECTION[donor][0]))))
        return 0
    args.output_root.mkdir(parents=True, exist_ok=True)
    reports_root = args.output_root.parent / "reports/image_conversion"
    reports_root.mkdir(parents=True, exist_ok=True)
    for donor in args.donors:
        filename, stem = SELECTION[donor]
        process(args.source_root / filename, args.output_root / (stem + ".ome.tiff"),
                reports_root / (stem + ".json"), args.threads)
        write_identity_audit(reports_root, args.donors)
    audit = write_identity_audit(reports_root, args.donors)
    if not audit["duplicate_check_passed"]:
        log("Conversion finished, but duplicate donor-labelled inputs prevent cohort acceptance; see cohort_identity_audit.json")
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
