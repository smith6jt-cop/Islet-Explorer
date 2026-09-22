#!/usr/bin/env python3
"""Inventory the six corrected image inputs without promoting a release.

The original conversion directory contains a rejected 6539 duplicate. Select
the recovered acquisition explicitly, verify current file hashes, and retain
the initial duplicate audit as historical evidence.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

from convert_images import REPO, SELECTION, process, sha256, write_json


ROOT = REPO / "data/development/expanded-20260911"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(path.read_text())


def record(path):
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256(path)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--with-preparation", action="store_true",
                        help="Also join the two completed preparation audits; fail if either is absent")
    args = parser.parse_args()
    root = args.root.resolve()
    recovered = root / "recovery"
    match_path = root / "reports/image_conversion/source_recovery_6090.json"
    reconstruction_path = recovered / "source_reconstruction_6539.json"
    project_path = recovered / "recovered_project_provenance.json"
    matched, rebuilt, project = map(read_json, [match_path, reconstruction_path, project_path])
    require(matched["all_planes_identical"] is True and len(matched["planes"]) == 36,
            "6090 source was not reconciled against all 36 original planes")
    require(rebuilt["status"] == "validated" and len(rebuilt["planes"]) == 36,
            "6539 source reconstruction is incomplete")
    require(project["status"] == "source_reconciled"
            and project["source_match_report_sha256"] == sha256(match_path)
            and project["reconstruction_report_sha256"] == sha256(reconstruction_path),
            "Recovered processing project has different source reconciliation evidence")
    rows = []
    for case, (filename, stem) in SELECTION.items():
        owner = recovered if case == "6539" else root
        report_path = owner / "reports/image_conversion" / f"{stem}.json"
        report = read_json(report_path)
        source, output = Path(report["source"]), Path(report["output"])
        require(source.name == filename and output.name == f"{stem}.ome.tiff"
                and report["status"] == "validated", f"Incorrect conversion identity for {case}")
        require(output.is_relative_to(owner / "images"), f"Unexpected viewer output for {case}")
        if case in {"6090", "6539"}:
            require(str(source) == project["cases"][case]["source"]
                    and report["source_sha256"] == project["cases"][case]["source_sha256"],
                    f"Source reconciliation does not match the {case} viewer")
        # The converter verifies current source, TIFF and all sidecar hashes
        # against its completed pixel-validation checkpoint before returning.
        process(source, output, report_path, threads=2)
        pixels = {entry["channel_name"]: entry["sha256_uint16_little_endian"]
                  for entry in report["validation"]["base_channel_pixel_sha256"]}
        require(len(pixels) == report["shape_cyx"][0], f"Incomplete pixel signature for {case}")
        signature = hashlib.sha256(json.dumps(pixels, sort_keys=True).encode()).hexdigest()
        rows.append({"case_id": case, "donor_source": "HPAP" if case == "78" else "nPOD",
                     "source": str(source), "source_sha256": report["source_sha256"],
                     "viewer_tiff": str(output), "viewer_tiff_bytes": output.stat().st_size,
                     "viewer_tiff_sha256": report["output_sha256"],
                     "sidecars": [record(output.parent / name) for name in report["sidecar_sha256"]],
                     "shape_cyx": report["shape_cyx"], "pixel_size_um": report["pixel_size_um"],
                     "pixel_signature_sha256": signature,
                     "validated_pixels_all_levels": report["validation"]["checked_pixels_all_levels"],
                     "conversion_report": record(report_path)})
    for key in ["case_id", "source_sha256", "pixel_signature_sha256"]:
        require(len({row[key] for row in rows}) == 6, f"Corrected inventory still contains duplicate {key}")
    summary = {"recorded_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "status": "six_sources_reconciled_and_viewer_hashes_verified", "release_ready": False,
               "script_sha256": sha256(Path(__file__)), "cases": rows,
               "reconciliation_evidence": list(map(record, [match_path, reconstruction_path, project_path])),
               "supersedes_for_source_selection": str(root / "reports/image_conversion/cohort_identity_audit.json"),
               "excluded_duplicate_viewer": str(root / "images/ND_6539.ome.tiff"),
               "totals": {"new_images": 6, "new_npod_images": 5, "new_hpap_images": 1,
                          "target_total_images": 21, "target_npod_images": 20, "target_hpap_images": 1,
                          "viewer_tiff_bytes": sum(row["viewer_tiff_bytes"] for row in rows),
                          "viewer_sidecar_logical_bytes": sum(p["bytes"] for row in rows for p in row["sidecars"]),
                          "validated_pixels_all_levels": sum(row["validated_pixels_all_levels"] for row in rows)},
               "limits": ["Source attribution is supported by the local raw acquisition folders and channel records.",
                          "Original 6539 annotations were drawn on duplicated 6090 pixels; their geometry was not transferred.",
                          "No phenotype assignment, joint refit, biological acceptance or runtime release is implied.",
                          "Stage viewer inputs from the six paths in this inventory, never the entire development images directory."]}
    if args.with_preparation:
        paths = [base / "reports/preparation/summary.json" for base in [root, recovered]]
        preparations = list(map(read_json, paths))
        cases = [case for prep in preparations for case in prep["cases"]]
        require(len(cases) == 6 and {case["case_id"] for case in cases} == set(SELECTION),
                "Preparation audits do not cover the six corrected cases exactly once")
        by_id = {row["case_id"]: row for row in rows}
        for case in cases:
            expected = by_id[case["case_id"]]
            require(case["donor_source"] == expected["donor_source"]
                    and case["viewer_tiff_bytes"] == expected["viewer_tiff_bytes"], "Preparation identity mismatch")
            conversions = [item for item in case["evidence"] if item["path"].startswith("reports/image_conversion/")]
            require(len(conversions) == 1 and conversions[0]["sha256"] == expected["conversion_report"]["sha256"],
                    "Preparation and inventory refer to different image conversions")
        summary["preparation"] = {
            "audits": list(map(record, paths)),
            "cases": [{k: v for k, v in case.items() if k not in {"evidence", "raw_h5ad"}} for case in cases],
            "totals": {key: sum(case[key] for case in cases) for key in
                       ["cell_objects", "fully_measured_cells", "unmeasured_cell_objects", "islets",
                        "composition_rows", "marker_rows", "target_rows"]}}
    destination = root / "reports/image_conversion/reconciled_inventory.json"
    write_json(destination, summary)
    print(destination, flush=True)
    print(json.dumps(summary["totals"]), flush=True)
    if args.with_preparation:
        print(json.dumps(summary["preparation"]["totals"]), flush=True)


if __name__ == "__main__":
    main()
