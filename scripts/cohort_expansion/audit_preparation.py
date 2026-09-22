#!/usr/bin/env python3
"""Verify completed offline preparation and write a candidate evidence inventory.

This checks the selected acquisitions without clearing identity-held inputs,
choosing a missing-measurement policy, or asserting scientific/release acceptance.
"""
import argparse
import datetime
import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape

from prepare_analysis import BASELINE, COHORT, DEFAULT_ROOT, QUANTIFICATION_PIXEL_UM, read_cohort, sha256


STAGES = ["normalize", "tissue", "sam_prepare", "sam", "sam_import", "cells", "markers", "targets", "export", "qc"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(path.read_text())


def file_record(path, root):
    return {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}


def validate_regions(annotations, composition, markers, targets, expected_islets):
    islets = annotations.loc[annotations["Classification"] == "Islet"].copy()
    require(len(islets) == expected_islets and not islets["Name"].isna().any()
            and not islets["Name"].duplicated().any(), "Invalid islet annotation identities/count")
    names = set(islets["Name"])
    require(not composition["name"].duplicated().any() and set(composition["name"]) == names,
            "Composition does not join one-to-one to islet annotations")
    annotation_counts = islets.set_index("Name")["Islet: cells total"].sort_index()
    composition_counts = composition.set_index("name")["cells_total"].sort_index()
    require(np.array_equal(annotation_counts.to_numpy(), composition_counts.to_numpy()),
            "Composition cell counts differ from the islet annotations")
    for table, key, type_key in [(markers, "marker", "region_type"), (targets, "class", "type")]:
        require(not table.duplicated(["region", key]).any(), "Duplicate regional measurement key")
        cores = table.loc[table[type_key] == "islet_core", "region"]
        require(set(cores) == {f"{name}_core" for name in names}, "Regional core coverage differs from islets")
        regional = table.loc[table[type_key].str.startswith("islet_"), "region"]
        base = regional.str.extract(r"^(Islet_[1-9][0-9]*)_(?:core|band|union)$", expand=False)
        require(not base.isna().any() and set(base).issubset(names), "Regional measurement has an unknown islet")
    return {"islets": len(islets), "composition_rows": len(composition),
            "marker_rows": len(markers), "target_rows": len(targets),
            "annotation_classes": {str(k): int(v) for k, v in annotations["Classification"].value_counts().items()}}


def validate_cell_roundtrip(data, source, marker_names, row, conversion):
    require(data.var_names.tolist() == marker_names, "Raw import marker order differs from baseline")
    require("phenotype" not in data.obs and data.uns["analysis_stage"] == "raw_measurements_unphenotyped",
            "Unexpected phenotyped/non-raw data in preparation inventory")
    require(data.obs_names.is_unique and set(data.obs["donor_id"].astype(str)) == {row["case_id"]}
            and set(data.obs["donor_source"].astype(str)) == {row["donor_source"]}, "Raw import donor identity mismatch")
    marker_columns = [f"Cell: {marker}: Mean" for marker in marker_names]
    columns = ["Object ID", "Centroid X µm", "Centroid Y µm"] + marker_columns
    start = 0
    with pd.read_csv(source, sep="\t", usecols=columns, chunksize=25000, float_precision="round_trip") as chunks:
        for frame in chunks:
            stop = start + len(frame)
            require(stop <= data.n_obs, "Raw import is truncated")
            expected = frame[marker_columns].to_numpy(dtype=np.float32)
            raw = np.asarray(data.X[start:stop])
            require(np.array_equal(raw, expected, equal_nan=True), "Raw H5AD values differ from the TSV")
            require(np.allclose(data.layers["log"][start:stop], np.log1p(expected), rtol=1e-6, atol=0, equal_nan=True),
                    "Log layer differs from log1p of the TSV")
            complete = np.isfinite(raw).all(axis=1)
            require(np.array_equal(complete, data.obs["marker_measurements_complete"].iloc[start:stop].to_numpy()),
                    "Missing-measurement flags differ from the raw matrix")
            require(np.isnan(raw[~complete]).all() and (raw[complete] >= 0).all(), "Partial/invalid marker measurements")
            require(data.obs["Object ID"].iloc[start:stop].tolist() == frame["Object ID"].tolist(),
                    "Raw import object identities/order differ from the TSV")
            for source_key, output_key, dimension in [("Centroid X µm", "X_centroid", 2), ("Centroid Y µm", "Y_centroid", 1)]:
                values = data.obs[output_key].iloc[start:stop].to_numpy()
                require(np.array_equal(values, frame[source_key].to_numpy()), "Raw import changed centroid coordinates")
                require(np.isfinite(values).all() and (values >= 0).all()
                        and (values <= conversion["shape_cyx"][dimension] * QUANTIFICATION_PIXEL_UM).all(),
                        "Cell centroid falls outside its source image")
            start = stop
    require(start == data.n_obs, "Raw import has extra rows")


def polygon_boundary_contains(vertices, x, y):
    """Match QuPath 0.7 PolygonROI/WindingTest at an exact boundary point.

    PolygonROI uses a nonzero winding number rather than JTS strict containment.
    Keep the original coordinates: buffering or nudging would change membership.
    """
    winding = 0
    for (x0, y0), (x1, y1) in zip(vertices, vertices[1:]):
        side = (x1 - x0) * (y - y0) - (x - x0) * (y1 - y0)
        if y0 <= y < y1 and side > 0:
            winding += 1
        elif y1 <= y < y0 and side < 0:
            winding -= 1
    return winding != 0


def validate_spatial_counts(data, geometry_path, composition):
    """Independently reproduce QuPath's core counts from unrounded centroids."""
    with gzip.open(geometry_path, "rt") as stream:
        features = json.load(stream)["features"]
    cores = [feature for feature in features
             if feature["properties"].get("classification", {}).get("name") == "Islet"]
    expected = composition.set_index("name")["cells_total"]
    require(len(cores) == len(expected)
            and {feature["properties"]["name"] for feature in cores} == set(expected.index),
            "GeoJSON islet identities differ from composition")
    x = data.obs["X_centroid"].to_numpy() / QUANTIFICATION_PIXEL_UM
    y = data.obs["Y_centroid"].to_numpy() / QUANTIFICATION_PIXEL_UM
    require(np.isfinite(x).all() and np.isfinite(y).all(), "Invalid spatial-audit coordinates")
    index = shapely.STRtree(shapely.points(x, y))
    membership_counts = np.zeros(data.n_obs, dtype=np.int32)
    boundary_checks = []
    for feature in cores:
        name = feature["properties"]["name"]
        geometry = shape(feature["geometry"])
        require(geometry.is_valid and not geometry.is_empty, f"Invalid islet geometry: {name}")
        require(geometry.geom_type == "Polygon" and not geometry.interiors,
                f"Spatial audit requires the workflow's simple PolygonROI cores: {name}")
        members = index.query(geometry, predicate="contains")
        boundary_members = []
        for i in index.query(geometry.boundary, predicate="intersects"):
            included = polygon_boundary_contains(list(geometry.exterior.coords), x[i], y[i])
            boundary_checks.append({"islet": name, "object_id": str(data.obs_names[i]),
                                    "x_px": float(x[i]), "y_px": float(y[i]), "included": included})
            if included:
                boundary_members.append(i)
        members = np.concatenate([members, np.asarray(boundary_members, dtype=np.int64)])
        require(len(members) == int(expected[name]),
                f"Spatial core count differs for {name}: {len(members)} versus {int(expected[name])}")
        membership_counts[members] += 1
    return {"all_core_counts_equal": True, "core_regions_checked": len(cores),
            "core_cell_membership_rows": int(membership_counts.sum()),
            "unique_cells_in_cores": int(np.count_nonzero(membership_counts)),
            "cells_in_multiple_cores": int(np.count_nonzero(membership_counts > 1)),
            "predicate": "PolygonROI.contains: strict interior plus QuPath WindingTest for exact boundary ties",
            "exact_boundary_checks": boundary_checks, "shapely_version": shapely.__version__,
            "limitation": "Overlapping core memberships are counted explicitly; no single-islet assignment is invented."}


def audit_case(root, row, marker_names):
    cid = row["case_id"]
    folder = root / "exports" / cid.zfill(4)
    log_folder = root / "logs" / cid.zfill(4)
    records = []
    for stage in STAGES:
        checkpoint = log_folder / f"{stage}.done.json"
        done = read_json(checkpoint)
        scripts = [Path(value) for value in done["command"] if value.endswith((".groovy", ".py"))]
        require(done["case"] == cid.zfill(4) and done["stage"] == stage and len(scripts) == 1
                and sha256(scripts[0]) == done["script_sha256"], f"Case {cid}: stale {stage} checkpoint")
        records.append(file_record(checkpoint, root))
    imported_path = root / "analysis/raw_cells" / f"{cid}.json"
    imported = read_json(imported_path)
    raw_path = imported_path.with_suffix(".h5ad")
    require(sha256(raw_path) == imported["output_sha256"], f"Case {cid}: raw import hash mismatch")
    manifest = read_json(folder / "export_manifest.json")
    require(manifest.get("numeric_encoding") == "java_double_to_string", "Rounded measurement export is not current")
    segmentation = read_json(folder / "segmentation_provenance.json")
    scope = read_json(folder / "structural_scope_provenance.json")
    for document in [manifest, segmentation]:
        require(document["source_filename"] == row["source_filename"], "Export/segmentation source mismatch")
    require(segmentation["model_version"] == "0.1.1", "Unexpected InstanSeg model")
    require(scope["scope"] == "Tissue annotations" and scope["cell_uuids_preserved"] is True
            and scope["cell_parents_preserved"] is True, "Structural measurement scope/identity check failed")
    counts = [manifest["cell_count"], segmentation["cell_count"], imported["cells"],
              scope["cell_count_before"], scope["cell_count_after"]]
    require(len(set(counts)) == 1, "Cell count changed across preparation stages")
    for name, key in [("cells.tsv", "cells_sha256"), ("annotations.tsv", "annotations_sha256")]:
        record = file_record(folder / name, root)
        require(record["sha256"] == manifest[key], f"{name} differs from the completed manifest")
        if name == "cells.tsv":
            require(record["sha256"] == imported["inputs"][key], "Raw import consumed a different export")
        records.append(record)
    require(sha256(folder / "export_manifest.json") == imported["inputs"]["export_manifest_sha256"],
            "Export manifest changed after raw import")
    conversions = []
    for path in (root / "reports/image_conversion").glob("*.json"):
        content = read_json(path)
        if isinstance(content, dict) and Path(content.get("source", "")).name == row["source_filename"]:
            conversions.append((path, content))
    require(len(conversions) == 1, "Missing/ambiguous viewer conversion record")
    conversion_path, conversion = conversions[0]
    require(conversion["status"] == "validated"
            and sha256(conversion_path) == imported["inputs"]["conversion_report_sha256"], "Viewer conversion record changed")
    annotations = pd.read_csv(folder / "annotations.tsv", sep="\t", dtype={"Name": str}, low_memory=False)
    require(len(annotations) == manifest["annotation_count"], "Annotation count mismatch")
    tables = [pd.read_csv(folder / f"islet_{kind}_{cid.zfill(4)}.tsv", sep="\t")
              for kind in ["composition", "markers", "targets"]]
    regions = validate_regions(annotations, *tables, segmentation["islet_count"])
    data = ad.read_h5ad(raw_path, backed="r")
    try:
        validate_cell_roundtrip(data, folder / "cells.tsv", marker_names, row, conversion)
        spatial = validate_spatial_counts(data, folder / "annotations.geojson.gz", tables[0])
        require(data.n_obs == counts[0], "Raw H5AD count mismatch")
        unmeasured = data.obs.loc[~data.obs["marker_measurements_complete"]]
        require(len(unmeasured) == imported["unmeasured_cell_objects"], "Unmeasured object count mismatch")
        review_columns = ["Object ID", "X_centroid", "Y_centroid", "Cell Area", "Nucleus Area"]
        review = unmeasured[review_columns].copy()
        review.insert(0, "case_id", cid)
        review["reason"] = "all_31_marker_means_missing_preserved_for_review"
        total, unmeasured_count = data.n_obs, len(unmeasured)
    finally:
        data.file.close()
    for name in ["export_manifest.json", "segmentation_provenance.json", "structural_scope_provenance.json",
                 "annotations.geojson.gz", "qc/crop.json", "qc/objects.geojson", "qc/DAPI_INS_GCG.png"]:
        records.append(file_record(folder / name, root))
    for kind in ["composition", "markers", "targets"]:
        records.append(file_record(folder / f"islet_{kind}_{cid.zfill(4)}.tsv", root))
    records += [file_record(imported_path, root), file_record(conversion_path, root),
                file_record(root / "reports/segmentation" / f"{cid.zfill(4)}_cell_review.png", root)]
    return {"case_id": cid, "donor_source": row["donor_source"], "source_filename": row["source_filename"],
            "cell_objects": total, "fully_measured_cells": total - unmeasured_count,
            "unmeasured_cell_objects": unmeasured_count, **regions,
            "spatial_count_validation": spatial,
            "quantification_pixel_size_um": QUANTIFICATION_PIXEL_UM,
            "viewer_pixel_size_um": conversion["pixel_size_um"],
            "viewer_tiff_bytes": Path(conversion["output"]).stat().st_size,
            "raw_h5ad": {"path": str(raw_path.relative_to(root)), "bytes": raw_path.stat().st_size,
                         "sha256": imported["output_sha256"]}, "evidence": records}, review


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--cohort", type=Path, default=COHORT)
    parser.add_argument("--baseline", type=Path, default=BASELINE)
    args = parser.parse_args()
    root = args.root.resolve()
    cohort = read_cohort(args.cohort)
    baseline = ad.read_h5ad(args.baseline, backed="r")
    marker_names = baseline.var_names.tolist()
    baseline.file.close()
    held = [row["case_id"] for row in cohort if row["inclusion_status"] == "identity_hold"]
    for cid in held:
        require(not (root / "analysis/raw_cells" / f"{cid}.h5ad").exists(), "Identity-held donor entered raw import")
    results, reviews = [], []
    for row in cohort:
        if row["inclusion_status"] != "selected":
            continue
        result, review = audit_case(root, row, marker_names)
        results.append(result)
        reviews.append(review)
        print(f"Validated {row['case_id']}: {result['cell_objects']:,} objects, {result['islets']:,} islets", flush=True)
    destination = root / "reports/preparation"
    destination.mkdir(parents=True, exist_ok=True)
    review_path = destination / "unmeasured_cell_objects.csv"
    pd.concat(reviews, ignore_index=True).to_csv(review_path, index=False)
    summary = {"recorded_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "status": "offline_preparation_validated_with_scientific_holds", "release_ready": False,
               "cohort_sha256": sha256(args.cohort), "audit_script_sha256": sha256(__file__),
               "selected_cases": [result["case_id"] for result in results], "identity_held_cases": held,
               "markers": marker_names, "cases": results,
               "totals": {key: sum(result[key] for result in results)
                          for key in ["cell_objects", "fully_measured_cells", "unmeasured_cell_objects", "islets",
                                      "composition_rows", "marker_rows", "target_rows", "viewer_tiff_bytes"]},
               "candidate_logical_file_bytes": {name: sum(p.stat().st_size for p in (root / name).rglob("*") if p.is_file())
                                                for name in ["images", "exports", "analysis", "qupath", "logs", "metadata"]},
               "unmeasured_inventory": file_record(review_path, root),
               "remaining": (["join_separate_recovered_case_audit_for_final_cohort"] if held else []) + ["review_segmentation_and_unmeasured_objects",
                             "recover_original_phenotype_method_or_review_replacement", "joint_refit_and_derived_app_data",
                             "candidate_release_evidence_and_capacity", "UFRC_assigned_host_and_acceptance"],
               "limits": ["Raw archive retains all objects and NaNs; no missing value imputation or cell exclusion was applied.",
                          "Regional joins use islet names; QuPath Parent remains the original export field.",
                          "Viewer pixel-validation reports are linked; this audit does not rerun whole-image pixel validation.",
                          "Staging bytes include invalidated development exports and are not a PubApps runtime payload estimate.",
                          "This inventory does not establish biological, disclosure, or deployment acceptance."]}
    (destination / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary["totals"]), flush=True)


if __name__ == "__main__":
    main()
