#!/usr/bin/env python3
"""Record the classifier-defined islet eligibility contract.

An islet candidate is created by the fixed endocrine pixel classifier.  Its
eligibility therefore depends only on the classifier geometry being clear of
the image edge and containing at least one canonical core cell.  Cell typing
is downstream evidence and must never revoke the islet's classifier identity.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import shape

from common import sha256, write_json
from reprocess_candidate import class_name, feature_list


CLASSIFIER_DISPLAY_CLASSES = {"Endocrine pixel candidate", "Islet"}


def materialize(root: Path) -> dict:
    root = Path(root).resolve(strict=True)
    cohort = pd.read_csv(root / "cohort.csv", dtype=str)
    islet_path = root / "analysis/islets.parquet"
    membership_path = root / "analysis/memberships.parquet"
    islets = pd.read_parquet(islet_path)
    memberships = pd.read_parquet(membership_path)
    islets["case_id"] = islets.case_id.astype(str)
    memberships["case_id"] = memberships.case_id.astype(str)
    if islets.islet_id.duplicated().any():
        raise ValueError("Classifier eligibility requires unique islet IDs")
    if set(islets.case_id) != set(cohort.case_id):
        raise ValueError("Classifier eligibility does not match the manifest cohort")

    edge_rows = []
    source_islets_sha256 = sha256(islet_path)
    source_memberships_sha256 = sha256(membership_path)
    sources = {str(islet_path): source_islets_sha256, str(membership_path): source_memberships_sha256}
    core20 = memberships.loc[
        memberships.compartment.eq("core")
        & (memberships.band_um.eq(20) if "band_um" in memberships else True)
    ].copy()
    if core20.duplicated(["case_id", "object_id", "islet_id"]).any():
        raise ValueError("Duplicate classifier-core cell membership")

    native_core_parts = []
    membership_differences = []
    measured_core_counts = {}
    for row in cohort.itertuples(index=False):
        cid = str(row.case_id)
        source_path = Path(row.boundary_source).resolve(strict=True)
        native_root = Path(row.native_core_root).resolve(strict=True)
        native_path = native_root / "cores.geojson.gz"
        native_manifest_path = native_root / "manifest.json"
        native_membership_path = native_root / "native_core_memberships.csv"
        adaptive_manifest_path = Path(row.boundary_candidates).resolve(strict=True).parent / "prediction_manifest.json"
        source = json.loads(source_path.read_text())
        native_manifest = json.loads(native_manifest_path.read_text())
        adaptive_manifest = json.loads(adaptive_manifest_path.read_text())
        if sha256(native_path) != native_manifest["cores_sha256"]:
            raise ValueError(f"{cid}: native classifier cores changed")
        if sha256(native_membership_path) != native_manifest["memberships_sha256"]:
            raise ValueError(f"{cid}: native classifier membership changed")
        edge_margin = float(adaptive_manifest.get("boundary_config", {}).get("edge_margin_um", np.nan))
        if not np.isfinite(edge_margin) or edge_margin <= 0:
            raise ValueError(f"{cid}: classifier boundary lacks an image-edge exclusion")
        height, width = map(float, source["image_shape_cyx"][-2:])
        pixel_size = float(source["pixel_size_um"])
        if min(height, width, pixel_size) <= 0:
            raise ValueError(f"{cid}: invalid source dimensions or calibration")

        features = feature_list(native_path)
        names = [str(feature.get("properties", {}).get("name", "")) for feature in features]
        if len(names) != len(set(names)) or any(not name for name in names):
            raise ValueError(f"{cid}: native classifier cores need unique names")
        expected = set(islets.loc[islets.case_id.eq(cid), "islet_key"].astype(str))
        if set(names) != expected:
            raise ValueError(f"{cid}: classifier geometry differs from analysis islets")
        if any(class_name(feature) not in CLASSIFIER_DISPLAY_CLASSES for feature in features):
            raise ValueError(f"{cid}: an unexpected display class is mixed into classifier cores")
        if any("INS/GCG/SST foreground" not in str(feature.get("properties", {}).get("metadata", {}).get("Segmentation method", ""))
               for feature in features):
            raise ValueError(f"{cid}: classifier marker provenance is missing")
        model_sha = str(adaptive_manifest.get("model_sha256", ""))
        if len(model_sha) != 64 or any(str(feature.get("properties", {}).get("metadata", {}).get("Model SHA256", "")) != model_sha
                                       for feature in features):
            raise ValueError(f"{cid}: promoted cores do not match the recorded endocrine classifier model")

        native = pd.read_csv(native_membership_path, dtype=str)
        native["case_id"] = native.case_id.astype(str)
        native["islet_id"] = native.case_id + "_" + native.islet_key
        native_core_parts.append(native[["case_id", "object_id", "islet_id"]].assign(compartment="core", band_um=20))
        current = core20.loc[core20.case_id.eq(cid), ["case_id", "object_id", "islet_id"]]
        native_keys = set(map(tuple, native[["case_id", "object_id", "islet_id"]].to_numpy()))
        current_keys = set(map(tuple, current[["case_id", "object_id", "islet_id"]].astype(str).to_numpy()))
        membership_differences.append({"case_id": cid, "added_from_native": len(native_keys-current_keys),
                                       "removed_against_native": len(current_keys-native_keys)})
        corrected = pd.read_parquet(root / "corrected" / f"{cid}-cells.parquet",
                                    columns=["object_id", "marker_measurements_complete"])
        corrected["object_id"] = corrected.object_id.astype(str)
        joined_native = native.merge(corrected, on="object_id", how="left", validate="many_to_one")
        if joined_native.marker_measurements_complete.isna().any():
            raise ValueError(f"{cid}: native classifier membership references a missing canonical cell")
        measured_core_counts.update(joined_native.groupby("islet_id").marker_measurements_complete.sum().astype(int).to_dict())

        for feature in features:
            geometry = shape(feature["geometry"])
            min_x, min_y, max_x, max_y = geometry.bounds
            clearance_um = min(min_x, min_y, width - max_x, height - max_y) * pixel_size
            edge_rows.append({
                "case_id": cid,
                "islet_key": str(feature["properties"]["name"]),
                "classifier_defined_core": True,
                "classifier_class": class_name(feature),
                "classifier_model_sha256": model_sha,
                "classifier_edge_margin_um": edge_margin,
                "image_edge_clearance_um": float(clearance_um),
                "image_edge_contact": bool(clearance_um <= 0),
            })
        for path in (source_path, native_path, native_manifest_path, native_membership_path, adaptive_manifest_path):
            sources[str(path)] = sha256(path)

    edge = pd.DataFrame(edge_rows)
    if edge.duplicated(["case_id", "islet_key"]).any() or len(edge) != len(islets):
        raise ValueError("Classifier eligibility audit did not produce one row per islet")
    native_core = pd.concat(native_core_parts, ignore_index=True)
    rebuilt_memberships = pd.concat([memberships.loc[~memberships.compartment.eq("core")], native_core], ignore_index=True)
    if rebuilt_memberships.duplicated(["case_id", "object_id", "islet_id", "compartment", "band_um"]).any():
        raise ValueError("Native classifier membership introduced duplicate compartment rows")
    membership_temporary = membership_path.with_suffix(".parquet.tmp")
    rebuilt_memberships.to_parquet(membership_temporary, index=False, compression="zstd")
    membership_temporary.replace(membership_path)

    replace = [column for column in edge if column not in {"case_id", "islet_key"}]
    enriched = islets.drop(columns=replace, errors="ignore").merge(
        edge, on=["case_id", "islet_key"], how="left", validate="one_to_one"
    )
    if enriched.classifier_defined_core.isna().any():
        raise ValueError("Classifier provenance is missing for an analysis islet")
    counts = native_core.groupby("islet_id").size()
    enriched["classifier_core_cells"] = enriched.islet_id.map(counts).fillna(0).astype(int)
    enriched["total_cells_core"] = enriched.classifier_core_cells
    enriched["measured_cells_core"] = enriched.islet_id.map(measured_core_counts).fillna(0).astype(int)
    enriched["classifier_eligible"] = (~enriched.image_edge_contact) & enriched.classifier_core_cells.gt(0)
    enriched["classifier_exclusion_reason"] = np.select(
        [enriched.image_edge_contact, enriched.classifier_core_cells.eq(0)],
        ["image_edge_contact", "no_core_cells_in_classifier_area"],
        default="",
    )

    temporary = islet_path.with_suffix(".parquet.tmp")
    enriched.to_parquet(temporary, index=False, compression="zstd")
    temporary.replace(islet_path)
    by_donor = {
        str(case_id): int(number)
        for case_id, number in enriched.loc[enriched.classifier_eligible].groupby("case_id").size().items()
    }
    report = {
        "schema_version": 1,
        "definition": {
            "all_of": [
                "the endocrine pixel-classifier core has no image-edge contact",
                "at least one canonical cell centroid is inside the classifier-defined core",
            ],
            "downstream_typing_can_revoke_eligibility": False,
            "endocrine_denominator": "all canonical core cells inside the classifier-defined area, counted once",
        },
        "islets": int(len(enriched)),
        "eligible_islets": int(enriched.classifier_eligible.sum()),
        "image_edge_islets": int(enriched.image_edge_contact.sum()),
        "empty_classifier_cores": int(enriched.classifier_core_cells.eq(0).sum()),
        "minimum_image_edge_clearance_um": float(enriched.image_edge_clearance_um.min()),
        "eligible_islets_by_donor": by_donor,
        "native_membership_reconciliation": membership_differences,
        "sources": sources,
        "source_islets_sha256": source_islets_sha256,
        "source_memberships_sha256": source_memberships_sha256,
        "memberships_sha256": sha256(membership_path),
        "islets_sha256": sha256(islet_path),
    }
    write_json(root / "analysis/classifier_eligibility_manifest.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    print(json.dumps(materialize(parser.parse_args().input_root), indent=2))
