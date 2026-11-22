"""Compare v25 tissue/endocrine residuals with the native outlines and original classifier mask."""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import Point, shape

sys.path.insert(0, str(Path(__file__).resolve().parent / "typing"))
from endocrine_object_audit import MARKERS, original_mask_support


def audit(args):
    project = json.loads((args.project / "project.qpproj").read_text())
    replay = json.loads(args.replay.read_text())
    records = {item["image"]: item for item in replay["images"]}
    frames = []
    for entry in project["images"]:
        tag = entry["imageName"]
        source = args.fit / f"{tag}.mask.npy"
        (digest,) = [
            v
            for k, v in records[tag]["source_sha256"].items()
            if k.endswith("/" + source.name)
        ]
        with source.open("rb") as stream:
            if hashlib.file_digest(stream, "sha256").hexdigest() != digest:
                raise ValueError(f"{tag}: mask differs from frozen classifier replay")
        mask = np.load(source, mmap_mode="r")
        metadata = json.loads(
            (args.project / "data" / str(entry["entryID"]) / "server.json").read_text()
        )["metadata"]
        cal = metadata["pixelCalibration"]
        sx, sy = (cal[key]["value"] for key in ("pixelWidth", "pixelHeight"))
        assert np.isclose(sx, sy) and mask.shape == (
            (metadata["height"] + 1) // 2,
            (metadata["width"] + 1) // 2,
        )
        with (args.geometry / f"{tag}.jsonl").open() as stream:
            native = [json.loads(line) for line in stream]
        ann = [r for r in native if r["kind"] == "annotation"]
        polys = [shape(r["geometry"]) for r in ann]
        tree = shapely.STRtree(polys)
        with (args.fit.parent / "script1_outlines" / f"{tag}.geojson").open() as stream:
            source_outlines = json.load(stream)["features"]
        original = [
            shape(f["geometry"])
            for f in source_outlines
            if f.get("properties", {}).get("classification", {}).get("name") == "Islet"
        ]
        original_tree = shapely.STRtree(original)
        cells = [r for r in native if r["kind"] == "cell"]
        cols = ["object_id", "start_class", "n_positive", "label", "status"] + [
            p + m for p in ("pos__", "canary__", "fold__") for m in MARKERS
        ]
        identity = pd.read_parquet(
            args.identity / f"{tag}.parquet",
            columns=cols,
            filters=[
                ("label", "==", "Tissue"),
                ("status", "==", "endocrine_marker_only"),
            ],
        ).set_index("object_id", verify_integrity=True)
        if set(identity.index) != {r["id"] for r in cells}:
            raise ValueError(f"{tag}: native endocrine residuals differ")
        ctx = pd.read_parquet(
            args.typing / "context" / f"{tag}.parquet",
            columns=["object_id", "in_islet", "endocrine_object_id"],
        ).set_index("object_id", verify_integrity=True)
        coverage = pd.read_parquet(
            args.typing / "coverage" / f"{tag}.parquet",
            columns=["object_id", *["covint__" + m for m in MARKERS]],
        ).set_index("object_id", verify_integrity=True)
        table = []
        for r in cells:
            pt = Point(r["x"], r["y"])
            cell = shape(r["geometry"])
            near = tree.nearest(pt)
            members = tree.query(pt, predicate="intersects")
            record = {
                "image": tag,
                "object_id": r["id"],
                "native_parent": r["parent"],
                "x_px": r["x"],
                "y_px": r["y"],
                "nearest_id": ann[near]["id"],
                "nearest_class": ann[near]["name"],
                "centroid_distance_um": pt.distance(polys[near]) * sx,
                "native_containing_objects": len(members),
                "overlaps_outline": bool(len(tree.query(cell, predicate="intersects"))),
                "classifier_centroid": bool(mask[int(r["y"] // 2), int(r["x"] // 2)]),
                "classifier_pixels_in_cell": original_mask_support(r["geometry"], mask),
            }
            final_ids = tree.query(cell, predicate="intersects")
            fraction, best = max(
                ((cell.intersection(polys[i]).area / cell.area, i) for i in final_ids),
                default=(0.0, None),
            )
            record["max_final_outline_cell_fraction"] = fraction
            record["overlap_object_id"] = ann[best]["id"] if best is not None else None
            record["overlap_object_class"] = (
                ann[best]["name"] if best is not None else None
            )
            record["centroid_inside_cell"] = cell.covers(pt)
            record["original_outline_centroid"] = bool(
                len(original_tree.query(pt, predicate="intersects"))
            )
            intersect = original_tree.query(cell, predicate="intersects")
            record["max_original_outline_cell_fraction"] = max(
                (cell.intersection(original[i]).area / cell.area for i in intersect),
                default=0.0,
            )
            table.append(record)
        frame = (
            pd.DataFrame(table)
            .set_index("object_id", verify_integrity=True)
            .join(identity)
            .join(coverage.loc[identity.index])
            .join(ctx.loc[identity.index])
        )
        frames.append(frame.reset_index())
    pd.concat(frames, ignore_index=True).to_parquet(args.output, index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in (
        "project",
        "replay",
        "fit",
        "geometry",
        "identity",
        "typing",
        "output",
    ):
        p.add_argument("--" + name, type=Path, required=True)
    audit(p.parse_args())
