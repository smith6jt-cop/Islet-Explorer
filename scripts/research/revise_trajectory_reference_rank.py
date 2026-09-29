#!/usr/bin/env python3
"""Create a fresh schema-4 candidate with fixed evenly spaced display ranks.

The fitted diffusion graph, raw pseudotime, components, typing, measurements,
and images are preserved byte-for-byte.  Only the historical
``pseudotime_percentile`` display coordinate and the metadata/model fields that
define it are replaced.  Exact raw-coordinate ties are ordered by SHA-256 of
the immutable islet identifier, independently of donor or biological values.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from common import require_fresh_output, sha256, write_json
from trajectory import fixed_reference_rank


ROOT_FILES_TO_REBUILD = {
    "candidate.json", "revision.json", "manifest.json", "capacity.json",
    "candidate_validation.json", "image_channels.json",
}
SKIP_ROOTS = {"preview", "evidence"}


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=f".{path.name}.", delete=False) as stream:
        temporary = Path(stream.name)
        frame.to_csv(stream, index=False)
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", delete=False) as stream:
        temporary = Path(stream.name)
    try:
        frame.to_parquet(temporary, index=False)
        temporary.replace(path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def atomic_npz(values: dict[str, np.ndarray], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", delete=False) as stream:
        temporary = Path(stream.name)
        np.savez_compressed(stream, **values)
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def coordinate_identity(frame: pd.DataFrame) -> str:
    digest = hashlib.sha256()
    for row in frame[["islet_id", "pseudotime_raw"]].itertuples(index=False):
        digest.update(str(row.islet_id).encode())
        digest.update(b"\0")
        digest.update(np.asarray([row.pseudotime_raw], dtype="<f8").tobytes())
    return digest.hexdigest()


def clone(source: Path, target: Path) -> tuple[int, int]:
    source=source.resolve(strict=True)
    target=require_fresh_output(target,[source])
    linked_files = linked_bytes = 0
    paths=[p for p in sorted(source.rglob("*"))
           if p.relative_to(source).as_posix() not in ROOT_FILES_TO_REBUILD
           and p.relative_to(source).parts[0] not in SKIP_ROOTS]
    for path in paths:
        if path.is_symlink() and not path.is_file():
            raise ValueError(f"Directory or broken symlink cannot be staged: {path}")
    target.mkdir(parents=True)
    for path in paths:
        relative = path.relative_to(source)
        if relative.as_posix() in ROOT_FILES_TO_REBUILD or relative.parts[0] in SKIP_ROOTS:
            continue
        destination = target / relative
        if path.is_dir():
            destination.mkdir(parents=True, exist_ok=True)
        elif path.is_file():
            destination.parent.mkdir(parents=True, exist_ok=True)
            if relative.parts[0] in {"raw", "images"}:
                os.link(path.resolve(strict=True), destination)
                linked_files += 1
                linked_bytes += path.stat().st_size
            else:
                shutil.copy2(path, destination)
    return linked_files, linked_bytes


def revise(source: Path, target: Path) -> dict:
    source = source.resolve(strict=True)
    target = require_fresh_output(target, [source])
    identity = json.loads((source / "candidate.json").read_text())
    manifest = json.loads((source / "manifest.json").read_text())
    if str(identity.get("schema_version")) != "4" or str(manifest.get("schema_version")) != "4":
        raise ValueError("Schema-4 source candidate required")
    if manifest.get("expected_images") != 21:
        raise ValueError("The source is not the explicit 21-donor candidate")

    csv = pd.read_csv(source / "trajectory/coordinates.csv", dtype={"case_id": str})
    parquet = pd.read_parquet(source / "trajectory/coordinates.parquet")
    if "case_id" in parquet:
        parquet["case_id"] = parquet.case_id.astype(str)
    csv_compare, parquet_compare = csv.copy(), parquet.copy()
    for column in parquet.columns:
        if pd.api.types.is_object_dtype(parquet[column]) or pd.api.types.is_string_dtype(parquet[column]):
            csv_compare[column] = csv_compare[column].replace("", np.nan)
            parquet_compare[column] = parquet_compare[column].replace("", np.nan)
    pd.testing.assert_frame_equal(csv_compare, parquet_compare, check_dtype=False)
    coordinates = parquet
    if coordinates.islet_id.duplicated().any():
        raise ValueError("Trajectory islet identifiers must be unique")
    finite = coordinates.pseudotime_raw.notna()
    if not finite.any():
        raise ValueError("No supported primary coordinates to rank")

    before_identity = coordinate_identity(coordinates)
    ranks = fixed_reference_rank(coordinates.pseudotime_raw.to_numpy(float), coordinates.islet_id.astype(str).to_numpy())
    revised = coordinates.copy()
    revised["pseudotime_percentile"] = ranks
    if not revised.loc[~finite, "pseudotime_percentile"].isna().all():
        raise AssertionError("Unsupported components received display ranks")
    expected = (np.arange(finite.sum(), dtype=float) + .5) / finite.sum()
    if not np.allclose(np.sort(revised.loc[finite, "pseudotime_percentile"]), expected, rtol=0, atol=1e-12):
        raise AssertionError("Reference rank is not the evenly spaced midpoint grid")
    unchanged = [column for column in coordinates if column != "pseudotime_percentile"]
    pd.testing.assert_frame_equal(coordinates[unchanged], revised[unchanged], check_dtype=False)

    linked_files, linked_bytes = clone(source, target)
    atomic_csv(revised, target / "trajectory/coordinates.csv")
    atomic_parquet(revised, target / "trajectory/coordinates.parquet")

    with np.load(source / "trajectory/model.npz", allow_pickle=True) as saved:
        model = {name: saved[name] for name in saved.files}
    model["reference_raw"] = np.sort(revised.loc[finite, "pseudotime_raw"].to_numpy(float))
    model["reference_ids"] = revised.loc[finite, "islet_id"].astype(str).to_numpy()
    model["reference_id_raw"] = revised.loc[finite, "pseudotime_raw"].to_numpy(float)
    model["reference_id_percentile"] = revised.loc[finite, "pseudotime_percentile"].to_numpy(float)
    atomic_npz(model, target / "trajectory/model.npz")

    raw_counts = revised.loc[finite, "pseudotime_raw"].value_counts()
    specification = json.loads((source / "trajectory/specification.json").read_text())
    specification.update(
        reference_percentile="fixed ordinal display rank: raw pseudotime order with exact ties spread by SHA-256(islet_id); midpoint ranks (i+0.5)/N",
        reference_percentile_version=2,
        reference_percentile_role="display and descriptive binning only; inference uses pseudotime_raw and preserves exact ties",
        reference_rank_count=int(finite.sum()),
        raw_pseudotime_unique=int(revised.loc[finite, "pseudotime_raw"].nunique()),
        largest_raw_pseudotime_tie=int(raw_counts.max()),
        runner_sha256=sha256(Path(__file__).with_name("trajectory.py")),
    )
    write_json(target / "trajectory/specification.json", specification)

    validation = json.loads((source / "trajectory/validation.json").read_text())
    validation["reference_rank"] = {
        "version": 2,
        "supported_islets": int(finite.sum()),
        "unique_display_ranks": int(revised.loc[finite, "pseudotime_percentile"].nunique()),
        "raw_coordinate_unique_values": int(revised.loc[finite, "pseudotime_raw"].nunique()),
        "largest_exact_raw_tie": int(raw_counts.max()),
        "within_tie_order": "SHA-256(islet_id), for display only and without biological meaning",
        "inference_coordinate": "pseudotime_raw",
    }
    write_json(target / "trajectory/validation.json", validation)

    identity.update(
        revision_of=str(source),
        scientific_status="classifier_eligibility_trajectory_reference_rank_v2",
        scientific_acceptance=False,
        ufrc_accepted=False,
        research_results_in_app=False,
    )
    write_json(target / "candidate.json", identity)
    after_identity = coordinate_identity(revised)
    if before_identity != after_identity:
        raise AssertionError("Raw pseudotime or islet alignment changed")
    provenance = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "parent": str(source),
        "source_manifest_sha256": sha256(source / "manifest.json"),
        "schema_version": 4,
        "kind": "trajectory_reference_rank_v2",
        "status": "awaiting_finalization_and_validation",
        "linked_files": linked_files,
        "linked_logical_bytes": linked_bytes,
        "raw_pseudotime_preserved": True,
        "raw_coordinate_identity_sha256": before_identity,
        "supported_islets": int(finite.sum()),
        "display_rank_method": specification["reference_percentile"],
        "display_rank_role": specification["reference_percentile_role"],
        "largest_exact_raw_tie": int(raw_counts.max()),
        "source_coordinates_sha256": sha256(source / "trajectory/coordinates.parquet"),
        "revised_coordinates_sha256": sha256(target / "trajectory/coordinates.parquet"),
    }
    write_json(target / "revision.json", provenance)
    write_json(target / "provenance/trajectory_reference_rank_v2.json", provenance)
    return provenance


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    arguments = parser.parse_args()
    print(json.dumps(revise(arguments.input_root, arguments.output_root), indent=2))
