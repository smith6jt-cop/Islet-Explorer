#!/usr/bin/env python3
"""Stage a fresh schema-4 revision with immutable inputs hardlinked.

Trajectory, runtime summaries, eligibility records and release manifests are
omitted because they must be regenerated. Only raw data and images are linked
read-only by policy. Derived files are independent copies, so later ordinary
writes cannot change their parent revision.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

from common import require_fresh_output, sha256, write_json


ROOT_FILES_TO_REBUILD = {
    "candidate.json", "revision.json", "manifest.json", "capacity.json",
    "candidate_validation.json", "image_channels.json",
}
ANALYSIS_TO_REBUILD = {
    "analysis/counts.parquet",
    "analysis/marker_summaries.parquet",
    "analysis/islet_embedding.parquet",
    "analysis/islet_embedding_manifest.json",
    "analysis/classifier_eligibility_manifest.json",
}


def omitted(relative: Path) -> bool:
    text = relative.as_posix()
    if text in ROOT_FILES_TO_REBUILD or text in ANALYSIS_TO_REBUILD:
        return True
    if relative.parts[0] in {"preview", "research", "summary", "trajectory"}:
        return True
    return len(relative.parts) >= 2 and relative.parts[:2] == ("analysis", "composition_parts")


def stage(source: Path, target: Path) -> dict:
    source = source.resolve(strict=True)
    target = require_fresh_output(target, [source])
    identity = json.loads((source / "candidate.json").read_text())
    manifest = json.loads((source / "manifest.json").read_text())
    if str(identity.get("schema_version")) != "4" or str(manifest.get("schema_version")) != "4":
        raise ValueError("Schema-4 parent candidate required")
    if manifest.get("expected_images") != 21:
        raise ValueError("The parent is not the explicit 21-donor candidate")

    paths = [p for p in sorted(source.rglob("*")) if not omitted(p.relative_to(source))]
    for path in paths:
        if path.is_symlink() and (not path.is_file() or path.is_dir()):
            raise ValueError(f"Directory or broken symlink cannot be staged: {path}")
    target.mkdir(parents=True)
    linked_files = 0
    linked_bytes = 0
    for path in paths:
        relative = path.relative_to(source)
        if omitted(relative):
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

    candidate = dict(identity)
    candidate.update(
        schema_version="4",
        revision_of=str(source),
        scientific_status="classifier_eligibility_and_trajectory_rebuild",
        scientific_acceptance=False,
        ufrc_accepted=False,
        research_results_in_app=False,
    )
    write_json(target / "candidate.json", candidate)
    (target / "provenance").mkdir(exist_ok=True)
    write_json(target / "provenance/parent_release.json", {
        "parent": str(source),
        "parent_candidate_sha256": sha256(source / "candidate.json"),
        "parent_manifest_sha256": sha256(source / "manifest.json"),
    })
    revision = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "parent": str(source),
        "schema_version": 4,
        "status": "reprocessing",
        "islet_eligibility": "no image-edge contact and at least one canonical cell inside the endocrine pixel-classifier core",
        "endocrine_denominator": "all canonical cells inside the endocrine pixel-classifier core, counted once",
        "downstream_typing_can_revoke_islet_eligibility": False,
        "linked_files": linked_files,
        "linked_logical_bytes": linked_bytes,
        "raw_inputs_preserved": True,
    }
    write_json(target / "revision.json", revision)
    return revision


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(stage(args.input_root, args.output_root), indent=2))
