#!/usr/bin/env python3
"""Materialize supported broad lineages as selectable schema-4 identities.

This applies the owner-approved current-data policy to rows already marked
``lineage_supported`` by the fitted broad ensemble.  Candidate, unsupported,
and unreviewed rows remain unresolved.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import write_json


LINEAGE_TYPES = {
    "Endocrine": "Hormone-unresolved endocrine cell",
    "Epithelial": "Epithelial lineage",
    "Immune": "Immune lineage",
    "Vascular": "Vascular lineage",
    "Neural": "Neural lineage",
    "Stromal": "Stromal lineage",
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    source = args.input_root / "typing"
    destination = args.output_root / "typing"
    (destination / "cells").mkdir(parents=True, exist_ok=True)

    totals = {}
    for path in sorted((source / "cells").glob("*.parquet")):
        frame = pd.read_parquet(path)
        supported = frame.classification_status.eq("lineage_supported") & frame.lineage.isin(LINEAGE_TYPES)
        frame.loc[supported, "base_type"] = frame.loc[supported, "lineage"].map(LINEAGE_TYPES)
        frame.loc[supported, "classification_status"] = "validated"
        frame.loc[supported, "is_endocrine"] = frame.loc[supported, "lineage"].eq("Endocrine")
        for name, count in frame.loc[supported, "base_type"].value_counts().items():
            totals[name] = totals.get(name, 0) + int(count)
        frame.to_parquet(destination / "cells" / path.name, index=False)

    catalog = json.loads((source / "catalog.json").read_text())
    existing = {entry["name"] for entry in catalog["base_types"]}
    for lineage, name in LINEAGE_TYPES.items():
        if name not in existing:
            catalog["base_types"].append({
                "name": name, "lineage": lineage,
                "description": "Supported broad-lineage assignment retained when no reviewed finer identity is available.",
                "marker_dependencies": [], "computational_features": catalog.get("features", []),
            })
    catalog["base_types"] = sorted(catalog["base_types"], key=lambda row: row["name"])
    catalog["base_type_names"] = [row["name"] for row in catalog["base_types"]]
    write_json(destination / "catalog.json", catalog)

    application = json.loads((source / "application.json").read_text())
    application["lineage_promotion_policy"] = "owner-approved lineage_supported_v1"
    application["promoted_lineage_counts"] = totals
    write_json(destination / "application.json", application)
    print(json.dumps(totals, indent=2))


if __name__ == "__main__":
    main()
