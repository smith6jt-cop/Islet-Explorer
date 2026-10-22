#!/usr/bin/env python3
"""Inventory runtime schemas/counts; never serialize donor rows or embedded text.

Use an interpreter with h5py and tifffile (plus pyarrow for a built release):
  python scripts/inventory_release_data.py --source-root data --image-root PATH \
    --parity-report /restricted/evidence/parity.json --output-prefix PATH
Add --release-root PATH after building the immutable release. Embedded TIFF/XML/
JSON free text still requires data-owner inspection; field names cannot establish
that the files contain no identifiers.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import subprocess
import sys
from urllib.parse import quote
import xml.etree.ElementTree as ET


def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def csv_schema(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = csv.reader(stream)
        fields = next(rows)
        return {"fields": fields, "rows": sum(1 for _ in rows)}


def group_csv_schemas(paths: list[Path]) -> dict:
    groups: dict[tuple[str, ...], dict] = {}
    for path in paths:
        value = csv_schema(path)
        key = tuple(value["fields"])
        group = groups.setdefault(key, {"fields": value["fields"], "files": 0, "rows": 0})
        group["files"] += 1
        group["rows"] += value["rows"]
    return {"files": len(paths), "schemas": list(groups.values())}


def h5ad_schema(path: Path) -> dict:
    import h5py
    with h5py.File(path, "r") as handle:
        slots = {}
        for name, value in handle.items():
            entry = {"kind": "group" if isinstance(value, h5py.Group) else "dataset",
                     "attribute_names": sorted(value.attrs.keys())}
            if isinstance(value, h5py.Group):
                entry["field_names"] = sorted(value.keys())
                if name in {"obs", "var"}:
                    index = value.attrs.get("_index", "_index")
                    if isinstance(index, bytes):
                        index = index.decode()
                    entry["rows"] = len(value[index])
                    entry["fields"] = {
                        key: {"encoding": str(child.attrs.get("encoding-type", "")),
                              "dtype": str(child.dtype) if isinstance(child, h5py.Dataset) else "categorical/group"}
                        for key, child in value.items()
                    }
            else:
                entry.update(shape=list(value.shape), dtype=str(value.dtype))
            slots[name] = entry
    return {"bytes": path.stat().st_size, "sha256": sha256(path), "slots": slots,
            "values_included": False}


def xml_field_names(text: str) -> dict:
    root = ET.fromstring(text)
    return {"element_names": sorted({node.tag.rsplit("}", 1)[-1] for node in root.iter()}),
            "attribute_names": sorted({key.rsplit("}", 1)[-1] for node in root.iter() for key in node.attrib}),
            "text_values_included": False}


def json_field_names(value: object) -> dict:
    # Numeric/object-instance keys can themselves contain scientific IDs. Only
    # identifier-shaped schema keys are retained; values are never serialized.
    fields: set[str] = set()
    types: set[str] = set()
    dynamic = False
    def visit(item: object) -> None:
        nonlocal dynamic
        types.add(type(item).__name__)
        if isinstance(item, dict):
            for key, child in item.items():
                if re.fullmatch(r"[A-Za-z_]+", key):
                    fields.add(key)
                else:
                    dynamic = True
                visit(child)
        elif isinstance(item, list):
            for child in item:
                visit(child)
    visit(value)
    return {"field_names": sorted(fields), "value_types": sorted(types),
            "dynamic_keys_omitted": dynamic, "values_included": False}


def image_inventory(root: Path) -> dict:
    import tifffile
    paths = sorted(path for path in root.rglob("*") if path.is_file())
    formats = Counter(path.suffix.lower() or "no_extension" for path in paths)
    tag_names: set[str] = set()
    ome_elements: set[str] = set()
    ome_attributes: set[str] = set()
    text_tags: set[str] = set()
    tiff_count = 0
    page_count = 0
    ome_count = 0
    sidecars: list[dict] = []
    for path in paths:
        suffix = path.suffix.lower()
        if suffix in {".tif", ".tiff"}:
            with tifffile.TiffFile(path) as image:
                tiff_count += 1
                page_count += len(image.pages)
                for page in image.pages:
                    for tag in page.tags.values():
                        tag_names.add(tag.name)
                        if int(tag.dtype) == 2:  # TIFF ASCII; record names only.
                            text_tags.add(tag.name)
                metadata = image.ome_metadata
                if metadata:
                    ome_count += 1
                    schema = xml_field_names(metadata)
                    ome_elements.update(schema["element_names"])
                    ome_attributes.update(schema["attribute_names"])
        elif suffix == ".xml":
            sidecars.append({"format": "xml", **xml_field_names(path.read_text(encoding="utf-8"))})
        elif suffix == ".json":
            sidecars.append({"format": "json", **json_field_names(json.loads(path.read_text(encoding="utf-8")))})
    return {"files": len(paths), "bytes": sum(path.stat().st_size for path in paths),
            "format_counts": dict(formats), "tiff_files": tiff_count, "tiff_pages": page_count,
            "ome_tiff_files": ome_count, "tiff_tag_names": sorted(tag_names),
            "embedded_text_tag_names": sorted(text_tags), "ome_element_names": sorted(ome_elements),
            "ome_attribute_names": sorted(ome_attributes), "sidecar_schemas": sidecars,
            "pixel_arrays_read": False, "embedded_text_values_included": False,
            "owner_inspection_required": "Review embedded text and sidecar values for identifiers; this report inventories names and counts only."}


def parquet_inventory(paths: list[Path]) -> dict:
    import pyarrow.parquet as pq
    groups = {}
    for path in paths:
        meta = pq.ParquetFile(path)
        fields = [(field.name, str(field.type)) for field in meta.schema_arrow]
        key = tuple(fields)
        group = groups.setdefault(key, {"fields": dict(fields), "files": 0, "rows": 0})
        group["files"] += 1
        group["rows"] += meta.metadata.num_rows
    return {"files": len(paths), "schemas": list(groups.values())}


def geometry_inventory(paths: list[Path]) -> dict:
    groups = {}
    for path in paths:
        uri = "file:" + quote(str(path.resolve()), safe="/") + "?mode=ro"
        with sqlite3.connect(uri, uri=True) as connection:
            fields = [(row[1], row[2]) for row in connection.execute("PRAGMA table_info(segmentation)")]
            rows = connection.execute("SELECT count(*) FROM segmentation").fetchone()[0]
            indexed = connection.execute("SELECT count(*) FROM sqlite_master WHERE name LIKE 'rtree_segmentation_%'").fetchone()[0] > 0
        key = tuple(fields)
        group = groups.setdefault(key, {"fields": dict(fields), "files": 0, "features": 0, "indexed_files": 0})
        group["files"] += 1
        group["features"] += rows
        group["indexed_files"] += int(indexed)
    return {"files": len(paths), "schemas": list(groups.values())}


def download_scope(comp_fields: list[str]) -> dict:
    neighborhood = [field for field in comp_fields if re.match(r"^(peri_prop_|peri_count_|enrich_z_|min_dist_|leiden_[0-9])", field)
                    or field in {"immune_frac_peri", "immune_frac_core", "immune_ratio", "cd8_to_macro_ratio",
                                 "tcell_density_peri", "leiden_umap_1", "leiden_umap_2"}]
    return {
        "basis": "Reviewer CSV handlers at the recorded source commit; fields may be omitted for unavailable data or current selection.",
        "plot_summary": ["donor_status", "diam_bin", "diam_mid", "n", "mean", "sd", "se", "median", "q1", "q3", "y", "ymin", "ymax"],
        "plot_distribution": ["Case ID", "islet_key", "donor_status", "value", "islet_diam_um", "total_cells_core", "total_cells_peri", "is_outlier", "outlier_score"],
        "trajectory_heatmap": ["marker", "pseudotime_mid", "bin_mean", "z_clamped", "n_islets_in_bin"],
        "trajectory_umap": ["case_id", "islet_key", "combined_islet_id", "donor_status", "umap_1", "umap_2", "pseudotime", "pt_rank", "feature_name", "feature_value", "total_cells", "islet_diam_um"],
        "trajectory_paga": ["islet_id", "case_id", "donor_status", "leiden", "fa_1", "fa_2", "pseudotime"],
        "spatial_neighborhood": ["Case ID", "Donor Status", "islet_key", "islet_diam_um", "total_cells_core", "total_cells_peri", *neighborhood],
        "spatial_enrichment": {
            "common_fields": ["phenotype", "donor_status", "source_column", "n_donors"],
            "median_fields": ["median_of_donor_means", "q25", "q75"],
            "mean_fields": ["mean_of_donor_means", "mean_minus_sem", "mean_plus_sem"],
        },
        "spatial_infiltration": ["level", "phenotype", "case_id", "donor_status", "islet_key", "enrich_z", "is_outlier", "fence_dist_iqr"],
        "spatial_distance": ["metric", "case_id", "islet_key", "donor_status", "distance_um", "is_outlier", "z_log1p"],
        "csv_comment_headers": "Method, selected scientific filters/features, bin/range settings, row counts, outlier rules and creation time; values excluded here.",
        "png_exports": "Rendered Plot/Trajectory/Spatial/segmentation/drilldown views, including plotted scientific labels; browser Plotly modebar also permits PNG exports.",
        "reviewer_disabled": ["Statistics", "Inflammation", "AI Assistant"],
    }


def build_inventory(args: argparse.Namespace) -> dict:
    repo = Path(__file__).resolve().parent.parent
    source = args.source_root.resolve()
    report = {"schema_version": 1, "created_at_utc": datetime.now(timezone.utc).isoformat(),
              "mode": "release_inventory" if args.release_root else "source_inventory_release_pending",
              "identity": {"git_commit": subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
                           "git_dirty": bool(subprocess.check_output(["git", "-C", str(repo), "status", "--porcelain"], text=True).strip())},
              "source": {}, "release_artifacts": {}}
    for name in ("islet_explorer.h5ad", "adata_ins_root.h5ad"):
        report["source"][name] = h5ad_schema(source / name)
    # A uniform 20-donor root (data schema 5) carries uniform_release.json, the object counts and the
    # protein phenotype rules, and may lack phenotype_rules.csv; the frozen root (schema 1) is the reverse.
    uniform = (source / "uniform_release.json").is_file()
    csv_sources = ["cell_distance.csv", "signed_distance_kde.csv", "barlow_islet_features.csv",
                   "islet_spatial_lookup.csv", "donor_metadata.csv"]
    if uniform:
        csv_sources += ["endocrine_objects_by_donor.csv", "protein_phenotype_rules.csv"]
        if (source / "phenotype_rules.csv").is_file():
            csv_sources.append("phenotype_rules.csv")
        with (source / "uniform_release.json").open(encoding="utf-8") as stream:
            release_keys = sorted(json.load(stream).keys())
        # Key names only; the file's values name run paths and counts, which stay out of the inventory.
        report["source"]["uniform_release.json"] = {"fields": release_keys,
                                                    "sha256": sha256(source / "uniform_release.json")}
    else:
        csv_sources.append("phenotype_rules.csv")
    for name in csv_sources:
        report["source"][name] = {**csv_schema(source / name), "sha256": sha256(source / name)}
    for name in ("cells", "donors"):
        report["source"][name] = group_csv_schemas(sorted((source / name).glob("*.csv")))
    report["images"] = image_inventory(args.image_root.resolve())
    comp_fields = []
    if args.parity_report:
        parity = json.loads(args.parity_report.read_text(encoding="utf-8"))
        if parity["identity"]["canonical_h5ad_sha256"] != report["source"]["islet_explorer.h5ad"]["sha256"]:
            raise ValueError("parity_source_identity_mismatch")
        if parity["identity"]["donor_metadata_sha256"] != report["source"]["donor_metadata.csv"]["sha256"]:
            raise ValueError("parity_metadata_identity_mismatch")
        report["prepared_summaries"] = parity["summaries"]
        report["identity"]["parity_report_sha256"] = sha256(args.parity_report)
        comp_fields = parity["summaries"]["comp"]["fields"]
    if args.release_root:
        release = args.release_root.resolve()
        manifest_path = release / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        report["identity"]["release_manifest_sha256"] = sha256(manifest_path)
        report["identity"]["release_git_commit"] = manifest["git_commit"]
        for name, artifact in manifest["artifacts"].items():
            path = (release / artifact["relative_path"]).resolve()
            if not path.is_relative_to(release):
                raise ValueError("unsafe_manifest_path")
            entry = {key: artifact[key] for key in ("type", "relative_path", "bytes", "sha256", "row_count", "file_count", "columns") if key in artifact}
            if path.suffix == ".parquet":
                entry["actual_inventory"] = parquet_inventory([path])
                if name == "comp":
                    comp_fields = list(entry["actual_inventory"]["schemas"][0]["fields"])
            elif name == "donor_tissue":
                entry["actual_inventory"] = parquet_inventory(sorted(path.rglob("*.parquet")))
            elif name == "geometry":
                entry["actual_inventory"] = geometry_inventory(sorted(path.glob("*.gpkg")))
            elif name == "cells":
                entry["actual_inventory"] = group_csv_schemas(sorted(path.glob("*.csv")))
            elif name == "images":
                entry["actual_inventory"] = image_inventory(path)
            elif path.suffix == ".h5ad":
                entry["actual_inventory"] = h5ad_schema(path)
            elif path.suffix == ".csv":
                entry["actual_inventory"] = csv_schema(path)
            report["release_artifacts"][name] = entry
    report["download_scope"] = download_scope(comp_fields)
    report["download_scope"]["dynamic_neighborhood_fields_resolved"] = bool(comp_fields)
    report["limitations"] = ["No scientific row values or embedded free text are emitted.",
                             "Owner inspection of TIFF/XML/JSON text is still required; a schema inventory is not de-identification certification.",
                             "Use the separate parity/hash report to establish release integrity; this inventory does not replace it."]
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=Path("data"))
    parser.add_argument("--image-root", type=Path, default=Path("app/shiny_app/www/local_images"))
    parser.add_argument("--release-root", type=Path)
    parser.add_argument("--parity-report", type=Path)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    paths = [Path(str(args.output_prefix) + suffix) for suffix in (".json", ".md")]
    if any(path.exists() for path in paths):
        parser.error("Evidence outputs already exist; choose a new prefix.")
    try:
        report = build_inventory(args)
        args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        paths[0].write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        lines = ["# Runtime data inventory", "", f"Mode: {report['mode']}", "",
                 f"Source commit: {report['identity']['git_commit']}", "",
                 "The accompanying JSON lists exact field names, H5AD keys, format counts, TIFF/sidecar metadata fields and reviewer download scope.", "",
                 "| Source table family | Rows |", "|---|---|"]
        for name, item in report["source"].items():
            if "rows" in item:
                lines.append(f"| {name} | {item['rows']} |")
            elif "schemas" in item:
                lines.append(f"| {name} | {sum(group['rows'] for group in item['schemas'])} |")
        lines += ["", f"Image files: {report['images']['files']}; TIFF files: {report['images']['tiff_files']}; pixel arrays read: no.", "",
                  "Embedded free text and sidecar values require data-owner inspection. Their contents are not printed in this inventory."]
        paths[1].write_text("\n".join(lines) + "\n", encoding="utf-8")
    except Exception:
        print("Inventory failed; inspect dependencies and inputs in the restricted environment. Raw errors suppressed.", file=sys.stderr)
        return 1
    print("Runtime data inventory completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
