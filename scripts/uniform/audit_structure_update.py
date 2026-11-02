#!/usr/bin/env python3
"""Audit the invariants of a Lymphatic-only update using actual run artifacts.

No scientific inputs are modified. UUIDs, REDSEA values, QC, calibration files,
every positive/canary/fold field, and non-Lymphatic context must remain exact.
Derived identity fields may change only inside the union of the old/new
Lymphatic centroid memberships or positive cell-area overlaps. Every changed
label must involve Endothelium Lymphatic. The named endocrine set, label,
compartment and lineage are protected everywhere.

This is a technical change-scope gate, not a biological accuracy assessment.
Exit 0 means the declared invariants passed; exit 1 means they did not. Missing
files, UUID/schema mismatches and incomplete calibration fail closed.

Example:
  python scripts/uniform/audit_structure_update.py --baseline OLD_RUN \
      --candidate NEW_RUN --output FRESH_AUDIT_DIRECTORY --workers 2
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

NAMED_ENDOCRINE = frozenset(("Beta cell", "Alpha cell", "Delta cell", "CHGA+ hormone-negative"))
LY_CONTEXT = frozenset(("in_Lymphatic", "crossed_frac_lymphatic", "object_context", "context"))
IDENTITY_DERIVED = frozenset((
    "label", "named_by", "compartment", "compartment_basis", "compartment_status",
    "n_compartments_canary", "n_compartments_plain", "compartment_depth",
    "runner_up_compartment", "runner_up_depth", "lineage", "lineage_basis",
    "lineage_status", "subtype_status", "status", "tags", "confirmed", "confirming_markers",
))
ENDOCRINE_PROTECTED = frozenset(("label", "compartment", "lineage"))
CALIBRATION_REQUIRED = frozenset((
    "levels.csv", "levels_differentiation.csv", "canary_levels.csv", "display_levels.csv",
))
CONTEXT_REQUIRED = frozenset((
    "qc_analysis_eligible", "x_um", "y_um", "in_islet", "islet_num",
    "endocrine_object_id", "endocrine_object_class", "in_Lymphatic", "crossed_frac_lymphatic",
    "in_Blood vessel", "crossed_frac_vessel", "in_Nerve", "crossed_frac_nerve",
    "object_context", "context",
))
IDENTITY_REQUIRED = frozenset(("label", "compartment", "lineage", "start_class"))


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, default=str, allow_nan=False) + "\n")


def stat_identity(path: Path) -> dict:
    st = path.stat()
    return {"size": st.st_size, "device": st.st_dev, "inode": st.st_ino,
            "mtime_ns": st.st_mtime_ns, "ctime_ns": st.st_ctime_ns}


def fingerprint(path: Path, cache: dict) -> dict:
    before = stat_identity(path)
    key = tuple(before.items())
    if key not in cache:
        h = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(8 << 20), b""):
                h.update(block)
        cache[key] = h.hexdigest()
    if stat_identity(path) != before:
        raise ValueError(f"Input changed during hashing: {path}")
    return {"path": str(path.resolve()), "sha256": cache[key], **before}


def changed(a: pd.Series, b: pd.Series) -> np.ndarray:
    """Exact value equality, with matching missing values equal; no rounding."""
    return ~(a.eq(b) | (a.isna() & b.isna())).fillna(False).to_numpy(dtype=bool)


def read_table(path: Path, columns=None) -> pd.DataFrame:
    frame = pd.read_parquet(path, columns=columns)
    if "object_id" not in frame or frame.object_id.isna().any() or frame.object_id.duplicated().any():
        raise ValueError(f"Missing, null or duplicate object_id: {path}")
    return frame.set_index("object_id")


def field_signatures(matrix: np.ndarray, names: list[str]) -> np.ndarray:
    """Factor repeated change patterns before making human-readable field lists."""
    packed = np.packbits(matrix, axis=1, bitorder="little")
    unique, inverse = np.unique(packed, axis=0, return_inverse=True)
    bits = np.unpackbits(unique, axis=1, count=len(names), bitorder="little")
    labels = np.array([";".join(n for n, hit in zip(names, row) if hit) for row in bits], dtype=object)
    return labels[inverse]


def calibration_names(root: Path) -> set[str]:
    return {p.name for p in (root / "typing").glob("*.csv")
            if "level" in p.name or p.name == "empty_area_nulls.csv"}


def audit_calibration(baseline: Path, candidate: Path, output: Path) -> list[dict]:
    records, cache = [], {}
    names = CALIBRATION_REQUIRED | calibration_names(baseline) | calibration_names(candidate)
    for name in sorted(names):
        pa, pb = (r / "typing" / name for r in (baseline, candidate))
        rec = {"file": name, "baseline_present": pa.is_file(), "candidate_present": pb.is_file()}
        if not pa.is_file() or not pb.is_file():
            rec.update(passed=False, reason="missing_calibration_file")
            records.append(rec)
            continue
        rec["baseline"], rec["candidate"] = fingerprint(pa, cache), fingerprint(pb, cache)
        rec["passed"] = rec["baseline"]["sha256"] == rec["candidate"]["sha256"]
        if not rec["passed"]:
            rec["reason"] = "calibration_bytes_changed"
            a, b = (pd.read_csv(p, keep_default_na=False) for p in (pa, pb))
            rec["baseline_rows"], rec["candidate_rows"] = len(a), len(b)
            rec["columns_equal"] = set(a.columns) == set(b.columns)
            keys = [k for k in ("image", "marker") if k in a and k in b]
            if keys and not a.duplicated(keys).any() and not b.duplicated(keys).any():
                a, b = a.set_index(keys).sort_index(), b.set_index(keys).sort_index()
                rec["keys_equal"] = a.index.equals(b.index)
                if rec["keys_equal"] and rec["columns_equal"]:
                    b = b.reindex(columns=a.columns)
                    differences = {c: changed(a[c], b[c]) for c in a}
                    rec["changed_values_by_column"] = {c: int(d.sum()) for c, d in differences.items()}
                    mask = np.column_stack(list(differences.values())).any(axis=1) if differences else np.zeros(len(a), bool)
                    rec["changed_value_rows"] = int(mask.sum())
                    a.add_suffix("_baseline").join(b.add_suffix("_candidate")).loc[mask].to_csv(
                        output / "calibration" / name)
        records.append(rec)
    return records


def aligned_pair(a: pd.DataFrame, b: pd.DataFrame, kind: str, issues: list, out: Path):
    missing, added = a.index.difference(b.index), b.index.difference(a.index)
    if len(missing) or len(added):
        issues.append({"kind": kind, "reason": "uuid_universe_changed",
                       "missing": len(missing), "added": len(added)})
        pd.DataFrame({"object_id": list(missing) + list(added),
                      "change": ["missing_from_candidate"] * len(missing) + ["added_in_candidate"] * len(added)}).to_parquet(
            out / f"{kind}_uuid_changes.parquet", index=False)
    if set(a.columns) != set(b.columns):
        issues.append({"kind": kind, "reason": "schema_columns_changed",
                       "removed": sorted(set(a) - set(b)), "added": sorted(set(b) - set(a))})
    shared = a.index[a.index.isin(b.index)]
    columns = [c for c in a if c in b]
    a, b = a.loc[shared, columns], b.reindex(shared)[columns]
    dtypes = [c for c in columns if a[c].dtype != b[c].dtype]
    if dtypes:
        issues.append({"kind": kind, "reason": "column_dtypes_changed", "fields": dtypes})
    return a, b


def audit_image(baseline: Path, candidate: Path, output: Path, image: str,
                component_scope: bool = False) -> dict:
    dest = output / "images" / image
    dest.mkdir()
    result = {"image": image, "passed": False, "issues": [], "inputs": {}, "field_changes": []}
    cache = {}
    try:
        tables = {}
        for kind in ("context", "identity", "redsea"):
            pa, pb = (r / "typing" / kind / f"{image}.parquet" for r in (baseline, candidate))
            fa, fb = fingerprint(pa, cache), fingerprint(pb, cache)
            result["inputs"][kind] = {"baseline": fa, "candidate": fb,
                                      "byte_identical": fa["sha256"] == fb["sha256"]}
            if kind == "redsea" and fa["sha256"] == fb["sha256"]:
                # Equal bytes prove every stored measurement/flag equal. Read IDs/QC
                # as well to verify the cross-table universe, even on this fast path.
                cols = ["object_id", "qc_analysis_eligible"]
                a = read_table(pa, cols)
                b = a.copy(deep=False)
                result["redsea_all_columns"] = [c for c in pq.read_schema(pa).names if c != "object_id"]
            else:
                a = read_table(pa)
                b = a.copy(deep=False) if pa.resolve() == pb.resolve() else read_table(pb)
            tables[kind] = (a, b)

        for kind, required in (("context", CONTEXT_REQUIRED), ("identity", IDENTITY_REQUIRED),
                               ("redsea", {"qc_analysis_eligible"})):
            for version, table in zip(("baseline", "candidate"), tables[kind]):
                if required - set(table):
                    raise ValueError(f"{version} {kind}: required fields absent: {sorted(required - set(table))}")
        for version, table in zip(("baseline", "candidate"), tables["identity"]):
            marker_sets = [{c.removeprefix(prefix) for c in table if c.startswith(prefix)}
                           for prefix in ("pos__", "canary__", "fold__")]
            if not marker_sets[0] or not all(s == marker_sets[0] for s in marker_sets):
                raise ValueError(f"{version}: incomplete positive/canary/fold marker fields")
            if not {"INS", "GCG", "SST", "CHGA"} <= marker_sets[0]:
                raise ValueError(f"{version}: endocrine marker fields are incomplete")
        for version, ci in (("baseline", 0), ("candidate", 1)):
            ctx, ident, rs = (tables[k][ci] for k in ("context", "identity", "redsea"))
            if not ctx.index.sort_values().equals(rs.index.sort_values()):
                result["issues"].append({"reason": "context_redsea_uuid_mismatch", "version": version})
            elif changed(ctx.qc_analysis_eligible, rs.reindex(ctx.index).qc_analysis_eligible).any():
                result["issues"].append({"reason": "context_redsea_qc_mismatch", "version": version})
            if ctx.qc_analysis_eligible.isna().any():
                raise ValueError(f"{version}: missing QC eligibility")
            eligible = ctx.index[ctx.qc_analysis_eligible.astype(bool)]
            if not eligible.sort_values().equals(ident.index.sort_values()):
                result["issues"].append({"reason": "identity_not_exactly_qc_eligible_uuid_universe", "version": version,
                                         "eligible": len(eligible), "identity_rows": len(ident)})

        ca, cb = aligned_pair(*tables["context"], "context", result["issues"], dest)
        if ca.in_Lymphatic.isna().any() or cb.in_Lymphatic.isna().any():
            raise ValueError("Missing Lymphatic centroid membership")
        footprint = (ca.in_Lymphatic.astype(bool) | cb.in_Lymphatic.astype(bool)
                     | ca.crossed_frac_lymphatic.gt(0) | cb.crossed_frac_lymphatic.gt(0))
        geometry_changed = changed(ca.in_Lymphatic, cb.in_Lymphatic) | changed(ca.crossed_frac_lymphatic, cb.crossed_frac_lymphatic)
        protected = tables["identity"][0].label.isin(NAMED_ENDOCRINE)
        candidate_endocrine = tables["identity"][1].label.isin(NAMED_ENDOCRINE)
        result["counts"] = {"context_baseline": len(tables["context"][0]), "context_candidate": len(tables["context"][1]),
                            "identity_baseline": len(tables["identity"][0]), "identity_candidate": len(tables["identity"][1]),
                            "redsea_baseline": len(tables["redsea"][0]), "redsea_candidate": len(tables["redsea"][1]),
                            "lymphatic_footprint_union": int(footprint.sum()),
                            "lymphatic_membership_or_overlap_changed": int(geometry_changed.sum()),
                            "baseline_named_endocrine": int(protected.sum())}

        for kind in ("context", "identity", "redsea"):
            a, b = (ca, cb) if kind == "context" else aligned_pair(*tables[kind], kind, result["issues"], dest)
            allowed_area = footprint.reindex(a.index, fill_value=False).to_numpy(bool)
            endo = protected.reindex(a.index, fill_value=False).to_numpy(bool)
            endocrine_union = endo | candidate_endocrine.reindex(a.index, fill_value=False).to_numpy(bool)
            columns = list(a)
            if kind == "redsea" and result["inputs"][kind]["byte_identical"]:
                columns = result["redsea_all_columns"]
                differences = np.zeros((len(a), len(columns)), dtype=bool)
            else:
                differences = np.column_stack([changed(a[c], b[c]) for c in columns])
            unexpected = differences.copy()
            for ci, column in enumerate(columns):
                allowed = kind == "context" and column in LY_CONTEXT or kind == "identity" and column in IDENTITY_DERIVED
                # A separately validated pair of component-selection runs also
                # retains the natural decisions before selection. Apply the same
                # local change restrictions to those diagnostic decisions. Never
                # permit a general natural__ prefix or change to source evidence.
                natural_field = column.removeprefix("natural__") if column.startswith("natural__") else None
                if component_scope and kind == "identity":
                    allowed |= natural_field in IDENTITY_DERIVED or column == "component_decision_source"
                if allowed:
                    unexpected[:, ci] &= ~allowed_area
                if kind == "identity" and column in ENDOCRINE_PROTECTED:
                    unexpected[:, ci] |= differences[:, ci] & endocrine_union
                if kind == "identity" and (column == "label" or component_scope and column == "natural__label"):
                    involves_lymphatic = (a[column].eq("Endothelium Lymphatic")
                                           | b[column].eq("Endothelium Lymphatic")).to_numpy(bool)
                    unexpected[:, ci] |= differences[:, ci] & ~involves_lymphatic
                    result["counts"][f"{column}_changes_neither_lymphatic"] = int(
                        (differences[:, ci] & ~involves_lymphatic).sum())
                    changed_geometry = pd.Series(geometry_changed, index=ca.index).reindex(a.index, fill_value=False).to_numpy(bool)
                    unexpected[:, ci] |= differences[:, ci] & ~changed_geometry
                    result["counts"][f"{column}_changes_with_unchanged_lymphatic_geometry"] = int(
                        (differences[:, ci] & ~changed_geometry).sum())
                # No marker flag/fold exception, including PDPN. Unrecognized
                # columns remain immutable rather than implicitly authorized.
                d, bad = differences[:, ci], unexpected[:, ci]
                record = {"image": image, "table": kind, "field": column, "changed": int(d.sum()),
                          "inside_lymphatic_footprint": int((d & allowed_area).sum()),
                          "outside_lymphatic_footprint": int((d & ~allowed_area).sum()),
                          "baseline_named_endocrine": int((d & endo).sum()), "unexpected": int(bad.sum())}
                if kind == "identity" and (column.startswith(("pos__", "canary__")) or column == "Ki67_pos"):
                    record["false_to_true"] = int((d & a[column].eq(0).to_numpy() & b[column].eq(1).to_numpy()).sum())
                    record["true_to_false"] = int((d & a[column].eq(1).to_numpy() & b[column].eq(0).to_numpy()).sum())
                result["field_changes"].append(record)
            any_change, any_bad = differences.any(axis=1), unexpected.any(axis=1)
            result["counts"][f"{kind}_changed_cells"] = int(any_change.sum())
            result["counts"][f"{kind}_unexpected_cells"] = int(any_bad.sum())
            if any_change.any():
                evidence = pd.DataFrame({"object_id": a.index[any_change],
                                         "inside_lymphatic_footprint": allowed_area[any_change],
                                         "baseline_named_endocrine": endo[any_change],
                                         "changed_fields": field_signatures(differences[any_change], columns),
                                         "unexpected_fields": field_signatures(unexpected[any_change], columns)})
                if kind == "identity":
                    for field in ("label", "compartment", "lineage"):
                        evidence[f"{field}_baseline"] = a.loc[any_change, field].to_numpy()
                        evidence[f"{field}_candidate"] = b.loc[any_change, field].to_numpy()
                    transitions = pd.DataFrame({"baseline": a.label, "candidate": b.label,
                                                "inside_lymphatic_footprint": allowed_area,
                                                "baseline_named_endocrine": endo})
                    transitions.loc[changed(a.label, b.label)].value_counts(dropna=False).rename("cells").reset_index().to_csv(
                        dest / "label_transitions.csv", index=False)
                evidence.to_parquet(dest / f"{kind}_changes.parquet", index=False, compression="zstd")
            if any_bad.any():
                result["issues"].append({"reason": "unexpected_field_changes", "table": kind,
                                         "cells": int(any_bad.sum()),
                                         "fields": [c for c, n in zip(columns, unexpected.sum(axis=0)) if n]})
        for pair in result["inputs"].values():
            for version in ("baseline", "candidate"):
                rec = pair[version]
                current = stat_identity(Path(rec["path"]))
                if current != {k: rec[k] for k in current}:
                    raise ValueError(f"Input changed during audit: {rec['path']}")
        result["passed"] = not result["issues"]
    except Exception as exc:
        result["issues"].append({"reason": "input_or_audit_error", "error": f"{type(exc).__name__}: {exc}"})
    write_json(dest / "summary.json", result)
    return result


def validate_component_pair(baseline: Path, candidate: Path) -> dict | None:
    """Natural decision exceptions require two complete, matching component contracts."""
    present = [(r / "typing/component_selection.json").exists() for r in (baseline, candidate)]
    if not any(present):
        return None
    if not all(present):
        raise ValueError("Both sides must have the same explicit component-selection scope")
    sys.path.insert(0, str(Path(__file__).resolve().parent / "typing"))
    import component_selection as CS
    bindings = [CS.validate(r) for r in (baseline, candidate)]
    audits = [CS.validate_audit(r) for r in (baseline, candidate)]
    a, b = bindings
    for field in ("scope", "endocrine_source", "nonendocrine_source", "endocrine_markers", "endocrine_labels"):
        if a[field] != b[field]:
            raise ValueError(f"Component control and candidate differ in {field}")
    if set(a["images"]) != set(b["images"]):
        raise ValueError("Component control and candidate have different donor inventories")
    keys = {k for k in a["files"] if k.startswith(("endocrine/", "nonendocrine/", "registry/", "code/"))}
    if keys != {k for k in b["files"] if k.startswith(("endocrine/", "nonendocrine/", "registry/", "code/"))}:
        raise ValueError("Component control and candidate have different source inventories")
    for key in keys:
        if a["files"][key]["sha256"] != b["files"][key]["sha256"]:
            raise ValueError(f"Component control and candidate source differs: {key}")
    return {"schema": a["schema"], "validated": True,
            "endocrine_source": a["endocrine_source"], "nonendocrine_source": a["nonendocrine_source"],
            "manifest_sha256": [CS.sha(r / CS.MANIFEST) for r in (baseline, candidate)],
            "audit_sha256": [CS.sha(r / CS.AUDIT) for r in (baseline, candidate)],
            "selected_named_endocrine_cells": [x["preserved_named_endocrine_cells"] for x in audits]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--images", nargs="+", default=["all"], help="all, or explicit image names (limited audit)")
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    baseline, candidate, output = args.baseline.resolve(), args.candidate.resolve(), args.output.resolve()
    if args.workers < 1:
        parser.error("workers must be positive")
    if output.exists():
        parser.error("output already exists; choose a fresh audit directory")
    if any(output == r or output.is_relative_to(r) for r in (baseline, candidate)):
        parser.error("audit output must be outside both scientific input roots")
    output.mkdir(parents=True)
    (output / "images").mkdir()
    (output / "calibration").mkdir()
    started = datetime.now(timezone.utc).isoformat()
    script_source = fingerprint(Path(__file__), {})
    issues = []
    component_pair = None
    try:
        component_pair = validate_component_pair(baseline, candidate)
    except Exception as exc:
        issues.append({"reason": "component_contract_invalid", "error": f"{type(exc).__name__}: {exc}"})
    inventories = {}
    for version, root in (("baseline", baseline), ("candidate", candidate)):
        inventories[version] = {k: sorted(p.stem for p in (root / "typing" / k).glob("*.parquet"))
                                for k in ("context", "identity", "redsea")}
    available = set().union(*(set(x) for inv in inventories.values() for x in inv.values()))
    images = sorted(available if args.images == ["all"] else set(args.images))
    if any(not x or Path(x).name != x or x in (".", "..") for x in images):
        parser.error("images must be individual image names, not paths")
    selected_inventory = {version: {kind: sorted(set(tags) & set(images))
                                    for kind, tags in inv.items()}
                          for version, inv in inventories.items()}
    for version, inv in selected_inventory.items():
        if len({tuple(v) for v in inv.values()}) != 1:
            issues.append({"reason": "within_run_selected_donor_inventory_mismatch", "version": version})
    if selected_inventory["baseline"] != selected_inventory["candidate"]:
        issues.append({"reason": "baseline_candidate_selected_donor_inventory_mismatch"})
    if not images:
        issues.append({"reason": "no_images"})
    if set(images) - available:
        issues.append({"reason": "requested_images_missing", "images": sorted(set(images) - available)})
    try:
        calibration = audit_calibration(baseline, candidate, output)
    except Exception as exc:
        calibration = []
        issues.append({"reason": "calibration_audit_error", "error": f"{type(exc).__name__}: {exc}"})
    if any(not x["passed"] for x in calibration):
        issues.append({"reason": "calibration_not_immutable",
                       "files": [x["file"] for x in calibration if not x["passed"]]})
    results = []
    with ProcessPoolExecutor(args.workers) as pool:
        futures = [pool.submit(audit_image, baseline, candidate, output, image, component_pair is not None)
                   for image in images]
        for image, future in zip(images, futures):
            result = future.result()
            results.append(result)
            print(f"{image}: {'PASS' if result['passed'] else 'FAIL'}", flush=True)
    fields = [x for r in results for x in r["field_changes"]]
    if fields:
        frame = pd.DataFrame(fields)
        frame.to_csv(output / "field_changes_by_image.csv", index=False)
        frame.groupby(["table", "field"], dropna=False).sum(numeric_only=True).reset_index().to_csv(
            output / "field_changes_totals.csv", index=False)
    totals = Counter()
    for result in results:
        totals.update(result.get("counts", {}))
    for record in calibration:
        for version in ("baseline", "candidate"):
            if version in record:
                prior = record[version]
                current = stat_identity(Path(prior["path"]))
                if current != {k: prior[k] for k in current}:
                    issues.append({"reason": "calibration_input_changed_during_audit",
                                   "version": version, "file": record["file"]})
    if fingerprint(Path(__file__), {})["sha256"] != script_source["sha256"]:
        issues.append({"reason": "audit_source_changed_during_execution"})
    passed = not issues and bool(results) and all(r["passed"] for r in results)
    summary = {"schema": "uniform.structure_update_invariant_audit/1", "passed": passed,
               "target_structure": "Lymphatic", "baseline": str(baseline), "candidate": str(candidate),
               "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "script_sha256": script_source["sha256"],
               "component_selection_pair": component_pair,
               "complete_donor_inventory_audited": set(images) == available,
               "images": images, "donor_inventory": inventories,
               "selected_donor_inventory": selected_inventory, "issues": issues,
               "calibration": calibration, "totals": dict(totals),
               "per_image": [{k: v for k, v in r.items() if k not in ("inputs", "field_changes")} for r in results],
               "rules": {"lymphatic_footprint": "old/new in_Lymphatic OR old/new crossed_frac_lymphatic > 0",
                         "calibration": "all required and discovered level/null CSV files must be byte-identical",
                         "markers": "all pos__/canary__/fold__ fields immutable everywhere, including PDPN",
                         "redsea": "every measurement and QC field immutable everywhere",
                         "non_target_context": "every context field except the four named Lymphatic fields immutable",
                         "named_endocrine": "named endocrine set, label, compartment and lineage immutable everywhere",
                         "expected_identity_changes": "only listed derived identity fields inside the Lymphatic footprint; every changed label must involve Endothelium Lymphatic",
                         "unknown_fields": "immutable"},
               "limitations": ["A passing change-scope audit is not biological accuracy or scientific acceptance.",
                               "The footprint uses existing centroid/area-overlap fields; no new geometry is inferred.",
                               "An explicit image subset is a limited audit, even when its invariants pass."]}
    write_json(output / "summary.json", summary)
    print(f"{'PASS' if passed else 'FAIL'}: {len(images)} images; report {output / 'summary.json'}", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
