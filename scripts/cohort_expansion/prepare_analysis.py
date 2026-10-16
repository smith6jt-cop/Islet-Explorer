#!/usr/bin/env python3
"""Import completed QuPath measurements into isolated, unphenotyped AnnData.

This input is for offline analysis only. It deliberately has no `phenotype`
column: QuPath endocrine/immune classifications are not the 21-class phenotype
method used by the frozen application dataset.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import sys

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from excluded_donors import find_excluded


REPO = Path(__file__).resolve().parents[2]
BASELINE = REPO / "single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad"
DEFAULT_ROOT = REPO / "data/development/expanded-20260911"
COHORT = Path(__file__).with_name("cohort.csv")
QUANTIFICATION_PIXEL_UM = 0.3774
DISTANCES = {
    "Islet": "Islet", "IsletExpanded": "Peri-Islet", "Tissue": "Tissue",
    "Nerve": "Nerve", "Lymphatic": "Lymphatic", "Capillary": "Capillary",
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def case_id(value):
    from decimal import Decimal, InvalidOperation
    try:
        number = Decimal(str(value).strip())
    except InvalidOperation as exc:
        raise ValueError(f"Invalid numeric case ID: {value!r}") from exc
    if not number.is_finite() or number <= 0 or number != number.to_integral_value():
        raise ValueError(f"Invalid numeric case ID: {value!r}")
    return str(int(number))


def read_cohort(path):
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    seen = set()
    for row in rows:
        row["case_id"] = case_id(row["case_id"])
        if row["case_id"] in seen:
            raise ValueError("Duplicate case ID in cohort")
        seen.add(row["case_id"])
        if row["donor_source"] not in {"nPOD", "HPAP"}:
            raise ValueError("Explicit nPOD or HPAP provenance is required")
        if row["inclusion_status"] not in {"selected", "identity_hold", "excluded"}:
            raise ValueError("Unrecognized cohort inclusion status")
        if find_excluded([row["case_id"]]) and row["inclusion_status"] != "excluded":
            raise ValueError(f"Case {row['case_id']} is in excluded_donors.csv and must be marked excluded")
    return rows


def donor_metadata(path, selected):
    path = Path(path)
    frame = pd.read_excel(path) if path.suffix.lower() == ".xlsx" else pd.read_csv(path)
    frame = frame.loc[frame["Case ID"].notna()].copy()
    frame["case_id"] = frame["Case ID"].map(case_id)
    if frame["case_id"].duplicated().any():
        raise ValueError("Donor key contains duplicate case IDs")
    frame = frame.set_index("case_id")
    output = {}
    for row in selected:
        cid = row["case_id"]
        if cid not in frame.index:
            raise ValueError(f"No donor-key record for {cid}")
        record = frame.loc[cid]
        if record[["Age", "Gender", "Donor Status"]].isna().any():
            raise ValueError(f"Incomplete demographics for {cid}")
        status = str(record["Donor Status"]).strip()
        status = {"AAb+": "Aab+", "AAb": "Aab+", "Aab": "Aab+"}.get(status, status)
        if status not in {"ND", "Aab+", "T1D"}:
            raise ValueError(f"Unsupported donor status for {cid}: {status!r}")
        age = float(record["Age"])
        if not np.isfinite(age) or age < 0:
            raise ValueError(f"Invalid age for {cid}")
        output[cid] = {"Age": age, "Gender": str(record["Gender"]).strip(),
                       "Donor Status": status, "donor_source": row["donor_source"]}
    return output


def convert_frame(frame, markers, row, metadata, retain_unmeasured=False):
    required = ["Image", "Object ID", "Parent", "Centroid X µm", "Centroid Y µm"]
    marker_columns = [f"Cell: {marker}: Mean" for marker in markers]
    missing = set(required + marker_columns) - set(frame.columns)
    if missing:
        raise ValueError(f"Missing measurements: {sorted(missing)}")
    if frame.empty or frame["Object ID"].isna().any() or frame["Object ID"].duplicated().any():
        raise ValueError("Cell export is empty or has missing/duplicate object IDs")
    images = set(frame["Image"].dropna().map(lambda value: Path(str(value)).name))
    if frame["Image"].isna().any() or images != {row["source_filename"]}:
        raise ValueError(f"Cell export image identity mismatch: {images}")
    raw = frame[marker_columns].apply(pd.to_numeric, errors="raise").to_numpy(dtype=np.float32)
    unmeasured = np.isnan(raw).all(axis=1)
    if (not np.isfinite(raw[~unmeasured]).all() or (raw < 0).any()
            or (unmeasured.any() and not retain_unmeasured)):
        raise ValueError("Marker means must be finite and nonnegative; missing markers are not imputed")
    if unmeasured.all():
        raise ValueError("No cells have complete marker measurements")
    obs = frame.drop(columns=marker_columns).copy()
    # Some tiny InstanSeg objects have no intensity measurements at all. The
    # explicit archival mode preserves their IDs and NaNs for review; it does
    # not choose an exclusion policy or make this raw matrix fit for training.
    obs["marker_measurements_complete"] = ~unmeasured
    obs["Object ID"] = obs["Object ID"].astype(str)
    for source, target in [("Centroid X µm", "X_centroid"), ("Centroid Y µm", "Y_centroid")]:
        # The canonical H5AD stores micrometres. The viewer divides by its
        # pixel calibration when positioning these points over the TIFF.
        coordinates = pd.to_numeric(obs.pop(source), errors="raise")
        if not np.isfinite(coordinates).all() or (coordinates < 0).any():
            raise ValueError("Invalid cell centroid")
        obs[target] = coordinates
    for compartment in ("Cell", "Nucleus"):
        for measure, unit in [("Area", " µm^2"), ("Length", " µm"), ("Circularity", ""), ("Solidity", "")]:
            source = f"{compartment}: {measure}{unit}"
            if source in obs:
                obs.rename(columns={source: f"{compartment} {measure}"}, inplace=True)
    for source, target in DISTANCES.items():
        name = f"Signed distance to annotation {source} µm"
        if name in obs:
            obs.rename(columns={name: f"Dist to Closest {target}"}, inplace=True)
    obs["imageid"] = row["case_id"]
    obs["donor_id"] = row["case_id"]
    for key, value in metadata.items():
        obs[key] = value
    obs["Parent"] = obs["Parent"].fillna("").astype(str)
    obs.index = pd.Index(row["case_id"] + ":" + obs["Object ID"], name="cell_id")
    for name in obs.select_dtypes(include=["object"]).columns:
        obs[name] = obs[name].fillna("").astype(str)
    result = ad.AnnData(X=raw, obs=obs, var=pd.DataFrame(index=pd.Index(markers, name="marker")))
    result.layers["log"] = np.log1p(raw)
    result.uns["analysis_stage"] = "raw_measurements_unphenotyped"
    result.uns["quantification_pixel_size_um"] = QUANTIFICATION_PIXEL_UM
    result.uns["centroid_coordinate_units"] = "micrometres"
    result.uns["instanseg_model_version"] = "0.1.1"
    result.uns["missing_marker_policy"] = "retain_wholly_unmeasured" if retain_unmeasured else "reject"
    return result


def import_case(root, row, markers, metadata, cohort_path, key_path, retain_unmeasured=False):
    cid = row["case_id"]
    if row["inclusion_status"] != "selected":
        raise ValueError(f"Case {cid} is {row['inclusion_status']}; source identity must be resolved first")
    folder = root / "exports" / cid.zfill(4)
    source = folder / "cells.tsv"
    completed = root / "logs" / cid.zfill(4) / "export.done.json"
    if not completed.exists():
        raise ValueError(f"Case {cid}: QuPath export has not completed")
    export_manifest_path = folder / "export_manifest.json"
    if not export_manifest_path.exists():
        raise ValueError(f"Case {cid}: missing full-tissue export manifest")
    exported = json.loads(export_manifest_path.read_text())
    if case_id(exported.get("case_id")) != cid or exported.get("source_filename") != row["source_filename"]:
        raise ValueError("Export manifest donor/image identity mismatch")
    if exported.get("coordinate_units") != "micrometres":
        raise ValueError("Export must declare micrometre coordinates")
    if exported.get("numeric_encoding") != "java_double_to_string":
        raise ValueError("Full-precision QuPath numeric export is required before raw import")
    for key in ["pixel_size_um_x", "pixel_size_um_y"]:
        if not np.isclose(float(exported.get(key, np.nan)), QUANTIFICATION_PIXEL_UM, rtol=0, atol=1e-9):
            raise ValueError("Export calibration differs from the original analysis convention")
    cells_hash = sha256(source)
    if exported.get("cells_sha256") != cells_hash:
        raise ValueError("Cell export hash differs from its completed manifest")
    # Converter reports use the viewer basename; match the explicit source file.
    reports = []
    for path in (root / "reports/image_conversion").glob("*.json"):
        report = json.loads(path.read_text())
        if isinstance(report, dict) and Path(report.get("source", "")).name == row["source_filename"]:
            reports.append((path, report))
    if len(reports) != 1:
        raise ValueError(f"Expected one image conversion report for {cid}, found {len(reports)}")
    conversion_path, conversion = reports[0]
    if conversion.get("status") != "validated":
        raise ValueError(f"Image conversion has not validated for {cid}")
    output = root / "analysis/raw_cells" / f"{cid}.h5ad"
    report_path = output.with_suffix(".json")
    marker_signature = hashlib.sha256(json.dumps(markers).encode()).hexdigest()
    fingerprints = {"cells_sha256": cells_hash, "cohort_sha256": sha256(cohort_path),
                    "donor_key_sha256": sha256(key_path), "script_sha256": sha256(__file__),
                    "conversion_report_sha256": sha256(conversion_path),
                    "export_manifest_sha256": sha256(export_manifest_path),
                    "marker_order_sha256": marker_signature,
                    "retain_unmeasured": bool(retain_unmeasured)}
    if output.exists():
        if not report_path.exists():
            raise ValueError(f"Uncheckpointed output exists: {output}")
        previous = json.loads(report_path.read_text())
        if previous.get("inputs") != fingerprints or previous.get("output_sha256") != sha256(output):
            raise ValueError(f"Input or output changed for {cid}; use a new candidate directory")
        return previous
    frame = pd.read_csv(source, sep="\t", float_precision="round_trip")
    if len(frame) != exported.get("cell_count"):
        raise ValueError("Cell export row count does not match the full-tissue cell count")
    data = convert_frame(frame, markers, row, metadata, retain_unmeasured=retain_unmeasured)
    data.uns["source_files"] = fingerprints
    data.uns["source_numeric_encoding"] = exported["numeric_encoding"]
    data.uns["marker_storage_dtype"] = str(data.X.dtype)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".partial.h5ad")
    data.write_h5ad(temporary, compression="gzip")
    temporary.replace(output)
    report = {"case_id": cid, "status": "raw_import_complete", "cells": data.n_obs,
              "fully_measured_cells": int(data.obs["marker_measurements_complete"].sum()),
              "unmeasured_cell_objects": int((~data.obs["marker_measurements_complete"]).sum()),
              "unmeasured_object_ids": data.obs.loc[~data.obs["marker_measurements_complete"], "Object ID"].tolist(),
              "measurement_review_status": "pending" if not data.obs["marker_measurements_complete"].all() else "complete",
              "markers": data.var_names.tolist(), "inputs": fingerprints,
              "source_numeric_encoding": exported["numeric_encoding"],
              "marker_storage_dtype": str(data.X.dtype),
              "output": str(output), "output_sha256": sha256(output),
              "phenotyping_status": "pending_original_rule_recovery_or_reviewed_replacement"}
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--cohort", type=Path, default=COHORT)
    parser.add_argument("--donor-key", type=Path, default=Path("/data/panc_CODEX/results/CODEX_Pancreas_Donors.xlsx"))
    parser.add_argument("--baseline", type=Path, default=BASELINE)
    parser.add_argument("--cases", nargs="+", required=True)
    parser.add_argument("--retain-unmeasured", action="store_true",
                        help="Archive wholly unmeasured cell objects as flagged NaN rows for review; never impute or drop them")
    args = parser.parse_args()
    requested = {case_id(value) for value in args.cases}
    rows = [row for row in read_cohort(args.cohort) if row["case_id"] in requested]
    if {row["case_id"] for row in rows} != requested:
        parser.error("Requested case is absent from cohort")
    baseline = ad.read_h5ad(args.baseline, backed="r")
    markers = baseline.var_names.tolist()
    baseline.file.close()
    metadata = donor_metadata(args.donor_key, rows)
    for row in rows:
        report = import_case(args.root.resolve(), row, markers, metadata[row["case_id"]], args.cohort,
                             args.donor_key, retain_unmeasured=args.retain_unmeasured)
        print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
