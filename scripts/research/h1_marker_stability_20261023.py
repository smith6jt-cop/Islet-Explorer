#!/usr/bin/env python3
"""Post-fit point/whole-donor stability checks, motivated by secondary signals.

No new inferential model selection: all supported secondary slots receive the
same leave-one-donor and threshold checks. Bootstrap is limited to the three
small-object secondary signals with q<.05, and clearly labeled post-fit.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from h1_immune_markers_size_20261023 import marker_frame, fit_conditional
from h1_immune_types_size_20261023 import input_data, write_json

ROOT = Path("output/h1_immune_size_20261023")


def subset(base, table, marker, threshold, size):
    a = marker_frame(base, table, marker, threshold)
    a = a[(a.equivalent_diameter == size) & a.base_eligible & a.marker_complete & (a.type_nearby > 0)].copy()
    z = a.groupby("case_id").n_beta_core.sum()
    return a[~a.case_id.isin(z[z == 0].index)].copy()


def main():
    out = ROOT/"marker_models"
    estimates = pd.read_csv(out/"estimates.csv")
    e = estimates[(estimates.family == "conditional_broad") & estimates.rr.notna()]
    write_json(out/"postfit_stability_plan.json", {"reason": "Unexpected secondary conditional-on-broad signals; primary marker family had no discoveries", "checks": "LODO point-only and threshold0.05/0.20 point-only for every supported conditional slot; 999 whole-donor resamples for each secondary q<.05 slot", "bootstrap_seed": 20261023, "limits": "Post-fit diagnostics, percentile intervals are exploratory and not new adjusted tests"})
    base = input_data(Path("output/h1_20261022/h1_object_features.csv"))
    table = pd.read_parquet(ROOT/"object_marker_counts.parquet")
    table["case_id"] = table.case_id.astype(str)
    leave, thresholds, boots, summaries = [], [], [], []
    rng = np.random.default_rng(20261023)
    for er in e.itertuples():
        d = subset(base, table, er.marker, .1, er.size_bin)
        for donor in sorted(d.case_id.unique()):
            row = {"marker": er.marker, "size_bin": er.size_bin, "omitted_donor": donor, "rr": np.nan, "status": "FAILED"}
            try:
                point, _ = fit_conditional(d[d.case_id != donor], infer=False)
                row.update(point, status="FITTED_POINT_ONLY")
            except Exception as exc:
                row["reason"] = str(exc)
            leave.append(row)
        for threshold in [.05, .20]:
            sub = subset(base, table, er.marker, threshold, er.size_bin)
            row = {"marker": er.marker, "size_bin": er.size_bin, "threshold": threshold, "objects": len(sub), "donors": sub.case_id.nunique(), "informative_donors": int(sub.groupby("case_id").type_allocation.nunique().gt(1).sum()), "core_positive_objects": int((sub.type_core > 0).sum()), "rr": np.nan, "status": "UNSUPPORTED"}
            if row["objects"] >= 100 and row["informative_donors"] >= 6 and row["core_positive_objects"] >= 30:
                try:
                    point, _ = fit_conditional(sub, infer=False)
                    row.update(point, status="FITTED_POINT_ONLY")
                except Exception as exc:
                    row["reason"] = str(exc)
            thresholds.append(row)
        if er.q_bh < .05:
            donor_parts = {donor: part for donor, part in d.groupby("case_id")}
            donors = np.array(sorted(donor_parts))
            values = []
            for b in range(999):
                drawn = rng.choice(donors, len(donors), replace=True)
                sub = pd.concat([donor_parts[donor].assign(case_id=f"draw_{i:02d}") for i, donor in enumerate(drawn)], ignore_index=True)
                row = {"marker": er.marker, "size_bin": er.size_bin, "replicate": b, "rr": np.nan, "status": "FAILED"}
                try:
                    point, _ = fit_conditional(sub, infer=False)
                    row.update(point, status="FITTED_POINT_ONLY")
                    values.append(point["rr"])
                except Exception as exc:
                    row["reason"] = str(exc)
                boots.append(row)
            summaries.append({"marker": er.marker, "size_bin": er.size_bin, "attempts": 999, "successful": len(values), "rr_original": er.rr, "percentile_lower": float(np.quantile(values, .025)), "percentile_upper": float(np.quantile(values, .975))})
            print(json.dumps(summaries[-1]), flush=True)
    pd.DataFrame(leave).to_csv(out/"conditional_leave_one_donor_out.csv", index=False)
    pd.DataFrame(thresholds).to_csv(out/"conditional_threshold_point_checks.csv", index=False)
    pd.DataFrame(boots).to_csv(out/"conditional_donor_bootstrap.csv", index=False)
    pd.DataFrame(summaries).to_csv(out/"conditional_bootstrap_summary.csv", index=False)


if __name__ == "__main__":
    main()
