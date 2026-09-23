"""Compare browser trajectory exports to the frozen H5AD and composition data.

Run after deploy/browser/check-frozen-trajectory.mjs. All five exports are
required; missing files and numerical failures produce a nonzero exit status.
"""

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


MARKER_EXPORTS = {
    "rank-ins": "INS",
    "raw-ins": "INS",
    "rank-gcg": "GCG",
    "raw-gcg": "GCG",
}
TOLERANCE = 1e-10


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", default="/tmp/islet-restored-trajectory-20260923")
    parser.add_argument("--data", type=Path,
                        default=Path(__file__).resolve().parents[2] / "data")
    args = parser.parse_args()
    report = {"checks": []}

    def finish():
        report["passed"] = bool(report["checks"]) and all(
            check["passed"] for check in report["checks"])
        result = json.dumps(report, indent=2, default=int)
        Path(args.prefix + ".exports.json").write_text(result + "\n")
        print(result)
        return 0 if report["passed"] else 1

    files = {label: Path(f"{args.prefix}.{label}.csv")
             for label in [*MARKER_EXPORTS, "ratio-raw"]}
    missing = [str(path) for path in files.values() if not path.is_file()]
    report["checks"].append({"name": "all_five_browser_exports_present",
                             "passed": not missing, "missing": missing})
    if missing:
        return finish()

    try:
        h5ad = ad.read_h5ad(args.data / "adata_ins_root.h5ad")
        lookup = pd.Series(np.arange(h5ad.n_obs),
                           index=h5ad.obs["combined_islet_id"].astype(str))
        report["trajectory_h5ad_islets"] = h5ad.n_obs
        report["trajectory_h5ad_donors"] = h5ad.obs["imageid"].nunique()
        for label, marker in MARKER_EXPORTS.items():
            exported = pd.read_csv(files[label], comment="#")
            index = lookup.loc[exported.combined_islet_id].to_numpy()
            umap = h5ad.obsm["X_umap"][index]
            feature = h5ad[index, [marker]].X
            if hasattr(feature, "toarray"):
                feature = feature.toarray()
            feature = np.asarray(feature).reshape(-1)
            umap_error = float(np.max(np.abs(
                umap - exported[["umap_1", "umap_2"]].to_numpy())))
            feature_error = float(np.max(np.abs(
                feature - exported.feature_value.to_numpy())))
            report["checks"].append({
                "name": label + "_umap_and_marker_match_frozen_h5ad",
                "passed": umap_error < TOLERANCE and feature_error < TOLERANCE,
                "rows": len(exported), "max_umap_error": umap_error,
                "max_expression_error": feature_error,
            })

        exported = pd.read_csv(files["ratio-raw"], comment="#")
        features = pd.read_csv(args.data / "barlow_islet_features.csv").set_index("islet_id")
        types = exported.feature_name.unique()
        errors = []
        for cell_type in types:
            selected = exported.loc[exported.feature_name == cell_type]
            source = features.loc[selected.combined_islet_id]
            numerator = (source["n_core_" + cell_type].fillna(0)
                         + source["n_peri_" + cell_type].fillna(0))
            expected = np.log1p(numerator / source.n_endocrine_core.clip(lower=1))
            errors.extend(np.abs(selected.feature_value.to_numpy() - expected.to_numpy()))
        max_error = float(np.max(errors)) if errors else float("inf")
        report["checks"].append({
            "name": "ratio_matches_core_plus_peri_over_floored_endocrine_count",
            "passed": max_error < TOLERANCE, "max_error": max_error,
            "rows": len(exported), "cell_types": list(types),
        })
    except Exception as error:
        report["checks"].append({"name": "export_audit_completed", "passed": False,
                                 "error": str(error)})
    return finish()


if __name__ == "__main__":
    raise SystemExit(main())
