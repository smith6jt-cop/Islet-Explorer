#!/usr/bin/env python3
"""Read-only comparison of retired-inclusive and retained-only live DPT inputs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "app/shiny_app/python/paga_pseudotime.py"


def audit(data_root):
    original = set(pd.read_csv(ROOT / "data/barlow_islet_features.csv").case_id.astype(str))
    h5ad = ad.read_h5ad(data_root / "adata_ins_root.h5ad")
    ids = pd.Index(h5ad.obs.combined_islet_id.astype(str))
    insulin = np.asarray(h5ad.X[:, list(h5ad.var_names).index("INS")]).ravel()
    output = {"data_root": str(data_root), "module_sha256": hashlib.sha256(MODULE.read_bytes()).hexdigest(),
              "features_sha256": hashlib.sha256((data_root / "barlow_islet_features.csv").read_bytes()).hexdigest(),
              "fits": {}}
    for policy in ("all_categories_before", "retained_categories_after"):
        spec = importlib.util.spec_from_file_location(policy, MODULE)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module._FEATURES_CSV = str(data_root / "barlow_islet_features.csv")
        if policy == "all_categories_before":
            # Reproduce the old feature selector; all other fitting steps match.
            module.composition_features = lambda cols: [c for c in cols if c.startswith(("core_", "peri_"))]
        for lo, hi in ((10, 1902), (25, 1902), (25, 200)):
            fit = module.compute_pseudotime(lo, hi)
            if "error" in fit:
                raise RuntimeError(fit["error"])
            ctx = module._context(lo, hi)
            selected = ctx["sub"].copy()
            indices = ids.get_indexer(selected.islet_id)
            if np.any(indices < 0):
                raise ValueError("Trajectory IDs absent from the marker artifact")
            selected["pt"] = fit["pseudotime"]
            selected["INS"] = insulin[indices]
            selected["added"] = ~selected.case_id.astype(str).isin(original)
            _, components = connected_components(ctx["A"].obsp["connectivities"])
            finite = selected[np.isfinite(selected.pt)]
            tail = finite[finite.pt >= finite.pt.quantile(.9)]
            between = sum(len(g) * (g.pt.mean() - finite.pt.mean()) ** 2
                          for _, g in finite.groupby("added"))
            total = np.sum((finite.pt - finite.pt.mean()) ** 2)
            deciles = pd.qcut(finite.pt, 10, labels=False, duplicates="drop")
            key = f"{policy}:{lo}-{hi}"
            output["fits"][key] = dict(
                selected=len(selected), finite=len(finite),
                components=np.bincount(components).tolist(),
                feature_columns=ctx["features"], root=selected.iloc[ctx["root"]].islet_id,
                insulin_spearman=float(spearmanr(finite.pt, finite.INS).statistic),
                beta_ratio_spearman=float(spearmanr(finite.pt, finite["core_Beta cell"]).statistic),
                cohort_eta_squared=float(between / total),
                high_decile_islets=len(tail), high_decile_added=int(tail.added.sum()),
                insulin_decile_means=finite.groupby(deciles).INS.mean().tolist())
            print(key, "n=", len(finite), "INS rho=", round(output["fits"][key]["insulin_spearman"], 4),
                  "cohort eta2=", round(between / total, 4), flush=True)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.data_root.resolve(strict=True))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.output.chmod(0o600)
