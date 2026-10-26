#!/usr/bin/env python3
"""Read-only, ID-linked source nuclear Ki67 audit for the accepted v10 cohort.

No threshold is trained and no Ki67-positive or proliferating cell is called.
Only selected columns are read from each source parquet; whole matrices are never
materialised. The source and app data are never changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import shapely

sys.path.insert(0, str(Path(__file__).resolve().parent))
from h1_features_20261022 import geometry, IMMUNE, sha256, finite_json

PHENOTYPES = ["CD8a Tcell", "CD4 Tcell", "T cell", "B cell", "Macrophage", "Immune"]
COMPARTMENTS = ["paired_core", "unpaired_core", "paired_peri_exclusive",
                "unpaired_peri_exclusive", "other_tissue", "outside_tissue"]
READOUTS = ["ki67_redsea_nuclear", "ki67_aligned_nuclear", "ki67_common_scaled",
            "ki67_donor_immune_midrank", "ki67_qupath_nuclear_mean",
            "ki67_log1p_nucleus_minus_cytoplasm"]


def record(path, columns=None):
    p = Path(path)
    return {"path": str(p.resolve()), "bytes": p.stat().st_size,
            "sha256": sha256(p), "columns_read": columns}


def stats(f, keys):
    r = dict(keys, n_cells=len(f))
    for col in READOUTS:
        v = f[col].to_numpy(float)
        v = v[np.isfinite(v)]
        r[col + "_n"] = len(v)
        r[col + "_missing"] = len(f) - len(v)
        for q, label in [(0.1, "q10"), (0.25, "q25"), (0.5, "median"), (0.75, "q75"), (0.9, "q90")]:
            r[col + "_" + label] = float(np.quantile(v, q)) if len(v) else np.nan
        r[col + "_mean"] = float(v.mean()) if len(v) else np.nan
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-root", type=Path, default=Path("data/development/uniform-20donor-20261015"))
    ap.add_argument("--app-root", type=Path, default=Path("data/development/uniform-20donor-20261022-v10"))
    ap.add_argument("--h1-root", type=Path, default=Path("output/h1_20261022"))
    ap.add_argument("--output-root", type=Path, default=Path("output/h1_immune_size_20261023/ki67"))
    args = ap.parse_args()
    source, app, h1, out = [p.resolve() for p in [args.source_root, args.app_root, args.h1_root, args.output_root]]
    if any(out == p or out.is_relative_to(p) for p in [source, app]):
        raise ValueError("Audit output must be outside source/app data")
    (out / "immune_cells_by_donor").mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(app / "donor_metadata.csv")
    manifest = json.loads((source / "project_manifest.json").read_text())
    entries = {int(e["case_id"]): e for e in manifest["entries"]}
    if len(meta) != 20 or set(meta.case_id) != set(entries) or 78 in set(meta.case_id):
        raise ValueError("Not the accepted complete 20-donor cohort")
    pairs = pd.read_csv(h1 / "h1_object_features.csv", usecols=["case_id", "islet_id"])
    scale = pd.read_parquet(source / "typing/scale.parquet").loc["Ki67"]
    shifts_all = json.loads((source / "typing/shifts.json").read_text())
    shifts = shifts_all["shifts"]["Ki67"]
    background = pd.read_csv(source / "typing/background_check.csv")
    background = background.loc[background.marker.eq("Ki67")].iloc[0].to_dict()
    sources = [record(p) for p in [source / "project_manifest.json", source / "typing/scale.parquet",
                                  source / "typing/shifts.json", source / "typing/background_check.csv",
                                  app / "donor_metadata.csv", h1 / "h1_object_features.csv",
                                  Path(__file__), Path(__file__).with_name("h1_features_20261022.py"),
                                  Path("scripts/uniform/typing/codex_typing.py"),
                                  Path("scripts/uniform/typing/align_scale.py"),
                                  Path("scripts/uniform/typing/refine.py")]]
    donor_rows, summary_rows, qc_rows, dfs = [], [], [], []
    started = time.monotonic()
    for md in meta.sort_values("case_id").to_dict("records"):
        case = int(md["case_id"])
        tag = entries[case]["image_name"]
        served_path = app / "cell_types" / f"{case}.parquet"
        served_cols = ["object_id", "phenotype", "cell_region", "islet_name", "x_um", "y_um"]
        served = pd.read_parquet(served_path, columns=served_cols)
        f = served.loc[served.phenotype.isin(IMMUNE)].copy()
        if served.object_id.duplicated().any():
            raise ValueError(f"Duplicate served ID {case}")
        f["case_id"], f["donor_status"] = case, entries[case]["donor_status"]
        f["image_name"], f["acquisition_group"] = tag, entries[case]["cohort"]
        al_path = source / "typing/aligned" / f"{tag}.parquet"
        al_cols = ["object_id", "Ki67", "qc_analysis_eligible", "qc_estimation_eligible", "qc_reasons"]
        al = pd.read_parquet(al_path, columns=al_cols).rename(columns={"Ki67": "ki67_aligned_nuclear"})
        re_path = source / "typing/redsea" / f"{tag}.parquet"
        re_cols = ["object_id", "Nucleus__Ki67", "Cytoplasm__Ki67", "Membrane__Ki67", "Cell__Ki67",
                   "Nucleus__n_px", "Cytoplasm__n_px", "qc_flags"]
        re = pd.read_parquet(re_path, columns=re_cols).rename(columns={
            "Nucleus__Ki67": "ki67_redsea_nuclear", "Cytoplasm__Ki67": "ki67_redsea_cytoplasm",
            "Membrane__Ki67": "ki67_redsea_membrane", "Cell__Ki67": "ki67_redsea_whole_cell"})
        raw_path = source / "cells" / f"{tag}.parquet"
        raw_cols = ["Object ID", "Nucleus: Ki67: Mean", "Cytoplasm: Ki67: Mean", "Cell: Ki67: Mean",
                    "Nucleus: Area µm^2", "Cell: Area µm^2"]
        raw = pd.read_parquet(raw_path, columns=raw_cols).rename(columns={
            "Object ID": "object_id", "Nucleus: Ki67: Mean": "ki67_qupath_nuclear_mean",
            "Cytoplasm: Ki67: Mean": "ki67_qupath_cytoplasm_mean", "Cell: Ki67: Mean": "ki67_qupath_whole_cell_mean",
            "Nucleus: Area µm^2": "nucleus_area_um2", "Cell: Area µm^2": "cell_area_um2"})
        # Identity joins are strict, never row-order joins. Indicators establish coverage.
        for part, label in [(al, "aligned"), (re, "redsea"), (raw, "qupath")]:
            f = f.merge(part, on="object_id", how="left", validate="one_to_one", indicator="_link")
            if not f._link.eq("both").all():
                raise ValueError(f"Missing {label} IDs in {case}")
            f = f.drop(columns="_link")
        if not f.qc_analysis_eligible.all():
            raise ValueError(f"Served immune cell analysis-ineligible in {case}")
        factor = shifts["scale"].get(tag)
        if not isinstance(factor, (int, float)) or factor <= 0:
            raise ValueError(f"No positive alignment factor for {tag}")
        expected = (f.ki67_redsea_nuclear.to_numpy(float) * factor).astype(np.float32)
        actual = f.ki67_aligned_nuclear.to_numpy(np.float32)
        if not np.array_equal(expected, actual, equal_nan=True):
            raise ValueError(f"Ki67 source alignment mismatch in {case}")
        f["ki67_common_scaled"] = np.clip((actual - scale["floor"]) / scale["span"], 0, 1)
        # Use un-clipped source nuclear values, so floor/top censoring cannot create artificial ties.
        valid = np.isfinite(f.ki67_redsea_nuclear)
        nref = int(valid.sum())
        f["ki67_donor_immune_midrank"] = np.nan
        f.loc[valid, "ki67_donor_immune_midrank"] = (f.loc[valid, "ki67_redsea_nuclear"].rank(method="average") - 0.5) / nref
        f["ki67_rank_reference_n"] = nref
        f["ki67_log1p_nucleus_minus_cytoplasm"] = np.log1p(f.ki67_redsea_nuclear) - np.log1p(f.ki67_redsea_cytoplasm)
        cores, _, tissue, _, geom_path = geometry(app, case)
        f["in_tissue_geometry"] = shapely.covers(tissue, shapely.points(f[["x_um", "y_um"]].to_numpy(float)))
        paired_names = set(pairs.loc[pairs.case_id.eq(case), "islet_id"].str.removeprefix(f"{case}_"))
        f["paired_core"] = f.cell_region.eq("core") & f.islet_name.isin(paired_names)
        peri_path = h1 / "peri_cells_by_donor" / f"{case}.parquet"
        peri_cols = ["object_id", "membership_multiplicity", "exclusive_owner", "exclusive_owner_in_paired_table"]
        peri = pd.read_parquet(peri_path, columns=peri_cols)
        f = f.merge(peri, on="object_id", how="left", validate="one_to_one")
        f["membership_multiplicity"] = f.membership_multiplicity.fillna(0).astype(int)
        f["exclusive_owner"] = f.exclusive_owner.fillna("")
        f["exclusive_owner_in_paired_table"] = f.exclusive_owner_in_paired_table.fillna(False).astype(bool)
        f["audit_compartment"] = np.select([
            ~f.in_tissue_geometry, f.paired_core, f.cell_region.eq("core"),
            f.membership_multiplicity.gt(0) & f.exclusive_owner_in_paired_table,
            f.membership_multiplicity.gt(0)],
            ["outside_tissue", "paired_core", "unpaired_core", "paired_peri_exclusive", "unpaired_peri_exclusive"],
            default="other_tissue")
        f.to_parquet(out / "immune_cells_by_donor" / f"{case}.parquet", index=False)
        dfs.append(f)
        donor_rows.append(stats(f, {"case_id": case, "donor_status": entries[case]["donor_status"],
                                    "acquisition_group": entries[case]["cohort"], "population": "all_served_curated_immune"}))
        donor_rows.append(stats(f.loc[f.in_tissue_geometry], {"case_id": case, "donor_status": entries[case]["donor_status"],
                                                           "acquisition_group": entries[case]["cohort"], "population": "tissue_curated_immune"}))
        for phenotype in ["All curated immune"] + PHENOTYPES:
            for comp in ["all_served", "in_tissue"] + COMPARTMENTS:
                keep = np.ones(len(f), bool) if phenotype == "All curated immune" else f.phenotype.eq(phenotype).to_numpy()
                if comp == "in_tissue":
                    keep &= f.in_tissue_geometry.to_numpy()
                elif comp != "all_served":
                    keep &= f.audit_compartment.eq(comp).to_numpy()
                summary_rows.append(stats(f.loc[keep], {"case_id": case, "donor_status": entries[case]["donor_status"],
                                                      "phenotype": phenotype, "audit_compartment": comp}))
        for pop, keep in [("all_curated_immune", np.ones(len(f), bool)),
                          ("upper_rank_decile_descriptive", f.ki67_donor_immune_midrank.ge(0.9).to_numpy()),
                          ("lower_rank_half_descriptive", f.ki67_donor_immune_midrank.le(0.5).to_numpy())]:
            z = f.loc[keep]
            q = stats(z, {"case_id": case, "donor_status": entries[case]["donor_status"], "population": pop})
            q.update(nucleus_signal_gt_cytoplasm_fraction=float(z.ki67_redsea_nuclear.gt(z.ki67_redsea_cytoplasm).mean()) if len(z) else np.nan,
                     nuclear_cytoplasm_spearman=z.ki67_redsea_nuclear.corr(z.ki67_redsea_cytoplasm, method="spearman") if len(z)>2 else np.nan,
                     qupath_redsea_nucleus_spearman=z.ki67_qupath_nuclear_mean.corr(z.ki67_redsea_nuclear, method="spearman") if len(z)>2 else np.nan,
                     n_nucleus_px_lt_10=int(z["Nucleus__n_px"].lt(10).sum()),
                     n_cytoplasm_px_zero=int(z["Cytoplasm__n_px"].eq(0).sum()),
                     n_common_scale_zero=int(z.ki67_common_scaled.eq(0).sum()),
                     n_common_scale_one=int(z.ki67_common_scaled.eq(1).sum()),
                     alignment_factor=float(factor))
            qc_rows.append(q)
        for p, cols in [(served_path, served_cols), (al_path, al_cols), (re_path, re_cols),
                        (raw_path, raw_cols), (peri_path, peri_cols), (geom_path, None)]:
            sources.append(record(p, cols))
        print(f"{case}: {len(f)} curated immune; {nref} finite source nuclear Ki67; aligned identity verified", flush=True)
    all_cells = pd.concat(dfs, ignore_index=True)
    all_cells.to_parquet(out / "ki67_immune_cells.parquet", index=False)
    pd.DataFrame(donor_rows).to_csv(out / "ki67_donor_summary.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(out / "ki67_phenotype_compartment_summary.csv", index=False)
    pd.DataFrame(qc_rows).to_csv(out / "ki67_localization_qc.csv", index=False)
    pd.DataFrame({"column": all_cells.columns, "dtype": all_cells.dtypes.astype(str).tolist()}).to_csv(out / "ki67_cell_schema.csv", index=False)
    audit = {
        "audit_date": "2026-10-23", "source_root": str(source), "app_root": str(app),
        "accepted_cohort_cases": sorted(meta.case_id.astype(int).tolist()), "excluded_case_78_present": False,
        "n_curated_immune": len(all_cells), "n_in_tissue": int(all_cells.in_tissue_geometry.sum()),
        "n_finite_nuclear_ki67": int(np.isfinite(all_cells.ki67_redsea_nuclear).sum()),
        "n_source_join_missing": 0, "n_alignment_identity_mismatch": 0,
        "nucleus_compartment": "Nucleus; raw source image units; no spillover compensation for Ki67",
        "source_common_scale": scale.to_dict(), "alignment": shifts,
        "existing_background_check": background,
        "rank_definition": "Within each donor, all current served curated immune cells with finite REDSEA nuclear Ki67 define the reference: (average ascending rank - 0.5)/n. Ranks computed before clipping. Larger=relatively greater nuclear Ki67 signal within that donor. Not a cross-donor absolute proliferation scale.",
        "summary_definition": "Every donor × six phenotypes plus all immune × all served/in tissue and six exclusive audit compartments, including n=0 cells. Empty measurements stay missing. Peri exclusivity uses the verified prior H1 nearest-boundary allocation over all Islets/Clusters. The cell table can be joined by case_id/object_id to overlapping object memberships for primary H1 analysis.",
        "assay_status": "Continuous exploratory nuclear-signal endpoint available on all 20 donors; no validated Ki67-positive threshold or proliferation label.",
        "limitations": [
            "Numerical alignment and a pooled floor/top scale are available; these are technical transformations, not calibrated molecule abundance or validated cross-donor proliferation comparisons.",
            "Ki67 is absent from the current common served marker panel. This readout is a source-level extension joined to current IDs, without modifying app data.",
            "Current curated immune IDs condition on lineage-marker gating and expression/neighbour refinement. Ki67 is among REFINE_MARKERS; outcome-dependent selection cannot be excluded. Require an independent Ki67-omitted typing sensitivity before mechanistic inference.",
            "Existing background_check selects ECAD/panCK-high epithelial cells as nominal negatives. Epithelial cells can proliferate; that check cannot establish Ki67 specificity, sensitivity, or a biological negative distribution.",
            "Nuclear/cytoplasmic comparisons and upper-rank-decile summaries describe localization only. They cannot establish a validated positive class; a donor rank always has a high tail even if every cell is background.",
            "Raw QuPath nuclear means and raster-recomputed REDSEA nuclear means differ in mask rasterization and are not interchangeable. The alignment identity is explicitly checked against the latter.",
            "A genuine proliferation claim needs direct nuclear localization review with DAPI and Ki67, positive/negative control evidence and a prespecified validated threshold; no such assay validation is established by this audit.",
            "Donors, rather than cells or objects, are the independent units. Within-donor ranks remove donor intensity shifts and cannot measure global donor-wide elevation. Small phenotype/compartment groups are descriptive and must retain their sample sizes.",
            "Current analysis eligibility excludes absent nuclei; other geometry QC is flagged but not automatically excluded. Localization results should be checked against nucleus pixel count and QC flags."
        ],
        "elapsed_seconds": round(time.monotonic()-started, 2), "source_provenance": sources,
        "outputs": {p.name: {"sha256": sha256(p), "bytes": p.stat().st_size} for p in out.glob("*.parquet")}
    }
    (out / "audit.json").write_text(json.dumps(finite_json(audit), indent=2, allow_nan=False) + "\n")
    # A concise human-readable audit complements the complete machine-readable provenance.
    qc = pd.DataFrame(qc_rows)
    low, high = qc.loc[qc.population.eq("lower_rank_half_descriptive")], qc.loc[qc.population.eq("upper_rank_decile_descriptive")]
    text = f"""# Source Ki67 audit — 23 October 2026

All 20 accepted donors have ID-linked source nuclear Ki67. The audit joins {len(all_cells):,} current curated immune cells ({int(all_cells.in_tissue_geometry.sum()):,} within tissue geometry), with no missing source IDs, {int(all_cells.ki67_redsea_nuclear.isna().sum()):,} missing nuclear Ki67 measurements, and exact float32 agreement between aligned Ki67 and source nuclear Ki67 × the recorded donor factor. Donor 78 is excluded.

Source `typing/aligned` Ki67 is the nuclear channel, prior to common scaling. The source `typing/redsea` nuclear channel passes without spillover compensation; raw QuPath nuclear means are kept separately because rasterized masks differ. The common floor is {scale['floor']:.6f}, top {scale['top']:.6f}. Alignment is `{shifts['provenance']}` across {shifts['n_donors']} donors; mode SD falls from {shifts['mode_sd_before']:.6f} to {shifts['mode_sd_after']:.6f} log. This establishes numerical availability and technical concordance, not biological calibration.

`ki67_donor_immune_midrank` = (average ascending rank − 0.5)/n among every finite current curated immune cell of the same donor, before clipping. This continuous endpoint can compare immune localization within donors. It cannot test whether all immune cells in one donor have globally increased Ki67. No threshold-positive or proliferating-cell calls are provided.

Across donors, the median nuclear-minus-cytoplasmic log1p signal ranges {low['ki67_log1p_nucleus_minus_cytoplasm_median'].min():.3f}–{low['ki67_log1p_nucleus_minus_cytoplasm_median'].max():.3f} in the lower rank half and {high['ki67_log1p_nucleus_minus_cytoplasm_median'].min():.3f}–{high['ki67_log1p_nucleus_minus_cytoplasm_median'].max():.3f} in the upper rank decile. These descriptive localization comparisons require image/control review before any positive class can be inferred.

`ki67_immune_cells.parquet` gives persistent `case_id`/`object_id` joins, phenotype, source region/name, coordinates, tissue membership, exclusive peri owner, source signals, common-scale signal, rank, nucleus geometry, and QC. `ki67_phenotype_compartment_summary.csv` enumerates each donor × six immune phenotypes plus all immune × all served/in tissue/exclusive compartments, retaining empty groups and missingness. Primary overlapping peri analyses should join this table's cell IDs to the full object membership table; the exclusive peri summaries here are a sensitivity.

Material limits:

- Current immune selection uses lineage gates and expression/neighbour refinement, whose marker list includes Ki67. This is a conditional analysis of current curated immune cells. Independent typing with Ki67 omitted is needed to address selection dependence.
- The existing background check reports Ki67 SD 0.087 log after alignment, but its ECAD/panCK-high epithelial pool is not a validated biological negative for a proliferation marker.
- Donor ranks have a high tail even under background-only signal. Neither rank tails, localization ratios nor the pooled scale validate proliferation. Direct DAPI/Ki67 image review, true controls and validated thresholds remain necessary.
- Absent nuclei are excluded upstream. Other geometry issues are flagged; pixel-count and QC sensitivity remain necessary. Small groups stay descriptive; donors are the inferential units.

The reproducible script is `scripts/research/h1_ki67_audit_20261023.py`. Full source paths, read-column lists, file SHA256 hashes, alignment values, definitions and limitations are in `audit.json`. Source and served data were read only.
"""
    (out / "audit.md").write_text(text)
    print(f"Complete: {len(all_cells)} immune cells; {out}", flush=True)


if __name__ == "__main__":
    main()
