#!/usr/bin/env python3
"""H3 observed compartment densities using current app capillaries and corrected cells.

No pipeline changes, new cell labels, manual-review filters, or simulated nulls.
Primary width 10 um; 5 and 20 um are prespecified sensitivity measurements.
Local method precedents are recorded in the output's method_basis.json.
"""
from __future__ import annotations

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/h3-20261025-mpl")
import gzip
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import shapely
from shapely.affinity import scale
from shapely.geometry import shape
from statsmodels.stats.multitest import multipletests

import h1_features_20261022 as hf

REPO = Path(__file__).resolve().parents[2]
APP = REPO / "data/development/uniform-20donor-20261022-v10"
CORRECTED = REPO / "output/h1_edge_corrected_20261024"
OUT = REPO / "output/h1_h4_literature_methods_20261025/h3"
GROUPS = ["ND", "Aab+", "T1D"]
WIDTHS = [5, 10, 20]
EPS = 1e-9


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def polygon_parts(g):
    if g.geom_type == "Polygon":
        return [g] if g.area > EPS else []
    return [p for part in getattr(g, "geoms", []) for p in polygon_parts(part)]


def read_current_geometry(case):
    path = APP / "gson" / f"{case:04d}.geojson.gz"
    with gzip.open(path, "rt") as f:
        obj = json.load(f)
    cores, vessels = {}, []
    for feature in obj["features"]:
        p = feature.get("properties") or {}
        label = (p.get("classification") or {}).get("name")
        if label not in {"Capillary", "Islet"}:
            continue
        g = shape(feature["geometry"])
        g = shapely.make_valid(scale(g, xfact=hf.PIXEL_UM, yfact=hf.PIXEL_UM,
                                    origin=(0, 0)))
        if label == "Capillary":
            vessels.extend(polygon_parts(g))
        else:
            cores[p["name"]] = g
    return cores, vessels, path


def crossing_geometry(core, vessels, tree):
    crossings = []
    for k in tree.query(core, predicate="intersects"):
        vessel = vessels[k]
        if vessel.intersection(core).area > EPS and vessel.difference(core).area > EPS:
            cross = vessel.intersection(core.boundary)
            if not cross.is_empty:
                crossings.append(cross)
    return shapely.union_all(crossings), len(crossings)


def measure_compartments(core, crossing, width, xy):
    band = core.intersection(core.boundary.buffer(width))
    near_unclipped = core.intersection(crossing.buffer(width))
    # Polygonal circle approximation can leave tiny slivers outside an otherwise
    # identical-width boundary buffer. Clip to the shared comparison band.
    near = near_unclipped.intersection(band)
    reference = band.difference(near)
    if near.area <= EPS or reference.area <= EPS:
        return None
    pts = shapely.points(np.asarray(xy).reshape(-1, 2))
    is_near = shapely.covers(near, pts)
    is_reference = shapely.covers(reference, pts) & ~is_near
    return dict(near_area_um2=float(near.area), reference_area_um2=float(reference.area),
                band_area_um2=float(band.area), n_near=int(is_near.sum()),
                n_reference=int(is_reference.sum()), n_band=int((is_near | is_reference).sum()),
                area_partition_error_um2=float(abs(near.area + reference.area - band.area)),
                buffer_approximation_clip_um2=float(near_unclipped.area - near.area))


def summarize(donors):
    rows = []
    for width in WIDTHS:
        for group in ["Pooled", *GROUPS]:
            d = donors.loc[donors.width_um.eq(width)]
            if group != "Pooled":
                d = d.loc[d.donor_status.eq(group)]
            delta = d.near_density_mm2 - d.reference_density_mm2
            if len(d) and np.any(delta != 0):
                test = stats.wilcoxon(delta, zero_method="wilcox", alternative="two-sided",
                                      correction=False, method="auto")
                w, p = float(test.statistic), float(test.pvalue)
            else:
                w, p = 0.0, 1.0
            rows.append(dict(width_um=width, group=group, donors=len(d),
                objects=int(d.objects.sum()), n_near=int(d.n_near.sum()),
                n_reference=int(d.n_reference.sum()), donors_near_higher=int((delta > 0).sum()),
                donors_near_lower=int((delta < 0).sum()), donors_equal=int((delta == 0).sum()),
                median_near_density_mm2=float(d.near_density_mm2.median()),
                median_reference_density_mm2=float(d.reference_density_mm2.median()),
                median_difference_mm2=float(delta.median()),
                mean_difference_mm2=float(delta.mean()), wilcoxon_statistic=w, p_value=p))
    result = pd.DataFrame(rows)
    result["p_fdr_four_groups"] = result.groupby("width_um").p_value.transform(
        lambda p: multipletests(p, method="fdr_bh")[1])
    return result


def make_plot(donors, summary):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.8), sharey=True)
    colors = {"ND": "#3976a8", "Aab+": "#a87924", "T1D": "#bd5151"}
    for ax, group in zip(axes, GROUPS):
        d = donors.loc[donors.width_um.eq(10) & donors.donor_status.eq(group)]
        s = summary.loc[summary.width_um.eq(10) & summary.group.eq(group)].iloc[0]
        for r in d.itertuples():
            ax.plot([0, 1], [r.reference_density_mm2, r.near_density_mm2],
                    color=colors[group], alpha=.65, lw=1.2, marker="o", ms=5)
        ax.set_xticks([0, 1], ["Remaining\nboundary band", "Near capillary\ncrossings"])
        ax.set_xlim(-.25, 1.25)
        ax.set_title(f"{group}: {len(d)} donors\nPaired Wilcoxon p = {s.p_value:.3g}")
        ax.spines[["right", "top"]].set_visible(False)
        ax.grid(axis="y", color="#e4e4e4", linewidth=.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Immune cells per mm²")
    pooled = summary.loc[summary.width_um.eq(10) & summary.group.eq("Pooled")].iloc[0]
    fig.suptitle("H3: immune density at capillary crossings of the islet boundary\n"
                 f"All {int(pooled.donors)} donors: paired Wilcoxon p = {pooled.p_value:.4f}", fontsize=14)
    fig.text(.5, .025, "10 µm inside the core boundary; paired points are donor totals / compartment area.",
             ha="center", fontsize=10)
    fig.tight_layout(rect=[0, .06, 1, .91])
    fig.savefig(OUT / "h3_paired_donor_densities.png", dpi=200)
    plt.close(fig)


def process_case(item):
    case, sub = item
    objects, hashes = [], {}
    case = int(case)
    cores, vessels, gp = read_current_geometry(case)
    tree = shapely.STRtree(vessels)
    mp = CORRECTED / "object_memberships_by_donor" / f"{case}.parquet"
    cp = CORRECTED / "analysis_cells_by_donor" / f"{case}.parquet"
    for p in [gp, mp, cp]:
        hashes[str(p)] = sha(p)
    cells = pd.read_parquet(cp, columns=["object_id", "x_um", "y_um", "edge_analysis_eligible"])
    mem = pd.read_parquet(mp)
    immune = mem.loc[mem.cell_region.eq("core") & mem.phenotype.isin(hf.IMMUNE)]
    immune = immune.merge(cells, on="object_id", validate="many_to_one")
    immune = immune.loc[immune.edge_analysis_eligible]
    if immune.duplicated(["islet_id", "object_id"]).any():
        raise ValueError("Duplicate physical-cell membership")
    grouped = {key: z[["x_um", "y_um"]].to_numpy() for key, z in immune.groupby("islet_id")}
    for r in sub.itertuples():
        core = cores[r.base_islet_id]
        xy = grouped.get(r.islet_id, np.empty((0, 2)))
        if len(xy) != r.n_immune_core:
            raise ValueError(f"Core immune count mismatch: {r.islet_id}")
        crossing, ncross = crossing_geometry(core, vessels, tree)
        for width in WIDTHS:
            row = dict(case_id=case, donor_status=r.donor_status, islet_id=r.islet_id,
                       object_class=r.object_class, width_um=width, core_immune=len(xy),
                       crossing_components=ncross)
            if ncross == 0:
                objects.append(dict(row, status="no_capillary_boundary_crossing"))
                continue
            result = measure_compartments(core, crossing, width, xy)
            if result is None:
                objects.append(dict(row, status="no_reference_or_near_area"))
            else:
                objects.append(dict(row, status="included", **result))
    input_audit = dict(case_id=case, donor_status=sub.donor_status.iloc[0],
        primary_objects=len(sub), core_immune=int(sub.n_immune_core.sum()),
        current_capillary_components=len(vessels), app_geometry=str(gp))
    print(f"H3 donor {case}: {len(sub)} objects, {len(vessels)} current capillary components", flush=True)
    pd.DataFrame(objects).to_csv(OUT / f"donor_{case}_objects.csv", index=False)
    return objects, input_audit, hashes


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    primary_path = CORRECTED / "models/primary_analysis_objects.csv"
    primary = pd.read_csv(primary_path)
    if len(primary) != 2937 or primary.case_id.nunique() != 20 or 78 in primary.case_id.values:
        raise ValueError("Unexpected corrected H1 population")
    hashes = {str(primary_path): sha(primary_path)}
    objects, input_audit = [], []
    with ProcessPoolExecutor(max_workers=4) as executor:
        for records, audit, donor_hashes in executor.map(process_case, primary.groupby("case_id")):
            objects.extend(records)
            input_audit.append(audit)
            hashes.update(donor_hashes)
    obj = pd.DataFrame(objects)
    obj.to_csv(OUT / "h3_object_compartments.csv", index=False)
    included = obj.loc[obj.status.eq("included")]
    donors = included.groupby(["case_id", "donor_status", "width_um"]).agg(
        objects=("islet_id", "size"), core_immune=("core_immune", "sum"),
        n_near=("n_near", "sum"), n_reference=("n_reference", "sum"),
        near_area_um2=("near_area_um2", "sum"), reference_area_um2=("reference_area_um2", "sum")).reset_index()
    for compartment in ["near", "reference"]:
        donors[f"{compartment}_density_mm2"] = donors[f"n_{compartment}"] * 1e6 / donors[f"{compartment}_area_um2"]
    donors["difference_density_mm2"] = donors.near_density_mm2 - donors.reference_density_mm2
    donors.to_csv(OUT / "h3_donor_densities.csv", index=False)
    summary = summarize(donors)
    summary.to_csv(OUT / "h3_summary.csv", index=False)
    pd.DataFrame(input_audit).to_csv(OUT / "h3_input_audit.csv", index=False)
    obj.groupby(["donor_status", "width_um", "status"]).size().rename("objects").reset_index().to_csv(
        OUT / "h3_object_eligibility.csv", index=False)
    (OUT / "source_hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
    method_basis = dict(primary_width_um=10, sensitivity_widths_um=[5, 20],
        near="core intersect buffer(capillary component intersect core boundary, width)",
        reference="core intersect buffer(core boundary, width), minus near",
        crossing="Current app capillary polygon has positive area both inside and outside core",
        paired_unit="Donor; sum counts / sum areas within each compartment",
        inclusion="Corrected H1 primary objects with both nonzero comparison areas; zero-immune objects included",
        test="Two-sided Wilcoxon signed-rank, scipy method=auto; FDR over four group tests per width",
        source_method_inventory=str(OUT.parent / "method_inventory.json"),
        sources=[dict(pdf=str(REPO / "literature/Barlow2025.pdf"), pages=[10, 14], method="paired donor comparisons and connected donor points"),
                 dict(pdf=str(REPO / "literature/Richardson2023.pdf"), pages=[3], method="structure density normalized by ROI area"),
                 dict(pdf=str(REPO / "literature/Heide.pdf"), pages=[30], method="geometric intersections and positive-area fractions"),
                 dict(pdf=str(REPO / "literature/Bodenmiller2025.pdf"), pages=[17, 27], method="distance-defined spatial compartments")])
    (OUT / "method_basis.json").write_text(json.dumps(method_basis, indent=2) + "\n")
    make_plot(donors, summary)
    validation = dict(primary_objects=len(primary), donors=primary.case_id.nunique(),
        primary_core_immune=int(primary.n_immune_core.sum()), donor78_excluded=True,
        all_input_hashes_unchanged=all(sha(p) == h for p, h in hashes.items()),
        largest_area_partition_error_um2=float(included.area_partition_error_um2.max()),
        largest_buffer_approximation_clip_um2=float(included.buffer_approximation_clip_um2.max()),
        count_partition_valid=bool((included.n_near + included.n_reference == included.n_band).all()),
        counts_bounded_by_original_core=bool((included.n_band <= included.core_immune).all()))
    (OUT / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    run()
