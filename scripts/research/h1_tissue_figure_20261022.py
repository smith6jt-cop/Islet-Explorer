#!/usr/bin/env python3
"""Auditable central matched-pair tissue maps for H1; no outcome-tail selection.

Run after the matching results are finalized. It selects the donor with most
matched pairs (at least5), then the pair closest to that donor's median Beta
density difference; ties use match distance and IDs. A second row applies the
same rule among ND donors. Source cells use persistent identities, not row joins.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "islet-explorer-h1-matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
import shapely
from shapely import STRtree
from shapely.affinity import translate
from shapely.geometry.polygon import orient

import h1_features_20261022 as hf

PALETTE = {"Beta": "#0072B2", "Other endocrine": "#E69F00",
           "Immune": "#CC3311", "Residual Exocrine": "#7D8790"}


def central_pair(pairs, status=None):
    d = pairs if status is None else pairs.loc[pairs.donor_status.eq(status)]
    donors = d.groupby("case_id").size().rename("n_pairs").reset_index()
    donors = donors.loc[donors.n_pairs.ge(5)].sort_values(["n_pairs", "case_id"], ascending=[False, True])
    if donors.empty:
        raise ValueError(f"No donor with at least5 pairs in {status or 'all donors'}")
    case = int(donors.iloc[0].case_id)
    candidates = d.loc[d.case_id.eq(case)].copy()
    outcome = "beta_density_difference_cells_per_1000um2"
    median = float(candidates[outcome].median())
    candidates["distance_from_donor_median"] = (candidates[outcome] - median).abs()
    candidates = candidates.sort_values(["distance_from_donor_median", "match_distance_squared",
                                         "treated_islet_id", "control_islet_id"], kind="stable")
    chosen = candidates.iloc[0].to_dict()
    chosen["selection_donor_pair_count"] = len(candidates)
    chosen["selection_donor_median_beta_density_difference"] = median
    chosen["selection_scope"] = status or "all donors"
    chosen["number_of_pairs_tied_at_median_distance"] = int(candidates.distance_from_donor_median.eq(
        chosen["distance_from_donor_median"]).sum())
    return chosen


def polygons(g):
    if g.geom_type == "Polygon":
        yield g
    elif hasattr(g, "geoms"):
        for part in g.geoms:
            yield from polygons(part)


def polygon_path(poly):
    poly = orient(poly, sign=1)
    verts, codes = [], []
    for ring in [poly.exterior, *poly.interiors]:
        xy = np.asarray(ring.coords)
        verts.extend(xy)
        codes.extend([MplPath.MOVETO] + [MplPath.LINETO] * (len(xy) - 2) + [MplPath.CLOSEPOLY])
    return MplPath(np.asarray(verts), codes)


def draw_shape(ax, geom, facecolor, edgecolor, linewidth, linestyle="solid", zorder=1):
    for poly in polygons(geom):
        ax.add_patch(PathPatch(polygon_path(poly), facecolor=facecolor, edgecolor=edgecolor,
                              linewidth=linewidth, linestyle=linestyle, zorder=zorder))


def cell_group(phenotype):
    if phenotype == "Beta cell":
        return "Beta"
    if phenotype in hf.IMMUNE:
        return "Immune"
    if phenotype in hf.ENDOCRINE:
        return "Other endocrine"
    if phenotype == "Exocrine":
        return "Residual Exocrine"
    raise ValueError(f"Unexpected phenotype {phenotype}")


def donor_data(root, case):
    cores, expanded, _, repairs, geometry_path = hf.geometry(root, case)
    cell_path = root / "cell_types" / f"{case}.parquet"
    cells = pd.read_parquet(cell_path)
    if cells.object_id.duplicated().any():
        raise ValueError("Persistent cell identities duplicated")
    points = shapely.points(cells[["x_um", "y_um"]].to_numpy())
    return dict(cores=cores, expanded=expanded, cells=cells, points=points, tree=STRtree(points),
                repairs=repairs, geometry_path=geometry_path, cell_path=cell_path)


def object_data(bundle, feature):
    name = feature.base_islet_id
    cls, core = bundle["cores"][name]
    if cls != feature.object_class:
        raise ValueError("Object class differs between source and analysis")
    expanded = bundle["expanded"][name]
    band = expanded.difference(core)
    cells = bundle["cells"]
    ci = np.flatnonzero(cells.cell_region.eq("core") & cells.islet_name.eq(name))
    pi = bundle["tree"].query(band, predicate="covers")
    pi = pi[~cells.iloc[pi].cell_region.eq("core").to_numpy()]
    core_beta = int(cells.iloc[ci].phenotype.eq("Beta cell").sum())
    core_immune = int(cells.iloc[ci].phenotype.isin(hf.IMMUNE).sum())
    peri_immune = int(cells.iloc[pi].phenotype.isin(hf.IMMUNE).sum())
    if (len(ci), len(pi), core_beta, core_immune, peri_immune) != (
            feature.n_cells_core, feature.n_cells_peri, feature.n_beta_core,
            feature.n_immune_core, feature.n_immune_peri):
        raise ValueError(f"Figure source counts disagree for {feature.islet_id}")
    ids = np.r_[ci, pi]
    if len(np.unique(cells.iloc[ids].object_id)) != len(ids):
        raise ValueError("Core+peri set repeats a persistent cell")
    bounds = expanded.bounds
    cx, cy = core.centroid.x, core.centroid.y
    reach = max(abs(bounds[0] - cx), abs(bounds[2] - cx), abs(bounds[1] - cy), abs(bounds[3] - cy))
    required_width = max(100, math.ceil((2 * reach + 30) / 25) * 25)
    return dict(feature=feature, core=core, band=band, expanded=expanded, ids=ids,
                ci=ci, pi=pi, center=(cx, cy), required_width_um=required_width)


def panel(ax, bundle, obj, width, letter, core_present):
    f, cells = obj["feature"], bundle["cells"]
    cx, cy = obj["center"]
    coords = cells[["x_um", "y_um"]].to_numpy() - [cx, cy]
    visible = (np.abs(coords[:, 0]) <= width / 2) & (np.abs(coords[:, 1]) <= width / 2)
    outside = np.flatnonzero(visible & ~np.isin(np.arange(len(cells)), obj["ids"]))
    ax.scatter(coords[outside, 0], coords[outside, 1], s=8, color="#C8CDD1", alpha=.40,
               linewidths=0, zorder=2, rasterized=False)
    band = translate(obj["band"], -cx, -cy)
    core = translate(obj["core"], -cx, -cy)
    draw_shape(ax, band, "#F1F4F6", "none", 0)
    draw_shape(ax, translate(obj["expanded"], -cx, -cy), "none", "#6B7780", .9, "dashed", 3)
    draw_shape(ax, core, "white", "#1E2B35", 1.35, zorder=4)
    chosen = cells.iloc[obj["ids"]]
    groups = chosen.phenotype.map(cell_group).to_numpy()
    for group, color in PALETTE.items():
        ids = obj["ids"][groups == group]
        ax.scatter(coords[ids, 0], coords[ids, 1], s=24, color=color, edgecolor="white", linewidth=.35,
                   zorder=6 if group == "Immune" else 5, alpha=.95)
    ax.set_xlim(-width / 2, width / 2)
    ax.set_ylim(width / 2, -width / 2)
    ax.set_aspect("equal")
    ax.set_xlabel("µm from core centroid", fontsize=9)
    ax.set_ylabel("µm from core centroid", fontsize=9)
    ax.tick_params(labelsize=8, length=3)
    for spine in ax.spines.values():
        spine.set_color("#BCC5CB")
    state = "Core immune present" if core_present else "Core immune absent"
    ax.set_title(f"{letter}   {state}\n{f.islet_id}", loc="left", fontsize=11, pad=9)
    beta_density = f.n_beta_core / f.core_area_um2 * 1000
    core_unit, peri_unit = ("cell" if f.n_cells_core == 1 else "cells"), ("cell" if f.n_cells_peri == 1 else "cells")
    text = (f"Core: {f.n_cells_core} {core_unit}  |  β {f.n_beta_core}  |  immune {f.n_immune_core}\n"
            f"Peri: {f.n_cells_peri} {peri_unit}  |  immune {f.n_immune_peri}\n"
            f"β density: {beta_density:.2f} cells / 1,000 µm²")
    ax.text(0, -.24, text, transform=ax.transAxes, va="top", fontsize=9, color="#24343F", linespacing=1.45)
    bar = 25 if width < 200 else (50 if width < 500 else 100)
    x, y = -width / 2 + width * .07, width / 2 - width * .09
    ax.plot([x, x + bar], [y, y], color="#172731", lw=2.5, solid_capstyle="butt", zorder=9)
    ax.text(x + bar / 2, y - width * .025, f"{bar} µm", fontsize=8, ha="center", va="bottom", zorder=9)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--matches", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    root, out = args.data_root.resolve(strict=True), args.output_root.resolve()
    out.mkdir(parents=True, exist_ok=True)
    features = pd.read_csv(args.features).set_index("islet_id", drop=False)
    pairs = pd.read_csv(args.matches)
    if 78 in set(pairs.case_id) or pairs.empty or pairs.pair_id.duplicated().any():
        raise ValueError("Invalid matched-pair cohort")
    protocol_path = args.matches.parent / "protocol.json"
    protocol = json.loads(protocol_path.read_text())
    feature_hash = hf.sha256(args.features)
    if protocol["feature_sha256"] != feature_hash:
        raise ValueError("Matched results reference a different feature table")
    selected = [central_pair(pairs), central_pair(pairs, "ND")]
    if selected[0]["case_id"] == selected[1]["case_id"]:
        raise ValueError("Overall and ND examples coincide; clarify second-row protocol")
    bundles, objects = {}, []
    for pair in selected:
        case = int(pair["case_id"])
        bundles[case] = donor_data(root, case)
        row = []
        for key in ["treated_islet_id", "control_islet_id"]:
            f = features.loc[pair[key]]
            if f.case_id != case or f.object_class != pair["object_class"]:
                raise ValueError("Matched identities differ from source feature table")
            row.append(object_data(bundles[case], f))
        objects.append(row)
    widths = [max(x["required_width_um"] for x in row) for row in objects]
    if max(widths) / min(widths) <= 1.8:
        widths = [max(widths)] * 2
        scale_note = f"Every panel has the same {widths[0]} µm physical width."
    else:
        scale_note = f"Matched-pair widths: top {widths[0]} µm; bottom {widths[1]} µm. Scale is shared within each pair."
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "savefig.facecolor": "white"})
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 11.5))
    fig.subplots_adjust(left=.08, right=.98, top=.79, bottom=.18, hspace=.92, wspace=.25)
    fig.suptitle("Central matched tissue examples", x=.08, y=.982, ha="left", fontsize=20, weight="bold", color="#183340")
    fig.text(.08, .948, "Fixed donor-and-median selection rule • cell centroids with actual core and 20 µm peri boundaries",
             fontsize=10.3, color="#4D5F6A")
    legend = [Line2D([0], [0], marker="o", color="none", markerfacecolor=c, markeredgecolor="white",
                      markersize=7, label=l) for l, c in PALETTE.items()]
    fig.legend(handles=legend, loc="upper left", bbox_to_anchor=(.073, .925), ncol=4, frameon=False, fontsize=10,
               columnspacing=2.0, handletextpad=.4)
    letters = [["A", "B"], ["C", "D"]]
    for r, pair in enumerate(selected):
        case = int(pair["case_id"])
        for c in range(2):
            panel(axes[r, c], bundles[case], objects[r][c], widths[r], letters[r][c], c == 0)
        label = (f"Donor {case} ({pair['donor_status']}) • {pair['object_class']} • "
                 f"central pair among {pair['selection_donor_pair_count']} matches")
        axes[r, 0].text(0, 1.36, label, transform=axes[r, 0].transAxes, fontsize=11.4, weight="bold", color="#183340")
    fig.text(.08, .065, "Solid: core boundary    Dashed: expanded 20 µm edge    Faint gray: nearby cells outside the focal core/peri region",
             fontsize=9, color="#52636E")
    fig.text(.08, .043, scale_note, fontsize=9, color="#52636E")
    fig.text(.08, .020, "These descriptive examples do not estimate the pooled association or establish a causal effect.",
             fontsize=9, color="#52636E")
    png, svg = out / "figure5_tissue_example.png", out / "figure5_tissue_example.svg"
    fig.savefig(png, dpi=250)
    fig.savefig(svg)
    plt.close(fig)
    metadata = dict(selection_rule={"donor": "largest matched-pair count among donors with>=5; tie lowest donor ID",
                                   "pair": "closest to donor median Beta-density difference; tie lowest squared match distance then treated/control IDs",
                                   "rows": ["overall donor", "ND donor"],
                                   "scale": "one width across all panels if geometry-derived pair widths differ<=1.8-fold, otherwise one width per matched pair"},
                    selected_pairs=selected, widths_um=widths, scale_note=scale_note,
                    source_sha256={str(args.features.resolve()): feature_hash, str(args.matches.resolve()): hf.sha256(args.matches),
                                   str(protocol_path.resolve()): hf.sha256(protocol_path),
                                   **{str(v[k]): hf.sha256(v[k]) for v in bundles.values() for k in ["geometry_path", "cell_path"]}},
                    plotted_objects=[dict(islet_id=o["feature"].islet_id,
                                          centroid_x_um=o["center"][0], centroid_y_um=o["center"][1],
                                          n_core=len(o["ci"]), n_peri=len(o["pi"]),
                                          n_beta_core=int(o["feature"].n_beta_core),
                                          n_immune_core=int(o["feature"].n_immune_core),
                                          n_immune_peri=int(o["feature"].n_immune_peri),
                                          core_area_um2=float(o["feature"].core_area_um2)) for row in objects for o in row],
                    script_sha256=hf.sha256(Path(__file__)),
                    helper_script_sha256=hf.sha256(Path(hf.__file__)),
                    software={"matplotlib": matplotlib.__version__, "shapely": shapely.__version__,
                              "numpy": np.__version__, "pandas": pd.__version__},
                    output_sha256={str(p): hf.sha256(p) for p in [png, svg]})
    hf.write_json(out / "figure5_selection.json", metadata)
    print(json.dumps(hf.finite_json({"selected_pair_ids": [p["pair_id"] for p in selected],
                                    "widths_um": widths, "png": str(png), "svg": str(svg)}), indent=2))


if __name__ == "__main__":
    main()
