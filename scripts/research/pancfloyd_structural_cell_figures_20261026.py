#!/usr/bin/env python3
"""Restore the reference structural-cell panels from the existing rule's recovery.

The recovery is analysis-only, with unchanged accepted cell universe and edges.
This plotting script neither assigns cells nor changes the current served data.
"""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/pancfloyd_20donor_20261026/structural_cells"
SOURCE = OUT / "object_structure_counts.csv"
FIG = OUT / "figures"
DATA = OUT / "data"
GROUPS = ["ND", "Aab+", "T1D"]
COLORS = dict(zip(GROUPS, ["#3366CC", "#F07618", "#008D28"]))
MIDS = np.arange(25, 350, 50)
PANELS = [
    ("slide03_image01", "vessel_core_percent", "Blood vessel", "core", "percent", "Blood-vessel cells · core"),
    ("slide03_image02", "vessel_peri_percent", "Blood vessel", "peri", "percent", "Blood-vessel cells · peri"),
    ("slide03_image03", "vessel_core_density", "Blood vessel", "core", "density", "Blood-vessel cell density · core"),
    ("slide03_image04", "vessel_peri_density", "Blood vessel", "peri", "density", "Blood-vessel cell density · peri"),
    ("slide04_image04", "nerve_core_percent", "Nerve", "core", "percent", "Neural cells · core"),
    ("slide04_image05", "nerve_peri_percent", "Nerve", "peri", "percent", "Neural cells · peri"),
    ("slide09_image01", "lymphatic_core_percent", "Lymphatic", "core", "percent", "Lymphatic cells · core"),
    ("slide09_image02", "lymphatic_peri_percent", "Lymphatic", "peri", "percent", "Lymphatic cells · peri"),
]


def digest(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(8*1024*1024), b""):
            h.update(b)
    return h.hexdigest()


def group_counts(frame):
    return "; ".join(f"{g}: {len(x):,} objects from {x.case_id.nunique()} donors"
                     for g in GROUPS for x in [frame[frame.donor_status.eq(g)]])


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13,
                         "axes.titlesize": 17, "axes.labelsize": 14,
                         "svg.fonttype": "none", "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    source = pd.read_csv(SOURCE)
    required = ["islet_id", "case_id", "donor_status", "object_class", "region", "structure",
                "n_recovered_cells", "n_accepted_cells", "region_area_um2", "core_area_um2",
                "core_complete", "neighborhood_complete_20", "count_fraction_percent", "density_cells_mm2"]
    assert set(required).issubset(source.columns), sorted(set(required)-set(source.columns))
    assert source.islet_id.nunique() == 4906 and len(source) == 4906*2*3
    assert source.case_id.nunique() == 20 and not source.case_id.eq(78).any()
    assert not source.duplicated(["islet_id", "region", "structure"]).any()
    assert source.groupby("case_id").donor_status.nunique().eq(1).all()
    assert source.n_recovered_cells.ge(0).all()
    assert source.n_recovered_cells.le(source.n_accepted_cells).all()
    sums = source.groupby(["islet_id", "region"]).agg(structural=("n_recovered_cells", "sum"),
                                                    accepted=("n_accepted_cells", "first"))
    assert sums.structural.le(sums.accepted).all()
    assert source.groupby(["islet_id", "region"]).n_accepted_cells.nunique().eq(1).all()
    denominator = source.n_accepted_cells.replace(0, np.nan)
    assert np.allclose(source.count_fraction_percent, 100*source.n_recovered_cells/denominator, equal_nan=True)
    assert np.allclose(source.density_cells_mm2, 1e6*source.n_recovered_cells/source.region_area_um2.replace(0,np.nan), equal_nan=True)
    assert source.count_fraction_percent.dropna().between(0,100).all()
    definitions = (
        "These are analysis-only structural-cell calls recovered from the existing rule and retained current inputs. "
        "The served v10 labels were produced with STRUCTURE_CELL_CALLS=False; the served labels and source data have not been overwritten. "
        "For this recovery, the cell centroid must lie inside an accepted structure annotation, more than 50% of the cell's total polygon area must overlap that structure, "
        "the structure's subpanel measurement must clear its existing reference level, and its creating marker must be the cell's highest expressed identity marker. "
        "The current-panel creating markers are CD31 for Blood vessel, B3TUBB for Nerve and PDPN for Lymphatic. "
        "The existing highest-marker comparison uses aligned measurements before per-marker scaling. "
        "Structural labels are mutually exclusive, with Blood vessel > Nerve > Lymphatic precedence. "
        "The original app phenotype is retained separately; recovered structural calls are not added to phenotype totals or substituted into endocrine/immune compositions. "
        "Cell/structure overlap is measured against the union of existing same-class annotations so overlapping fragments are counted only once. "
        "Accepted edge-cell exclusions and complete-core flags are retained. Islets and EndocrineClusters are included; Singles are not analysis objects. "
        "Cell identities and per-object memberships are persistent-ID joins. "
    )
    manifest, coverage = [], []
    for panel, metric, structure, region, measure, title in PANELS:
        d = source[source.structure.eq(structure) & source.region.eq(region)].copy()
        d["numerator"] = d.n_recovered_cells
        d["denominator"] = d.n_accepted_cells if measure == "percent" else d.region_area_um2
        d["value"] = d.count_fraction_percent if measure == "percent" else d.density_cells_mm2
        d["units"] = "% of accepted cells in region" if measure == "percent" else "cells/mm2 of assessed region"
        d["diameter_um"] = 2*np.sqrt(d.core_area_um2/np.pi)
        d["bin_lo_um"] = np.floor(d.diameter_um/50)*50
        d["diam_mid"] = d.bin_lo_um+25
        d["diam_bin"] = d.bin_lo_um.map(lambda x:f"[{int(x)}, {int(x+50)})")
        d["display_range"] = d.diameter_um.ge(0)&d.diameter_um.lt(350)
        d["eligible"] = d.core_complete & (d.neighborhood_complete_20 if region == "peri" else True)
        d["eligible"] &= d.denominator.gt(0)&np.isfinite(d.value)
        d["included"] = d.eligible&d.display_range
        used = d[d.included]
        assert used.case_id.nunique() == 20
        donor = used.groupby(["case_id", "donor_status", "diam_bin", "diam_mid"], observed=True).agg(
            value=("value", "mean"), n_objects=("value", "size")).reset_index()
        summary = used.groupby(["donor_status", "diam_bin", "diam_mid"], observed=True).agg(
            mean=("value", "mean"), sd=("value", "std"), n_objects=("value", "size"),
            n_donors=("case_id", "nunique")).reset_index()
        summary["lower_sd"] = summary["mean"]-summary.sd
        summary["upper_sd"] = summary["mean"]+summary.sd
        d.to_csv(DATA/f"{metric}_object_values.csv", index=False)
        donor.to_csv(DATA/f"{metric}_donor_bin_means.csv", index=False)
        summary.to_csv(DATA/f"{metric}_curve_values.csv", index=False)
        fig, ax = plt.subplots(figsize=(7.4, 5.4))
        for group in GROUPS:
            s = summary[summary.donor_status.eq(group)].sort_values("diam_mid")
            ax.errorbar(s.diam_mid, s["mean"], yerr=s.sd, color=COLORS[group], marker="o", markersize=5,
                        capsize=3, linewidth=1.7, elinewidth=1.1, label=group)
        noun = {"Blood vessel":"Blood-vessel", "Nerve":"Neural", "Lymphatic":"Lymphatic"}[structure]
        ylabel = f"{noun} cells / all {region} cells (%)" if measure == "percent" else f"{noun} cells / {region} area (cells/mm²)"
        ax.set(title=title, xlabel="Endocrine-object diameter (µm)", ylabel=ylabel, xlim=(0,350))
        ax.set_xticks(MIDS)
        ax.grid(axis="y", color="#E5E5E5", linewidth=.65)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, loc="best")
        for ext in ["png", "pdf", "svg"]:
            fig.savefig(FIG/f"{panel}.{ext}", dpi=250, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        method = (f"The numerator is the recovered {structure} cell count in the {region}. " +
                  (f"The denominator is every accepted analysed cell in that {region}; the ratio is multiplied by 100. Objects with zero accepted cells have undefined percentages and are not plotted. "
                   if measure == "percent" else
                   f"Density is the recovered cell count divided by the assessed {region} area in mm². Zero counts, including zero-cell objects with a positive assessed area, are retained. "))
        if region == "peri":
            method += ("Peri refers to the saved 20 µm band outside the core, clipped to Tissue and excluding all Islet/EndocrineCluster cores, matching the current app. "
                       "The complete 20 µm neighborhood flag is required. Cells are unique within an object's region; overlapping peri bands may include a cell in more than one object's separate value. ")
        outside = int((d.eligible&~d.display_range).sum())
        legend = (definitions+method+
                  "Equivalent diameter is 2√(core area/π). The reference-matched 50 µm bins are [0,50), [50,100), …, [300,350), displayed at 25–325 µm. "
                  f"{outside} otherwise eligible objects at ≥350 µm remain in the source CSV but are outside the plotted diameter range. "+
                  group_counts(used)+". The full cohort contains 7 ND, 7 Aab+ and 6 T1D donors. "
                  "Points show arithmetic object means within each disease-state/diameter bin and error bars are ±1 sample SD across objects, matching the existing app mean±SD plot. "
                  "Negative lower SD limits are mathematical intervals, not negative measurements. Missing bins are left missing; zero counts are retained when their denominator is defined. "
                  "These descriptive panels contain no hypothesis test or significance annotation. Exact object values, donor-bin means, plotted means/SD and per-bin donor/object counts are exported. "
                  "Passthrough overlap fractions remain separately identified measurements and are not substituted for structural-cell counts. "
                  "Method sources: scripts/uniform/typing/classify.py::structure_calls; scripts/uniform/typing/codex_typing.py; "
                  "Senior-v3 core/phenocycler/scripts/typing/subpanel.py::_structure_calls; existing app diameter binning and mean±SD display.")
        (FIG/f"{panel}_legend.txt").write_text(legend+"\n")
        manifest.append({"panel":panel, "status":"recreated", "metric":metric, "structure":structure,
                         "region":region, "endpoint":"structural-cell percentage" if measure == "percent" else "structural-cell density",
                         "x":"diam_mid", "y":"mean", "error":"sd", "units":d.units.iloc[0],
                         "aggregation":"object mean ± sample SD within disease and diameter bin",
                         "object_source_csv":f"data/{metric}_object_values.csv",
                         "donor_source_csv":f"data/{metric}_donor_bin_means.csv",
                         "curve_source_csv":f"data/{metric}_curve_values.csv", "recovery":"existing rule reapplied offline; served labels unchanged"})
        coverage.append({"panel":panel, "metric":metric, "geometric_objects":len(d),
                         "eligible_objects":int(d.eligible.sum()), "plotted_objects":len(used),
                         "n_donors":used.case_id.nunique(), "zero_count_objects_plotted":int(used.n_recovered_cells.eq(0).sum()),
                         "zero_accepted_cell_objects_plotted":int(used.n_accepted_cells.eq(0).sum()),
                         "undefined_percent_objects":int(d.n_accepted_cells.eq(0).sum()),
                         "outside_display":outside})
    pd.DataFrame(manifest).to_csv(OUT/"panel_manifest.csv", index=False)
    pd.DataFrame(coverage).to_csv(DATA/"plot_coverage.csv", index=False)
    (OUT/"panel_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    passthrough_source = OUT/"object_passthrough_summaries.csv"
    passthrough = pd.read_csv(passthrough_source)
    assert len(passthrough) == 4906*2*3
    assert not passthrough.duplicated(["islet_id","region","structure"]).any()
    passthrough["mean_cell_area_overlap_percent_among_measured_passthrough"] = 100*passthrough.mean_passthrough_fraction_among_measured
    passthrough["percent_accepted_cells_with_measured_passthrough"] = 100*passthrough.n_passthrough_measured/passthrough.n_accepted_cells.replace(0,np.nan)
    passthrough["percent_centroid_inside_cells_with_measured_passthrough"] = 100*passthrough.n_passthrough_measured/passthrough.n_centroid_inside.replace(0,np.nan)
    for c in ["mean_cell_area_overlap_percent_among_measured_passthrough", "percent_accepted_cells_with_measured_passthrough", "percent_centroid_inside_cells_with_measured_passthrough"]:
        assert passthrough[c].dropna().between(-1e-8,100+1e-6).all()
    passthrough.to_csv(DATA/"passthrough_percentage_object_values.csv",index=False)
    (DATA/"passthrough_percentages_methods.txt").write_text(
        "Passthrough remains a separate measurement from structural-cell identity. For a given structure, a measured passthrough occurs when the centroid lies inside the structure but the cell does not receive that structure's own call. Its fraction is union(cell polygon ∩ same-class structure annotations) / total cell polygon area; multiplying by100 gives percent cell area. Values remain NA outside that condition or where no overlap measurement is available.\n\n"
        "mean_cell_area_overlap_percent_among_measured_passthrough:100 × the arithmetic mean of measured cell-area overlap fractions among that object's region-specific passthrough cells. Denominator:n_passthrough_measured; undefined if zero.\n"
        "percent_accepted_cells_with_measured_passthrough:100 × n_passthrough_measured / n_accepted_cells in that object's core or peri region; undefined if zero accepted cells. This is a cell frequency, not cell-area coverage.\n"
        "percent_centroid_inside_cells_with_measured_passthrough:100 × n_passthrough_measured / n_centroid_inside in that object's region; undefined if no centroid-inside cells. This is a cell frequency, not cell-area coverage.\n\n"
        "All4,906 geometric objects, both regions and all three structures are exported with eligibility flags; no disease summary, hypothesis test or companion plot is assigned. Per-cell fractions are retained by the recovery producer under per_cell/. Peri memberships are unique within each object and may recur between overlapping bands.\n")
    sources = [SOURCE,passthrough_source,OUT/"analysis_specification.json",Path(__file__).resolve(),ROOT/"scripts/uniform/typing/classify.py",ROOT/"scripts/uniform/typing/codex_typing.py"]
    summary = {"source_sha256":{str(p):digest(p) for p in sources}, "figures":8, "donors":20,
               "objects":4906,"interpretation":"none assigned","statistics":"descriptive object mean ± SD; no tests",
               "served_data_changed":False,"area_assets_changed":False,"quantitative_INS":False,
               "files":{str(p.relative_to(OUT)):digest(p) for folder in [FIG,DATA] for p in sorted(folder.glob("*")) if p.is_file()}}
    (OUT/"figure_manifest.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(pd.DataFrame(coverage).to_string(index=False),flush=True)
    print("Eight structural-cell panels complete",flush=True)


if __name__ == "__main__":
    main()
