#!/usr/bin/env python3
"""Recreate quantitative reference-deck panels with accepted 20-donor counts/areas.

No cell typing, structure classification, embeddings or signal thresholds are fit.
Per-bin tests call the existing app's R helpers, with its multiplicity families.
"""
from pathlib import Path
import hashlib
import json
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/pancfloyd_20donor_20261026/slides1_4"
FIG = OUT / "figures"
DATA = OUT / "data"
GROUPS = ["ND", "Aab+", "T1D"]
COLORS = dict(zip(GROUPS, ["#3366CC", "#F07618", "#008D28"]))
MIDS = np.arange(25, 350, 50)
COUNT_PATH = ROOT / "output/h1_edge_corrected_20261024/h1_object_features.parquet"
CENSUS_PATH = ROOT / "output/h5_h10_original_plan_20261026/h8/full_object_census.parquet"
TARGET_PATH = ROOT / "tmp/ufrc-release-20261022/artifacts/release-d55680a/data/summary/targets.parquet"
UNION_PATH = DATA / "structure_union_area_values.csv"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13,
                     "axes.titlesize": 17, "axes.labelsize": 14,
                     "svg.fonttype": "none", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})


def savefig(fig, panel):
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(FIG / f"{panel}.{ext}", dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def bin_data(df):
    df = df.copy()
    df["diameter_um"] = 2 * np.sqrt(df.core_area_um2 / np.pi)
    df["bin_lo_um"] = np.floor(df.diameter_um / 50) * 50
    df["diam_mid"] = df.bin_lo_um + 25
    df["diam_bin"] = df.bin_lo_um.map(lambda x: f"[{int(x)}, {int(x+50)})")
    df["display_range"] = (df.diameter_um >= 0) & (df.diameter_um < 350)
    df["eligible"] = df.eligible & np.isfinite(df.value) & (df.denominator > 0)
    df["included"] = df.eligible & df.display_range
    return df


def summarize(metric, df):
    assert not df.islet_id.duplicated().any(), f"{metric}: duplicate object rows"
    finite=df.value[np.isfinite(df.value)]
    assert finite.between(-1e-7,100+1e-6).all(), f"{metric}: physical percentage outside 0–100%"
    df.to_csv(DATA / f"{metric}_object_values.csv", index=False)
    used = df[df.included].copy()
    donor = used.groupby(["case_id", "donor_status", "diam_bin", "diam_mid"], observed=True).agg(
        value=("value", "mean"), n_objects=("value", "size")).reset_index()
    donor.to_csv(DATA / f"{metric}_donor_bin_means.csv", index=False)
    summary = used.groupby(["donor_status", "diam_bin", "diam_mid"], observed=True).agg(
        mean=("value", "mean"), sd=("value", "std"), n_objects=("value", "size"),
        n_donors=("case_id", "nunique")).reset_index()
    summary["lower_sd"] = summary["mean"] - summary.sd
    summary["upper_sd"] = summary["mean"] + summary.sd
    summary.to_csv(DATA / f"{metric}_curve_values.csv", index=False)
    return used, donor, summary


def count_text(df):
    return "; ".join(f"{g}: {len(x):,} objects from {x.case_id.nunique()} donors"
                     for g in GROUPS for x in [df[df.donor_status == g]])


def base_legend(df, definition, peri):
    used = df[df.included]
    outside = int((df.eligible & ~df.display_range).sum())
    return (definition + " Accepted uniform-20donor-20261022-v10 phenotypes and existing app structure annotations were retained. "
            "Endocrine objects include Islets and EndocrineClusters; Singles are not included. "
            "Core completeness follows the corrected acquisition-edge flags. "
            + ("The full 20 µm neighborhood must additionally be complete. " if peri else "")
            + "Diameter is 2√(core area/π). Bins are [0,50), [50,100), …, [300,350) µm, shown at their midpoints; "
            f"{outside} otherwise eligible objects at ≥350 µm are outside this reference-matched display range and remain in the object CSV. "
            + count_text(used) + ". ")


def curve(panel, metric, title, ylabel, df, summary, definition, peri=False):
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    for g in GROUPS:
        x = summary[summary.donor_status == g].sort_values("diam_mid")
        ax.errorbar(x.diam_mid, x["mean"], yerr=x.sd, color=COLORS[g], marker="o",
                    markersize=5, capsize=3, linewidth=1.7, elinewidth=1.1, label=g)
    ax.set(title=title, xlabel="Endocrine-object diameter (µm)", ylabel=ylabel, xlim=(0,350))
    ax.set_xticks(MIDS)
    ax.grid(axis="y", color="#E5E5E5", linewidth=.65)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="best")
    savefig(fig, panel)
    legend = base_legend(df, definition, peri) + (
        "Points are arithmetic means across objects in each disease-state/diameter bin; error bars are ±1 SD across objects, "
        "matching the reference curve and the app's mean±SD display option. Negative lower SD limits are mathematical intervals, not negative measurements. "
        "No hypothesis test is encoded on this descriptive curve. Exact object values, donor-bin means and plotted summaries are exported separately.")
    (FIG / f"{panel}_legend.txt").write_text(legend + "\n")


def run_app_tests(metrics):
    # Source the established app code rather than implementing a new test recipe.
    r_script = r'''
    suppressPackageStartupMessages({library(dplyr); library(tibble)})
    args <- commandArgs(trailingOnly=TRUE)
    root <- args[1]; folder <- args[2]
    source(file.path(root, "app/shiny_app/R/utils_stats.R"))
    for (metric in args[-c(1,2)]) {
      d <- read.csv(file.path(folder, paste0(metric,"_object_values.csv")), check.names=FALSE)
      d <- d[tolower(as.character(d$included)) == "true", ]; d$`Case ID` <- d$case_id
      d$donor_status <- factor(d$donor_status, levels=c("ND","Aab+","T1D"))
      normality <- normality_tests(aggregate_to_donor(d, "mean"))
      parametric <- if (all(is.na(normality$is_normal))) FALSE else isTRUE(all(normality$is_normal,na.rm=TRUE))
      normality$selected_pairwise <- if(parametric) "Welch t test" else "Wilcoxon rank-sum test"
      a <- per_bin_donor_anova(d,"diam_bin","donor_status","value")
      a$q_anova <- p.adjust(a$p_anova,"BH")
      k <- per_bin_donor_kendall(d,"diam_bin","donor_status","value")
      k$q_kendall <- p.adjust(k$p_kendall,"BH")
      p <- per_bin_donor_pairwise(d,"diam_bin","donor_status","value",parametric=parametric)
      p <- p %>% group_by(bin) %>% mutate(q_within_bin=p.adjust(p_value,"BH")) %>% ungroup()
      p$protected_bin <- p$mid %in% a$mid[is.finite(a$q_anova) & a$q_anova < .05]
      p$test <- if(parametric) "Welch t test" else "Wilcoxon rank-sum test"
      write.csv(normality,file.path(folder,paste0(metric,"_normality.csv")),row.names=FALSE)
      write.csv(a,file.path(folder,paste0(metric,"_anova.csv")),row.names=FALSE)
      write.csv(k,file.path(folder,paste0(metric,"_kendall.csv")),row.names=FALSE)
      write.csv(p,file.path(folder,paste0(metric,"_pairwise.csv")),row.names=FALSE)
    }
    '''
    script = OUT / "app_bin_tests.R"
    script.write_text(r_script)
    process = subprocess.run(["Rscript", str(script), str(ROOT), str(DATA), *metrics],
                             text=True, capture_output=True)
    (OUT / "statistics_run.log").write_text(process.stdout + process.stderr)
    process.check_returncode()


def auc_plot(panel, metric, title, df, summary):
    rows = []; full_rows = []
    shared_midpoints = sorted(set.intersection(*[set(summary.loc[summary.donor_status == g, "diam_mid"]) for g in GROUPS]))
    for g in GROUPS:
        x = summary[summary.donor_status == g].sort_values("diam_mid")
        full_rows.append({"donor_status":g, "auc_percent_um":np.trapezoid(x["mean"], x.diam_mid),
                          "first_mid_um":x.diam_mid.min(), "last_mid_um":x.diam_mid.max(),
                          "n_bins":len(x), "n_objects":int(x.n_objects.sum())})
        x = x[x.diam_mid.isin(shared_midpoints)]
        rows.append({"donor_status":g, "auc_percent_um":np.trapezoid(x["mean"], x.diam_mid),
                     "first_mid_um":x.diam_mid.min(), "last_mid_um":x.diam_mid.max(),
                     "n_bins":len(x), "n_objects":int(x.n_objects.sum())})
    values = pd.DataFrame(rows)
    values.to_csv(DATA / f"{metric}_auc_values.csv", index=False)
    pd.DataFrame(full_rows).to_csv(DATA / f"{metric}_auc_full_available_range.csv", index=False)
    fig, ax = plt.subplots(figsize=(5.8,5.4))
    ax.bar(GROUPS, values.auc_percent_um, color=[COLORS[g] for g in GROUPS], width=.62)
    ax.set(title=title, ylabel="Area under fraction–diameter curve (%·µm)")
    ax.grid(axis="y", color="#E5E5E5", linewidth=.65); ax.set_axisbelow(True)
    savefig(fig,panel)
    spans = "; ".join(f"{r.donor_status}: {r.first_mid_um:g}–{r.last_mid_um:g} µm, {r.n_bins} bins" for r in values.itertuples())
    legend = ("Bars show the trapezoidal integral of the disease-state object-mean fraction curve over the same bin midpoints observed in all three disease states, "
              "using the app's trapezoidal calculation on this common domain. Values use accepted typed-cell counts as a percentage of all analysed core cells; they do not quantify hormone signal. "
              + spans + ". " + count_text(df[df.included & df.diam_mid.isin(shared_midpoints)]) + ". "
              "This is a descriptive integral of each group curve, not an individual-donor distribution; no uncertainty or significance test is attached. "
              "Missing bins are not filled. The main AUC CSV records the common integration limits and values; a companion CSV preserves the full available-range app integral for each group. The corresponding curve legend specifies edge eligibility and diameter range.")
    (FIG / f"{panel}_legend.txt").write_text(legend+"\n")


def heatmap(panel, metric, title, df):
    a = pd.read_csv(DATA / f"{metric}_anova.csv")
    k = pd.read_csv(DATA / f"{metric}_kendall.csv")
    p = pd.read_csv(DATA / f"{metric}_pairwise.csv")
    norm = pd.read_csv(DATA / f"{metric}_normality.csv")
    rows = ["ANOVA", "ND vs Aab+", "ND vs T1D", "Aab+ vs T1D", "Kendall τ"]
    q = np.full((5,len(MIDS)),np.nan)
    for j,mid in enumerate(MIDS):
        ar=a[a.mid==mid]; kr=k[k.mid==mid]
        if len(ar): q[0,j]=ar.q_anova.iloc[0]
        if len(kr): q[4,j]=kr.q_kendall.iloc[0]
        for i,pair in enumerate(rows[1:4],1):
            z=p[(p.mid==mid)&(p.pair==pair)&p.protected_bin]
            if len(z): q[i,j]=z.q_within_bin.iloc[0]
    plot_rows=[]
    for i,row in enumerate(rows):
        for j,mid in enumerate(MIDS):
            plot_rows.append({"comparison":row,"diam_mid":mid,"displayed_q":q[i,j],
                              "minus_log10_q":-np.log10(q[i,j]) if np.isfinite(q[i,j]) else np.nan})
    pd.DataFrame(plot_rows).to_csv(DATA/f"{metric}_heatmap_values.csv",index=False)
    value=-np.log10(q)
    cmap=plt.get_cmap("magma").copy();cmap.set_bad("#E5E5E5")
    fig,ax=plt.subplots(figsize=(8.7,4.7))
    im=ax.imshow(np.ma.masked_invalid(value),aspect="auto",cmap=cmap,norm=Normalize(0,max(2,float(np.nanmax(value)))))
    for i in range(5):
        for j in range(len(MIDS)):
            if np.isfinite(q[i,j]) and q[i,j]<.05:
                star="***" if q[i,j]<.001 else "**" if q[i,j]<.01 else "*"
                ax.text(j,i,star,ha="center",va="center",fontsize=17,color="white" if value[i,j]<np.nanmax(value)*.65 else "black")
    ax.set(xticks=np.arange(len(MIDS)),xticklabels=MIDS,yticks=np.arange(5),yticklabels=rows,
           xlabel="Endocrine-object diameter bin midpoint (µm)",title=title)
    cb=fig.colorbar(im,ax=ax,pad=.025);cb.set_label("−log₁₀(q)")
    savefig(fig,panel)
    test=norm.selected_pairwise.iloc[0]
    legend=("Each test uses one arithmetic object mean per donor within a 50 µm diameter bin. "
            "At least two disease groups must be present, with at least two donors in every present group, following the existing app coverage rule. "
            "The ANOVA row is an ordinary one-way ANOVA of donor-bin means. The Kendall row tests ordered disease codes ND=0, Aab+=1, T1D=2 against donor-bin means, using the asymptotic Kendall test. "
            "Benjamini–Hochberg correction is applied separately across eligible diameter bins for the ANOVA and Kendall families. "
            f"Pairwise rows use {test}, selected by the app's Shapiro–Wilk rule on whole-range donor means (all available group tests p>0.05 selects Welch). "
            "Wilcoxon tests, when selected, use exact p values when untied and the tie-corrected approximation with continuity correction otherwise. "
            "Pairwise rows are displayed only where the bin's BH-adjusted ANOVA q<0.05; their three comparisons are BH-adjusted within that bin. "
            "Grey cells indicate unavailable tests or pairwise tests withheld by this prespecified omnibus display rule, not q=1. "
            "Asterisks indicate *q<0.05, **q<0.01, ***q<0.001. All raw p values, adjusted q values, normality decisions and donor-bin values are exported. "
            + count_text(df[df.included])+". Fractions use typed-cell counts divided by all analysed core cells. The curve legend specifies edge eligibility and display range.")
    (FIG/f"{panel}_legend.txt").write_text(legend+"\n")


def main():
    FIG.mkdir(parents=True,exist_ok=True); DATA.mkdir(parents=True,exist_ok=True)
    counts=pd.read_parquet(COUNT_PATH); census=pd.read_parquet(CENSUS_PATH); targets=pd.read_csv(UNION_PATH)
    assert counts.case_id.nunique()==20 and 78 not in counts.case_id.values
    assert not counts.islet_id.duplicated().any()
    immune=["CD8a_Tcell","CD4_Tcell","T_cell","B_cell","Macrophage","Immune"]
    assert np.array_equal(counts[[f"n_core_{x}" for x in immune]].sum(axis=1),counts.n_immune_core)
    assert np.array_equal(counts[[f"n_peri_{x}" for x in immune]].sum(axis=1),counts.n_immune_peri)
    specs={
      "alpha_core":("slide02_image01","α-cell fraction","α cells / all core cells (%)",["Alpha_cell"],"core"),
      "beta_core":("slide02_image02","β-cell fraction","β cells / all core cells (%)",["Beta_cell"],"core"),
      "delta_core":("slide08_image01","δ-cell fraction","δ cells / all core cells (%)",["Delta_cell"],"core"),
      "cd8_core_peri":("slide05_image06","CD8⁺ T-cell fraction · core + peri","CD8⁺ T cells / all core + peri cells (%)",["CD8a_Tcell"],"core_peri"),
      "immune_core_peri":("slide10_image01","Immune-cell fraction · core + peri","Immune cells / all core + peri cells (%)",immune,"core_peri"),
      "tcell_core_peri":("slide10_image02","T-cell fraction · core + peri","T cells / all core + peri cells (%)",["CD8a_Tcell","CD4_Tcell","T_cell"],"core_peri"),
    }
    metrics={}; manifest=[]; coverage=[]
    common=["islet_id","case_id","donor_status","object_class","core_area_um2","core_complete","neighborhood_complete_20"]
    for metric,(panel,title,ylabel,types,region) in specs.items():
        d=counts[common].copy(); peri=region=="core_peri"
        d["region"]=region
        d["numerator"]=counts[[f"n_core_{x}" for x in types]].sum(axis=1)
        d["denominator"]=counts.n_cells_core
        if peri:
            d["numerator"]+=counts[[f"n_peri_{x}" for x in types]].sum(axis=1)
            d["denominator"]+=counts.n_cells_peri
        d["value"]=100*d.numerator/d.denominator
        d["eligible"]=d.core_complete & (d.neighborhood_complete_20 if peri else True)
        d=bin_data(d); used,donor,summary=summarize(metric,d);metrics[metric]=(d,summary)
        definition=f"The numerator is the mutually exclusive accepted {', '.join(types)} label count; the denominator is all analysed {('core + 20 µm peri' if peri else 'core')} cells. Values are percentages. "
        if peri: definition+="The peri band excludes endocrine cores and retains existing per-object membership; shared peri cells can belong to more than one object's band. "
        curve(panel,metric,title,ylabel,d,summary,definition,peri)
        manifest.append({"panel":panel,"status":"recreated","metric":metric,"endpoint":"typed-cell fraction","replacement":"count-based accepted phenotype; no quantitative hormone signal"})
        coverage.append({"metric":metric,"source_objects":len(d),"edge_eligible_objects":int(d.eligible.sum()),"plotted_objects":len(used),"donors":used.case_id.nunique(),"outside_display":int((d.eligible&~d.display_range).sum())})
    geometry=census[census.object_class.isin(["Islet","EndocrineCluster"])][common].copy()
    assert len(geometry)==4906 and not geometry.islet_id.duplicated().any()
    areas=[("Capillary","core","slide03_image05"),("Capillary","peri","slide03_image06"),
           ("Nerve","core","slide04_image03"),("Nerve","peri","slide04_image06"),
           ("Lymphatic","core","slide09_image03"),("Lymphatic","peri","slide09_image04")]
    for kind,region,panel in areas:
        metric=f"{kind.lower()}_{region}_area";peri=region=="peri"
        t=targets[(targets["class"]==kind)&(targets.type==("islet_band" if peri else "islet_core"))]
        assert len(t)==4906 and not t.islet_id.duplicated().any()
        d=geometry.merge(t[["islet_id","area_um2","region_um2"]],on="islet_id",validate="one_to_one")
        d["region"]=region;d["numerator"]=d.area_um2;d["denominator"]=d.region_um2
        d["value"]=100*d.numerator/d.denominator
        d["eligible"]=d.core_complete & (d.neighborhood_complete_20 if peri else True)
        d=bin_data(d); used,donor,summary=summarize(metric,d);metrics[metric]=(d,summary)
        definition=(f"{kind} area is the union of the existing app's {kind.lower()} annotations, clipped to the assessed {'20 µm peri-band' if peri else 'endocrine-core'} region, divided by that region's area and multiplied by 100. Overlapping annotations contribute physical area once; no values are clamped. "
                    "This is a structure-area measurement, not a typed-cell percentage or cell density. Geometric objects with an area denominator remain eligible even if they have no paired cell-count row. ")
        if peri: definition += "The exact saved app peri definition is retained: tissue-clipped expanded area minus Islet/Cluster cores; Single-cell endocrine objects remain in this app region. "
        curve(panel,metric,f"{kind} area · {region}",f"{kind} area / {region} area (%)",d,summary,definition,peri)
        manifest.append({"panel":panel,"status":"recreated","metric":metric,"endpoint":"structure union area percent","replacement":"existing current app structures; overlapping annotation area counted once"})
        coverage.append({"metric":metric,"source_objects":len(d),"edge_eligible_objects":int(d.eligible.sum()),"plotted_objects":len(used),"donors":used.case_id.nunique(),"outside_display":int((d.eligible&~d.display_range).sum())})
    run_app_tests(["alpha_core","beta_core","delta_core"])
    for metric,ap,hp,cell in [("beta_core","slide02_image03","slide02_image05","β"),("alpha_core","slide02_image04","slide02_image06","α"),("delta_core","slide08_image02","slide08_image03","δ")]:
        d,s=metrics[metric]
        auc_plot(ap,metric,f"{cell}-cell fraction across diameters",d,s)
        heatmap(hp,metric,f"{cell}-cell fraction across diameter bins",d)
        manifest.extend([{"panel":ap,"status":"recreated","metric":metric,"endpoint":"descriptive fraction-curve AUC"},
                         {"panel":hp,"status":"recreated","metric":metric,"endpoint":"existing app donor-bin statistics"}])
    structural_cell_panels={"slide03_image01":"slide03_image05","slide03_image02":"slide03_image06",
                            "slide03_image03":"slide03_image05","slide03_image04":"slide03_image06",
                            "slide04_image04":"slide04_image03","slide04_image05":"slide04_image06",
                            "slide09_image01":"slide09_image03","slide09_image02":"slide09_image04"}
    for panel,area_panel in structural_cell_panels.items():
        manifest.append({"panel":panel,"status":"superseded_by_structural_cell_panel","current_panel":panel,
                         "current_manifest":"../structural_cells/panel_manifest.csv","separate_area_panel":area_panel,
                         "reason":"Structural object density/area was an app display choice, not removal of the base structural-cell rule. The original cell panel is recreated separately from retained inputs using the established centroid, cell-area containment and marker/highest-marker steps; passthrough measurements are retained separately. Structural area remains a distinct endpoint."})
    for panel in ["slide03_image07","slide03_image08"]:
        manifest.append({"panel":panel,"status":"no_current_equivalent","reason":"Reference SMA+ arteriole typed-cell percent is not an accepted current cell label or current served marker-state call. No new marker gate or structure model was introduced."})
    pd.DataFrame(coverage).to_csv(DATA/"coverage.csv",index=False)
    pd.DataFrame(manifest).to_csv(OUT/"panel_manifest.csv",index=False)
    specification={"data_root":"uniform-20donor-20261022-v10","donors":20,"display_diameter_um":[0,350],"bin_width_um":50,
      "curves":"object arithmetic mean ± sample SD, existing app mean_sd option",
      "inference":"existing app donor-bin ANOVA, Kendall and protected pairwise helpers; no biological interpretation",
      "quantitative_INS":False,"count_labels_retained":True,"input_files":[str(COUNT_PATH),str(CENSUS_PATH),str(UNION_PATH)],
      "structure_area":"Union of existing accepted annotations clipped to current regions, not sum of overlapping annotations",
      "auc_primary_domain_um":[25,225],"auc_full_available_range_companion_retained":True,
      "source_code":["app/shiny_app/R/data_loading.R::bin_islet_sizes","app/shiny_app/R/utils_stats.R::summary_stats",
                     "app/shiny_app/R/utils_stats.R::per_bin_donor_anova","app/shiny_app/R/utils_stats.R::per_bin_donor_kendall",
                     "app/shiny_app/R/utils_stats.R::per_bin_donor_pairwise","app/shiny_app/R/mod_statistics_server.R::AUC"],
      "native_figure_count":18,"figures":[p.stem for p in sorted(FIG.glob("*.png"))]}
    (OUT/"analysis_specification.json").write_text(json.dumps(specification,indent=2)+"\n")
    hashes={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.rglob("*")) if p.is_file() and p.name!="file_manifest.json"}
    (OUT/"file_manifest.json").write_text(json.dumps(hashes,indent=2)+"\n")
    print(pd.DataFrame(coverage).to_string(index=False));print(f"Saved {len(list(FIG.glob('*.png')))} native figures to {FIG}")


if __name__=="__main__":
    main()
