#!/usr/bin/env python3
"""Compare existing donor capillary-crossing density differences across disease states."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/h3-disease-20261026-mpl")
import hashlib
import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[2]
SOURCE = REPO / "output/h1_h4_literature_methods_20261025/h3/h3_donor_densities.csv"
OUT = REPO / "output/h1_h4_disease_comparisons_20261026/h3"
GROUPS = ["ND", "Aab+", "T1D"]
COLORS = ["#4477AA", "#CC6633", "#228833"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source_hash = sha(SOURCE)
    all_widths = pd.read_csv(SOURCE)
    d = all_widths.loc[all_widths.width_um.eq(10)].copy()
    assert len(d) == 20 and d.case_id.nunique() == 20 and 78 not in d.case_id.values
    assert d.groupby("donor_status").size().to_dict() == {"ND": 7, "Aab+": 7, "T1D": 6}
    assert np.allclose(d.difference_density_mm2, d.near_density_mm2 - d.reference_density_mm2)
    d.to_csv(OUT / "h3_primary_donor_values.csv", index=False)

    samples = [d.loc[d.donor_status.eq(g), "difference_density_mm2"].to_numpy() for g in GROUPS]
    kw = stats.kruskal(*samples)
    pd.DataFrame([dict(endpoint="Near-minus-reference immune density at 10 um",
        test="Kruskal-Wallis", donors=20, groups=3, statistic=float(kw.statistic),
        degrees_of_freedom=2, p_value=float(kw.pvalue))]).to_csv(OUT / "h3_global_test.csv", index=False)
    rows = []
    for i, j in itertools.combinations(range(3), 2):
        # Standard rule shared across hypotheses: exact if untied, otherwise
        # tie-corrected normal approximation with continuity correction.
        combined = np.concatenate([samples[i], samples[j]])
        exact_available = len(np.unique(combined)) == len(combined)
        method = "exact" if exact_available else "asymptotic"
        result = stats.mannwhitneyu(samples[i], samples[j], alternative="two-sided",
                                   method=method, use_continuity=True)
        rows.append(dict(group_1=GROUPS[i], group_2=GROUPS[j], donors_1=len(samples[i]),
            donors_2=len(samples[j]), statistic_U=float(result.statistic), p_value=float(result.pvalue),
            method="Exact two-sided Wilcoxon rank-sum" if exact_available else "Tie-corrected normal approximation with continuity correction",
            median_difference_group2_minus_group1=float(np.median(samples[j]) - np.median(samples[i]))))
    comparisons = pd.DataFrame(rows)
    comparisons["p_bh"] = multipletests(comparisons.p_value, method="fdr_bh")[1]
    comparisons.to_csv(OUT / "h3_pairwise_disease_tests.csv", index=False)

    context = all_widths.groupby(["width_um", "donor_status"]).agg(
        donors=("case_id", "size"), objects=("objects", "sum"),
        median_difference_mm2=("difference_density_mm2", "median"),
        mean_difference_mm2=("difference_density_mm2", "mean"),
        q25_difference_mm2=("difference_density_mm2", lambda x: x.quantile(.25)),
        q75_difference_mm2=("difference_density_mm2", lambda x: x.quantile(.75))).reset_index()
    context.to_csv(OUT / "h3_width_context_descriptive.csv", index=False)

    fig, ax = plt.subplots(figsize=(5.5, 4.8))
    boxes = ax.boxplot(samples, positions=np.arange(3), widths=.48, patch_artist=True,
        showfliers=False, medianprops=dict(color="#222222", linewidth=1.5),
        whiskerprops=dict(color="#555555", linewidth=1.1), capprops=dict(color="#555555"))
    for box, color in zip(boxes["boxes"], COLORS):
        box.set_facecolor(color)
        box.set_alpha(.18)
        box.set_edgecolor(color)
    for i, (values, color) in enumerate(zip(samples, COLORS)):
        # Fixed symmetric jitter displays every donor without altering group summaries.
        jitter = np.linspace(-.09, .09, len(values))
        ax.scatter(i + jitter, values, s=45, color=color, edgecolor="white", linewidth=.6, zorder=3)
    ax.axhline(0, color="#888888", linewidth=.8, linestyle="--", zorder=0)
    ax.set_xticks(np.arange(3), GROUPS)
    ax.set_xlabel("Disease state")
    ax.set_ylabel("Immune-density difference\nnear crossing − reference (cells/mm²)")
    ax.spines[["top", "right"]].set_visible(False)
    span = max(d.difference_density_mm2.max() - d.difference_density_mm2.min(), 1)
    top = d.difference_density_mm2.max()
    for r in comparisons.loc[comparisons.p_bh.lt(.05)].itertuples():
        i, j = GROUPS.index(r.group_1), GROUPS.index(r.group_2)
        base, height = top + .10 * span, .035 * span
        ax.plot([i, i, j, j], [base, base + height, base + height, base], color="#222222", lw=1)
        star = "***" if r.p_bh < .001 else "**" if r.p_bh < .01 else "*"
        ax.text((i + j) / 2, base + height, star, ha="center", va="bottom", fontsize=13)
        top += .16 * span
    ax.set_ylim(d.difference_density_mm2.min() - .12 * span, top + .16 * span)
    fig.tight_layout()
    fig.savefig(OUT / "h3_density_difference_by_disease.png", dpi=220)
    fig.savefig(OUT / "h3_density_difference_by_disease.svg")
    plt.close(fig)

    pair_text = "; ".join(f"{r.group_1} versus {r.group_2}, adjusted p={r.p_bh:.3f}"
                          for r in comparisons.itertuples())
    legend = (
        "Immune-cell localization near capillary boundary crossings across disease states. "
        "Each point represents one donor: ND, n=7; Aab+, n=7; T1D, n=6. "
        "The plotted value is immune-cell density near capillary crossings minus density in the "
        "remaining peripheral core region. Both regions are restricted to 10 µm inside the endocrine "
        "core boundary; the near region lies within 10 µm of a capillary–boundary intersection. "
        "For each donor, counts and areas were summed across eligible islets/clusters before calculating "
        "each region's density. Zero denotes equal density; positive values indicate greater density near "
        "crossings. Boxes show the median and interquartile range; whiskers extend to the most extreme "
        "values within 1.5 times the interquartile range. All donor values are displayed. "
        f"The disease groups were compared by Kruskal–Wallis test (H={kw.statistic:.3f}, df=2, "
        f"p={kw.pvalue:.3f}), followed by two-sided Wilcoxon rank-sum comparisons with "
        f"Benjamini–Hochberg adjustment across three comparisons: {pair_text}. "
        "Rank-sum p values use the exact calculation for untied values and a tie-corrected "
        "normal approximation with continuity correction otherwise. "
        "Brackets and asterisks indicate adjusted p<.05 (*), p<.01 (**), or p<.001 (***); "
        "none of the pairwise comparisons meets this criterion."
    )
    (OUT / "figure_legend.txt").write_text(legend + "\n")
    report = ["Capillary-crossing immune localization: comparison across disease states", "",
        "Primary question: does the donor-level near-minus-reference density difference at 10 µm differ between ND, Aab+, and T1D?",
        "Existing measurements were reused; no geometry, source cells, or data pipeline was rebuilt.", "",
        f"Kruskal-Wallis: H={kw.statistic:.6f}, df=2, p={kw.pvalue:.6f}.", "",
        "Group summaries (cells/mm²):"]
    for g in GROUPS:
        r = context.loc[context.width_um.eq(10) & context.donor_status.eq(g)].iloc[0]
        report.append(f"{g}: n={int(r.donors)}, median={r.median_difference_mm2:.3f}, "
                      f"IQR={r.q25_difference_mm2:.3f} to {r.q75_difference_mm2:.3f}, "
                      f"mean={r.mean_difference_mm2:.3f}.")
    report += ["", "Pairwise Wilcoxon rank-sum comparisons (two-sided; BH over three comparisons):"]
    for r in comparisons.itertuples():
        report.append(f"{r.group_1} vs {r.group_2}: U={r.statistic_U:.1f}, p={r.p_value:.6f}, adjusted p={r.p_bh:.6f}; {r.method}.")
    report += ["", "Conclusion: these data do not establish a difference in capillary-associated immune localization across disease states. "
        "The previously observed pooled near-versus-reference association answers a different question and does not establish a disease-group difference.", "",
        "The 5 and 20 µm measurements are retained as descriptive context only; no new group tests were run at those widths. "
        "The endpoint is a difference between observed spatial compartments, not a conditional random spatial null.", "",
        "Method precedents: literature/Richardson2023.pdf, PDF p8 (donor-group Kruskal-Wallis); "
        "literature/Bodenmiller2025.pdf, PDF p28 (donor-aggregated Wilcoxon rank-sum and FDR adjustment). "
        "Boxplots with individual donor points are also described in Bodenmiller2025, PDF pp17-18."]
    (OUT / "RESULTS.txt").write_text("\n".join(report) + "\n")
    provenance = dict(source=str(SOURCE), source_sha256=source_hash,
        primary_width_um=10, donor_unit="One near-minus-reference immune-density difference per donor",
        global_test="Kruskal-Wallis", pairwise_test="Two-sided Wilcoxon rank-sum / Mann-Whitney U",
        rank_sum_settings="exact for untied combined samples; tie-corrected normal approximation with continuity correction otherwise",
        multiplicity="Benjamini-Hochberg over three pairwise disease comparisons",
        source_unchanged=sha(SOURCE) == source_hash, donor78_excluded=True,
        analysis_script=str(Path(__file__).resolve()), analysis_script_sha256=sha(__file__))
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    main()
