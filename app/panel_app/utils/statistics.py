"""
Statistical analysis utilities for Islet Explorer Panel App.
Provides ANOVA, pairwise comparisons, and AUC analysis.
"""

import numpy as np
import pandas as pd
from typing import Dict
from scipy import stats
from scipy.stats import f_oneway, kruskal, ttest_ind
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)


class StatisticalAnalyzer:
    """Performs statistical analyses on islet data."""

    def __init__(self, alpha: float = 0.05):
        self.alpha = alpha
        self.results = {}

    def run_anova(
        self,
        df: pd.DataFrame,
        value_col: str,
        group_col: str = "Donor Status",
        remove_outliers: bool = False,
        outlier_threshold: float = 3.0
    ) -> Dict:
        """Run one-way ANOVA with optional outlier removal."""

        analysis_df = df.dropna(subset=[value_col, group_col]).copy()

        if analysis_df.empty:
            return {"error": "No valid data for analysis"}

        # Outlier removal
        if remove_outliers:
            mean = analysis_df[value_col].mean()
            std = analysis_df[value_col].std()
            mask = np.abs(analysis_df[value_col] - mean) <= outlier_threshold * std
            n_removed = (~mask).sum()
            analysis_df = analysis_df[mask]
        else:
            n_removed = 0

        # Get groups
        groups = analysis_df.groupby(group_col)[value_col].apply(list).to_dict()

        if len(groups) < 2:
            return {"error": "Need at least 2 groups for ANOVA"}

        # Check assumptions
        group_sizes = {k: len(v) for k, v in groups.items()}
        min_size = min(group_sizes.values())

        # Normality test (Shapiro-Wilk) for each group
        normality_tests = {}
        for name, values in groups.items():
            if len(values) >= 3:
                stat, p = stats.shapiro(values)
                normality_tests[name] = {"statistic": stat, "p_value": p, "normal": p > 0.05}
            else:
                normality_tests[name] = {"statistic": None, "p_value": None, "normal": None}

        # Levene's test for homogeneity of variances
        if min_size >= 2:
            levene_stat, levene_p = stats.levene(*groups.values())
            equal_variances = levene_p > 0.05
        else:
            levene_stat, levene_p = None, None
            equal_variances = None

        # Decide between parametric and non-parametric
        all_normal = all(t.get("normal", False) for t in normality_tests.values() if t["normal"] is not None)
        use_parametric = all_normal and equal_variances

        # Run ANOVA or Kruskal-Wallis
        if use_parametric:
            f_stat, p_value = f_oneway(*groups.values())
            test_name = "One-way ANOVA"
        else:
            f_stat, p_value = kruskal(*groups.values())
            test_name = "Kruskal-Wallis H-test"

        # Effect size (eta-squared for ANOVA)
        if use_parametric:
            ss_between = sum(len(g) * (np.mean(g) - analysis_df[value_col].mean())**2 for g in groups.values())
            ss_total = sum((analysis_df[value_col] - analysis_df[value_col].mean())**2)
            eta_squared = ss_between / ss_total if ss_total > 0 else 0
        else:
            # Epsilon-squared for Kruskal-Wallis
            n = len(analysis_df)
            eta_squared = (f_stat - len(groups) + 1) / (n - len(groups)) if n > len(groups) else 0

        result = {
            "test_name": test_name,
            "statistic": float(f_stat),
            "p_value": float(p_value),
            "significant": p_value < self.alpha,
            "effect_size": float(eta_squared),
            "effect_size_name": "η²" if use_parametric else "ε²",
            "n_total": len(analysis_df),
            "n_groups": len(groups),
            "group_sizes": group_sizes,
            "outliers_removed": n_removed,
            "normality_tests": normality_tests,
            "levene_test": {"statistic": levene_stat, "p_value": levene_p, "equal_variances": equal_variances},
            "parametric": use_parametric,
            "alpha": self.alpha
        }

        return result

    def run_pairwise_comparisons(
        self,
        df: pd.DataFrame,
        value_col: str,
        group_col: str = "Donor Status",
        method: str = "tukey"  # "tukey", "bonferroni", "holm"
    ) -> pd.DataFrame:
        """Run pairwise comparisons between groups."""

        analysis_df = df.dropna(subset=[value_col, group_col]).copy()

        if analysis_df.empty:
            return pd.DataFrame()

        groups = analysis_df[group_col].unique()

        if len(groups) < 2:
            return pd.DataFrame()

        if method == "tukey":
            # Tukey HSD
            tukey = pairwise_tukeyhsd(
                analysis_df[value_col],
                analysis_df[group_col],
                alpha=self.alpha
            )

            results = pd.DataFrame({
                "group1": tukey.groupsunique[tukey.pairindices[0]],
                "group2": tukey.groupsunique[tukey.pairindices[1]],
                "mean_diff": tukey.meandiffs,
                "p_adj": tukey.pvalues,
                "lower_ci": tukey.confint[:, 0],
                "upper_ci": tukey.confint[:, 1],
                "reject": tukey.reject
            })
        else:
            # Manual pairwise t-tests with correction
            comparisons = []
            p_values = []

            for i, g1 in enumerate(groups):
                for g2 in groups[i+1:]:
                    v1 = analysis_df[analysis_df[group_col] == g1][value_col]
                    v2 = analysis_df[analysis_df[group_col] == g2][value_col]

                    stat, p = ttest_ind(v1, v2)
                    comparisons.append({
                        "group1": g1,
                        "group2": g2,
                        "mean_diff": v1.mean() - v2.mean(),
                        "t_statistic": stat,
                        "p_value": p
                    })
                    p_values.append(p)

            results = pd.DataFrame(comparisons)

            # Apply correction
            if method == "bonferroni":
                results["p_adj"] = np.minimum(np.array(p_values) * len(p_values), 1.0)
            elif method == "holm":
                from statsmodels.stats.multitest import multipletests
                _, p_adj, _, _ = multipletests(p_values, method="holm")
                results["p_adj"] = p_adj
            else:
                # Benjamini-Hochberg
                from statsmodels.stats.multitest import multipletests
                _, p_adj, _, _ = multipletests(p_values, method="fdr_bh")
                results["p_adj"] = p_adj

            results["reject"] = results["p_adj"] < self.alpha

        return results

    def compute_descriptive_stats(
        self,
        df: pd.DataFrame,
        value_col: str,
        group_col: str = "Donor Status"
    ) -> pd.DataFrame:
        """Compute descriptive statistics by group."""

        analysis_df = df.dropna(subset=[value_col, group_col])

        stats_df = analysis_df.groupby(group_col)[value_col].agg([
            ("n", "count"),
            ("mean", "mean"),
            ("std", "std"),
            ("se", lambda x: x.std() / np.sqrt(len(x))),
            ("median", "median"),
            ("q25", lambda x: x.quantile(0.25)),
            ("q75", lambda x: x.quantile(0.75)),
            ("min", "min"),
            ("max", "max")
        ]).reset_index()

        return stats_df

    def compute_auc(
        self,
        df: pd.DataFrame,
        x_col: str,
        y_col: str,
        group_col: str = "Donor Status"
    ) -> pd.DataFrame:
        """Compute area under curve for each group using trapezoidal integration."""

        results = []

        for group in df[group_col].unique():
            group_df = df[df[group_col] == group].dropna(subset=[x_col, y_col])

            if len(group_df) < 2:
                continue

            # Sort by x value
            group_df = group_df.sort_values(x_col)

            # Trapezoidal integration
            auc = np.trapz(group_df[y_col], group_df[x_col])

            # Normalized AUC (per unit x)
            x_range = group_df[x_col].max() - group_df[x_col].min()
            auc_normalized = auc / x_range if x_range > 0 else 0

            results.append({
                "group": group,
                "auc": auc,
                "auc_normalized": auc_normalized,
                "n_points": len(group_df),
                "x_range": x_range
            })

        return pd.DataFrame(results)


def format_p_value(p: float) -> str:
    """Format p-value for display."""
    if p < 0.001:
        return "< 0.001"
    elif p < 0.01:
        return f"{p:.3f}"
    elif p < 0.05:
        return f"{p:.3f}"
    else:
        return f"{p:.3f}"


def significance_stars(p: float) -> str:
    """Return significance stars for p-value."""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def create_results_table(anova_result: Dict) -> pd.DataFrame:
    """Create formatted results table from ANOVA results."""

    if "error" in anova_result:
        return pd.DataFrame({"Error": [anova_result["error"]]})

    rows = [
        {"Metric": "Test", "Value": anova_result["test_name"]},
        {"Metric": "Statistic", "Value": f"{anova_result['statistic']:.4f}"},
        {"Metric": "P-value", "Value": format_p_value(anova_result["p_value"])},
        {"Metric": "Significance", "Value": significance_stars(anova_result["p_value"])},
        {"Metric": f"Effect Size ({anova_result['effect_size_name']})",
         "Value": f"{anova_result['effect_size']:.4f}"},
        {"Metric": "Sample Size", "Value": str(anova_result["n_total"])},
        {"Metric": "Groups", "Value": str(anova_result["n_groups"])},
        {"Metric": "Outliers Removed", "Value": str(anova_result["outliers_removed"])},
        {"Metric": "Parametric Test", "Value": "Yes" if anova_result["parametric"] else "No"},
    ]

    return pd.DataFrame(rows)
