"""
Islet Explorer Panel Application

A comprehensive spatial omics visualization tool for analyzing
pancreatic islet data from multiplex imaging studies.

Run with: panel serve app.py --address 0.0.0.0 --port 8080
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import panel as pn
import param
from tornado.web import StaticFileHandler

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from utils.data_loader import IsletDataLoader, REGION_TYPES
from utils.plotting import (
    create_scatter_plot, create_distribution_plot, create_umap_plot
)
from utils.statistics import StatisticalAnalyzer, create_results_table, format_p_value
from components.ai_assistant import AIAssistantPanel
from components.image_viewer import ImageViewerPanel

# Initialize Panel extensions
pn.extension(
    "plotly",
    "tabulator",
    sizing_mode="stretch_width",
    notifications=True
)

# Application version
APP_VERSION = "panel-v1.0.0"

# CSS Styling
CSS = """
:root {
    --primary-color: #2c5aa0;
    --primary-dark: #1e3a72;
    --primary-light: #66b3ff;
    --background: #f8f9fa;
    --card-bg: #ffffff;
    --text-primary: #1f2937;
    --text-secondary: #6b7280;
    --border-color: #e5e7eb;
}

.bk-root {
    font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
}

.header-gradient {
    background: linear-gradient(135deg, var(--primary-dark) 0%, var(--primary-color) 100%);
    color: white;
    padding: 15px 20px;
    border-radius: 8px;
    margin-bottom: 15px;
}

.card {
    background: var(--card-bg);
    border-radius: 8px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    padding: 15px;
    margin-bottom: 15px;
}

.sidebar {
    background: #f1f5f9;
    border-radius: 8px;
    padding: 15px;
}

.stat-value {
    font-size: 24px;
    font-weight: 600;
    color: var(--primary-color);
}

.stat-label {
    font-size: 12px;
    color: var(--text-secondary);
    text-transform: uppercase;
}
"""

pn.config.raw_css.append(CSS)


class IsletExplorerApp(param.Parameterized):
    """Main application class for Islet Explorer."""

    # Data parameters
    mode = param.Selector(default="Markers", objects=["Markers", "Targets", "Composition"])
    region = param.Selector(default="islet_core", objects=REGION_TYPES)
    metric = param.Selector(default=None, objects=[])

    # Filter parameters
    donor_nd = param.Boolean(default=True, label="ND")
    donor_aab = param.Boolean(default=True, label="Aab+")
    donor_t1d = param.Boolean(default=True, label="T1D")

    # Autoantibody filters
    filter_gada = param.Boolean(default=True, label="GADA")
    filter_ia2a = param.Boolean(default=True, label="IA2A")
    filter_znt8a = param.Boolean(default=True, label="ZnT8A")
    filter_iaa = param.Boolean(default=True, label="IAA")
    filter_miaa = param.Boolean(default=True, label="mIAA")

    # Visualization parameters
    stat_type = param.Selector(
        default="mean_se",
        objects=["mean_se", "mean_sd", "median_iqr"],
        label="Error Bars"
    )
    point_size = param.Number(default=8, bounds=(2, 20), label="Point Size")
    opacity = param.Number(default=0.7, bounds=(0.1, 1.0), label="Opacity")
    log_x = param.Boolean(default=False, label="Log X-axis")
    log_y = param.Boolean(default=False, label="Log Y-axis")
    show_trend = param.Boolean(default=False, label="Show Trend Line")
    trend_type = param.Selector(default="lowess", objects=["lowess", "ols"], label="Trend Type")

    # Islet size parameters
    diameter_min = param.Number(default=0, bounds=(0, 500), label="Min Diameter (µm)")
    diameter_max = param.Number(default=500, bounds=(0, 500), label="Max Diameter (µm)")
    bin_width = param.Number(default=25, bounds=(5, 100), label="Bin Width")

    # Statistics parameters
    alpha = param.Selector(default=0.05, objects=[0.05, 0.01, 0.001], label="α Level")
    remove_outliers = param.Boolean(default=False, label="Remove Outliers (>3 SD)")

    def __init__(self, data_path=None, **params):
        super().__init__(**params)

        # Initialize data loader
        self.data_loader = IsletDataLoader(data_path)
        self._data_loaded = False
        self._load_data()

        # Initialize statistics analyzer
        self.stats_analyzer = StatisticalAnalyzer(alpha=self.alpha)

        # Create AI assistant
        self.ai_assistant = AIAssistantPanel(
            context="islet_analysis",
            on_plot_context=self._get_plot_context
        )

        # Create image viewer
        self.image_viewer = ImageViewerPanel()

        # Setup watchers
        self._setup_watchers()

    def _load_data(self):
        """Load data from Excel file."""
        try:
            self.data_loader.load_all()
            self._data_loaded = True
            self._update_metric_options()
            pn.state.notifications.success("Data loaded successfully", duration=3000)
        except FileNotFoundError as e:
            self._data_loaded = False
            pn.state.notifications.error(str(e), duration=5000)
        except Exception as e:
            self._data_loaded = False
            pn.state.notifications.error(f"Error loading data: {e}", duration=5000)

    def _update_metric_options(self):
        """Update metric options based on selected mode."""
        if not self._data_loaded:
            return

        if self.mode == "Markers":
            metrics = self.data_loader.get_markers()
            if not metrics:
                metrics = ["pos_frac"]
        elif self.mode == "Targets":
            df = self.data_loader.targets
            metrics = [c for c in ["area_density", "count", "area_um2"] if c in df.columns]
        else:  # Composition
            metrics = ["Ins_any_frac", "Glu_any_frac", "Stt_any_frac"]

        self.param.metric.objects = metrics
        if metrics:
            self.metric = metrics[0]

    def _setup_watchers(self):
        """Setup parameter watchers."""
        self.param.watch(lambda e: self._update_metric_options(), ["mode"])

    def _get_selected_donors(self):
        """Get list of selected donor statuses."""
        donors = []
        if self.donor_nd:
            donors.append("ND")
        if self.donor_aab:
            donors.append("Aab+")
        if self.donor_t1d:
            donors.append("T1D")
        return donors

    def _get_filtered_data(self):
        """Get filtered data based on current selections."""
        if not self._data_loaded:
            return pd.DataFrame()

        # Get dataframe based on mode
        if self.mode == "Markers":
            df = self.data_loader.markers.copy()
        elif self.mode == "Targets":
            df = self.data_loader.targets.copy()
        else:
            df = self.data_loader.composition.copy()

        # Filter by donor status
        donors = self._get_selected_donors()
        if donors and "Donor Status" in df.columns:
            df = df[df["Donor Status"].isin(donors)]

        # Filter by region
        if "region" in df.columns:
            df = df[df["region"].str.contains(self.region, case=False, na=False)]

        # Filter by marker (for Markers mode)
        if self.mode == "Markers" and self.metric and "marker" in df.columns:
            df = df[df["marker"] == self.metric]

        return df

    def _get_plot_context(self):
        """Get current plot context for AI assistant."""
        df = self._get_filtered_data()
        if df.empty:
            return None

        context = f"""
Current analysis:
- Mode: {self.mode}
- Region: {self.region}
- Metric: {self.metric}
- Donor groups: {', '.join(self._get_selected_donors())}
- Sample size: {len(df)} observations
- Groups: {df.get('Donor Status', pd.Series()).value_counts().to_dict()}
"""
        return context

    # =========================================================================
    # UI Components
    # =========================================================================

    @pn.depends("mode", "region", "metric")
    def sidebar_filters(self):
        """Create sidebar filter controls."""
        mode_select = pn.widgets.RadioButtonGroup(
            name="Mode",
            options=["Markers", "Targets", "Composition"],
            value=self.mode,
            button_type="primary"
        )
        mode_select.link(self, value="mode")

        region_select = pn.widgets.Select(
            name="Region",
            options=REGION_TYPES,
            value=self.region,
            width=200
        )
        region_select.link(self, value="region")

        metric_select = pn.widgets.Select(
            name="Metric",
            options=self.param.metric.objects,
            value=self.metric,
            width=200
        )
        metric_select.link(self, value="metric")

        # Donor status filters
        donor_filters = pn.Column(
            pn.pane.Markdown("**Donor Status**"),
            pn.widgets.Checkbox.from_param(self.param.donor_nd),
            pn.widgets.Checkbox.from_param(self.param.donor_aab),
            pn.widgets.Checkbox.from_param(self.param.donor_t1d),
        )

        # Autoantibody filters (only show for Aab+ selected)
        aab_filters = pn.Column(
            pn.pane.Markdown("**Autoantibodies (Aab+ only)**"),
            pn.widgets.Checkbox.from_param(self.param.filter_gada),
            pn.widgets.Checkbox.from_param(self.param.filter_ia2a),
            pn.widgets.Checkbox.from_param(self.param.filter_znt8a),
            pn.widgets.Checkbox.from_param(self.param.filter_iaa),
            pn.widgets.Checkbox.from_param(self.param.filter_miaa),
            visible=self.donor_aab
        )

        # Visualization options
        viz_options = pn.Column(
            pn.pane.Markdown("**Visualization**"),
            pn.widgets.Select.from_param(self.param.stat_type, width=150),
            pn.widgets.IntSlider.from_param(self.param.point_size, width=150),
            pn.widgets.FloatSlider.from_param(self.param.opacity, width=150),
            pn.widgets.Checkbox.from_param(self.param.log_x),
            pn.widgets.Checkbox.from_param(self.param.log_y),
            pn.widgets.Checkbox.from_param(self.param.show_trend),
        )

        return pn.Column(
            pn.pane.Markdown("## Filters", styles={"margin-bottom": "10px"}),
            mode_select,
            pn.Spacer(height=10),
            region_select,
            metric_select,
            pn.Spacer(height=15),
            donor_filters,
            aab_filters,
            pn.Spacer(height=15),
            viz_options,
            css_classes=["sidebar"],
            width=250
        )

    @pn.depends(
        "mode", "region", "metric", "donor_nd", "donor_aab", "donor_t1d",
        "point_size", "opacity", "log_x", "log_y", "show_trend", "trend_type"
    )
    def main_scatter_plot(self):
        """Create main scatter plot."""
        df = self._get_filtered_data()

        if df.empty:
            return pn.pane.Markdown("No data available. Check filters or data file.")

        # Determine x and y columns
        x_col = "islet_diameter_um" if "islet_diameter_um" in df.columns else df.columns[0]
        y_col = self.metric if self.metric in df.columns else "pos_frac"

        if y_col not in df.columns:
            y_col = [c for c in df.columns if "frac" in c.lower() or "density" in c.lower()]
            y_col = y_col[0] if y_col else df.select_dtypes(include=[np.number]).columns[0]

        fig = create_scatter_plot(
            df,
            x=x_col,
            y=y_col,
            color="Donor Status",
            title=f"{self.metric or y_col} by Islet Size",
            x_title="Islet Diameter (µm)",
            y_title=self.metric or y_col,
            log_x=self.log_x,
            log_y=self.log_y,
            point_size=self.point_size,
            opacity=self.opacity,
            trendline=self.trend_type if self.show_trend else None,
            height=450
        )

        return pn.pane.Plotly(fig, sizing_mode="stretch_width")

    @pn.depends("mode", "region", "metric", "donor_nd", "donor_aab", "donor_t1d", "stat_type")
    def distribution_plot(self):
        """Create distribution comparison plot."""
        df = self._get_filtered_data()

        if df.empty:
            return pn.pane.Markdown("No data available.")

        y_col = self.metric if self.metric in df.columns else "pos_frac"
        if y_col not in df.columns:
            y_col = df.select_dtypes(include=[np.number]).columns[0] if len(df.select_dtypes(include=[np.number]).columns) > 0 else None

        if y_col is None:
            return pn.pane.Markdown("No numeric columns for distribution plot.")

        fig = create_distribution_plot(
            df,
            x=y_col,
            group="Donor Status",
            plot_type="box",
            title="Distribution by Donor Status",
            show_points=True,
            height=400
        )

        return pn.pane.Plotly(fig, sizing_mode="stretch_width")

    @pn.depends("mode", "region", "metric", "donor_nd", "donor_aab", "donor_t1d")
    def summary_table(self):
        """Create summary statistics table."""
        df = self._get_filtered_data()

        if df.empty:
            return pn.pane.Markdown("No data available.")

        y_col = self.metric if self.metric in df.columns else "pos_frac"
        if y_col not in df.columns:
            return pn.pane.Markdown("Metric not found in data.")

        stats_df = df.groupby("Donor Status")[y_col].agg([
            ("N", "count"),
            ("Mean", "mean"),
            ("SD", "std"),
            ("Median", "median"),
            ("Min", "min"),
            ("Max", "max")
        ]).round(4).reset_index()

        return pn.widgets.Tabulator(
            stats_df,
            sizing_mode="stretch_width",
            height=150,
            show_index=False
        )

    def download_summary(self):
        """Create download button for summary data."""
        def get_csv():
            df = self._get_filtered_data()
            return df.to_csv(index=False)

        return pn.widgets.FileDownload(
            callback=get_csv,
            filename="islet_summary.csv",
            button_type="primary",
            label="Download Data"
        )

    # =========================================================================
    # Trajectory Tab
    # =========================================================================

    @pn.depends("donor_nd", "donor_aab", "donor_t1d", "point_size", "opacity")
    def trajectory_scatter(self):
        """Create trajectory scatterplot."""
        df = self._get_filtered_data()

        if df.empty or "islet_diameter_um" not in df.columns:
            return pn.pane.Markdown("No trajectory data available.")

        y_col = self.metric if self.metric in df.columns else "pos_frac"
        if y_col not in df.columns:
            y_col = df.select_dtypes(include=[np.number]).columns[0] if len(df.select_dtypes(include=[np.number]).columns) > 0 else None

        if y_col is None:
            return pn.pane.Markdown("No numeric data for trajectory.")

        fig = create_scatter_plot(
            df,
            x="islet_diameter_um",
            y=y_col,
            color="Donor Status",
            title="Trajectory Analysis",
            point_size=self.point_size,
            opacity=self.opacity,
            height=400
        )

        return pn.pane.Plotly(fig, sizing_mode="stretch_width")

    @pn.depends("donor_nd", "donor_aab", "donor_t1d")
    def umap_by_status(self):
        """Create UMAP colored by donor status."""
        # For demo, create synthetic UMAP data
        df = self._get_filtered_data()

        if df.empty:
            return pn.pane.Markdown("No UMAP data available.")

        # Check if UMAP columns exist, otherwise create placeholder
        if "UMAP1" not in df.columns:
            # Create placeholder message
            return pn.Column(
                pn.pane.Markdown("### UMAP by Donor Status"),
                pn.pane.HTML("""
                <div style="
                    display: flex;
                    align-items: center;
                    justify-content: center;
                    height: 350px;
                    background: #f3f4f6;
                    border-radius: 8px;
                    color: #6b7280;
                ">
                    <div style="text-align: center;">
                        <p>UMAP coordinates not found in data.</p>
                        <p style="font-size: 12px;">Run trajectory analysis to generate UMAP.</p>
                    </div>
                </div>
                """)
            )

        fig = create_umap_plot(
            df,
            color="Donor Status",
            title="UMAP by Donor Status",
            point_size=self.point_size,
            opacity=self.opacity,
            height=400
        )

        return pn.pane.Plotly(fig, sizing_mode="stretch_width")

    @pn.depends("metric", "donor_nd", "donor_aab", "donor_t1d")
    def umap_by_feature(self):
        """Create UMAP colored by selected feature."""
        df = self._get_filtered_data()

        if df.empty or "UMAP1" not in df.columns:
            return pn.Column(
                pn.pane.Markdown("### UMAP by Feature"),
                pn.pane.HTML("""
                <div style="
                    display: flex;
                    align-items: center;
                    justify-content: center;
                    height: 350px;
                    background: #f3f4f6;
                    border-radius: 8px;
                    color: #6b7280;
                ">
                    <p>UMAP coordinates not found in data.</p>
                </div>
                """)
            )

        color_col = self.metric if self.metric in df.columns else "pos_frac"

        fig = create_umap_plot(
            df,
            color=color_col,
            title=f"UMAP by {color_col}",
            point_size=self.point_size,
            opacity=self.opacity,
            height=400
        )

        return pn.pane.Plotly(fig, sizing_mode="stretch_width")

    @pn.depends("mode", "region", "metric")
    def trajectory_heatmap(self):
        """Create trajectory heatmap."""
        df = self._get_filtered_data()

        if df.empty:
            return pn.pane.Markdown("No data for heatmap.")

        # For now, return placeholder
        return pn.pane.HTML("""
        <div style="
            display: flex;
            align-items: center;
            justify-content: center;
            height: 200px;
            background: linear-gradient(90deg, #3b82f6, #8b5cf6, #ec4899);
            border-radius: 8px;
            color: white;
        ">
            <p>Feature expression heatmap will appear here when trajectory data is loaded.</p>
        </div>
        """)

    # =========================================================================
    # Statistics Tab
    # =========================================================================

    @pn.depends("mode", "region", "metric", "alpha", "remove_outliers",
                "donor_nd", "donor_aab", "donor_t1d")
    def anova_results(self):
        """Run ANOVA and display results."""
        df = self._get_filtered_data()

        if df.empty:
            return pn.pane.Markdown("No data for statistical analysis.")

        y_col = self.metric if self.metric in df.columns else "pos_frac"
        if y_col not in df.columns:
            return pn.pane.Markdown("Selected metric not in data.")

        self.stats_analyzer.alpha = self.alpha
        result = self.stats_analyzer.run_anova(
            df, y_col,
            group_col="Donor Status",
            remove_outliers=self.remove_outliers
        )

        results_df = create_results_table(result)

        return pn.Column(
            pn.pane.Markdown("### ANOVA Results"),
            pn.widgets.Tabulator(
                results_df,
                sizing_mode="stretch_width",
                height=250,
                show_index=False
            )
        )

    @pn.depends("mode", "region", "metric", "alpha", "donor_nd", "donor_aab", "donor_t1d")
    def pairwise_results(self):
        """Run pairwise comparisons and display results."""
        df = self._get_filtered_data()

        if df.empty:
            return pn.pane.Markdown("No data for pairwise comparisons.")

        y_col = self.metric if self.metric in df.columns else "pos_frac"
        if y_col not in df.columns:
            return pn.pane.Markdown("Selected metric not in data.")

        pairwise_df = self.stats_analyzer.run_pairwise_comparisons(
            df, y_col, group_col="Donor Status"
        )

        if pairwise_df.empty:
            return pn.pane.Markdown("Not enough groups for pairwise comparisons.")

        # Format p-values
        pairwise_df["p_adj"] = pairwise_df["p_adj"].apply(format_p_value)
        pairwise_df["Significant"] = pairwise_df["reject"].map({True: "Yes", False: "No"})

        display_df = pairwise_df[["group1", "group2", "mean_diff", "p_adj", "Significant"]].copy()
        display_df.columns = ["Group 1", "Group 2", "Mean Diff", "Adj. P-value", "Significant"]
        display_df["Mean Diff"] = display_df["Mean Diff"].round(4)

        return pn.Column(
            pn.pane.Markdown("### Pairwise Comparisons (Tukey HSD)"),
            pn.widgets.Tabulator(
                display_df,
                sizing_mode="stretch_width",
                height=150,
                show_index=False
            )
        )

    @pn.depends("mode", "region", "metric", "donor_nd", "donor_aab", "donor_t1d")
    def stats_distribution_plot(self):
        """Create distribution plot for statistics tab."""
        df = self._get_filtered_data()

        if df.empty:
            return pn.pane.Markdown("No data available.")

        y_col = self.metric if self.metric in df.columns else "pos_frac"
        if y_col not in df.columns:
            return pn.pane.Markdown("Metric not found.")

        fig = create_distribution_plot(
            df, x=y_col, group="Donor Status",
            plot_type="violin",
            title="Distribution by Donor Status",
            show_points=True,
            height=350
        )

        return pn.pane.Plotly(fig, sizing_mode="stretch_width")

    # =========================================================================
    # Layout Assembly
    # =========================================================================

    def create_header(self):
        """Create application header."""
        return pn.Row(
            pn.pane.Markdown(
                "# Islet Explorer",
                styles={"margin": "0", "color": "white"}
            ),
            pn.Spacer(),
            pn.pane.Markdown(
                f"v{APP_VERSION} | Spatial Omics Analysis Platform",
                styles={"margin": "0", "color": "rgba(255,255,255,0.8)", "font-size": "14px"}
            ),
            css_classes=["header-gradient"],
            sizing_mode="stretch_width"
        )

    def create_plot_tab(self):
        """Create the main Plot tab."""
        # Left sidebar with filters
        sidebar = pn.Column(
            self.sidebar_filters,
            self.download_summary,
            width=260,
            sizing_mode="stretch_height"
        )

        # Main content area
        scatter_card = pn.Column(
            pn.pane.Markdown("### Islet Size Distribution"),
            self.main_scatter_plot,
            css_classes=["card"]
        )

        dist_card = pn.Column(
            pn.pane.Markdown("### Distribution Comparison"),
            self.distribution_plot,
            css_classes=["card"]
        )

        summary_card = pn.Column(
            pn.pane.Markdown("### Summary Statistics"),
            self.summary_table,
            css_classes=["card"]
        )

        # AI Assistant
        ai_panel = pn.Column(
            self.ai_assistant.get_panel(),
            sizing_mode="stretch_both",
            min_width=300
        )

        # Arrange in columns
        main_content = pn.Row(
            pn.Column(scatter_card, summary_card, sizing_mode="stretch_both"),
            pn.Column(dist_card, sizing_mode="stretch_both"),
            ai_panel,
            sizing_mode="stretch_both"
        )

        return pn.Row(
            sidebar,
            main_content,
            sizing_mode="stretch_both"
        )

    def create_trajectory_tab(self):
        """Create the Trajectory tab."""
        # Left panel - scatter and heatmap
        left_panel = pn.Column(
            pn.pane.Markdown("### Feature Trajectory"),
            self.trajectory_scatter,
            pn.Spacer(height=15),
            self.trajectory_heatmap,
            sizing_mode="stretch_both"
        )

        # Right panel - UMAPs
        right_panel = pn.Column(
            self.umap_by_status,
            pn.Spacer(height=15),
            self.umap_by_feature,
            sizing_mode="stretch_both"
        )

        return pn.Row(
            left_panel,
            right_panel,
            sizing_mode="stretch_both"
        )

    def create_viewer_tab(self):
        """Create the Image Viewer tab."""
        return self.image_viewer.get_panel()

    def create_statistics_tab(self):
        """Create the Statistics tab."""
        # Controls
        controls = pn.Row(
            pn.widgets.Select.from_param(self.param.alpha, width=100),
            pn.widgets.Checkbox.from_param(self.param.remove_outliers),
            sizing_mode="stretch_width"
        )

        # Results section
        results = pn.Row(
            pn.Column(
                self.anova_results,
                self.pairwise_results,
                sizing_mode="stretch_both"
            ),
            pn.Column(
                self.stats_distribution_plot,
                sizing_mode="stretch_both"
            ),
            sizing_mode="stretch_both"
        )

        # Method information
        info_text = """
### Statistical Methods

**ANOVA / Kruskal-Wallis:**
- Parametric one-way ANOVA is used when data is normally distributed and variances are equal
- Non-parametric Kruskal-Wallis H-test is used otherwise
- Normality assessed via Shapiro-Wilk test; homogeneity via Levene's test

**Pairwise Comparisons:**
- Tukey HSD for parametric data
- Controls family-wise error rate

**Effect Size:**
- η² (eta-squared) for ANOVA: small=0.01, medium=0.06, large=0.14
- ε² (epsilon-squared) for Kruskal-Wallis

**Outlier Detection:**
- Points beyond 3 standard deviations from mean
"""

        info = pn.pane.Markdown(info_text, styles={"font-size": "13px"})

        return pn.Column(
            pn.pane.Markdown("## Statistical Analysis"),
            controls,
            pn.Spacer(height=15),
            results,
            pn.Spacer(height=15),
            pn.Accordion(("Method Information", info), active=[]),
            sizing_mode="stretch_both"
        )

    def create_app(self):
        """Create the complete application."""
        header = self.create_header()

        tabs = pn.Tabs(
            ("Plot", self.create_plot_tab()),
            ("Trajectory", self.create_trajectory_tab()),
            ("Viewer", self.create_viewer_tab()),
            ("Statistics", self.create_statistics_tab()),
            dynamic=True,
            sizing_mode="stretch_both"
        )

        return pn.Column(
            header,
            tabs,
            sizing_mode="stretch_both",
            styles={"background": "#f8f9fa", "min-height": "100vh"}
        )


def create_app():
    """Factory function to create the app."""
    # Check for data path environment variable
    data_path = os.environ.get("ISLET_DATA_PATH")

    app = IsletExplorerApp(data_path=data_path)
    return app.create_app()


# Build static directories for image serving
def get_static_dirs():
    """Get static directory mappings."""
    static_dirs = {}
    # Check for images directory
    env_root = os.environ.get("LOCAL_IMAGE_ROOT", "")
    if env_root and Path(env_root).exists():
        static_dirs["images"] = env_root
    else:
        local_images = Path(__file__).parent / "local_images"
        if local_images.exists():
            static_dirs["images"] = str(local_images)

    # Check for avivator directory
    avivator_dir = Path(__file__).parent / "assets" / "avivator"
    if avivator_dir.exists():
        static_dirs["avivator"] = str(avivator_dir)

    return static_dirs

STATIC_DIRS = get_static_dirs()


# Custom static handler with CORS support for AVIVATOR
class CORSStaticHandler(StaticFileHandler):
    """Static file handler with CORS headers for AVIVATOR."""

    def set_default_headers(self):
        self.set_header("Access-Control-Allow-Origin", "*")
        self.set_header("Access-Control-Allow-Methods", "GET, OPTIONS")
        self.set_header("Access-Control-Allow-Headers", "Range, Accept-Encoding")
        self.set_header("Access-Control-Expose-Headers", "Content-Length, Content-Range, Accept-Ranges")
        self.set_header("Accept-Ranges", "bytes")

    def options(self, *args):
        self.set_status(204)
        self.finish()


def modify_doc(doc):
    """Add CORS routes to the Bokeh document."""
    pass


# For panel serve - add CORS routes
if pn.state.served:
    create_app().servable(title="Islet Explorer")

    # Add static handlers for images and avivator
    def add_static_routes(server):
        image_dir = STATIC_DIRS.get("images")
        avivator_dir = STATIC_DIRS.get("avivator")
        handlers = []
        if image_dir:
            handlers.append((r"/images/(.*)", CORSStaticHandler, {"path": image_dir}))
            print(f"[CORS] Added handler for /images/ -> {image_dir}")
        if avivator_dir:
            handlers.append((r"/avivator/(.*)", StaticFileHandler, {"path": avivator_dir}))
            print(f"[Static] Added handler for /avivator/ -> {avivator_dir}")
        if handlers:
            server._tornado.add_handlers(r".*", handlers)

    pn.state.on_session_created(lambda ctx: None)  # Placeholder

# For direct execution
if __name__ == "__main__":

    app = create_app()

    # Custom server setup with CORS for images and local AVIVATOR
    extra_patterns = []

    # Serve images with CORS
    image_dir = STATIC_DIRS.get("images")
    if image_dir:
        extra_patterns.append((r"/images/(.*)", CORSStaticHandler, {"path": image_dir}))
        print(f"[CORS] Serving images from: {image_dir}")

    # Serve local AVIVATOR (same origin = no CORS issues)
    avivator_dir = STATIC_DIRS.get("avivator")
    if avivator_dir:
        extra_patterns.append((r"/avivator/(.*)", StaticFileHandler, {"path": avivator_dir}))
        print(f"[Static] Serving AVIVATOR from: {avivator_dir}")

    pn.serve(
        app,
        port=8080,
        show=False,
        title="Islet Explorer",
        extra_patterns=extra_patterns,
        allow_websocket_origin=["*"]
    )
