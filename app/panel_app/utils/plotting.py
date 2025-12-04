"""
Plotting utilities for Islet Explorer Panel App.
Provides interactive plots using Plotly and HoloViews.
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import Dict, List, Optional, Tuple

# Donor status colors
DONOR_COLORS = {
    "ND": "#1f77b4",
    "Aab+": "#ff7f0e",
    "T1D": "#d62728"
}

# Default plot styling
PLOT_TEMPLATE = "plotly_white"


def create_scatter_plot(
    df: pd.DataFrame,
    x: str,
    y: str,
    color: str = "Donor Status",
    size: Optional[str] = None,
    title: str = "",
    x_title: str = "",
    y_title: str = "",
    log_x: bool = False,
    log_y: bool = False,
    point_size: float = 8,
    opacity: float = 0.7,
    color_map: Optional[Dict] = None,
    trendline: Optional[str] = None,  # "ols", "lowess", None
    height: int = 500
) -> go.Figure:
    """Create an interactive scatter plot with optional trendline."""

    if color_map is None:
        color_map = DONOR_COLORS

    # Filter out NaN values
    plot_df = df.dropna(subset=[x, y])

    if plot_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16)
        )
        fig.update_layout(height=height, template=PLOT_TEMPLATE, title=title)
        return fig

    # Create scatter plot
    fig = px.scatter(
        plot_df,
        x=x,
        y=y,
        color=color if color in plot_df.columns else None,
        size=size if size and size in plot_df.columns else None,
        color_discrete_map=color_map if color in plot_df.columns else None,
        opacity=opacity,
        trendline=trendline,
        log_x=log_x,
        log_y=log_y,
        template=PLOT_TEMPLATE,
        title=title,
        height=height
    )

    # Update marker size if no size column specified
    if not size:
        fig.update_traces(marker=dict(size=point_size))

    # Update axis labels
    fig.update_xaxes(title_text=x_title or x)
    fig.update_yaxes(title_text=y_title or y)

    # Update layout
    fig.update_layout(
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        margin=dict(l=60, r=20, t=60, b=60)
    )

    return fig


def create_distribution_plot(
    df: pd.DataFrame,
    x: str,
    group: str = "Donor Status",
    plot_type: str = "box",  # "box", "violin", "histogram"
    title: str = "",
    x_title: str = "",
    y_title: str = "",
    color_map: Optional[Dict] = None,
    show_points: bool = True,
    height: int = 400
) -> go.Figure:
    """Create distribution comparison plot (box, violin, or histogram)."""

    if color_map is None:
        color_map = DONOR_COLORS

    plot_df = df.dropna(subset=[x])

    if plot_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        fig.update_layout(height=height, template=PLOT_TEMPLATE)
        return fig

    if plot_type == "box":
        fig = px.box(
            plot_df,
            x=group if group in plot_df.columns else None,
            y=x,
            color=group if group in plot_df.columns else None,
            color_discrete_map=color_map,
            points="all" if show_points else False,
            template=PLOT_TEMPLATE,
            title=title,
            height=height
        )
    elif plot_type == "violin":
        fig = px.violin(
            plot_df,
            x=group if group in plot_df.columns else None,
            y=x,
            color=group if group in plot_df.columns else None,
            color_discrete_map=color_map,
            points="all" if show_points else False,
            box=True,
            template=PLOT_TEMPLATE,
            title=title,
            height=height
        )
    else:  # histogram
        fig = px.histogram(
            plot_df,
            x=x,
            color=group if group in plot_df.columns else None,
            color_discrete_map=color_map,
            barmode="overlay",
            opacity=0.7,
            template=PLOT_TEMPLATE,
            title=title,
            height=height
        )

    fig.update_xaxes(title_text=x_title or group)
    fig.update_yaxes(title_text=y_title or x)
    fig.update_layout(margin=dict(l=60, r=20, t=60, b=60))

    return fig


def create_heatmap(
    df: pd.DataFrame,
    x: str,
    y: str,
    z: str,
    title: str = "",
    colorscale: str = "Viridis",
    height: int = 400
) -> go.Figure:
    """Create a heatmap from pivoted data."""

    # Pivot data for heatmap
    try:
        pivot_df = df.pivot_table(values=z, index=y, columns=x, aggfunc="mean")
    except Exception:
        pivot_df = pd.DataFrame()

    if pivot_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No data available for heatmap",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        fig.update_layout(height=height, template=PLOT_TEMPLATE)
        return fig

    fig = px.imshow(
        pivot_df,
        color_continuous_scale=colorscale,
        title=title,
        template=PLOT_TEMPLATE,
        height=height,
        aspect="auto"
    )

    fig.update_layout(margin=dict(l=60, r=20, t=60, b=60))
    return fig


def create_umap_plot(
    df: pd.DataFrame,
    umap1: str = "UMAP1",
    umap2: str = "UMAP2",
    color: str = "Donor Status",
    title: str = "UMAP",
    point_size: float = 6,
    opacity: float = 0.8,
    color_map: Optional[Dict] = None,
    height: int = 450
) -> go.Figure:
    """Create UMAP scatter plot."""

    if color_map is None and color == "Donor Status":
        color_map = DONOR_COLORS

    plot_df = df.dropna(subset=[umap1, umap2])

    if plot_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No UMAP data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        fig.update_layout(height=height, template=PLOT_TEMPLATE, title=title)
        return fig

    # Check if color column is numeric or categorical
    is_numeric = pd.api.types.is_numeric_dtype(plot_df[color]) if color in plot_df.columns else False

    if is_numeric:
        fig = px.scatter(
            plot_df,
            x=umap1,
            y=umap2,
            color=color,
            color_continuous_scale="Viridis",
            opacity=opacity,
            template=PLOT_TEMPLATE,
            title=title,
            height=height
        )
    else:
        fig = px.scatter(
            plot_df,
            x=umap1,
            y=umap2,
            color=color if color in plot_df.columns else None,
            color_discrete_map=color_map if color in plot_df.columns else None,
            opacity=opacity,
            template=PLOT_TEMPLATE,
            title=title,
            height=height
        )

    fig.update_traces(marker=dict(size=point_size))
    fig.update_xaxes(title_text="UMAP 1", showgrid=False)
    fig.update_yaxes(title_text="UMAP 2", showgrid=False)
    fig.update_layout(
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
        margin=dict(l=40, r=20, t=60, b=40)
    )

    return fig


def create_bar_plot(
    df: pd.DataFrame,
    x: str,
    y: str,
    error_y: Optional[str] = None,
    color: str = "Donor Status",
    title: str = "",
    x_title: str = "",
    y_title: str = "",
    color_map: Optional[Dict] = None,
    height: int = 400
) -> go.Figure:
    """Create grouped bar plot with optional error bars."""

    if color_map is None:
        color_map = DONOR_COLORS

    fig = px.bar(
        df,
        x=x,
        y=y,
        error_y=error_y if error_y and error_y in df.columns else None,
        color=color if color in df.columns else None,
        color_discrete_map=color_map if color in df.columns else None,
        barmode="group",
        template=PLOT_TEMPLATE,
        title=title,
        height=height
    )

    fig.update_xaxes(title_text=x_title or x)
    fig.update_yaxes(title_text=y_title or y)
    fig.update_layout(margin=dict(l=60, r=20, t=60, b=60))

    return fig


def create_line_plot(
    df: pd.DataFrame,
    x: str,
    y: str,
    color: str = "Donor Status",
    error_band: Optional[str] = None,
    title: str = "",
    x_title: str = "",
    y_title: str = "",
    color_map: Optional[Dict] = None,
    height: int = 400
) -> go.Figure:
    """Create line plot with optional confidence bands."""

    if color_map is None:
        color_map = DONOR_COLORS

    fig = px.line(
        df,
        x=x,
        y=y,
        color=color if color in df.columns else None,
        color_discrete_map=color_map if color in df.columns else None,
        markers=True,
        template=PLOT_TEMPLATE,
        title=title,
        height=height
    )

    fig.update_xaxes(title_text=x_title or x)
    fig.update_yaxes(title_text=y_title or y)
    fig.update_layout(margin=dict(l=60, r=20, t=60, b=60))

    return fig


def create_islet_size_scatter(
    df: pd.DataFrame,
    x_col: str = "islet_diameter_um",
    y_col: str = "pos_frac",
    color_col: str = "Donor Status",
    diameter_range: Tuple[float, float] = (0, 500),
    bin_width: float = 25,
    show_individuals: bool = True,
    show_trend: bool = False,
    trend_type: str = "lowess",  # "lowess" or "ols"
    point_size: float = 8,
    opacity: float = 0.7,
    title: str = "Islet Size Distribution",
    height: int = 500
) -> go.Figure:
    """Create islet size distribution scatter plot with binning support."""

    plot_df = df.dropna(subset=[x_col, y_col]).copy()

    # Apply diameter range filter
    plot_df = plot_df[
        (plot_df[x_col] >= diameter_range[0]) &
        (plot_df[x_col] <= diameter_range[1])
    ]

    if plot_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No data in selected range",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        fig.update_layout(height=height, template=PLOT_TEMPLATE, title=title)
        return fig

    # Create bin column if binning
    if bin_width > 0:
        plot_df["diameter_bin"] = (plot_df[x_col] // bin_width) * bin_width + bin_width / 2

    fig = create_scatter_plot(
        plot_df,
        x=x_col,
        y=y_col,
        color=color_col,
        title=title,
        x_title="Islet Diameter (µm)",
        y_title=y_col,
        point_size=point_size,
        opacity=opacity,
        trendline=trend_type if show_trend else None,
        height=height
    )

    return fig
