
import io
import os
import re
import math
import zipfile
import tempfile
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from PIL import Image
import streamlit as st

st.set_page_config(page_title="Insulitis Explorer", layout="wide")

# Header image
try:
    with open("assets/header.png", "rb") as _fh:
        st.image(_fh.read(), use_container_width=True)
except Exception:
    pass

st.title("Insulitis Explorer")
st.markdown("""
Explore across donors:
- Islet size distributions with counts or cumulative area
- Composition overlays (INS/GCG/SST) aligned to islet size
- Peri-Islet composition and targets vs islet size
- Binned heatmaps by donor and islet size
""")

# ------------------------- Helpers -------------------------

@st.cache_data
def load_excel(file_bytes: bytes, sheet_name: Optional[str] = None) -> pd.DataFrame:
    if sheet_name is not None:
        return pd.read_excel(io.BytesIO(file_bytes), sheet_name=sheet_name)
    xls = pd.ExcelFile(io.BytesIO(file_bytes))
    return pd.read_excel(io.BytesIO(file_bytes), sheet_name=xls.sheet_names[0])

def coalesce_area_columns(df: pd.DataFrame, region: str) -> Tuple[np.ndarray, List[str]]:
    exp1 = f"islet_area_{region}_um2"
    exp2 = f"{region}__islet_area_{region}_um2"
    if exp1 in df.columns:
        return pd.to_numeric(df[exp1], errors="coerce").to_numpy(), [exp1]
    if exp2 in df.columns:
        return pd.to_numeric(df[exp2], errors="coerce").to_numpy(), [exp2]
    cands = [f"{region}__region_Nerve", f"{region}__region_Capillary", f"{region}__region_Lymphatic"]
    present = [c for c in cands if c in df.columns]
    if present:
        M = df[present].apply(pd.to_numeric, errors="coerce")
        area = M.bfill(axis=1).iloc[:,0].to_numpy()  # first non-NA across present columns
        return area, present
    return np.full(len(df), np.nan), []

def compute_diameter_um(area_um2: np.ndarray) -> np.ndarray:
    area_um2 = np.where(np.isfinite(area_um2) & (area_um2 > 0), area_um2, np.nan)
    return 2.0 * np.sqrt(area_um2 / math.pi)

def parse_mid_from_cutlabel(lbl: str) -> float:
    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(lbl) or "")
    if len(nums) >= 2:
        return (float(nums[0]) + float(nums[1])) / 2.0
    return float("nan")

def compute_log_bin_edges(log_vals: np.ndarray, binwidth: float) -> np.ndarray:
    finite = np.isfinite(log_vals)
    if not finite.any():
        return np.array([])
    lo = math.floor(np.nanmin(log_vals[finite]))
    hi = math.ceil(np.nanmax(log_vals[finite]))
    return np.arange(lo, hi + binwidth, binwidth)

def assign_bins(log_vals: np.ndarray, edges: np.ndarray) -> np.ndarray:
    if edges.size < 2:
        return np.full_like(log_vals, fill_value=-1, dtype=int)
    digitized = np.digitize(log_vals, edges) - 1
    valid = (digitized >= 0) & (digitized < (len(edges) - 1))
    out = np.where(valid, digitized, -1)
    return out

def format_um(x: float) -> str:
    try:
        if x >= 100:
            return f"{x:.0f}"
        elif x >= 10:
            return f"{x:.1f}"
        else:
            return f"{x:.2f}"
    except Exception:
        return str(x)

def normalize_remote_url(url: str) -> str:
    """Attempt to coerce OneDrive/SharePoint share URLs into direct download/view links.
    This is best-effort; manifests should ideally contain directly fetchable image URLs.
    """
    if not isinstance(url, str):
        return url
    u = url.strip()
    if ("1drv.ms" in u) or ("onedrive.live.com" in u) or ("sharepoint.com" in u):
        # Append download=1 if not present
        if "download=1" not in u:
            sep = "&" if ("?" in u) else "?"
            u = f"{u}{sep}download=1"
    return u

# Marker configuration (expected columns like: core__pct_INS, band__dens_CD163, etc.)
ENDO_MARKERS = ["INS", "GCG", "SST"]
IMMUNE_MARKERS = ["CD3E", "CD68", "CD163", "LGALS3"]

PREFERRED_MEASURES = ["pct", "frac", "dens", "count"]
MEASURE_SYNONYMS: Dict[str, List[str]] = {
    "pct": ["pct", "percent", "percentage", "perc", "pctg"],
    "frac": ["frac", "fraction", "prop", "proportion"],
    "dens": ["dens", "density", "per_area", "cells_per_um2", "cells_per_mm2"],
    "count": ["count", "n", "num", "cells"]
}

def _normalize(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())

def find_marker_columns_for(df: pd.DataFrame, region: str, markers: List[str]) -> Dict[str, Dict[str, str]]:
    """Return mapping: measure -> {marker -> column_name} using flexible, case-insensitive matching.
    Accepts columns containing region, a measure synonym, and the marker (ignoring punctuation/case).
    """
    out: Dict[str, Dict[str, str]] = {}
    region_n = _normalize(region)
    # Precompute normalized columns
    cols = list(df.columns)
    cols_norm = [(_normalize(c), c) for c in cols]
    # Build reverse lookup for marker normalized names
    marker_norms = {m: _normalize(m) for m in markers}
    # Build measure synonym lookup
    syn2canon = {}
    for canon, syns in MEASURE_SYNONYMS.items():
        for s in syns:
            syn2canon[_normalize(s)] = canon

    for norm, original in cols_norm:
        # Check region token present
        if region_n not in norm:
            continue
        # Detect canonical measure via any synonym token
        canon_measure = None
        for syn_norm, canon in syn2canon.items():
            if syn_norm and syn_norm in norm:
                canon_measure = canon
                break
        if canon_measure is None:
            continue
        # Detect marker
        for m, mn in marker_norms.items():
            if mn and mn in norm:
                out.setdefault(canon_measure, {})
                # Prefer first seen column per (measure, marker)
                out[canon_measure].setdefault(m, original)
    return out

def compute_marker_bin_means(df: pd.DataFrame, donor_col: str, diam_um: np.ndarray, binwidth: float,
                             marker_col: str, max_um: float = 300.0) -> pd.DataFrame:
    vals = pd.to_numeric(df[marker_col], errors='coerce')
    # Convert to percentage if it's a fraction-like column
    colnorm = _normalize(marker_col)
    if '__frac_' in colnorm or colnorm.endswith('frac'):
        vals = vals * 100.0
    finite = np.isfinite(diam_um) & np.isfinite(vals.to_numpy())
    d = diam_um[finite]
    v = vals.to_numpy()[finite]
    g = df.loc[finite, donor_col].astype(str).to_numpy()
    if d.size == 0:
        return pd.DataFrame(columns=['donor_status','diam_mid_um','mean_pct'])
    # Build bins with 300+ bucket
    logs = np.log10(d)
    lo = math.floor(np.nanmin(logs))
    max_log = math.log10(max_um)
    hi_needed = max(max_log, math.ceil(np.nanmax(logs)))
    edges = np.arange(lo, hi_needed + binwidth, binwidth)
    bin_idx = np.digitize(logs, edges) - 1
    plus_bin = len(edges) - 1
    over = (logs >= edges[-1])
    bin_idx[over] = plus_bin
    valid = (bin_idx >= 0) & (bin_idx <= plus_bin)
    dfb = pd.DataFrame({
        'donor_status': g[valid],
        'bin': bin_idx[valid],
        'val': v[valid],
    })
    agg = dfb.groupby(['donor_status','bin'])['val'].mean().reset_index(name='mean_pct')
    mids = (edges[:-1] + edges[1:]) / 2.0
    agg['diam_mid_um'] = agg['bin'].map(lambda i: (10**mids[int(i)] if i < len(mids) else max_um))
    return agg

def try_merge_targets_from_zip(df: pd.DataFrame, zip_path: str = "data/results.zip") -> pd.DataFrame:
    """If Peri-Islet target columns are missing, attempt to load them from results.zip targets TSVs and merge by islet_key.
    Produces columns like band__dens_Capillary / band__count_Capillary, etc.
    """
    need = [
        'band__dens_Capillary','band__dens_Nerve','band__dens_Lymphatic',
        'band__count_Capillary','band__count_Nerve','band__count_Lymphatic'
    ]
    if all(c in df.columns for c in need):
        return df
    try:
        import zipfile
        import pandas as pd
        import numpy as np
        with zipfile.ZipFile(zip_path) as zf:
            names = [n for n in zf.namelist() if n.lower().endswith('.tsv') and ('targets' in n.lower())]
            if not names:
                return df
            rows = []
            for name in names:
                with zf.open(name) as fh:
                    tdf = pd.read_csv(fh, sep='\t')
                    if {'region','class'}.issubset(set(tdf.columns)):
                        # extract islet key & region label
                        def _parse(s):
                            m = re.match(r'^(Islet_\d+)_(core|band|union)$', str(s))
                            return (m.group(1), m.group(2)) if m else (None, None)
                        tdf['islet_key'] = tdf['region'].apply(lambda s: _parse(s)[0])
                        tdf['region_label'] = tdf['region'].apply(lambda s: _parse(s)[1])
                        tdf = tdf[tdf['region_label'] == 'band']
                        if tdf.empty:
                            continue
                        # pivot dens & count
                        out = pd.DataFrame()
                        if 'area_density' in tdf.columns:
                            p = tdf.pivot_table(index='islet_key', columns='class', values='area_density', aggfunc='first')
                            p.columns = [f"band__dens_{c}" for c in p.columns]
                            out = pd.concat([out, p], axis=1)
                        if 'count' in tdf.columns:
                            p = tdf.pivot_table(index='islet_key', columns='class', values='count', aggfunc='first')
                            p.columns = [f"band__count_{c}" for c in p.columns]
                            out = pd.concat([out, p], axis=1)
                        if not out.empty:
                            out.reset_index(inplace=True)
                            rows.append(out)
            if not rows:
                return df
            merged = rows[0]
            for r in rows[1:]:
                merged = pd.merge(merged, r, on='islet_key', how='outer')
            # ensure df has islet_key to merge on
            if 'islet_key' not in df.columns and 'islet_id' in df.columns:
                try:
                    df = df.copy()
                    df['islet_key'] = 'Islet_' + df['islet_id'].astype(str)
                except Exception:
                    pass
            if 'islet_key' in df.columns:
                df = df.merge(merged, on='islet_key', how='left')
    except Exception:
        return df
    return df

def build_bin_summary(df: pd.DataFrame, donor_col: str, log_diam: np.ndarray, area: np.ndarray, binwidth: float, max_um: float = 300.0) -> pd.DataFrame:
    finite = np.isfinite(log_diam)
    x = log_diam[finite]
    a = area[finite]
    g = df.loc[finite, donor_col].astype(str).fillna("NA").to_numpy()
    if x.size == 0:
        return pd.DataFrame(columns=["donor_status","bin","x_mid","n_islets","total_area_um2"])

    lo = math.floor(np.nanmin(x))
    # Create edges up to 300 µm, and one extra open bin for >= 300
    max_log = math.log10(max_um)
    edges = np.arange(lo, max_log, binwidth)
    if len(edges) == 0 or edges[-1] < max_log:
        edges = np.append(edges, max_log)
    # assign bins
    digitized = np.digitize(x, edges) - 1
    # Values >= last edge go to a special plus bin
    plus_bin = len(edges) - 1
    over_mask = (log_diam[finite] >= max_log)
    digitized[over_mask] = plus_bin
    valid = (digitized >= 0) & (digitized <= plus_bin)

    tbl = pd.DataFrame({
        "bin": digitized[valid],
        "donor_status": pd.Series(g[valid], dtype="category"),
        "area": a[valid]
    })
    counts = tbl.groupby(["donor_status","bin"]).size().reset_index(name="n_islets")
    areas  = tbl.groupby(["donor_status","bin"])["area"].sum().reset_index(name="total_area_um2")
    out = pd.merge(counts, areas, on=["donor_status","bin"], how="outer").sort_values(["donor_status","bin"])

    mids = (edges[:-1] + edges[1:]) / 2.0
    out["x_mid"] = out["bin"].map(lambda i: (mids[int(i)] if (0 <= int(i) < len(mids)) else max_log))
    out["diam_mid_um"] = out["x_mid"].map(lambda v: float(10**v) if pd.notna(v) else max_um)
    out["diam_lo_um"] = out["bin"].map(lambda i: (float(10**edges[int(i)]) if (0 <= int(i) < len(edges)-1) else max_um))
    out["diam_hi_um"] = out["bin"].map(lambda i: (float(10**edges[int(i)+1]) if (0 <= int(i) < len(edges)-1) else np.inf))
    # reorder columns
    out = out[["donor_status","bin","x_mid","diam_mid_um","diam_lo_um","diam_hi_um","n_islets","total_area_um2"]]
    return out

def _sort_groups_preferred(groups: List[str]) -> List[str]:
    preferred = ["ND", "Aab+", "T1D"]
    # Map lower -> original for stable, case-insensitive matching
    lower_to_original = {str(g): str(g) for g in groups}
    # Build ordered list using the exact strings if present (case-insensitive)
    ordered = []
    used = set()
    for p in preferred:
        # find first group that matches p ignoring case
        match = next((g for g in groups if str(g).lower() == p.lower()), None)
        if match is not None:
            ordered.append(match)
            used.add(match)
    # append any remaining groups in original order
    ordered.extend([g for g in groups if g not in used])
    return ordered

def _preferred_order_from_df(df: pd.DataFrame, donor_col: str) -> List[str]:
    try:
        groups = df[donor_col].astype(str).dropna().unique().tolist()
    except Exception:
        return ["ND","Aab+","T1D"]
    return _sort_groups_preferred(groups)

def make_hist_and_area_plot(summary: pd.DataFrame, binwidth: float, title: Optional[str],
                            show_counts: bool = True, show_area: bool = True, show_error: bool = False,
                            overlay_marker_means: Optional[Dict[str, pd.DataFrame]] = None,
                            region_label_for_title: Optional[str] = None,
                            left_metric: Optional[str] = None,
                            overlay_fill: bool = False,
                            overlay_normalize: bool = False):
    if summary.empty:
        st.warning("No valid diameter values to plot.")
        return
    groups = summary["donor_status"].astype(str).unique().tolist()
    groups = _sort_groups_preferred(groups)
    # Prepare cumulative area by donor across increasing diameter
    # If composition overlays are present and left_metric == 'area', area is drawn on left y-axis.
    overlay_present = bool(overlay_marker_means)
    if overlay_present and left_metric:
        if left_metric == 'counts':
            show_counts, show_area = True, False
        elif left_metric == 'area':
            show_counts, show_area = False, True

    fig = go.Figure()
    # Modern color mapping (consistent across traces)
    color_map = {"ND": "#636EFA", "Aab+": "#EF553B", "T1D": "#00CC96"}
    # bars (counts) — x is diameter mid (µm) on a linear axis with real-diameter labels
    if show_counts:
        n_groups = max(1, len(groups))
        for gi, grp in enumerate(groups):
            sub = summary[summary["donor_status"].astype(str) == grp].copy()
            # base bin width per row
            base_w = (sub["diam_hi_um"] - sub["diam_lo_um"]).replace([np.inf, -np.inf], np.nan)
            base_w = base_w.fillna(method='ffill').fillna(method='bfill')
            # Each group's bar width is a fraction of the bin width
            bar_w = base_w * (0.9 / n_groups)
            # Compute per-group centers to avoid overlap
            sub["_x_center"] = sub["diam_lo_um"] + (gi + 0.5) * (base_w / n_groups)
            err = None
            if show_error:
                err = dict(type='data', array=np.sqrt(sub["n_islets"].fillna(0)).to_numpy(), visible=True)
            fig.add_trace(go.Bar(
                x=sub["_x_center"],
                y=sub["n_islets"],
                name=f"{grp} (count)",
                marker_color=color_map.get(grp, None),
                width=bar_w,
                offsetgroup=str(grp),
                opacity=0.9,
                error_y=err
            ))
    # lines (cumulative area)
    if show_area:
        for grp in groups:
            sub = summary[summary["donor_status"].astype(str) == grp]
            sub = sub.sort_values("diam_mid_um")
            cum = sub["total_area_um2"].fillna(0).cumsum()
            fig.add_trace(go.Scatter(
                x=sub["diam_mid_um"],
                y=cum,
                name=f"{grp} (cum area)",
                mode="lines+markers",
                line=dict(color=color_map.get(grp, None), width=2),
                yaxis=("y" if (overlay_present and left_metric == 'area') else "y2")
            ))

    # Overlay marker means (percentage) per bin per donor group
    overlay_present = overlay_marker_means is not None and len(overlay_marker_means) > 0
    max_overlay_val = 0.0
    if overlay_present:
        # Helper to convert hex to rgba with alpha
        def _hex_to_rgba(hex_color: str, alpha: float) -> str:
            hex_color = hex_color.lstrip('#')
            r = int(hex_color[0:2], 16)
            g = int(hex_color[2:4], 16)
            b = int(hex_color[4:6], 16)
            return f"rgba({r},{g},{b},{alpha})"

        for marker_name, dfm in overlay_marker_means.items():
            # Optional normalization per marker
            marker_max = float(dfm['mean_pct'].max(skipna=True)) if overlay_normalize else 1.0
            marker_max = marker_max if marker_max not in [0.0, float('nan')] else 1.0
            for grp in groups:
                subm = dfm[dfm["donor_status"].astype(str) == grp].sort_values("diam_mid_um")
                if subm.empty:
                    continue
                series_vals = subm["mean_pct"] / marker_max if overlay_normalize else subm["mean_pct"]
                local_max = float(series_vals.max(skipna=True)) if series_vals.size else 0.0
                if local_max > max_overlay_val:
                    max_overlay_val = local_max
                line_color = color_map.get(grp, None)
                fig.add_trace(go.Scatter(
                    x=subm["diam_mid_um"],
                    y=series_vals,
                    name=f"{marker_name} — {grp}",
                    mode="lines",
                    line=dict(color=line_color, width=3),
                    yaxis='y2',
                    fill=('tozeroy' if overlay_fill else None),
                    fillcolor=(_hex_to_rgba(line_color, 0.15) if (overlay_fill and line_color) else None)
                ))

    final_title = title
    if final_title is None and region_label_for_title:
        final_title = f"Distribution: Counts + Cumulative Area ({region_label_for_title})"
    # Axis titles and ranges depend on overlays
    y_title = ""
    y2_title = None
    if overlay_present:
        y2_title = "Composition" + (" (normalized)" if overlay_normalize else " (%)")
        if left_metric == 'area':
            y_title = "Cumulative islet area (µm²)"
        elif left_metric == 'counts':
            y_title = "Islet count per bin"
        else:
            y_title = "Islet count per bin" if show_counts else ("Cumulative islet area (µm²)" if show_area else "")
    else:
        y_title = "Islet count per bin"
        y2_title = "Cumulative islet area (µm²)"

    # Build x ticks including explicit 300+
    xtick_vals = [50, 100, 150, 200, 250, 300]
    xtick_text = ["50", "100", "150", "200", "250", "300+"]

    # Compute y2 range for composition to zoom into smaller values
    y2_range = None
    if overlay_present:
        try:
            ymax = max(5.0, max_overlay_val * 1.25)
            y2_range = [0, ymax]
        except Exception:
            y2_range = None

    fig.update_layout(
        title=final_title,
        barmode="group",
        xaxis=dict(title="Estimated islet diameter (µm)", tickmode="array", tickvals=xtick_vals, ticktext=xtick_text,
                   title_font=dict(size=28), tickfont=dict(size=22)),
        yaxis=dict(title=y_title, title_font=dict(size=28), tickfont=dict(size=22)),
        yaxis2=dict(title=y2_title, overlaying="y", side="right", tickformat=",.1f", range=y2_range,
                    title_font=dict(size=28), tickfont=dict(size=22)),
        legend=dict(orientation="h", yanchor="bottom", y=-0.22, xanchor="left", x=0.0, font=dict(size=21)),
        margin=dict(l=80, r=80, t=90, b=100),
        height=700,
        bargap=0.0,
        bargroupgap=0.08,
        font=dict(size=25)
    )
    # Force x-axis to 300+ max visual range
    try:
        xmin = float(np.nanmin(summary['diam_lo_um']))
        fig.update_xaxes(range=[xmin, 310])
    except Exception:
        pass
    st.plotly_chart(fig, use_container_width=True)

def find_column(df: pd.DataFrame, candidates: List[str], fallback_msg: str) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    st.error(fallback_msg)
    return None

def detect_density_columns(df: pd.DataFrame, region: str, cls: str, measure_prefix: str="dens") -> List[str]:
    patt = f"{region}__{measure_prefix}_{cls}"
    return [c for c in df.columns if c == patt]

# ------------------------- Sidebar: Data -------------------------

def _clear_and_rerun():
    try:
        st.cache_data.clear()
    except Exception:
        pass
    # Clear session state so UI tweaks reapply
    for k in list(st.session_state.keys()):
        del st.session_state[k]
    try:
        # Streamlit >= 1.25
        st.rerun()
    except Exception:
        try:
            # Backward compatibility
            st.experimental_rerun()
        except Exception:
            st.warning("Please manually rerun the app; rerun() not available in this Streamlit version.")

st.sidebar.button("Clear cache and rerun", on_click=_clear_and_rerun)

# Load consolidated data from workspace (no uploader)
data_path = "consolidated_streamlined.xlsx"
try:
    xls = pd.ExcelFile(data_path)
    use_sheet = "Islets" if "Islets" in xls.sheet_names else xls.sheet_names[0]
    df = pd.read_excel(data_path, sheet_name=use_sheet)
except Exception:
    st.error("Could not find consolidated_streamlined.xlsx. Please run the compiler to generate it.")
    st.stop()

# Optional structure mapping (e.g., TissueSquare -> Acinar)
if 'structure' in df.columns:
    df['_structure'] = df['structure'].astype(str)
    df['_structure'] = df['_structure'].replace({r'(?i)tissuesquare': 'Acinar'}, regex=True)
else:
    df['_structure'] = 'Islet'

# Column selections
donor_col = "donor_status" if "donor_status" in df.columns else st.sidebar.selectbox("Select donor status column", options=df.columns.tolist())
# Size uses Islet (core) region always
region = 'core'
binwidth = st.sidebar.slider("Bin width (log10 scale)", min_value=0.02, max_value=0.20, value=0.05, step=0.01)

# removed scatter aesthetics per request
marker_size = 5
marker_opacity = 0.45

# Filter donor_status
donor_levels = sorted([str(x) for x in pd.Series(df[donor_col]).dropna().unique()])
chosen = st.sidebar.multiselect("Filter donor_status", options=donor_levels, default=donor_levels)
df = df[df[donor_col].astype(str).isin(chosen)].copy()

# Structure filter removed if only a single value present
if "_structure" in df.columns:
    structures = sorted(df["_structure"].astype(str).unique().tolist())
    if len(structures) > 1:
        chosen_struct = st.sidebar.multiselect("Structure", options=structures, default=structures)
        df = df[df["_structure"].astype(str).isin(chosen_struct)].copy()

# Compute area & diameter
area, used_cols = coalesce_area_columns(df, region=region)
diam = compute_diameter_um(area)
log_diam = np.log10(diam)

## removed column detection & status per request

# ------------------------- Islet size filter -------------------------
finite_d = np.isfinite(diam)
if finite_d.any():
    dmin = float(np.nanmin(diam[finite_d]))
    dmax = float(np.nanmax(diam[finite_d]))
    sel_min, sel_max = st.slider(
        "Filter by estimated diameter (µm)",
        min_value=float(max(0.0, math.floor(dmin))),
        max_value=float(math.ceil(dmax)),
        value=(float(max(0.0, math.floor(dmin))), float(math.ceil(dmax))),
        step=1.0,
    )
    keep = finite_d & (diam >= sel_min) & (diam <= sel_max)
else:
    keep = np.zeros_like(diam, dtype=bool)

# Apply filter consistently
df = df.loc[keep].copy()
area = area[keep]
diam = diam[keep]
log_diam = np.log10(diam)

# ------------------------- Main: Distribution + overlays -------------------------

with st.sidebar.expander("Islets", expanded=True):
    opt_counts = st.checkbox("Show counts", value=True)
    opt_area = st.checkbox("Show cumulative area", value=True)
    opt_err = st.checkbox("Show error bars (counts)", value=False)
    st.caption("Islet composition (Islet region):")
    opt_ins = st.checkbox("INS_pos", value=False)
    opt_gcg = st.checkbox("GCG_pos", value=False)
    opt_sst = st.checkbox("SST_pos", value=False)
    opt_fill_islet = st.checkbox("Fill under overlays", value=False)
    opt_norm_islet = st.checkbox("Normalize overlays (per-marker)", value=False)
    left_metric_overlay = 'counts'
    if opt_ins or opt_gcg or opt_sst:
        left_metric_overlay = st.radio("Left y-axis (with composition):", options=['counts','area'], index=0, horizontal=True)
        if left_metric_overlay == 'counts':
            opt_counts, opt_area = True, False
        else:
            opt_counts, opt_area = False, True

summary = build_bin_summary(df, donor_col=donor_col, log_diam=log_diam, area=area, binwidth=binwidth)

# Prepare marker overlays (Islet/core)
overlay = {}
endo_maps = find_marker_columns_for(df, region='core', markers=ENDO_MARKERS)
def _add_overlay_if(name_key, marker_name):
    if name_key and marker_name in endo_maps.get(name_key, {}):
        col = endo_maps[name_key][marker_name]
        overlay[marker_name] = compute_marker_bin_means(df, donor_col, diam, binwidth, col)

if opt_ins:
    key = 'pct' if 'pct' in endo_maps else ('frac' if 'frac' in endo_maps else None)
    _add_overlay_if(key, 'INS')
if opt_gcg:
    key = 'pct' if 'pct' in endo_maps else ('frac' if 'frac' in endo_maps else None)
    _add_overlay_if(key, 'GCG')
if opt_sst:
    key = 'pct' if 'pct' in endo_maps else ('frac' if 'frac' in endo_maps else None)
    _add_overlay_if(key, 'SST')

make_hist_and_area_plot(
    summary,
    binwidth=binwidth,
    title="Distribution: Counts + Cumulative Area (Islet)",
    show_counts=opt_counts,
    show_area=opt_area,
    show_error=opt_err,
    overlay_marker_means=overlay,
    region_label_for_title="Islet",
    left_metric=left_metric_overlay,
    overlay_fill=opt_fill_islet,
    overlay_normalize=opt_norm_islet
)


# ------------------------- Peri-Islet distribution with overlays -------------------------

st.markdown("\n")
st.subheader("Peri-Islet Composition by Islet Size")
with st.sidebar.expander("Peri-Islet", expanded=True):
    opt_cd3e = st.checkbox("CD3E_pos (Peri-Islet)", value=False)
    opt_cd68 = st.checkbox("CD68_pos (Peri-Islet)", value=False)
    opt_cd163 = st.checkbox("CD163_pos (Peri-Islet)", value=False)
    opt_lgals3 = st.checkbox("LGALS3_pos (Peri-Islet)", value=False)
    st.caption("Targets (Peri-Islet) percentages:")
    opt_cap_pct = st.checkbox("Capillary %", value=False)
    opt_ner_pct = st.checkbox("Nerve %", value=False)
    opt_lym_pct = st.checkbox("Lymphatic %", value=False)
    opt_fill_band = st.checkbox("Fill under overlays (Peri-Islet)", value=False)
    opt_norm_band = st.checkbox("Normalize overlays (per-marker, Peri-Islet)", value=False)

overlay_band = {}
df = try_merge_targets_from_zip(df)
imm_maps_band = find_marker_columns_for(df, region='band', markers=IMMUNE_MARKERS)
def _add_overlay_band_if(marker_name):
    key = 'pct' if 'pct' in imm_maps_band else ('frac' if 'frac' in imm_maps_band else None)
    if key and marker_name in imm_maps_band.get(key, {}):
        col = imm_maps_band[key][marker_name]
        overlay_band[marker_name] = compute_marker_bin_means(df, donor_col, diam, binwidth, col)

if opt_cd3e: _add_overlay_band_if('CD3E')
if opt_cd68: _add_overlay_band_if('CD68')
if opt_cd163: _add_overlay_band_if('CD163')
if opt_lgals3: _add_overlay_band_if('LGALS3')

"""
Derive Peri-Islet target percentages robustly from densities (preferred) or counts (fallback).
"""
den_cols = ['band__dens_Capillary','band__dens_Nerve','band__dens_Lymphatic']
cnt_cols = ['band__count_Capillary','band__count_Nerve','band__count_Lymphatic']

def _safe_pct(num, den):
    den = den.replace(0, np.nan)
    return (num / den) * 100.0

if all(c in df.columns for c in den_cols) or all(c in df.columns for c in cnt_cols):
    if all(c in df.columns for c in den_cols):
        cap = pd.to_numeric(df['band__dens_Capillary'], errors='coerce')
        ner = pd.to_numeric(df['band__dens_Nerve'], errors='coerce')
        lym = pd.to_numeric(df['band__dens_Lymphatic'], errors='coerce')
    else:
        cap = pd.to_numeric(df['band__count_Capillary'], errors='coerce')
        ner = pd.to_numeric(df['band__count_Nerve'], errors='coerce')
        lym = pd.to_numeric(df['band__count_Lymphatic'], errors='coerce')
    denom = cap + ner + lym
    if opt_cap_pct:
        df['_band_pct_Capillary'] = _safe_pct(cap, denom)
        overlay_band['Capillary'] = compute_marker_bin_means(df, donor_col, diam, binwidth, '_band_pct_Capillary')
    if opt_ner_pct:
        df['_band_pct_Nerve'] = _safe_pct(ner, denom)
        overlay_band['Nerve'] = compute_marker_bin_means(df, donor_col, diam, binwidth, '_band_pct_Nerve')
    if opt_lym_pct:
        df['_band_pct_Lymphatic'] = _safe_pct(lym, denom)
        overlay_band['Lymphatic'] = compute_marker_bin_means(df, donor_col, diam, binwidth, '_band_pct_Lymphatic')

make_hist_and_area_plot(
    summary,
    binwidth=binwidth,
    title="Peri-Islet Composition by Islet Size",
    show_counts=False,
    show_area=False,
    show_error=False,
    overlay_marker_means=overlay_band,
    region_label_for_title="Peri-Islet",
    overlay_fill=opt_fill_band,
    overlay_normalize=opt_norm_band
)

# Download per-bin summary
csv_bytes = summary.to_csv(index=False).encode("utf-8")
st.download_button("Download per-bin summary CSV", data=csv_bytes, file_name="islet_size_bin_summary.csv", mime="text/csv")

# Facet option
st.checkbox("Facet by donor_status (small multiples)", value=False, key="facet_flag")
if st.session_state["facet_flag"] and not summary.empty:
    # Enforce preferred donor order if present
    cat_order = _sort_groups_preferred(summary["donor_status"].astype(str).unique().tolist())
    fig_facet = px.bar(
        summary,
        x="diam_mid_um",
        y="n_islets",
        facet_col="donor_status",
        facet_col_wrap=3,
        title="Counts by donor_status (faceted)",
        opacity=0.8,
        category_orders={"donor_status": cat_order},
    )
    fig_facet.update_xaxes(title_text="Estimated diameter (µm)")
    st.plotly_chart(fig_facet, use_container_width=True)

# ------------------------- Alternative: Binned Heatmap -------------------------

st.subheader("Binned Heatmap by Donor and Islet Size")

# Build 300+ bins
if np.isfinite(log_diam).any():
    lo = math.floor(np.nanmin(log_diam[np.isfinite(log_diam)]))
    max_log = math.log10(300.0)
    edges = np.arange(lo, max_log, binwidth)
    if len(edges) == 0 or edges[-1] < max_log:
        edges = np.append(edges, max_log)
else:
    edges = np.array([])

if edges.size >= 2 and len(df) > 0:
    mids = (edges[:-1] + edges[1:]) / 2.0
    mids_um = 10**mids
    bins = np.digitize(log_diam, edges) - 1
    plus_bin = len(edges) - 1
    bins = np.where(np.isfinite(log_diam) & (log_diam >= max_log), plus_bin, bins)
    df_heat = pd.DataFrame({
        "donor": df[donor_col].astype(str).to_numpy(),
        "bin": bins,
        "diam_label": np.where(bins >= 0,
                                 np.where(bins == plus_bin, '300+', mids_um[np.clip(bins, 0, len(mids_um)-1)].round(1).astype(str)),
                                 '')
    })

    # Build metric options: counts + mean percentages across regions/markers
    heat_metric_label_to_col = {"count": None}
    # Endocrine (Islet only)
    endo_maps_hm = find_marker_columns_for(df, region='core', markers=ENDO_MARKERS)
    meas_e = 'pct' if 'pct' in endo_maps_hm else ('frac' if 'frac' in endo_maps_hm else None)
    if meas_e:
        for m, col in endo_maps_hm[meas_e].items():
            heat_metric_label_to_col[f"mean {m} (percentage, Islet)"] = col
    # Immune/myeloid (all regions)
    for rgn in ['core','band','union']:
        maps = find_marker_columns_for(df, region=rgn, markers=IMMUNE_MARKERS)
        meas = 'pct' if 'pct' in maps else ('frac' if 'frac' in maps else None)
        if not meas:
            continue
        friendly = {"core":"Islet","band":"Peri-Islet","union":"Islet+20um"}.get(rgn, rgn)
        for m, col in maps[meas].items():
            heat_metric_label_to_col[f"mean {m} (percentage, {friendly})"] = col

    metric_choice = st.selectbox("Heatmap metric", options=list(heat_metric_label_to_col.keys()), index=0,
                                 help="Counts per bin or mean percentage per bin")
    show_values = st.checkbox("Show value labels", value=False)

    if metric_choice == 'count':
        agg = df_heat[df_heat["bin"] >= 0].groupby(["donor","bin","diam_label"]).size().reset_index(name="val")
    else:
        col = heat_metric_label_to_col[metric_choice]
        vals = pd.to_numeric(df[col], errors='coerce')
        if '__frac_' in _normalize(col) or _normalize(col).endswith('frac'):
            vals = vals * 100.0
        temp = df_heat.copy()
        temp['metric'] = vals
        temp = temp.replace([np.inf, -np.inf], np.nan)
        agg = temp[temp["bin"] >= 0].groupby(["donor","bin","diam_label"]).agg(val=("metric","mean")).reset_index()

    cat_order = _sort_groups_preferred(sorted(agg["donor"].unique().tolist()))
    agg["donor"] = pd.Categorical(agg["donor"], categories=cat_order, ordered=True)
    agg = agg.sort_values(["donor","bin"]) 
    mat = agg.pivot_table(index="donor", columns="diam_label", values="val", aggfunc="first").fillna(0)
    fig_hm = px.imshow(mat, aspect="auto", color_continuous_scale="Viridis", origin="lower", text_auto=show_values)
    if show_values:
        fig_hm.update_traces(text=mat.round(1), texttemplate="%{text}", textangle=0, textfont=dict(size=10))
    fig_hm.update_layout(margin=dict(l=60, r=60, t=80, b=80), font=dict(size=18))
    fig_hm.update_xaxes(title_text="Estimated diameter (µm)")
    fig_hm.update_yaxes(title_text="Donor status")
    st.plotly_chart(fig_hm, use_container_width=True)
else:
    st.info("Not enough data to create a heatmap.")

## Removed Object Density and Cell Type Analysis sections per request
## Donor aggregates removed per request

# ------------------------- Islet Image Browser -------------------------

st.header("🖼️ Islet Image Browser")
st.markdown("Choose a source for thumbnails: upload a ZIP or provide a manifest CSV with remote URLs.")

thumb_source = st.radio("Thumbnail source", ["ZIP upload", "Manifest CSV (remote)"])
img_zip = None
manifest_df = None
if thumb_source == "ZIP upload":
    img_zip = st.file_uploader("Upload ZIP of islet images", type=["zip"], key="zipuploader")
else:
    manifest_file = st.file_uploader(
        "Upload manifest CSV (columns: donor_id,islet_id,image_url)",
        type=["csv"],
        key="manifestuploader",
    )
    if manifest_file is not None:
        try:
            manifest_df = pd.read_csv(manifest_file)
        except Exception:
            st.error("Could not read manifest CSV. Ensure it has donor_id, islet_id, image_url columns.")
# Heuristics for ID columns
donor_id_col = None
islet_id_col = None
for col in df.columns:
    lc = col.lower()
    if donor_id_col is None and (lc == "donor_id" or lc.endswith("_donor_id") or lc == "donor"):
        donor_id_col = col
    if islet_id_col is None and (lc == "islet_id" or lc.endswith("_islet_id") or "islet_id" in lc):
        islet_id_col = col

donor_id_col = st.selectbox("Select donor ID column", options=[donor_id_col] + [c for c in df.columns if c != donor_id_col]) if donor_id_col else st.selectbox("Select donor ID column", options=df.columns.tolist())
islet_id_col = st.selectbox("Select islet ID column", options=[islet_id_col] + [c for c in df.columns if c != islet_id_col]) if islet_id_col else st.selectbox("Select islet ID column", options=df.columns.tolist())

donors_list = sorted([str(x) for x in pd.Series(df[donor_id_col]).dropna().unique()])
selected_donor = st.selectbox("Donor", options=donors_list)
islets_for_donor = sorted([str(x) for x in df[df[donor_id_col].astype(str) == selected_donor][islet_id_col].dropna().unique()])
selected_islet = st.selectbox("Islet ID", options=islets_for_donor)

def try_extract_zip(file) -> Optional[str]:
    if file is None:
        return None
    tdir = tempfile.mkdtemp(prefix="islet_imgs_")
    with zipfile.ZipFile(file) as zf:
        zf.extractall(tdir)
    return tdir

img_dir = try_extract_zip(img_zip)
if img_dir:
    st.caption("Images extracted to a temporary folder.")

def find_images_for(donor: str, islet: str, root: Optional[str]) -> List[str]:
    if not root or not os.path.isdir(root):
        return []
    hits = []
    pat_donor = re.escape(str(donor))
    pat_islet = re.escape(str(islet))
    for base, _, files in os.walk(root):
        for fn in files:
            if fn.lower().endswith((".png",".jpg",".jpeg",".tif",".tiff",".bmp")):
                f = fn.lower()
                if re.search(pat_donor, f) and re.search(pat_islet, f):
                    hits.append(os.path.join(base, fn))
    return hits

img_paths = find_images_for(selected_donor, selected_islet, img_dir)
def show_local_images(paths: List[str]):
    if not paths:
        return False
    st.success(f"Found {len(paths)} local image(s).")
    ncols = min(3, len(paths))
    cols = st.columns(ncols)
    for i, p in enumerate(paths):
        with cols[i % ncols]:
            try:
                im = Image.open(p)
                st.image(im, caption=os.path.basename(p), use_container_width=True)
            except Exception:
                st.warning(f"Could not open image: {os.path.basename(p)}")
    return True

def show_remote_images(manifest: pd.DataFrame, donor: str, islet: str):
    if manifest is None or manifest.empty:
        return False
    req_cols = {"donor_id","islet_id","image_url"}
    if not req_cols.issubset(set(c.lower() for c in manifest.columns)):
        st.error("Manifest must include donor_id, islet_id, image_url columns (case-insensitive).")
        return False
    # normalize columns
    cols_map = {c.lower(): c for c in manifest.columns}
    sub = manifest[(manifest[cols_map["donor_id"]].astype(str) == str(donor)) & (manifest[cols_map["islet_id"]].astype(str) == str(islet))]
    urls = [normalize_remote_url(u) for u in sub[cols_map["image_url"]].dropna().astype(str).tolist()]
    if not urls:
        return False
    st.success(f"Found {len(urls)} remote image(s).")
    ncols = min(3, len(urls))
    cols = st.columns(ncols)
    for i, u in enumerate(urls):
        with cols[i % ncols]:
            try:
                st.image(u, caption=os.path.basename(u), use_container_width=True)
            except Exception:
                st.warning(f"Could not display image URL: {u}")
    return True

displayed = False
if thumb_source == "ZIP upload":
    displayed = show_local_images(img_paths)
elif thumb_source == "Manifest CSV (remote)":
    displayed = show_remote_images(manifest_df, selected_donor, selected_islet)

if not displayed:
    if thumb_source == "ZIP upload":
        st.info("No matching local images found. Ensure filenames contain both donor and islet IDs.")
    else:
        st.info("No matching entries in manifest yet. Verify donor/islet IDs and URLs.")

# ------------------------- Gallery by size bin (optional) -------------------------

st.subheader("Gallery by Size Bin (across donors)")
use_bin_gallery = st.checkbox("Show gallery by size bin", value=False)
if use_bin_gallery:
    # Compute edges based on current filtered data and selected binwidth
    edges = compute_log_bin_edges(log_diam, binwidth)
    if edges.size < 2:
        st.info("Not enough data to build size bins.")
    else:
        # Build human-readable bin labels in µm
        bin_labels = []
        for i in range(len(edges)-1):
            lo_um = 10**edges[i]
            hi_um = 10**edges[i+1]
            bin_labels.append(f"[{format_um(lo_um)}–{format_um(hi_um)}) µm")
        sel_bin_idx = st.selectbox("Select size bin", options=list(range(len(bin_labels))), format_func=lambda i: bin_labels[i])
        donor_choices = st.multiselect("Restrict donors", options=donors_list, default=donors_list)
        max_imgs = st.slider("Max images to show", min_value=6, max_value=100, value=30, step=2)

        # Assign bins for each islet
        bins = assign_bins(log_diam, edges)
        mask = (bins == sel_bin_idx) & df[donor_id_col].astype(str).isin(donor_choices)
        df_bin = df.loc[mask, [donor_id_col, islet_id_col]].dropna().astype(str)
        pairs = list(df_bin.drop_duplicates().itertuples(index=False, name=None))

        st.caption(f"Found {len(pairs)} islets in selected bin across chosen donors.")

        # Resolve images for pairs
        resolved_paths = []
        resolved_urls = []
        if thumb_source == "ZIP upload":
            if img_dir:
                for donor_val, islet_val in pairs:
                    hits = find_images_for(donor_val, islet_val, img_dir)
                    for h in hits:
                        resolved_paths.append(h)
                        if len(resolved_paths) >= max_imgs:
                            break
                    if len(resolved_paths) >= max_imgs:
                        break
        else:
            if manifest_df is not None and not manifest_df.empty:
                cols_map = {c.lower(): c for c in manifest_df.columns}
                need = set((str(d), str(i)) for d, i in pairs)
                for _, row in manifest_df.iterrows():
                    try:
                        d = str(row[cols_map.get("donor_id")])
                        i = str(row[cols_map.get("islet_id")])
                        if (d, i) in need:
                            u = str(row[cols_map.get("image_url")])
                            u = normalize_remote_url(u)
                            resolved_urls.append(u)
                            if len(resolved_urls) >= max_imgs:
                                break
                    except Exception:
                        continue

        # Display gallery
        if thumb_source == "ZIP upload" and resolved_paths:
            ncols = min(4, len(resolved_paths))
            cols = st.columns(ncols)
            for idx, pth in enumerate(resolved_paths):
                with cols[idx % ncols]:
                    try:
                        im = Image.open(pth)
                        st.image(im, caption=os.path.basename(pth), use_container_width=True)
                    except Exception:
                        st.warning(f"Could not open image: {os.path.basename(pth)}")
        elif thumb_source == "Manifest CSV (remote)" and resolved_urls:
            ncols = min(4, len(resolved_urls))
            cols = st.columns(ncols)
            for idx, url in enumerate(resolved_urls):
                with cols[idx % ncols]:
                    try:
                        st.image(url, caption=os.path.basename(url), use_container_width=True)
                    except Exception:
                        st.warning(f"Could not display image URL: {url}")
        else:
            st.info("No images resolved for this bin selection.")
