"""
Data loading utilities for Islet Explorer Panel App.
Handles loading and preprocessing of master_results.xlsx data.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Optional, Tuple


# Default data path
DEFAULT_DATA_PATH = Path(__file__).parent.parent.parent.parent / "data" / "master_results.xlsx"

# Donor status colors (matching R Shiny app)
DONOR_STATUS_COLORS = {
    "ND": "#1f77b4",      # blue
    "Aab+": "#ff7f0e",    # orange
    "T1D": "#d62728"      # red
}

# Autoantibody columns
AAB_COLUMNS = ["GADA", "IA2A", "ZnT8A", "IAA", "mIAA"]

# Region types
REGION_TYPES = ["islet_core", "islet_band", "islet_union"]


class IsletDataLoader:
    """Loads and manages islet data from Excel workbook."""

    def __init__(self, data_path: Optional[str] = None):
        self.data_path = Path(data_path) if data_path else DEFAULT_DATA_PATH
        self._markers_df: Optional[pd.DataFrame] = None
        self._targets_df: Optional[pd.DataFrame] = None
        self._composition_df: Optional[pd.DataFrame] = None
        self._prepared_data: Optional[Dict] = None

    def validate_file(self) -> Tuple[bool, str]:
        """Check if data file exists and is readable."""
        if not self.data_path.exists():
            return False, f"Data file not found: {self.data_path}"
        try:
            # Quick validation by reading sheet names
            xl = pd.ExcelFile(self.data_path)
            required_sheets = ["Islet_Markers", "Islet_Targets", "Islet_Composition"]
            missing = [s for s in required_sheets if s not in xl.sheet_names]
            if missing:
                return False, f"Missing required sheets: {missing}"
            return True, "Data file validated successfully"
        except Exception as e:
            return False, f"Error validating data file: {e}"

    def load_all(self, force_reload: bool = False) -> Dict[str, pd.DataFrame]:
        """Load all data sheets from Excel file."""
        if self._prepared_data is not None and not force_reload:
            return self._prepared_data

        valid, msg = self.validate_file()
        if not valid:
            raise FileNotFoundError(msg)

        # Load each sheet
        self._markers_df = self._load_markers()
        self._targets_df = self._load_targets()
        self._composition_df = self._load_composition()

        # Prepare combined data
        self._prepared_data = self._prepare_data()
        return self._prepared_data

    def _load_markers(self) -> pd.DataFrame:
        """Load and preprocess Islet_Markers sheet."""
        df = pd.read_excel(self.data_path, sheet_name="Islet_Markers")

        # Try to load LGALS3 sheet and merge if exists
        try:
            lgals3_df = pd.read_excel(self.data_path, sheet_name="LGALS3")
            df = pd.concat([df, lgals3_df], ignore_index=True)
        except Exception:
            pass  # LGALS3 sheet may not exist

        # Standardize column names
        df.columns = df.columns.str.strip()

        # Ensure required columns exist
        required = ["Case ID", "Donor Status", "region", "marker"]
        for col in required:
            if col not in df.columns:
                # Try case-insensitive match
                for c in df.columns:
                    if c.lower() == col.lower():
                        df.rename(columns={c: col}, inplace=True)
                        break

        return df

    def _load_targets(self) -> pd.DataFrame:
        """Load and preprocess Islet_Targets sheet."""
        df = pd.read_excel(self.data_path, sheet_name="Islet_Targets")
        df.columns = df.columns.str.strip()
        return df

    def _load_composition(self) -> pd.DataFrame:
        """Load and preprocess Islet_Composition sheet."""
        df = pd.read_excel(self.data_path, sheet_name="Islet_Composition")
        df.columns = df.columns.str.strip()

        # Compute composition fractions if not present
        if "cells_total" in df.columns:
            for col in ["Ins_any", "Glu_any", "Stt_any"]:
                frac_col = f"{col}_frac"
                if col in df.columns and frac_col not in df.columns:
                    df[frac_col] = (df[col] / df["cells_total"]) * 100

        return df

    def _prepare_data(self) -> Dict[str, pd.DataFrame]:
        """Prepare combined data with computed metrics."""
        markers = self._markers_df.copy()
        targets = self._targets_df.copy()
        composition = self._composition_df.copy()

        # Extract islet key from region (e.g., "Islet_200_core" -> "Islet_200")
        if "region" in targets.columns:
            targets["islet_key"] = targets["region"].str.replace(r"_(core|band|union)$", "", regex=True)

        # Compute islet diameter from core area
        if "region_um2" in targets.columns:
            core_areas = targets[targets["region"].str.contains("core", case=False, na=False)]
            if not core_areas.empty:
                core_areas = core_areas.copy()
                core_areas["islet_diameter_um"] = 2 * np.sqrt(core_areas["region_um2"] / np.pi)
                # Merge diameter back to targets
                diameter_map = core_areas.set_index(["Case ID", "islet_key"])["islet_diameter_um"].to_dict()
                targets["islet_diameter_um"] = targets.apply(
                    lambda r: diameter_map.get((r.get("Case ID"), r.get("islet_key")), np.nan), axis=1
                )

        return {
            "markers": markers,
            "targets": targets,
            "composition": composition
        }

    @property
    def markers(self) -> pd.DataFrame:
        """Get markers dataframe."""
        if self._prepared_data is None:
            self.load_all()
        return self._prepared_data["markers"]

    @property
    def targets(self) -> pd.DataFrame:
        """Get targets dataframe."""
        if self._prepared_data is None:
            self.load_all()
        return self._prepared_data["targets"]

    @property
    def composition(self) -> pd.DataFrame:
        """Get composition dataframe."""
        if self._prepared_data is None:
            self.load_all()
        return self._prepared_data["composition"]

    def get_unique_values(self, column: str, df_name: str = "markers") -> list:
        """Get unique values for a column in specified dataframe."""
        df = getattr(self, df_name, None)
        if df is None or column not in df.columns:
            return []
        return sorted(df[column].dropna().unique().tolist())

    def get_donor_statuses(self) -> list:
        """Get unique donor status values."""
        return self.get_unique_values("Donor Status", "markers")

    def get_markers(self) -> list:
        """Get unique marker names."""
        return self.get_unique_values("marker", "markers")

    def get_regions(self) -> list:
        """Get unique region types."""
        # Return standard region types
        return REGION_TYPES

    def filter_data(
        self,
        df_name: str = "markers",
        donor_statuses: Optional[list] = None,
        region: Optional[str] = None,
        marker: Optional[str] = None,
        aab_filters: Optional[Dict[str, bool]] = None
    ) -> pd.DataFrame:
        """Filter dataframe by various criteria."""
        df = getattr(self, df_name).copy()

        if donor_statuses:
            df = df[df["Donor Status"].isin(donor_statuses)]

        if region and "region" in df.columns:
            df = df[df["region"].str.contains(region, case=False, na=False)]

        if marker and "marker" in df.columns:
            df = df[df["marker"] == marker]

        # Apply autoantibody filters for Aab+ donors
        if aab_filters:
            for aab, include in aab_filters.items():
                if aab in df.columns and not include:
                    # Exclude donors positive for this autoantibody
                    df = df[~((df["Donor Status"] == "Aab+") & df[aab])]

        return df


def compute_summary_stats(
    df: pd.DataFrame,
    value_col: str,
    group_col: str = "Donor Status",
    stat_type: str = "mean_se"
) -> pd.DataFrame:
    """Compute summary statistics grouped by a column."""

    def mean_se(x):
        return pd.Series({
            "mean": x.mean(),
            "se": x.std() / np.sqrt(len(x)) if len(x) > 1 else 0,
            "n": len(x)
        })

    def mean_sd(x):
        return pd.Series({
            "mean": x.mean(),
            "sd": x.std(),
            "n": len(x)
        })

    def median_iqr(x):
        return pd.Series({
            "median": x.median(),
            "q25": x.quantile(0.25),
            "q75": x.quantile(0.75),
            "n": len(x)
        })

    stat_funcs = {
        "mean_se": mean_se,
        "mean_sd": mean_sd,
        "median_iqr": median_iqr
    }

    func = stat_funcs.get(stat_type, mean_se)
    return df.groupby(group_col)[value_col].apply(func).unstack()


def detect_outliers(df: pd.DataFrame, value_col: str, n_sd: float = 3.0) -> pd.Series:
    """Detect outliers beyond n standard deviations."""
    mean = df[value_col].mean()
    std = df[value_col].std()
    return (np.abs(df[value_col] - mean) > n_sd * std)
