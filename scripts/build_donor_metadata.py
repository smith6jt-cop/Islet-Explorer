#!/usr/bin/env python3
"""
Derive per-donor clinical metadata for the 15 CODEX donors from an nPOD export.

The nPOD "All-Download" spreadsheet carries ~741 donors and 47 columns, including
free-text admission courses, clinical histories, medication lists and HLA types.
Only 15 donors and a handful of numeric/categorical fields are relevant here, so
this script extracts those into a small, reviewable CSV that CAN be committed.
The source spreadsheet is deliberately left out of git (`*.xlsx` in .gitignore).

What this adds over the app's existing donor metadata (age + sex, read from the
H5AD obs):
  - t1d_duration_years  (nPOD `diabetes_hx_years`) — disease duration
  - BMI, HbA1c_percent, C_peptide_ng_mL
  - the four AAb results as logicals, plus n_aab
  - age_group (<13 / >=13), the binning used by Steenbuck et al. 2025 (Bodenmiller
    lab, bioRxiv 2025.03.05.641526) for their age-endotype covariate

Two things worth knowing about this cohort, both surfaced by this script:

  * `t1d_duration_years` is NOT a usable regressor. All five T1D donors are
    recent-onset (0, 0, 0, 0.58, 0 years) and non-diabetic donors have no value at
    all. It is a cohort DESCRIPTOR — it establishes that we have no long-duration
    T1D donors — not a covariate. Do not fit a model on it.
  * Donor 112 appears in nPOD as `HDL112` (89 donors carry that prefix), and its
    nPOD age (3.72) disagrees with the app's H5AD age (3.00). The other 14 match
    exactly. This script warns rather than resolving it: the H5AD `age` stays the
    canonical modelling variable so nothing existing shifts, and `age_years_npod`
    rides alongside so the discrepancy stays visible.

Outputs:
  - data/donor_metadata.csv  (15 rows)

Usage:
    python scripts/build_donor_metadata.py
    python scripts/build_donor_metadata.py --source path/to/All-Download_*.xlsx
"""

import argparse
import glob
import os
import sys

import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEFAULT_OUT = os.path.join(REPO_ROOT, "data", "donor_metadata.csv")

# The 15 CODEX donors. Keys are the app's `Case ID`; values are the nPOD case_id
# where it differs (donor 112 is filed under the HDL prefix).
CODEX_DONORS = {
    "112": "HDL112",
    "6356": "6356", "6479": "6479", "6516": "6516", "6548": "6548",
    "6450": "6450", "6505": "6505", "6521": "6521", "6538": "6538", "6549": "6549",
    "6533": "6533", "6534": "6534", "6550": "6550", "6551": "6551", "6563": "6563",
}

# The app's `age` (H5AD obs), for the discrepancy check only.
H5AD_AGE = {
    "112": 3.00, "6356": 1.58, "6479": 21.67, "6516": 20.75, "6548": 20.24,
    "6450": 22.00, "6505": 20.59, "6521": 19.77, "6538": 19.14, "6549": 4.22,
    "6533": 3.75, "6534": 4.19, "6550": 25.06, "6551": 20.70, "6563": 14.56,
}

AAB_COLS = {
    "GADA_Result": "AAb_GADA",
    "IA_2A_Result": "AAb_IA2A",
    "mIAA_Result": "AAb_mIAA",
    "ZnT8A_Result": "AAb_ZnT8A",
}

# Steenbuck et al. pooled donors at 13 years. That cut is data-supported here
# rather than arbitrary: our ages are bimodal with a gap from 4.22 to 14.56.
AGE_GROUP_CUT = 13.0


def find_source(explicit=None):
    if explicit:
        if not os.path.exists(explicit):
            sys.exit(f"[donor-meta] source not found: {explicit}")
        return explicit
    hits = sorted(glob.glob(os.path.join(REPO_ROOT, "All-Download_From_nPOD_*.xlsx")))
    if not hits:
        sys.exit("[donor-meta] no All-Download_From_nPOD_*.xlsx in repo root; pass --source")
    return hits[-1]


def parse_aab(value, case_id, column):
    """'Positive'/'Negative' -> bool. Anything else is a data error, not a silent NA."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return False
    v = str(value).strip().lower()
    if v.startswith("pos"):
        return True
    if v.startswith("neg"):
        return False
    sys.exit(f"[donor-meta] donor {case_id}: unparseable {column} = {value!r}")


def parse_c_peptide(value):
    """
    C_peptide_ng_mL is a text column and one donor is left-censored ('<0.02').
    Keep the numeric bound and flag it, rather than dropping the donor or
    silently treating the bound as an observed value.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return (float("nan"), False)
    s = str(value).strip()
    if s.startswith("<"):
        return (float(s.lstrip("<").strip()), True)
    try:
        return (float(s), False)
    except ValueError:
        return (float("nan"), False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=None,
                    help="nPOD All-Download xlsx (default: newest in repo root)")
    ap.add_argument("--output", default=DEFAULT_OUT)
    args = ap.parse_args()

    src = find_source(args.source)
    print(f"[donor-meta] source: {src}")
    raw = pd.read_excel(src, sheet_name="data", dtype={"case_id": str})
    raw["case_id"] = raw["case_id"].astype(str).str.strip()
    print(f"[donor-meta] read {len(raw)} donors x {raw.shape[1]} columns")

    npod_to_app = {v: k for k, v in CODEX_DONORS.items()}
    sub = raw[raw["case_id"].isin(npod_to_app)].copy()

    missing = set(CODEX_DONORS.values()) - set(sub["case_id"])
    if missing:
        sys.exit(f"[donor-meta] not found in source: {sorted(missing)}")
    if len(sub) != len(CODEX_DONORS):
        sys.exit(f"[donor-meta] expected {len(CODEX_DONORS)} rows, matched {len(sub)}")

    rows = []
    for _, r in sub.iterrows():
        app_id = npod_to_app[r["case_id"]]
        c_pep, c_pep_censored = parse_c_peptide(r.get("C_peptide_ng_mL"))
        aabs = {out: parse_aab(r.get(src_col), app_id, src_col)
                for src_col, out in AAB_COLS.items()}
        age = float(r["age_years"])

        # nPOD leaves diabetes_hx_years blank for non-diabetic donors. That is
        # "not applicable", not zero -- keep it NA so it can never be averaged in.
        dur = r.get("diabetes_hx_years")
        dur = float(dur) if pd.notna(dur) else float("nan")

        rows.append({
            "case_id": app_id,
            "npod_case_id": r["case_id"],
            "donor_type": r["donor_type"],
            "age_years_npod": age,
            "age_group": "<13" if age < AGE_GROUP_CUT else ">=13",
            "t1d_duration_years": dur,
            "sex": r["sex"],
            "BMI": float(r["BMI"]) if pd.notna(r["BMI"]) else float("nan"),
            "HbA1c_percent": float(r["HbA1c_percent"]) if pd.notna(r["HbA1c_percent"]) else float("nan"),
            "C_peptide_ng_mL": c_pep,
            "C_peptide_censored": c_pep_censored,
            **aabs,
            "n_aab": sum(aabs.values()),
        })

    out = pd.DataFrame(rows).sort_values(["donor_type", "age_years_npod"])
    out = out.reset_index(drop=True)

    # ---- Loud checks on the columns that will feed models or tables ----
    for col in ["age_years_npod", "sex", "BMI", "HbA1c_percent", "C_peptide_ng_mL"]:
        bad = out.loc[out[col].isna(), "case_id"].tolist()
        if bad:
            sys.exit(f"[donor-meta] unexpected missing {col} for: {bad}")

    # Age discrepancy vs the H5AD: warn, never overwrite. The H5AD `age` stays
    # canonical for modelling so existing results do not shift underneath us.
    for _, r in out.iterrows():
        h5 = H5AD_AGE.get(r["case_id"])
        if h5 is not None and abs(h5 - r["age_years_npod"]) > 0.01:
            print(f"[donor-meta] WARNING donor {r['case_id']}: age differs -- "
                  f"H5AD {h5} vs nPOD {r['age_years_npod']} "
                  f"(H5AD value remains canonical for modelling)")

    # Disease duration is a descriptor, not a regressor -- say so at build time.
    dur = out["t1d_duration_years"].dropna()
    if len(dur):
        print(f"[donor-meta] t1d_duration_years: n={len(dur)} "
              f"range {dur.min()}-{dur.max()} y, unique {sorted(dur.unique())}")
        if dur.nunique() < 3:
            print("[donor-meta] NOTE: duration has too little variance to model; "
                  "it is a cohort descriptor (all T1D donors are recent-onset).")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"[donor-meta] wrote {args.output} ({len(out)} rows x {out.shape[1]} cols)")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
