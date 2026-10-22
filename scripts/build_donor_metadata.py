#!/usr/bin/env python3
"""
Derive donor metadata from an explicit cohort/key or the original nPOD cohort.

With --cohort, read case_id and donor_source (nPOD or HPAP) from a CSV and
demographics from --donor-key. Optional --source adds nPOD clinical fields only
to nPOD donors; unavailable clinical values remain null. Keep the candidate
output separate from the original CSV. inclusion_status is retained, not used
to silently remove identity-held donors or certify scientific inclusion.

Without --cohort, the original command, cohort, and output schema are retained.

With --uniform-cohort (scripts/uniform/cohort.csv), every donor of the 20-donor uniform cohort is read from the nPOD
sheet, all in the original schema plus the extra donor fields (owner, 2026-10-22: the nPOD sheet is the source for
every donor field, donor 112's age included). scripts/uniform/assemble_uniform.py --donor-metadata takes the age
from it.

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
    python scripts/build_donor_metadata.py --cohort cohort.csv \
        --donor-key CODEX_Pancreas_Donors.xlsx --output candidate/donor_metadata.csv
"""

import argparse
import glob
import hashlib
import math
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from excluded_donors import assert_no_excluded

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


def normalize_case_id(value):
    """Normalize spreadsheet integers and HPAP labels without truncating fractions."""
    if value is None or isinstance(value, bool) or pd.isna(value):
        raise ValueError("missing or invalid donor case ID")
    match = re.fullmatch(r"(?:HPAP-|HDL)?([0-9]+)(?:\.0+)?", str(value).strip(), re.I)
    if not match:
        raise ValueError(f"invalid donor case ID: {value!r}")
    return str(int(match.group(1)))


def read_table(path, sheet_name=0):
    if str(path).lower().endswith(".csv"):
        return pd.read_csv(path, dtype=str)
    return pd.read_excel(path, sheet_name=sheet_name)


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def optional_number(value, label):
    if value is None or pd.isna(value) or str(value).strip() == "":
        return float("nan")
    try:
        number = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"invalid numeric {label}: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"nonfinite {label}")
    return number


def read_cohort(path):
    cohort = pd.read_csv(path, dtype=str)
    required = {"case_id", "donor_source"}
    if not required.issubset(cohort.columns) or cohort.empty:
        raise ValueError("cohort must contain rows with case_id and donor_source")
    original_ids = cohort["case_id"].copy()
    cohort["case_id"] = cohort["case_id"].map(normalize_case_id)
    cohort["donor_source"] = cohort["donor_source"].str.strip().str.upper()
    if not cohort["donor_source"].isin(["NPOD", "HPAP"]).all():
        raise ValueError("donor_source must be explicitly nPOD or HPAP")
    for raw_id, source in zip(original_ids, cohort["donor_source"]):
        if str(raw_id).strip().upper().startswith("HPAP-") and source != "HPAP":
            raise ValueError("an HPAP case ID cannot be assigned to nPOD")
        if str(raw_id).strip().upper().startswith("HDL") and source != "NPOD":
            raise ValueError("an nPOD HDL case ID cannot be assigned to HPAP")
    if cohort["case_id"].duplicated().any():
        raise ValueError("duplicate normalized case_id in cohort")
    if "inclusion_status" not in cohort:
        cohort["inclusion_status"] = "selected"
    if not cohort["inclusion_status"].isin(["selected", "identity_hold"]).all():
        raise ValueError("inclusion_status must be selected or identity_hold; no rows are silently filtered")
    assert_no_excluded(cohort["case_id"], "Donor metadata cohort")
    return cohort


def build_cohort_metadata(cohort_path, donor_key_path, npod_source=None):
    """Build only the explicitly supplied cohort, keeping source identity distinct."""
    cohort = read_cohort(cohort_path)
    key = read_table(donor_key_path)
    required = {"Case ID", "Donor Status", "Age", "Gender"}
    if not required.issubset(key.columns):
        raise ValueError(f"donor key is missing columns: {sorted(required - set(key.columns))}")
    key = key.loc[key["Case ID"].notna()].copy()
    key["case_id"] = key["Case ID"].map(normalize_case_id)
    key = key.loc[key["case_id"].isin(cohort["case_id"])].copy()
    if key["case_id"].duplicated().any():
        raise ValueError("duplicate normalized case_id in donor key")
    missing = set(cohort["case_id"]) - set(key["case_id"])
    if missing:
        raise ValueError(f"cohort donors missing from donor key: {sorted(missing)}")
    key = key.set_index("case_id")

    npod = None
    if npod_source:
        npod = read_table(npod_source, sheet_name="data")
        if "case_id" not in npod:
            raise ValueError("nPOD export is missing case_id")
        # Preserve nPOD's HDL namespace. Numeric spreadsheet values normalize
        # without allowing HPAP-078 to become a lookup of nPOD donor 78.
        npod["case_id"] = npod["case_id"].map(
            lambda value: re.sub(r"\.0+$", "", str(value).strip()).upper()
        )

    rows = []
    key_hash = file_sha256(donor_key_path)
    source_hash = file_sha256(npod_source) if npod_source else None
    for _, selected in cohort.iterrows():
        app_id = selected["case_id"]
        source = "nPOD" if selected["donor_source"] == "NPOD" else "HPAP"
        demographics = key.loc[app_id]
        if str(demographics["Case ID"]).strip().upper().startswith("HPAP-") and source != "HPAP":
            raise ValueError("an HPAP donor-key row cannot supply nPOD demographics")
        if str(demographics["Case ID"]).strip().upper().startswith("HDL") and source != "nPOD":
            raise ValueError("an nPOD HDL donor-key row cannot supply HPAP demographics")
        age = optional_number(demographics["Age"], f"donor {app_id} Age")
        if pd.isna(age) or age < 0:
            raise ValueError(f"donor {app_id} requires a nonnegative Age")
        for field in ["Gender", "Donor Status"]:
            if pd.isna(demographics[field]) or not str(demographics[field]).strip():
                raise ValueError(f"donor {app_id} is missing {field}")
        source_id = CODEX_DONORS.get(app_id, app_id) if source == "nPOD" else f"HPAP-{int(app_id):03d}"
        row = {
            "case_id": app_id,
            "npod_case_id": source_id if source == "nPOD" else None,
            "donor_type": str(demographics["Donor Status"]).strip(),
            "age_years_npod": float("nan"),
            "age_group": "<13" if age < AGE_GROUP_CUT else ">=13",
            "t1d_duration_years": float("nan"),
            "sex": str(demographics["Gender"]).strip(),
            "BMI": float("nan"),
            "HbA1c_percent": float("nan"),
            "C_peptide_ng_mL": float("nan"),
            "C_peptide_censored": None,
            **{field: None for field in AAB_COLS.values()},
            "n_aab": None,
            "donor_source": source,
            "source_donor_id": source_id,
            "inclusion_status": selected["inclusion_status"],
            "age_years_key": age,
            "metadata_source": "donor_key_only",
            "donor_key_sha256": key_hash,
            "npod_source_sha256": None,
        }
        if source == "nPOD" and npod is not None:
            matches = npod.loc[npod["case_id"] == source_id]
            if len(matches) != 1:
                raise ValueError(f"expected exactly one nPOD row for {source_id}, found {len(matches)}")
            clinical = matches.iloc[0]
            for output, column in {
                "age_years_npod": "age_years", "t1d_duration_years": "diabetes_hx_years",
                "BMI": "BMI", "HbA1c_percent": "HbA1c_percent",
            }.items():
                row[output] = optional_number(clinical.get(column), f"donor {app_id} {column}")
            c_peptide = clinical.get("C_peptide_ng_mL")
            if c_peptide is not None and pd.notna(c_peptide):
                row["C_peptide_ng_mL"], row["C_peptide_censored"] = parse_c_peptide(c_peptide)
                if not math.isfinite(row["C_peptide_ng_mL"]):
                    raise ValueError(f"invalid C_peptide_ng_mL for donor {app_id}")
            for column, output in AAB_COLS.items():
                value = clinical.get(column)
                if value is not None and pd.notna(value):
                    row[output] = parse_aab(value, app_id, column)
            known = [row[field] for field in AAB_COLS.values()]
            row["n_aab"] = sum(known) if all(value is not None for value in known) else None
            row["metadata_source"] = "donor_key+npod_export"
            row["npod_source_sha256"] = source_hash
        rows.append(row)
    return pd.DataFrame(rows)


# Extra donor fields (owner, 2026-10-22), each read from the nPOD sheet. Numbers keep nPOD's units; a value nPOD
# records that cannot be true (a negative transit time) is left out with a warning.
EXTRA_NUMERIC = {
    "height_cm": "height_cm", "weight_kg": "weight_kg", "BMI_percentile": "BMI_percentile",
    "pancreas_weight_g": "pancreas_weight_grams", "admission_glucose_mg_dL": "admission_glucose_mg_dL",
    "peak_glucose_mg_dL": "peak_glucose_mg_dL", "ICU_days": "ICU_time_days",
    "transit_time_min": "transit_time_minutes", "RIN": "RIN",
}
STATUS_SHORT = {"No Diabetes": "ND", "Autoantibody Positive": "Aab+", "Type 1 Diabetes": "T1D"}


def leading_number(value):
    """'4.32 days' -> 4.32; blank -> NaN."""
    if value is None or (isinstance(value, float) and pd.isna(value)) or str(value).strip() == "":
        return float("nan")
    match = re.match(r"\s*(-?[0-9]+(?:\.[0-9]+)?)", str(value))
    if not match:
        raise ValueError(f"not a number: {value!r}")
    return float(match.group(1))


def build_uniform_metadata(cohort_path, source):
    """Every donor of the uniform cohort from the nPOD sheet: the original columns plus EXTRA_NUMERIC and four
    categories (race_ethnicity, recovery_type, cause_of_death, ABO_group)."""
    cohort = pd.read_csv(cohort_path, dtype=str)
    ids = [normalize_case_id(v) for v in cohort["case_id"]]
    assert_no_excluded(ids, "Uniform cohort")
    raw = pd.read_excel(source, sheet_name="data", dtype={"case_id": str})
    raw["app_id"] = raw["case_id"].astype(str).str.strip().map(
        lambda v: str(int(re.sub(r"^HDL", "", v))) if re.fullmatch(r"(HDL)?[0-9]+(\.0+)?", v) else None)
    rows = []
    for app_id, status in zip(ids, cohort["donor_status"]):
        hits = raw.loc[raw["app_id"] == app_id]
        if len(hits) != 1:
            sys.exit(f"[donor-meta] expected one nPOD row for donor {app_id}, found {len(hits)}")
        r = hits.iloc[0]
        if STATUS_SHORT.get(r["donor_type"]) != status:
            sys.exit(f"[donor-meta] donor {app_id}: cohort status {status}, nPOD {r['donor_type']}")
        c_pep, c_pep_censored = parse_c_peptide(r.get("C_peptide_ng_mL"))
        aabs = {out: parse_aab(r.get(col), app_id, col) for col, out in AAB_COLS.items()}
        age = float(r["age_years"])
        dur = r.get("diabetes_hx_years")
        row = {
            "case_id": app_id, "npod_case_id": r["case_id"], "donor_type": r["donor_type"],
            "age_years_npod": age, "age_group": "<13" if age < AGE_GROUP_CUT else ">=13",
            # Blank for donors without diabetes: not applicable, never zero.
            "t1d_duration_years": float(dur) if pd.notna(dur) else float("nan"),
            "sex": r["sex"], "BMI": optional_number(r.get("BMI"), "BMI"),
            # nPOD records HbA1c 0.0 for two donors (6090, 6232) whose notes give no value: not measured.
            "HbA1c_percent": hba1c if (hba1c := optional_number(r.get("HbA1c_percent"), "HbA1c")) > 0 else float("nan"),
            "C_peptide_ng_mL": c_pep, "C_peptide_censored": c_pep_censored,
            **aabs, "n_aab": sum(aabs.values()),
        }
        for out, col in EXTRA_NUMERIC.items():
            value = leading_number(r.get(col))
            if out == "transit_time_min" and value < 0:
                print(f"[donor-meta] WARNING donor {app_id}: transit time {value} min is impossible; left out")
                value = float("nan")
            row[out] = value
        row["race_ethnicity"] = r.get("race_ethnicity") if pd.notna(r.get("race_ethnicity")) else None
        # nPOD's own labels ("Organ Donor", "DCD Organ Donor"; DCD: donation after circulatory death).
        row["recovery_type"] = r.get("case_recovery_type") if pd.notna(r.get("case_recovery_type")) else None
        row["cause_of_death"] = r.get("cause_of_death") if pd.notna(r.get("cause_of_death")) else None
        abo = re.match(r"\s*(AB|A|B|O)", str(r.get("ABO_blood_type") or ""))
        row["ABO_group"] = abo.group(1) if abo else None
        rows.append(row)
    out = pd.DataFrame(rows)
    for col in ["age_years_npod", "sex", "BMI", "C_peptide_ng_mL"]:
        bad = out.loc[out[col].isna(), "case_id"].tolist()
        if bad:
            sys.exit(f"[donor-meta] unexpected missing {col} for: {bad}")
    for _, r in out.iterrows():
        h5 = H5AD_AGE.get(r["case_id"])
        if h5 is not None and abs(h5 - r["age_years_npod"]) > 0.01:
            print(f"[donor-meta] donor {r['case_id']}: H5AD age {h5}, nPOD {r['age_years_npod']} (nPOD is used)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=None,
                    help="nPOD All-Download XLSX (original default: newest in repo root; with --cohort: optional, CSV also accepted)")
    ap.add_argument("--output", default=DEFAULT_OUT)
    ap.add_argument("--cohort", help="CSV with case_id, donor_source, and optional inclusion_status")
    ap.add_argument("--donor-key", help="CSV/XLSX with Case ID, Donor Status, Age, and Gender")
    ap.add_argument("--uniform-cohort", help="scripts/uniform/cohort.csv: every donor from the nPOD sheet")
    args = ap.parse_args()

    if args.uniform_cohort:
        if os.path.realpath(args.output) == os.path.realpath(DEFAULT_OUT):
            ap.error("--uniform-cohort requires a separate --output; the original donor_metadata.csv is preserved")
        src = find_source(args.source)
        out = build_uniform_metadata(args.uniform_cohort, src)
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        out.to_csv(args.output, index=False)
        print(f"[donor-meta] wrote {args.output} ({len(out)} rows x {out.shape[1]} cols) from {os.path.basename(src)}")
        return

    if args.cohort:
        if not args.donor_key:
            ap.error("--cohort requires --donor-key")
        protected = [DEFAULT_OUT, args.cohort, args.donor_key]
        if args.source:
            protected.append(args.source)
        if os.path.realpath(args.output) in {os.path.realpath(path) for path in protected}:
            ap.error("--cohort requires a separate --output; preserve inputs and the original donor_metadata.csv")
        try:
            out = build_cohort_metadata(args.cohort, args.donor_key, args.source)
        except (ValueError, OSError) as exc:
            ap.error(str(exc))
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        out.to_csv(args.output, index=False)
        held = int((out["inclusion_status"] == "identity_hold").sum())
        print(f"[donor-meta] wrote {len(out)} explicit cohort rows; {held} remain on identity hold")
        print("[donor-meta] metadata does not establish image attribution or scientific inclusion")
        return
    if args.donor_key:
        ap.error("--donor-key requires --cohort")

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
