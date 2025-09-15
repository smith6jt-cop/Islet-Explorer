import os
import re
import glob
import pandas as pd
from pathlib import Path


def _parse_aab_flags(aabs_value: str) -> dict:
    """Parse an autoantibody string like "GADA, ZnT8A" into boolean flags.
    Recognizes: GADA, IA2A, ZnT8A, IAA, mIAA (case-insensitive).
    """
    flags = {
        'AAb_GADA': False,
        'AAb_IA2A': False,
        'AAb_ZnT8A': False,
        'AAb_IAA': False,
        'AAb_mIAA': False,
    }
    if not isinstance(aabs_value, str) or not aabs_value.strip():
        return flags
    parts = [p.strip() for p in aabs_value.split(',') if p.strip()]
    norm = set()
    for p in parts:
        token = re.sub(r'[^A-Za-z0-9]', '', p).upper()
        if token:
            norm.add(token)
    # Set flags based on normalized tokens
    if 'GADA' in norm:
        flags['AAb_GADA'] = True
    if 'IA2A' in norm:
        flags['AAb_IA2A'] = True
    # ZnT8A can appear as ZNT8A after normalization
    if 'ZNT8A' in norm or 'ZNT8' in norm:
        flags['AAb_ZnT8A'] = True
    # Distinguish mIAA vs IAA
    if 'MIAA' in norm:
        flags['AAb_mIAA'] = True
    if 'IAA' in norm:
        flags['AAb_IAA'] = True
    return flags


def load_donor_key(key_path: str) -> pd.DataFrame:
    key = pd.read_excel(key_path, sheet_name=0)
    # Normalize Case ID as zero-padded 4-digit string for joining
    key['Case ID'] = key['Case ID'].astype(str).str.zfill(4)
    # Expand autoantibody string column into logical indicator columns
    if 'Aabs' in key.columns:
        flags_df = key['Aabs'].apply(_parse_aab_flags).apply(pd.Series)
        for c in flags_df.columns:
            key[c] = flags_df[c].astype(bool)
    return key


def parse_case_id_from_filename(path: str) -> str:
    m = re.search(r"_(\d+)\.tsv$", os.path.basename(path))
    if not m:
        raise ValueError(f"Could not parse Case ID from filename: {path}")
    return m.group(1).zfill(4)


def group_result_files(results_dir: Path) -> dict:
    groups = {}
    for tsv in sorted(results_dir.glob('*.tsv')):
        base = tsv.name
        # Group by prefix before the last underscore + ID
        # e.g., 'islet_composition_6505.tsv' -> 'islet_composition'
        prefix = re.sub(r"_\d+\.tsv$", "", base)
        groups.setdefault(prefix, []).append(tsv)
    return groups


def build_sheet_dataframe(files: list[Path], key_df: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        case_id = parse_case_id_from_filename(str(f))
        df.insert(0, 'Case ID', case_id)
        df.insert(1, 'source_file', f.name)
        frames.append(df)

    if not frames:
        return pd.DataFrame()

    combined = pd.concat(frames, ignore_index=True)
    # Merge donor metadata on 'Case ID'
    merged = combined.merge(key_df, on='Case ID', how='left')

    # Reorder columns: donor metadata first if present
    donor_cols = ['Case ID', 'Donor Status', 'Age', 'Gender', 'Aabs',
                  'AAb_GADA', 'AAb_IA2A', 'AAb_ZnT8A', 'AAb_IAA', 'AAb_mIAA']
    existing_donor_cols = [c for c in donor_cols if c in merged.columns]

    other_cols = [c for c in merged.columns if c not in existing_donor_cols]
    ordered_cols = existing_donor_cols + other_cols
    merged = merged[ordered_cols]
    return merged


def write_master_excel(results_dir: str, key_path: str, out_path: str) -> list[tuple[str, list[str]]]:
    results_dir = Path(results_dir)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    key_df = load_donor_key(key_path)
    groups = group_result_files(results_dir)

    # Map prefixes to readable sheet names if desired
    sheet_name_map = {
        'islet_composition': 'Islet_Composition',
        'islet_markers': 'Islet_Markers',
        'islet_targets': 'Islet_Targets',
        'LGALS3': 'LGALS3',
    }

    with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
        sheet_infos: list[tuple[str, list[str]]] = []
        for prefix, files in sorted(groups.items()):
            # Derive grouping key (prefix without ID suffix already computed)
            # Normalize to a concise sheet name
            base_prefix = prefix
            sheet_name = sheet_name_map.get(base_prefix, base_prefix[:31])  # Excel sheet name limit 31

            df_sheet = build_sheet_dataframe(files, key_df)
            # Ensure at least an empty sheet exists
            if df_sheet.empty:
                df_sheet = pd.DataFrame({'Case ID': [], 'source_file': []})

            df_sheet.to_excel(writer, sheet_name=sheet_name, index=False)
            sheet_infos.append((sheet_name, list(df_sheet.columns)))

    return sheet_infos


if __name__ == '__main__':
    # Defaults tuned to repository structure
    results_dir = 'data/results'
    key_path = 'CODEX_Pancreas_Donors.xlsx'
    out_path = 'data/master_results.xlsx'

    infos = write_master_excel(results_dir, key_path, out_path)
    print(f"Wrote master workbook: {out_path}")
    for name, cols in infos:
        print(f"Sheet: {name}")
        print("Columns:", ", ".join(cols))
