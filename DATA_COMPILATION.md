Compiling donor spreadsheets into app-ready data

Overview

- Use `scripts/compile_donors.py` to combine per-donor CSV/XLSX files into a single `consolidated_streamlined.xlsx` with an `Islets` sheet.
- The script auto-detects core measurements and common marker columns; you can override detection with a mapping YAML if needed.

Inputs supported

- Directory or `.zip` containing per-donor files (CSV or XLSX). If an XLSX has multiple sheets, the script uses `Islets` if present, otherwise the first sheet.

What the script extracts

- `donor_id`, `islet_id`, `donor_status` (inferred if necessary from file names or columns).
- Area columns (prefers `islet_area_core_um2`, also considers `islet_area_union_um2` if found).
- Endocrine markers in `core`: INS, GCG, SST (any of pct/frac/dens/count).
- Immune/myeloid markers in `core`, `band`, `union`: CD3E, CD68, CD163, LGALS3 (any of pct/frac/dens/count).
- The script keeps detected marker columns as-is, preserving names the app can read.

Column name flexibility

- The script is case-insensitive and punctuation-tolerant. It recognizes measure synonyms (e.g., `percentage`, `density`, `cells_per_mm2`). It matches any token order as long as region, measure, and marker appear in the name.

Usage

- From repo root:

  python scripts/compile_donors.py --input donors.zip --output consolidated_streamlined.xlsx

  python scripts/compile_donors.py --input /path/to/donors_dir --output consolidated_streamlined.xlsx

Donor status mapping

- To apply donor_status from an external spreadsheet (e.g., `CODEX_Pancreas_Donors.xlsx`), pass:

  python scripts/compile_donors.py --input donors.zip --donor-map-xlsx CODEX_Pancreas_Donors.xlsx --output consolidated_streamlined.xlsx

- The mapping file should include columns resembling `Case ID` and `Donor Status` (case-insensitive). IDs are normalized to strings (e.g., 6356.0 → "6356").

Optional mapping

- If auto-detection needs help, provide a YAML mapping file. Example skeleton:

  donor_id: [donor_id, DonorID]
  islet_id: [islet_id, Islet]
  donor_status: [donor_status, Status]
  area_core: [islet_area_core_um2, core_area_um2]
  area_union: [islet_area_union_um2]
  markers:
    core:
      pct:
        INS: core__pct_INS
        GCG: core__pct_GCG
        SST: core__pct_SST
      dens:
        CD3E: core__dens_CD3E
        CD68: core__dens_CD68
        CD163: core__dens_CD163
        LGALS3: core__dens_LGALS3
    band:
      dens:
        CD163: band__dens_CD163
    union:
      dens:
        CD68: union__dens_CD68

- Mapping is not required; the script works heuristically without it.

Output

- Writes `consolidated_streamlined.xlsx` with sheet `Islets` in the repo root (or the path you specify).

Next

- Load the output into the app (sidebar file upload) to explore distributions, heatmaps, and thumbnails.
