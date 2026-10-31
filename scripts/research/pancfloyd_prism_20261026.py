#!/usr/bin/env python3
"""Export explicit PancFloyd plot-data schemas to native Prism, CSV and TSV.

No scientific tests, model fitting, data-pipeline rebuilding or diagnostic-series
promotion. Rerun after the figure producers update their source tables.
"""
from pathlib import Path
from datetime import datetime, timezone
import importlib.util
import json
import re
import shutil
import hashlib
import numpy as np
import pandas as pd
from lxml import etree as ET

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "output/pancfloyd_20donor_20261026"
OUT = BASE / "prism"
GROUPS = ["ND", "Aab+", "T1D"]
INDEX = []
SOURCES = {}

# Reuse the already validated registration and GraphPad-schema writer.
spec = importlib.util.spec_from_file_location("validated_prism_export", ROOT / "scripts/research/export_prism_20261026.py")
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)
P.OUT = OUT
P.TABLES = []
P.INS_NOTICE = OUT / "NO_INS_NOTICE_PRESENT"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    path = Path(path)
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return pd.read_csv(path, dtype={"case_id": str, "islet_id": str})


def safe(x):
    return re.sub(r"[^a-zA-Z0-9]+", "_", str(x)).strip("_")


def labels(d, extra=None):
    if "islet_id" in d:
        z = d.case_id.astype(str) + " | " + d.islet_id.astype(str)
    elif "case_id" in d:
        z = "Donor " + d.case_id.astype(str)
    else:
        z = pd.Series(np.arange(1, len(d)+1).astype(str), index=d.index)
    if extra is not None:
        z = z + " | " + d[extra].astype(str)
    return z.to_numpy()


def reg(key, title, frame, numeric, description, source, *, x=None, panels="", kind="measurements", units="", table_type=None, endpoint_status="recreated"):
    if any("diagnostic_beta" in str(v) for v in [key, source]):
        raise ValueError("The diagnostic Beta smoother is not a figure series")
    assert key not in {s["key"] for s in P.TABLES}, key
    frame = frame.reset_index(drop=True)
    assert not frame.iloc[:, 0].duplicated().any(), (key, "duplicate row labels")
    if x:
        assert np.isfinite(pd.to_numeric(frame[x])).all(), (key, x)
    for name in numeric:
        vals = pd.to_numeric(frame[name])
        assert not np.isinf(vals).any(), (key, name)
    P.register(key, title, frame, numeric, description, source, x=x, table_type=table_type)
    INDEX.append({"table_key": key, "title": title, "reference_panels": panels,
                  "endpoint_status": endpoint_status,
                  "row_role": kind, "units": units, "rows": len(frame),
                  "x_field": x or "", "y_fields": ";".join(numeric),
                  "source_csv": str(Path(source).relative_to(ROOT)),
                  "source_sha256": sha(source), "description": description})


def xy_disease(key, title, d, x, y, source, desc, *, panels="", units="", row_label=None, extra_numeric=(), endpoint_status="recreated"):
    d = d.reset_index(drop=True)
    z = pd.DataFrame({"Row_ID": labels(d) if row_label is None else np.asarray(row_label), x: d[x]})
    nums = []
    for g in GROUPS:
        for v in [y, *extra_numeric]:
            name = g if v == y and not extra_numeric else f"{g} | {v}"
            z[name] = np.where(d.donor_status.eq(g), d[v], np.nan)
            nums.append(name)
    reg(key, title, z, nums, desc, source, x=x, panels=panels, units=units, endpoint_status=endpoint_status)


def donor_columns(key, title, d, y, source, desc, panels="", units=""):
    assert not d.case_id.duplicated().any()
    z = pd.DataFrame({"Donor_ID": d.case_id.to_numpy()})
    for g in GROUPS:
        z[g] = np.where(d.donor_status.eq(g), d[y], np.nan)
    reg(key, title, z, GROUPS, desc, source, panels=panels, units=units, kind="independent donor observations")


def matrix(key, title, d, row, col, value, source, desc, *, panels="", units="", row_order=None, col_order=None):
    z = d.pivot(index=row, columns=col, values=value)
    if row_order is not None:
        z = z.reindex(row_order)
    if col_order is not None:
        z = z.reindex(columns=col_order)
    z.columns = [str(c) for c in z.columns]
    numeric = z.columns.tolist()
    z = z.reset_index()
    reg(key, title, z, numeric, desc, source, panels=panels, units=units, kind="plotted matrix; not replicate observations", table_type="TwoWay")


def metric_tables(module):
    folder = BASE / module / "data"
    pm = read(BASE / module / "panel_manifest.csv")
    accepted_status = ["recreated", "explicit_area_endpoint_replacement", "restored_existing_cell_rule", "restored_structure_cell_endpoint"]
    for path in sorted(folder.glob("*_object_values.csv")):
        metric = path.name.removesuffix("_object_values.csv")
        records = pm.loc[pm.metric.eq(metric) & pm.status.isin(accepted_status)]
        panels = ";".join(records["panel"].astype(str))
        if not panels:
            continue
        endpoint_status = ";".join(records.status.drop_duplicates().astype(str))
        endpoint_note = ("Panel status: " + endpoint_status + ". " +
                         " ".join(records.endpoint.dropna().drop_duplicates().astype(str)) + ". ")
        if "replacement" in records:
            endpoint_note += " ".join(records.replacement.dropna().drop_duplicates().astype(str)) + ". "
        all_d = read(path)
        d = all_d[all_d.included.astype(str).str.lower().eq("true")].copy()
        unit = ("cells/mm²" if "density" in metric else "structure area / region area (%)" if "area" in metric
                else "structure-qualified cells / accepted cells (%)" if module=="structural_cells"
                else "typed cells / analysed cells (%)")
        if "units" in records and records.units.notna().any():
            unit = str(records.units.dropna().iloc[0])
        title = metric.replace("_", " ")
        desc = (endpoint_note + "Contributing eligible endocrine-object values within the reference diameter display range. "
                f"Row labels retain donor and object IDs. X is equivalent diameter in µm; Y units: {unit}. "
                "Objects from one donor are not independent biological replicates; published figure curves use object-bin summaries. "
                "The source CSV retains excluded objects, denominators and eligibility; this table includes only included=True rows.")
        xy_disease(metric+"_object_measurements", title+" — object values", d, "diameter_um", "value", path, desc,
                   panels=panels, units=unit, endpoint_status=endpoint_status)
        source = folder / f"{metric}_donor_bin_means.csv"
        if source.exists():
            q = read(source)
            q["Donor_series"] = q.donor_status + " | donor " + q.case_id
            wide = q.pivot(index="diam_mid", columns="Donor_series", values="value").sort_index()
            ys = wide.columns.tolist()
            wide = wide.reset_index()
            wide.insert(0, "Diameter_bin", wide.diam_mid.map(lambda x: f"{x-25:g}–{x+25:g} µm"))
            reg(metric+"_donor_bins", title+" — donor-bin means", wide, ys,
                endpoint_note+"One series per named donor; arithmetic mean of eligible object values within each 50-µm diameter bin. "
                "Blank donor/bin cells are unobserved bins, not zero. Repeated bins are not independent donors.", source,
                x="diam_mid", panels=panels, units=unit, kind="donor-bin observations", endpoint_status=endpoint_status)
        source = folder / f"{metric}_curve_values.csv"
        if source.exists():
            q = read(source).sort_values(["donor_status", "diam_mid"])
            row = q.donor_status + " | diameter " + q.diam_mid.astype(str)
            xy_disease(metric+"_plotted_curve", title+" — plotted means and SD limits", q, "diam_mid", "mean", source,
                       endpoint_note+"Exact plotted object-bin arithmetic means and lower/upper one-SD limits. Bounds are not additional observations. "
                       "The intervals are SD, not SEM or confidence intervals; group curve points are summaries, not donor replicates.",
                       panels=panels, units=unit, row_label=row, extra_numeric=["lower_sd", "upper_sd"], endpoint_status=endpoint_status)
        source = folder / f"{metric}_auc_values.csv"
        if source.exists():
            q = read(source).set_index("donor_status").reindex(GROUPS).reset_index()
            q["Disease_position"] = np.arange(1, len(q)+1)
            z = q[["donor_status", "Disease_position", "auc_percent_um"]]
            reg(metric+"_descriptive_auc", title+" — descriptive curve AUC", z, ["auc_percent_um"],
                "One trapezoidal integral of each disease group's curve on the common reported domain. "
                "Units percent·µm. These three group summaries are not individual donor values; no uncertainty or Prism analysis is embedded.",
                source, x="Disease_position", panels=panels, units="%·µm", kind="descriptive group summary")
        source = folder / f"{metric}_heatmap_values.csv"
        if source.exists():
            q = read(source)
            order = ["ANOVA", "ND vs Aab+", "ND vs T1D", "Aab+ vs T1D", "Kendall τ"]
            for value in ["minus_log10_q", "displayed_q"]:
                matrix(metric+"_"+value, title+" — "+value, q, "comparison", "diam_mid", value, source,
                       "Previously computed plotted statistics, supplied as numeric heatmap values only. "
                       "No new test is run by this exporter and no Prism analysis is embedded. "
                       "Blank protected pairwise/absent-bin cells remain blank; see the figure's original legend for multiplicity definitions.",
                       panels=panels, units="−log10(q)" if value=="minus_log10_q" else "adjusted p value", row_order=order)


def neighborhoods():
    folder = BASE / "neighborhoods/data"
    pm = read(BASE / "neighborhoods/panel_manifest.csv")
    summary_path = folder / "enrichment_group_summary.csv"
    summaries = read(summary_path)
    for row in pm.itertuples():
        key = row.image_id
        source = folder / f"{row.figure_stem}_plot_data.csv"
        if not source.exists():
            continue
        d = read(source)
        desc = f"Reference {row.image_id}. {row.aggregation}. Units: {row.units}."
        if "donor_mean_peri_enrichment_z" in d:
            y = "donor_mean_peri_enrichment_z"
            donor_columns(key+"_donors", row.figure_stem.replace("_", " "), d, y, source,
                          desc+" One observation per donor; disease columns are independent groups, not paired rows.", key, row.units)
            xy_disease(key+"_plotted_points", "Enrichment donor points — "+row.target, d, "plot_x", y, source,
                       desc+" X is a horizontal display offset only; donor IDs are row labels.", panels=key, units=row.units)
            q = summaries[summaries.panel.eq(row.figure_stem)].copy().reset_index(drop=True)
            q["lower_SEM"] = q["mean"]-q["sem"]
            q["upper_SEM"] = q["mean"]+q["sem"]
            q["Disease_position"] = q.donor_status.map(dict(zip(GROUPS, range(3))))
            z = q[["donor_status", "Disease_position", "mean", "sem", "lower_SEM", "upper_SEM"]]
            reg(key+"_summary", "Enrichment means and SEM — "+row.target, z, ["mean", "sem", "lower_SEM", "upper_SEM"],
                desc+" Plotted group mean ±SEM; summary columns and error limits are not replicate observations.", summary_path,
                x="Disease_position", panels=key, units=row.units, kind="descriptive group summary")
        elif "plotted_value_per_um" in d:
            row_ids = d.donor_status + " | bin " + d.bin_index.astype(str)
            for x, tag in [("center_symlog", "plotted_axis"), ("center_distance_um", "physical_distance")]:
                xy_disease(key+"_"+tag, row.figure_stem.replace("_", " ")+" — "+tag, d, x, "plotted_value_per_um", source,
                           desc+" Histogram bins are not biological replicates. The transformed x coordinate is only for matching the displayed axis; "
                           "the companion table preserves physical signed distance in µm. Exact source bin edges, widths and counts remain in the source CSV.",
                           panels=key, units=row.units, row_label=row_ids)
        elif "age_years_npod" in d:
            donor_columns(key+"_donors", "Donor age", d, "age_years_npod", source, desc, key, "years")
            xy_disease(key+"_plotted_points", "Donor age — plotted points", d, "plot_x", "age_years_npod", source,
                       desc+" X is only a horizontal display offset.", panels=key, units="years")
            src = folder / "slide07_image02_donor_age_summary.csv"
            if src.exists():
                q = read(src)
                q["Disease_position"] = q.donor_status.map(dict(zip(GROUPS, range(3))))
                z = q[["donor_status", "Disease_position", "mean_age_years", "sd_age_years"]]
                reg(key+"_summary", "Donor age — mean and SD", z, ["mean_age_years", "sd_age_years"],
                    "Exact age mean and sample SD shown in the donor-age figure. Summary fields are not individual donors.", src,
                    x="Disease_position", panels=key, units="years", kind="descriptive group summary")
        else:
            y = row.y
            assert y in d and "islet_id" in d, (source, y)
            z = pd.DataFrame({"Row_ID":labels(d), "plot_x":d.plot_x})
            ys = []
            for group in GROUPS:
                for donor in sorted(d.loc[d.donor_status.eq(group),"case_id"].unique(),key=int):
                    name = f"{group} | donor {donor}"
                    z[name] = np.where(d.case_id.eq(donor), d[y], np.nan); ys.append(name)
            reg(key+"_plotted_points", row.figure_stem.replace("_", " "), z, ys,
                desc+" One endocrine object per row and one named series per donor, matching the figure's donor colors. "
                "Donor and object IDs are retained. X is a display offset, not a biological predictor. "
                "Objects are not independent biological replicates.", source, x="plot_x", panels=key, units=row.units)


def trajectory():
    folder = BASE / "trajectory/data"
    pm = read(BASE / "trajectory/panel_mapping.csv")
    source = folder / "plot_object_representation.csv"
    d = read(source)
    desc = ("Saved current-app composition UMAP coordinates; corrected core and 20-µm neighborhood display subset. "
            "Full app fit: 1,215 objects; plotted subset: 1,126 objects, all 20 donors. No refit or reranking. "
            "Row labels preserve donor/object identity. Objects are not independent donors.")
    xy_disease("slide05_image01_umap_disease", "Composition UMAP — disease", d, "umap_1", "umap_2", source, desc,
               panels="slide05_image01", units="UMAP coordinates")
    z = pd.DataFrame({"Row_ID": labels(d), "umap_1": d.umap_1})
    clusters = sorted(d["leiden_0.5"].unique(), key=int)
    ys = []
    for cl in clusters:
        name = f"Leiden 0.5 cluster {cl}"
        z[name] = np.where(d["leiden_0.5"].eq(cl), d.umap_2, np.nan); ys.append(name)
    reg("slide05_image01_umap_leiden", "Composition UMAP — Leiden 0.5", z, ys, desc, source,
        x="umap_1", panels="slide05_image01", units="UMAP coordinates")
    source = folder / "cluster_composition_means.csv"
    q = read(source)
    present = q.groupby("phenotype").mean_core_fraction.max()
    q = q[q.phenotype.isin(present[present>0].index)]
    matrix("slide05_image02_cluster_means", "Mean core phenotype composition by cluster", q, "cluster", "phenotype", "mean_core_fraction", source,
           "Exact stacked-bar means of object-level core fractions. Denominator:all analysed core cells. "
           "Retired generic phenotypes are omitted without renormalization; bar sums may be below1. Cluster means are not donor replicates.",
           panels="slide05_image02", units="fraction of all core cells")
    source = folder / "cluster_composition_objects.csv"
    q = read(source)
    cols = [x for x in q if x.startswith("prop_") and q[x].max()>0]
    z = pd.DataFrame({"Row_ID": labels(q)})
    for col in cols:
        z[col] = q[col]
    reg("slide05_image02_contributing_objects", "Cluster-composition contributing object fractions", z, cols,
        "Each row identifies a donor and object. These are the values averaged in the cluster-composition bar chart; "
        "cluster assignments are in the source CSV. Object rows are not independent donors.", source,
        panels="slide05_image02", units="fraction of all core cells", kind="contributing object measurements")
    for row in pm.itertuples():
        if "_points.csv" not in row.data_csv:
            continue
        key = f"slide{row.reference_slide:02d}_image{row.reference_image:02d}"
        panel_title = "Beta-cell fraction along pseudotime" if key == "slide05_image03" else row.reference_panel
        sources = row.data_csv.split(";")
        source = folder / sources[0]
        q = read(source)
        for x,y,tag in [("x_plot","y_plot","plotted_points"),("pt_rank","value","measured_values")]:
            xy_disease(key+"_"+tag, panel_title+" — "+tag, q, x, y, source,
                       f"{row.change_note}. {row.units}. Exact object values or display coordinates as named; donor/object IDs retained. "
                       "Raw rank comes from the full 1,215-object fit, not a rerank of displayed objects. "
                       "Beta fraction has horizontal jitter only and no curve. Marker values use the served core-cell mean scale. "
                       "Object points are not independent donors.", panels=key, units=row.units)
        if len(sources)>1:
            source = folder / sources[1]
            assert "beta" not in source.name.lower()
            q = read(source)
            row_ids = q.donor_status + " | grid " + q.groupby("donor_status").cumcount().add(1).astype(str)
            xy_disease(key+"_loess", panel_title+" — descriptive app LOESS", q, "x", "y", source,
                       "Exact existing-app descriptive weighted LOESS coordinates and pointwise 95% band. "
                       "Bounds are not additional observations; they are not donor-replicated disease-comparison intervals. "
                       "No test or model is run by this exporter.", panels=key, units=row.units, row_label=row_ids,
                       extra_numeric=["lower","upper"])
    source = folder / "immune_heatmap_values.csv"
    q = read(source)
    q["Bin_label"] = q.apply(lambda r:f"bin{int(r['bin']):02d} | rank {r.midpoint:.6g}",axis=1)
    order = ["CD68","CD20","CD3e","CD4","CD8a"]
    for value in ["z_plot","bin_mean"]:
        matrix("slide05_image07_"+value, "Immune pseudotime heatmap — "+value, q, "marker", "Bin_label", value, source,
               "The 20 equal-width bins use the full fitted rank domain and only edge-eligible displayed objects for means. "
               "z_plot is standardized across the 20 marker-bin means using sample SD and clipped to ±2.5; "
               "bin_mean is the unstandardized served core-marker mean. Bins are not donor replicates.",
               panels="slide05_image07", units="z-score" if value=="z_plot" else "served core marker mean 0–1", row_order=order)


def main():
    structural_manifest = BASE / "structural_cells/panel_manifest.csv"
    if not structural_manifest.exists():
        raise RuntimeError("The restored structural-cell figure inputs are still being generated; rerun after their panel manifest is saved.")
    structural_panels = pd.read_csv(structural_manifest)
    expected = {"slide03_image01","slide03_image02","slide03_image03","slide03_image04",
                "slide04_image04","slide04_image05","slide09_image01","slide09_image02"}
    assert expected.issubset(set(structural_panels.panel)), "A restored structural-cell panel is missing from its manifest"
    for metric in structural_panels.loc[structural_panels.panel.isin(expected),"metric"].drop_duplicates():
        for suffix in ["object_values","donor_bin_means","curve_values"]:
            assert (BASE/"structural_cells/data"/f"{metric}_{suffix}.csv").exists(), f"Incomplete structural-cell export: {metric}_{suffix}"
    for folder in ["csv", "txt", "support"]:
        (OUT/folder).mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT/"output/prism_export_20261026/support/Prism8XMLSchema.xsd", OUT/"support/Prism8XMLSchema.xsd")
    previous_keys = set()
    if (OUT/"table_index.csv").exists():
        previous_keys = set(pd.read_csv(OUT/"table_index.csv").table_key)
    metric_tables("slides1_4"); metric_tables("structural_cells"); neighborhoods(); trajectory()
    restored_covered = set(p for t in INDEX if "/structural_cells/" in t["source_csv"]
                          for p in t["reference_panels"].split(";") if p)
    assert expected.issubset(restored_covered), "A restored cell panel was not exported; check source manifest status and metric"
    replacement_covered = set(p for t in INDEX if "explicit_area_endpoint_replacement" in t["endpoint_status"]
                             for p in t["reference_panels"].split(";") if p)
    assert {"slide03_image07", "slide03_image08"}.issubset(replacement_covered), "An explicit SMA area replacement was not exported"
    new_keys = {t["key"] for t in P.TABLES}
    for key in previous_keys-new_keys:
        for folder,suffix in [("csv","csv"),("txt","txt")]:
            stale = OUT/folder/f"{key}.{suffix}"
            if stale.exists():
                stale.unlink()
    validation = P.prism_file()
    original = OUT / "Islet_Explorer_plot_data.pzfx"
    xml = ET.parse(str(original))
    ns = {"p": P.NS}
    xml.xpath("/p:GraphPadPrismFile/p:Info/p:Title", namespaces=ns)[0].text = "PancFloyd figures — current 20-donor plot data"
    xml.xpath("/p:GraphPadPrismFile/p:Info/p:Notes", namespaces=ns)[0].text = (
        "Native Prism data tables only. The matching PNG/PDF/SVG figures and exact method legends are supplied outside this file. "
        "Donors are biological replicates. Object points, histogram bins, matrix entries, curve summaries and fitted curve coordinates "
        "must not be treated as independent donor observations. Row and series labels retain donor/object identities where applicable. "
        "Missing cells remain blank, not zero. No Prism analyses, model fits or graph layouts are embedded. "
        "The diagnostic Beta LOESS curve is excluded; the replacement Beta fraction panel contains measured label fractions only.")
    native = OUT / "PancFloyd_20donor_plot_data.pzfx"
    xml.write(str(native), encoding="UTF-8", xml_declaration=True, pretty_print=True)
    ET.XMLSchema(ET.parse(str(OUT/"support/Prism8XMLSchema_plain_cells.xsd"))).assertValid(ET.parse(str(native)))
    original.unlink()
    index = pd.DataFrame(INDEX)
    index.to_csv(OUT/"table_index.csv", index=False)
    index.to_csv(OUT/"table_index.txt", index=False, sep="\t")
    # Every source actually read is bound; report new source changes on rerun.
    validation.update({"exporter": str(Path(__file__).relative_to(ROOT)), "exporter_sha256":sha(__file__),
                       "created_utc":datetime.now(timezone.utc).isoformat(),
                       "native_file":native.name,"native_sha256":sha(native),
                       "source_files":SOURCES,"no_diagnostic_beta_loess":True,
                       "restored_structural_cell_panels":sorted(expected),
                       "explicit_area_endpoint_replacement_tables":[t["table_key"] for t in INDEX if "explicit_area_endpoint_replacement" in t["endpoint_status"]],
                       "no_new_scientific_test":True,
                       "csv_and_tsv_table_counts":len(P.TABLES),
                       "plotted_panels_covered":sorted(set(p for t in INDEX for p in t["reference_panels"].split(";") if p)),
                       "files":{str(p.relative_to(OUT)):sha(p) for p in sorted(OUT.rglob("*"))
                                if p.is_file() and p.name not in ["validation.json","README.md"]}})
    for table in P.TABLES:
        # CSV roundtrip is checked at17-digit precision, including identifiers.
        q = pd.read_csv(OUT/"csv"/(table["key"]+".csv"), dtype={table["frame"].columns[0]:str}, float_precision="round_trip")
        assert q.iloc[:,0].astype(str).tolist()==table["frame"].iloc[:,0].astype(str).tolist()
        for c in table["numeric"] + ([table["x"]] if table["x"] else []):
            assert np.array_equal(q[c].to_numpy(dtype=float),table["frame"][c].to_numpy(dtype=float),equal_nan=True),(table["key"],c)
    validation["csv_numeric_and_identifier_roundtrip_exact"] = True
    assert all(sha(ROOT/path)==expected for path,expected in SOURCES.items()), "A source changed during export; rerun for a consistent snapshot"
    validation["source_snapshot_unchanged_during_export"] = True
    (OUT/"validation.json").write_text(json.dumps(validation,indent=2)+"\n")
    (OUT/"README.md").write_text(
        "# PancFloyd current 20-donor Prism tables\n\n"
        f"`PancFloyd_20donor_plot_data.pzfx` contains {len(P.TABLES)} explicitly defined native Prism tables. "
        "`table_index.csv` maps each table to its source, units, row role and reference panels. "
        "Matching CSV and tab-delimited `.txt` tables are in `csv/` and `txt/`.\n\n"
        "The file contains data tables, without embedded Prism analyses or graph layouts. Use the supplied publication plots and legends to reconstruct graph formatting. "
        "Donor observations, object points, donor-bin means, group summaries, heatmap matrices and plotted curve coordinates are separated and labeled. "
        "Do not treat objects, bins, bounds or curve predictions as additional donors. Blank values are unobserved, not zero. "
        "The Beta fraction diagnostic LOESS is excluded. No quantitative INS endpoint is introduced.\n\n"
        "`validation.json` records GraphPad schema validation, exact numeric/row-label roundtrips, source hashes and covered panels. "
        "The published Prism8 schema has an overlapping optional wildcard that libxml2 cannot compile; the same previously validated stricter plain-cell subset is used. "
        "The native file has not been opened in the Prism application on this host.\n\n"
        "Rerun `/home/smith6jt/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 "
        "scripts/research/pancfloyd_prism_20261026.py` after figure producers update their files. "
        "Structure-qualified cell percentage/density panels are taken from `structural_cells/panel_manifest.csv` and remain separate from annotation-union area. "
        "SMA raster-area tables retain `explicit_area_endpoint_replacement` in the table index and descriptions; they do not represent arteriole cell percentages. "
        "Metrics are discovered through the relevant recreated/restored or explicit-area-replacement rows in the two source manifests. "
        "The exporter does not rerun any scientific analysis.\n")
    print(json.dumps({"native_file":str(native),"tables":len(P.TABLES),"covered_panels":len(validation["plotted_panels_covered"]),
                      "numeric_values":validation["numeric_values"],"validated":True},indent=2))


if __name__=="__main__":
    main()
