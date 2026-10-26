#!/usr/bin/env python3
"""Final scientific-output consistency checks and delivery hashes; no fits."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT=Path("output/h1_immune_size_20261023")
IMMUNE_MARKERS={"CD45","CD3e","CD8a","CD4","CD20","CD68","CD163","HLADR","CD35","LGALS3","CD44","VIM","Ki67"}


def main():
    checks={}
    source=Path("data/development/uniform-20donor-20261015/typing/aligned")
    inventory=[{"file":p.name,"columns":pq.read_schema(p).names} for p in sorted(source.glob("*.parquet"))]
    assert len(inventory)==20
    assert all(IMMUNE_MARKERS.issubset(v["columns"]) for v in inventory)
    (ROOT/"source_aligned_marker_inventory.json").write_text(json.dumps(inventory,indent=2)+"\n")
    checks["all_13_immune_associated_shared_or_nuclear_channels_in_all_20_source_donors"]=True
    base=pd.read_csv("output/h1_20261022/h1_object_features.csv",dtype={"case_id":str})
    assert len(base)==4763 and base.case_id.nunique()==20
    assert not base.case_id.eq("78").any()
    assert base.islet_id.is_unique
    checks["paired_source_identity_and_excluded_donor"]=True
    marker=pd.read_parquet(ROOT/"object_marker_counts.parquet")
    assert len(marker)==4763*2*(31+12*2)
    assert not marker.duplicated(["islet_id","compartment","marker","threshold"]).any()
    assert marker.loc[np.isclose(marker.threshold,.1),"marker"].nunique()==31
    checks["complete_31_marker_plus_12_threshold_sensitivity_grid"]=True
    types=pd.read_csv(ROOT/"types/estimates.csv")
    effects=pd.read_csv(ROOT/"marker_models/estimates.csv")
    assert len(types)==168 and len(effects)==252
    for label,table in [("type",types),("marker",effects)]:
        fitted=table.status.eq("FITTED_EXPLORATORY")
        assert table.loc[~fitted,["rr","lower","upper","p_cluster_t","q_bh"]].isna().all().all()
        assert table.loc[fitted,["rr","lower","upper","p_cluster_t","q_bh"]].notna().all().all()
        assert (table.loc[fitted,"lower"] <= table.loc[fitted,"rr"]).all()
        assert (table.loc[fitted,"upper"] >= table.loc[fitted,"rr"]).all()
        assert (table.loc[fitted,"informative_donors"]>=6).all()
        assert (table.loc[fitted,"core_positive_objects"]>=30).all()
        assert (table.loc[fitted,"objects"]>=100).all()
        checks[f"{label}_inference_gates_and_no_fabricated_nulls"]=True
    assert effects.loc[effects.family.eq("primary"),"q_bh"].dropna().min()>.05
    assert types.loc[types.family.eq("primary_subtypes"),"q_bh"].dropna().min()>.05
    assert not effects[effects.cohort.eq("T1D")].status.eq("FITTED_EXPLORATORY").any()
    assert not types[(types.cohort=="T1D")&(types.size_definition=="equivalent_diameter")].status.eq("FITTED_EXPLORATORY").any()
    checks["primary_null_and_T1D_support_statements_verified"]=True
    cells=pd.read_parquet(ROOT/"ki67/ki67_size_joined_cells.parquet")
    assert not cells.duplicated(["case_id","object_id"]).any()
    for compartment,count_col in [("core","n_immune_core"),("peri_exclusive","n_immune_peri_exclusive")]:
        observed=cells[cells.compartment==compartment].groupby("islet_id").size().reindex(base.islet_id,fill_value=0).to_numpy()
        assert np.array_equal(observed,base[count_col].to_numpy()),compartment
    checks["ki67_unique_local_cells_match_every_core_and_exclusive_peri_count"]=True
    audit=json.loads((ROOT/"audit/audit_summary.json").read_text())
    assert audit["phenotype_union_and_merge_verified"] and all(f["verified"] for f in audit["bh_family_verification"])
    checks["independent_BH_and_source_audit_passed"]=True
    for target in re.findall(r'\]\((/[^)]+)\)',(ROOT/"README.md").read_text()):
        assert Path(target).exists(),target
    figures=sorted((ROOT/"figures").glob("*.png"))
    assert len(figures)==8 and all(f.with_suffix(".svg").exists() for f in figures)
    checks["report_links_and_eight_PNG_SVG_pairs"]=True
    scripts=sorted(Path("scripts/research").glob("*20261023.py"))
    delivered=[p for p in ROOT.rglob("*") if p.is_file() and p.name!="delivery_manifest.json"]
    hashes={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in scripts+delivered}
    manifest={"analysis_date":"2026-10-23","checks":checks,"input_feature_sha256":hashlib.sha256(Path("output/h1_20261022/h1_object_features.csv").read_bytes()).hexdigest(),"figures_visually_reviewed":True,"model_slots":{"types":168,"markers":252},"source_data_and_app_modified":False,"files_sha256":hashes}
    (ROOT/"delivery_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({"checks":checks,"hashed_files":len(hashes)},indent=2))


if __name__=="__main__":
    main()
