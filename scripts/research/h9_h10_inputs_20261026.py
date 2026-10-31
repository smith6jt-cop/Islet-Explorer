#!/usr/bin/env python3
"""Assemble existing H9/H10 donor measurements without fitting clinical outcomes."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import shapely
from h1_features_20261022 import geometry, sha256

REPO = Path(__file__).resolve().parents[2]
APP = REPO / 'data/development/uniform-20donor-20261022-v10'
EDGE = REPO / 'output/h1_edge_corrected_20261024'
OUT = REPO / 'output/h5_h10_original_plan_20261026'


def main():
    for name in ('h9','h10'):
        (OUT/name).mkdir(parents=True, exist_ok=True)
    clinical_path = APP/'donor_metadata.csv'
    h2_path = REPO/'output/h1_h4_literature_methods_20261025/h2/h2_analysis_objects.csv'
    slopes_path = REPO/'output/h1_h4_disease_comparisons_20261026/h2/primary_donor_coefficients.csv'
    clinical = pd.read_csv(clinical_path)
    objects = pd.read_csv(h2_path)
    slopes = pd.read_csv(slopes_path)
    assert len(clinical)==20 and clinical.case_id.nunique()==20 and 78 not in clinical.case_id.values
    assert len(objects)==2937 and objects.islet_id.is_unique
    quantity = []
    sources = {str(p):sha256(p) for p in (clinical_path,h2_path,slopes_path)}
    for case in sorted(clinical.case_id):
        p = EDGE/'analysis_cells_by_donor'/f'{case}.parquet'
        cells = pd.read_parquet(p)
        assert cells.object_id.is_unique
        beta = cells[cells.edge_analysis_eligible & cells.phenotype.eq('Beta cell')]
        _,_,tissue,_,gp = geometry(APP,int(case))
        inside = shapely.covers(tissue,shapely.points(beta[['x_um','y_um']].to_numpy()))
        n = int(inside.sum()); area = tissue.area/1e6
        quantity.append(dict(case_id=case,tissue_beta_count=n,tissue_area_mm2=area,
            tissue_beta_density_mm2=n/area,beta_outside_tissue=int((~inside).sum()),
            beta_in_incomplete_objects=int((inside & beta.primary_compartment.eq('excluded_object').to_numpy()).sum())))
        sources[str(p)] = sha256(p); sources[str(gp)] = sha256(gp)
    quantity = pd.DataFrame(quantity)
    summaries = []
    for case,d in objects.groupby('case_id'):
        row = dict(case_id=case,donor_status=d.donor_status.iloc[0],immune_exposed_objects=len(d),
            core_immune_cells=int(d.n_immune_core.sum()),core_immune_positive_objects=int(d.n_immune_core.gt(0).sum()),
            mean_boundary_coverage=d.coliv_coverage_2um.mean(),
            mean_immune_localization=d.immune_localization_fraction.mean(),
            tissue_immune_density_mm2=d.tissue_immune_density_mm2.iloc[0])
        for col in ('coliv_coverage_2um','n_immune_nearby','field_density'):
            for label,q in [('min',0),('q10',.1),('q25',.25),('median',.5),('q75',.75),('q90',.9),('max',1)]:
                row[f'{col}_{label}'] = d[col].quantile(q)
        summaries.append(row)
    donors = clinical.merge(quantity,on='case_id',validate='one_to_one').merge(pd.DataFrame(summaries),on='case_id',validate='one_to_one')
    donors = donors.merge(slopes[['case_id','standardized_beta','within_donor_se']],on='case_id',validate='one_to_one')
    donors = donors.rename(columns={'standardized_beta':'coverage_core_immune_slope','within_donor_se':'coverage_core_immune_slope_se'})
    donors.to_csv(OUT/'h9/donor_inputs.csv',index=False)
    donors.to_csv(OUT/'h10/donor_inputs.csv',index=False)
    quantity.to_csv(OUT/'h9/tissue_beta_quantity.csv',index=False)
    quantity.to_csv(OUT/'h10/tissue_beta_quantity.csv',index=False)
    objects.to_csv(OUT/'h9/analysis_objects_existing_h2.csv',index=False)
    specification = dict(
        status='Data assembled; censored prediction fitting awaits method approval',
        outcome='Natural log C-peptide; retain <0.02 left-censoring and score probability below bound',
        baseline=['tissue_beta_density_mm2','age_years_npod','HbA1c_percent'],
        score='0.5 * (training-z(mean_boundary_coverage) - training-z(mean_immune_localization))',
        score_direction='Greater boundary coverage and lower core localization indicate greater support/lower access',
        score_selection='Two existing measurement components available in all 20 donors; equal weights fixed before clinical fitting; no association-based selection',
        aggregation='Arithmetic mean across the existing H1/H2 immune-exposed eligible objects per donor',
        preprocessing='Training donors only: median imputation, scaling, score calculation and predictive variance; fixed regularization',
        planned_validation='Leave one donor out; compare augmented-minus-baseline predictive log score; disease-conditioned sensitivity',
        quantity='Unique edge-eligible typed Beta cells inside saved Tissue polygons / Tissue area; valid cells in incomplete endocrine objects remain in tissue census',
        censoring_metadata='scripts/build_donor_metadata.py parse_c_peptide preserves numeric bound and censor flag; no full assay censoring rule provided',
        source_files=sources,donors=20,donor78_excluded=True,
        HbA1c_missing=donors.loc[donors.HbA1c_percent.isna(),'case_id'].tolist(),
        censored_C_peptide=donors.loc[donors.C_peptide_censored,['case_id','C_peptide_ng_mL']].to_dict('records'))
    (OUT/'h10/input_specification.json').write_text(json.dumps(specification,indent=2)+'\n')
    (OUT/'h9/source_provenance.json').write_text(json.dumps(dict(source_files=sources,quantity_definition=specification['quantity']),indent=2)+'\n')
    print(donors[['case_id','donor_status','tissue_beta_density_mm2','coverage_core_immune_slope','mean_boundary_coverage','mean_immune_localization','HbA1c_percent','C_peptide_ng_mL','C_peptide_censored']].to_string(index=False))


if __name__=='__main__':
    main()
