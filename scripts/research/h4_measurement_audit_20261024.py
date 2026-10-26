#!/usr/bin/env python3
"""Audit H4 on current corrected objects; no surrogate hypothesis test."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import h1_analysis_20261022 as h1
from h234_existing_data_20261024 import OUT, CORRECTED, SRC, sha


def run():
    d=h1.eligible(pd.read_parquet(CORRECTED/'h1_object_features.parquet'))
    audit=pd.read_csv(OUT/'structural_evidence_audit.csv',dtype={'case_id':str})
    duct=audit.loc[audit.structure.eq('Duct')].set_index('case_id')
    rows=[];hashes={};unique_counts={}
    for case,g in d.groupby('case_id'):
        membership_path=CORRECTED/'object_memberships_by_donor'/f'{case}.parquet'
        hashes[str(membership_path)]=sha(membership_path)
        membership=pd.read_parquet(membership_path)
        unique_counts[case]=membership.loc[membership.islet_id.isin(g.islet_id)&
            membership.phenotype.eq('Macrophage'),'object_id'].nunique()
        p=Path(f'output/h1_immune_size_20261023/marker_features_by_donor/{case}/marker_profiles.parquet')
        hashes[str(p)]=sha(p)
        profiles=pd.read_parquet(p)
        profiles=profiles.loc[profiles.islet_id.isin(g.islet_id)&profiles.marker.eq('LGALS3')]
        for rec in g.to_dict('records'):
            obj=profiles.loc[profiles.islet_id.eq(rec['islet_id'])]
            for compartment in ['core','peri']:
                for pop,expected in [('Macrophage',rec[f'n_{compartment}_Macrophage']),
                                     ('IMMUNE_UNION',rec[f'n_immune_{compartment}']),
                                     ('ALL',rec[f'n_cells_{compartment}'])]:
                    x=obj.loc[obj.compartment.eq(compartment)&obj.phenotype.eq(pop)]
                    if expected:
                        if len(x)!=1 or x.iloc[0].n_cells!=expected:
                            raise ValueError('Marker population does not match corrected cells')
                        r=x.iloc[0]
                        nmeasured=int(r.n_measured);signal=float(r['sum']);npos=int(r.n_positive) if nmeasured else None
                    else:
                        if len(x):raise ValueError('Unexpected nonempty profile')
                        nmeasured=0;signal=0.;npos=0
                    rows.append(dict(case_id=case,donor_status=rec['donor_status'],
                        islet_id=rec['islet_id'],compartment=compartment,population=pop,
                        cells=expected,measured=nmeasured,missing=expected-nmeasured,
                        sum_measured_cell_means=signal,n_scaled_ge_0p10=npos,
                        complete_signal=nmeasured==expected,
                        epithelial_mask_accepted=bool(duct.loc[case,'validated'])))
    r=pd.DataFrame(rows);r.to_csv(OUT/'h4_lgals3_measurement_audit.csv',index=False)
    z=r.loc[r.population.eq('Macrophage')].groupby(['donor_status','case_id']).agg(
        objects=('islet_id','nunique'),macrophage_memberships=('cells','sum'),
        lgals3_measured=('measured','sum'),lgals3_missing=('missing','sum'),
        epithelial_mask_accepted=('epithelial_mask_accepted','first')).reset_index()
    z.to_csv(OUT/'h4_donor_readiness.csv',index=False)
    z['unique_nearby_macrophages']=z.case_id.map(unique_counts)
    z.to_csv(OUT/'h4_donor_readiness.csv',index=False)
    for p in [SRC/'typing/background_check.csv',SRC/'typing/scale.csv',
              Path('scripts/uniform/README.md'),Path('scripts/uniform/typing/background_check.py')]:
        hashes[str(p)]=sha(p)
    (OUT/'h4_source_hashes.json').write_text(json.dumps(hashes,indent=2))
    print(z.to_string(index=False))


if __name__=='__main__':run()
