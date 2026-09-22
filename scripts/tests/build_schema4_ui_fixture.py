#!/usr/bin/env python3
"""Build explicitly synthetic schema-4 trajectory data for browser regression."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'research'))
from common import MARKERS,write_json

def build(root):
    root=Path(root);root.mkdir(parents=True,exist_ok=False)
    for folder in ['trajectory','typing','images']:(root/folder).mkdir()
    write_json(root/'candidate.json',dict(schema_version='4',synthetic=True))
    donors=['78']+[str(i) for i in range(1001,1021)]
    metadata=pd.DataFrame({'case_id':donors,'donor_status':np.repeat(['ND','Aab+','T1D'],7)})
    metadata.to_csv(root/'donor_metadata.csv',index=False);metadata.to_csv(root/'cohort.csv',index=False)
    names=['Beta cell','CHGA-supported hormone-unresolved endocrine']
    state='Beta cell: CD44 high'
    write_json(root/'typing/catalog.json',dict(schema_version=1,base_types=[dict(name=name) for name in names],states=[dict(name=state,parent='Beta cell')],endocrine_population='Endocrine union'))
    write_json(root/'typing/specification.json',dict(synthetic=True))
    coords=[];counts=[];denominators=[];summaries=[]
    for index,donor in enumerate(donors):
        for j in range(2):
            key=f'SyntheticIslet_{j+1}';iid=f'{donor}_{key}';e=20
            coords.append(dict(case_id=donor,islet_id=iid,islet_key=key,donor_status=metadata.donor_status.iloc[index],measured_cells=20+index,
                area_um2=100+10*index,eligible=True,supported=True,landmark=True,component=0,pseudotime_raw=(index*2+j)/42,pseudotime_percentile=(index*2+j+.5)/42,DC1=index,DC2=j))
            denominators.append(dict(case_id=donor,islet_id=iid,endocrine_count=e,validated=True))
            for comp in ['core','peri']:
                for typ,num in [(names[0],3),(names[1],1),(state,2)]:
                    counts.append(dict(case_id=donor,islet_id=iid,phenotype=typ,compartment=comp,count=num,validated=True))
                for pop in ['All cells','Endocrine union',*names,state]:
                    n=(20 if comp=='core' else 5) if pop=='All cells' else (4 if comp=='core' else 2)
                    for m,marker in enumerate(MARKERS):
                        available=marker!='CD45'
                        summaries.append(dict(case_id=donor,islet_id=iid,compartment=comp,population=pop,marker=marker,n_cells=n,
                          mean_log1p=(m+1)/31+index/100+(1 if comp=='peri' else 0),positive_count=1 if available else np.nan,positive_available=available))
    graph_layout=[]
    for row in coords:
        secondary=row['case_id']==donors[-1]
        graph_layout.append(dict(case_id=row['case_id'],islet_id=row['islet_id'],component=int(secondary),
            graph_x=row['DC1'],graph_y=row['DC2']+10*secondary,graph_landmark=True,layout_kind='synthetic componentwise spectral graph layout'))
        if secondary:row.update(component=1,supported=False,pseudotime_raw=np.nan,pseudotime_percentile=np.nan,DC1=np.nan,DC2=np.nan)
    for name,records in [('coordinates',coords),('counts',counts),('endocrine_counts',denominators),('marker_summaries',summaries),('graph_layout',graph_layout)]:
        table=pd.DataFrame(records);table.to_csv(root/'trajectory'/f'{name}.csv',index=False);table.to_parquet(root/'trajectory'/f'{name}.parquet',index=False)
    pairs=[(a,b) for a,b in zip(coords[:-1],coords[1:]) if a['component']==b['component']]
    pd.DataFrame({'from':[a['islet_id'] for a,b in pairs],'to':[b['islet_id'] for a,b in pairs]}).to_csv(root/'trajectory/edges.csv',index=False)
    write_json(root/'trajectory/marker_panel.json',dict(markers=MARKERS))
    write_json(root/'trajectory/validation.json',dict(status='synthetic_browser_fixture',coordinates_available=True,interpretation='Synthetic test data; not a scientific result.',median_lodo_absolute_spearman=.9,supported_primary_fraction=1))
    return root

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-root',type=Path,required=True)
    print(build(parser.parse_args().output_root))
