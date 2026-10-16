import gzip
import json
from pathlib import Path
import sys
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box, mapping
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from assemble_protein_expansion import process_donor, append_h5ad, validate_review, append_csv, donor_regions, neighborhood, ENDOCRINE
from research.common import MARKERS


def test_review_required_before_processing(tmp_path):
    (tmp_path/'preparation.json').write_text(json.dumps({'donors':[]}))
    (tmp_path/'manifest.json').write_text(json.dumps({'donors':[]}))
    p=tmp_path/'decision.json';p.write_text('{}')
    with pytest.raises(ValueError,match='owner acceptance'):
        validate_review(tmp_path,tmp_path,p)


def test_region_counts_and_append_preserve_originals(tmp_path):
    # Two measured core cells, one measured peri cell and an unmeasured core
    # object: raw QuPath and QC phenotype denominators must remain distinct.
    xy=np.array([[10,10],[20,20],[60,20],[25,25]],dtype=float)
    index=['6232:a','6232:b','6232:c','6232:d']
    obs=pd.DataFrame({'Object ID':['a','b','c','d'],'X_centroid':xy[:,0]*.3774,
        'Y_centroid':xy[:,1]*.3774,'Cell Area':[10.]*4,'Nucleus Area':[5.]*4,
        'Parent':['old']*4,'Age':[20.]*4,'Gender':['F']*4,'imageid':[6232]*4},index=index)
    raw_values=np.full((4,31),1000.);raw_values[-1]=np.nan
    raw=ad.AnnData(raw_values,obs=obs,var=pd.DataFrame(index=MARKERS))
    scaled=np.full((3,31),.1);scaled[:,MARKERS.index('INS')]=.8
    cells=ad.AnnData(scaled,obs=obs.iloc[:3].copy(),var=raw.var.copy())
    cells.raw=raw[:3].copy();cells.layers['log']=np.log1p(raw_values[:3])
    features=[dict(type='Feature',geometry=mapping(box(0,0,40,40)),properties=dict(name='Islet_1',classification=dict(name='Islet'))),
              dict(type='Feature',geometry=mapping(box(-100,-100,200,200)),properties=dict(name='tissue',classification=dict(name='Tissue')))]
    geometry=tmp_path/'geometry.gz'
    with gzip.open(geometry,'wt') as stream:json.dump(dict(type='FeatureCollection',features=features),stream)
    centroid=tmp_path/'centroids.tsv'
    pd.DataFrame(dict(object_id=['a','b','c','d'],x_source_px=xy[:,0],y_source_px=xy[:,1])).to_csv(centroid,sep='\t',index=False)
    for d in ['cells','donors','gson','additions']:(tmp_path/d).mkdir()
    snapshot=dict(geometry=str(geometry),centroids=str(centroid),pixel_size_um=.3774)
    meta={'Case ID':'6232','Donor Status':'ND','Age':20.,'Gender':'F',
          'AAb_GADA':np.nan,'AAb_IA2A':np.nan,'AAb_ZnT8A':np.nan,'AAb_mIAA':np.nan}
    frequencies={k:0.2 for k in [*ENDOCRINE,'Unknown']}
    result,tables=process_donor('6232',cells,raw,snapshot,meta,frequencies,tmp_path,['Cell Area','Nucleus Area'])
    assert result.n_obs==1
    assert result.obs.total_cells.iloc[0]==2
    assert tables['composition'][0]['cells_total']==3
    assert tables['composition'][0]['unmeasured_cells']==1
    assert tables['features'][0]['n_endocrine_core']==2
    assert result.obs['peri_count_Beta_cell'].iloc[0]==1
    assert set(pd.read_csv(tmp_path/'cells/6232_Islet_1.csv').cell_region)=={'core','peri'}
    old=result.copy();old.obs_names=['1_Islet_1'];old.obs['imageid']=1
    old.obs['islet_id']='1_Islet_1';old.obs['combined_islet_id']='1_Islet_1'
    old.obs.index.name='islet_id'
    old.obsm['X_scVI_mean']=np.ones((1,10));old.layers['measurement']=old.X.copy()
    joined=append_h5ad(old,result)
    np.testing.assert_array_equal(joined.X[:1],old.X)
    np.testing.assert_array_equal(joined.obsm['phenotype_counts'][:1],old.obsm['phenotype_counts'])
    assert np.isnan(joined.obsm['X_scVI_mean'][1]).all()
    joined.write_h5ad(tmp_path/'joined.h5ad')
    # Coordinate-roundoff tolerance must never admit a meaningful overlap.
    features.append(dict(type='Feature',geometry=mapping(box(20,0,50,40)),
                         properties=dict(name='Islet_2',classification=dict(name='Islet'))))
    with gzip.open(geometry,'wt') as stream:json.dump(dict(type='FeatureCollection',features=features),stream)
    with pytest.raises(ValueError,match='Overlapping reviewed cores'):
        donor_regions(snapshot,raw,pd.read_csv(centroid,sep='\t').set_index('object_id'))


def test_append_csv_preserves_existing_bytes(tmp_path):
    path=tmp_path/'data.csv';original=b'case_id,value\n1,0.1234567890123456789\n'
    path.write_bytes(original);append_csv(path,[dict(case_id=2,value=.2)])
    assert path.read_bytes().startswith(original)


def test_immune_distance_uses_nearest_cell_like_original_pipeline():
    labels=np.array(['Beta cell','CD4 Tcell','Macrophage'])
    coords=np.array([[0.,0.],[3.,4.],[6.,8.]])
    row=neighborhood(labels,np.array([0]),np.array([1,2]),coords,{'Beta cell':.5,'CD4 Tcell':.25,'Macrophage':.25})
    assert row['min_dist_immune_mean']==5.
