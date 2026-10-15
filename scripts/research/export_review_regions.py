#!/usr/bin/env python3
"""Export exact analysis compartments as separate QuPath review overlays."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import shapely
from build_candidate import geometry_data
from common import sha256, write_json


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input-root', required=True, type=Path)
    parser.add_argument('--output-root', required=True, type=Path)
    args=parser.parse_args()
    root=args.input_root.resolve(strict=True);out=args.output_root.resolve()
    if out.exists():raise ValueError('Use a fresh review overlay folder')
    out.mkdir(parents=True)
    cohort=pd.read_csv(root/'source_manifest.csv',dtype=str)
    counts=pd.read_csv(root/'trajectory/counts.csv',dtype={'case_id':str})
    areas=counts[['islet_id','compartment','valid_area_um2']].drop_duplicates()
    if areas.duplicated(['islet_id','compartment']).any():raise ValueError('Inconsistent compartment areas')
    areas=areas.set_index(['islet_id','compartment']).valid_area_um2
    records=[]
    for donor in cohort.itertuples(index=False):
        cores,geoms,_=geometry_data(donor.geometry)
        tissue=geoms['Tissue'].difference(geoms.get('Exclude',shapely.GeometryCollection()))
        union=shapely.union_all([g for _,g in cores if g.is_valid and not g.is_empty])
        features=[]
        for name,g in cores:
            valid=g.is_valid and not g.is_empty and g.area>0
            regions={'core':g.intersection(tissue) if valid else shapely.GeometryCollection(),
                     'peri':g.buffer(20/.3774).difference(union).intersection(tissue) if valid else shapely.GeometryCollection()}
            for compartment,region in regions.items():
                identity=donor.case_id+'_'+name;area=region.area*.3774**2
                if not np.isclose(area,areas.loc[(identity,compartment)],atol=1e-6,rtol=1e-9):
                    raise ValueError(f'Review overlay area differs from count artifact: {identity}/{compartment}')
                records.append(dict(islet_id=identity,case_id=donor.case_id,compartment=compartment,
                                    valid_area_um2=area,geometry_valid=valid,area_matches_counts=True))
                if region.is_empty:continue
                features.append(dict(type='Feature',geometry=shapely.geometry.mapping(region),properties=dict(
                    objectType='annotation',name=f'Analysis review {compartment}: {name}',
                    classification=dict(name=f'QC analysis {compartment}',color=[40,150,220] if compartment=='core' else [190,80,200]),
                    measurements={'Analysis: valid area um2':area,'Analysis: band width um':0 if compartment=='core' else 20},
                    source_islet_id=identity,identity_review='pending; this overlay does not validate the source islet')))
        write_json(out/f'{donor.case_id}-analysis-core-peri20.geojson',dict(type='FeatureCollection',features=features))
        print(donor.case_id,len(features),'review overlays; areas match count artifact',flush=True)
    pd.DataFrame(records).to_csv(out/'area_verification.csv',index=False)
    write_json(out/'manifest.json',dict(source_candidate=str(root),coordinate_units='full-resolution source image pixels',
        micrometres_per_pixel=.3774,source_counts_sha256=sha256(root/'trajectory/counts.csv'),
        donors=cohort.case_id.tolist(),compartments_verified=len(records),review_project_modified=False,
        script_sha256=sha256(__file__),files={p.name:sha256(p) for p in sorted(out.glob('*.geojson'))}))
    (out/'README.md').write_text('''# Exact analysis compartments for QuPath review

Each donor file contains the saved islet core clipped to valid tissue and the
20-µm peri-islet band clipped to tissue and excluding all neighboring cores.
Areas are independently checked against every corresponding count-artifact row.
These are review annotations, not accepted biological islets or new cell objects.

Open the matching full-resolution image in a separate QuPath review entry and
run `importObjectsFromFile('/absolute/path/78-analysis-core-peri20.geojson')`,
substituting the matching donor and actual folder. Import once to avoid duplicates.
Existing upstream IsletExpanded objects use a different region definition.
All coordinates use source image pixels at 0.3774 µm/pixel; use the corrected
6539 image. Current review-project objects and cells have not been changed.
''')


if __name__=='__main__':main()
