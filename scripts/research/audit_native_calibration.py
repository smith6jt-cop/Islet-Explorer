#!/usr/bin/env python3
"""Reconcile native analysis calibration with independently encoded OME units."""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import pandas as pd
import tifffile
from common import write_json,sha256

UNITS={'µm':1,'um':1,'mm':1000,'nm':.001,'m':1e6}


def audit(root):
    cohort=pd.read_csv(root/'cohort.csv',dtype=str);rows=[]
    for row in cohort.to_dict('records'):
        cid=row['case_id'];native=json.loads((root/'native_cells'/cid/'manifest.json').read_text())
        path=Path(native['source_native']);project=path.parents[2]/'project.qpproj'
        document=json.loads(project.read_text())
        entry=next(x for x in document['images'] if int(x['entryID'])==int(path.parent.name))
        calibration=entry['serverBuilder'].get('metadata',{}).get('pixelCalibration',{})
        x,y=calibration['pixelWidth'],calibration['pixelHeight']
        source=json.loads(Path(row['boundary_source']).read_text());declared=source['pixel_size_um']
        if x['unit']!='µm' or y['unit']!='µm' or x['value']!=declared or y['value']!=declared:
            raise ValueError(f'{cid}: native analysis calibration differs from explicit source specification')
        with tifffile.TiffFile(row['image']) as image:
            xml=image.ome_metadata
            pixels=next(e for e in ET.fromstring(xml).iter() if e.tag.rsplit('}',1)[-1]=='Pixels')
        ome_x=float(pixels.attrib['PhysicalSizeX'])*UNITS[pixels.attrib.get('PhysicalSizeXUnit','µm')]
        ome_y=float(pixels.attrib['PhysicalSizeY'])*UNITS[pixels.attrib.get('PhysicalSizeYUnit','µm')]
        if max(abs(ome_x-declared),abs(ome_y-declared))/declared>.001:
            raise ValueError(f'{cid}: OME/native difference exceeds rounding tolerance; investigate before analysis')
        rows.append(dict(case_id=cid,analysis_pixel_size_um=declared,native_qupath_x_um=x['value'],
                         native_qupath_y_um=y['value'],native_project=str(project),native_project_sha256=sha256(project),
                         ome_x_um=ome_x,ome_y_um=ome_y,ome_relative_difference=(ome_x-declared)/declared,
                         ome_metadata_sha256=hashlib.sha256(xml.encode()).hexdigest(),
                         authority='Explicit native QuPath calibration used for source CellArea and geometry; acquisition OME rounding recorded separately'))
    write_json(root/'preprocessing/calibration.json',dict(status='all_manifest_native_analysis_calibration_matches',
               cohort_sha256=sha256(root/'cohort.csv'),donors=rows,code_sha256=sha256(__file__)))
    pd.DataFrame(rows).to_csv(root/'preprocessing/calibration.csv',index=False)
    print(f'{len(rows)} donors: native analysis calibration matches; OME differences recorded')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-root',type=Path,required=True)
    audit(p.parse_args().input_root.resolve(strict=True))
