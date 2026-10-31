#!/usr/bin/env python3
"""Recreate the reference field maps from the accepted 20-donor app data."""
from pathlib import Path
import gzip
import json
import re
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/pancfloyd-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import shapely
from shapely.geometry import shape
from shapely.affinity import scale
from h1_tissue_figure_20261022 import draw_shape

ROOT = Path(__file__).resolve().parents[2]
APP = ROOT/'data/development/uniform-20donor-20261022-v10'
OUT = ROOT/'output/pancfloyd_20donor_20261026/context'
PX = .3774
for sub in ('figures', 'data'):
    (OUT/sub).mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10, 'axes.titlesize':11,
                     'axes.titleweight':'bold', 'svg.fonttype':'none', 'pdf.fonttype':42})
palette = dict(re.findall(r'"([^"]+)"\s*=\s*"(#[A-Fa-f0-9]+)"',
    (ROOT/'app/shiny_app/R/drilldown_helpers.R').read_text().split('PHENOTYPE_COLORS_HIGH_CONTRAST')[0]))
with gzip.open(APP/'gson/6533.geojson.gz','rt') as f:
    geo = json.load(f)
features = []
for f in geo['features']:
    p=f.get('properties') or {}
    cl=(p.get('classification') or {}).get('name')
    g=shapely.make_valid(scale(shape(f['geometry']),xfact=PX,yfact=PX,origin=(0,0)))
    features.append((p.get('name',''),cl,g))
# Use the same physical field shown in the supplied figure. Object names changed
# during the accepted reconstruction; report the current name rather than reusing199.
anchor = shapely.Point(10100*PX,5950*PX)
cores = [(n,g) for n,cl,g in features if cl=='Islet']
name,core = min(cores,key=lambda x:x[1].distance(anchor))
assert name=='Islet_79'
region = shapely.box(9650*PX,5450*PX,10500*PX,6450*PX)
structure_colors={'IsletExpanded':'#1DBFC3','Islet':'#2784B8',
                  'Lymphatic':'#25A03E','Capillary':'#D8232A','Nerve':'#C832D0'}
structure_labels={'IsletExpanded':'20 µm outer boundary','Islet':'Endocrine object',
                  'Lymphatic':'Lymphatic','Capillary':'Capillary','Nerve':'Nerve'}
def finish(fig,stem):
    for suffix in ('png','pdf','svg'):
        fig.savefig(OUT/'figures'/f'{stem}.{suffix}',dpi=300,facecolor='white',bbox_inches='tight')
    plt.close(fig)
fig,ax=plt.subplots(figsize=(4.2,4.7))
saved=[]
for n,cl,g in features:
    if cl not in structure_colors or not g.intersects(region): continue
    g=g.intersection(region)
    draw_shape(ax,g,'none' if cl in ('Islet','IsletExpanded') else structure_colors[cl],
               structure_colors[cl],1.2 if cl in ('Islet','IsletExpanded') else .3)
    saved.append({'type':'Feature','properties':{'name':n,'classification':cl},'geometry':shapely.geometry.mapping(g)})
ax.set(xlim=(region.bounds[0],region.bounds[2]),ylim=(region.bounds[1],region.bounds[3]),
       xlabel='X position (µm)',ylabel='Y position (µm)',title=f'{name} · Donor 6533')
ax.set_aspect('equal');ax.spines[['top','right']].set_visible(False)
finish(fig,'slide01_image02')
(OUT/'data/slide01_structure_field.geojson').write_text(json.dumps({'type':'FeatureCollection','features':saved}))
fig,ax=plt.subplots(figsize=(2.9,2.5));ax.axis('off')
handles=[Line2D([0],[0],color=structure_colors[k],lw=2,label=structure_labels[k]) for k in structure_colors]
ax.legend(handles=handles,title='Structures',loc='upper left',frameon=False,fontsize=10,title_fontsize=11)
finish(fig,'slide01_image09')
cells=pd.read_parquet(ROOT/'output/h1_edge_corrected_20261024/analysis_cells_by_donor/6533.parquet')
cells=cells[cells.edge_analysis_eligible].copy()
points=shapely.points(cells[['x_um','y_um']].to_numpy())
cells=cells.loc[shapely.covers(core.buffer(20),points)].copy()
cells['in_core']=shapely.covers(core,shapely.points(cells[['x_um','y_um']].to_numpy()))
cells.to_csv(OUT/'data/slide01_image03_cells.csv',index=False)
fig,ax=plt.subplots(figsize=(4.4,4.0))
for ph,d in cells.groupby('phenotype',sort=True):
    ax.scatter(d.x_um,d.y_um,s=24,c=palette.get(ph,'#CCCCCC'),edgecolors='#555555',linewidths=.12,label=ph)
draw_shape(ax,core,'none','#2784B8',1.2,zorder=4)
draw_shape(ax,core.buffer(20),'none','#1DBFC3',1.0,zorder=3)
ax.set(xlabel='X position (µm)',ylabel='Y position (µm)',title=f'{name} · Cells')
ax.set_aspect('equal');ax.spines[['top','right']].set_visible(False)
finish(fig,'slide01_image03')
fig,ax=plt.subplots(figsize=(2.8,4));ax.axis('off')
handles=[Line2D([0],[0],marker='o',linestyle='none',color=palette.get(p,'#CCCCCC'),label=p,markersize=6)
         for p in sorted(cells.phenotype.unique())]
ax.legend(handles=handles,title='Cell phenotype',loc='upper left',frameon=False,fontsize=10,title_fontsize=11)
finish(fig,'slide01_image04')
features=pd.read_parquet(ROOT/'output/h1_edge_corrected_20261024/h1_object_features.parquet')
donors=features[['case_id','donor_status']].drop_duplicates().sort_values(['donor_status','case_id'])
assert len(donors)==20 and not donors.case_id.astype(str).eq('78').any()
donors.to_csv(OUT/'data/cohort.csv',index=False)
manifest={'dataset':str(APP),'donors':20,'donors_by_group':donors.groupby('donor_status',observed=True).size().to_dict(),
 'source_reference_region_pixels':[9650,5450,10500,6450], 'current_object':'6533_Islet_79',
 'old_reference_label':'6533_Islet_199','object_mapping':'Same physical field selected by source-figure coordinates; current reconstructed object name used.',
 'pixel_size_um':PX,'cells_in_display':len(cells),'continuous_INS_used':False,
 'source_micrographs':'Retained unchanged from supplied deck as illustrative images, without quantitative INS comparison or inferred donor identity.'}
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
(OUT/'legends.md').write_text('''# Context figure legends\n\nFigure 1. Accepted cohort and spatial organization. The cohort contains 20 donors (7 ND, 7 Aab+, 6 T1D). The structure and cell maps use the same physical field in donor 6533 as the reference illustration, with accepted current geometry and cell labels. The reconstructed object in that field is Islet_79; the historical figure called its earlier annotation Islet_199. Coordinates use micrometres (0.3774 µm per source pixel). The outer outline lies 20 µm from the endocrine-object boundary. Cell points represent unique edge-eligible cells within that expanded object. The schematic and image montage are retained from the source presentation. These illustrative panels contain no hypothesis test.\n\nFigure 4 images. The two source micrographs and their channel labels are retained unchanged from the supplied deck. No quantitative INS comparison is made from their display brightness. Acquisition identities and scale-bar lengths are not added where the source does not provide them.\n''')
print(json.dumps(manifest,indent=2))
