#!/usr/bin/env python3
"""Create a fresh full native QuPath review project from final identity tables.

Every detection receives its final display class or exclusion reason. The copy retains
Tissue/endocrine geometry and adds the four endpoint structure classes. Assignment CSVs,
logs and the completion manifest live in a separate fresh directory.

    python scripts/uniform/review_project.py --run <final-run> \\
      --project <new-project> --assignments <new-assignment-directory>
"""
from __future__ import annotations

import argparse
import colorsys
import csv
import json
import os
import shutil
import hashlib
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent / 'qupath-review-20261027'
QUPATH = '/home/smith6jt/QuPath/bin/QuPath'
SUP = {'+': '⁺', '-': '⁻'}
FIXED = {'Beta cell': 'Endocrine Beta', 'Alpha cell': 'Endocrine Alpha', 'Delta cell': 'Endocrine Delta',
         'CHGA+ hormone-negative': 'Endocrine CHGA⁺ hormone-negative', 'Endocrine': 'Endocrine',
         'CD8a Tcell': 'T Cell CD8a⁺', 'CD4 Tcell': 'T Cell CD4⁺', 'T cell': 'T Cell', 'B cell': 'B Cell',
         'Immune': 'Immune', 'CD45+': 'CD45⁺', 'Hematopoietic CD45+': 'Hematopoietic CD45⁺', 'Hematopoietic CD45-': 'Hematopoietic CD45⁻', 'CD45+ (top marker)': 'CD45⁺ (top marker)',
         'Exocrine': 'Exocrine'}
COLORS = {'T Cell CD8a⁺': '#E41A1C', 'T Cell CD4⁺': '#FF7F00', 'T Cell': '#F781BF', 'B Cell': '#FFFF33',
          'Immune': '#1B9E77', 'CD45⁺': '#1B9E77', 'Hematopoietic CD45⁺': '#1B9E77',
          'Hematopoietic CD45⁻': '#74C476', 'Endocrine Beta': '#9ECAE1', 'Endocrine Alpha': '#A1D99B',
          'Endocrine Delta': '#BCBDDC', 'Endocrine CHGA⁺ hormone-negative': '#8D6E63', 'Endocrine': '#BDBDBD',
          'Exocrine': '#E0E0E0', 'Excluded - cell QC': '#555555', 'Excluded - exocrine inside islet border': '#A0522D',
          'Excluded - hormone under 2x runner-up': '#9932CC', 'Excluded - single without usable hormone': '#6A5ACD',
          'Excluded - ColIV top marker': '#8B7D6B', 'Excluded - lineage ambiguous': '#B8860B',
          'Excluded - doublet suspected': '#FF1493',
          'CD3e\u207a': '#C9A3A3', 'CD20\u207a': '#C9C2A3', 'CD68\u207a': '#A3C9C2', 'CD163\u207a': '#A3B5C9',
          'ECAD\u207a': '#B5A3C9', 'panCK\u207a': '#C9A3BD',
          'B3TUBB⁺': '#B8A47E', 'CD31⁺': '#C98B8B', 'CD34⁺': '#D4A5A5', 'CD4⁺': '#9FB7C9',
          'CD8a⁺': '#8FA8D0', 'GCG⁺': '#B5C99A', 'INS⁺': '#A9C4D9', 'LGALS3⁺': '#C3B1D6',
          'PDPN⁺': '#D9B38C', 'SMA⁺': '#A3A3C2', 'SST⁺': '#C2C2A3', 'VIM⁺': '#B7C9B7',
          'CD45⁺ (top marker)': '#8CC2B5',
          'Endothelium Blood Vessel': '#D7301F', 'Endothelium Lymphatic': '#41B6C4', 'Neural': '#FEC44F',
          # final typing run (2026-10-28): cells no marker names keep their parent object's class; other marker names
          'Tissue cell - no marker above level': '#3A3A3A', 'Islet cell - no marker above level': '#5A4A6A',
          'EndocrineCluster cell - no marker above level': '#5A4A6A', 'EndocrineSingle cell - no marker above level': '#5A4A6A',
          'Tissue cell - endocrine marker only': '#7A6A8A', 'Islet cell - endocrine marker only': '#9A7AAA',
          'EndocrineCluster cell - endocrine marker only': '#9A7AAA', 'EndocrineSingle cell - endocrine marker only': '#9A7AAA',
          'BCatenin⁺': '#CDBE9A', 'CD35⁺': '#C2D69A', 'GAP43⁺': '#C8B068', 'HLADR⁺': '#7FB3A0', 'PGP9.5⁺': '#B08FC0',
          'CHGA⁺': '#B0A0C8', 'Islet cell - no hormone above level': '#6E8FA6',
          'EndocrineCluster cell - no hormone above level': '#6E8FA6', 'EndocrineSingle cell - no hormone above level': '#6E8FA6'}
START_STATUS = {'no_marker_above_background': 'no marker above level', 'endocrine_marker_only': 'endocrine marker only',
                'no_hormone_above_level': 'no hormone above level',
                # no cell is named by a marker (owner 2026-10-28, "Get rid of naming by top marker logic")
                'no_lineage_named': 'no lineage named', 'shared_marker_only': 'shared marker only',
                'no_endocrine_marker_above_level': 'no endocrine marker above level',
                'tied_highest_hormones': 'tied highest hormones'}
for _start in ('Tissue', 'Islet', 'EndocrineCluster', 'EndocrineSingle'):
    for _st, _col in (('no lineage named', '#6B5B3A'), ('shared marker only', '#4A5A4A'),
                      ('no endocrine marker above level', '#6E8FA6')):
        COLORS.setdefault(f'{_start} cell - {_st}', _col)
# classes earlier typing imports added; a class of this kind that no cell carries now is taken off the project's list
TYPING_CLASS = ('Macrophage ', 'Excluded - ', 'Conflict ', 'T Cell', 'B Cell', 'Endocrine ', 'Exocrine', 'Immune',
                'Doublet suspected', 'Endothelial', 'Stromal', 'Neural', 'Hematopoietic')
IDENTITY_META = ['start_class', 'compartment', 'compartment_status', 'lineage', 'lineage_status', 'subtype_status',
                 'status', 'named_by', 'tags', 'confirmed', 'confirming_markers', 'established_label']


def identity_tables(run: Path, t: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(typed rows, excluded rows) of one image from the final typing run's typing/identity (not assembled)."""
    d = pd.read_parquet(run / 'typing' / 'identity' / f'{t}.parquet', columns=['object_id', 'label', 'excluded_reason']
                        + IDENTITY_META)
    typed = d[d['label'] != 'Excluded'].copy()
    name = typed['label'].map(display)
    kept = (typed['status'] != '') & (typed['label'] == typed['start_class'])   # the cell kept its start class
    name[kept] = typed.loc[kept, 'start_class'] + ' cell - ' + typed.loc[kept, 'status'].map(START_STATUS)
    if name.isna().any():
        raise SystemExit(f'{t}: unmapped status')
    typed['confirmed'] = typed['confirmed'].map(lambda v: '' if pd.isna(v) else ('yes' if v else 'no'))
    out = pd.DataFrame({'image': t, 'object_id': typed['object_id'], 'phenotype': name, 'phenotype_label': typed['label']})
    for c in IDENTITY_META:
        out[c] = typed[c].astype(str).str.replace(',', ';')
    ex = d[d['label'] == 'Excluded'][['object_id', 'excluded_reason']]
    qc = pd.read_parquet(run / 'typing' / 'redsea' / f'{t}.parquet', columns=['object_id', 'qc_analysis_eligible', 'qc_reasons'])
    qc = qc[~qc['qc_analysis_eligible']]
    cls = ex['excluded_reason'].map(EXCLUDED)
    if cls.isna().any():
        raise SystemExit(f'{t}: unmapped exclusion reasons {sorted(set(ex.loc[cls.isna(), "excluded_reason"]))}')
    excl = pd.concat([pd.DataFrame({'object_id': ex['object_id'], 'class': cls, 'reason': ex['excluded_reason']}),
                      pd.DataFrame({'object_id': qc['object_id'], 'class': 'Excluded - cell QC',
                                    'reason': qc['qc_reasons'].astype(str).str.replace(',', ';')})])
    return out, excl
EXCLUDED = {'exocrine_inside_islet_border': 'Excluded - exocrine inside islet border',
            'hormone_unresolved_margin': 'Excluded - hormone under 2x runner-up',
            'endocrine_single_hormone_margin_below_2x': 'Excluded - hormone under 2x runner-up',
            'endocrine_single_no_usable_hormone': 'Excluded - single without usable hormone',
            'top_marker_ColIV': 'Excluded - ColIV top marker',
            'lineage_ambiguous': 'Excluded - lineage ambiguous', 'doublet_suspected': 'Excluded - doublet suspected'}
COLS = ['object_id', 'phenotype', 'phenotype_label', 'cell_class', 'refined_class', 'immune_definition',
        'cell_region', 'islet_name', 'typing_route']


def display(label: str) -> str:
    if label in FIXED:
        return FIXED[label]
    if 'Macrophage_' in label:
        head, tail = label.split('Macrophage_')
        return 'Macrophage ' + ''.join(SUP.get(c, c) for c in head + tail)
    if label.endswith('+') and ' ' not in label:
        return label[:-1] + '⁺'
    return label


def color(name: str) -> str:
    if name in COLORS:
        return COLORS[name]
    if name.startswith('Conflict '):                          # canary conflicts (final typing run, 2026-10-28)
        return '#B8860B'
    if name == 'Doublet suspected':
        return '#FF1493'
    if name.startswith('Macrophage '):                       # as import_review.groovy colors them
        p, m = SUP['+'], SUP['-']
        hue = 0.58 if f'CD68{p}CD163{p}' in name else (0.50 if f'CD68{p}' in name else 0.75)
        k = (f'LGALS3{m}' in name) + (f'CD4{m}' in name)
        r, g, b = colorsys.hsv_to_rgb(hue, 1 - 0.3 * k, 1 - 0.15 * k)
        return '#%02X%02X%02X' % (round(r * 255), round(g * 255), round(b * 255))
    return '#808080'


def argb(hexcolor: str) -> int:
    v = (0xFF << 24) | int(hexcolor[1:], 16)
    return v - (1 << 32)


def main() -> int:
    """Create a fresh project so an open older review cannot overwrite this import."""
    ap = argparse.ArgumentParser(description='Build a fresh full QuPath review from final typing')
    ap.add_argument('--run',type=Path,required=True)
    ap.add_argument('--project',type=Path,required=True)
    ap.add_argument('--assignments',type=Path,required=True)
    ap.add_argument('--workers',type=int,default=2)
    a=ap.parse_args();run=a.run.resolve();project=a.project.resolve();out=a.assignments.resolve()
    if project.exists() or out.exists():
        raise FileExistsError('Use fresh project and assignment paths; existing review projects are never overwritten')
    rows=list(csv.DictReader((REPO/'scripts/uniform/cohort.csv').open()))
    # Build and reconcile every assignment before making the native project copy.
    tables={}; all_classes=set(); summary=[]
    for row in rows:
        tag=row['image_name'];typed,excluded=identity_tables(run,tag)
        ids=pd.read_parquet(run/'cells'/f'{tag}.parquet',columns=['Object ID'])['Object ID']
        assigned=pd.concat([typed.object_id,excluded.object_id],ignore_index=True)
        if assigned.duplicated().any() or set(assigned)!=set(ids):
            raise ValueError(tag+': assignments do not cover every detection exactly once')
        tables[tag]=(typed,excluded);all_classes |= set(typed.phenotype)|set(excluded['class'])
        summary.append({'image':tag,'detections':len(ids),'analyzed':len(typed),'excluded':len(excluded)})
    project.parent.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True)
    shutil.copytree(run/'qupath',project)
    (out/'structures').mkdir();(out/'logs').mkdir()
    structure_colors={'Capillary':[255,45,45],'Nerve':[255,215,0],'Lymphatic':[230,153,77],'Collagen IV':[216,45,189]}
    classes_path=project/'classifiers'/'classes.json'
    data=json.loads(classes_path.read_text()) if classes_path.exists() else {'pathClasses':[]}
    kept=[c for c in data['pathClasses'] if c['name'] in ('Tissue','Islet','EndocrineCluster','EndocrineSingle','pLN')]
    classes={c['name']:c for c in kept}
    classes.update({c:{'name':c,'color':argb(color(c))} for c in sorted(all_classes) if c not in classes})
    classes.update({c:{'name':c,'color':argb('#%02x%02x%02x'%tuple(rgb))} for c,rgb in structure_colors.items()})
    data['pathClasses']=list(classes.values())
    classes_path.parent.mkdir(exist_ok=True);classes_path.write_text(json.dumps(data,indent=2,ensure_ascii=False))
    for tag,(typed,excluded) in tables.items():
        typed.to_csv(out/f'{tag}.csv',index=False);excluded.to_csv(out/f'{tag}.excluded.csv',index=False)
        for cls,rgb in structure_colors.items():
            src=run/'structures'/cls/f'{cls}__{tag}.geojson';obj=json.loads(src.read_text())
            for f in obj['features']:
                p=f.setdefault('properties',{});p['classification']={'name':cls,'color':rgb};p['isLocked']=True
            (out/'structures'/src.name).write_text(json.dumps(obj))
    pd.DataFrame(summary).to_csv(out/'summary_by_image.csv',index=False)
    script=Path(__file__).parent/'qupath'/'import_final_review.groovy'
    def import_one(tag):
        log=out/'logs'/f'{tag}.log'
        with log.open('w') as fh:
            rc=subprocess.run([QUPATH,'script','-p',str(project/'project.qpproj'),'-i',tag,'-a',str(out),str(script)],
                stdout=fh,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,
                env={**os.environ,'JAVA_TOOL_OPTIONS':'-Xmx64g -XX:ActiveProcessorCount=12'}).returncode
        text=log.read_text();ok=rc==0 and 'in neither 0' in text and f'saved {tag}' in text
        print(tag,'review imported' if ok else 'REVIEW FAILED',flush=True)
        return {'image':tag,'pass':ok,'log':str(log)}
    with ThreadPoolExecutor(a.workers) as pool:results=list(pool.map(import_one,tables))
    record={'run':str(run),'project':str(project),'images':results,'status':'complete' if all(r['pass'] for r in results) else 'failed'}
    (out/'review_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    (project/'REVIEW.txt').write_text('Current final cell typing, all detections accounted for; Tissue/endocrine objects and four structure classes.\nSource run: '+str(run)+'\nAssignments: '+str(out)+'\n')
    return 0 if record['status']=='complete' else 1

if __name__=='__main__':sys.exit(main())
