#!/usr/bin/env python3
"""Hash the new runtime and original input provenance; no deployment attestation."""
import argparse
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq
from common import sha256, write_json, require_cohort

def artifact(root,relative,startup=True):
    path=root/relative
    if path.is_dir():
        files=sorted(p for p in path.rglob('*') if p.is_file())
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=4) as pool:
            records=list(pool.map(lambda p:dict(path=str(p.relative_to(path)),bytes=p.stat().st_size,sha256=sha256(p)),files))
        import hashlib
        digest=hashlib.sha256('\n'.join(r['path']+'  '+r['sha256'] for r in records).encode()).hexdigest()
        return dict(relative_path=relative,type='tree',sha256=digest,bytes=sum(r['bytes'] for r in records),files=records,startup_verify=startup)
    result=dict(relative_path=relative,type='file',sha256=sha256(path),bytes=path.stat().st_size,startup_verify=startup)
    if relative.endswith('.parquet'):result['row_count']=pq.read_metadata(path).num_rows
    return result

def startup_verify_artifact(name, schema):
    """Keep startup checks bounded while retaining creation-time full hashes."""
    large_immutable = {'images', 'geometry', 'cells', 'islet_boundaries', 'provenance'}
    if str(schema) == '4':
        large_immutable.update({'typing', 'trajectory', 'spatial_tissue'})
    return name not in large_immutable

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-root',required=True);p.add_argument('--schema',choices=['2','3','4'],default='3');p.add_argument('--add-cell-artifact',action='store_true');a=p.parse_args();root=Path(a.input_root).resolve(strict=True)
    if a.add_cell_artifact:
        # Append a newly completed optional artifact to this candidate's already
        # verified manifest. Other inputs must remain immutable; this is not a
        # general refresh operation and never attests to changed existing files.
        manifest=json.loads((root/'manifest.json').read_text())
        if 'cells' in manifest['artifacts']:raise ValueError('Cell artifact is already declared')
        manifest['artifacts']['cells']=artifact(root,'cells',False)
        write_json(root/'manifest.json',manifest)
        capacity=json.loads((root/'capacity.json').read_text())
        capacity['runtime_bytes']+=manifest['artifacts']['cells']['bytes']
        capacity['recommended_two_releases_bytes']=2*capacity['runtime_bytes']
        capacity['candidate_logical_bytes']=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
        write_json(root/'capacity.json',capacity)
        return
    cohort=pd.read_csv(root/'cohort.csv',dtype=str)
    if len(cohort)!=21:raise ValueError('Wrong cohort')
    for path in [('trajectory/composition_features.csv' if a.schema=='4' else 'trajectory/features.csv'),'trajectory/coordinates.csv','donor_metadata.csv']:
        require_cohort(pd.read_csv(root/path,dtype={'case_id':str}),cohort.case_id)
    mapping={'core_area':'summary/core_area.parquet','targets_all':'summary/targets.parquet','markers_all':'summary/markers.parquet',
             'comp':'summary/islets.parquet','combined_pairs':'summary/combined_pairs.parquet','geometry':'geometry',
             'islet_spatial_lookup':'islet_spatial_lookup.csv','phenotype_rules':'phenotype_rules.csv','donor_metadata':'donor_metadata.csv',
             'images':'images','cohort':'cohort.csv','trajectory':'trajectory','cells':'cells'}
    # Exact source vertices support outline downloads without another rounding step.
    if (root/'islet_boundaries').exists():mapping['islet_boundaries']='islet_boundaries'
    if (root/'provenance').exists():mapping['provenance']='provenance'
    if a.schema=='2':mapping['research']='research'
    elif a.schema=='4':
        import tifffile
        import xml.etree.ElementTree as ET
        channel_records={}
        for image in sorted((root/'images').glob('*.ome.tif*')):
            with tifffile.TiffFile(image) as tif:
                if not tif.ome_metadata:raise ValueError('Embedded OME metadata required: '+image.name)
                xml=ET.fromstring(tif.ome_metadata)
                channels=[node.attrib.get('Name','') for node in xml.findall('.//{*}Pixels/{*}Channel')]
                if len(channels)!=len(set(channels)) or any(not name for name in channels):raise ValueError('Invalid OME channel panel: '+image.name)
                pixels=xml.find('.//{*}Pixels')
                if len(channels)!=int(pixels.attrib['SizeC']):raise ValueError('OME channel count mismatch: '+image.name)
                units={'µm':1.,'um':1.,'micrometer':1.,'mm':1000.,'nm':.001}
                row=cohort.loc[cohort.image.map(lambda p:Path(p).name).eq(image.name)].iloc[0]
                calibration=json.loads((root/'boundaries/sources'/str(row.case_id)/'source.json').read_text())
                channel_records[image.name]=dict(channel_names=channels,source='embedded OME XML',
                    pixel_size_x_um=float(pixels.attrib['PhysicalSizeX'])*units[pixels.attrib.get('PhysicalSizeXUnit','µm')],
                    pixel_size_y_um=float(pixels.attrib['PhysicalSizeY'])*units[pixels.attrib.get('PhysicalSizeYUnit','µm')],
                    analysis_pixel_size_um=calibration['pixel_size_um'],
                    acquired_ome_pixel_size_um=float(pixels.attrib['PhysicalSizeX'])*units[pixels.attrib.get('PhysicalSizeXUnit','µm')])
        write_json(root/'image_channels.json',dict(schema_version=1,images=channel_records))
        mapping['image_channels']='image_channels.json'
        if (root/'donors').is_dir():mapping['spatial_tissue']='donors'
        mapping['islet_embedding']='analysis/islet_embedding.parquet'
        mapping['islet_embedding_manifest']='analysis/islet_embedding_manifest.json'
        mapping['classifier_eligibility']='analysis/classifier_eligibility_manifest.json'
        mapping.pop('phenotype_rules')
        mapping['typing']='typing'
        files=['typing/catalog.json','typing/specification.json','analysis/classifier_eligibility_manifest.json','trajectory/feature_manifest.json',
               'trajectory/count_manifest.json','trajectory/marker_panel.json','trajectory/validation.json']
        for path in files:
            if not (root/path).is_file():raise ValueError('Missing schema-4 artifact: '+path)
        for name in ['composition_features','marker_summaries','counts','endocrine_counts','coordinates','edges','graph_layout']:
            if not any((root/'trajectory'/(name+suffix)).is_file() for suffix in ['.parquet','.csv']):
                raise ValueError('Missing schema-4 trajectory table: '+name)
        if (root/'research').exists():raise ValueError('Research reports belong outside the schema-4 runtime')
        # Large immutable tables/models get full verification at release creation
        # and deployment; startup checks the small contracts and fixed coordinates.
        mapping.update(typing_catalog='typing/catalog.json',typing_specification='typing/specification.json',
                       trajectory_feature_manifest='trajectory/feature_manifest.json',
                       trajectory_count_manifest='trajectory/count_manifest.json',trajectory_validation='trajectory/validation.json',
                       trajectory_coordinates='trajectory/coordinates.csv',trajectory_graph_layout='trajectory/graph_layout.csv')
    else:
        for name in ['counts.csv','endocrine_counts.csv','phenotype_availability.csv','marker_panel.json','count_manifest.json']:
            if not (root/'trajectory'/name).is_file():raise ValueError('Missing schema-3 trajectory artifact: '+name)
        if (root/'research').exists():raise ValueError('Research reports belong outside the schema-3 runtime')
    artifacts={}
    for name,rel in mapping.items():
        print('Hashing',name,flush=True)
        artifacts[name]=artifact(root,rel,startup_verify_artifact(name,a.schema))
    sources=pd.read_csv(root/'source_manifest.csv',dtype=str)
    paths=sorted(set(sources.cells)|set(sources.geometry)|set(sources.metadata))
    evidence=[dict(path=str(p),sha256=sha256(p),bytes=Path(p).stat().st_size) for p in paths]
    unavailable={'legacy_spatial_cell_distances':'not rebuilt for the expanded cohort; no frozen fallback'}
    if a.schema=='2':unavailable.update(phenotype_ratios='upstream gate validation incomplete; no certified denominator',hypotheses='per-test status in research/results.csv')
    else:
        counts=json.loads((root/'trajectory/count_manifest.json').read_text())
        if a.schema=='4':
            current=pd.read_parquet(root/'trajectory/counts.parquet') if (root/'trajectory/counts.parquet').exists() else pd.read_csv(root/'trajectory/counts.csv')
            availability=current['validated'] if 'validated' in current else current['available']
            if not availability.fillna(False).all():unavailable['phenotype_ratios']='see trajectory/counts for unsupported phenotype numerators and unresolved cells; classifier-core denominators remain independent of downstream typing'
            validation=json.loads((root/'trajectory/validation.json').read_text())
            if validation.get('status') == 'unavailable' or validation.get('coordinates_available') is False:
                unavailable['pseudotime']='outcome-blind fit withheld because the primary typed core/peri eligibility gate is unavailable for donors listed in trajectory/validation.json'
        elif counts['status']!='available':unavailable['phenotype_ratios']='see trajectory/phenotype_availability.csv and endocrine_counts.csv for donor/type-specific missing evidence'
        if a.schema!='4' and (root/'trajectory/quality_review.json').exists():
            quality=json.loads((root/'trajectory/quality_review.json').read_text())
            if quality['status']=='awaiting_reviewed_refit':
                unavailable['pseudotime']='Previous fit withdrawn after geometry correction; remaining donor boundary reviews and a new offline fit/validation are required'
                unavailable['islet_identity']='HPAP-078 owner-approved; remaining donor islet identity reviews pending'
            else:unavailable['islet_identity']='candidate islet identity requires review; existing coordinates are provisional'
    manifest=dict(schema_version=a.schema,created_at_utc=datetime.now(timezone.utc).isoformat(),builder='scripts/research/finalize_candidate.py',
                  scientific_acceptance=False,ufrc_accepted=False,expected_images=21,artifacts=artifacts,canonical_sources=evidence,
                  unavailable_analyses=unavailable,
                  research_results_in_app=False,
                  materialization_mode='hardlink raw inputs/images; newly written derived data; never mutate linked inputs')
    if a.schema=='4':
        manifest['measurement_definitions']=dict(
            trajectory_features='log1p_validated_base_type_count_per_classifier_defined_core_cell',
            spatial_embedding_features='same_typed_classifier_core_normalized_composition_as_primary_trajectory',
            islet_eligibility='no_image_edge_contact_and_at_least_one_canonical_classifier_core_cell',
            endocrine_denominator='all_canonical_cells_inside_endocrine_pixel_classifier_core_counted_once',
            heatmap_positive='log1p_marker_positive_count_per_classifier_defined_core_cell',
            heatmap_expression='mean_log1p_corrected_intensity_within_typed_population',
            heatmap_bins='equal_islet_mean_within_donor_then_equal_donor_mean')
        manifest['provenance_hashes']={name:artifacts[name]['sha256'] for name in
            ['typing_specification','typing_catalog','trajectory_feature_manifest']}
        manifest['provenance_hashes']['geometry']=artifacts['geometry']['sha256']
        manifest['provenance_hashes']['typing_models']=artifacts['typing']['sha256']
        manifest['provenance_hashes']['islet_embedding_manifest']=artifacts['islet_embedding_manifest']['sha256']
        manifest['provenance_hashes']['classifier_eligibility']=artifacts['classifier_eligibility']['sha256']
    write_json(root/'manifest.json',manifest)
    runtime_bytes=sum(x['bytes'] for x in artifacts.values())
    build_report = root/'build_report.json'
    analysis_seconds = 0.0
    if build_report.exists():
        analysis_seconds = sum(d.get('seconds', 0.0) for d in json.loads(build_report.read_text()).get('donor_reports', []))
    usage=dict(runtime_bytes=runtime_bytes,images_bytes=artifacts['images']['bytes'],
               candidate_logical_bytes=sum(p.stat().st_size for p in root.rglob('*') if p.is_file()),
               recommended_two_releases_bytes=2*runtime_bytes,
               cpu=platform.processor(),trajectory=json.loads((root/'trajectory/validation.json').read_text()),
               analysis_seconds=analysis_seconds,
               note='Logical sizes include hardlinked source files; production copy needs full runtime bytes. No image/container acceptance is implied.')
    write_json(root/'capacity.json',usage)

if __name__=='__main__':main()
