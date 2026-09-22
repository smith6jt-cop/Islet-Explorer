"""Independent schema-4 validation of native cell universes and typed endpoints."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from common import MARKERS, require_cohort, sha256, write_json
from build_runtime_schema4 import read_table, flags


def validate_structural_areas(frame):
    keys=['case_id','islet_id','compartment','structure'] + (['band_um'] if 'band_um' in frame else [])
    if frame.duplicated(keys).any():raise ValueError('Duplicate structural compartments')
    available=flags(frame.available)
    if frame.loc[~available,'positive_area_um2'].notna().any():
        raise ValueError('Unavailable mask has a numeric structural area')
    positive=frame.loc[available,'positive_area_um2']; valid=frame.loc[available,'valid_area_um2']
    if not (np.isfinite(positive)&np.isfinite(valid)&(positive>=0)&(valid>0)&(positive<=valid+1e-6)).all():
        raise ValueError('Invalid physical structural area or denominator')


def validate_count_arithmetic(counts, endocrine, base_types):
    if counts.duplicated(['case_id','islet_id','compartment','phenotype']).any():raise ValueError('Duplicate typed counts')
    if endocrine.islet_id.duplicated().any():raise ValueError('Duplicate endocrine denominators')
    available=flags(counts.validated if 'validated' in counts else counts.available)
    if counts.loc[~available,'count'].notna().any():raise ValueError('Unavailable typed counts must remain missing')
    if (counts.loc[available,'count']<0).any():raise ValueError('Negative typed count')
    den=endocrine.set_index('islet_id').endocrine_count
    expected=counts.islet_id.map(den).to_numpy(float)
    raw=np.divide(counts['count'].to_numpy(float),expected,out=np.full(len(counts),np.nan),where=np.isfinite(expected)&(expected>0))
    if not np.allclose(raw,counts.raw_ratio,equal_nan=True):raise ValueError('Typed numerator/denominator ratio mismatch')
    if not np.allclose(np.log1p(raw),counts.plotted_value,equal_nan=True):raise ValueError('Typed plotted ratio mismatch')
    base=counts[counts.phenotype.isin(base_types)]
    for _,group in base.groupby(['case_id','islet_id','compartment']):
        if len(group)!=len(base_types):raise ValueError('Incomplete base-identity count panel')
        if group['count'].notna().all() and not np.isclose(group['count'].sum(),group.n_typed.iloc[0]):
            raise ValueError('Base identities are not a disjoint census of typed cells')


def assert_parquet_csv_equal(parquet_path, csv_path):
    """Compare mirror tables while treating blank strings as serialized NA."""
    left=pd.read_parquet(parquet_path)
    right=pd.read_csv(csv_path,dtype={'case_id':str})
    for column in left.columns:
        if pd.api.types.is_object_dtype(left[column]) or pd.api.types.is_string_dtype(left[column]):
            left[column]=left[column].replace('',np.nan)
            right[column]=right[column].replace('',np.nan)
    pd.testing.assert_frame_equal(left,right,check_dtype=False)


def validate(root, report_root=None):
    root=Path(root).resolve(strict=True)
    cohort=pd.read_csv(root/'cohort.csv',dtype={'case_id':str})
    if len(cohort)!=21 or cohort.donor_status.value_counts().to_dict()!={'ND':7,'Aab+':7,'T1D':7}:
        raise ValueError('Expected exactly seven donors in each status group')
    if cohort.loc[cohort.donor_source.eq('HPAP'),'case_id'].tolist()!=['78']:raise ValueError('HPAP provenance mismatch')
    actual={p.name for p in (root/'images').glob('*.ome.tif*')}
    if actual!={Path(p).name for p in cohort.image}:raise ValueError('Image cohort mismatch')
    inputs=json.loads((root/'preprocessing/inputs.json').read_text())
    records={str(row['case_id']):row for row in inputs['donors']}
    if set(records)!=set(cohort.case_id):raise ValueError('Native source universe omits or mixes donors')
    islets=pd.read_parquet(root/'analysis/islets.parquet');islets['case_id']=islets.case_id.astype(str)
    membership=pd.read_parquet(root/'analysis/memberships.parquet');membership['case_id']=membership.case_id.astype(str)
    if 'band_um' in membership:membership=membership[membership.band_um.eq(20)]
    if membership.duplicated(['case_id','object_id','islet_id','compartment']).any():raise ValueError('Duplicate region membership')
    if set(membership.islet_id)-set(islets.islet_id):raise ValueError('Membership references obsolete islet geometry')
    require_cohort(islets,cohort.case_id)
    catalog=json.loads((root/'typing/catalog.json').read_text())
    base=[x if isinstance(x,str) else x['name'] for x in catalog['base_types']]
    endocrine_components=catalog['endocrine_components']
    if len(endocrine_components)!=4:raise ValueError('Expanded endocrine union must declare four identity components')
    counts=read_table(root,'trajectory/counts');counts['case_id']=counts.case_id.astype(str)
    endocrine=read_table(root,'trajectory/endocrine_counts');endocrine['case_id']=endocrine.case_id.astype(str)
    features=read_table(root,'trajectory/composition_features');features['case_id']=features.case_id.astype(str)
    coordinates=read_table(root,'trajectory/coordinates');coordinates['case_id']=coordinates.case_id.astype(str)
    for frame in [counts,endocrine,features,coordinates]:
        require_cohort(frame,cohort.case_id)
        if set(frame.islet_id)!=set(islets.islet_id):raise ValueError('Derived table does not match current islet geometry')
    validate_count_arithmetic(counts,endocrine,base)
    totals=dict(raw_objects=0,archival_objects=0,incomplete_raw_objects=0,wholly_unmeasured=0,corrected_unavailable=0,validated_cells=0)
    census=[]
    for cid in cohort.case_id:
        record=records[cid]
        raw=pd.read_parquet(root/'raw'/f'{cid}-cells.parquet')
        corrected=pd.read_parquet(root/'corrected'/f'{cid}-cells.parquet')
        typed=pd.read_parquet(root/'typing/cells'/f'{cid}.parquet')
        archived=pd.read_parquet(root/'provenance/archival_raw'/f'{cid}-cells.parquet',columns=['object_id'])
        if len(raw)!=record['canonical_rows'] or len(archived)!=record['archived_rows']:raise ValueError(f'{cid}: native universe changed')
        for frame in [raw,corrected,typed]:
            if frame.object_id.duplicated().any() or set(frame.case_id.astype(str))!={cid}:raise ValueError(f'{cid}: duplicate or cross-donor cells')
        if set(raw.object_id)!=set(corrected.object_id) or set(raw.object_id)!=set(typed.object_id):raise ValueError(f'{cid}: corrected/typed objects were discarded')
        finite=np.isfinite(raw[MARKERS]).all(axis=1)
        if not np.array_equal(finite,flags(raw.marker_measurements_complete)):raise ValueError(f'{cid}: raw completeness differs')
        if int((~finite).sum())!=record['incomplete_raw_cells']:raise ValueError(f'{cid}: raw exclusion accounting differs')
        measured=np.isfinite(corrected[MARKERS]).all(axis=1)
        if not np.array_equal(measured,flags(corrected.marker_measurements_complete)):raise ValueError(f'{cid}: correction completeness differs')
        if (corrected.loc[measured,MARKERS]<0).any().any():raise ValueError(f'{cid}: negative corrected intensity')
        if corrected.loc[~measured,'measurement_exclusion_reason'].fillna('').eq('').any():raise ValueError(f'{cid}: unexplained measurement exclusion')
        typed=typed.merge(corrected[['object_id','marker_measurements_complete']],on='object_id',suffixes=('','_corrected'),validate='one_to_one')
        eligible=typed.classification_status.eq('validated') & flags(typed.projection_supported)
        for column in ['segmentation_quality_excluded','technical_region_excluded']:
            if column in typed:eligible &= ~flags(typed[column])
        if set(typed.loc[eligible,'base_type'])-set(base):raise ValueError(f'{cid}: unknown validated identity')
        if not flags(typed.loc[eligible,'marker_measurements_complete']).all():raise ValueError(f'{cid}: unmeasured cell was classified')
        typed['validated_identity']=eligible
        typed['endocrine_member']=eligible & flags(typed.is_endocrine)
        if not (typed.loc[eligible,'endocrine_member']==typed.loc[eligible,'base_type'].isin(endocrine_components)).all():
            raise ValueError(f'{cid}: endocrine union is not the four declared identities')
        members=membership[membership.case_id.eq(cid)]
        joined=members.merge(typed[['object_id','base_type','validated_identity','endocrine_member']],on='object_id',how='left',validate='many_to_one',indicator=True)
        if not joined._merge.eq('both').all():raise ValueError(f'{cid}: missing typed membership')
        # The endocrine pixel classifier defines the core.  Every canonical
        # core member enters the denominator exactly once; downstream identity
        # support only controls phenotype numerators.
        expected=joined[joined.compartment.eq('core')].groupby('islet_id').size()
        observed=endocrine[endocrine.case_id.eq(cid)]
        available=observed.endocrine_count.notna()
        if not np.allclose(observed.loc[available,'endocrine_count'],observed.loc[available,'islet_id'].map(expected).fillna(0)):
            raise ValueError(f'{cid}: classifier-core denominator differs from unique cell membership')
        expected_counts=joined[joined.validated_identity].groupby(['islet_id','compartment','base_type']).size()
        current=counts[counts.case_id.eq(cid)&counts.phenotype.isin(base)]
        for row in current[current['count'].notna()].itertuples():
            if row.count!=expected_counts.get((row.islet_id,row.compartment,row.phenotype),0):raise ValueError(f'{cid}: phenotype census differs')
        totals['raw_objects']+=len(raw);totals['archival_objects']+=len(archived)
        totals['incomplete_raw_objects']+=int((~finite).sum());totals['wholly_unmeasured']+=int((~np.isfinite(raw[MARKERS]).any(axis=1)).sum())
        totals['corrected_unavailable']+=int((~measured).sum());totals['validated_cells']+=int(eligible.sum())
        census.append({'case_id':cid,'canonical_cells':len(raw),'validated_cells':int(eligible.sum()),'islets':int(islets.case_id.eq(cid).sum())})
    manifest=json.loads((root/'trajectory/feature_manifest.json').read_text())
    columns=manifest['columns']
    if manifest.get('representation')!='typed_composition' or not manifest.get('typing_required'):raise ValueError('Trajectory is not downstream of typed composition')
    eligibility_manifest=json.loads((root/'analysis/classifier_eligibility_manifest.json').read_text())
    eligibility=manifest.get('eligibility',{})
    if (manifest.get('minimum_classifier_core_cells')!=1 or len(eligibility.get('all_of',[]))!=2
            or eligibility.get('downstream_typing_can_revoke_eligibility') is not False):
        raise ValueError('Classifier-defined trajectory eligibility rules are not recorded')
    if eligibility_manifest.get('definition',{}).get('downstream_typing_can_revoke_eligibility') is not False:
        raise ValueError('Classifier eligibility was overridden by downstream typing')
    if eligibility_manifest.get('islets_sha256')!=sha256(root/'analysis/islets.parquet'):
        raise ValueError('Classifier eligibility geometry provenance changed')
    classifier_columns={'classifier_defined_core','image_edge_contact','image_edge_clearance_um',
                        'classifier_core_cells','classifier_eligible','classifier_exclusion_reason'}
    if not classifier_columns.issubset(islets):raise ValueError('Analysis islets lack classifier eligibility fields')
    expected_islet_eligibility=(~flags(islets.image_edge_contact)) & islets.classifier_core_cells.gt(0)
    if not flags(islets.classifier_defined_core).all() or not np.array_equal(flags(islets.classifier_eligible),expected_islet_eligibility):
        raise ValueError('Analysis islet eligibility differs from classifier geometry and core-cell membership')
    if eligibility_manifest.get('eligible_islets')!=int(expected_islet_eligibility.sum()):
        raise ValueError('Classifier eligibility manifest census differs')
    eligible_by_donor=features.loc[features.eligible].groupby('case_id').size()
    if set(eligible_by_donor.index.astype(str))!=set(cohort.case_id) or int(eligible_by_donor.min())<2:
        raise ValueError('Every donor must contribute at least two eligible islets')
    recorded_by_donor={str(key):int(value) for key,value in manifest.get('eligible_islets_by_donor',{}).items()}
    observed_by_donor={str(key):int(value) for key,value in eligible_by_donor.items()}
    if recorded_by_donor!=observed_by_donor:raise ValueError('Recorded donor eligibility census differs')
    if set(columns)!={f'{c}__count__{t}' for c in ['core','peri'] for t in base}:raise ValueError('Unexpected trajectory feature matrix')
    primary_columns=manifest.get('primary_columns',[c for c in columns if c.startswith('core__')])
    if not set(primary_columns).issubset(columns):raise ValueError('Primary trajectory columns are not in the feature manifest')
    expected_eligibility=(~flags(features.image_edge_contact)) & features.raw_cells.gt(0)
    if not np.array_equal(flags(features.eligible),expected_eligibility):
        raise ValueError('Eligibility must be only no image-edge contact plus at least one classifier-core cell')
    if not np.array_equal(features.endocrine_count.to_numpy(float),features.raw_cells.to_numpy(float)):
        raise ValueError('Classifier-core cell census is not the endocrine denominator')
    supported=flags(features.eligible)&flags(features.trajectory_feature_supported)
    observed_support=flags(features.eligible)&np.isfinite(features[primary_columns]).all(axis=1)
    if not np.array_equal(supported,observed_support):raise ValueError('Trajectory feature support differs from available typed composition')
    if manifest.get('feature_supported_islets')!=int(supported.sum()):raise ValueError('Feature-supported islet census differs')
    if manifest['features_sha256']!=sha256(root/'trajectory/composition_features.parquet'):raise ValueError('Trajectory feature provenance changed')
    embedding=pd.read_parquet(root/'analysis/islet_embedding.parquet');embedding['case_id']=embedding.case_id.astype(str)
    embedding_manifest=json.loads((root/'analysis/islet_embedding_manifest.json').read_text())
    expected_embedding_ids=set(features.loc[supported,'islet_id'])
    if set(embedding.islet_id)!=expected_embedding_ids or embedding.islet_id.duplicated().any():raise ValueError('Spatial embedding does not cover each feature-supported eligible islet exactly once')
    require_cohort(embedding,cohort.case_id)
    if embedding_manifest.get('representation')!='typed_composition' or not embedding_manifest.get('outcomes_excluded'):
        raise ValueError('Spatial embedding is not outcome-blind typed composition')
    if not set(embedding_manifest.get('feature_columns',[])).issubset(primary_columns):raise ValueError('Spatial embedding uses undeclared trajectory features')
    if 'whole-islet marker averages' not in embedding_manifest.get('excluded_feature_families',[]):raise ValueError('Spatial embedding does not exclude whole-islet marker averages')
    if embedding_manifest.get('embedding_sha256')!=sha256(root/'analysis/islet_embedding.parquet'):raise ValueError('Spatial embedding provenance changed')
    if (embedding_manifest.get('classifier_eligible_islets')!=int(flags(features.eligible).sum())
            or embedding_manifest.get('eligible_islets')!=int(flags(features.eligible).sum())
            or embedding_manifest.get('feature_supported_islets')!=len(embedding)
            or embedding_manifest.get('embedded_islets')!=len(embedding)
            or embedding_manifest.get('donors')!=21):raise ValueError('Spatial embedding manifest census mismatch')
    embed_columns=['leiden_umap_1','leiden_umap_2','leiden_0.3','leiden_0.5','leiden_0.8','leiden_1.0']
    if not set(embed_columns).issubset(embedding) or not np.isfinite(embedding[['leiden_umap_1','leiden_umap_2']]).all().all():raise ValueError('Spatial UMAP or Leiden outputs are incomplete')
    markers=read_table(root,'trajectory/marker_summaries')
    if set(markers.marker)!=set(MARKERS):raise ValueError('Incomplete heatmap panel')
    if markers.duplicated(['case_id','islet_id','compartment','population','marker']).any():raise ValueError('Duplicate marker summary')
    if markers.loc[~flags(markers.positive_available),'positive_count'].notna().any():raise ValueError('Unavailable positivity was replaced by a number')
    structural=pd.read_parquet(root/'analysis/spatial_measurements.parquet');validate_structural_areas(structural)
    regions=pd.read_parquet(root/'analysis/regions.parquet')
    if set(regions.band_um)!={10,20,30}:raise ValueError('Missing spatial sensitivity bands')
    status=json.loads((root/'trajectory/validation.json').read_text())
    supported=flags(coordinates.supported)&coordinates.pseudotime_raw.notna()
    if not supported.any():raise ValueError('All-missing trajectory does not satisfy local delivery acceptance')
    if status.get('bootstrap_completed')!=200 or status.get('lodo_completed')!=21:raise ValueError('Incomplete trajectory validation')
    if coordinates.islet_id.duplicated().any():raise ValueError('Duplicate trajectory coordinates')
    for name in ['coordinates','composition_features','counts','endocrine_counts','edges']:
        parquet_path=root/'trajectory'/f'{name}.parquet';csv=parquet_path.with_suffix('.csv')
        if parquet_path.exists() and csv.exists():assert_parquet_csv_equal(parquet_path,csv)
    eligible_ids=set(features.loc[features.eligible,'islet_id'])
    covered_ids=set(membership.loc[membership.band_um.eq(20),'islet_id']) if 'band_um' in membership else set(membership.islet_id)
    if not eligible_ids.issubset(covered_ids):raise ValueError('Eligible islet drill-down membership missing')
    cell_files={p.name.removesuffix('-cells.parquet'):p for p in (root/'cells').glob('*-cells.parquet')}
    if set(cell_files)!=set(cohort.case_id):raise ValueError('Per-donor cell partitions do not match cohort')
    for cid in cohort.case_id:
        # ``schema.names`` flattens list-valued Arrow fields to their child
        # name (for example, phenotype_states becomes ``element``).  Validate
        # the top-level columns that DuckDB and pandas expose to the app.
        cell_columns=set(pq.ParquetFile(cell_files[cid]).schema_arrow.names)
        typing_columns=set(pq.ParquetFile(root/'typing/cells'/f'{cid}.parquet').schema_arrow.names)
        if not set(MARKERS+['object_id','X_centroid','Y_centroid']).issubset(cell_columns):raise ValueError(f'{cid}: cell partition lacks coordinates or markers')
        if not {'object_id','base_type','phenotype_state','phenotype_states','classification_status'}.issubset(typing_columns):raise ValueError(f'{cid}: typing partition lacks drill-down identities')
    if (root/'research').exists():raise ValueError('Research reports must remain outside runtime')
    result=dict(schema_version=4,donors=21,shared_markers=31,islets=len(islets),eligible=int(features.eligible.sum()),
                feature_supported_islets=int(flags(features.trajectory_feature_supported).sum()),embedded_islets=len(embedding),supported_islets=int(supported.sum()),
                trajectory_status=status['status'],native_cell_census=census,**totals,passed=True,scientific_acceptance=False,
                canonical_count_reconciliation=True,research_not_in_runtime=True,validation_script_sha256=sha256(__file__))
    if report_root is not None:
        reports=Path(report_root).resolve()
        if reports.is_relative_to(root):raise ValueError('Reports are inside runtime')
        results=pd.read_csv(reports/'results.csv')
        if set(results.id)!={f'H{i}' for i in range(1,11)}:raise ValueError('Primary result slots differ')
        if results.loc[results.status.eq('unavailable'),'p_value'].notna().any():raise ValueError('Unavailable test has a p-value')
        result['estimated_primary_tests']=int(results.p_value.notna().sum())
    write_json(root/'candidate_validation.json',result)
    return result
