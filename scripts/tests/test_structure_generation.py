"""Explicit structure selection is immutable, reviewed and tied to the cell lineage."""
import json
import os
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'uniform'))
import structure_generation as S
import collagen_structures as C


def write(p, value):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(value))
    return {'path': str(p), 'sha256': S._sha(p)}


@pytest.fixture
def bound(tmp_path, monkeypatch):
    run = tmp_path / 'run'; root = tmp_path / 'bundle'; tags = ['one', 'two']
    S._CACHE.clear()
    cohort = tmp_path / 'cohort.csv'; cohort.write_text('image_name\none\ntwo\n')
    monkeypatch.setattr(S, 'COHORT', cohort)
    monkeypatch.setattr(S, '_cohort', lambda: tags)
    calls = []
    def lineage(inv, current_run, tag, pinned):
        calls.append(tag)
        actual = json.loads((run / 'lineage' / f'{tag}.json').read_text())
        if pinned != actual: raise ValueError('cell lineage changed')
        inv.add(run / 'lineage' / f'{tag}.json'); inv.records(actual)
        return actual
    monkeypatch.setattr(S, '_lineage', lineage)
    monkeypatch.setattr(C, 'validate_record', lambda *a, **k: {})
    method = write(root / 'method.py', 'common method')
    edge = write(root / 'old_edge.groovy', 'historical edge predicate')
    restoration = {'schema': 'lymphatic_bounded_restoration/1', 'source_run': str(run),
                   'baseline_run': str(tmp_path / 'baseline'), 'method': {
                       'name': 'common rule', 'code_path': method['path'], 'code_sha256': method['sha256'],
                       'original_edge_script_path': edge['path'], 'original_edge_script_sha256': edge['sha256']},
                   'images': {}}
    wrapper = {'schema': 'uniform.structure_generation/1', 'cell_run': str(run),
               'cohort_sha256': S._sha(cohort), 'images': {}}
    for tag in tags:
        tissue = write(run / 'inputs/tissue' / f'{tag}.geojson', {'tissue': tag})
        lineage = {'files': {'tissue': tissue}, 'filter_manifest_sha256': 'filter-' + tag}
        write(run / 'lineage' / f'{tag}.json', lineage)
        old = write(tmp_path / 'baseline/structures/Lymphatic' / f'Lymphatic__{tag}.geojson', {'old': tag})
        old_tissue = write(tmp_path / 'baseline/inputs/tissue' / f'{tag}.geojson', {'old tissue': tag})
        item = {'cell_inputs': lineage, 'classes': {}, 'reviews': {}}
        for cls in S.CLASSES:
            current = write(run / 'structures' / cls / f'{cls}__{tag}.geojson', {'class': cls, 'image': tag})
            item['classes'][cls] = (write(root / 'analysis' / f'{tag}.geojson', {'new lymphatic': tag})
                                    if cls == 'Lymphatic' else current)
        selected = item['classes']['Lymphatic']
        classifier = {'path': str(run / 'structures/Lymphatic' / f'Lymphatic__{tag}.geojson'),
                      'sha256': S._sha(run / 'structures/Lymphatic' / f'Lymphatic__{tag}.geojson')}
        q = {'image': tag, 'method': 'common rule', 'output_path': selected['path'], 'output_sha256': selected['sha256'],
             'inputs': {'tissue': tissue, 'old_tissue': old_tissue, 'baseline': old, 'current_classifier': classifier},
             'checks': {'passed': True, 'tolerance_um2': 1e-7,
              'original_layer_bytes_equal': True, 'source_files_unchanged': True,
              **{name: 0. for name in ('baseline_missing_inside_final_tissue_um2', 'outside_final_tissue_um2',
                  'added_in_old_tissue_outside_selected_edge_components_um2', 'unauthorized_addition_um2',
                  'reconstructed_edge_geometry_missing_um2', 'delta_inference_missing_um2')}}}
        item['quantitative'] = write(root / 'records' / f'{tag}.json', q)
        restored = {'classes': {'Lymphatic': selected}, 'record': item['quantitative'], 'tissue': tissue,
                    'old_tissue': old_tissue, 'current_classifier': classifier,
                    'original_layer': {**old, 'source_path': old['path'], 'source_sha256': old['sha256']}}
        for cls in ('Nerve', 'Capillary', 'Lymphatic'):
            folder = root / 'reviews' / cls
            views = [f'{tag}_whole.png', *[f'{tag}_tile{i:02}.png' for i in range(1, 7)]]
            hashes = {v: write(folder / v, {'native view': v})['sha256'] for v in views}
            tiles = write(folder / f'{tag}_tiles.json', [{'tile': i, 'x_px': i * 800, 'y_px': 0,
                            'width_px': 800, 'height_px': 800} for i in range(1, 7)])
            checks = write(folder.parent / 'checks' / f'{tag}.csv', {'image': tag})
            review = {'image': tag, 'class': cls, 'status': 'acceptable_with_notes',
                      'structure_sha256': item['classes'][cls]['sha256'], 'tissue_sha256': tissue['sha256'],
                      'views_reviewed': views, 'file_sha256': hashes,
                      'tiles_metadata_sha256': tiles['sha256'], 'checks_sha256': checks['sha256']}
            if cls == 'Lymphatic':
                review.update(quantitative_sha256=item['quantitative']['sha256'], tile_manifest_path=tiles['path'],
                              tile_manifest_sha256=tiles['sha256'])
                meta = {k: review[k] for k in ('image', 'structure_sha256', 'tissue_sha256',
                        'quantitative_sha256', 'tile_manifest_sha256', 'file_sha256')}
                meta['source_inputs'] = q['inputs']
                restored['renders'] = write(folder / f'{tag}.metadata.json', meta)
                review['metadata_sha256'] = restored['renders']['sha256']
            item['reviews'][cls] = write(folder / f'{tag}.review.json', review)
        item['collagen_sidecar'] = write(C.checkpoint_path(run, tag), {'output_sha256': item['classes']['Collagen IV']['sha256']})
        wrapper['images'][tag] = item; restoration['images'][tag] = restored
    def publish():
        wrapper['restoration'] = write(root / 'restoration.json', restoration)
        pointer = write(root / 'manifest.json', wrapper)
        write(run / S.BINDING, {'schema': 'uniform.structure_generation_binding/1',
                               'cell_run': str(run), 'manifest': pointer})
    publish()
    return run, wrapper, restoration, publish, calls


def test_legacy_and_broken_binding_never_falls_back(tmp_path):
    assert S.validate(tmp_path) is None
    assert S.provenance(tmp_path) is None
    assert S.path(tmp_path, 'one', 'Lymphatic') == tmp_path / 'structures/Lymphatic/Lymphatic__one.geojson'
    (tmp_path / S.BINDING).symlink_to(tmp_path / 'missing.json')
    with pytest.raises(FileNotFoundError): S.path(tmp_path, 'one', 'Lymphatic')


def test_directory_entry_guard_detects_addition_with_unchanged_stat(tmp_path, monkeypatch):
    folder = tmp_path / 'raw_evidence'; folder.mkdir()
    inv = S._Inventory(); inv.directory(folder)
    initial = folder.stat(); original_stat = Path.stat
    monkeypatch.setattr(Path, 'stat', lambda p, *a, **k: initial if p == folder else original_stat(p, *a, **k))
    (folder / 'new_evidence.json').write_text('{}')
    assert not inv.unchanged()


def test_selection_cohort_proof_is_cached_with_file_guards(bound):
    run, wrapper, _, _, calls = bound
    result = S.validate(run)
    assert set(result) == {'binding', 'manifest', 'provenance', 'paths', 'files', 'identities'}
    assert calls == ['one', 'two']
    for _ in range(3):
        assert S.path(run, 'one', 'Lymphatic') == Path(wrapper['images']['one']['classes']['Lymphatic']['path'])
        assert S.provenance(run) == result['provenance']
    assert calls == ['one', 'two']
    assert not any('Collagen IV' in str(x) and 'review' in str(x) for x in result['files'])
    with pytest.raises(ValueError, match='outside'): S.validate(run, ['three'])


def test_cache_rejects_changed_view_even_if_mtime_restored(bound):
    run, wrapper, _, _, _ = bound; S.validate(run)
    p = Path(wrapper['images']['one']['reviews']['Lymphatic']['path']).with_name('one_tile01.png')
    before = p.stat(); p.write_text('changed actual view'); os.utime(p, ns=(before.st_atime_ns, before.st_mtime_ns))
    with pytest.raises(ValueError, match='hash mismatch'): S.validate(run)


@pytest.mark.parametrize('kind', ['Nerve', 'Capillary', 'Collagen IV'])
def test_only_lymphatic_can_change(bound, kind):
    run, wrapper, _, publish, _ = bound
    wrapper['images']['one']['classes'][kind] = write(run.parent / f'new-{kind}.geojson', {'replacement': kind})
    publish()
    with pytest.raises(ValueError, match='hash mismatch'): S.validate(run)


@pytest.mark.parametrize('field,value', [('status', 'needs_review'), ('image', 'wrong'),
    ('structure_sha256', 'stale'), ('views_reviewed', ['one_whole.png', 'one_tile01.png'])])
def test_requires_complete_actual_current_review(bound, field, value):
    run, wrapper, _, publish, _ = bound
    pointer = wrapper['images']['one']['reviews']['Lymphatic']; p = Path(pointer['path'])
    review = json.loads(p.read_text()); review[field] = value
    wrapper['images']['one']['reviews']['Lymphatic'] = write(p, review); publish()
    with pytest.raises(ValueError): S.validate(run)


@pytest.mark.parametrize('field,value', [('unauthorized_addition_um2', 0.01),
    ('tolerance_um2', 1.), ('outside_final_tissue_um2', -1.), ('outside_final_tissue_um2', float('nan')),
    ('original_layer_bytes_equal', False)])
def test_quantitative_invariants_are_not_self_waived(bound, field, value):
    run, wrapper, restoration, publish, _ = bound
    pointer = wrapper['images']['one']['quantitative']; p = Path(pointer['path'])
    q = json.loads(p.read_text()); q['checks'][field] = value
    pointer = write(p, q); wrapper['images']['one']['quantitative'] = pointer
    restoration['images']['one']['record'] = pointer; publish()
    with pytest.raises(ValueError, match='tolerance|quantitative|geometry invariant'): S.validate(run)


def test_cell_lineage_and_method_changes_invalidate_cached_proof(bound):
    run, wrapper, restoration, publish, _ = bound; S.validate(run)
    wrapper['images']['one']['cell_inputs']['filter_manifest_sha256'] = 'other-filter'
    publish()
    with pytest.raises(ValueError, match='cell lineage'): S.validate(run)


def test_method_code_and_source_files_are_guarded(bound):
    run, _, restoration, _, _ = bound; S.validate(run)
    Path(restoration['method']['code_path']).write_text('edited method')
    with pytest.raises(ValueError, match='hash mismatch'): S.provenance(run)


def test_superseded_actual_review_archive_is_pinned(bound):
    run, wrapper, _, publish, _ = bound
    pointer = wrapper['images']['one']['reviews']['Lymphatic']; p = Path(pointer['path'])
    review = json.loads(p.read_text())
    prior = p.with_name('one.initial_review.json')
    review['supersedes'] = write(prior, {'image': 'one', 'class': 'Lymphatic', 'status': 'acceptable'})
    wrapper['images']['one']['reviews']['Lymphatic'] = write(p, review); publish()
    assert str(prior) in S.validate(run)['files']
    prior.write_text('changed archival evidence')
    with pytest.raises(ValueError, match='hash mismatch'): S.validate(run)


def test_exact_cohort_required(bound):
    run, wrapper, _, publish, _ = bound
    wrapper['images'].pop('two'); publish()
    with pytest.raises(ValueError, match='exactly the cohort'): S.validate(run)


@pytest.mark.parametrize('defect', [None, 'missing_actual_extra', 'wrong_supplement_geometry', 'duplicate_windows'])
def test_supplemental_native_fields_have_exact_ancestry_and_distinct_coverage(bound, defect):
    run, wrapper, restoration, publish, _ = bound
    item = wrapper['images']['one']; restored = restoration['images']['one']
    review_path = Path(item['reviews']['Lymphatic']['path']); review = json.loads(review_path.read_text())
    original_path = Path(review['tile_manifest_path']); positions = json.loads(original_path.read_text())
    # Preserve the first six records, including their old duplicate; add a new distinct field.
    positions[-1]['x_px'] = positions[-2]['x_px']
    original = write(original_path, positions)
    meta_path = Path(restored['renders']['path']); meta = json.loads(meta_path.read_text())
    meta['tile_manifest_sha256'] = original['sha256']
    original_meta = write(meta_path, meta); restored['renders'] = original_meta
    review['metadata_sha256'] = original_meta['sha256']
    extra = 'one_tile07.png'; extra_hash = write(review_path.parent / extra, {'new native field': 7})['sha256']
    extra_position = {'tile': 7, 'x_px': 4000 if defect == 'duplicate_windows' else 8000,
                      'y_px': 0, 'width_px': 800, 'height_px': 800}
    supplement = {'schema': 'bounded_lymphatic_supplemental_native/1', 'image': 'one',
                  'source_inputs': meta['source_inputs'], 'original_manifest': original, 'original_metadata': original_meta,
                  **{k: review[k] for k in ('structure_sha256', 'tissue_sha256', 'quantitative_sha256')},
                  'file_sha256': {extra: extra_hash}, 'tiles': [extra_position]}
    if defect == 'wrong_supplement_geometry': supplement['structure_sha256'] = 'unrelated'
    supplemental_pointer = write(review_path.parent / 'one.supplemental_metadata.json', supplement)
    combined = {'schema': 'bounded_lymphatic_combined_native_windows/1', 'image': 'one',
                'original_manifest': original, 'original_metadata': original_meta,
                'supplemental_metadata': supplemental_pointer,
                **{k: review[k] for k in ('structure_sha256', 'tissue_sha256', 'quantitative_sha256')},
                'tiles': [*positions, extra_position], 'distinct_windows': 6, 'separated_windows': 6}
    combined_pointer = write(review_path.parent / 'one_tiles_combined.json', combined)
    review['tile_manifest_path'] = combined_pointer['path']; review['tile_manifest_sha256'] = combined_pointer['sha256']
    review['supplemental_metadata_sha256'] = supplemental_pointer['sha256']
    if defect != 'missing_actual_extra':
        review['views_reviewed'].append(extra); review['file_sha256'][extra] = extra_hash
    item['reviews']['Lymphatic'] = write(review_path, review); publish()
    if defect is None:
        assert S.validate(run)['manifest']['images']['one']['reviews']['Lymphatic'] == item['reviews']['Lymphatic']
    else:
        with pytest.raises(ValueError): S.validate(run)


def test_endpoint_reads_selected_mask_under_stable_logical_key(tmp_path, monkeypatch):
    import numpy as np
    import pandas as pd
    from shapely.geometry import box, mapping
    import structure_metrics as M
    def feature(g): return {'geometry': mapping(g), 'properties': {}}
    write(tmp_path / 'inputs/tissue/image.geojson', {'features': [feature(box(0, 0, 10, 10))]})
    write(tmp_path / 'islets/image.geojson', {'features': []})
    (tmp_path / 'cells').mkdir()
    pd.DataFrame({'Object ID': ['a'], 'Centroid X µm': [1.], 'Centroid Y µm': [1.]}).to_parquet(tmp_path / 'cells/image.parquet')
    selected = tmp_path / 'bundle/selected.geojson'
    write(selected, {'features': [feature(box(0, 0, 2, 10))]})
    for cls in M.CLASSES:
        write(tmp_path / 'structures' / cls / f'{cls}__image.geojson', {'features': [feature(box(0, 0, 10, 10))]})
    legacy = S.path
    monkeypatch.setattr(M.SG, 'path', lambda run, tag, cls: selected if cls == 'Lymphatic' else legacy(run, tag, cls))
    monkeypatch.setattr(M.SG, 'provenance', lambda run: {'manifest_sha256': 'chosen'})
    result = M.one(tmp_path, {'image_name': 'image', 'case_id': '1'})
    lymph = next(r for r in result['records'] if r['class'] == 'Lymphatic')
    assert np.isclose(lymph['area_density'], .2)
    assert result['inputs']['structures/Lymphatic/Lymphatic__image.geojson'] == S._sha(selected)
    assert result['structure_generation'] == {'manifest_sha256': 'chosen'}


def test_cell_context_uses_selected_geometry_and_records_generation(tmp_path, monkeypatch):
    import pandas as pd
    from shapely.geometry import box, mapping
    import cell_context as CC
    cells = pd.DataFrame({'Object ID': ['a', 'b'], 'Centroid X µm': [2. * CC.CT.PIXEL_UM, 8. * CC.CT.PIXEL_UM],
                          'Centroid Y µm': [5. * CC.CT.PIXEL_UM] * 2,
                          'Nucleus: Max diameter µm': [4., 4.], 'Nucleus: Min diameter µm': [3., 3.]})
    for folder in ('cells', 'typing/redsea', 'context'):
        (tmp_path / folder).mkdir(parents=True)
    cells.to_parquet(tmp_path / 'cells/image.parquet')
    pd.DataFrame({'object_id': ['a', 'b'], 'qc_analysis_eligible': [True, True]}).to_parquet(tmp_path / 'typing/redsea/image.parquet')
    pd.DataFrame({'Object ID': ['a', 'b'], 'Parent': ['Image', 'Image']}).to_parquet(tmp_path / 'context/image.parquet')
    write(tmp_path / 'islets/image.geojson', {'features': []})
    write(tmp_path / 'cell_geojson/image.geojson', {'features': [
        {'id': 'a', 'geometry': mapping(box(1, 4, 3, 6))}, {'id': 'b', 'geometry': mapping(box(7, 4, 9, 6))}]})
    selected = tmp_path / 'bundle/selected.geojson'
    write(selected, {'features': [{'geometry': mapping(box(0, 0, 4, 10))}]})
    for cls in ('Nerve', 'Capillary', 'Lymphatic'):
        write(tmp_path / 'structures' / cls / f'{cls}__image.geojson',
              {'features': [{'geometry': mapping(box(0, 0, 10, 10))}]} if cls == 'Lymphatic' else {'features': []})
    generation = {'manifest_sha256': 'selected'}
    legacy = S.path
    monkeypatch.setattr(CC.SG, 'path', lambda run, tag, cls: selected if cls == 'Lymphatic' else legacy(run, tag, cls))
    monkeypatch.setattr(CC.SG, 'provenance', lambda run: generation)
    record = CC.process(tmp_path, 'image')
    result = pd.read_parquet(tmp_path / 'typing/context/image.parquet')
    assert result['in_Lymphatic'].tolist() == [True, False]
    assert result['crossed_frac_lymphatic'].tolist() == [1., 0.]
    assert record['structure_generation'] == generation
    assert json.loads((tmp_path / 'typing/context/image.json').read_text())['structure_generation'] == generation
    assert pd.read_parquet(tmp_path / 'context/image.parquet')['Parent'].tolist() == ['Image', 'Image']


@pytest.fixture
def previous_bound(bound, monkeypatch):
    run, wrapper, selection, publish, calls = bound
    baseline = Path(selection['baseline_run'])
    monkeypatch.setattr(S, 'PREVIOUS_BASELINE', baseline)
    instruction = write(run.parent / 'user_instruction.json', {
        'decision': 'reuse_exact_previous_lymphatic_objects', 'cell_run': str(run),
        'baseline_run': str(baseline), 'user_instruction': 'Pause regeneration and use previous lymphatic objects.',
        'supersedes_reconstruction_gate': True})
    selection['schema'] = 'lymphatic_exact_previous_reuse/1'
    selection['user_instruction'] = instruction
    selection.pop('method')
    for tag, item in wrapper['images'].items():
        old = selection['images'][tag]['original_layer']
        selected = {'path': old['source_path'], 'sha256': old['source_sha256']}
        selection['images'][tag] = {'classes': {'Lymphatic': selected},
                                    'tissue': selection['images'][tag]['tissue']}
        item['classes']['Lymphatic'] = selected
        item.pop('quantitative'); item['reviews'].pop('Lymphatic')
    publish()
    return run, wrapper, selection, publish, calls


def test_previous_objects_are_explicit_and_do_not_claim_new_review(previous_bound):
    run, wrapper, selection, _, calls = previous_bound
    result = S.validate(run)
    assert result['provenance']['lymphatic_selection'] == 'exact_previous_objects'
    assert result['provenance']['user_instruction_sha256'] == selection['user_instruction']['sha256']
    assert S.path(run, 'one', 'Lymphatic') == Path(selection['baseline_run']) / 'structures/Lymphatic/Lymphatic__one.geojson'
    assert 'Lymphatic' not in wrapper['images']['one']['reviews']
    assert str(Path(selection['user_instruction']['path'])) in result['files']
    assert calls == ['one', 'two']
    assert S.provenance(run) == result['provenance']
    assert calls == ['one', 'two']


@pytest.mark.parametrize('defect', ['wrong_baseline', 'replacement', 'no_override', 'wrong_tissue', 'changed_nerve'])
def test_previous_objects_branch_is_narrow_and_fail_closed(previous_bound, defect):
    run, wrapper, selection, publish, _ = previous_bound
    if defect == 'wrong_baseline': selection['baseline_run'] = str(run)
    elif defect == 'replacement':
        wrong = write(run.parent / 'new_lymphatic.geojson', {'not': 'original'})
        wrapper['images']['one']['classes']['Lymphatic'] = wrong
        selection['images']['one']['classes']['Lymphatic'] = wrong
    elif defect == 'no_override':
        p = Path(selection['user_instruction']['path']); value = json.loads(p.read_text())
        value['supersedes_reconstruction_gate'] = False
        selection['user_instruction'] = write(p, value)
    elif defect == 'wrong_tissue':
        selection['images']['one']['tissue'] = selection['images']['two']['tissue']
    elif defect == 'changed_nerve':
        wrapper['images']['one']['classes']['Nerve'] = write(run.parent / 'changed_nerve.geojson', {'changed': True})
    publish()
    with pytest.raises(ValueError): S.validate(run)


def test_previous_original_and_instruction_bytes_remain_guarded(previous_bound):
    run, _, selection, _, _ = previous_bound
    S.validate(run)
    p = Path(selection['images']['one']['classes']['Lymphatic']['path'])
    p.write_text('changed baseline bytes')
    with pytest.raises(ValueError, match='hash mismatch'): S.validate(run)


DG_METHOD = 'despeckled_detection_mean_despeckled_pdpn_at_or_above_image_level/1'


def polygons_geojson(order):
    squares = [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]], [[2, 0], [3, 0], [3, 1], [2, 1], [2, 0]]]
    return {'type': 'FeatureCollection', 'features': [
        {'type': 'Feature', 'id': f'uuid-{order}-{i}', 'geometry': {'type': 'Polygon', 'coordinates': [squares[i]]}}
        for i in (order if isinstance(order, list) else [0, 1])]}


@pytest.fixture
def despeckled_bound(bound):
    """bound, turned into the owner's DG selection (2026-11-01): despeckled detection, then the PDPN level gate."""
    run, wrapper, selection, publish, calls = bound
    base = run.parent / 'typing_20261030'
    classifier = write(run.parent / 'qpi/classifiers/pixel_classifiers/Lymphatics.json', {'classifier': 'Lymphatics'})
    detection = {'driver': write(base / 'lymphatic_despeckle_run.py', 'driver'),
                 'groovy': write(base / 'lymphatics_despeckled.groovy', 'groovy'),
                 'despeckle_driver': write(base / 'despeckle_pdpn.py', 'despeckle'), 'classifier': classifier,
                 'despeckle_code_sha256': 'a' * 64}
    reviewed = write(run / 'structures/Lymphatic/Lymphatic__one.geojson', polygons_geojson([0, 1]))
    original = write(base / 'runs/original/one/Lymphatic.geojson', polygons_geojson([1, 0]))
    detection['validation'] = {'one': {'result': write(base / 'runs/original/one/result.json',
                                                       {'mode': 'original', 'image': 'one', 'output': original}),
                                       'reviewed_rerun_objects': reviewed}}
    entries = []
    instruction = write(base / 'gate/owner_instruction.json', {'decision': 'despeckled_detection_level_gate',
                                                               'cell_run': str(run)})
    images = {}
    for i, (tag, item) in enumerate(wrapper['images'].items(), start=1):
        source = run.parent / 'images' / f'{tag}.tif'
        source.parent.mkdir(parents=True, exist_ok=True); source.write_bytes(b'tif ' + tag.encode())
        entries.append({'entry_id': i, 'image_name': tag, 'image_path': str(source)})
        server = write(run / 'qupath/data' / str(i) / 'server.json', {'uri': str(source)})
        raw = write(base / 'lymphatic_despeckle/pdpn' / f'{tag}.u16', 'despeckled ' + tag)
        st = source.stat()
        pdpn = write(base / 'lymphatic_despeckle/pdpn' / f'{tag}.json', {
            'image': tag, 'channel': 'PDPN', 'path': raw['path'], 'sha256': raw['sha256'],
            'despeckle_code_sha256': 'a' * 64, 'source': {'path': str(source), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns}})
        tissue = item['cell_inputs']['files']['tissue']
        objects = write(base / 'lymphatic_despeckle/runs/despeckled' / tag / 'Lymphatic.geojson', polygons_geojson([0, 1]))
        config = write(base / 'lymphatic_despeckle/runs/despeckled' / f'{tag}.config.json', {'image': tag})
        result = write(base / 'lymphatic_despeckle/runs/despeckled' / tag / 'result.json', {
            'schema': 'uniform.lymphatics_despeckled_result/1', 'image': tag, 'mode': 'despeckled',
            'inputs': {server['path']: server['sha256'], classifier['path']: classifier['sha256'],
                       tissue['path']: tissue['sha256'], raw['path']: raw['sha256']},
            'output': objects, 'config': config})
        det = {'result': result, 'config': config, 'pdpn_record': pdpn, 'pdpn_raw': raw}
        previous = write(run.parent / 'previous' / f'Lymphatic__{tag}.geojson', {'previous': tag})
        selected = item['classes']['Lymphatic']
        item['quantitative'] = write(base / 'gate/records' / f'{tag}.json', {
            'image': tag, 'method': DG_METHOD, 'output_path': selected['path'], 'output_sha256': selected['sha256'],
            'inputs': {'tissue': tissue, 'source_objects': objects, 'previous_objects': previous, 'detection': det},
            'checks': {'passed': True}})
        review_path = Path(item['reviews']['Lymphatic']['path'])
        review = json.loads(review_path.read_text()); review['quantitative_sha256'] = item['quantitative']['sha256']
        item['reviews']['Lymphatic'] = write(review_path, review)
        images[tag] = {'tissue': tissue, 'source_objects': objects, 'previous_objects': previous, 'detection': det,
                       'pdpn_level_raw': 1.0, 'classes': {'Lymphatic': selected}, 'record': item['quantitative'],
                       'review': item['reviews']['Lymphatic']}
    write(run / 'project_manifest.json', {'entries': entries})
    baseline = selection['baseline_run']
    selection.clear()
    selection.update({'schema': 'lymphatic_despeckled_level_gate/1', 'source_run': str(run), 'baseline_run': baseline,
                      'method': {'name': DG_METHOD, **{k: v for k, v in zip(('code_path', 'code_sha256'),
                                 write(base / 'lymphatic_despeckle_gate.py', 'builder').values())},
                                 'detection': detection},
                      'owner_instruction': instruction, 'levels': write(base / 'gate/levels_v19.csv', 'levels'),
                      'images': images})
    publish()
    return run, wrapper, selection, publish, calls


def test_despeckled_gate_validates_and_names_the_selection(despeckled_bound):
    run, wrapper, selection, _, calls = despeckled_bound
    result = S.validate(run)
    assert result['provenance']['lymphatic_selection'] == 'despeckled_level_gate'
    assert result['provenance']['owner_instruction_sha256'] == selection['owner_instruction']['sha256']
    one = selection['images']['one']['detection']
    for pointer in (one['result'], one['pdpn_record'], one['pdpn_raw'], selection['method']['detection']['groovy'],
                    selection['method']['detection']['classifier']):
        assert pointer['path'] in result['files']
    assert any(p.endswith('images/one.tif') for p in result['identities'])
    assert S.path(run, 'one', 'Lymphatic') == Path(wrapper['images']['one']['classes']['Lymphatic']['path'])
    assert calls == ['one', 'two']


def rewrite(pointer, change):
    p = Path(pointer['path']); value = json.loads(p.read_text()); change(value)
    return write(p, value)


DG_DEFECTS = {'decision': 'lacks the owner decision', 'mode': 'does not produce the gate source',
              'output': 'does not produce the gate source', 'tissue': 'read another Tissue',
              'classifier': 'read another Tissue', 'pdpn': 'not this image', 'source_image': 'not this image',
              'validation': 'does not reproduce the reviewed re-run objects', 'record_inputs': 'different inputs',
              'server': 'read another Tissue'}


@pytest.mark.parametrize('defect', sorted(DG_DEFECTS))
def test_despeckled_gate_is_fail_closed(despeckled_bound, defect):
    run, wrapper, selection, publish, _ = despeckled_bound
    one = selection['images']['one']; det = one['detection']
    if defect == 'decision':
        selection['owner_instruction'] = rewrite(selection['owner_instruction'], lambda v: v.update(decision='pdpn_level_object_gate'))
    elif defect == 'mode':
        det['result'] = rewrite(det['result'], lambda v: v.update(mode='original'))
    elif defect == 'output':
        det['result'] = rewrite(det['result'], lambda v: v.update(output=selection['images']['two']['source_objects']))
    elif defect == 'tissue':
        det['result'] = rewrite(det['result'], lambda v: v['inputs'].update({one['tissue']['path']: 'b' * 64}))
    elif defect == 'classifier':
        path = selection['method']['detection']['classifier']['path']
        det['result'] = rewrite(det['result'], lambda v: v['inputs'].update({path: 'c' * 64}))
    elif defect == 'pdpn':
        det['pdpn_record'] = rewrite(det['pdpn_record'], lambda v: v.update(sha256='d' * 64))
    elif defect == 'source_image':
        det['pdpn_record'] = rewrite(det['pdpn_record'], lambda v: v['source'].update(path=str(run / 'other.tif')))
    elif defect == 'validation':
        v = selection['method']['detection']['validation']['one']
        out = json.loads(Path(v['result']['path']).read_text())['output']
        changed = polygons_geojson([0]); write(Path(out['path']), changed)
        v['result'] = rewrite(v['result'], lambda r: r.update(output=write(Path(out['path']), changed)))
    elif defect == 'record_inputs':
        wrapper['images']['one']['quantitative'] = rewrite(
            wrapper['images']['one']['quantitative'], lambda q: q['inputs'].update(detection={**det, 'config': det['result']}))
        one['record'] = wrapper['images']['one']['quantitative']
    elif defect == 'server':
        det['result'] = rewrite(det['result'], lambda v: v['inputs'].update({str(run.parent / 'x/qupath/data/9/server.json'): 'e' * 64}))
    if defect in ('mode', 'output', 'tissue', 'classifier', 'pdpn', 'source_image', 'server', 'record_inputs'):
        q = wrapper['images']['one']['quantitative']
        if defect != 'record_inputs':
            wrapper['images']['one']['quantitative'] = rewrite(q, lambda value: value['inputs'].update(detection=det))
            one['record'] = wrapper['images']['one']['quantitative']
    publish()
    with pytest.raises(ValueError, match=DG_DEFECTS[defect]): S.validate(run)


def test_despeckled_source_image_identity_is_guarded(despeckled_bound):
    run, _, _, _, _ = despeckled_bound
    S.validate(run)
    source = run.parent / 'images/one.tif'
    source.write_bytes(b'changed acquisition bytes')
    S._CACHE.clear()
    with pytest.raises(ValueError, match='identity changed'): S.validate(run)
