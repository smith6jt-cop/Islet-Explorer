"""Preservation and restart tests for the additive Tissue augmentation path."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd
import pytest

UNIFORM = Path(__file__).resolve().parents[1] / 'uniform'
sys.path.insert(0, str(UNIFORM))
import augment_tissue as AU


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def test_augmentation_has_no_destructive_or_islet_steps():
    steps = [step for stage in AU.STAGES for step in AU.stage_steps(stage)]
    assert [s.id for s in steps] == ['augment_detect', 'augment_measure', 'augment_3', 'augment_8',
                                    'augment_C_cells', 'augment_C_geojson', 'augment_J_cells',
                                    'augment_J_context', 'augment_J_annotations', 'augment_J_islets', 'augment_verify']
    assert not any(s.script.endswith(('import_tissue.groovy', '0_cell_detection.groovy', '1_islet_objects.groovy',
                                     '1b_islet_outlines_to_cells.groovy')) for s in steps)


def test_source_cannot_be_a_linked_output(tmp_path):
    source = tmp_path / 'source'
    source.write_bytes(b'qpdata')
    linked = tmp_path / 'hardlink'
    linked.hardlink_to(source)
    with pytest.raises(ValueError, match='independent'):
        AU.check_independent(source, linked)
    independent = tmp_path / 'copy'
    shutil.copy2(source, independent)
    AU.check_independent(source, independent)


@pytest.fixture
def exports(tmp_path):
    source, out = tmp_path / 'source', tmp_path / 'out'
    def cell(oid, x):
        return {'type': 'Feature', 'id': oid,
                'geometry': {'type': 'Polygon', 'coordinates': [[[x, 0], [x + 1, 0], [x + 1, 1], [x, 0]]]},
                'nucleusGeometry': {'type': 'Point', 'coordinates': [x + .2, .2]}}
    obj = {'type': 'Feature', 'id': 'islet-id', 'geometry': cell('unused', 0)['geometry'],
           'properties': {'name': 'Islet_1', 'classification': {'name': 'Islet'}, 'measurements': {'Islet: UID': 1}}}
    for run, ids in [(source, ['a']), (out, ['a', 'b'])]:
        for folder in ('cells', 'final_cells', 'context'):
            (run / folder).mkdir(parents=True)
            pd.DataFrame({'Object ID': ids, 'Cell: DAPI: Mean': [12.] * len(ids),
                          'Cell: Area µm^2': [1.] * len(ids), 'Centroid X µm': list(range(len(ids)))}
                         ).to_parquet(run / folder / 'image.parquet', index=False)
        write_json(run / 'cell_geojson/image.geojson', {'features': [cell(oid, i) for i, oid in enumerate(ids)]})
        write_json(run / 'islets/image.geojson', {'features': [obj]})
    return source, out


def test_export_verification_preserves_values_geometry_and_ids(exports):
    source, out = exports
    result = AU.verify_exports(source, out, 'image', {'b'})
    assert result['original_cells'] == 1 and result['new_cells'] == 1
    table = pd.read_parquet(out / 'cells/image.parquet')
    table.loc[0, 'Cell: DAPI: Mean'] = 11.
    table.to_parquet(out / 'cells/image.parquet')
    with pytest.raises(ValueError, match='original measurement changed'):
        AU.verify_exports(source, out, 'image', {'b'})


def test_export_verification_rejects_changed_endocrine_geometry(exports):
    source, out = exports
    obj = AU.read_json(out / 'islets/image.geojson')
    obj['features'][0]['geometry']['coordinates'][0][1][0] = 2
    write_json(out / 'islets/image.geojson', obj)
    with pytest.raises(ValueError, match='endocrine'):
        AU.verify_exports(source, out, 'image', {'b'})


def test_export_verification_rejects_missing_new_id(exports):
    source, out = exports
    with pytest.raises(ValueError, match='accepted additions'):
        AU.verify_exports(source, out, 'image', {'missing'})


def test_checkpoint_resume_and_saved_stage_recovery(tmp_path, monkeypatch):
    source, out, outlines = tmp_path / 'source', tmp_path / 'out', tmp_path / 'outlines'
    for run in (source, out):
        write_json(run / 'qupath/project.qpproj', {'images': [{'imageName': 'image', 'entryID': 1}]})
        (run / 'qupath/data/1').mkdir(parents=True)
        (run / 'qupath/data/1/data.qpdata').write_bytes(b'original')
    write_json(outlines / 'image.geojson', {'features': []})
    runner = AU.RU.Runner(out, [{'image_name': 'image'}], UNIFORM / 'cohort.csv')
    calls = []
    def fake_run(tag, step):
        calls.append(step.id)
        cfg = AU.read_json(out / 'augmentation/image/config.json')
        write_json(Path(cfg['inventory']), {'cells': {'a': 'hash'}, 'annotations': {}})
        write_json(Path(cfg['new_ids']), ['b'])
        target = AU.qpdata(out, tag)
        target.write_bytes(b'augmented')
        write_json(runner.record(tag, step.id), {'status': 'pass',
                   'script': {'source_sha256': AU.RU.sha256(UNIFORM / step.script)},
                   'qpdata': {'after': AU.RU.stat_of(target)},
                   'outputs': [{'path': cfg['inventory'], 'exists': True}, {'path': cfg['new_ids'], 'exists': True}]})
        return True
    monkeypatch.setattr(runner, 'run_qupath', fake_run)
    initial = AU.process(runner, source, outlines, 'image', ['detect'])
    assert initial['inputs']['detection_context_um'] == 30.0
    assert AU.read_json(out / 'augmentation/image/config.json')['detection_context_um'] == 30.0
    AU.process(runner, source, outlines, 'image', ['detect'])
    assert calls == ['augment_detect']
    with pytest.raises(ValueError, match='augmentation inputs changed'):
        AU.process(runner, source, outlines, 'image', ['detect'], detection_context_um=40.)
    # Simulate interruption after QuPath saved its successful provenance but before Python checkpointed it.
    state_path = out / 'augmentation/image/state.json'
    interrupted = dict(initial)
    interrupted['last_qpdata_sha256'] = AU.RU.sha256(AU.qpdata(source, 'image'))
    interrupted['in_progress'] = {'id': 'augment_detect',
                                  'fingerprint': initial['steps']['augment_detect']['fingerprint']}
    interrupted['steps'] = {}
    write_json(state_path, interrupted)
    AU.process(runner, source, outlines, 'image', ['detect'])
    assert calls == ['augment_detect']
    AU.qpdata(out, 'image').write_bytes(b'external modification')
    with pytest.raises(ValueError, match='qpdata changed'):
        AU.process(runner, source, outlines, 'image', ['detect'])


def test_failed_stage_can_retry_without_losing_baseline(tmp_path, monkeypatch):
    source, out, outlines = tmp_path / 'source', tmp_path / 'out', tmp_path / 'outlines'
    for run in (source, out):
        write_json(run / 'qupath/project.qpproj', {'images': [{'imageName': 'image', 'entryID': 1}]})
        (run / 'qupath/data/1').mkdir(parents=True)
        (run / 'qupath/data/1/data.qpdata').write_bytes(b'original')
    write_json(outlines / 'image.geojson', {'features': []})
    runner = AU.RU.Runner(out, [{'image_name': 'image'}], UNIFORM / 'cohort.csv')
    monkeypatch.setattr(runner, 'run_qupath', lambda *args: False)
    for _ in range(2):
        with pytest.raises(RuntimeError, match='failed'):
            AU.process(runner, source, outlines, 'image', ['detect'])
        assert AU.qpdata(out, 'image').read_bytes() == b'original'


def test_groovy_compiles_and_real_qupath_object_invariants(tmp_path):
    jars = Path('/home/smith6jt/QuPath/lib/app')
    qupath = Path('/home/smith6jt/QuPath/bin/QuPath')
    if not jars.is_dir() or not qupath.is_file():
        pytest.skip('QuPath Java installation is unavailable')
    fixture = r'''
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.objects.hierarchy.PathObjectHierarchy
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.GeometryTools
import org.locationtech.jts.geom.prep.PreparedGeometryFactory
def shell = new GroovyShell()
args.each { shell.parse(new File(it)) }
def helper = new GroovyShell().evaluate(new File(args[0]))
def detection = shell.getClassLoader().loadClass('AugmentationDetectionGeometry')
def plane = ImagePlane.getDefaultPlane()
def h = new PathObjectHierarchy()
def t = PathObjects.createAnnotationObject(ROIs.createRectangleROI(0, 0, 100, 100, plane), PathClass.fromString('Tissue'))
t.getMeasurementList().put('Area µm^2', 10000.0)
def islet = PathObjects.createAnnotationObject(ROIs.createRectangleROI(2, 2, 20, 20, plane), PathClass.fromString('Islet'))
islet.setName('Islet_1')
def cell = PathObjects.createCellObject(ROIs.createRectangleROI(3, 3, 4, 4, plane),
    ROIs.createRectangleROI(4, 4, 2, 2, plane), PathClass.fromString('Beta cell'), null)
cell.getMeasurementList().put('Cell: DAPI: Mean', 12.0)
h.addObject(t); t.addChildObject(islet); islet.addChildObject(cell)
def before = helper.snapshot(h)
t.setROI(ROIs.createRectangleROI(0, 0, 200, 200, plane))
cell.getMeasurementList().put('Signed distance to annotation with Tissue µm', -3.0)
assert helper.verify(h, before).old_cells == 1
cell.getMeasurementList().put('Cell: DAPI: Mean', 13.0)
boolean rejected = false
try { helper.verify(h, before) } catch (IllegalStateException e) { rejected = true }
assert rejected
cell.getMeasurementList().put('Cell: DAPI: Mean', 12.0)
islet.setROI(ROIs.createRectangleROI(2, 2, 21, 20, plane))
rejected = false
try { helper.verify(h, before) } catch (IllegalStateException e) { rejected = true }
assert rejected
// A full cell straddling the old/new boundary stays whole when detection has context.
def finalROI = ROIs.createRectangleROI(0, 0, 100, 100, plane)
def oldROI = ROIs.createRectangleROI(0, 0, 50, 100, plane)
def target = finalROI.getGeometry(), delta = target.difference(oldROI.getGeometry())
def padded = detection.detectionRegion(delta, 30.0, 1.0, 140, 140)
assert padded.isValid()
assert padded.getEnvelopeInternal().getMinX() >= 0 && padded.getEnvelopeInternal().getMinY() >= 0
assert padded.getEnvelopeInternal().getMaxX() <= 140 && padded.getEnvelopeInternal().getMaxY() <= 140
assert padded.covers(ROIs.createRectangleROI(105, 50, 1, 1, plane).getGeometry()) // not clipped to final Tissue
def addedPrepared = PreparedGeometryFactory.prepare(delta)
def finalPrepared = PreparedGeometryFactory.prepare(target)
def borderCell = PathObjects.createCellObject(ROIs.createRectangleROI(46, 40, 10, 10, plane),
    ROIs.createRectangleROI(48, 42, 6, 6, plane), null, null)
def full = borderCell.getROI().getGeometry()
assert full.intersection(delta).getArea() == 60.0
assert full.intersection(padded).equalsTopo(full) && full.getArea() == 100.0
assert detection.rejectionReason(borderCell, addedPrepared, finalPrepared) == null
// A cell centered in OLD tissue would yield a spurious new fragment if clipped to delta.
def oldCentered = PathObjects.createCellObject(ROIs.createRectangleROI(42, 60, 10, 10, plane),
    ROIs.createRectangleROI(44, 62, 6, 6, plane), null, null)
def clipped = oldCentered.getROI().getGeometry().intersection(delta)
def fragment = PathObjects.createCellObject(GeometryTools.geometryToROI(clipped, plane), null, null, null)
assert detection.rejectionReason(fragment, addedPrepared, finalPrepared) == null
assert detection.rejectionReason(oldCentered, addedPrepared, finalPrepared) == 'outside_added'
assert oldCentered.getROI().getGeometry().intersection(padded).equalsTopo(oldCentered.getROI().getGeometry())
// Positive-area nucleus containment is independent of whole-cell containment.
def badNucleus = PathObjects.createCellObject(ROIs.createRectangleROI(92, 40, 5, 5, plane),
    ROIs.createRectangleROI(98, 40, 5, 5, plane), null, null)
assert detection.rejectionReason(badNucleus, addedPrepared, finalPrepared) == 'nucleus_outside_final'
def badCell = PathObjects.createCellObject(ROIs.createRectangleROI(97, 40, 6, 5, plane), null, null, null)
assert detection.rejectionReason(badCell, addedPrepared, finalPrepared) == 'cell_outside_final'
// Context never alters original cells: add the complete accepted cell beside a fixed old one.
def h2 = new PathObjectHierarchy()
def tissue2 = PathObjects.createAnnotationObject(oldROI, PathClass.fromString('Tissue'))
def fixed = PathObjects.createCellObject(ROIs.createRectangleROI(10, 10, 8, 8, plane),
    ROIs.createRectangleROI(12, 12, 4, 4, plane), PathClass.fromString('Beta cell'), null)
fixed.getMeasurementList().put('Cell: DAPI: Mean', 12.0)
h2.addObject(tissue2); tissue2.addChildObject(fixed)
def saved2 = helper.snapshot(h2)
assert !fixed.getROI().getGeometry().intersects(borderCell.getROI().getGeometry())
tissue2.setROI(finalROI); tissue2.addChildObject(borderCell)
assert helper.verify(h2, saved2).old_cells == 1
assert fixed.getROI().getArea() == 64.0 && fixed.getNucleusROI().getArea() == 16.0
assert fixed.getMeasurementList().get('Cell: DAPI: Mean') == 12.0
println('preservation checks passed')
println('padded detection geometry checks passed')
'''
    script = tmp_path / 'fixture.groovy'
    script.write_text(fixture)
    paths = [UNIFORM / 'qupath' / name for name in ('augmentation_invariants.groovy', 'augment_tissue.groovy',
             'measure_augmented_cells.groovy', 'verify_tissue_augmentation.groovy', 'refresh_augmented_context.groovy')]
    command = [str(qupath), 'script']
    for path in paths:
        command.extend(['-a', str(path)])
    command.append(str(script))
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'preservation checks passed' in result.stdout
    assert 'padded detection geometry checks passed' in result.stdout
