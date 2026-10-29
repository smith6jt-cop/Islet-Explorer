"""Independent samples and prediction-free native QuPath review artifacts."""
from pathlib import Path
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'uniform'))
import blind_review_project as B


@pytest.fixture
def population():
    return pd.DataFrame([{'image': f'image-{i % 3}', 'object_id': f'{label}-{i}', 'label': label}
                         for label, count in [('Common prediction', 100), ('Rare prediction', 7)]
                         for i in range(count)])


def test_sample_is_reproducible_stratified_and_excludes_prior_reviews(population):
    excluded = set(population.iloc[:4][['image', 'object_id']].itertuples(index=False, name=None))
    first, strata = B.sample_by_label(population, excluded, per_label=20, seed=1)
    second, _ = B.sample_by_label(population.sample(frac=1, random_state=5), excluded, per_label=20, seed=1)
    pd.testing.assert_frame_equal(first, second)
    assert first.groupby('label').size().to_dict() == {'Common prediction': 20, 'Rare prediction': 7}
    assert not (set(first[['image', 'object_id']].itertuples(index=False, name=None)) & excluded)
    assert not first[['image', 'object_id']].duplicated().any()
    assert first.sample_id.is_unique and first.blind_uuid.is_unique
    assert first.label.tolist() != sorted(first.label)
    assert strata.set_index('label').loc['Common prediction', 'inverse_sampling_weight'] == 96 / 20
    assert np.isclose(strata.eligible_prevalence.sum(), 1)
    # The inverse weights recover the eligible population size for each stratum.
    assert first.groupby('label').inverse_sampling_weight.sum().to_dict() == {
        'Common prediction': 96., 'Rare prediction': 7.}


def test_duplicate_or_missing_native_ids_are_rejected(population):
    with pytest.raises(ValueError, match='duplicate'):
        B.sample_by_label(pd.concat([population, population.iloc[:1]]), set())
    for column in ('image', 'object_id', 'label'):
        bad = population.copy(); bad.loc[0, column] = ''
        with pytest.raises(ValueError, match='missing'):
            B.sample_by_label(bad, set())
    with pytest.raises(ValueError, match='empty'):
        B.sample_by_label(population.iloc[:0], set())


def test_public_table_has_no_predictions_native_ids_or_weights(population):
    key, _ = B.sample_by_label(population, set())
    key['review_image'] = 'Review image 001'
    key['x_um'] = 1.5; key['y_um'] = 2.5
    public = B.public_table(key)
    assert list(public) == B.PUBLIC_COLUMNS
    assert public[['reviewer_label', 'reviewer_confidence', 'notes']].eq('').all().all()
    text = public.to_csv(index=False)
    assert not any(label in text for label in population.label.unique())
    assert not any(oid in text for oid in key.object_id)
    assert set(public.image) == {'Review image 001'}


def test_project_copies_only_selected_raw_server_definitions(tmp_path):
    raw = {'builderType': 'uri', 'uri': 'file:///raw-image.tif'}
    source = {'version': '0.7.0', 'uri': 'file:///original/project.qpproj',
              'metadata': {'prediction': 'should disappear'}, 'images': [
                  {'entryID': 5, 'imageName': 'original', 'serverBuilder': raw,
                   'metadata': {'prediction': 'Beta cell'}, 'description': 'Beta cell'},
                  {'entryID': 6, 'imageName': 'unused', 'serverBuilder': raw}]}
    result = B.minimal_project(source, tmp_path / 'project.qpproj', {'original': 'Review image 001'})
    assert result['images'] == [{'entryID': 5, 'imageName': 'Review image 001', 'serverBuilder': raw}]
    assert 'prediction' not in json.dumps(result)
    result['images'][0]['serverBuilder']['uri'] = 'changed'
    assert raw['uri'] == 'file:///raw-image.tif'


def test_exclusion_manifest_deduplicates_and_requires_ids(tmp_path):
    path = tmp_path / 'previous.csv'
    path.write_text('image,object_id,label\nimage-1,one,private\nimage-1,one,private\n')
    assert B.exclusions(path) == {('image-1', 'one')}
    path.write_text('image,object_id\nimage-1,\n')
    with pytest.raises(ValueError, match='missing'):
        B.exclusions(path)


@pytest.mark.parametrize('project,key_dir', [('review', 'review/key'), ('review/key', 'review'),
                                           ('run/review', 'key'), ('review', 'run/key'), ('run', 'key')])
def test_private_key_and_source_are_separate(tmp_path, project, key_dir):
    with pytest.raises(ValueError):
        B.ensure_separate(tmp_path / project, tmp_path / key_dir, tmp_path / 'run')


def test_qupath_sanitizer_preserves_geometry_but_removes_predictions(tmp_path):
    qupath = Path(B.QUPATH)
    if not qupath.is_file():
        pytest.skip('QuPath is unavailable')
    assignments = tmp_path / 'assignments.json'
    assignments.write_text(json.dumps([{'object_id': '00000000-0000-0000-0000-000000000001',
                                      'sample_id': 'Review-000001',
                                      'blind_uuid': '11111111-1111-1111-1111-111111111111'}]))
    fixture = tmp_path / 'fixture.groovy'
    fixture.write_text(r'''
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClass
import qupath.lib.objects.hierarchy.PathObjectHierarchy
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
import qupath.lib.images.ImageData
import qupath.lib.images.servers.WrappedBufferedImageServer
import qupath.lib.io.PathIO
import qupath.lib.plugins.workflow.DefaultScriptableWorkflowStep
import java.awt.image.BufferedImage

def rect = { x, y, w, h -> ROIs.createRectangleROI(x, y, w, h, ImagePlane.getDefaultPlane()) }
def hierarchy = new PathObjectHierarchy()
def tissue = PathObjects.createAnnotationObject(rect(0,0,100,100), PathClass.fromString('Tissue'))
def islet = PathObjects.createAnnotationObject(rect(5,5,40,40), PathClass.fromString('Islet'))
def collagen = PathObjects.createAnnotationObject(rect(10,10,40,40), PathClass.fromString('Collagen IV'))
def cell = PathObjects.createCellObject(rect(12,12,8,8), rect(14,14,3,3), PathClass.fromString('SECRET_PREDICTION'), null)
def unused = PathObjects.createCellObject(rect(60,60,8,8), rect(62,62,3,3), PathClass.fromString('SECRET_UNUSED'), null)
cell.setID(UUID.fromString('00000000-0000-0000-0000-000000000001'))
cell.setName('SECRET_NAME')
cell.setPathClass(cell.getPathClass(), 0.95)
cell.setColor(0xFF0000)
cell.getMetadata().put('phenotype_label','SECRET_METADATA')
cell.getMeasurementList().put('SECRET_PREDICTED_PROBABILITY',0.95)
tissue.setName('SECRET_TISSUE_NAME')
tissue.setDescription('SECRET_DESCRIPTION')
tissue.getMetadata().put('prediction','SECRET_TISSUE_METADATA')
hierarchy.addObject(tissue)
tissue.addChildObjects([islet,collagen,unused])
islet.addChildObject(cell)
hierarchy.getRootObject().getMetadata().put('prediction','SECRET_ROOT_METADATA')
def server = new WrappedBufferedImageServer('raw-image',new BufferedImage(100,100,BufferedImage.TYPE_BYTE_GRAY))
def imageData = new ImageData(server,hierarchy,ImageData.ImageType.FLUORESCENCE)
imageData.setProperty('prediction','SECRET_PROPERTY')
imageData.getHistoryWorkflow().addStep(new DefaultScriptableWorkflowStep('SECRET_WORKFLOW','predict()'))
def beforeCell = cell.getROI().getGeometry().toText()
def beforeNucleus = cell.getNucleusROI().getGeometry().toText()
def beforeTissue = tissue.getROI().getGeometry().toText()
def destination = new File(args[3])
def entry = new Expando(getImageName: { -> 'Review image 001' }, saveImageData: { data -> PathIO.writeImageData(destination,data) })
def local = new Binding([args: [args[1],args[2]] as String[],getCurrentHierarchy:{ -> hierarchy },
                        getCurrentImageData:{ -> imageData },getProjectEntry:{ -> entry }])
new GroovyShell(local).evaluate(new File(args[0]))
def reloaded = PathIO.readImageData(destination,server)
def objects = reloaded.getHierarchy().getDetectionObjects()
assert objects.size() == 1
def result = objects[0]
assert result.getID().toString() == '11111111-1111-1111-1111-111111111111'
assert result.getName() == 'Review-000001'
assert result.getPathClass().toString() == 'Review cell'
assert Double.isNaN(result.getClassProbability())
assert result.getMetadata().isEmpty() && result.getMeasurementList().isEmpty()
assert result.getROI().getGeometry().toText() == beforeCell
assert result.getNucleusROI().getGeometry().toText() == beforeNucleus
def annotations = reloaded.getHierarchy().getAnnotationObjects()
assert annotations.size() == 1 && annotations[0].getPathClass().toString() == 'Tissue'
assert annotations[0].getROI().getGeometry().toText() == beforeTissue
assert annotations[0].getMetadata().isEmpty() && annotations[0].getName() == null
assert annotations[0].getDescription() == null
assert reloaded.getProperties().isEmpty()
// QuPath recreates its harmless "Set image type" step while loading; old prediction steps must be gone.
assert reloaded.getHistoryWorkflow().getSteps().every { it.getName() == 'Set image type' }
assert !new String(destination.bytes, 'ISO-8859-1').contains('SECRET_')
reloaded.close()
imageData.close()
println('blind saved/reloaded geometry and metadata checks passed')
''')
    report = tmp_path / 'report.json'
    command = [str(qupath), 'script', '-a', str(B.SCRIPT), '-a', str(assignments),
               '-a', str(report), '-a', str(tmp_path / 'sanitized.qpdata'), str(fixture)]
    result = subprocess.run(command, capture_output=True, text=True, timeout=60,
                            env={**os.environ, 'JAVA_TOOL_OPTIONS': '-Xmx1g -XX:ActiveProcessorCount=2'})
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'blind saved/reloaded geometry and metadata checks passed' in result.stdout
    assert json.loads(report.read_text())['passed'] is True
