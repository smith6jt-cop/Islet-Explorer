#!/usr/bin/env python3
"""Stage the selected source project and compatible copies of its analysis scripts."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil

SOURCE = Path('/data/QuPath_Projects/QuPath_project_IsletExplorer')
SELECTED = {38: '6539', 39: '6090', 40: '6520', 42: '0078', 43: '6232', 44: '6558'}


def compatible(text):
    # QuPath 0.7 replaced the legacy write methods with MeasurementList.put.
    text = re.sub(r"if \(ml\.metaClass\.respondsTo\(ml,\s*'addMeasurement',\s*String,\s*double\)\)",
                  "if (ml.metaClass.respondsTo(ml, 'put', String, double))", text)
    text = re.sub(r"if \(ml\.metaClass\.respondsTo\(ml,\s*'putMeasurement',\s*String,\s*double\)\)",
                  "if (ml.metaClass.respondsTo(ml, 'put', String, double))", text)
    text = text.replace('.addMeasurement(', '.put(').replace('.putMeasurement(', '.put(').replace('.setMeasurement(', '.put(')
    # Explicit comparator works across Groovy versions; its ordering is unchanged.
    text = text.replace(".sort { a -> [a.getROI().getBoundsX(), a.getROI().getBoundsY()] }",
                        ".sort { a, b -> (a.getROI().getBoundsX() <=> b.getROI().getBoundsX()) ?: (a.getROI().getBoundsY() <=> b.getROI().getBoundsY()) }")
    return text


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, default=Path('data/development/expanded-20260911'))
    args = p.parse_args()
    root = args.root.resolve()
    project = root / 'qupath'
    project.mkdir(parents=True, exist_ok=True)
    source_json = json.loads((SOURCE / 'project.qpproj').read_text())
    source_json['images'] = [e for e in source_json['images'] if e['entryID'] in SELECTED]
    def fix_uris(obj):
        if isinstance(obj, dict):
            for key, value in obj.items():
                if key == 'uri' and isinstance(value, str) and '/Images/' in value:
                    obj[key] = value.replace('/Images/', '/Images/New/')
                else:
                    fix_uris(value)
        elif isinstance(obj, list):
            for value in obj:
                fix_uris(value)
    fix_uris(source_json)
    source_json['uri'] = (project / 'project.qpproj').as_uri()
    if not (project / 'project.qpproj').exists():
        (project / 'project.qpproj').write_text(json.dumps(source_json, indent=2) + '\n')
    shutil.copytree(SOURCE / 'classifiers', project / 'classifiers', dirs_exist_ok=True)
    for entry in source_json['images']:
        data_dir = project / 'data' / str(entry['entryID'])
        data_dir.mkdir(parents=True, exist_ok=True)
        for name in ('data.qpdata', 'summary.json'):
            original = SOURCE / 'data' / str(entry['entryID']) / name
            if original.exists() and not (data_dir / name).exists():
                shutil.copy2(original, data_dir / name)
    stage = project / 'scripts'
    stage.mkdir(exist_ok=True)
    provenance = []
    composition_source = SOURCE / 'scripts/3_Islet_Counter.groovy'
    provenance.append({'role': 'composition writer retained from preceding counter',
                       'source': str(composition_source),
                       'sha256': hashlib.sha256(composition_source.read_bytes()).hexdigest()})
    for name in ('1_Panc_Objects.groovy', '3_Islet_Counter_v2.groovy', '4_Islet_Cells.groovy', '5_Islet_Objects.groovy'):
        original = SOURCE / 'scripts' / name
        content = original.read_text()
        if name.startswith('1_'):
            content = content.replace("== 'Islets'", "== 'Islet'")
        if name.startswith('4_'):
            scope_setup = '''
// Restore the Tissue selection that is not serialized between CLI stages.
def scopeCellsBefore = getCellObjects()
def scopeIdsBefore = scopeCellsBefore.collect { it.getID() }.toSet()
def scopeParentsBefore = scopeCellsBefore.collectEntries { [(it.getID()):it.getParent()?.getID()] }
def scopeGeneratedClasses = ['Nerve', 'Capillary', 'Lymphatic', 'TissueSample'] as Set
def scopePrevious = getAnnotationObjects().findAll { scopeGeneratedClasses.contains(it.getPathClass()?.toString()) }
if (!scopePrevious.isEmpty()) removeObjects(scopePrevious, true)
if (getCellObjects().collect { it.getID() }.toSet() != scopeIdsBefore)
    throw new IllegalStateException('Cell identities changed while clearing structural annotations')
if (getCellObjects().collectEntries { [(it.getID()):it.getParent()?.getID()] } != scopeParentsBefore)
    throw new IllegalStateException('Cell parents changed while clearing structural annotations')
scopeCellsBefore.each { cell ->
    def ml = cell.getMeasurementList()
    def stale = ml.getNames().findAll { it.startsWith('Distance to detection ') || it.startsWith('Signed distance to annotation ') }
    if (!stale.isEmpty()) ml.removeAll(stale as String[])
}
selectObjectsByClassification('Tissue')
if (getSelectedObjects().isEmpty()) throw new IllegalStateException('No Tissue region for structural classifiers')
'''
            content = content.replace('// --- Your detections & housekeeping ---', scope_setup + '\n// --- Your detections & housekeeping ---')
            for classifier in ('Nerve', 'Capillary', 'Lymphatics'):
                content = content.replace(f'createAnnotationsFromPixelClassifier("{classifier}"',
                                          f'selectObjectsByClassification("Tissue")\ncreateAnnotationsFromPixelClassifier("{classifier}"')
            content += '''
def scopeCellsAfter = getCellObjects()
if (scopeCellsAfter.collect { it.getID() }.toSet() != scopeIdsBefore)
    throw new IllegalStateException('Cell identities changed during structural measurements')
if (scopeCellsAfter.collectEntries { [(it.getID()):it.getParent()?.getID()] } != scopeParentsBefore)
    throw new IllegalStateException('Cell parents changed during structural measurements')
new File(System.getenv('ISLET_EXPORT_DIR'), 'structural_scope_provenance.json').text =
    qupath.lib.io.GsonTools.getInstance(true).toJson([scope:'Tissue annotations',
        classifiers:['Nerve', 'Capillary', 'Lymphatics'],
        previous_generated_annotations_removed:scopePrevious.size(),
        cell_count_before:scopeCellsBefore.size(), cell_count_after:scopeCellsAfter.size(),
        cell_uuids_preserved:true, cell_parents_preserved:true,
        finished_at:java.time.Instant.now().toString()])
'''
        if name.startswith('3_'):
            # v2 repeats the endocrine block and omits the composition writer.
            start = content.index('// ───────────────────── COLOR CELLS & ENSURE CONSISTENCY')
            second = content.index('// ───────────────────── COLOR CELLS & ENSURE CONSISTENCY', start + 1)
            content = content[:start] + content[second:]
            old = (SOURCE / 'scripts/3_Islet_Counter.groovy').read_text()
            composition = old[old.index('// ───────────────────── ISLET COMPOSITION & TSV'):]
            content += '\n' + composition
            content = content.replace('/home/smith6jt/QuPath/v0.6/downloaded/fluorescence_nuclei_and_cells-0.1.0',
                                      '/home/smith6jt/QuPath/v0.7/instanseg/downloaded/fluorescence_nuclei_and_cells-0.1.1')
            content = content.replace("String  instansegDevice    = 'gpu'", "String  instansegDevice    = System.getenv('ISLET_INSTANSEG_DEVICE') ?: 'cpu'")
            content = content.replace('int     instansegThreads   = 32', "int     instansegThreads   = Integer.parseInt(System.getenv('ISLET_INSTANSEG_THREADS') ?: '32')")
            content = content.replace('System.setProperty("ai.djl.offline", "false")',
                                      'System.setProperty("ai.djl.offline", "true")\n'
                                      'System.setProperty("PYTORCH_LIBRARY_PATH", System.getenv("PYTORCH_LIBRARY_PATH") ?: "/home/smith6jt/.djl.ai/pytorch/2.5.1-cpu-linux-x86_64")\n'
                                      'System.setProperty("ai.djl.pytorch.num_interop_threads", "1")\n'
                                      'System.setProperty("ai.djl.pytorch.num_threads", "1")')
            content = content.replace('if (allCells2.isEmpty()) { println "No cell objects found; skipping coloring & composition."; return }',
                                      'if (allCells2.isEmpty()) throw new IllegalStateException("InstanSeg returned no tissue cells")')
            content += '''
def segmentationProvenance = [model:'fluorescence_nuclei_and_cells', model_version:'0.1.1',
    model_sha256:'636833d26b5afcc03bde1357b3e669cf015e91dd1edda863aa51c19eee986562',
    extension_version:'0.1.7', pytorch_version:ai.djl.engine.Engine.getInstance().getVersion(),
    device:instansegDevice, tile_workers:instansegThreads, tile_size:512, tile_padding:32,
    source_filename:getProjectEntry().getImageName(), segmentation_region:'Tissue',
    cell_count:getCellObjects().size(), islet_count:islets.size(),
    calibration_um:server.getPixelCalibration().getAveragedPixelSizeMicrons()]
new File(System.getenv('ISLET_EXPORT_DIR'), 'segmentation_provenance.json').text =
    qupath.lib.io.GsonTools.getInstance(true).toJson(segmentationProvenance)
'''
        content = compatible(content)
        if not (stage / name).exists() or (stage / name).read_text() != content:
            (stage / name).write_text(content)
        provenance.append({'source': str(original), 'sha256': hashlib.sha256(original.read_bytes()).hexdigest(),
                           'staged': str(stage / name), 'staged_sha256': hashlib.sha256(content.encode()).hexdigest()})
    (project / 'processing_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(project / 'project.qpproj')


if __name__ == '__main__':
    main()
