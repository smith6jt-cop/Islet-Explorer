"""Resolve an explicitly reviewed structure generation without changing frozen cell runs.

An absent binding preserves the legacy layout. A present binding is fail-closed.
Only Lymphatic may change; the selected Nerve, Capillary and Collagen IV bytes must
equal the frozen cell run. Native historical annotations remain ancestry, while
the selected derived masks are used by context, endpoints and new review copies.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
import stat
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'typing'))
import codex_typing as CT

CLASSES = ('Nerve', 'Capillary', 'Lymphatic', 'Collagen IV')
BINDING = 'structure_generation.json'
COHORT = HERE / 'cohort.csv'
PREVIOUS_BASELINE = HERE.parents[1] / 'data/development/uniform-20donor-20261027'
_CACHE = {}
_HASH_CACHE = {}


def _stat(target):
    p = Path(target)
    s = p.stat()
    value = {'path': str(p.resolve()), 'size': s.st_size, 'mtime_ns': s.st_mtime_ns,
             'ctime_ns': s.st_ctime_ns, 'inode': s.st_ino}
    if stat.S_ISDIR(s.st_mode):
        value['entries_sha256'] = hashlib.sha256(json.dumps(sorted(x.name for x in p.iterdir()),
            ensure_ascii=True, separators=(',', ':')).encode()).hexdigest()
    return value


def _sha(target):
    before = _stat(target); key = tuple(before.values())
    if key not in _HASH_CACHE:
        h = hashlib.sha256()
        with Path(target).open('rb') as stream:
            for block in iter(lambda: stream.read(1 << 20), b''): h.update(block)
        if before != _stat(target): raise ValueError(f'Structure input changed while hashing: {target}')
        _HASH_CACHE[key] = h.hexdigest()
    return _HASH_CACHE[key]


_LOADED_SHA = _sha(__file__)


def _cohort():
    tags = [r['image_name'] for r in CT.cohort_rows()]
    if len(tags) != len(set(tags)): raise ValueError('Duplicate cohort image')
    return tags


class _Inventory:
    def __init__(self):
        self.files = {}; self.identities = {}

    def add(self, target, expected=None):
        p = Path(target)
        if not p.is_absolute(): raise ValueError(f'Structure generation requires an absolute path: {p}')
        before = _stat(p); digest = _sha(p)
        if before != _stat(p): raise ValueError(f'Structure input changed: {p}')
        if expected is not None and digest != expected: raise ValueError(f'Structure generation hash mismatch: {p}')
        self.files[str(p)] = {**before, 'sha256': digest}
        return p

    def pointer(self, item):
        if not isinstance(item, dict) or not isinstance(item.get('sha256'), str):
            raise ValueError('Missing hash-pinned structure-generation pointer')
        return self.add(item['path'], item['sha256'])

    def json(self, target):
        p = self.add(Path(target).absolute())
        value = json.loads(p.read_text())
        if self.files[str(p)] != {**_stat(p), 'sha256': _sha(p)}:
            raise ValueError(f'Structure record changed during reading: {p}')
        return value

    def records(self, value, relative_root=None):
        """Retain all file-backed ancestry already verified by the underlying validators."""
        if isinstance(value, list):
            for v in value: self.records(v, relative_root)
        elif isinstance(value, dict):
            if 'path' in value and 'sha256' in value:
                self.pointer(value)
            elif 'relative_path' in value and 'sha256' in value and relative_root is not None:
                self.add(relative_root / value['relative_path'], value['sha256'])
            elif all(k in value for k in ('path', 'size', 'mtime_ns')):
                actual = _stat(value['path'])
                if any(actual[k] != value[k] for k in actual if k in value):
                    raise ValueError(f'Source image identity changed: {value["path"]}')
                self.identities[value['path']] = _stat(value['path'])
            for k, v in value.items():
                if isinstance(k, str) and k.startswith('/') and isinstance(v, str) and re.fullmatch('[0-9a-f]{64}', v):
                    self.add(k, v)
                self.records(v, relative_root)

    def unchanged(self):
        for p, item in {**self.files, **self.identities}.items():
            current = _stat(p)
            if current != {k: item[k] for k in current}: return False
        return True

    def directory(self, target):
        """Guard file-set predicates (raw evidence globs and optional C inputs)."""
        p = Path(target)
        while not p.exists(): p = p.parent
        self.identities[str(p)] = _stat(p)


def _review(inv, pointer, tag, cls, geometry_sha, tissue_sha, quantitative_sha=None):
    p = inv.pointer(pointer); r = inv.json(p)
    if r.get('supersedes') is not None:
        prior = inv.json(inv.pointer(r['supersedes']))
        if prior.get('image') != tag or prior.get('class', cls) != cls:
            raise ValueError(f'{tag} {cls}: superseded review names another image/class')
    if r.get('image') != tag or r.get('class', cls) != cls:
        raise ValueError(f'{tag} {cls}: review names another image/class')
    if r.get('verdict', r.get('status')) not in ('acceptable', 'acceptable_with_notes', 'good'):
        raise ValueError(f'{tag} {cls}: actual review is not accepted')
    if r.get('structure_sha256') != geometry_sha or r.get('tissue_sha256') != tissue_sha:
        raise ValueError(f'{tag} {cls}: review geometry/Tissue is stale')
    views = r.get('views_reviewed', [])
    if len(views) != len(set(views)) or not any('whole' in Path(v).stem for v in views):
        raise ValueError(f'{tag} {cls}: missing whole-image actual review')
    native = [v for v in views if re.search(r'(?:tile|zoom)_?\d+', Path(v).stem)]
    if len(native) < 6: raise ValueError(f'{tag} {cls}: fewer than six actual native fields')
    for v in views:
        target = Path(v) if Path(v).is_absolute() else p.parent / v
        if v not in r.get('file_sha256', {}): raise ValueError(f'{tag} {cls}: unpinned reviewed view {v}')
        inv.add(target, r['file_sha256'][v])
    if quantitative_sha is not None:
        if r.get('quantitative_sha256') != quantitative_sha:
            raise ValueError(f'{tag}: Lymphatic quantitative review is stale')
        tile_path = inv.add(Path(r.get('tile_manifest_path', p.parent / f'{tag}_tiles.json')), r['tile_manifest_sha256'])
        document = inv.json(tile_path)
        positions = document
        if isinstance(document, dict):
            if document.get('schema') != 'bounded_lymphatic_combined_native_windows/1' or document.get('image') != tag:
                raise ValueError(f'{tag}: unsupported combined native-field manifest')
            original = inv.json(inv.pointer(document['original_manifest']))
            original_meta = inv.json(inv.pointer(document['original_metadata']))
            supplement = inv.json(inv.pointer(document['supplemental_metadata']))
            if supplement.get('schema') != 'bounded_lymphatic_supplemental_native/1' or supplement.get('image') != tag:
                raise ValueError(f'{tag}: unsupported supplemental native fields')
            if (document['original_metadata']['sha256'] != r.get('metadata_sha256')
                    or document['original_manifest']['sha256'] != original_meta.get('tile_manifest_sha256')
                    or supplement.get('source_inputs') != original_meta.get('source_inputs')
                    or supplement.get('original_manifest') != document['original_manifest']
                    or supplement.get('original_metadata') != document['original_metadata']
                    or r.get('supplemental_metadata_sha256') != document['supplemental_metadata']['sha256']):
                raise ValueError(f'{tag}: supplemental view ancestry differs')
            for field in ('structure_sha256', 'tissue_sha256', 'quantitative_sha256'):
                if any(obj.get(field) != r.get(field) for obj in (document, original_meta, supplement)):
                    raise ValueError(f'{tag}: supplemental geometry differs: {field}')
            inv.records(document); inv.records(supplement); inv.records(original_meta)
            positions = document['tiles']
            if positions != original + supplement['tiles']:
                raise ValueError(f'{tag}: combined manifest rewrites original fields')
            known_views = {**original_meta['file_sha256'], **supplement['file_sha256']}
            if set(original_meta['file_sha256']) & set(supplement['file_sha256']):
                raise ValueError(f'{tag}: supplemental view overwrites an original')
            if any(known_views.get(v) != r['file_sha256'][v] for v in views):
                raise ValueError(f'{tag}: actual supplemental view is not render-provenanced')
        required = {f'{tag}_tile{t["tile"]:02}.png' for t in positions}
        if len(required) != len(positions) or len(required) < 6 or not required <= {Path(v).name for v in native}:
            raise ValueError(f'{tag}: incomplete actual review of declared native fields')
        windows = {(t['x_px'], t['y_px'], t['width_px'], t['height_px']) for t in positions}
        if len(windows) < 6:
            raise ValueError(f'{tag}: fewer than six distinct actual native windows')
        if isinstance(document, dict):
            separated = []
            for t in positions:
                if all((t['x_px']-p['x_px'])**2 + (t['y_px']-p['y_px'])**2 >= 400**2 for p in separated):
                    separated.append(t)
            if (len(separated) < 6 or document.get('distinct_windows') != len(windows)
                    or document.get('separated_windows') != len(separated)):
                raise ValueError(f'{tag}: supplemental spatial coverage differs')
    else:
        inv.add(p.parent.parent / 'checks' / f'{tag}.csv', r['checks_sha256'])
        if 'tiles_metadata_sha256' in r:
            inv.add(p.parent / f'{tag}_tiles.json', r['tiles_metadata_sha256'])
    return r


def _lineage(inv, run, tag, pinned):
    import new_detection_filter as NF
    current = NF.accepted_inputs(run, tag)
    if current is None or current != pinned:
        raise ValueError(f'{tag}: structure generation cell lineage changed')
    inv.records(current)
    m = inv.json(Path(current['filter_manifest_path']))
    inv.json(run / 'filter/generation.json')
    inv.records(m, run)
    raw = Path(current['source_run'])
    for folder in (raw / 'augmentation' / tag, raw / 'provenance' / tag, raw / 'logs' / tag):
        inv.directory(folder)
    # Reviewed decision ancestry and all referenced evidence are part of the cache guard.
    decisions = inv.json(Path(m['decision_source']))
    inv.records(decisions['images'][tag])
    reuse = NF.structure_reuse_current(run, tag)
    if reuse is None: raise ValueError(f'{tag}: missing unchanged native structure ancestry')
    inv.records(reuse)
    certificate = inv.json(Path(reuse['certificate']))
    inv.records(certificate, run)
    # The underlying reuse validator checks these code/classifier files. Pin them
    # explicitly as well so a cached proof cannot outlive a classifier edit.
    import redraw_structures as R
    for source in (Path(NF.__file__), Path(R.__file__), R.SCRIPT, Path(R.ST.__file__), Path(CT.__file__)):
        inv.add(source.resolve())
    for name in ('Nerve', 'Capillary', 'Lymphatics'): inv.add((R.CLASSIFIERS / f'{name}.json').resolve())
    threshold = Path(current['source_run']) / 'structures/thresholded.json'
    if threshold.exists(): inv.add(threshold)
    inv.directory(threshold.parent)
    return current


def _validate_previous(inv, run, manifest, selection, tags):
    """Honor explicit reuse of original objects without inventing a new review claim."""
    import collagen_structures as C
    if selection.get('source_run') != str(run) or set(selection.get('images', {})) != set(tags):
        raise ValueError('Previous Lymphatic selection has a different cell run/cohort')
    baseline = Path(selection['baseline_run']).resolve()
    if baseline != PREVIOUS_BASELINE.resolve():
        raise ValueError('Previous Lymphatic selection is not the original baseline')
    instruction = inv.json(inv.pointer(selection['user_instruction']))
    if (instruction.get('decision') != 'reuse_exact_previous_lymphatic_objects'
            or instruction.get('cell_run') != str(run)
            or instruction.get('baseline_run') != str(baseline)
            or instruction.get('supersedes_reconstruction_gate') is not True
            or not isinstance(instruction.get('user_instruction'), str)
            or not instruction['user_instruction'].strip()):
        raise ValueError('Previous Lymphatic selection lacks the explicit user instruction')
    paths = {}
    for tag in tags:
        item = manifest['images'][tag]; chosen = selection['images'][tag]
        current = _lineage(inv, run, tag, item['cell_inputs'])
        inv.records(chosen)
        expected_tissue = run / 'inputs/tissue' / f'{tag}.geojson'
        if (Path(chosen['tissue']['path']).resolve() != expected_tissue.resolve()
                or chosen['tissue']['sha256'] != current['files']['tissue']['sha256']):
            raise ValueError(f'{tag}: previous Lymphatic selection uses another Tissue lineage')
        if set(item['classes']) != set(CLASSES):
            raise ValueError(f'{tag}: incomplete selected structure inventory')
        paths[tag] = {}
        for cls in CLASSES:
            selected = item['classes'][cls]; selected_path = inv.pointer(selected)
            paths[tag][cls] = str(selected_path)
            if cls == 'Lymphatic':
                expected = baseline / 'structures/Lymphatic' / f'Lymphatic__{tag}.geojson'
                if selected != chosen['classes']['Lymphatic'] or selected_path.resolve() != expected.resolve():
                    raise ValueError(f'{tag}: selected Lymphatic is not the exact previous object file')
                inv.add(expected, selected['sha256'])
                # The user's explicit reuse supersedes regeneration review. Do not
                # turn held reconstruction reviews into acceptance of these originals.
                continue
            inv.add(run / 'structures' / cls / f'{cls}__{tag}.geojson', selected['sha256'])
            if cls != 'Collagen IV':
                _review(inv, item['reviews'][cls], tag, cls, selected['sha256'], current['files']['tissue']['sha256'])
        sidecar = inv.json(inv.pointer(item['collagen_sidecar']))
        inv.add(C.checkpoint_path(run, tag), item['collagen_sidecar']['sha256'])
        C.validate_record(run, tag, sidecar)
        inv.records(sidecar); inv.directory(run / 'structures' / C.PREVIOUS)
        for source in (Path(C.__file__), Path(C.ST.__file__), Path(C.RS.__file__), Path(CT.__file__)):
            inv.add(source.resolve())
    return paths


def _validate_pdpn_gate(inv, run, manifest, gate, tags):
    """Owner 2026-10-30 (option 1): the reviewed re-run Lymphatic objects whose mean PDPN reaches the image's own
    PDPN level. The gate record proves the subset; the agent's actual review pins the gated geometry and views."""
    import collagen_structures as C
    if gate.get('source_run') != str(run) or set(gate.get('images', {})) != set(tags):
        raise ValueError('PDPN-gated Lymphatic selection has a different cell run/cohort')
    instruction = inv.json(inv.pointer(gate['owner_instruction']))
    if instruction.get('decision') != 'pdpn_level_object_gate' or instruction.get('cell_run') != str(run):
        raise ValueError('PDPN-gated Lymphatic selection lacks the owner decision')
    method = gate['method']
    inv.add(Path(method['code_path']), method['code_sha256'])
    inv.pointer(gate['levels'])
    paths = {}
    for tag in tags:
        item = manifest['images'][tag]; chosen = gate['images'][tag]
        current = _lineage(inv, run, tag, item['cell_inputs'])
        inv.records(chosen)
        if (Path(chosen['tissue']['path']).resolve() != (run / 'inputs/tissue' / f'{tag}.geojson').resolve()
                or chosen['tissue']['sha256'] != current['files']['tissue']['sha256']):
            raise ValueError(f'{tag}: PDPN-gated Lymphatic uses another Tissue lineage')
        if Path(chosen['source_objects']['path']).resolve() != (run / 'structures/Lymphatic' / f'Lymphatic__{tag}.geojson').resolve():
            raise ValueError(f'{tag}: PDPN gate source is not this run\'s re-run objects')
        if item.get('quantitative') != chosen['record']:
            raise ValueError(f'{tag}: PDPN gate record pointer differs')
        q = inv.json(inv.pointer(item['quantitative']))
        if q.get('image') != tag or q.get('method') != method['name'] or q.get('checks', {}).get('passed') is not True:
            raise ValueError(f'{tag}: PDPN gate record is missing or failed')
        if any(q['inputs'].get(k) != chosen[k] for k in ('tissue', 'source_objects', 'previous_objects')):
            raise ValueError(f'{tag}: PDPN gate record has different inputs')
        if set(item['classes']) != set(CLASSES): raise ValueError(f'{tag}: incomplete selected structure inventory')
        paths[tag] = {}
        for cls in CLASSES:
            selected = item['classes'][cls]; selected_path = inv.pointer(selected)
            paths[tag][cls] = str(selected_path)
            if cls == 'Lymphatic':
                if (selected != chosen['classes'][cls] or q.get('output_sha256') != selected['sha256']
                        or Path(q['output_path']).resolve() != selected_path.resolve()):
                    raise ValueError(f'{tag}: PDPN gate output differs from selected Lymphatic')
            else:
                inv.add(run / 'structures' / cls / f'{cls}__{tag}.geojson', selected['sha256'])
            if cls != 'Collagen IV':
                _review(inv, item['reviews'][cls], tag, cls, selected['sha256'], current['files']['tissue']['sha256'],
                        item['quantitative']['sha256'] if cls == 'Lymphatic' else None)
        sidecar = inv.json(inv.pointer(item['collagen_sidecar']))
        inv.add(C.checkpoint_path(run, tag), item['collagen_sidecar']['sha256'])
        C.validate_record(run, tag, sidecar)
        inv.records(sidecar); inv.directory(run / 'structures' / C.PREVIOUS)
        for source in (Path(C.__file__), Path(C.ST.__file__), Path(C.RS.__file__), Path(CT.__file__)):
            inv.add(source.resolve())
    return paths




def _normalized_polygons(path):
    """Geometry of a QuPath export independent of object order and IDs (lymphatic_despeckle_run.normalized)."""
    import shapely
    from shapely.geometry import shape
    return sorted(shapely.to_wkb(shapely.normalize(shape(f['geometry']))) for f in json.loads(Path(path).read_text())['features'])

def _validate_despeckled_gate(inv, run, manifest, gate, tags):
    """Owner 2026-11-01 ("Do dg"): the stored Lymphatics classifier run on despeckled PDPN (owner 2026-10-31, "the
    lymphatic object detection needs the despeckle"), then option 1's PDPN level gate. Each image's detection result
    pins the run's server, the classifier, the current Tissue and the despeckled PDPN, whose record names the source
    image; the original-mode reproduction of the reviewed re-run objects validates the detection script; the gate
    record proves the subset; the agent's review pins the gated geometry and views."""
    import collagen_structures as C
    if gate.get('source_run') != str(run) or set(gate.get('images', {})) != set(tags):
        raise ValueError('Despeckled Lymphatic selection has a different cell run/cohort')
    instruction = inv.json(inv.pointer(gate['owner_instruction']))
    if instruction.get('decision') != 'despeckled_detection_level_gate' or instruction.get('cell_run') != str(run):
        raise ValueError('Despeckled Lymphatic selection lacks the owner decision')
    method = gate['method']; detection = method['detection']
    inv.add(Path(method['code_path']), method['code_sha256'])
    inv.records(detection)
    if not detection.get('validation'): raise ValueError('Despeckled detection lacks its original-mode validation')
    for tag, v in detection['validation'].items():
        result = inv.json(inv.pointer(v['result']))
        produced = inv.pointer(result['output']); reviewed = inv.pointer(v['reviewed_rerun_objects'])
        if (result.get('mode') != 'original' or result.get('image') != tag
                or reviewed.resolve() != (run / 'structures/Lymphatic' / f'Lymphatic__{tag}.geojson').resolve()
                or _normalized_polygons(produced) != _normalized_polygons(reviewed)):
            raise ValueError(f'{tag}: the detection script does not reproduce the reviewed re-run objects')
    inv.pointer(gate['levels'])
    classifier = detection['classifier']
    sources = {e['image_name']: e['image_path'] for e in inv.json(run / 'project_manifest.json')['entries']}
    paths = {}
    for tag in tags:
        item = manifest['images'][tag]; chosen = gate['images'][tag]
        current = _lineage(inv, run, tag, item['cell_inputs'])
        inv.records(chosen)
        tissue_sha = current['files']['tissue']['sha256']
        if (Path(chosen['tissue']['path']).resolve() != (run / 'inputs/tissue' / f'{tag}.geojson').resolve()
                or chosen['tissue']['sha256'] != tissue_sha):
            raise ValueError(f'{tag}: despeckled Lymphatic uses another Tissue lineage')
        det = chosen['detection']
        result = inv.json(inv.pointer(det['result'])); pdpn = inv.json(inv.pointer(det['pdpn_record']))
        if (result.get('schema') != 'uniform.lymphatics_despeckled_result/1' or result.get('image') != tag
                or result.get('mode') != 'despeckled' or result.get('output') != chosen['source_objects']):
            raise ValueError(f'{tag}: the despeckled detection result does not produce the gate source')
        read = result.get('inputs', {})
        servers = [k for k in read if k.endswith('/server.json')]
        if (read.get(chosen['tissue']['path']) != tissue_sha or read.get(classifier['path']) != classifier['sha256']
                or read.get(det['pdpn_raw']['path']) != det['pdpn_raw']['sha256'] or len(servers) != 1
                or Path(servers[0]).resolve().parents[2] != (run / 'qupath').resolve()):
            raise ValueError(f'{tag}: the despeckled detection read another Tissue, classifier, PDPN or project')
        if (pdpn.get('image') != tag or pdpn.get('channel') != 'PDPN' or pdpn.get('path') != det['pdpn_raw']['path']
                or pdpn.get('sha256') != det['pdpn_raw']['sha256']
                or pdpn.get('despeckle_code_sha256') != detection['despeckle_code_sha256']
                or pdpn.get('source', {}).get('path') != sources.get(tag)):
            raise ValueError(f'{tag}: despeckled PDPN is not this image\'s despeckled source channel')
        inv.records(result); inv.records(pdpn)
        if item.get('quantitative') != chosen['record']:
            raise ValueError(f'{tag}: despeckled Lymphatic gate record pointer differs')
        q = inv.json(inv.pointer(item['quantitative']))
        if q.get('image') != tag or q.get('method') != method['name'] or q.get('checks', {}).get('passed') is not True:
            raise ValueError(f'{tag}: despeckled Lymphatic gate record is missing or failed')
        if any(q['inputs'].get(k) != chosen[k] for k in ('tissue', 'source_objects', 'previous_objects', 'detection')):
            raise ValueError(f'{tag}: despeckled Lymphatic gate record has different inputs')
        if set(item['classes']) != set(CLASSES): raise ValueError(f'{tag}: incomplete selected structure inventory')
        paths[tag] = {}
        for cls in CLASSES:
            selected = item['classes'][cls]; selected_path = inv.pointer(selected)
            paths[tag][cls] = str(selected_path)
            if cls == 'Lymphatic':
                if (selected != chosen['classes'][cls] or q.get('output_sha256') != selected['sha256']
                        or Path(q['output_path']).resolve() != selected_path.resolve()):
                    raise ValueError(f'{tag}: despeckled Lymphatic gate output differs from selected Lymphatic')
            else:
                inv.add(run / 'structures' / cls / f'{cls}__{tag}.geojson', selected['sha256'])
            if cls != 'Collagen IV':
                _review(inv, item['reviews'][cls], tag, cls, selected['sha256'], tissue_sha,
                        item['quantitative']['sha256'] if cls == 'Lymphatic' else None)
        sidecar = inv.json(inv.pointer(item['collagen_sidecar']))
        inv.add(C.checkpoint_path(run, tag), item['collagen_sidecar']['sha256'])
        C.validate_record(run, tag, sidecar)
        inv.records(sidecar); inv.directory(run / 'structures' / C.PREVIOUS)
        for source in (Path(C.__file__), Path(C.ST.__file__), Path(C.RS.__file__), Path(CT.__file__)):
            inv.add(source.resolve())
    return paths

def _result(inv, run, binding_path, binding, manifest_path, manifest, paths, extra=None):
    if not inv.unchanged(): raise ValueError('Structure generation inputs changed during cohort validation')
    provenance = {'schema': 'uniform.structure_generation_provenance/1', 'cell_run': str(run),
                  'binding_path': str(binding_path), 'binding_sha256': inv.files[str(binding_path)]['sha256'],
                  'manifest_path': str(manifest_path), 'manifest_sha256': binding['manifest']['sha256'],
                  'restoration_manifest_sha256': manifest['restoration']['sha256'], 'resolver_sha256': _LOADED_SHA,
                  **(extra or {})}
    return {'binding': binding, 'manifest': manifest, 'provenance': provenance, 'paths': paths,
            'files': inv.files, 'identities': inv.identities, '_inventory': inv}


def _validate(run):
    import collagen_structures as C
    inv = _Inventory()
    inv.add(Path(__file__).resolve(), _LOADED_SHA)
    binding_path = run / BINDING
    binding = inv.json(binding_path)
    if binding.get('schema') != 'uniform.structure_generation_binding/1' or binding.get('cell_run') != str(run):
        raise ValueError('Unsupported or misplaced structure-generation binding')
    manifest_path = inv.pointer(binding['manifest']); manifest = inv.json(manifest_path)
    if manifest.get('schema') != 'uniform.structure_generation/1' or manifest.get('cell_run') != str(run):
        raise ValueError('Unsupported or misplaced structure-generation manifest')
    inv.add(COHORT.resolve(), manifest['cohort_sha256'])
    tags = _cohort()
    if set(manifest.get('images', {})) != set(tags): raise ValueError('Structure generation must cover exactly the cohort')
    restoration_path = inv.pointer(manifest['restoration']); restoration = inv.json(restoration_path)
    if restoration.get('schema') == 'lymphatic_exact_previous_reuse/1':
        paths = _validate_previous(inv, run, manifest, restoration, tags)
        return _result(inv, run, binding_path, binding, manifest_path, manifest, paths,
                       {'lymphatic_selection': 'exact_previous_objects',
                        'user_instruction_sha256': restoration['user_instruction']['sha256']})
    if restoration.get('schema') == 'lymphatic_pdpn_level_gate/1':
        paths = _validate_pdpn_gate(inv, run, manifest, restoration, tags)
        return _result(inv, run, binding_path, binding, manifest_path, manifest, paths,
                       {'lymphatic_selection': 'pdpn_level_gate',
                        'owner_instruction_sha256': restoration['owner_instruction']['sha256']})
    if restoration.get('schema') == 'lymphatic_despeckled_level_gate/1':
        paths = _validate_despeckled_gate(inv, run, manifest, restoration, tags)
        return _result(inv, run, binding_path, binding, manifest_path, manifest, paths,
                       {'lymphatic_selection': 'despeckled_level_gate',
                        'owner_instruction_sha256': restoration['owner_instruction']['sha256']})
    if restoration.get('schema') != 'lymphatic_bounded_restoration/1' or restoration.get('source_run') != str(run):
        raise ValueError('Unsupported or misplaced Lymphatic restoration')
    if set(restoration.get('images', {})) != set(tags): raise ValueError('Restoration cohort differs')
    method = restoration['method']
    inv.records(method)
    for stem in ('code', 'original_edge_script'): inv.add(Path(method[stem + '_path']), method[stem + '_sha256'])
    paths = {}
    residues = ('baseline_missing_inside_final_tissue_um2', 'outside_final_tissue_um2',
                'added_in_old_tissue_outside_selected_edge_components_um2', 'unauthorized_addition_um2',
                'reconstructed_edge_geometry_missing_um2', 'delta_inference_missing_um2')
    for tag in tags:
        item = manifest['images'][tag]; restored = restoration['images'][tag]
        current = _lineage(inv, run, tag, item['cell_inputs'])
        inv.records(restored)
        if restored['tissue']['sha256'] != current['files']['tissue']['sha256']:
            raise ValueError(f'{tag}: restoration uses another Tissue')
        baseline = Path(restoration['baseline_run'])
        for pointer, expected in ((restored['tissue'], run / 'inputs/tissue' / f'{tag}.geojson'),
                                  (restored['current_classifier'], run / 'structures/Lymphatic' / f'Lymphatic__{tag}.geojson'),
                                  (restored['old_tissue'], baseline / 'inputs/tissue' / f'{tag}.geojson')):
            if Path(pointer['path']).resolve() != expected.resolve():
                raise ValueError(f'{tag}: restoration source path names another input')
        if Path(restored['original_layer']['source_path']).resolve() != (baseline / 'structures/Lymphatic' / f'Lymphatic__{tag}.geojson').resolve():
            raise ValueError(f'{tag}: original layer is not the declared baseline')
        if restored['original_layer']['sha256'] != restored['original_layer']['source_sha256']:
            raise ValueError(f'{tag}: native original layer bytes were changed')
        inv.add(Path(restored['original_layer']['source_path']), restored['original_layer']['source_sha256'])
        if item['quantitative'] != restored['record']:
            raise ValueError(f'{tag}: quantitative restoration pointer differs')
        q = inv.json(inv.pointer(item['quantitative'])); checks = q['checks']
        if q.get('method') != method['name']:
            raise ValueError(f'{tag}: quantitative proof uses another restoration method')
        for name in ('tissue', 'old_tissue', 'current_classifier'):
            if q['inputs'].get(name) != restored[name]:
                raise ValueError(f'{tag}: quantitative proof has different {name} inputs')
        if q['inputs'].get('baseline') != {'path': restored['original_layer']['source_path'],
                                         'sha256': restored['original_layer']['source_sha256']}:
            raise ValueError(f'{tag}: quantitative proof has a different original layer')
        tol = checks.get('tolerance_um2')
        if not isinstance(tol, (float, int)) or not math.isfinite(tol) or not 0 <= tol <= 1e-7:
            raise ValueError(f'{tag}: invalid quantitative tolerance')
        if any(checks.get(k) is not True for k in ('passed', 'original_layer_bytes_equal', 'source_files_unchanged')):
            raise ValueError(f'{tag}: restoration quantitative proof failed')
        if any(not isinstance(checks.get(k), (float, int)) or not math.isfinite(checks[k])
               or not 0 <= checks[k] <= tol for k in residues):
            raise ValueError(f'{tag}: restoration geometry invariant failed')
        inv.records(q)
        if set(item['classes']) != set(CLASSES): raise ValueError(f'{tag}: incomplete selected structure inventory')
        paths[tag] = {}
        for cls in CLASSES:
            selected = item['classes'][cls]; selected_path = inv.pointer(selected)
            paths[tag][cls] = str(selected_path)
            if cls == 'Lymphatic':
                if selected != restored['classes'][cls] or q.get('image') != tag or q.get('output_sha256') != selected['sha256'] or Path(q['output_path']) != selected_path:
                    raise ValueError(f'{tag}: restoration output differs from selected Lymphatic')
            else:
                inv.add(run / 'structures' / cls / f'{cls}__{tag}.geojson', selected['sha256'])
            if cls != 'Collagen IV':
                review = _review(inv, item['reviews'][cls], tag, cls, selected['sha256'], current['files']['tissue']['sha256'],
                                 item['quantitative']['sha256'] if cls == 'Lymphatic' else None)
                if cls == 'Lymphatic':
                    metadata = inv.json(inv.pointer(restored['renders']))
                    if metadata.get('image') != tag or metadata.get('source_inputs') != q['inputs']:
                        raise ValueError(f'{tag}: Lymphatic render source inputs differ')
                    if review.get('metadata_sha256') != restored['renders']['sha256']:
                        raise ValueError(f'{tag}: Lymphatic render provenance is stale')
                    for field in ('structure_sha256', 'tissue_sha256', 'quantitative_sha256'):
                        if metadata.get(field) != review.get(field):
                            raise ValueError(f'{tag}: Lymphatic render/review mismatch: {field}')
                    tile_path = Path(review.get('tile_manifest_path', Path(item['reviews'][cls]['path']).parent / f'{tag}_tiles.json'))
                    tile_document = inv.json(tile_path)
                    original_tile_sha = (tile_document['original_manifest']['sha256'] if isinstance(tile_document, dict)
                                         else review['tile_manifest_sha256'])
                    if metadata.get('tile_manifest_sha256') != original_tile_sha:
                        raise ValueError(f'{tag}: original native field provenance differs')
                    if any(metadata['file_sha256'][v] != review['file_sha256'].get(v) for v in metadata['file_sha256']):
                        raise ValueError(f'{tag}: Lymphatic render views differ from actual review')
                    inv.records(metadata)
        sidecar_path = inv.pointer(item['collagen_sidecar'])
        inv.add(C.checkpoint_path(run, tag), item['collagen_sidecar']['sha256'])
        sidecar = inv.json(sidecar_path)
        C.validate_record(run, tag, sidecar)
        # C provenance has explicit input paths and code hashes; retain them all.
        inv.records(sidecar)
        inv.directory(run / 'structures' / C.PREVIOUS)
        for source in (Path(C.__file__), Path(C.ST.__file__), Path(C.RS.__file__), Path(CT.__file__)):
            inv.add(source.resolve())
    return _result(inv, run, binding_path, binding, manifest_path, manifest, paths)


def validate(run, tags=None):
    """Prove the complete bound cohort once; invalidate cache on any tracked-file change.

    Call before and after a stage. Path resolution reuses this proof with stat
    checks, never repeated full-cohort hashing. No binding returns None.
    """
    run = Path(run).resolve(); binding = run / BINDING
    if _sha(__file__) != _LOADED_SHA: raise ValueError('Structure resolver changed since import; restart')
    if not binding.exists() and not binding.is_symlink():
        _CACHE.pop(str(run), None)
        return None
    cached = _CACHE.get(str(run))
    try:
        current = cached is not None and cached['_inventory'].unchanged()
    except OSError:
        current = False
    if not current:
        _CACHE.pop(str(run), None)
        cached = _validate(run); _CACHE[str(run)] = cached
    if tags is not None and (len(tags) != len(set(tags)) or set(tags) - set(cached['paths'])):
        raise ValueError('Requested images are outside the bound structure cohort')
    return {k: v for k, v in cached.items() if k != '_inventory'}


def provenance(run):
    result = validate(run)
    return None if result is None else result['provenance']


def path(run, tag, cls):
    if cls not in CLASSES: raise ValueError(f'Unknown structure class: {cls}')
    result = validate(run, [tag])
    return Path(run) / 'structures' / cls / f'{cls}__{tag}.geojson' if result is None else Path(result['paths'][tag][cls])
