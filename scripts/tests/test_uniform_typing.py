"""Tests for the final typing chain (final typing plan of 2026-10-28, Step 6): the registry, the canary status logic, and
identity.py's rules on a small synthetic image built with the run's own registry, plus a check that no module of the
chain names a marker in code (owner 2026-10-28: marker roles are data, panel- and tissue-specific)."""
from __future__ import annotations

import ast
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

TYPING = Path(__file__).resolve().parents[1] / 'uniform' / 'typing'
sys.path.insert(0, str(TYPING))
import registry as RG  # noqa: E402
import identity as ID  # noqa: E402
import levels as LV  # noqa: E402
import area_coverage as AC  # noqa: E402

CHAIN = ('registry.py', 'levels.py', 'canary.py', 'identity.py', 'acceptance_checks.py')


# -- the registry ----------------------------------------------------------------------------------------------------------
def test_registry_typing_markers_are_on_every_image_with_one_allowed_measurement():
    reg = RG.load()
    for m in reg.measured_markers():
        assert reg.table.loc[m, 'on_every_image'] == 'yes', m
        assert reg.measurement(m) in RG.COMPARTMENT_COLUMN, m
        for c in reg.candidates(m):
            assert c in RG.COMPARTMENT_COLUMN, (m, c)
    for m in reg.table.index:
        if reg.table.loc[m, 'on_every_image'] != 'yes':
            assert reg.role(m) == 'never_used', m
    assert reg.role('ColIV') == 'none'                      # owner 2026-10-28: excluded from typing completely


def test_labels_file_names_registry_lineages_and_markers():
    reg = RG.load()
    lab = ID.labels_table(reg)
    lineages = {f'{lin}:{s}' if s else lin for m in reg.typing_markers() for lin, s in reg.lineage_entries(m)}
    lineages |= {lin for m in reg.typing_markers() for lin, _ in reg.lineage_entries(m)}
    for r in lab.itertuples():
        if r.kind == 'structure':            # a structure object's cell type, named by the object's marker
            assert r.named_by in reg.typing_markers() and r.structure, r
            continue
        assert r.lineage in lineages, r.lineage
        for m in RG.split(r.profile.replace('|', ';')):
            assert m in reg.table.index, m


def test_negative_quantile():
    assert LV.negative_quantile('negative_p99') == 99.0
    assert LV.negative_quantile('negative_p999') == pytest.approx(99.9)
    assert LV.negative_quantile('peak_dense') is None


def test_no_module_of_the_chain_names_a_marker():
    reg = RG.load()
    names = set(reg.table.index) | set(reg.table['channel'])
    for f in CHAIN:
        tree = ast.parse((TYPING / f).read_text())
        docstrings = {id(n.body[0].value) for n in ast.walk(tree)
                      if isinstance(n, (ast.Module, ast.FunctionDef, ast.ClassDef)) and n.body
                      and isinstance(n.body[0], ast.Expr) and isinstance(n.body[0].value, ast.Constant)}
        found = sorted({n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)
                        and id(n) not in docstrings and n.value in names})
        assert not found, f'{f} names markers in code: {found}'


# -- the canary status (Xenium notebook 02 cell 22, margin between groups) ------------------------------------------------
def test_group_status():
    det = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 0], [1, 1, 1], [1, 1, 0], [1, 1, 0], [1, 1, 0]], bool)
    depth = np.array([[0, 0, 0], [2, 0, 0], [3, 1, 0], [2, 2, 0], [3, 1, 1], [2, 2, 0], [2, 2, 0], [1, 1, 0]])
    # equal depth: the group both brighter and covering more of the cell is placed (owner 2026-10-29)
    bright = np.array([[0, 0, 0]] * 5 + [[1, 2, 0], [2, 1, 0], [1, 2, 0]], float)
    cover = np.array([[0, 0, 0]] * 5 + [[.2, .9, 0], [.2, .9, 0], [.2, .9, 0]], float)
    status, top, runner, td, rd, ng, settled = ID.group_status(det, depth, np.zeros(det.shape), 1, 2, bright, cover, 0.15)
    assert list(status) == ['none', 'single', 'resolved', 'lineage_ambiguous', 'doublet_suspected', 'resolved',
                            'lineage_ambiguous', 'lineage_ambiguous']   # brighter and covers more / split / under depth 2
    assert list(ng) == [0, 1, 2, 2, 3, 2, 2, 2] and top[2] == 0 and runner[2] == 1
    assert top[5] == 1 and runner[5] == 0 and list(settled) == [False] * 5 + [True, False, False]


@pytest.mark.parametrize('lead,brighter,has_interior,expected', [
    (0.149999, True, True, False), (0.15, True, True, True), (0.150001, True, True, True),
    (0.20, False, True, False), (0.20, True, False, False)])
def test_interior_tie_boundary_and_required_evidence(lead, brighter, has_interior, expected):
    bright = np.array([[2., 1.]]) if brighter else np.array([[1., 2.]])
    result = ID.group_status(np.ones((1, 2), bool), np.full((1, 2), 2), np.zeros((1, 2)), 1, 2,
                             bright, np.array([[0.20 + lead, 0.20]]), 0.15, np.array([has_interior]))
    assert bool(result[-1][0]) is expected
    assert result[0][0] == ('resolved' if expected else 'lineage_ambiguous')


# -- identity on a synthetic image --------------------------------------------------------------------------------------
LEVEL_RAW, CANARY_RAW, HIGH = 100.0, 1000.0, 2000.0
CELLS = [  # (expected label, start class, {marker: raw value}, established endocrine call); every other marker 0
    ('CD8a Tcell', None, {'CD3e': HIGH, 'CD45': HIGH, 'CD8a': HIGH}, None),
    ('CD4 Tcell', None, {'CD3e': HIGH, 'CD45': HIGH, 'CD8a': HIGH, 'CD4': 2 * HIGH}, None),   # both: the higher
    ('T cell', None, {'CD3e': HIGH, 'CD45': HIGH}, None),                                     # neither: interim status
    ('Hematopoietic CD45-', None, {'CD20': HIGH}, None),                                      # CD20 alone: never a B cell
    ('B cell', None, {'CD20': HIGH, 'CD45': HIGH}, None),
    ('CD68+CD163-Macrophage_LGALS3-CD4-', None, {'CD68': HIGH, 'CD45': HIGH}, None),
    ('Hematopoietic CD45+', None, {'CD45': HIGH}, None),
    ('Conflict T/Myeloid', None, {'CD163': HIGH, 'CD8a': 2 * HIGH}, None),                   # no CD8a+CD163+ cell type
    ('Exocrine', None, {'ECAD': HIGH, 'panCK': HIGH}, None),
    ('Exocrine', 'Islet', {'ECAD': HIGH}, None),                                              # exocrine inside the border (10-29)
    ('Beta cell', 'Islet', {'INS': HIGH, 'ECAD': HIGH}, None),                                # hormone, whole cell
    ('Beta cell', 'Islet', {'INS': HIGH, 'GCG': 0.8 * HIGH}, None),                           # the highest scaled hormone
    ('CHGA+ hormone-negative', 'EndocrineCluster', {'CHGA': HIGH}, None),
    ('EndocrineSingle', 'EndocrineSingle', {'CD31': HIGH}, None),                             # the object names the cell
    ('Alpha cell', 'EndocrineSingle', {'GCG': HIGH, 'CD31': HIGH}, None),                     # a verified single stays
    ('Exocrine', None, {'INS': HIGH, 'ECAD': HIGH}, None),                                    # no endocrine outside objects
    ('Beta cell', 'Islet', {'INS': HIGH}, None),                                              # hormone in an islet
    ('T cell', 'Islet', {'CD3e': HIGH, 'CD45': HIGH}, ('Immune', '')),                        # a T cell inside an islet
    ('Tissue', None, {}, None),                                                               # no marker above background
    ('Tissue', None, {'ColIV': 1e6}, None),                                                   # ColIV is never read
    ('Endothelial', None, {'CD31': HIGH}, None),                                              # lineage marker
    ('Tissue', None, {'BCatenin': HIGH}, None),                                               # only a shared marker
    ('Endothelium Blood Vessel', None, {'CD31': 2 * HIGH, 'SMA': HIGH}, None),                # in a vessel object
    ('Conflict Mesenchymal/Endothelial', None, {'CD31': HIGH, 'SMA': 2 * HIGH}, None),       # CD31 not on top: passes through
    ('Endothelial', None, {'CD31': HIGH}, None),                                              # under half its area inside
    ('Islet', 'Islet', {'INS': HIGH, 'GCG': HIGH}, None),                                     # exact tie: no subtype
    ('CD8a Tcell', None, {'CD3e': HIGH, 'CD45': HIGH, 'CD8a': HIGH, 'CD163': HIGH}, None),     # T leads Myeloid by one
    ('T cell', None, {'CD3e': 2 * HIGH, 'CD45': HIGH, 'CD163': HIGH}, None),                  # tie: brighter, covers more
    ('Conflict T/Myeloid', None, {'CD3e': 2 * HIGH, 'CD45': HIGH, 'CD163': HIGH}, None),      # tie: brighter, covers less
    ('Hematopoietic CD45-', None, {'CD4': HIGH, 'CD8a': HIGH}, None),                         # no lineage, CD45 negative
]
# area coverage (area_coverage.py): a positive marker covers 0.8 of the cell unless given here; index -> {marker: coverage}
COVERAGE = {27: {'CD3e': 0.9, 'CD163': 0.2}, 28: {'CD3e': 0.2, 'CD163': 0.9}}
# cells in a structure object: index -> (class, area fraction inside)
IN_STRUCTURE = {22: ('Blood vessel', 0.9), 23: ('Blood vessel', 0.9), 24: ('Blood vessel', 0.3)}


@pytest.fixture()
def synthetic_run(tmp_path):
    reg = RG.load()
    markers = reg.measured_markers()
    tag = 'SYN_0001'
    typ = tmp_path / 'typing'
    calls = tmp_path / ID.ESTABLISHED / 'calls'
    for sub in (typ / 'redsea', typ / 'context', calls):
        sub.mkdir(parents=True)
    n = len(CELLS)
    ids = [f'c{i}' for i in range(n)]
    rs = {'object_id': ids, 'qc_analysis_eligible': [True] * n}
    for m in reg.table.index:
        for meas in RG.COMPARTMENT_COLUMN:
            rs[reg.column(m, meas)] = [cell[2].get(m, 0.0) for cell in CELLS]
    pd.DataFrame(rs).to_parquet(typ / 'redsea' / f'{tag}.parquet')
    (typ / 'coverage').mkdir()
    cover = {m: [COVERAGE.get(i, {}).get(m, 0.8 if c[2].get(m, 0.0) > LEVEL_RAW else 0.0) for i, c in enumerate(CELLS)]
             for m in reg.table.index}
    coverage = pd.DataFrame({'object_id': ids, AC.INTERIOR_COUNT: 10,
                            **{ID.COV_PREFIX + m: v for m, v in cover.items()},
                            **{ID.COV_INTERIOR_PREFIX + m: v for m, v in cover.items()}})
    pd.DataFrame({'object_id': ids, 'x_um': np.arange(n) * 100.0, 'y_um': 0.0,
                  'context': ['islet' if c[1] else 'tissue' for c in CELLS],
                  'endocrine_object_class': [c[1] for c in CELLS],
                  **{col: [IN_STRUCTURE.get(i, ('', 0.0))[0] == s if k == 0 else
                           (IN_STRUCTURE[i][1] if IN_STRUCTURE.get(i, ('',))[0] == s else np.nan) for i in range(n)]
                     for s, cols in ID.STRUCTURE_COLUMNS.items() for k, col in enumerate(cols)}
                  }).to_parquet(typ / 'context' / f'{tag}.parquet')
    pd.DataFrame({'object_id': ids, 'cell_class': [(c[3] or ('Exocrine', ''))[0] for c in CELLS],
                  'excluded_reason': [(c[3] or ('', ''))[1] for c in CELLS]}).to_parquet(calls / f'{tag}.parquet')
    lv = pd.DataFrame({'image': tag, 'marker': markers, 'measurement': [reg.measurement(m) for m in markers],
                       'level_log': np.log1p(LEVEL_RAW)})
    lv.to_csv(typ / 'levels.csv', index=False)
    pd.DataFrame({'marker': markers, 'decision': 'use'}).to_csv(typ / 'levels_selection.csv', index=False)
    lv.assign(canary_log=np.log1p(CANARY_RAW))[['image', 'marker', 'canary_log']].to_csv(
        typ / 'canary_levels.csv', index=False)
    # Tiny provenance sources; identity tests do not render pixels from this synthetic image.
    image = tmp_path / 'image.tiff'
    image.write_bytes(b'synthetic-image')
    (tmp_path / 'project_manifest.json').write_text(json.dumps({'entries': [
        {'image_name': tag, 'entry_id': 1, 'image_path': str(image)}]}))
    server = tmp_path / 'qupath' / 'data' / '1'
    server.mkdir(parents=True)
    (server / 'server.json').write_text(json.dumps({'metadata': {'channels': [{'name': m} for m in markers]}}))
    (tmp_path / 'cell_geojson').mkdir()
    (tmp_path / 'cell_geojson' / f'{tag}.geojson').write_text('{}')
    AC.write_coverage(tmp_path, tag, coverage)
    return tmp_path, tag


def test_identity_rules(synthetic_run):
    run, tag = synthetic_run
    out, audit = ID.classify_image(run, tag)
    got = out['label'].tolist()
    want = [c[0] for c in CELLS]
    assert got == want, [(i, g, w) for i, (g, w) in enumerate(zip(got, want)) if g != w]
    assert out.iloc[2]['subtype_status'] == 'T_subtype_not_shown'
    assert out.iloc[11]['named_by'] == 'INS' and out.iloc[25]['status'] == 'tied_highest_hormones'
    assert out.iloc[13]['status'] == 'no_endocrine_marker_above_level'
    assert out.iloc[10]['named_by'] == 'INS' and out.iloc[14]['named_by'] == 'GCG'
    assert 'endocrine_outside_object' in out.iloc[15]['tags']
    assert out.iloc[18]['status'] == out.iloc[19]['status'] == 'no_marker_above_background'
    assert not any(str(x).startswith('ColIV') for x in got)
    # no cell is named by a marker (owner 2026-10-28): the start class stays, with the status saying why
    assert out.iloc[3]['status'] == 'no_lineage_named' and out.iloc[21]['status'] == 'shared_marker_only'
    assert out.iloc[22]['named_by'] == 'CD31' and 'passed_through_Blood_vessel' in out.iloc[23]['tags']
    assert not any(str(x)[-1:] in '+-' and ' ' not in str(x) and 'Macrophage_' not in str(x) for x in got)


@pytest.mark.parametrize('damage', ['missing_column', 'missing_id', 'extra_id', 'duplicate_id'])
def test_coverage_rejects_incomplete_or_nonunique_data(synthetic_run, damage):
    run, tag = synthetic_run
    frame = pd.read_parquet(run / 'typing' / 'coverage' / f'{tag}.parquet')
    ids = frame['object_id'].copy()
    markers = ID.used_markers(run, RG.load())
    if damage == 'missing_column':
        frame = frame.drop(columns=[AC.INTERIOR_PREFIX + markers[0]])
    elif damage == 'missing_id':
        frame = frame.iloc[1:]
    elif damage == 'extra_id':
        frame.loc[0, 'object_id'] = 'not-a-source-cell'
    else:
        frame.loc[0, 'object_id'] = frame.loc[1, 'object_id']
    with pytest.raises(ValueError, match='missing columns|object IDs'):
        AC.validate_frame(frame, ids, markers, tag)


@pytest.mark.parametrize('source', ['levels', 'geometry', 'image', 'provenance', 'parquet'])
def test_identity_rejects_missing_or_stale_coverage(synthetic_run, source):
    run, tag = synthetic_run
    if source == 'levels':
        path = run / 'typing' / 'levels.csv'
        frame = pd.read_csv(path)
        frame.loc[0, 'level_log'] += 0.1
        frame.to_csv(path, index=False)
    elif source == 'geometry':
        (run / 'cell_geojson' / f'{tag}.geojson').write_text('{"changed":true}')
    elif source == 'image':
        (run / 'image.tiff').write_bytes(b'changed-source-image')
    elif source == 'provenance':
        (run / 'typing' / 'coverage' / f'{tag}.json').unlink()
    else:
        path = run / 'typing' / 'coverage' / f'{tag}.parquet'
        frame = pd.read_parquet(path)
        frame.loc[0, AC.INTERIOR_COUNT] = 11
        frame.to_parquet(path, index=False)
    with pytest.raises(ValueError, match='stale coverage|provenance'):
        ID.classify_image(run, tag)


def test_no_interior_is_explicit_and_cannot_settle_identity_tie(synthetic_run):
    run, tag = synthetic_run
    path = run / 'typing' / 'coverage' / f'{tag}.parquet'
    frame = pd.read_parquet(path)
    before, _ = ID.classify_image(run, tag)
    resolved = before['tags'].str.contains('tie_brightness_coverage')
    assert resolved.any()
    frame.loc[resolved, AC.INTERIOR_COUNT] = 0
    frame.loc[resolved, frame.columns.str.startswith(AC.INTERIOR_PREFIX)] = np.nan
    AC.write_coverage(run, tag, frame)
    after, audit = ID.classify_image(run, tag)
    assert not after.loc[resolved, 'tags'].str.contains('tie_brightness_coverage').any()
    assert (after.loc[resolved, AC.INTERIOR_COUNT] == 0).all()
    assert audit['cells_without_interior'] == int(resolved.sum())


def test_coverage_reorders_by_exact_ids(synthetic_run):
    run, tag = synthetic_run
    frame = pd.read_parquet(run / 'typing' / 'coverage' / f'{tag}.parquet')
    AC.write_coverage(run, tag, frame.iloc[::-1])
    out, _ = ID.classify_image(run, tag)
    assert out['label'].tolist() == [c[0] for c in CELLS]


@pytest.fixture()
def accuracy_tiles():
    path = TYPING.parents[1] / 'research' / 'typing_accuracy_tiles_20261029.py'
    spec = importlib.util.spec_from_file_location('accuracy_tiles_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fresh_review_is_seeded_nonoverlapping_and_unique(accuracy_tiles, tmp_path):
    candidates = pd.DataFrame({'image': 'SYN', 'object_id': [f'c{i}' for i in range(60)],
                               'group': 'compartment_tie', 'pair': 'a over b'})
    # A cell can qualify for both a compartment and a lineage review; show it only once.
    candidates = pd.concat([candidates, candidates.assign(group='lineage_tie')], ignore_index=True)
    previous = tmp_path / 'old.csv'
    candidates.iloc[:10].to_csv(previous, index=False)
    fresh = accuracy_tiles.exclude_reviewed(candidates, [previous])
    first = accuracy_tiles.sample(fresh, seed=19)
    repeat = accuracy_tiles.sample(fresh, seed=19)
    assert first.equals(repeat)
    assert not first.duplicated(['image', 'object_id']).any()
    assert not set(first['object_id']) & set(candidates.iloc[:10]['object_id'])
    assert set(first['object_id']) != set(accuracy_tiles.sample(fresh, seed=20)['object_id'])


def test_blinded_review_hides_calls_and_keeps_separate_key(accuracy_tiles, tmp_path, monkeypatch):
    text = []
    class CaptureDraw:
        def text(self, position, value, **kwargs):
            text.append(value)
    monkeypatch.setattr(accuracy_tiles.ImageDraw, 'Draw', lambda image: CaptureDraw())
    record = {k: 'secret-value' for k in ('group', 'pair', 'image', 'object_id', 'label', 'status', 'tags',
                                         'compartment', 'runner_up_compartment', 'lineage')}
    record.update({'review_id': 'R0001', 'blind_color_swap': True, 'green': 'A', 'magenta': 'B',
                   'panel_names': ['composite', 'A', 'B', 'DAPI'],
                   'tiles': [np.zeros((16, 16, 3), np.uint8) for _ in range(4)], 'markers': {}})
    accuracy_tiles.sheets([record], tmp_path, blind=True)
    reviewer = pd.read_csv(tmp_path / 'blinded' / 'manifest.csv')
    key = pd.read_csv(tmp_path / 'key' / 'manifest.csv')
    assert set(reviewer.columns) == {'sheet', 'row', 'review_id', 'green', 'magenta', 'verdict', 'notes'}
    assert key.loc[0, 'label'] == 'secret-value'
    assert all('secret-value' not in label and 'call:' not in label for label in text)
    assert (tmp_path / 'blinded' / 'review_01.png').exists()
