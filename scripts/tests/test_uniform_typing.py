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
import tifffile

TYPING = Path(__file__).resolve().parents[1] / 'uniform' / 'typing'
sys.path.insert(0, str(TYPING))
import registry as RG  # noqa: E402
import identity as ID  # noqa: E402
import levels as LV  # noqa: E402
import area_coverage as AC  # noqa: E402
from interior_signal import compare_signals  # noqa: E402

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
        if r.kind == 'pattern':              # a marker-pattern cell type (owner 2026-11-06): named_by lists its markers
            assert r.named_by and all(m in reg.typing_markers() for m in RG.split(r.named_by)), r
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


# -- the canary status (Xenium notebook 02 cell 22, margin between groups; tie assignment owner 2026-11-18) ----------------
def test_group_status():
    det = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 0], [1, 1, 1], [1, 1, 0], [1, 1, 0], [1, 1, 0]], bool)
    depth = np.array([[0, 0, 0], [2, 0, 0], [3, 1, 0], [2, 2, 0], [3, 1, 1], [2, 2, 0], [2, 2, 0], [1, 1, 0]])
    # equal depth: more exclusive interior pixels assign the side; pixels beat brightness; depth-1 ties are assigned too
    bright = np.array([[0, 0, 0]] * 5 + [[1, 2, 0], [2, 1, 0], [1, 2, 0]], float)
    cover = np.array([[0, 0, 0]] * 5 + [[.2, .9, 0], [.2, .9, 0], [.2, .9, 0]], float)
    status, top, runner, td, rd, ng, settled, basis = ID.group_status(det, depth, np.zeros(det.shape), 1, 2,
                                                                      bright, cover)
    # r3: measured but every basis exactly tied -> the one remaining conflict; r4: three groups reduce to the top two
    # and the depth margin resolves; r6: pixels (side 1) beat brightness (side 0); r7: a depth-1 tie is assigned
    assert list(status) == ['none', 'single', 'resolved', 'lineage_ambiguous', 'resolved', 'resolved',
                            'resolved', 'resolved']
    assert list(ng) == [0, 1, 2, 2, 3, 2, 2, 2] and top[2] == 0 and runner[2] == 1
    assert top[4] == 0 and runner[4] in (1, 2)
    assert top[5] == 1 and top[6] == 1 and top[7] == 1
    assert list(settled) == [False] * 5 + [True, True, True]
    assert list(basis[5:]) == ['interior_pixels'] * 3 and basis[3] == ''


@pytest.mark.parametrize('cover_pair,bright_pair,score_pair,expected_top,expected_basis', [
    ((0.21, 0.20), (1., 2.), (0., 0.), 0, 'interior_pixels'),      # pixels beat brightness, no minimum lead
    ((0.20, 0.21), (2., 1.), (0., 0.), 1, 'interior_pixels'),
    ((0.20, 0.20), (1., 2.), (0., 0.), 1, 'interior_brightness'),  # exact pixel tie: the brighter median
    ((0.20, 0.20), (2., 2.), (3., 1.), 0, 'fold_ranking'),         # both tied: the fold ranking assigns
    ((np.nan, np.nan), (np.nan, np.nan), (3., 1.), 0, 'fold_ranking'),  # no interior evidence at all
])
def test_interior_tie_assignment(cover_pair, bright_pair, score_pair, expected_top, expected_basis):
    status, top, runner, *_, settled, basis = ID.group_status(
        np.ones((1, 2), bool), np.full((1, 2), 2), np.array([list(score_pair)]), 1, 2,
        np.array([list(bright_pair)]), np.array([list(cover_pair)]))
    assert status[0] == 'resolved' and bool(settled[0]) and top[0] == expected_top and basis[0] == expected_basis


def test_interior_tie_exact_tie_stays_conflict():
    status, top, runner, *_, settled, basis = ID.group_status(
        np.ones((1, 2), bool), np.full((1, 2), 2), np.zeros((1, 2)), 1, 2,
        np.array([[1., 1.]]), np.array([[0.2, 0.2]]))
    assert status[0] == 'lineage_ambiguous' and not settled[0] and basis[0] == ''


def test_interior_invalid_falls_back_to_fold():
    status, top, *_ , settled, basis = ID.group_status(
        np.ones((1, 2), bool), np.full((1, 2), 2), np.array([[1., 3.]]), 1, 2,
        np.array([[2., 1.]]), np.array([[0.9, 0.2]]), np.array([False]))
    assert status[0] == 'resolved' and top[0] == 1 and basis[0] == 'fold_ranking'


# -- identity on a synthetic image --------------------------------------------------------------------------------------
LEVEL_RAW, CANARY_RAW, HIGH = 100.0, 1000.0, 2000.0
CELLS = [  # (expected label, start class, {marker: raw value}, established endocrine call); every other marker 0
    ('CD8a Tcell', None, {'CD3e': HIGH, 'CD45': HIGH, 'CD8a': HIGH}, None),
    ('CD4 Tcell', None, {'CD3e': HIGH, 'CD45': HIGH, 'CD8a': HIGH, 'CD4': 2 * HIGH}, None),   # both: the higher (10-27; restored 11-06)
    ('Hematopoietic CD45+', None, {'CD3e': HIGH, 'CD45': HIGH}, None),                        # DN: regular workflow (11-06)
    ('Hematopoietic CD45-', None, {'CD20': HIGH}, None),                                      # CD20 alone: never a B cell
    ('B cell', None, {'CD20': HIGH, 'CD45': HIGH}, None),
    ('Macrophage_LGALS3-', None, {'CD68': HIGH, 'CD45': HIGH}, None),
    ('Hematopoietic CD45+', None, {'CD45': HIGH}, None),
    ('Hematopoietic CD45-', None, {'CD163': HIGH, 'CD8a': 2 * HIGH}, None),                  # tie assigned to T by interior pixels (11-18); T without CD45 is not named
    ('Acinar cell', None, {'ECAD': HIGH, 'panCK': HIGH}, None),                               # CK19- epithelial residual (11-06)
    ('Acinar cell', 'Islet', {'ECAD': HIGH}, None),                                           # exocrine inside the border (10-29)
    ('Beta cell', 'Islet', {'INS': HIGH, 'ECAD': HIGH}, None),                                # hormone, whole cell
    ('Beta cell', 'Islet', {'INS': HIGH, 'GCG': 0.8 * HIGH}, None),                           # the highest scaled hormone
    ('CHGA+ hormone-negative', 'EndocrineCluster', {'CHGA': HIGH}, None),
    ('Unclassified', 'EndocrineSingle', {'CD31': HIGH}, None),                                # unnamed in its object: explicit Unclassified (2026-11-05)
    ('Alpha cell', 'EndocrineSingle', {'GCG': HIGH, 'CD31': HIGH}, None),                     # a verified single stays
    ('Acinar cell', None, {'INS': HIGH, 'ECAD': HIGH}, None),                                 # no endocrine outside objects
    ('Beta cell', 'Islet', {'INS': HIGH}, None),                                              # hormone in an islet
    ('Hematopoietic CD45+', 'Islet', {'CD3e': HIGH, 'CD45': HIGH}, ('Immune', '')),           # DN in an islet: regular workflow
    ('Tissue', None, {}, None),                                                               # no marker above background
    ('Tissue', None, {'ColIV': 1e6}, None),                                                   # ColIV is never read
    ('Tissue', None, {'CD31': HIGH}, None),                                                   # no vessel: CD31 names nothing (10-30)
    ('Tissue', None, {'BCatenin': HIGH}, None),                                               # only a shared marker
    ('Endothelium Blood Vessel', None, {'CD31': 2 * HIGH, 'SMA': HIGH}, None),                # in a vessel object
    ('Stromal', None, {'CD31': HIGH, 'SMA': 2 * HIGH}, None),                                 # CD31 not on top: passes through
    ('Tissue', None, {'CD31': HIGH}, None),                                                   # under half its area inside
    ('Unclassified', 'Islet', {'INS': HIGH, 'GCG': HIGH}, None),                              # exact tie: no subtype, explicit Unclassified (2026-11-05)
    ('CD8a Tcell', None, {'CD3e': HIGH, 'CD45': HIGH, 'CD8a': HIGH, 'CD163': HIGH}, None),     # T leads Myeloid by one
    ('Hematopoietic CD45+', None, {'CD3e': 2 * HIGH, 'CD45': HIGH, 'CD163': HIGH}, None),     # tie won, DN: regular workflow
    ('Hematopoietic CD45+', None, {'CD3e': 2 * HIGH, 'CD45': HIGH, 'CD163': HIGH}, None),     # tie assigned to T by median brightness at a 0-0 interior-pixel tie (11-18); DN T takes the pan label
    ('Hematopoietic CD45-', None, {'CD4': HIGH, 'CD8a': HIGH}, None),                         # no lineage, CD45 negative
    ('Macrophage_LGALS3-', None, {'CD68': HIGH, 'CD45': HIGH, 'CD4': HIGH}, None),
    ('Macrophage_LGALS3+', None, {'CD163': HIGH, 'CD45': HIGH, 'LGALS3': HIGH}, None),
    ('Macrophage_LGALS3+', None, {'CD163': HIGH, 'CD45': HIGH, 'LGALS3': HIGH, 'CD4': HIGH}, None),
    ('Macrophage_LGALS3+', None, {'CD68': HIGH, 'CD163': HIGH, 'CD45': HIGH, 'LGALS3': HIGH}, None),
    ('Macrophage_LGALS3-', None, {'CD68': HIGH, 'CD163': HIGH, 'CD45': HIGH}, None),
    ('Hematopoietic CD45-', None, {'CD68': HIGH, 'LGALS3': HIGH}, None),
    ('Hematopoietic CD45+', None, {'CD45': HIGH, 'LGALS3': HIGH}, None),
    ('Hematopoietic CD45+', None, {'CD4': HIGH, 'CD45': HIGH}, None),
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
                       'level_log': np.log1p(LEVEL_RAW), 'level_raw': LEVEL_RAW})
    lv.to_csv(typ / 'levels.csv', index=False)
    # the reading competing groups are compared on (owner 2026-10-30): same values in every compartment here
    diff = [m for m in markers if reg.differentiation(m)]
    pd.DataFrame({'image': tag, 'marker': diff, 'measurement': [reg.differentiation(m) for m in diff],
                  'level_log': np.log1p(LEVEL_RAW)}).to_csv(typ / 'levels_differentiation.csv', index=False)
    pd.DataFrame({'marker': markers, 'decision': 'use'}).to_csv(typ / 'levels_selection.csv', index=False)
    lv.assign(canary_log=np.log1p(CANARY_RAW))[['image', 'marker', 'canary_log']].to_csv(
        typ / 'canary_levels.csv', index=False)
    # Cells occupy disjoint 9 by 6 rectangles: their eroded interiors contain 5 by 2 pixels.
    features = []
    channels = np.zeros((len(markers), 12, 12 * n), dtype=np.uint16)
    for i, (_, _, values, _) in enumerate(CELLS):
        x = 12 * i + 1
        coords = [[x,2],[x+8,2],[x+8,7],[x,7],[x,2]]
        features.append(dict(id=ids[i], geometry=dict(type='Polygon',coordinates=[coords])))
        for m, val in values.items():
            if m in markers:
                channels[markers.index(m), 2:8, x:x+9] = val
        if i in (27, 28):
            # Both markers clear the old cell-mean gate in every interior pixel.
            # In cell 27 T wins; cell 28 has an exact calibrated tie.
            channels[markers.index('CD3e'), 2:8, x:x+9] = 8000 if i == 27 else 0
            channels[markers.index('CD163'), 2:8, x:x+9] = 1000 if i == 27 else 0

    image = tmp_path / 'image.tiff'
    tifffile.imwrite(image, channels, photometric='minisblack', metadata={'axes':'CYX'})
    (tmp_path / 'project_manifest.json').write_text(json.dumps({'entries': [
        {'image_name': tag, 'entry_id': 1, 'image_path': str(image)}]}))
    server = tmp_path / 'qupath' / 'data' / '1'
    server.mkdir(parents=True)
    (server / 'server.json').write_text(json.dumps({'metadata': {'channels': [{'name': m} for m in markers]}}))
    (tmp_path / 'cell_geojson').mkdir()
    (tmp_path / 'cell_geojson' / f'{tag}.geojson').write_text(json.dumps(features))
    AC.write_coverage(tmp_path, tag, coverage)
    return tmp_path, tag


def test_exclusive_signal_counts_and_medians():
    # Counts from the three case-6533 interior masks. Both markers are positive everywhere.
    for n, sst in [(108, 86), (69, 55), (125, 123)]:
        left = np.ones((1, n))
        right = np.full((1, n), .5)
        right[:, :sst] = 2.
        assert compare_signals(left, right) == (n, n - sst, sst, 1., 2.)
    assert compare_signals([[1., 1., -1.]], [[1., 2., -2.]])[:3] == (3, 0, 1)
    assert compare_signals([[1., np.nan]], [[2., 0.]]) is None


def test_tied_hormones_reach_spatial_rule(synthetic_run):
    run, tag = synthetic_run
    # Keep the insulin/glucagon cell scores exactly tied and change only glucagon pixels in c25.
    image = run / 'image.tiff'
    data = tifffile.imread(image)
    markers = RG.load().measured_markers()
    x = 12 * 25 + 1
    data[markers.index('GCG'), 2:8, x:x+9] = 6000
    tifffile.imwrite(image, data, photometric='minisblack', metadata={'axes':'CYX'})
    frame = pd.read_parquet(run / 'typing/coverage' / f'{tag}.parquet')
    AC.write_coverage(run, tag, frame)
    out, _ = ID.classify_image(run, tag)
    row = out.set_index('object_id').loc['c25']
    assert row['label'] == 'Alpha cell' and 'hormone_tie_spatial_dominance' in row['tags']
    evidence = json.loads(row['spatial_tie_evidence'])['hormone']
    assert evidence['interior_pixels'] == 10
    assert sorted(evidence['winner_pixels']) == [0, 10]


def test_spatial_review_displays_the_actual_competing_groups(accuracy_tiles):
    comparison = {'markers': [['INS'], ['SST']], 'interior_pixels': 69,
                  'winner_pixels': [14, 55], 'median_scaled': [0.6, 1.1]}
    row = pd.Series({'group': 'hormone_tie', 'spatial_tie_evidence': json.dumps({'hormone': comparison})})
    assert accuracy_tiles.spatial_comparison(row) == comparison
    assert accuracy_tiles.shown_markers(row, {}, [], None) == (['SST'], ['INS'])


def test_identity_rules(synthetic_run):
    run, tag = synthetic_run
    out, audit = ID.classify_image(run, tag)
    got = out['label'].tolist()
    want = [c[0] for c in CELLS]
    assert got == want, [(i, g, w) for i, (g, w) in enumerate(zip(got, want)) if g != w]
    # owner 2026-11-06: a DN T-lineage cell is not named T cell; it takes the regular cascade (pan label)
    assert out.iloc[2]['label'] == 'Hematopoietic CD45+' and out.iloc[2]['subtype_status'] == 'T_dn_regular_workflow'
    assert out.iloc[11]['named_by'] == 'INS' and out.iloc[25]['status'] == 'tied_highest_hormones'
    assert out.iloc[5]['named_by'] == 'CD68;CD163'
    assert out.iloc[13]['status'] == 'no_endocrine_marker_above_level'
    assert out.iloc[10]['named_by'] == 'INS' and out.iloc[14]['named_by'] == 'GCG'
    assert 'endocrine_outside_object' in out.iloc[15]['tags']
    assert out.iloc[18]['status'] == out.iloc[19]['status'] == 'no_marker_above_background'
    assert not any(str(x).startswith('ColIV') for x in got)
    # no cell is named by a marker (owner 2026-10-28): the start class stays, with the status saying why
    assert out.iloc[3]['status'] == 'no_lineage_named' and out.iloc[21]['status'] == 'shared_marker_only'
    assert out.iloc[22]['named_by'] == 'CD31' and 'passed_through_Blood_vessel' in out.iloc[23]['tags']
    assert not any(str(x)[-1:] in '+-' and ' ' not in str(x) and 'Macrophage_' not in str(x) for x in got)


def test_macrophage_with_missing_gal3_retains_broad_identity(synthetic_run):
    run, tag = synthetic_run
    path = run / 'typing' / 'redsea' / f'{tag}.parquet'
    cells = pd.read_parquet(path)
    cells.loc[5, 'Cell__LGALS3'] = np.nan
    cells.to_parquet(path)
    result, _ = ID.classify_image(run, tag)
    assert result.iloc[5].label == 'Macrophage'
    assert result.iloc[5].lineage == 'Myeloid'
    assert result.iloc[5].named_by == 'CD68;CD163'


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


def test_no_interior_is_explicit_and_falls_back_to_fold(synthetic_run):
    run, tag = synthetic_run
    path = run / 'typing' / 'coverage' / f'{tag}.parquet'
    frame = pd.read_parquet(path)
    before, _ = ID.classify_image(run, tag)
    resolved = before['tags'].str.contains('tie_interior_')
    assert resolved.any()
    frame.loc[resolved, AC.INTERIOR_COUNT] = 0
    frame.loc[resolved, frame.columns.str.startswith(AC.INTERIOR_PREFIX)] = np.nan
    AC.write_coverage(run, tag, frame)
    after, audit = ID.classify_image(run, tag)
    # without interior pixels a tie is still assigned by the fold ranking; an exact fold tie (cell 30: CD68 and
    # CD4 both at HIGH) is the one case left as a conflict, and the absence stays audited
    assert not after.loc[resolved, 'tags'].str.contains('tie_interior_').any()
    fold = after.loc[resolved, 'tags'].str.contains('tie_fold_ranking')
    conflict = after.loc[resolved, 'label'].str.startswith('Conflict ')
    assert (fold | conflict).all() and fold.sum() == 3 and conflict.sum() == 1
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


def test_review_revision_matches_frozen_replay(synthetic_run):
    from review_revision import apply_decisions, reconsider
    run, tag = synthetic_run
    prior, _ = ID.classify_image(run, tag)
    candidate, audit = ID.classify_image(run, tag, review_revision=True)
    coverage = pd.read_parquet(run / 'typing' / 'coverage' / f'{tag}.parquet')
    context = pd.read_parquet(run / 'typing' / 'context' / f'{tag}.parquet')
    raw = pd.read_parquet(run / 'typing' / 'redsea' / f'{tag}.parquet').set_index('object_id')
    decisions = reconsider(prior, coverage, context, np.isfinite(raw.Cell__LGALS3))
    pd.testing.assert_frame_equal(candidate, apply_decisions(prior, decisions))
    assert decisions.disposition.value_counts().to_dict() == audit['review_revision']
