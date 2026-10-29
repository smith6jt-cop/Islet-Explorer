"""Tests for the final typing chain (final typing plan of 2026-10-28, Step 6): the registry, the canary status logic, and
identity.py's rules on a small synthetic image built with the run's own registry, plus a check that no module of the
chain names a marker in code (owner 2026-10-28: marker roles are data, panel- and tissue-specific)."""
from __future__ import annotations

import ast
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
    status, top, runner, td, rd, ng, settled = ID.group_status(det, depth, np.zeros(det.shape), 1, 2, bright, cover)
    assert list(status) == ['none', 'single', 'resolved', 'lineage_ambiguous', 'doublet_suspected', 'resolved',
                            'lineage_ambiguous', 'lineage_ambiguous']   # brighter and covers more / split / under depth 2
    assert list(ng) == [0, 1, 2, 2, 3, 2, 2, 2] and top[2] == 0 and runner[2] == 1
    assert top[5] == 1 and runner[5] == 0 and list(settled) == [False] * 5 + [True, False, False]


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
    pd.DataFrame({'object_id': ids, **{ID.COV_PREFIX + m: [COVERAGE.get(i, {}).get(m, 0.8 if c[2].get(m, 0.0) > LEVEL_RAW
                                                                                    else 0.0)
                                                           for i, c in enumerate(CELLS)] for m in reg.table.index}}
                 ).to_parquet(typ / 'coverage' / f'{tag}.parquet')
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
