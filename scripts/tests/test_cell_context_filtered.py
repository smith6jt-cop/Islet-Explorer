"""A filtered candidate must use certified current masks, never legacy fallback provenance."""
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1] / 'uniform'
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'typing'))
import new_detection_filter as NF
import cell_context as C


def setup_run(tmp_path):
    (tmp_path / 'filter').mkdir()
    (tmp_path / 'filter/generation.json').write_text('{}')
    (tmp_path / 'structures').mkdir()
    return [{'case_id': 'synthetic', 'image_name': 'A'}]


def test_current_certificate_keeps_reuse_ancestry_without_legacy_access(tmp_path, monkeypatch):
    rows = setup_run(tmp_path)
    reuse = {'certificate': '/cert/current.json', 'certificate_sha256': 'proof',
             'outputs': {cls: {'sha256': cls+'-mask', 'source': {'path': '/raw/'+cls}}
                         for cls in ('Nerve', 'Capillary', 'Lymphatic')}}
    monkeypatch.setattr(NF, 'structure_reuse_current', lambda *a: reuse)
    monkeypatch.setattr(C, 'STRUCTURAL', tmp_path / 'must-not-read-legacy')
    result = C.import_structures(tmp_path, rows)
    assert result['A']['Nerve']['source'] == '/cert/current.json'
    assert result['A']['Nerve']['raw_source'] == {'path': '/raw/Nerve'}
    assert result['A']['Nerve']['reuse_certificate_sha256'] == 'proof'


def test_missing_certificate_refuses_baseline_fallback(tmp_path, monkeypatch):
    rows = setup_run(tmp_path)
    monkeypatch.setattr(NF, 'structure_reuse_current', lambda *a: None)
    with pytest.raises(ValueError, match='lacks a current structure-reuse certificate'):
        C.import_structures(tmp_path, rows)
    assert not (tmp_path / 'structures/provenance.json').exists()
