"""A failed original-donor comparison must precede preview side effects."""
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'research'))
import serve_preview


def test_failed_parity_does_not_stop_or_create_preview(tmp_path, monkeypatch):
    candidate = tmp_path / 'candidate'
    candidate.mkdir()
    (candidate / 'candidate.json').write_text(json.dumps({'schema_version': '4'}))
    app = tmp_path / 'app'
    app.mkdir()
    baseline = tmp_path / 'baseline'
    baseline.mkdir()
    called = []
    monkeypatch.setattr(serve_preview, 'audit_original_parity',
                        lambda *args: {'passed': False})
    monkeypatch.setattr(serve_preview, 'stop_recorded_preview',
                        lambda *args: called.append('stopped'))
    monkeypatch.setattr(sys, 'argv', ['serve_preview.py',
        '--data-root', str(candidate), '--app-root', str(app),
        '--baseline-root', str(baseline), '--python', sys.executable,
        '--restart', '--replace-preview', str(candidate)])
    before = list(candidate.iterdir())
    with pytest.raises(ValueError, match='Original-donor parity failed'):
        serve_preview.main()
    assert not called
    assert list(candidate.iterdir()) == before


def test_guard_recomputes_comparison_instead_of_trusting_saved_report(monkeypatch):
    calls = []
    def compare(baseline, candidate):
        calls.append((baseline, candidate))
        return {'passed': len(calls) == 1}
    monkeypatch.setattr(serve_preview, 'audit_original_parity', compare)
    assert serve_preview.require_original_parity('original', 'expanded')['passed']
    with pytest.raises(ValueError, match='Original-donor parity failed'):
        serve_preview.require_original_parity('original', 'expanded')
    assert calls == [('original', 'expanded'), ('original', 'expanded')]
