"""Read-only regression checks against retained real classifier/geometry inputs.

No fabricated cells, masks, positivity values or scientific decisions are used.
Tests skip when this checkout does not contain the retained research evidence.
Temporary malformed metadata copies exercise rejection paths only.
"""
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/uniform/typing'))
import endocrine_object_audit as audit

BASE = ROOT / 'data/development/uniform-20donor-20261027'
RUN = ROOT / 'data/development/uniform-20donor-20261029-v19-qc-work'
REPLAY = ROOT / 'output/v20_classifier_support_20261102_r2/summary.json'
PROVENANCE = ROOT / 'output/v20_classifier_trace_20261102/classifier-provenance.json'
TAG = 'Aab_6538'
PINS = dict(classifier_run=BASE, classifier_replay=REPLAY,
            classifier_replay_sha256='476fc65e42bbd521a12f67e6ecb2c875726d3e8b8c0e0ab4f4cbc6db098ce877',
            geometry_reference=BASE, geometry_provenance=PROVENANCE,
            geometry_provenance_sha256='597603c1db91bc92effc3f8ece6f745b7f80d4fb7ec7384e5b9ae6761bde6646')


@pytest.fixture(scope='module', autouse=True)
def retained_inputs():
    required = [REPLAY, PROVENANCE, BASE / 'fit' / f'{TAG}.mask.npy',
                RUN / 'typing/redsea' / f'{TAG}.parquet', RUN / 'cell_geojson' / f'{TAG}.geojson']
    if not all(path.is_file() for path in required):
        pytest.skip('retained real classifier and cell evidence is unavailable')


@pytest.fixture(scope='module')
def reports(tmp_path_factory):
    root = tmp_path_factory.mktemp('real_endocrine_object_audit')
    strict = audit.audit(RUN, root / 'strict', images=[TAG])
    lineage = audit.audit(RUN, root / 'lineage', images=[TAG], **PINS)
    return root, strict, lineage


def test_preserved_corrected_cell_diagnostics_and_default_failure(reports):
    root, strict, lineage = reports
    assert strict['audit_mode'] == 'corrected_cell'
    assert strict['passes'] is False
    assert lineage['passes'] is True
    assert lineage['corrected_cell_diagnostic_passes'] is False
    for key in ('totals', 'by_class', 'by_image_and_class', 'checks'):
        assert strict[key] == lineage[key]
    old = pd.read_csv(root / 'strict/objects.csv')
    new = pd.read_csv(root / 'lineage/objects.csv')
    pd.testing.assert_frame_equal(old, new.loc[:, old.columns])
    assert lineage['totals']['objects_with_all_analysis_cells_negative'] == 10
    assert lineage['totals']['single_membership_anomalies'] == 1
    with pytest.raises(FileExistsError):
        audit.audit(RUN, root / 'lineage', images=[TAG], **PINS)


def test_real_disconnected_single_is_associated_without_rewriting_centroid(reports):
    root, _, report = reports
    objects = pd.read_csv(root / 'lineage/objects.csv').set_index('endocrine_object_id')
    row = objects.loc['a4281a63-be08-4a0a-91bd-efb5ed4407f6']
    assert row.single_associated_cell_id == 'd3a3179d-7791-4e76-a4a2-78fbc4975528'
    assert row.single_association_match_count == 1
    assert row.context_cells == 0
    assert row.single_membership_anomaly
    assert not row.single_association_matches_centroid_membership
    assert row.classifier_lineage_passes
    proof = report['classifier_lineage']['images'][0]
    assert proof['failures'] == []
    assert proof['objects'] == proof['objects_with_classifier_support'] == 234
    assert proof['singles'] == proof['singles_with_exact_cell_association'] == 128
    assert proof['singles_with_centroid_membership_difference'] == 1


def test_real_split_descendants_retain_ambiguous_source_sets(reports):
    root, _, report = reports
    objects = pd.read_csv(root / 'lineage/objects.csv').set_index('endocrine_object_id')
    for key in ('c7ba9741-e2e6-4772-8a68-f7079148fbc5', 'd3c038cb-30eb-4958-a421-b0afde92876d'):
        row = objects.loc[key]
        assert row.original_source_relation == 'split_descendant_spatial_support'
        assert not row.original_parent_unique
        parents = json.loads(row.original_parent_candidates)
        assert len(parents) == 2
        assert parents == json.loads(row.original_supported_parent_candidates)
        assert row.original_overlap_area_pixels > 0
        assert row.classifier_lineage_passes
    assert report['classifier_lineage']['totals']['split_descendants_with_multiple_parent_candidates'] == 2


def test_actual_receipt_pins_and_source_parameter_bytes_are_required(tmp_path):
    evidence = audit.load_lineage_evidence(**PINS, tags=[TAG])
    assert len(evidence['replay']) == len(evidence['geometry']) == 20
    wrong_pin = dict(PINS, classifier_replay_sha256='0' * 64)
    with pytest.raises(ValueError, match='SHA256 pin'):
        audit.load_lineage_evidence(**wrong_pin, tags=[TAG])
    incomplete = json.loads(REPLAY.read_text())
    incomplete['images'] = incomplete['images'][:-1]
    path = tmp_path / 'incomplete-replay.json'
    path.write_text(json.dumps(incomplete))
    with pytest.raises(ValueError, match='complete 20-image cohort'):
        audit.load_lineage_evidence(**dict(PINS, classifier_replay=path,
                                          classifier_replay_sha256=audit.sha256(path)), tags=[TAG])
    original = BASE / 'fit' / f'{TAG}.json'
    records = evidence['replay'][TAG]['source_sha256']
    assert audit.verify_recorded_source(original, f'fit/{TAG}.json', records, {}) == audit.sha256(original)
    changed = tmp_path / original.name
    changed.write_bytes(original.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='input changed'):
        audit.verify_recorded_source(changed, f'fit/{TAG}.json', records, {})
    with pytest.raises(ValueError, match='all six'):
        audit.audit(RUN, tmp_path / 'partial', images=[TAG], classifier_run=BASE)
    assert not (tmp_path / 'partial').exists()


def test_streaming_reader_matches_actual_export_and_rejects_truncation(tmp_path):
    path = BASE / 'script1_outlines' / f'{TAG}.geojson'
    assert list(audit.streamed_cell_features(path)) == audit.features(path)
    truncated = tmp_path / 'truncated.geojson'
    truncated.write_bytes(path.read_bytes()[:-3])
    with pytest.raises(ValueError, match='truncated|trailer'):
        list(audit.streamed_cell_features(truncated))
