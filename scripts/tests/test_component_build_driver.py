"""Engineering contracts for the sealed component build; no simulated biology."""
import importlib.util
import json
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    'component_build_driver', Path(__file__).resolve().parents[1] / 'hpg/build_component_candidate.py')
DRIVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DRIVER)


@pytest.fixture
def config(tmp_path):
    value = {'expected_commit': 'b' * 40, 'selection_reason': 'Explicit selected component baseline'}
    for name in DRIVER.INPUTS:
        path = tmp_path / 'sources' / name
        path.parent.mkdir(exist_ok=True)
        if name in DRIVER.PINNED_REPORTS or name == 'donor_metadata':
            path.write_text('{}\n')
        else:
            path.mkdir()
        value[name] = str(path)
    for name in DRIVER.OUTPUTS:
        value[name] = str(tmp_path / 'fresh' / name)
    for name in DRIVER.PINNED_REPORTS:
        value[name + '_sha256'] = DRIVER.sha256(Path(value[name]))
    return value


def test_explicit_valid_config_is_read_only(config):
    paths = DRIVER.validate_config(config)
    assert len(paths) == len(DRIVER.INPUTS + DRIVER.OUTPUTS)
    assert all(not paths[name].exists() for name in DRIVER.OUTPUTS)


@pytest.mark.parametrize('key', DRIVER.PINNED_REPORTS)
def test_tampered_lineage_receipt_fails_before_outputs(config, key):
    Path(config[key]).write_text('{"changed":true}\n')
    with pytest.raises(ValueError, match='receipt SHA256 differs'):
        DRIVER.validate_config(config)
    assert all(not Path(config[name]).exists() for name in DRIVER.OUTPUTS)


@pytest.mark.parametrize('key', ('classifier_run', 'classifier_replay', 'classifier_replay_sha256',
                                'geometry_reference', 'geometry_provenance', 'geometry_provenance_sha256'))
def test_all_lineage_inputs_are_required(config, key):
    del config[key]
    with pytest.raises(ValueError, match='Config must explicitly name'):
        DRIVER.validate_config(config)


@pytest.mark.parametrize('placement', ('same_source', 'within_source', 'above_source', 'within_output'))
def test_output_cannot_overlap_input_or_output(config, placement):
    source = Path(config['classifier_run'])
    config['candidate_run'] = str({
        'same_source': source, 'within_source': source / 'new-run',
        'above_source': source.parent, 'within_output': Path(config['report_root']) / 'run',
    }[placement])
    with pytest.raises(ValueError, match='already exists|must be separate'):
        DRIVER.validate_config(config)


def test_dangling_output_link_is_not_fresh(config):
    path = Path(config['candidate_run'])
    path.parent.mkdir()
    path.symlink_to(path.parent / 'missing')
    with pytest.raises(ValueError, match='already exists'):
        DRIVER.validate_config(config)


@pytest.fixture
def gate(tmp_path):
    path = tmp_path / 'summary.json'
    report = {'schema': 'uniform.endocrine_object_audit/2', 'passes': True,
              'audit_mode': 'classifier_lineage', 'run': '/explicit/candidate',
              'images': [f'image-{index:02}' for index in range(20)]}
    return path, report


def check(path):
    return DRIVER.checked_gate(path, schema='uniform.endocrine_object_audit/2', pass_key='passes',
                               expected={'audit_mode': 'classifier_lineage', 'run': '/explicit/candidate'})


def test_gate_records_exact_report_hash(gate):
    path, report = gate
    path.write_text(json.dumps(report))
    assert check(path)['sha256'] == DRIVER.sha256(path)


@pytest.mark.parametrize('field,value', [
    ('passes', False), ('passes', 1), ('schema', 'uniform.endocrine_object_audit/1'),
    ('audit_mode', 'corrected_cell'), ('run', '/old/candidate'),
    ('images', ['one'] * 20), ('images', ['one']),
])
def test_failed_stale_or_partial_report_cannot_pass(gate, field, value):
    path, report = gate
    report[field] = value
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        check(path)


def test_missing_report_cannot_pass(gate):
    path, _ = gate
    with pytest.raises(FileNotFoundError):
        check(path)
