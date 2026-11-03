"""Rejection controls using the actual completed release evidence, never science fixtures."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('composite', ROOT / 'scripts/compose_release_parity.py')
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)
EVIDENCE = ROOT / 'tmp/pubapps-component-release-20261102-v2/evidence'


class ActualEvidenceRejectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = Path(os.environ.get('ISLET_PARITY_REVALIDATION_REPORT',
                        str(EVIDENCE / 'release-parity-e5689f6-revalidated-r3.json')))
        if not cls.path.exists():
            raise unittest.SkipTest('Completed real release revalidation evidence is unavailable')
        cls.report = MOD.read(cls.path)
        cls.proof = MOD.load_reference(cls.report['validation_provenance']['evidence']['unchanged_inputs'])
        cls.manifest = MOD.load_reference(cls.report['validation_provenance']['evidence']['release_manifest'])
        cls.original = MOD.load_reference(cls.report['validation_provenance']['evidence']['original_full_report'])

    def reject_report(self, report):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder) / 'mutated-report.json'
            p.write_text(json.dumps(report))
            with self.assertRaises((ValueError, KeyError)):
                MOD.validate_report(p)

    def test_actual_complete_evidence_passes(self):
        self.assertTrue(MOD.validate_report(self.path)['passed'])

    def test_original_second_failure_rejected(self):
        value = copy.deepcopy(self.original)
        value['checks'][0]['status'] = 'failed'
        with self.assertRaises(ValueError):
            MOD.check_original(value)

    def test_changed_retained_check_rejected(self):
        value = copy.deepcopy(self.report)
        value['checks'][0]['details']['tables'] += 1
        self.reject_report(value)

    def test_altered_runtime_identity_rejected(self):
        value = copy.deepcopy(self.report)
        value['identity']['git_commit'] = MOD.VALIDATOR
        self.reject_report(value)

    def test_altered_validator_identity_rejected(self):
        value = copy.deepcopy(self.report)
        value['validation_provenance']['validator_commit'] = MOD.PRODUCER
        self.reject_report(value)

    def test_missing_or_mismatched_evidence_rejected(self):
        for action in ('missing', 'hash'):
            value = copy.deepcopy(self.report)
            evidence = value['validation_provenance']['evidence']
            if action == 'missing':
                del evidence['independent_uuid_audit']
            else:
                evidence['independent_uuid_audit']['sha256'] = '0' * 64
            self.reject_report(value)

    def test_omitted_proof_file_rejected(self):
        proof = copy.deepcopy(self.proof)
        del proof['files'][next(iter(proof['files']))]
        with self.assertRaises(ValueError):
            MOD.validate_inputs(proof, self.manifest)

    def test_omitted_roster_rejected(self):
        proof = copy.deepcopy(self.proof)
        del proof['rosters'][next(iter(proof['rosters']))]
        with self.assertRaises(ValueError):
            MOD.validate_inputs(proof, self.manifest)

    def test_image_pin_change_rejected(self):
        proof = copy.deepcopy(self.proof)
        path = next(p for p in proof['files'] if p.endswith('.ome.tiff'))
        proof['files'][path]['sha256'] = '0' * 64
        with self.assertRaises(ValueError):
            MOD.validate_inputs(proof, self.manifest)

    def test_missing_source_link_rejected(self):
        proof = copy.deepcopy(self.proof)
        del proof['source_links'][next(iter(proof['source_links']))]
        with self.assertRaises(ValueError):
            MOD.validate_inputs(proof, self.manifest)

    def test_altered_revalidation_metadata_rejected(self):
        value = copy.deepcopy(self.report)
        value['validation_provenance']['method'] = 'all checks rerun'
        self.reject_report(value)

    def reject_stripped_details(self, mode):
        evidence = self.report['validation_provenance']['evidence']
        source = MOD.load_reference(evidence['source_manifest'])
        correction = MOD.load_reference(evidence['corrected_state_report'])
        uuid = MOD.load_reference(evidence['independent_uuid_audit'])
        if mode == 'coordinates':
            for row in uuid['donor_rows']:
                row['coordinate_comparison'] = {}
        elif mode == 'uuid_distance':
            del uuid['distance_measurement']
        else:
            correction['distance_comparison'] = {'max_absolute_difference': 0}
        with self.assertRaises((ValueError, KeyError)):
            MOD.check_correction(correction, uuid, self.original, source, self.manifest)

    def test_stripped_coordinate_evidence_rejected(self):
        self.reject_stripped_details('coordinates')

    def test_stripped_uuid_distance_evidence_rejected(self):
        self.reject_stripped_details('uuid_distance')

    def test_stripped_corrected_distance_evidence_rejected(self):
        self.reject_stripped_details('corrected_distance')


if __name__ == '__main__':
    unittest.main()
