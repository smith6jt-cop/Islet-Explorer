"""A withdrawn candidate cannot be reused for a fresh native delivery."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cohort_expansion'))
from prepare_endocrine_candidates import refuse_withdrawn_candidate
from package_endocrine_review import package


class QualityWithdrawalTests(unittest.TestCase):
    def test_fresh_project_cannot_bypass_candidate_withdrawal(self):
        with tempfile.TemporaryDirectory() as tmp:
            candidate=Path(tmp)/'candidate'; candidate.mkdir()
            project=Path(tmp)/'fresh_project'; project.mkdir()
            (candidate/'quality_status.json').write_text(json.dumps({'status':'withdrawn_pending_boundary_correction'}))
            with self.assertRaisesRegex(ValueError,'withdrawal blocks delivery'):
                package(candidate,project)
            self.assertEqual(list(project.iterdir()),[])

    def test_project_withdrawal_also_blocks_delivery(self):
        with tempfile.TemporaryDirectory() as tmp:
            candidate=Path(tmp)/'candidate'; candidate.mkdir()
            project=Path(tmp)/'project'; project.mkdir()
            (project/'quality_status.json').write_text(json.dumps({'status':'failed_boundary_review'}))
            with self.assertRaisesRegex(ValueError,'withdrawal blocks delivery'):
                package(candidate,project)

    def test_new_candidate_is_not_scientifically_accepted_by_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            refuse_withdrawn_candidate(Path(tmp))
            self.assertEqual(list(Path(tmp).iterdir()),[])


if __name__ == '__main__': unittest.main()
