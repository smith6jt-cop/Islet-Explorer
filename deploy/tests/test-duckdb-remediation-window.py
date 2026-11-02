#!/usr/bin/env python3
"""Exercise the narrow, dated release decision independently of the live clock."""
import copy
import datetime
import importlib.util
import json
import pathlib
import unittest

deploy = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("decision", deploy / "evaluate-brotli-report.py")
decision = importlib.util.module_from_spec(spec)
spec.loader.exec_module(decision)


class RemediationWindowTests(unittest.TestCase):
    def setUp(self):
        self.today = datetime.date(2026, 11, 2)
        self.window = json.loads((deploy / "duckdb-remediation-window.json").read_text())
        self.report = {
            "schema": "islet.brotli_verification/2", "cve": "CVE-2026-111497",
            "result": "fail", "errors": [decision.DUCKDB_ERROR],
            "duckdb": {"version": "1.5.2", "result": 0, "error_code": -9,
                       "expected_error_code": -17, "library_sha256": "a" * 64},
        }

    def test_known_finding_is_visible_and_nonblocking_through_due_date(self):
        original = copy.deepcopy(self.report)
        for day in [self.today, datetime.date(2026, 11, 10)]:
            self.assertEqual(decision.evaluate(self.report, self.window, day), decision.DEFERRED_RESULT)
        self.assertEqual(self.report, original)

    def test_expired_or_not_yet_authorized_window_rejects(self):
        for day in [datetime.date(2026, 11, 1), datetime.date(2026, 11, 11)]:
            with self.assertRaises(ValueError):
                decision.evaluate(self.report, self.window, day)

    def test_missing_authorization_rejects(self):
        with self.assertRaises(ValueError):
            decision.evaluate(self.report, today=self.today)

    def test_broader_scope_or_later_deadline_rejects(self):
        for key, value in [("component", "all Brotli"), ("patch_by", "2026-12-31"),
                           ("cve", "another-CVE"), ("institutional_approval_claimed", True)]:
            window = dict(self.window, **{key: value})
            with self.subTest(key=key), self.assertRaises(ValueError):
                decision.evaluate(self.report, window, self.today)

    def test_other_errors_including_pillow_remain_blocking(self):
        for errors in [[], ["timeout"], [decision.DUCKDB_ERROR, "vendored Pillow Brotli"]]:
            report = dict(self.report, errors=errors)
            with self.subTest(errors=errors), self.assertRaises(ValueError):
                decision.evaluate(report, self.window, self.today)

    def test_unknown_or_partial_binary_evidence_rejects(self):
        for key, value in [("version", "1.6.0"), ("result", 3), ("error_code", -1),
                           ("expected_error_code", -9), ("library_sha256", "missing")]:
            report = copy.deepcopy(self.report)
            report["duckdb"][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                decision.evaluate(report, self.window, self.today)

    def test_fixed_binary_passes_after_deferral_expires(self):
        self.report.update(result="pass", errors=[])
        self.report["duckdb"]["error_code"] = -17
        self.assertEqual(decision.evaluate(self.report, today=datetime.date(2026, 12, 1)),
                         "pass-no-high-or-critical")


if __name__ == "__main__":
    unittest.main()
