#!/usr/bin/env python3
"""Focused regression tests for the single CRAN VIM/editor collision."""

import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "sbom_correction", Path(__file__).resolve().parents[1] / "correct-syft-vim-cpe.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def fixture():
    return {"bomFormat": "CycloneDX", "serialNumber": "preserved",
            "metadata": {"tools": {"components": [{"name": "syft", "version": "1.51.1"}]}},
            "components": [{"type": "library", "name": "VIM", "version": "6.2.6",
                            "purl": "pkg:cran/VIM@6.2.6", "cpe": MODULE.INCORRECT_CPE,
                            "bom-ref": "preserved-reference",
                            "licenses": [{"license": {"id": "GPL-2.0-or-later"}}],
                            "properties": [{"name": k, "value": v}
                                           for k, v in MODULE.REQUIRED_PROPERTIES.items()]}],
            "dependencies": [{"ref": "preserved-reference", "dependsOn": []}]}


class CorrectionTests(unittest.TestCase):
    def test_exact_collision_changes_only_cpe_and_preserves_original(self):
        original = fixture()
        before = copy.deepcopy(original)
        actual, changes = MODULE.correct_document(original)
        expected = copy.deepcopy(original)
        del expected["components"][0]["cpe"]
        self.assertEqual(actual, expected)
        self.assertEqual(original, before)
        self.assertEqual(changes, [{"json_pointer": "/components/0/cpe",
                                   "removed_value": MODULE.INCORRECT_CPE,
                                   "component_purl": "pkg:cran/VIM@6.2.6"}])

    def test_genuine_editor_and_unrelated_components_are_unchanged(self):
        original = fixture()
        original["components"] += [
            {"type": "application", "name": "vim", "version": "6.2.6",
             "purl": "pkg:deb/ubuntu/vim@6.2.6", "cpe": MODULE.INCORRECT_CPE},
            {"type": "library", "name": "other", "version": "1", "cpe": "other-cpe"},
        ]
        actual, changes = MODULE.correct_document(original)
        self.assertEqual(actual["components"][1:], original["components"][1:])
        self.assertEqual(len(changes), 1)

    def test_any_component_identity_change_is_unchanged(self):
        for key, value in {"type": "application", "name": "vim", "version": "6.2.7",
                           "purl": "pkg:cran/VIM@6.2.7", "cpe": "another-cpe"}.items():
            with self.subTest(field=key):
                original = fixture()
                original["components"][0][key] = value
                self.assertEqual(MODULE.correct_document(original), (original, []))

    def test_missing_changed_or_duplicate_required_properties_are_unchanged(self):
        for index in range(len(MODULE.REQUIRED_PROPERTIES)):
            for mode in ("missing", "changed", "duplicate"):
                with self.subTest(property=index, mode=mode):
                    original = fixture()
                    properties = original["components"][0]["properties"]
                    if mode == "missing":
                        properties.pop(index)
                    elif mode == "changed":
                        properties[index]["value"] = "different"
                    else:
                        properties.append(copy.deepcopy(properties[index]))
                    self.assertEqual(MODULE.correct_document(original), (original, []))

    def test_changed_scanner_version_is_unchanged(self):
        original = fixture()
        original["metadata"]["tools"]["components"][0]["version"] = "1.51.2"
        self.assertEqual(MODULE.correct_document(original), (original, []))

    def test_correction_is_idempotent(self):
        corrected, _ = MODULE.correct_document(fixture())
        self.assertEqual(MODULE.correct_document(corrected), (corrected, []))

    def test_original_bytes_and_hash_bound_record(self):
        with tempfile.TemporaryDirectory() as temp:
            source, target, record = [Path(temp) / name for name in ("raw.json", "scan.json", "record.json")]
            raw = json.dumps(fixture()).encode()
            source.write_bytes(raw)
            result = MODULE.write_correction(source, target, record)
            self.assertEqual(source.read_bytes(), raw)
            self.assertEqual(result["original_sbom_sha256"], hashlib.sha256(raw).hexdigest())
            self.assertEqual(result["scan_sbom_sha256"], hashlib.sha256(target.read_bytes()).hexdigest())
            self.assertEqual(json.loads(record.read_text()), result)
            with self.assertRaises(ValueError):
                MODULE.write_correction(source, target, record)

    def test_noop_preserves_original_formatting(self):
        with tempfile.TemporaryDirectory() as temp:
            source, target, record = [Path(temp) / name for name in ("raw.json", "scan.json", "record.json")]
            raw = b'{"bomFormat":"CycloneDX", "components":[]}\n'
            source.write_bytes(raw)
            result = MODULE.write_correction(source, target, record)
            self.assertEqual(target.read_bytes(), raw)
            self.assertEqual(result["correction_count"], 0)

    def test_invalid_document_fails(self):
        for document in ({}, {"bomFormat": "CycloneDX", "components": {}}):
            with self.assertRaises(ValueError):
                MODULE.correct_document(document)


if __name__ == "__main__":
    unittest.main()
