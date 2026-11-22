#!/usr/bin/env python3
"""Base-lock parsing and independent runtime identity regression tests."""
import pathlib
import runpy
import subprocess
import sys
import tempfile
import unittest

HELPER = pathlib.Path(__file__).resolve().parents[1] / "runtime-inputs.py"
MODULE = runpy.run_path(str(HELPER))
base_image_requirement = MODULE["base_image_requirement"]
runtime_fingerprint = MODULE["runtime_fingerprint"]
BASE = "docker.io/example/runtime@sha256:" + "a" * 64
OTHER = "docker.io/example/runtime@sha256:" + "b" * 64


class BaseImageRequirementTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = pathlib.Path(self.temp.name)
        (self.root / "deploy").mkdir()
        self.lock = self.root / "deploy/base-image.lock"
        self.lock.write_text(BASE + "\n")

    def test_normalizes_blank_lines_comments_and_outer_whitespace(self):
        self.lock.write_text("\n  # selected base\n\t" + BASE + "  \n\n# notes\n")
        self.assertEqual(base_image_requirement(self.root), BASE)

    def test_requires_exactly_one_digest_pinned_reference(self):
        for contents in ("", "# comment only\n", BASE + "\n" + OTHER,
                         BASE + "\n" + BASE, "example/runtime:latest",
                         "sha256:" + "a" * 64, BASE + " # trailing comment",
                         BASE[:-1], "example/runtime@sha512:" + "a" * 64,
                         "example/runtime@sha256:" + "A" * 64):
            with self.subTest(contents=contents):
                self.lock.write_text(contents)
                with self.assertRaises(ValueError):
                    base_image_requirement(self.root)

    def test_rejects_symlink_and_missing_lock(self):
        self.lock.unlink()
        with self.assertRaises(ValueError):
            base_image_requirement(self.root)
        target = self.root / "elsewhere.lock"
        target.write_text(BASE)
        self.lock.symlink_to(target)
        with self.assertRaises(ValueError):
            base_image_requirement(self.root)

    def test_git_selection_ignores_dirty_lock(self):
        subprocess.run(["git", "init", "-q", str(self.root)], check=True)
        subprocess.run(["git", "-C", str(self.root), "add", "."], check=True)
        subprocess.run(["git", "-C", str(self.root), "-c", "user.name=Test",
                        "-c", "user.email=test@example.invalid", "commit", "-qm", "base"], check=True)
        self.lock.write_text(OTHER)
        self.assertEqual(base_image_requirement(self.root, "HEAD"), BASE)
        self.assertEqual(base_image_requirement(self.root), OTHER)
        actual = subprocess.check_output([sys.executable, str(HELPER), "--base-image",
                                          str(self.root), "HEAD"], text=True).strip()
        self.assertEqual(actual, BASE)

    def test_base_identity_is_separate_from_unchanged_v1_fingerprint(self):
        (self.root / "deploy/Containerfile").write_text("ARG BASE_IMAGE\nFROM ${BASE_IMAGE}\n")
        (self.root / "deploy/containerignore").write_text("**\n")
        original = runtime_fingerprint(self.root)
        self.lock.write_text(OTHER)
        self.assertEqual(runtime_fingerprint(self.root), original)
        self.assertEqual(base_image_requirement(self.root), OTHER)
        # Existing CLI remains the fingerprint operation, not the base selector.
        actual = subprocess.check_output([sys.executable, str(HELPER), str(self.root)], text=True).strip()
        self.assertEqual(actual, original)


if __name__ == "__main__":
    unittest.main()
