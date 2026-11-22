#!/usr/bin/env python3
"""Local cache integrity and no-repeat-read regression tests."""
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location(
    "prepare_local_data", Path(__file__).resolve().parents[1] / "prepare-local-data.py"
)
cache = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cache)


class LocalDataCacheTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "source"
        self.target = self.root / "cache" / "dataset-one"
        entries = []
        for tree in cache.TREES:
            path = self.source / "data" / tree / "table.bin"
            path.parent.mkdir(parents=True)
            path.write_bytes(tree.encode())
            entries.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  data/{tree}/table.bin\n")
        self.image = self.source / "data/images/test.ome.tiff"
        self.image.parent.mkdir()
        self.image.write_bytes(b"microscopy bytes must never be opened by cache preparation")
        entries.append(f"{hashlib.sha256(self.image.read_bytes()).hexdigest()}  data/images/test.ome.tiff\n")
        (self.source / "MANIFEST.sha256").write_text("".join(entries))

    def prepare(self):
        with contextlib.redirect_stdout(io.StringIO()):
            cache.prepare(self.source, self.target)

    def test_initial_copy_checks_hashes_without_opening_microscopy(self):
        original = Path.open

        def guarded_open(path, *args, **kwargs):
            if path == self.image:
                raise AssertionError("cache preparation opened TIFF")
            return original(path, *args, **kwargs)

        with patch.object(Path, "open", guarded_open):
            self.prepare()
        self.assertTrue((self.target / ".verified.json").is_file())
        self.assertFalse((self.target / "data/images").exists())

    def test_hashing_does_not_require_python_311(self):
        path = self.root / "multiple-blocks.bin"
        content = bytes(range(256)) * 8193
        path.write_bytes(content)
        with patch.object(cache.hashlib, "file_digest", create=True,
                          side_effect=AssertionError("Python 3.11 API used")):
            self.assertEqual(cache.digest(path), hashlib.sha256(content).hexdigest())
            self.prepare()

    def test_reuse_only_hashes_small_manifest_and_does_not_copy(self):
        self.prepare()
        original = cache.digest

        def guarded_digest(path):
            self.assertEqual(path, self.source / "MANIFEST.sha256")
            return original(path)

        with patch.object(cache, "digest", guarded_digest), patch.object(
            cache.subprocess, "run", side_effect=AssertionError("reused cache copied again")
        ):
            self.prepare()

    def test_changed_cache_rejected_without_repair(self):
        self.prepare()
        path = self.target / "data/spatial/table.bin"
        path.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "Published local cache changed"):
            self.prepare()
        self.assertEqual(path.read_bytes(), b"changed")

    def test_extra_file_rejected(self):
        self.prepare()
        (self.target / "data/cells/extra.bin").write_bytes(b"extra")
        with self.assertRaisesRegex(ValueError, "Published local cache changed"):
            self.prepare()

    def test_manifest_change_invalidates_cache(self):
        self.prepare()
        with (self.source / "MANIFEST.sha256").open("a") as stream:
            stream.write("\n")
        with self.assertRaisesRegex(ValueError, "Published local cache changed"):
            self.prepare()

    def test_bad_source_hash_never_publishes_cache(self):
        (self.source / "data/geometry/table.bin").write_bytes(b"corruption")
        with self.assertRaisesRegex(ValueError, "do not match"):
            self.prepare()
        self.assertFalse(self.target.exists())
        self.assertEqual(list(self.target.parent.iterdir()), [])

    def test_symlink_in_cache_rejected(self):
        self.prepare()
        (self.target / "data/spatial/link").symlink_to(self.image)
        with self.assertRaisesRegex(ValueError, "symlinks"):
            self.prepare()

    def test_receipt_has_data_manifest_identity(self):
        self.prepare()
        receipt = json.loads((self.target / ".verified.json").read_text())
        self.assertEqual(receipt["manifest_sha256"], cache.digest(self.source / "MANIFEST.sha256"))
        self.assertEqual(set(receipt["inventory"]), {
            f"data/{tree}/table.bin" for tree in cache.TREES
        })


if __name__ == "__main__":
    unittest.main()
