#!/usr/bin/env python3
"""Exercise the split-release trust boundaries with real small artifact fixtures."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

SPEC = importlib.util.spec_from_file_location("release_artifacts", Path(__file__).resolve().parents[1] / "release-artifacts.py")
art = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(art)
IMAGE = "sha256:" + "a" * 64
COMMIT = "b" * 40


class Artifacts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.project = self.root / "project"
        self.project.mkdir()
        self.repo = self.root / "repo"
        self.repo.mkdir()
        subprocess.run(["git", "init", "-q", str(self.repo)], check=True)
        self.put(self.repo / "deploy/Containerfile", 'FROM example\nLABEL io.islet-explorer.runtime-layout="split-v1"\nCOPY deploy/support.txt /support\n')
        self.put(self.repo / "deploy/containerignore", "**\n!deploy/**\n")
        self.put(self.repo / "deploy/base-image.lock", "# fixture base\nexample@sha256:" + "d" * 64 + "\n")
        self.put(self.repo / "deploy/support.txt", "runtime support\n")
        self.put(self.repo / "app/shiny_app/app.R", "shinyApp(ui, server)\n")
        self.put(self.repo / "app/shiny_app/R/components.R", "card <- 1\n")
        self.put(self.repo / "app/docs/logo.png", "logo\n")
        self.put(self.repo / "scripts/verify_release_data.R", "print('verify')\n")
        (self.repo / "app/shiny_app/www").mkdir()
        (self.repo / "app/shiny_app/www/logo.png").symlink_to("../../docs/logo.png")
        self.commit = self.commit_repo()
        self.app = art.package_app(self.project, "app1", self.repo, self.commit)
        self.fingerprint = art.load_json(self.app / "app.json")["runtime_inputs_sha256"]
        self.source = self.project / "releases/source1"
        self.put(self.source / "data/manifest.json", json.dumps({"schema_version": 1, "git_commit": COMMIT, "git_dirty": False}))
        for i in range(15):
            self.put(self.source / f"data/images/{i}.ome.tiff", f"image {i}\n")
        self.put(self.source / "data/cells/example.parquet", "fake parquet\n")
        self.put(self.source / "release.env", f"MANIFEST_VERSION=1\nRELEASE=source1\nDATA_SCHEMA=1\nGIT_COMMIT={COMMIT}\n")
        art.write_manifest(self.source)
        self.data = art.register_data(self.project, "data1", "source1")
        self.scan = self.make_scan()
        self.runtime = art.register_runtime(self.project, "runtime1", self.scan, self.fingerprint)
        validation_path = self.validation("runtime")
        validation = art.load_json(validation_path)
        validation["checks"]["scientific"] = "pass"
        self.put(validation_path, json.dumps(validation))
        art.assemble(self.project, "baseline", "app1", "data1", "runtime1", validation_path, "runtime")

    def tearDown(self):
        for current, directories, _ in os.walk(self.root):
            Path(current).chmod(0o755)
        self.temp.cleanup()

    @staticmethod
    def put(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and not path.is_symlink():
            path.chmod(0o644)
        path.write_text(value)

    def commit_repo(self):
        subprocess.run(["git", "-C", str(self.repo), "add", "."], check=True)
        subprocess.run(["git", "-C", str(self.repo), "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                        "commit", "-qm", "fixture"], check=True)
        return art.git(self.repo, "rev-parse", "HEAD").decode().strip()

    def make_scan(self):
        root = self.root / "scan"
        root.mkdir()
        prefix = IMAGE.replace(":", "_")
        scan = {"IMAGE_ID": IMAGE, "IMAGE_REF": IMAGE, "RESULT": "pass-no-high-or-critical"}
        files = [("SBOM_SHA256", ".cyclonedx.json", "{}"),
                 ("RAW_SBOM_SHA256", ".syft.cyclonedx.json", "{}"),
                 ("SBOM_CORRECTIONS_SHA256", ".sbom-corrections.json", "{}"),
                 ("REPORT_SHA256", ".grype.json", '{"matches":[]}'),
                 ("R_PACKAGES_SHA256", ".R-packages.tsv", "package\tversion\n"),
                 ("PYTHON_PACKAGES_SHA256", ".python-packages.txt", "test==1\n"),
                 ("BROTLI_REPORT_SHA256", ".brotli.json", json.dumps({
                     "schema": "islet.brotli_verification/2", "cve": "CVE-2026-111497", "result": "pass", "errors": [],
                     "duckdb": {"version": "1.5.2", "result": 0, "error_code": -17, "library_sha256": "c" * 64}}))]
        for key, suffix, contents in files:
            self.put(root / (prefix + suffix), contents)
            scan[key] = art.digest(root / (prefix + suffix))
        self.put(root / (prefix + ".scan.txt"), "".join(f"{k}={v}\n" for k, v in scan.items()))
        self.put(root / "published-image.txt", f"LOCAL_IMAGE_ID={IMAGE}\nGIT_COMMIT={COMMIT}\nBASE_IMAGE=example@sha256:{'d'*64}\n")
        return root

    def validation(self, kind="app", app="app1", runtime="runtime1"):
        resolved = art.resolve_candidate(self.project, app, "data1", runtime)
        result = {"app_manifest_sha256": resolved["APP_MANIFEST_SHA256"], "image_ref": resolved["IMAGE_REF"],
                  "data_receipt_sha256": resolved["DATA_RECEIPT_SHA256"], "change_kind": kind,
                  "checks": {"app_tests": "pass", "smoke": "pass"}}
        path = self.root / "validation.json"
        self.put(path, json.dumps(result))
        return path

    def assemble(self, kind="app"):
        return art.assemble(self.project, "release1", "app1", "data1", "runtime1", self.validation(kind), kind, "baseline")

    def test_full_combination_independent_commits_and_readonly_bundle(self):
        release = self.assemble()
        resolved = art.resolve_release(release)
        self.assertEqual(resolved["APP_COMMIT"], self.commit)
        self.assertEqual(resolved["DATA_COMMIT"], COMMIT)
        self.assertEqual(resolved["RUNTIME_COMMIT"], COMMIT)
        self.assertEqual(resolved["DATA_ID"], "data1")
        self.assertEqual(resolved["APP_DIR"], str(self.app / "app"))
        self.assertEqual((self.app / "app/shiny_app/app.R").stat().st_mode & 0o222, 0)
        self.assertEqual((self.app / "app/shiny_app/www/logo.png").read_text(), "logo\n")
        self.assertTrue(resolved["DESCRIPTOR_SHA256"])

    def test_app_package_ignores_dirty_and_untracked_work(self):
        self.put(self.repo / "app/shiny_app/app.R", "uncommitted app\n")
        self.put(self.repo / "scripts/research/unrelated.R", "untracked research\n")
        app2 = art.package_app(self.project, "app2", self.repo, self.commit)
        self.assertEqual((app2 / "app/shiny_app/app.R").read_text(), "shinyApp(ui, server)\n")
        self.assertEqual((self.repo / "app/shiny_app/app.R").read_text(), "uncommitted app\n")
        self.assertFalse((app2 / "scripts/research").exists())

    def test_package_excludes_credentials_and_microscopy(self):
        self.put(self.repo / "app/shiny_app/.Renviron", "SECRET=fixture\n")
        self.put(self.repo / "app/shiny_app/www/local_images/private.ome.tiff", "pixels")
        self.put(self.repo / "app/shared/.env", "SECRET=fixture\n")
        self.put(self.repo / "app/docs/local_images/private.ome.tiff", "pixels")
        commit = self.commit_repo()
        app2 = art.package_app(self.project, "app2", self.repo, commit)
        self.assertFalse((app2 / "app/shiny_app/.Renviron").exists())
        self.assertFalse((app2 / "app/shared/.env").exists())
        self.assertFalse((app2 / "app/docs/local_images").exists())
        self.assertEqual(list((app2 / "app/shiny_app/www/local_images").iterdir()), [])

    def test_entire_app_tree_additions_and_removals_are_packaged(self):
        asset = "app/docs/help.svg"
        helper = "app/shared/helpers.R"
        self.put(self.repo / asset, "<svg><title>Component help</title></svg>\n")
        self.put(self.repo / helper, "card_title <- function(x) paste('Card', x)\n")
        added_commit = self.commit_repo()
        app2 = art.package_app(self.project, "app2", self.repo, added_commit)
        for relative in (asset, helper):
            self.assertEqual((app2 / relative).read_bytes(), (self.repo / relative).read_bytes())
            self.assertIn(relative, art.read_manifest(app2 / "MANIFEST.sha256"))
        self.assertEqual(art.resolve_candidate(self.project, "app2", "data1", "runtime1")["IMAGE_REF"], IMAGE)

        (self.repo / asset).unlink()
        (self.repo / helper).unlink()
        removed_commit = self.commit_repo()
        app3 = art.package_app(self.project, "app3", self.repo, removed_commit)
        for relative in (asset, helper):
            self.assertFalse((app3 / relative).exists())
            self.assertNotIn(relative, art.read_manifest(app3 / "MANIFEST.sha256"))
            self.assertTrue((app2 / relative).exists())
        self.assertEqual(art.resolve_candidate(self.project, "app3", "data1", "runtime1")["IMAGE_REF"], IMAGE)

    def test_arbitrary_app_changes_reuse_runtime(self):
        self.put(self.repo / "app/shiny_app/R/components.R", "new_card <- 2\n")
        self.put(self.repo / "app/shiny_app/www/behavior.js", "console.log('new behavior');\n")
        commit = self.commit_repo()
        app2 = art.package_app(self.project, "app2", self.repo, commit)
        self.assertEqual(art.verify_app(app2)["runtime_inputs_sha256"], self.fingerprint)
        release = art.assemble(self.project, "release2", "app2", "data1", "runtime1", self.validation(app="app2"), "app", "baseline")
        self.assertEqual(art.resolve_release(release)["APP_COMMIT"], commit)

    def test_runtime_dependency_change_fails_compatibility(self):
        self.put(self.repo / "deploy/support.txt", "changed runtime\n")
        commit = self.commit_repo()
        art.package_app(self.project, "app2", self.repo, commit)
        with self.assertRaisesRegex(ValueError, "runtime inputs differ"):
            art.resolve_candidate(self.project, "app2", "data1", "runtime1")

    def test_science_data_runtime_require_scientific_evidence(self):
        for kind in ("science", "data", "runtime"):
            with self.subTest(kind=kind):
                path = self.validation(kind)
                with self.assertRaisesRegex(ValueError, "scientific"):
                    art.assemble(self.project, "release-" + kind, "app1", "data1", "runtime1", path, kind, "baseline")
                validation = art.load_json(path)
                validation["checks"]["scientific"] = "pass"
                self.put(path, json.dumps(validation))
                release = art.assemble(self.project, "release-" + kind, "app1", "data1", "runtime1", path, kind, "baseline")
                self.assertEqual(art.resolve_release(release)["IMAGE_REF"], IMAGE)

    def test_validation_is_bound_to_exact_combination(self):
        path = self.validation()
        value = art.load_json(path)
        for key in ("app_manifest_sha256", "image_ref", "data_receipt_sha256", "change_kind"):
            changed = dict(value)
            changed[key] = "wrong"
            self.put(path, json.dumps(changed))
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "not bound"):
                art.assemble(self.project, "release1", "app1", "data1", "runtime1", path, "app", "baseline")
        self.assertFalse((self.project / "releases/release1").exists())

    def test_reuse_never_opens_microscopy_contents(self):
        release = self.assemble()
        original = Path.open
        def guarded(path, *args, **kwargs):
            if str(path).endswith(".ome.tiff"):
                raise AssertionError("Microscopy content read during app release")
            return original(path, *args, **kwargs)
        with mock.patch.object(Path, "open", guarded):
            art.resolve_release(release)

    def test_full_data_registration_checks_microscopy(self):
        self.put(self.source / "data/images/0.ome.tiff", "changed pixels\n")
        with self.assertRaisesRegex(ValueError, "Data checksum mismatch"):
            art.register_data(self.project, "data2", "source1")

    def test_reuse_rejects_changed_inventory(self):
        self.put(self.source / "data/images/0.ome.tiff", "changed pixels\n")
        with self.assertRaisesRegex(ValueError, "inventory changed"):
            art.verify_data(self.project, self.data)

    def test_reuse_rejects_new_files(self):
        self.put(self.source / "data/new.txt", "new")
        with self.assertRaisesRegex(ValueError, "inventory differs"):
            art.verify_data(self.project, self.data)

    def test_data_receipt_and_original_manifest_are_bound(self):
        self.put(self.source / "MANIFEST.sha256", (self.source / "MANIFEST.sha256").read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "manifest changed"):
            art.verify_data(self.project, self.data)

    def test_metadata_contents_always_checked(self):
        manifest = self.source / "data/manifest.json"
        previous = manifest.stat()
        text = manifest.read_text().replace('"schema_version": 1', '"schema_version": 5')
        self.put(manifest, text)
        os.utime(manifest, ns=(previous.st_atime_ns, previous.st_mtime_ns))
        with self.assertRaisesRegex(ValueError, "manifest hash mismatch"):
            art.verify_data(self.project, self.data)

    def test_app_tamper_and_extra_file_rejected(self):
        self.put(self.app / "app/shiny_app/app.R", "tampered\n")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            art.verify_app(self.app)

    def test_descriptor_tamper_rejected(self):
        release = self.assemble()
        self.put(release / "deployment.json", "{}\n")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            art.resolve_release(release)

    def test_security_evidence_tamper_rejected(self):
        self.put(self.runtime / "attestations" / (IMAGE.replace(":", "_") + ".grype.json"), '{"matches":[{}]}')
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            art.verify_runtime(self.runtime)

    def test_bad_runtime_security_rejected_before_publication(self):
        prefix = IMAGE.replace(":", "_")
        path = self.scan / (prefix + ".grype.json")
        self.put(path, '{"matches":[{"vulnerability":{"severity":"High"}}]}')
        scanpath = self.scan / (prefix + ".scan.txt")
        fields = art.fields(scanpath)
        fields["REPORT_SHA256"] = art.digest(path)
        self.put(scanpath, "".join(f"{k}={v}\n" for k, v in fields.items()))
        with self.assertRaisesRegex(ValueError, "High/Critical"):
            art.register_runtime(self.project, "runtime2", self.scan, self.fingerprint)
        self.assertFalse((self.project / "runtimes/runtime2").exists())

    def test_symlink_escape_rejected(self):
        self.put(self.repo / "app/shiny_app/www/escape", "placeholder")
        (self.repo / "app/shiny_app/www/escape").unlink()
        (self.repo / "app/shiny_app/www/escape").symlink_to("../../../../outside")
        commit = self.commit_repo()
        with self.assertRaisesRegex(ValueError, "escapes"):
            art.package_app(self.project, "app2", self.repo, commit)

    def test_path_traversal_and_duplicate_metadata_rejected(self):
        for value in ("../app", ".", "..", "bad name", "bad\nname"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                art.identifier(value)
        path = self.root / "duplicate.json"
        self.put(path, '{"id":"a","id":"b"}')
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            art.load_json(path)

    def test_runtime_published_fingerprint_must_match(self):
        published = self.scan / "published-image.txt"
        self.put(published, published.read_text() + f"RUNTIME_LAYOUT=split-v1\nRUNTIME_INPUTS_SHA256={'0'*64}\n")
        with self.assertRaisesRegex(ValueError, "fingerprint/layout"):
            art.register_runtime(self.project, "runtime2", self.scan, self.fingerprint)

    def test_registration_preserves_source_hardlinks_and_permissions(self):
        original = self.source / "data/images/0.ome.tiff"
        linked = self.root / "shared-image.ome.tiff"
        linked.hardlink_to(original)
        before = original.stat()
        art.register_data(self.project, "data2", "source1")
        after = original.stat()
        self.assertEqual((before.st_ino, before.st_mode, before.st_mtime_ns),
                         (after.st_ino, after.st_mode, after.st_mtime_ns))
        self.assertEqual(after.st_ino, linked.stat().st_ino)

    def test_descriptor_pins_app_manifest_even_if_manifest_regenerated(self):
        release = self.assemble()
        self.put(self.app / "app/shiny_app/app.R", "changed app\n")
        (self.app / "MANIFEST.sha256").chmod(0o644)
        art.write_manifest(self.app)
        with self.assertRaisesRegex(ValueError, "artifact identity changed"):
            art.resolve_release(release)

    def new_data_source(self):
        source = self.project / "releases/new-data"
        shutil.copytree(self.source / "data", source / "data")
        return source

    def test_new_data_seals_without_runtime_or_app(self):
        source = self.new_data_source()
        art.seal_data(self.project, "new-data", "1", 15)
        metadata = art.fields(source / "release.env")
        self.assertNotIn("IMAGE_REF", metadata)
        data2 = art.register_data(self.project, "data2", "new-data")
        receipt, actual = art.verify_data(self.project, data2)
        self.assertEqual(actual, source)
        self.assertEqual(receipt["git_commit"], COMMIT)

    def test_new_data_schema_and_count_conflicts_fail_before_sealing(self):
        source = self.new_data_source()
        for schema, count in (("5", 15), ("1", 14), ("1", 0)):
            with self.subTest(schema=schema, count=count), self.assertRaises(ValueError):
                art.seal_data(self.project, "new-data", schema, count)
        self.assertFalse((source / "release.env").exists())
        self.assertFalse((source / "MANIFEST.sha256").exists())

    def test_seal_data_cannot_overwrite_existing_metadata(self):
        self.new_data_source()
        art.seal_data(self.project, "new-data", "1", 15)
        with self.assertRaisesRegex(ValueError, "fresh source"):
            art.seal_data(self.project, "new-data", "1", 15)

    def test_app_release_requires_existing_split_baseline(self):
        with self.assertRaisesRegex(ValueError, "require --base-release"):
            art.assemble(self.project, "release1", "app1", "data1", "runtime1", self.validation(), "app")
        with self.assertRaises((ValueError, FileNotFoundError)):
            art.assemble(self.project, "release1", "app1", "data1", "runtime1", self.validation(), "app", "source1")

    def test_changed_runtime_cannot_be_claimed_app_only(self):
        # A new registered runtime artifact is independently versioned even if
        # this tiny fixture deliberately keeps the same image ID.
        art.register_runtime(self.project, "runtime2", self.scan, self.fingerprint)
        with self.assertRaisesRegex(ValueError, "baseline runtime and data"):
            art.assemble(self.project, "release1", "app1", "data1", "runtime2",
                         self.validation(runtime="runtime2"), "app", "baseline")

    def test_changed_data_cannot_be_claimed_app_only(self):
        data2 = art.register_data(self.project, "data2", "source1")
        validation_path = self.validation()
        validation = art.load_json(validation_path)
        validation["data_receipt_sha256"] = art.digest(data2 / "data.json")
        self.put(validation_path, json.dumps(validation))
        with self.assertRaisesRegex(ValueError, "baseline runtime and data"):
            art.assemble(self.project, "release1", "app1", "data2", "runtime1",
                         validation_path, "app", "baseline")

    def test_baseline_descriptor_is_bound_but_not_recursively_resolved(self):
        release = self.assemble()
        # Baseline metadata validation does not re-run its old runtime/app/data
        # or recurse through older ancestors. It checks only its pinned files.
        descriptor = art.load_json(release / "deployment.json")
        self.assertEqual(descriptor["base_release"], "baseline")
        self.assertEqual(descriptor["base_descriptor_sha256"],
                         art.digest(self.project / "releases/baseline/deployment.json"))
        baseline = self.project / "releases/baseline"
        path = baseline / "deployment.json"
        self.put(path, path.read_text() + "\n")
        (baseline / "MANIFEST.sha256").chmod(0o644)
        art.write_manifest(baseline)
        with self.assertRaisesRegex(ValueError, "Baseline deployment descriptor changed"):
            art.resolve_release(release)

    def test_invalid_legacy_baseline_never_reads_microscopy(self):
        validation = self.validation()
        original = Path.open
        def guarded(path, *args, **kwargs):
            if str(path).endswith(".ome.tiff"):
                raise AssertionError("Legacy baseline rejection read microscopy")
            return original(path, *args, **kwargs)
        with mock.patch.object(Path, "open", guarded), self.assertRaisesRegex(ValueError, "split deployment"):
            art.assemble(self.project, "release1", "app1", "data1", "runtime1", validation, "app", "source1")

    def test_recovery_app_check_rehashes_actual_contents(self):
        expected = art.digest(self.app / "MANIFEST.sha256")
        self.assertEqual(art.verify_app_identity(self.app, expected), self.app)
        self.put(self.app / "app/shiny_app/app.R", "modified since activation\n")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            art.verify_app_identity(self.app, expected)

    def test_recovery_app_check_rejects_changed_manifest_and_unsafe_path(self):
        with self.assertRaisesRegex(ValueError, "previously activated"):
            art.verify_app_identity(self.app, "0" * 64)
        linked = self.project / "apps/linked"
        linked.symlink_to(self.app)
        with self.assertRaisesRegex(ValueError, "symlinks"):
            art.verify_app_identity(linked, art.digest(self.app / "MANIFEST.sha256"))

    def test_registration_rejects_dirty_or_unproven_data_producer(self):
        path = self.source / "data/manifest.json"
        original = art.load_json(path)
        for dirty in (True, None):
            manifest = dict(original)
            if dirty is None:
                manifest.pop("git_dirty")
            else:
                manifest["git_dirty"] = dirty
            self.put(path, json.dumps(manifest))
            art.write_manifest(self.source)
            with self.subTest(dirty=dirty), self.assertRaisesRegex(ValueError, "clean committed"):
                art.register_data(self.project, "data2", "source1")

    def test_seal_rejects_dirty_or_unproven_data_producer(self):
        source = self.new_data_source()
        path = source / "data/manifest.json"
        original = art.load_json(path)
        for dirty in (True, None):
            manifest = dict(original)
            if dirty is None:
                manifest.pop("git_dirty")
            else:
                manifest["git_dirty"] = dirty
            self.put(path, json.dumps(manifest))
            with self.subTest(dirty=dirty), self.assertRaisesRegex(ValueError, "clean committed"):
                art.seal_data(self.project, "new-data", "1", 15)
        self.assertFalse((source / "release.env").exists())

    def test_base_lock_only_change_rejects_runtime_reuse(self):
        self.put(self.repo / "deploy/base-image.lock", "example@sha256:" + "e" * 64 + "\n")
        commit = self.commit_repo()
        app2 = art.package_app(self.project, "app2", self.repo, commit)
        self.assertEqual(art.verify_app(app2)["runtime_inputs_sha256"], self.fingerprint)
        with self.assertRaisesRegex(ValueError, "base-image requirement differs"):
            art.resolve_candidate(self.project, "app2", "data1", "runtime1")

    def test_runtime_base_requirement_is_bound_to_publication(self):
        path = self.runtime / "runtime.json"
        metadata = art.load_json(path)
        metadata["base_image"] = "example@sha256:" + "e" * 64
        self.put(path, json.dumps(metadata))
        (self.runtime / "MANIFEST.sha256").chmod(0o644)
        art.write_manifest(self.runtime)
        with self.assertRaisesRegex(ValueError, "bound published evidence"):
            art.verify_runtime(self.runtime)

    def test_resolution_includes_bound_base_image(self):
        resolved = art.resolve_candidate(self.project, "app1", "data1", "runtime1")
        self.assertEqual(resolved["RUNTIME_BASE_IMAGE"], "example@sha256:" + "d" * 64)

    def test_immutable_artifacts_cannot_be_overwritten(self):
        with self.assertRaisesRegex(ValueError, "overwrite"):
            art.package_app(self.project, "app1", self.repo, self.commit)
        self.assertEqual(art.verify_app(self.app)["git_commit"], self.commit)


if __name__ == "__main__":
    unittest.main()
