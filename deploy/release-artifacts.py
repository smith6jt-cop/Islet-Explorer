#!/usr/bin/env python3
"""Immutable app/runtime/data artifacts for rootless PubApps deployments.

Data registration hashes every data file once. Reuse verifies the receipt,
metadata and inventory only; it relies on published data remaining immutable.
All paths emitted by resolve are validated, and are data, never shell code.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import runpy
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile

SCRIPT_DIR = Path(__file__).resolve().parent
ID_RE = r"[A-Za-z0-9][A-Za-z0-9._-]*"
HASH_RE = r"[0-9a-f]{64}"
COMMIT_RE = r"[0-9a-f]{40}"
IMAGE_RE = r"(?:sha256:[0-9a-f]{64}|[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64})"
SCHEMA5_ARTIFACTS = {"core_area", "targets_all", "markers_all", "comp", "combined_pairs",
                     "cell_distance", "signed_distance_kde", "donor_tissue", "geometry",
                     "adata_ins_root", "barlow_islet_features", "islet_spatial_lookup",
                     "donor_metadata", "cells", "images", "uniform_release",
                     "endocrine_objects_by_donor", "protein_phenotype_rules", "embedding_metadata"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identifier(value):
    require(isinstance(value, str) and re.fullmatch(ID_RE, value), "Invalid artifact identifier")
    require(value not in (".", ".."), "Invalid artifact identifier")
    return value


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def load_json(path):
    def unique(items):
        result = {}
        for key, value in items:
            require(key not in result, f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=unique)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def fields(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        key, sep, value = line.partition("=")
        require(sep and key not in result, f"Malformed or duplicate field in {path}")
        result[key] = value
    return result


def safe_relative(name):
    require(isinstance(name, str) and name and "\\" not in name and
            not any(ord(c) < 32 for c in name), "Unsafe relative path")
    p = PurePosixPath(name)
    require(not p.is_absolute() and all(x not in (".", "..", "") for x in name.split("/")),
            f"Unsafe relative path: {name}")
    return name


def project_path(value):
    path = Path(value).absolute()
    require(re.fullmatch(r"/[A-Za-z0-9._/-]+", str(path)),
            "Project root must be absolute and contain no whitespace or shell metacharacters")
    require(path == path.resolve(), "Project root may not contain symlinks or traversal")
    return path


def artifact_path(project, category, name):
    path = project / category / identifier(name)
    require(path == path.resolve(), f"Artifact path may not contain symlinks: {path}")
    return path


def regular(path):
    require(path.is_file() and not path.is_symlink(), f"Expected regular file: {path}")
    return path


def tree_files(root):
    """List files and file symlinks without following directory symlinks."""
    result = {}
    require(root.is_dir() and not root.is_symlink(), f"Expected real directory: {root}")
    for current, directories, names in os.walk(root, followlinks=False):
        for name in directories:
            require(not (Path(current) / name).is_symlink(), "Directory symlinks are not permitted")
        for name in names:
            path = Path(current) / name
            relative = safe_relative(path.relative_to(root).as_posix())
            require(path.is_file(), f"Missing target or special file: {relative}")
            require(path.resolve().is_relative_to(root.resolve()), f"Symlink escapes artifact: {relative}")
            result[relative] = path
    return result


def read_manifest(path):
    regular(path)
    result = {}
    for line in path.read_text().splitlines():
        match = re.fullmatch(r"([0-9a-f]{64}) [ *](.+)", line)
        require(match is not None, f"Malformed checksum entry in {path}")
        name = safe_relative(match[2])
        require(name not in result, f"Duplicate checksum path: {name}")
        result[name] = match[1]
    require(result, f"Empty manifest: {path}")
    return result


def write_manifest(root):
    entries = tree_files(root)
    entries.pop("MANIFEST.sha256", None)
    (root / "MANIFEST.sha256").write_text("".join(
        f"{digest(path)}  {name}\n" for name, path in sorted(entries.items())))


def verify_manifest(root):
    expected = read_manifest(root / "MANIFEST.sha256")
    actual = tree_files(root)
    actual.pop("MANIFEST.sha256", None)
    require(set(expected) == set(actual), f"Artifact inventory changed: {root}")
    for name, checksum in expected.items():
        require(digest(actual[name]) == checksum, f"Artifact hash mismatch: {root / name}")
    return expected


def publish(project, category, name, build):
    target = artifact_path(project, category, name)
    require(not target.exists(), f"Refusing to overwrite immutable artifact: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{name}-", dir=target.parent))
    try:
        build(stage)
        write_manifest(stage)
        # Read-only is additional accident protection, not a substitute for hashes.
        for path in tree_files(stage).values():
            if not path.is_symlink():
                path.chmod(0o555 if path.stat().st_mode & 0o111 else 0o444)
        for current, directories, _ in os.walk(stage, topdown=False):
            for directory in directories:
                (Path(current) / directory).chmod(0o555)
        stage.chmod(0o555)
        require(not target.exists(), f"Refusing to overwrite immutable artifact: {target}")
        stage.rename(target)
    finally:
        if stage.exists():
            for current, _, _ in os.walk(stage):
                Path(current).chmod(0o755)
            shutil.rmtree(stage)
    return target


def forbidden(name):
    parts = PurePosixPath(name).parts
    leaf = parts[-1]
    return (any(x in (".git", "__pycache__", "node_modules", "local_images") for x in parts) or
            leaf in (".Renviron", ".Rhistory", "APP_LOADED_FLAG", "restart.txt") or
            leaf == ".env" or leaf.startswith(".env.") or
            leaf.endswith((".pem", ".key", ".bak", ".pyc", ".pyo", "~")))


def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args])


def package_app(project, name, repo, revision):
    commit = git(repo, "rev-parse", "--verify", revision + "^{commit}").decode().strip()
    require(re.fullmatch(COMMIT_RE, commit), "Expected a full Git commit")
    containerfile = git(repo, "show", f"{commit}:deploy/Containerfile").decode()
    require("io.islet-explorer.runtime-layout" in containerfile and "split-v1" in containerfile,
            "App source must use the split-v1 runtime layout")
    runtime_requirements = runpy.run_path(str(SCRIPT_DIR / "runtime-inputs.py"))
    fingerprint = runtime_requirements["runtime_fingerprint"](repo, commit)
    base_image = runtime_requirements["base_image_requirement"](repo, commit)
    archive = git(repo, "archive", "--format=tar", commit, "app",
                  "scripts/verify_release_data.R")

    def build(stage):
        links = {}
        with tarfile.open(fileobj=io.BytesIO(archive)) as stream:
            for member in stream:
                relative = member.name.rstrip("/")
                if not relative or forbidden(relative):
                    continue
                safe_relative(relative)
                path = stage / relative
                if member.isdir():
                    path.mkdir(parents=True, exist_ok=True)
                elif member.isfile():
                    path.parent.mkdir(parents=True, exist_ok=True)
                    require(not path.exists(), f"Duplicate archive path: {relative}")
                    with stream.extractfile(member) as source, path.open("wb") as output:
                        shutil.copyfileobj(source, output)
                    path.chmod(0o755 if member.mode & 0o111 else 0o644)
                elif member.issym():
                    # Create links only after every regular file has been extracted.
                    links[relative] = member.linkname
                else:
                    raise ValueError(f"Unsupported archive entry: {relative}")
        for relative, target in links.items():
            path = stage / relative
            require(not Path(target).is_absolute() and not any(ord(c) < 32 for c in target),
                    f"Unsafe symlink: {relative}")
            require((path.parent / target).resolve().is_relative_to(stage),
                    f"Symlink escapes app bundle: {relative}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.symlink_to(target)
        (stage / "app/shiny_app/www/local_images").mkdir(parents=True, exist_ok=True)
        regular(stage / "app/shiny_app/app.R")
        regular(stage / "scripts/verify_release_data.R")
        write_json(stage / "app.json", {"schema": "islet.app/1", "id": name,
                   "git_commit": commit, "runtime_inputs_sha256": fingerprint, "base_image": base_image, "symlinks": links})
    return publish(project, "apps", name, build)


def verify_app(root):
    verify_manifest(root)
    meta = load_json(regular(root / "app.json"))
    require(meta.get("schema") == "islet.app/1" and meta.get("id") == root.name,
            "Invalid app metadata")
    require(re.fullmatch(COMMIT_RE, meta.get("git_commit", "")) and
            re.fullmatch(HASH_RE, meta.get("runtime_inputs_sha256", "")) and
            re.fullmatch(r"[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}", meta.get("base_image", "")),
            "Invalid app provenance or missing committed base-image requirement")
    actual_links = {name: os.readlink(path) for name, path in tree_files(root).items() if path.is_symlink()}
    require(actual_links == meta.get("symlinks"), "App symlink targets changed")
    regular(root / "app/shiny_app/app.R")
    regular(root / "scripts/verify_release_data.R")
    require((root / "app/shiny_app/www/local_images").is_dir(), "Missing image mountpoint")
    return meta


def verify_app_identity(root, expected_manifest):
    require(root.parent.name == "apps", "App bundle must be under PROJECT/apps/ID")
    project = project_path(root.parent.parent)
    require(artifact_path(project, "apps", root.name) == root, "Invalid app bundle path")
    require(re.fullmatch(HASH_RE, expected_manifest), "Invalid expected app manifest hash")
    require(digest(regular(root / "MANIFEST.sha256")) == expected_manifest,
            "App manifest differs from the previously activated app")
    verify_app(root)
    return root


def data_identity(source):
    entries = read_manifest(source / "MANIFEST.sha256")
    require("release.env" in entries and "data/manifest.json" in entries, "Unbound data identity")
    require(digest(regular(source / "release.env")) == entries["release.env"], "Data release metadata changed")
    require(digest(regular(source / "data/manifest.json")) == entries["data/manifest.json"],
            "Data manifest hash mismatch")
    meta = fields(source / "release.env")
    require(meta.get("RELEASE") == source.name, "Source release ID mismatch")
    version = meta.get("MANIFEST_VERSION")
    require(version in ("1", "2"), "Register data from a version 1 or 2 data release")
    count = meta.get("EXPECTED_IMAGES", "15" if version == "1" else "")
    require(re.fullmatch(r"[1-9][0-9]*", count), "Missing or invalid expected image count")
    require(version != "1" or count == "15", "Legacy version 1 requires 15 images")
    manifest = load_json(source / "data/manifest.json")
    require(manifest.get("git_dirty") is False, "Data producer must be a clean committed revision (git_dirty=false)")
    schema = meta.get("DATA_SCHEMA")
    require(schema in ("1", "5"), "New split releases support accepted data schemas 1 and 5")
    require(str(manifest.get("schema_version")) == schema, "Data schema mismatch")
    commit = manifest.get("git_commit", "")
    require(re.fullmatch(COMMIT_RE, commit), "Data producer commit is missing or invalid")
    require(commit == meta.get("GIT_COMMIT"), "Legacy data producer and release commit differ")
    if schema == "5":
        require(manifest.get("expected_images") == int(count) and
                SCHEMA5_ARTIFACTS.issubset(manifest.get("artifacts", {})), "Incomplete schema-5 data manifest")
    return entries, schema, commit, int(count)


def data_inventory(source, entries, full=False):
    actual = {"data/" + name: path for name, path in tree_files(source / "data").items()}
    expected = {name: checksum for name, checksum in entries.items() if name.startswith("data/")}
    require(set(actual) == set(expected), "Data inventory differs from its original manifest")
    inventory = {}
    for name, path in sorted(actual.items()):
        if full:
            require(digest(path) == expected[name], f"Data checksum mismatch: {name}")
        status = path.stat()
        inventory[name] = {"size": status.st_size, "mtime_ns": status.st_mtime_ns,
                           "sha256": expected[name]}
        if path.is_symlink():
            inventory[name]["symlink"] = os.readlink(path)
    return inventory


def image_count(inventory):
    # Match legacy find -type f: sidecar symlinks are accepted, image symlinks are not.
    return sum(bool(re.search(r"[.]ome[.]tiff?$", name, re.I)) and "symlink" not in info
               for name, info in inventory.items())


def seal_data(project, source_name, schema, expected_images):
    """Seal a newly built data directory without selecting any runtime image."""
    source = artifact_path(project, "releases", source_name)
    require(source.is_dir(), "Stage the generated data under releases/ID/data first")
    require(schema in ("1", "5") and isinstance(expected_images, int) and expected_images > 0,
            "Expected an accepted schema and a positive approved image count")
    require(set(p.name for p in source.iterdir()) == {"data"},
            "seal-data requires a fresh source release containing only data/")
    manifest = load_json(regular(source / "data/manifest.json"))
    require(manifest.get("git_dirty") is False, "Data producer must be a clean committed revision (git_dirty=false)")
    require(str(manifest.get("schema_version")) == schema, "Requested schema differs from data manifest")
    commit = manifest.get("git_commit", "")
    require(re.fullmatch(COMMIT_RE, commit), "Data producer commit is missing or invalid")
    if schema == "5":
        require(manifest.get("expected_images") == expected_images and
                SCHEMA5_ARTIFACTS.issubset(manifest.get("artifacts", {})),
                "Incomplete schema-5 data manifest or conflicting approved image count")
    elif "expected_images" in manifest:
        require(manifest["expected_images"] == expected_images, "Conflicting approved image count")
    files = {"data/" + name: path for name, path in tree_files(source / "data").items()}
    require(sum(bool(re.search(r"[.]ome[.]tiff?$", name, re.I)) and not path.is_symlink()
                for name, path in files.items()) == expected_images, "Image count differs from approved count")
    for name in files:
        require(not forbidden(name), f"Forbidden filename in data release: {name}")
    metadata = (f"MANIFEST_VERSION=2\nRELEASE={source_name}\nDATA_SCHEMA={schema}\n"
                f"EXPECTED_IMAGES={expected_images}\nGIT_COMMIT={commit}\n")
    # Write both files privately, then publish without replacing existing files.
    # A crash before the second link leaves an unusable (not accepted) release.
    stage = Path(tempfile.mkdtemp(prefix=".seal-", dir=source))
    published = []
    try:
        (stage / "release.env").write_text(metadata)
        entries = {name: digest(path) for name, path in files.items()}
        entries["release.env"] = digest(stage / "release.env")
        (stage / "MANIFEST.sha256").write_text("".join(
            f"{checksum}  {name}\n" for name, checksum in sorted(entries.items())))
        for name in ("release.env", "MANIFEST.sha256"):
            (stage / name).chmod(0o444)
            os.link(stage / name, source / name)
            published.append(source / name)
    except BaseException:
        for path in published:
            path.unlink()
        raise
    finally:
        shutil.rmtree(stage)
    return source


def register_data(project, name, source_name):
    source = artifact_path(project, "releases", source_name)
    entries, schema, commit, count = data_identity(source)
    inventory = data_inventory(source, entries, full=True)
    require(image_count(inventory) == count, "Data image count differs from the declared cohort")
    # Bind metadata after reading contents; reject concurrent manifest edits.
    source_hash = digest(source / "MANIFEST.sha256")
    require(entries == read_manifest(source / "MANIFEST.sha256"), "Source manifest changed during registration")
    receipt = {"schema": "islet.data/1", "id": name, "source_release": source_name,
               "source_manifest_sha256": source_hash, "data_schema": schema, "git_commit": commit,
               "data_manifest_sha256": entries["data/manifest.json"], "expected_images": count,
               "verified_at": datetime.now(timezone.utc).isoformat(), "inventory": inventory}
    return publish(project, "datasets", name, lambda stage: write_json(stage / "data.json", receipt))


def verify_data(project, root):
    verify_manifest(root)
    receipt = load_json(regular(root / "data.json"))
    require(receipt.get("schema") == "islet.data/1" and receipt.get("id") == root.name, "Invalid data receipt")
    source = artifact_path(project, "releases", receipt["source_release"])
    require(digest(regular(source / "MANIFEST.sha256")) == receipt["source_manifest_sha256"],
            "Registered data manifest changed; full registration is required")
    entries, schema, commit, count = data_identity(source)
    require((schema, commit, count, entries["data/manifest.json"]) ==
            (receipt["data_schema"], receipt["git_commit"], receipt["expected_images"], receipt["data_manifest_sha256"]),
            "Registered data identity mismatch")
    inventory = data_inventory(source, entries)
    require(inventory == receipt["inventory"] and image_count(inventory) == count,
            "Registered data inventory changed; full registration is required")
    return receipt, source


def security_verify(root):
    runpy.run_path(str(SCRIPT_DIR / "verify-release-security.py"))["verify"](root)


def register_runtime(project, name, scan_dir, fingerprint):
    require(re.fullmatch(HASH_RE, fingerprint), "Invalid runtime input fingerprint")
    require(scan_dir.is_dir() and not scan_dir.is_symlink(), "Scan directory must be a real directory")
    published = fields(regular(scan_dir / "published-image.txt"))
    if "RUNTIME_INPUTS_SHA256" in published or "RUNTIME_LAYOUT" in published:
        require(published.get("RUNTIME_INPUTS_SHA256") == fingerprint and
                published.get("RUNTIME_LAYOUT") == "split-v1",
                "Published runtime fingerprint/layout differs from registration")
    local = published.get("LOCAL_IMAGE_ID")
    registry = published.get("REGISTRY_REF")
    published_id = published.get("IMAGE_ID")
    require((local and not registry and not published_id) or
            (not local and registry and published_id), "Ambiguous published image identity")
    image = registry or local
    commit = published.get("GIT_COMMIT", "")
    require(re.fullmatch(IMAGE_RE, image) and re.fullmatch(COMMIT_RE, commit), "Invalid runtime provenance")
    scans = list(scan_dir.glob("*.scan.txt"))
    require(len(scans) == 1, "Exactly one image scan is required")
    scan = fields(regular(scans[0]))
    require(scan.get("BROTLI_REPORT_SHA256"), "Current decoder evidence is required")
    base_image = published.get("BASE_IMAGE", "")
    require(re.fullmatch(r"[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}", base_image),
            "Runtime base image must be digest-pinned")

    def build(stage):
        evidence = stage / "attestations"
        evidence.mkdir()
        for path in scan_dir.iterdir():
            regular(path)
            safe_relative(path.name)
            require(not forbidden(path.name), "Forbidden filename in runtime evidence")
            shutil.copyfile(path, evidence / path.name)
        metadata = {"MANIFEST_VERSION": "3", "RELEASE": name, "IMAGE_REF": image,
                    "GIT_COMMIT": commit, "SECURITY_RESULT": scan.get("RESULT", "")}
        if "DUCKDB_PATCH_BY" in scan:
            metadata["DUCKDB_PATCH_BY"] = scan["DUCKDB_PATCH_BY"]
        (stage / "release.env").write_text("".join(f"{key}={value}\n" for key, value in metadata.items()))
        write_json(stage / "runtime.json", {"schema": "islet.runtime/1", "id": name,
                   "image_ref": image, "git_commit": commit, "runtime_inputs_sha256": fingerprint, "base_image": base_image})
        write_manifest(stage)
        security_verify(stage)
    return publish(project, "runtimes", name, build)


def verify_runtime(root):
    verify_manifest(root)
    meta = load_json(regular(root / "runtime.json"))
    require(meta.get("schema") == "islet.runtime/1" and meta.get("id") == root.name, "Invalid runtime metadata")
    require(re.fullmatch(IMAGE_RE, meta.get("image_ref", "")) and
            re.fullmatch(COMMIT_RE, meta.get("git_commit", "")) and
            re.fullmatch(HASH_RE, meta.get("runtime_inputs_sha256", "")) and
            re.fullmatch(r"[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}", meta.get("base_image", "")),
            "Invalid runtime identity or missing published base-image requirement")
    env = fields(root / "release.env")
    require(env.get("IMAGE_REF") == meta["image_ref"] and env.get("GIT_COMMIT") == meta["git_commit"],
            "Runtime metadata and security identity differ")
    require(env.get("SECURITY_RESULT"), "Current runtime security decision is required")
    published = fields(regular(root / "attestations/published-image.txt"))
    require(published.get("BASE_IMAGE") == meta["base_image"],
            "Runtime base-image requirement differs from its bound published evidence")
    if "RUNTIME_INPUTS_SHA256" in published or "RUNTIME_LAYOUT" in published:
        require(published.get("RUNTIME_INPUTS_SHA256") == meta["runtime_inputs_sha256"] and
                published.get("RUNTIME_LAYOUT") == "split-v1",
                "Published runtime fingerprint/layout differs from registered identity")
    security_verify(root)
    return meta


def resolve_candidate(project, app, data, runtime):
    app_root = artifact_path(project, "apps", app)
    data_root = artifact_path(project, "datasets", data)
    runtime_root = artifact_path(project, "runtimes", runtime)
    app_meta = verify_app(app_root)
    data_meta, source = verify_data(project, data_root)
    runtime_meta = verify_runtime(runtime_root)
    require(app_meta["runtime_inputs_sha256"] == runtime_meta["runtime_inputs_sha256"],
            "App runtime inputs differ; build and verify a compatible runtime image")
    require(app_meta["base_image"] == runtime_meta["base_image"],
            "App base-image requirement differs from the runtime; select a compatible verified runtime")
    return {"IMAGE_REF": runtime_meta["image_ref"], "RUNTIME_COMMIT": runtime_meta["git_commit"],
            "APP_COMMIT": app_meta["git_commit"], "APP_DIR": str(app_root / "app"),
            "APP_SCRIPTS_DIR": str(app_root / "scripts"), "DATA_DIR": str(source / "data"),
            "DATA_SCHEMA": data_meta["data_schema"], "DATA_COMMIT": data_meta["git_commit"],
            "DATA_ID": data, "DATA_RELEASE_DIR": str(source), "DATA_MANIFEST": str(source / "MANIFEST.sha256"),
            "APP_MANIFEST": str(app_root / "MANIFEST.sha256"),
            "APP_MANIFEST_SHA256": digest(app_root / "MANIFEST.sha256"),
            "DATA_RECEIPT_SHA256": digest(data_root / "data.json"),
            "RUNTIME_INPUT_SHA256": runtime_meta["runtime_inputs_sha256"],
            "RUNTIME_BASE_IMAGE": runtime_meta["base_image"], "DESCRIPTOR_SHA256": ""}


def check_validation(validation, resolved, kind):
    require(kind in ("app", "science", "data", "runtime"), "Unknown change kind")
    for key, expected in {"app_manifest_sha256": resolved["APP_MANIFEST_SHA256"],
                          "image_ref": resolved["IMAGE_REF"],
                          "data_receipt_sha256": resolved["DATA_RECEIPT_SHA256"],
                          "change_kind": kind}.items():
        require(validation.get(key) == expected, f"Validation is not bound to this combination: {key}")
    checks = validation.get("checks", {})
    required = ["app_tests", "smoke"] + (["scientific"] if kind != "app" else [])
    require(all(checks.get(key) == "pass" for key in required),
            "Missing passing validation: " + ", ".join(required))


def verify_descriptor_files(root):
    # Reject legacy/full data bundles before any content hashing. In particular,
    # an accidentally selected legacy baseline must never trigger a TIFF pass.
    expected = {"deployment.json", "validation.json", "MANIFEST.sha256"}
    require(root.is_dir() and {p.name for p in root.iterdir()} == expected,
            "Expected a split deployment descriptor; migrate legacy releases with change-kind runtime")
    for name in expected:
        regular(root / name)
    verify_manifest(root)


def check_baseline(project, kind, runtime_manifest, data_receipt, base_release, expected_hash=None):
    """Check one immutable predecessor descriptor, without recursive validation.

    Reuse eligibility is mechanical for runtime/data identities. Scientific app
    behavior remains an explicit classification backed by its test evidence.
    """
    if kind in ("app", "science"):
        require(base_release, "App/science releases require --base-release from a verified split deployment")
    if not base_release:
        require(expected_hash is None, "Unexpected baseline hash without a baseline release")
        return None
    root = artifact_path(project, "releases", base_release)
    verify_descriptor_files(root)
    path = regular(root / "deployment.json")
    descriptor = load_json(path)
    require(descriptor.get("schema") == "islet.deployment/1" and descriptor.get("id") == base_release,
            "Baseline must be a split deployment; migrate a legacy deployment with change-kind runtime")
    identity = digest(path)
    if expected_hash is not None:
        require(expected_hash == identity, "Baseline deployment descriptor changed")
    if kind in ("app", "science"):
        require(descriptor.get("runtime_manifest_sha256") == runtime_manifest and
                descriptor.get("data_receipt_sha256") == data_receipt,
                "App/science reuse requires the baseline runtime and data; use change-kind data or runtime with scientific validation")
    return identity


def assemble(project, name, app, data, runtime, validation_path, kind, base_release=None):
    resolved = resolve_candidate(project, app, data, runtime)
    validation = load_json(regular(validation_path))
    check_validation(validation, resolved, kind)
    runtime_manifest = digest(artifact_path(project, "runtimes", runtime) / "MANIFEST.sha256")
    baseline_hash = check_baseline(project, kind, runtime_manifest,
                                   resolved["DATA_RECEIPT_SHA256"], base_release)
    descriptor = {"schema": "islet.deployment/1", "id": name, "app": app, "data": data, "runtime": runtime,
                  "app_manifest_sha256": resolved["APP_MANIFEST_SHA256"],
                  "data_receipt_sha256": resolved["DATA_RECEIPT_SHA256"],
                  "runtime_manifest_sha256": runtime_manifest,
                  "change_kind": kind, "validation_sha256": digest(validation_path)}
    if base_release is not None:
        descriptor.update(base_release=base_release, base_descriptor_sha256=baseline_hash)

    def build(stage):
        shutil.copyfile(validation_path, stage / "validation.json")
        write_json(stage / "deployment.json", descriptor)
    return publish(project, "releases", name, build)


def resolve_release(release_dir):
    require(release_dir.parent.name == "releases", "Descriptor must be under PROJECT/releases/ID")
    project = project_path(release_dir.parent.parent)
    require(artifact_path(project, "releases", release_dir.name) == release_dir, "Invalid release path")
    verify_descriptor_files(release_dir)
    descriptor = load_json(regular(release_dir / "deployment.json"))
    require(descriptor.get("schema") == "islet.deployment/1" and descriptor.get("id") == release_dir.name,
            "Invalid deployment descriptor")
    resolved = resolve_candidate(project, descriptor["app"], descriptor["data"], descriptor["runtime"])
    require(descriptor["app_manifest_sha256"] == resolved["APP_MANIFEST_SHA256"] and
            descriptor["data_receipt_sha256"] == resolved["DATA_RECEIPT_SHA256"] and
            descriptor["runtime_manifest_sha256"] == digest(
                artifact_path(project, "runtimes", descriptor["runtime"]) / "MANIFEST.sha256"),
            "Deployment artifact identity changed")
    validation_path = release_dir / "validation.json"
    require(digest(validation_path) == descriptor["validation_sha256"], "Validation evidence changed")
    check_validation(load_json(validation_path), resolved, descriptor["change_kind"])
    if "base_release" in descriptor:
        require(re.fullmatch(HASH_RE, descriptor.get("base_descriptor_sha256", "")),
                "Missing baseline descriptor identity")
    check_baseline(project, descriptor["change_kind"], descriptor["runtime_manifest_sha256"],
                   descriptor["data_receipt_sha256"], descriptor.get("base_release"),
                   descriptor.get("base_descriptor_sha256"))
    resolved["DESCRIPTOR_SHA256"] = digest(release_dir / "deployment.json")
    return resolved


def emit(values):
    for key, value in values.items():
        require(not any(ord(c) < 32 for c in str(value)), "Unsafe resolved value")
        print(f"{key}={value}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command in ("package-app", "seal-data", "register-data", "register-runtime", "assemble", "resolve-candidate"):
        sub = commands.add_parser(command)
        sub.add_argument("--project-root", required=True, type=project_path)
        if command in ("package-app", "assemble", "resolve-candidate"):
            sub.add_argument("--app", required=True, type=identifier)
        if command in ("register-data", "assemble", "resolve-candidate"):
            sub.add_argument("--data", required=True, type=identifier)
        if command in ("register-runtime", "assemble", "resolve-candidate"):
            sub.add_argument("--runtime", required=True, type=identifier)
        if command == "package-app":
            sub.add_argument("--revision", required=True)
            sub.add_argument("--repo", type=Path, default=SCRIPT_DIR.parent)
        if command in ("register-data", "seal-data"):
            sub.add_argument("--source-release", required=True, type=identifier)
        if command == "seal-data":
            sub.add_argument("--data-schema", required=True, choices=("1", "5"))
            sub.add_argument("--expected-images", required=True, type=int)
        if command == "register-runtime":
            sub.add_argument("--scan-dir", required=True, type=Path)
            sub.add_argument("--runtime-input-sha256", required=True)
        if command == "assemble":
            sub.add_argument("--release", required=True, type=identifier)
            sub.add_argument("--base-release", type=identifier)
            sub.add_argument("--validation", required=True, type=Path)
            sub.add_argument("--change-kind", required=True, choices=("app", "science", "data", "runtime"))
    app_check = commands.add_parser("verify-app")
    app_check.add_argument("--app-dir", required=True, type=Path)
    app_check.add_argument("--expected-manifest", required=True)
    for command in ("verify", "resolve"):
        sub = commands.add_parser(command)
        sub.add_argument("--release-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "package-app":
        result = package_app(args.project_root, args.app, args.repo, args.revision)
    elif args.command == "seal-data":
        result = seal_data(args.project_root, args.source_release, args.data_schema, args.expected_images)
    elif args.command == "register-data":
        result = register_data(args.project_root, args.data, args.source_release)
    elif args.command == "register-runtime":
        result = register_runtime(args.project_root, args.runtime, args.scan_dir, args.runtime_input_sha256)
    elif args.command == "assemble":
        result = assemble(args.project_root, args.release, args.app, args.data, args.runtime, args.validation, args.change_kind, args.base_release)
    elif args.command == "verify-app":
        result = verify_app_identity(args.app_dir.absolute(), args.expected_manifest)
    elif args.command == "resolve-candidate":
        emit(resolve_candidate(args.project_root, args.app, args.data, args.runtime))
        return
    else:
        resolved = resolve_release(args.release_dir.absolute())
        if args.command == "resolve":
            emit(resolved)
        else:
            print(f"Verified deployment {args.release_dir.name}: app {resolved['APP_COMMIT']}, data {resolved['DATA_ID']}, runtime {resolved['IMAGE_REF']}")
        return
    print(result)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError, AttributeError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f"ERROR: {exc}")
