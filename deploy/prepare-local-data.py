#!/usr/bin/env python3
"""Create a data-keyed local cache once; reuse unchanged, verified copies.

The receipt is an operational integrity check, not protection against an owner
who can rewrite both files and receipt. Never repair a published cache in place.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

TREES = ("geometry", "cells", "spatial")


def digest(path):
    # PubApps currently uses Python 3.9; hashlib.file_digest needs Python 3.11.
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def inventory(root):
    result = {}
    for tree in TREES:
        base = root / "data" / tree
        if not base.is_dir() or base.is_symlink():
            raise ValueError(f"Missing or linked cache tree: {base}")
        for path in sorted(base.rglob("*")):
            if path.is_symlink():
                raise ValueError(f"Local cache may not contain symlinks: {path}")
            if path.is_file():
                stat = path.stat()
                result[path.relative_to(root).as_posix()] = [stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]
    return result


def prepare(source, target):
    manifest = source / "MANIFEST.sha256"
    identity = digest(manifest)
    if target.exists():
        if target.is_symlink():
            raise ValueError("Local data cache must not be a symlink")
        receipt = json.loads((target / ".verified.json").read_text())
        if receipt.get("manifest_sha256") != identity or receipt.get("inventory") != inventory(target):
            raise ValueError("Published local cache changed; use a new data ID, never repair an active cache in place")
        print(f"Reusing verified local data: {target}")
        return
    expected = {}
    for line in manifest.read_text().splitlines():
        match = re.fullmatch(r"([0-9a-f]{64}) [ *](data/(?:geometry|cells|spatial)/.+)", line)
        if match:
            name = match[2]
            if ".." in Path(name).parts or name in expected:
                raise ValueError("Unsafe or duplicate local-data manifest path")
            expected[name] = match[1]
    if not expected:
        raise ValueError("Manifest contains no local-data entries")
    target.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".cache-", dir=target.parent))
    try:
        for tree in TREES:
            destination = stage / "data" / tree
            destination.mkdir(parents=True)
            subprocess.run(["rsync", "-a", "--", str(source / "data" / tree) + "/", str(destination) + "/"], check=True)
        actual = inventory(stage)
        if set(actual) != set(expected) or any(digest(stage / name) != sha for name, sha in expected.items()):
            raise ValueError("Local copies do not match the verified data manifest")
        (stage / ".verified.json").write_text(json.dumps({"manifest_sha256": identity, "inventory": actual}, sort_keys=True) + "\n")
        os.rename(stage, target)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    print(f"Verified local data cache: {target}")


if __name__ == "__main__":
    try:
        if len(sys.argv) != 3:
            raise ValueError("Usage: prepare-local-data.py SOURCE_RELEASE_DIR CACHE_DIR")
        prepare(Path(sys.argv[1]).resolve(strict=True), Path(sys.argv[2]).absolute())
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f"ERROR: {exc}")
