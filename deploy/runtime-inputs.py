#!/usr/bin/env python3
"""Fingerprint the committed files consumed by the software-only runtime build.

Usage: runtime-inputs.py REPO_DIR [REVISION]
       runtime-inputs.py --base-image REPO_DIR [REVISION]
With REVISION, read Git objects without consulting or changing the worktree.
Without REVISION, read an isolated archive of that revision from disk.
"""
import hashlib
import pathlib
import re
import shlex
import stat
import subprocess
import sys


def base_image_requirement(repo, revision=None):
    """Read the selected source's one digest-pinned base requirement.

    This is deliberately separate from the v1 runtime file fingerprint so
    existing split runtimes retain their fingerprint identity. Both identities
    must match when selecting a runtime for an app bundle.
    """
    repo = pathlib.Path(repo)
    relative = "deploy/base-image.lock"
    if revision is not None:
        revision = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "--verify", revision + "^{commit}"],
            text=True,
        ).strip()
        entry = subprocess.check_output(
            ["git", "-C", str(repo), "ls-tree", revision, "--", relative], text=True
        ).strip()
        if entry.split(" ", 1)[0] not in ("100644", "100755"):
            raise ValueError("Base-image lock must be a committed regular file")
        contents = subprocess.check_output(
            ["git", "-C", str(repo), "show", f"{revision}:{relative}"], text=True
        )
    else:
        path = repo / relative
        if not path.is_file() or path.is_symlink():
            raise ValueError("Base-image lock must be a regular file")
        contents = path.read_text()
    entries = [line.strip() for line in contents.splitlines()
               if line.strip() and not line.lstrip().startswith("#")]
    if len(entries) != 1 or not re.fullmatch(r"[A-Za-z0-9._:/-]+@sha256:[0-9a-f]{64}", entries[0]):
        raise ValueError("Base-image lock must contain exactly one digest-pinned image reference")
    return entries[0]


def runtime_fingerprint(repo, revision=None):
    repo = pathlib.Path(repo)
    if revision is not None:
        revision = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "--verify", revision + "^{commit}"],
            text=True,
        ).strip()

    def read(relative):
        if revision is not None:
            entry = subprocess.check_output(
                ["git", "-C", str(repo), "ls-tree", revision, "--", relative], text=True
            ).strip()
            mode = entry.split(" ", 1)[0]
            if mode not in ("100644", "100755"):
                raise ValueError(f"Runtime input must be a committed regular file: {relative}")
            data = subprocess.check_output(["git", "-C", str(repo), "show", f"{revision}:{relative}"])
            return mode, data
        path = repo / relative
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"Runtime input must be a regular file: {relative}")
        mode = "100755" if path.stat().st_mode & stat.S_IXUSR else "100644"
        return mode, path.read_bytes()

    containerfile = read("deploy/Containerfile")[1].decode("utf-8")
    paths = {"deploy/Containerfile", "deploy/containerignore"}
    for line in re.sub(r"\\\r?\n", " ", containerfile).splitlines():
        tokens = shlex.split(line, comments=True)
        if not tokens or tokens[0].upper() not in ("COPY", "ADD"):
            continue
        if tokens[0].upper() == "ADD":
            raise ValueError("Runtime fingerprint does not permit ADD; use explicit local COPY inputs")
        if any(token.startswith("--from=") for token in tokens[1:]):
            continue
        operands = [token for token in tokens[1:] if not token.startswith("--")]
        if len(operands) < 2:
            raise ValueError(f"Unsupported COPY syntax: {line}")
        for relative in operands[:-1]:
            if not re.fullmatch(r"[A-Za-z0-9_.+/-]+", relative):
                raise ValueError(f"Use explicit runtime COPY paths: {relative}")
            path = pathlib.PurePosixPath(relative)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError(f"Unsafe runtime input path: {relative}")
            paths.add(relative)

    digest = hashlib.sha256(b"islet-runtime-inputs-v1\0")
    for relative in sorted(paths):
        mode, data = read(relative)
        digest.update(relative.encode("utf-8") + b"\0" + mode.encode("ascii") + b"\0")
        digest.update(hashlib.sha256(data).digest())
    return digest.hexdigest()


if __name__ == "__main__":
    try:
        args = sys.argv[1:]
        operation = runtime_fingerprint
        if args and args[0] == "--base-image":
            args = args[1:]
            operation = base_image_requirement
        if len(args) not in (1, 2):
            raise ValueError("Usage: runtime-inputs.py [--base-image] REPO_DIR [REVISION]")
        print(operation(args[0], args[1] if len(args) == 2 else None))
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f"ERROR: {exc}")
