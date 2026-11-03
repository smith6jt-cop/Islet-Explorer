#!/usr/bin/env python3
"""Backport the reviewed Brotli fix into the exact CRAN DuckDB 1.5.2 source.

The DuckDB port is an open, unmerged upstream PR, not an official CRAN release.
Preserve the R API/version and every unrelated source member. Patch the nested
archive that configure actually extracts, update its package MD5 entry, and
emit a receipt. Installation must rebuild from this archive without a prebuilt
object cache; the installed decoder regression remains a separate required gate.
"""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import io
import json
import lzma
import pathlib
import re
import struct
import subprocess
import tarfile
import tempfile


SOURCE_SHA256 = "55f98a336513a94354f328791c8b94521b9f2e918a57d7e9fe6c0cb6b9a38db6"
PATCH_SHA256 = "d00033af815fcb214daa9630d9da11b39a191dfab17b0f8cff367157af5806ca"
UPSTREAM_COMMIT = "0ed4268a4401abcb4472096340477cf4af7d1a2b"
BACKPORT_COMMIT = "da9806362874b9bd4b552a3af360b5d6347beb12"
PORT_COMMIT = "2b942af6319783a5cc3c21d76a8e560e191978b2"
NESTED = "duckdb/src/duckdb.tar.xz"
MD5_FILE = "duckdb/MD5"
FILES = {
    "third_party/brotli/dec/decode.cpp": (
        "27062f01a0f429f14ea0bc0587d9db83230e60f7725c7253817a6c3517eb08f4",
        "fff654fe532c55f0f70503b6af1a0ad7166b33edb68ff89c2eb021353b40595d",
    ),
    "third_party/brotli/include/brotli/decode.h": (
        "a39110bec706e846f51c05a6374282c4cec04b099ead97e425f049840f2a02b8",
        "781b0e6ece03bdc8800189983c9e636bd1ac550f5bab582c5c36960042a84155",
    ),
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def members(data: bytes) -> list[tuple[tarfile.TarInfo, bytes | None]]:
    """Read a hash-verified archive without extracting it onto the filesystem."""
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:*") as archive:
        result = []
        seen = set()
        for info in archive:
            require(info.name not in seen, f"Duplicate archive member: {info.name}")
            seen.add(info.name)
            require(info.isfile() or info.isdir(), f"Unexpected archive member type: {info.name}")
            require(not info.name.startswith("/") and ".." not in pathlib.PurePosixPath(info.name).parts,
                    f"Unsafe archive member: {info.name}")
            result.append((info, archive.extractfile(info).read() if info.isfile() else None))
        return result


def repack(records: list, replacements: dict[str, bytes]) -> bytes:
    """Keep member order, ownership, modes and original mtimes; change payloads only."""
    require(set(replacements) <= {info.name for info, _ in records}, "Replacement member is absent")
    target = io.BytesIO()
    with tarfile.open(fileobj=target, mode="w", format=tarfile.GNU_FORMAT) as archive:
        for original, data in records:
            info = copy.copy(original)
            info.pax_headers = dict(original.pax_headers)
            if info.name in replacements:
                data = replacements[info.name]
                info.size = len(data)
                if "size" in info.pax_headers:
                    info.pax_headers["size"] = str(info.size)
            archive.addfile(info, None if data is None else io.BytesIO(data))
    return target.getvalue()


def verify_repack(before: list, packed: bytes, expected: dict[str, bytes]) -> int:
    after = members(packed)
    require([i.name for i, _ in before] == [i.name for i, _ in after], "Repack changed member inventory/order")
    fields = ("name", "mode", "uid", "gid", "mtime", "type", "linkname", "uname", "gname", "devmajor", "devminor")
    for (a, raw), (b, new) in zip(before, after):
        require(all(getattr(a, k) == getattr(b, k) for k in fields), f"Repack changed metadata: {a.name}")
        require(new == expected.get(a.name, raw), f"Unexpected payload change: {a.name}")
    return len(after)


def build(source: bytes, patch: pathlib.Path) -> tuple[bytes, dict]:
    require(sha256(source) == SOURCE_SHA256, "Unknown or modified CRAN source; refusing to patch")
    patch_bytes = patch.read_bytes()
    require(sha256(patch_bytes) == PATCH_SHA256, "Patch SHA-256 mismatch; refusing to patch")
    outer = members(source)
    content = {info.name: raw for info, raw in outer}
    nested_original = content[NESTED]
    inner = members(nested_original)
    nested_content = {info.name: raw for info, raw in inner}
    replacements, file_receipts = {}, []
    with tempfile.TemporaryDirectory(prefix="islet-duckdb-brotli-") as temporary:
        stage = pathlib.Path(temporary)
        for name, (before, _) in FILES.items():
            raw = nested_content[f"duckdb/{name}"]
            require(sha256(raw) == before, f"Unexpected decoder source: {name}")
            destination = stage / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(raw)
        # Apply the verified bytes, not a potentially replaced external path.
        patch_result = subprocess.run(
            ["patch", "--batch", "--forward", "--fuzz=0", "-p1"],
            input=patch_bytes, cwd=stage, capture_output=True, check=True,
        )
        patch_log = (patch_result.stdout + patch_result.stderr).decode("utf-8", errors="replace")
        require(not re.search(r"offset|fuzz|FAILED|Reversed", patch_log, re.IGNORECASE),
                f"Patch did not apply exactly: {patch_log}")
        for name, (before, after) in FILES.items():
            raw = (stage / name).read_bytes()
            require(sha256(raw) == after, f"Unexpected patched decoder: {name}")
            replacements[f"duckdb/{name}"] = raw
            file_receipts.append({"path": name, "before_sha256": before, "after_sha256": after})
    inner_tar = repack(inner, replacements)
    inner_count = verify_repack(inner, inner_tar, replacements)
    nested_new = lzma.compress(inner_tar, format=lzma.FORMAT_XZ, check=lzma.CHECK_CRC64, preset=9)
    old_md5 = hashlib.md5(nested_original).hexdigest()
    new_md5 = hashlib.md5(nested_new).hexdigest()
    md5_text = content[MD5_FILE].decode("utf-8")
    old_entry = f"{old_md5} *src/duckdb.tar.xz"
    require(md5_text.splitlines().count(old_entry) == 1, "Package MD5 does not bind the original nested archive")
    new_entry = f"{new_md5} *src/duckdb.tar.xz"
    md5_new = md5_text.replace(old_entry, new_entry).encode("utf-8")
    outer_replacements = {NESTED: nested_new, MD5_FILE: md5_new}
    outer_tar = repack(outer, outer_replacements)
    outer_count = verify_repack(outer, outer_tar, outer_replacements)
    require(source[:3] == b"\x1f\x8b\x08", "Expected the pinned gzip source archive")
    original_mtime = struct.unpack("<I", source[4:8])[0]
    destination = io.BytesIO()
    with gzip.GzipFile(fileobj=destination, filename="", mode="wb", compresslevel=9, mtime=original_mtime) as archive:
        archive.write(outer_tar)
    patched = destination.getvalue()
    receipt = {
        "schema": "islet.duckdb_brotli_source_patch/1", "cve": "CVE-2026-111497",
        "package": "duckdb", "package_version": "1.5.2", "result": "source_patch_verified",
        "source_url": "https://cran.r-project.org/src/contrib/Archive/duckdb/duckdb_1.5.2.tar.gz",
        "source_sha256": SOURCE_SHA256, "patched_source_sha256": sha256(patched),
        "patch_sha256": PATCH_SHA256, "upstream_commit": UPSTREAM_COMMIT,
        "upstream_1_1_1_backport_commit": BACKPORT_COMMIT, "port_commit": PORT_COMMIT,
        "patch_url": f"https://github.com/duckdb/duckdb/commit/{PORT_COMMIT}.patch",
        "port_pr": "https://github.com/duckdb/duckdb/pull/27791",
        "port_status_at_review": "open_unmerged_2026-11-02",
        "files": file_receipts, "patch_application": patch_log.strip(),
        "nested_archive": {"path": NESTED, "before_sha256": sha256(nested_original),
                           "after_sha256": sha256(nested_new), "before_md5": old_md5, "after_md5": new_md5},
        "package_md5": {"before_sha256": sha256(content[MD5_FILE]), "after_sha256": sha256(md5_new),
                        "updated_entry": "src/duckdb.tar.xz"},
        "verification": {"inner_members": inner_count, "outer_members": outer_count,
                         "unrelated_member_payloads_and_metadata_unchanged": True,
                         "member_order_and_mtimes_preserved": True, "gzip_mtime": original_mtime},
        "helper_sha256": sha256(pathlib.Path(__file__).read_bytes()),
        "limits": "Source backport only; requires a fresh compiled decoder regression and valid-data checks.",
    }
    return patched, receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--receipt", type=pathlib.Path, required=True)
    parser.add_argument("--patch", type=pathlib.Path,
                        default=pathlib.Path(__file__).parent / "patches/duckdb-brotli-CVE-2026-111497.patch")
    args = parser.parse_args()
    require(len({p.resolve() for p in (args.source, args.output, args.receipt, args.patch)}) == 4,
            "Source, patch, output and receipt paths must be distinct")
    for path in (args.output, args.receipt):
        require(not path.exists() and not path.is_symlink(), f"Output must be fresh: {path}")
        require(path.parent.is_dir(), f"Output parent does not exist: {path.parent}")
    patched, receipt = build(args.source.read_bytes(), args.patch)
    receipt.update({"source_path": str(args.source.resolve()), "output_path": str(args.output.resolve()),
                    "patch_path": str(args.patch.resolve())})
    with args.output.open("xb") as handle:
        handle.write(patched)
    with args.receipt.open("x") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"result": receipt["result"], "patched_source_sha256": sha256(patched),
                      "receipt": str(args.receipt)}))


if __name__ == "__main__":
    main()
