#!/usr/bin/env python3
"""Check the image's Brotli copies, including Pillow libraries absent from Syft.

Run with the image's application Python. This gate deliberately requires our
source-built Pillow/system FreeType layout, not a wheel's vendored libraries.
Ubuntu backports and Node's patched release retain older Brotli version strings;
check their distributor package versions instead of rejecting those strings.
"""

import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys


DUCKDB_SOURCE_RECEIPT = pathlib.Path("/opt/islet-explorer-build/duckdb-brotli-source.json")
DUCKDB_SOURCE_PINS = {
    "schema": "islet.duckdb_brotli_source_patch/1",
    "cve": "CVE-2026-111497",
    "package": "duckdb",
    "package_version": "1.5.2",
    "result": "source_patch_verified",
    "source_sha256": "55f98a336513a94354f328791c8b94521b9f2e918a57d7e9fe6c0cb6b9a38db6",
    "patch_sha256": "d00033af815fcb214daa9630d9da11b39a191dfab17b0f8cff367157af5806ca",
    "upstream_commit": "0ed4268a4401abcb4472096340477cf4af7d1a2b",
    "upstream_1_1_1_backport_commit": "da9806362874b9bd4b552a3af360b5d6347beb12",
    "port_commit": "2b942af6319783a5cc3c21d76a8e560e191978b2",
}
DUCKDB_SOURCE_FILES = {
    "third_party/brotli/dec/decode.cpp": {
        "before_sha256": "27062f01a0f429f14ea0bc0587d9db83230e60f7725c7253817a6c3517eb08f4",
        "after_sha256": "fff654fe532c55f0f70503b6af1a0ad7166b33edb68ff89c2eb021353b40595d",
    },
    "third_party/brotli/include/brotli/decode.h": {
        "before_sha256": "a39110bec706e846f51c05a6374282c4cec04b099ead97e425f049840f2a02b8",
        "after_sha256": "781b0e6ece03bdc8800189983c9e636bd1ac550f5bab582c5c36960042a84155",
    },
}


def file_sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def unique_json_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON field: {key}")
        result[key] = value
    return result


def verify_duckdb_source(path, runtime_version, errors):
    """Validate the reviewed source receipt; the installed decoder is tested separately."""
    path = pathlib.Path(path)
    result = {"path": str(path), "result": "fail", "runtime_package_version": runtime_version}
    initial_errors = len(errors)
    if runtime_version != DUCKDB_SOURCE_PINS["package_version"]:
        errors.append(f"DuckDB {runtime_version} is not the pinned 1.5.2 source backport.")
    try:
        if path.is_symlink() or not path.is_file():
            raise ValueError("receipt must be a regular, non-symlink file")
        raw = path.read_bytes()
        result["receipt_sha256"] = hashlib.sha256(raw).hexdigest()
        receipt = json.loads(raw, object_pairs_hook=unique_json_object)
        if not isinstance(receipt, dict):
            raise ValueError("receipt must be a JSON object")
        result["receipt"] = receipt
        for key, expected in DUCKDB_SOURCE_PINS.items():
            if receipt.get(key) != expected:
                errors.append(f"DuckDB source receipt has an unverified {key}.")
        for key in ("patched_source_sha256", "helper_sha256"):
            if not isinstance(receipt.get(key), str) or not re.fullmatch(r"[0-9a-f]{64}", receipt[key]):
                errors.append(f"DuckDB source receipt is missing a valid {key}.")
        files = receipt.get("files")
        if (not isinstance(files, list) or len(files) != len(DUCKDB_SOURCE_FILES)
                or any(not isinstance(item, dict) for item in files)):
            raise ValueError("receipt must describe exactly the two reviewed decoder source files")
        paths = [item.get("path") for item in files]
        if any(not isinstance(name, str) for name in paths) or set(paths) != set(DUCKDB_SOURCE_FILES):
            raise ValueError("receipt decoder source inventory differs from the reviewed patch")
        for item in files:
            for key, expected in DUCKDB_SOURCE_FILES[item["path"]].items():
                if item.get(key) != expected:
                    errors.append(f"DuckDB source receipt has an unverified {item['path']} {key}.")
        verification = receipt.get("verification")
        if not isinstance(verification, dict):
            raise ValueError("receipt is missing the unchanged-source verification")
        for key in ("unrelated_member_payloads_and_metadata_unchanged", "member_order_and_mtimes_preserved"):
            if verification.get(key) is not True:
                errors.append(f"DuckDB source receipt did not verify {key}.")
    except (OSError, ValueError, TypeError) as exc:
        errors.append(f"DuckDB source provenance failed: {exc}")
    if len(errors) == initial_errors:
        result["result"] = "pass"
    return result


# Use the upstream regression vector, in an isolated subprocess with no data
# mounts/network (scan-image.sh). A version string cannot identify a backport.
# https://github.com/duckdb/duckdb/pull/27791
DUCKDB_PROBE = r'''
import ctypes as c, json, sys
c.CDLL(sys.argv[2], mode=c.RTLD_GLOBAL)
lib = c.CDLL(sys.argv[1])
u8 = c.c_uint8
ptr = c.POINTER(u8)
sizeptr = c.POINTER(c.c_size_t)
create = lib._ZN13duckdb_brotli27BrotliDecoderCreateInstanceEPFPvS0_mEPFvS0_S0_ES0_
create.argtypes = [c.c_void_p] * 3
create.restype = c.c_void_p
decode = lib._ZN13duckdb_brotli29BrotliDecoderDecompressStreamEPNS_24BrotliDecoderStateStructEPmPPKhS2_PPhS2_
decode.argtypes = [c.c_void_p, sizeptr, c.POINTER(ptr), sizeptr, c.POINTER(ptr), sizeptr]
decode.restype = c.c_int
error = lib._ZN13duckdb_brotli25BrotliDecoderGetErrorCodeEPKNS_24BrotliDecoderStateStructE
error.argtypes = [c.c_void_p]
error.restype = c.c_int
destroy = lib._ZN13duckdb_brotli28BrotliDecoderDestroyInstanceEPNS_24BrotliDecoderStateStructE
destroy.argtypes = [c.c_void_p]
destroy.restype = None
data = bytes.fromhex('200010414141fcffff07002228f00b60f7f41f000000')
source = (u8 * len(data)).from_buffer_copy(data)
target = (u8 * 65536)()
next_in, next_out = c.cast(source, ptr), c.cast(target, ptr)
available_in, available_out, total_out = c.c_size_t(len(data)), c.c_size_t(len(target)), c.c_size_t()
state = create(None, None, None)
if not state:
    raise RuntimeError('Brotli decoder allocation failed')
try:
    for _ in range(256):
        next_out, available_out = c.cast(target, ptr), c.c_size_t(len(target))
        result = decode(state, c.byref(available_in), c.byref(next_in), c.byref(available_out), c.byref(next_out), c.byref(total_out))
        if result != 3:  # BROTLI_DECODER_RESULT_NEEDS_MORE_OUTPUT
            break
    print(json.dumps({'result': result, 'error_code': error(state), 'expected_error_code': -17}))
finally:
    destroy(state)
'''


def command(*args):
    return subprocess.check_output(args, text=True).strip()


def verify():
    report = {"schema": "islet.brotli_verification/2", "cve": "CVE-2026-111497", "result": "fail", "errors": []}
    errors = report["errors"]
    os_release = pathlib.Path("/etc/os-release").read_text()
    if 'VERSION_ID="24.04"' not in os_release or '\nID=ubuntu\n' not in os_release:
        errors.append("Review Brotli's fixed package version for this base OS.")
    packages = {}
    for package in ("libbrotli1", "libbrotli-dev"):
        version = command("dpkg-query", "-W", "-f=${Version}", package)
        packages[package] = version
        if subprocess.call(["dpkg", "--compare-versions", version, "ge", "1.1.0-2ubuntu0.1"]):
            errors.append(f"{package} {version} predates Ubuntu's Brotli fix.")
    report["ubuntu_packages"] = packages

    node = json.loads(command("/opt/shiny-server/ext/node/bin/node", "-p", "JSON.stringify(process.versions)"))
    report["node"] = {key: node[key] for key in ("node", "brotli")}
    version = tuple(map(int, node["node"].split(".")))
    minimum = {22: (22, 23, 4), 24: (24, 22, 1), 26: (26, 12, 1)}.get(version[0])
    if minimum is None or version < minimum:
        errors.append(f"Node {node['node']} is not a verified patched release.")

    from PIL import __version__ as pillow_version, _imagingft

    report["pillow_version"] = pillow_version
    vendored = sorted(str(p) for p in pathlib.Path(sys.prefix).rglob("*brotli*") if ".so" in p.name or p.suffix == ".a")
    report["vendored_brotli"] = vendored
    if vendored:
        errors.append("Application Python contains a vendored Brotli library; rebuild Pillow from source.")
    linkage = command("ldd", _imagingft.__file__)
    report["pillow_freetype_linkage"] = linkage.splitlines()
    managed = {pathlib.Path(p).resolve() for p in command("dpkg-query", "-L", "libbrotli1").splitlines()}
    decoders = re.findall(r"libbrotlidec[^\s]* => (\S+)", linkage)
    if len(decoders) != 1 or pathlib.Path(decoders[0]).resolve() not in managed:
        errors.append("Pillow must link transitively to the Ubuntu-managed Brotli decoder.")
    if "not found" in linkage:
        errors.append("Pillow has an unresolved shared-library dependency.")

    # DuckDB statically vendors Brotli in its Parquet reader; it is absent from
    # package-only scans and cannot be repaired by upgrading libbrotli1.
    r_home = pathlib.Path(command("R", "RHOME"))
    duckdb_version, duckdb_library = command(
        "Rscript", "--vanilla", "-e",
        'cat(as.character(packageVersion("duckdb")), system.file("libs", "duckdb.so", package="duckdb"), sep="\\t")',
    ).split("\t")
    report["duckdb"] = {
        "version": duckdb_version,
        "library": duckdb_library,
        "library_sha256": file_sha256(duckdb_library),
        "source_provenance": verify_duckdb_source(DUCKDB_SOURCE_RECEIPT, duckdb_version, errors),
    }
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = str(r_home / "lib") + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
    probe = subprocess.run(
        [sys.executable, "-c", DUCKDB_PROBE, duckdb_library, str(r_home / "lib/libR.so")],
        env=env, text=True, capture_output=True, timeout=15, check=True,
    )
    report["duckdb"].update(json.loads(probe.stdout))
    if report["duckdb"]["result"] != 0 or report["duckdb"]["error_code"] != -17:
        errors.append("DuckDB's bundled Brotli fails the upstream block-switch regression (CVE-2026-111497).")
    if not errors:
        report["result"] = "pass"
    return report


if __name__ == "__main__":
    try:
        result = verify()
    except Exception as exc:
        result = {"schema": "islet.brotli_verification/2", "cve": "CVE-2026-111497", "result": "fail", "errors": [str(exc)]}
    print(json.dumps(result, indent=2))
    sys.exit(0 if result["result"] == "pass" else 1)
