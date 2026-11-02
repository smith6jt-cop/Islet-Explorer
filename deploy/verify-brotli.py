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


def file_sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


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
