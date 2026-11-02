#!/usr/bin/env python3
"""Recheck bounded security evidence on every full or quick release verification."""
import hashlib
import json
import pathlib
import re
import runpy
import sys

script_dir = pathlib.Path(__file__).resolve().parent
decision = runpy.run_path(str(script_dir / "evaluate-brotli-report.py"))


def fields(path):
    result = {}
    for line in path.read_text().splitlines():
        key, sep, value = line.partition("=")
        if not sep or key in result:
            raise ValueError(f"Malformed or duplicate field in {path.name}.")
        result[key] = value
    return result


def verify(root):
    entries = {}
    for line in (root / "MANIFEST.sha256").read_text().splitlines():
        match = re.fullmatch(r"([0-9a-f]{64}) [ *](.+)", line)
        if match:
            entries.setdefault(match[2], []).append(match[1])

    def bound_file(relative, expected=None):
        path = root / relative
        hashes = entries.get(relative, [])
        if (len(hashes) != 1 or not path.is_file() or path.is_symlink()
                or (root / "attestations").is_symlink()):
            raise ValueError(f"Missing, linked or unbound security evidence: {relative}")
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1 << 20), b""):
                digest.update(chunk)
        if digest.hexdigest() != hashes[0] or (expected is not None and expected != hashes[0]):
            raise ValueError(f"Security evidence hash mismatch: {relative}")
        return path

    meta = fields(bound_file("release.env"))
    scan_paths = {name for name in entries if name.startswith("attestations/") and name.endswith(".scan.txt")}
    scan_paths.update(p.relative_to(root).as_posix() for p in (root / "attestations").glob("*.scan.txt"))
    # Explicit historical compatibility: releases predating scan evidence have
    # no security decision/deferral artifacts. New manifests always require them.
    if not scan_paths and not meta.get("SECURITY_RESULT"):
        if meta.get("DUCKDB_PATCH_BY") or any("duckdb-deferral" in name or "brotli.json" in name for name in entries):
            raise ValueError("A release cannot discard its DuckDB decision and become legacy evidence.")
        return
    if len(scan_paths) != 1:
        raise ValueError("Exactly one manifest-bound image scan is required.")
    scan_path = next(iter(scan_paths))
    scan = fields(bound_file(scan_path))
    image_id = scan.get("IMAGE_ID", "")
    if not re.fullmatch(r"sha256:[0-9a-f]{64}", image_id) or scan.get("IMAGE_REF") != image_id:
        raise ValueError("Scan must identify its immutable image.")
    prefix = "attestations/" + image_id.replace(":", "_")
    if scan_path != prefix + ".scan.txt":
        raise ValueError("Scan filename/image mismatch.")
    result = scan.get("RESULT")
    if result not in ("pass-no-high-or-critical", decision["DEFERRED_RESULT"]):
        raise ValueError("Security scan did not permit release.")
    if meta.get("SECURITY_RESULT", "pass-no-high-or-critical") != result:
        raise ValueError("Release and scan security decisions differ.")
    for key, suffix in [("SBOM_SHA256", ".cyclonedx.json"), ("RAW_SBOM_SHA256", ".syft.cyclonedx.json"),
                        ("SBOM_CORRECTIONS_SHA256", ".sbom-corrections.json"), ("REPORT_SHA256", ".grype.json"),
                        ("R_PACKAGES_SHA256", ".R-packages.tsv"), ("PYTHON_PACKAGES_SHA256", ".python-packages.txt")]:
        if not re.fullmatch(r"[0-9a-f]{64}", scan.get(key, "")):
            raise ValueError(f"Missing scan evidence hash: {key}")
        bound_file(prefix + suffix, scan[key])
    grype = json.loads((root / (prefix + ".grype.json")).read_text())
    if not isinstance(grype.get("matches"), list) or any(
            m["vulnerability"]["severity"].lower() in ("high", "critical") for m in grype["matches"]):
        raise ValueError("High/Critical package finding or incomplete Grype report.")
    published = fields(bound_file("attestations/published-image.txt"))
    if (published.get("GIT_COMMIT") != meta.get("GIT_COMMIT")
            or published.get("LOCAL_IMAGE_ID", published.get("IMAGE_ID")) != image_id
            or published.get("REGISTRY_REF", published.get("LOCAL_IMAGE_ID")) != meta.get("IMAGE_REF")):
        raise ValueError("Published image, scan and release identities differ.")
    report_hash = scan.get("BROTLI_REPORT_SHA256")
    report_name = prefix + ".brotli.json"
    if not report_hash and report_name not in entries and not (root / report_name).exists():
        if meta.get("SECURITY_RESULT") or result != "pass-no-high-or-critical":
            raise ValueError("Current security decision requires actual decoder evidence.")
        return  # historical package-only receipt; never a dated DuckDB deferral
    if not re.fullmatch(r"[0-9a-f]{64}", report_hash or ""):
        raise ValueError("Missing decoder report hash.")
    report = json.loads(bound_file(report_name, report_hash).read_text())
    window = None
    if result == decision["DEFERRED_RESULT"]:
        window_hash = scan.get("DUCKDB_DEFERRAL_SHA256", "")
        if not re.fullmatch(r"[0-9a-f]{64}", window_hash):
            raise ValueError("Missing remediation-window hash.")
        window = json.loads(bound_file(prefix + ".duckdb-deferral.json", window_hash).read_text())
        if window != json.loads((script_dir / "duckdb-remediation-window.json").read_text()):
            raise ValueError("DuckDB authorization differs from the recorded owner direction.")
        if (scan.get("PACKAGE_RESULT") != "pass-no-high-or-critical"
                or scan.get("DUCKDB_PATCH_BY") != window["patch_by"]
                or meta.get("DUCKDB_PATCH_BY") != window["patch_by"]):
            raise ValueError("Missing or inconsistent DuckDB deadline/package result.")
    if decision["evaluate"](report, window) != result:
        raise ValueError("Measured decoder and release decision disagree.")


if __name__ == "__main__":
    try:
        if len(sys.argv) != 2:
            raise ValueError("Usage: verify-release-security.py RELEASE_DIR")
        verify(pathlib.Path(sys.argv[1]).resolve())
    except (ValueError, OSError, TypeError, KeyError, AttributeError) as exc:
        raise SystemExit(f"ERROR: {exc}")
