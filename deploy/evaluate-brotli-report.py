#!/usr/bin/env python3
"""Separate the measured decoder result from the dated release decision.

No general CVE suppression: only the known DuckDB 1.5.2 decoder result may be
deferred, and only through the patch date in the retained owner instruction.
"""
import datetime
import json
import pathlib
import re
import sys

DUCKDB_ERROR = "DuckDB's bundled Brotli fails the upstream block-switch regression (CVE-2026-111497)."
DEFERRAL_ID = "duckdb-cve-2026-111497-until-20261110"
DEFERRED_RESULT = "pass-with-duckdb-remediation-window"


def check_window(window, today=None):
    if (window.get("schema") != "islet.duckdb_remediation_window/1"
            or window.get("id") != DEFERRAL_ID
            or window.get("cve") != "CVE-2026-111497"
            or window.get("component") != "DuckDB 1.5.2 bundled Brotli"
            or window.get("duckdb_version") != "1.5.2"
            or window.get("patch_by") != "2026-11-10"
            or window.get("deadline_timezone") != "UTC"
            or window.get("authorization_date") != "2026-11-02"
            or not window.get("authorization")
            or window.get("institutional_approval_claimed") is not False):
        raise ValueError("Missing or unsupported DuckDB remediation-window authorization.")
    today = today or datetime.datetime.now(datetime.timezone.utc).date()
    if not datetime.date(2026, 11, 2) <= today <= datetime.date.fromisoformat(window["patch_by"]):
        raise ValueError("DuckDB remediation window is not current; patch deadline is 2026-11-10 UTC.")


def evaluate(report, window=None, today=None):
    duckdb = report.get("duckdb", {})
    if (report.get("schema") != "islet.brotli_verification/2"
            or report.get("cve") != "CVE-2026-111497"
            or duckdb.get("version") != "1.5.2"
            or duckdb.get("result") != 0
            or not re.fullmatch(r"[0-9a-f]{64}", str(duckdb.get("library_sha256", "")))):
        raise ValueError("A complete current DuckDB binary report is required.")
    if report.get("result") == "pass" and report.get("errors") == [] and duckdb.get("error_code") == -17:
        return "pass-no-high-or-critical"
    if (report.get("result") != "fail" or report.get("errors") != [DUCKDB_ERROR]
            or duckdb.get("error_code") != -9 or duckdb.get("expected_error_code") != -17):
        raise ValueError("Brotli checks failed outside the narrowly deferred DuckDB finding.")
    if window is None:
        raise ValueError("DuckDB is unresolved and no remediation-window record was supplied.")
    check_window(window, today)
    return DEFERRED_RESULT


if __name__ == "__main__":
    try:
        if len(sys.argv) not in (2, 3):
            raise ValueError("Usage: evaluate-brotli-report.py REPORT [DUCKDB_REMEDIATION_WINDOW]")
        report = json.loads(pathlib.Path(sys.argv[1]).read_text())
        window = json.loads(pathlib.Path(sys.argv[2]).read_text()) if len(sys.argv) == 3 else None
        print(evaluate(report, window))
    except (ValueError, OSError, TypeError, AttributeError) as exc:
        raise SystemExit(f"ERROR: {exc}")
