#!/usr/bin/env python3
"""Synthetic HTTP transport tests; the production CLI still requires HTTPS."""
import base64
from contextlib import contextmanager, redirect_stderr
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import http.client
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import threading
import unittest

SPEC = importlib.util.spec_from_file_location("reviewer_access", Path(__file__).resolve().parents[1] / "verify-reviewer-access.py")
verifier = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = verifier
SPEC.loader.exec_module(verifier)
PREFIX = b"II*\x00" + bytes(range(12))
PREFIX_SHA = hashlib.sha256(PREFIX).hexdigest()
IMAGE = "local_images/PRIVATE_FILENAME_CANARY.ome.tiff"
SIDECAR = "local_images/PRIVATE_FILENAME_CANARY.ome.xml"
AUTH = verifier.Credentials("PRIVATE_USERNAME_CANARY", "PRIVATE_PASSWORD_CANARY")


class LocalHttpTransport:
    def __init__(self, target):
        self.target = target
    def connect(self, target, timeout):
        # A redirect must never send credentials to a different origin.
        assert target == self.target
        return http.client.HTTPConnection(target.hostname, target.port, timeout=timeout)


@contextmanager
def synthetic_server(mode="pass"):
    requests = []
    closed_ignored_range = threading.Event()

    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"
        def log_message(self, *args):
            pass
        def send(self, status, body=b"", extra=None):
            self.send_response(status)
            for name, value in (extra or {}).items():
                self.send_header(name, value)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            if body:
                try:
                    self.wfile.write(body)
                except (BrokenPipeError, ConnectionResetError):
                    pass
        def do_GET(self):
            self.close_connection = True
            authenticated = self.headers.get("Authorization") == AUTH.header()
            requests.append((self.path, self.headers.get("Authorization")))
            if not authenticated:
                if mode == "public_alias" and "/img_data/" in self.path:
                    self.send(200, b"PRIVATE_RESPONSE_CANARY")
                else:
                    self.send(401, b"PRIVATE_ERROR_CANARY", {"WWW-Authenticate": 'Basic realm="PRIVATE_REALM_CANARY"'})
                return
            if self.path.endswith("websocket/"):
                key = self.headers["Sec-WebSocket-Key"]
                accept = base64.b64encode(hashlib.sha1(
                    (key + "258EAFA5-E914-47DA-95CA-C5AB0DC85B11").encode(),
                    usedforsecurity=False).digest()).decode()
                self.send(101, extra={"Upgrade": "websocket", "Connection": "Upgrade",
                          "Sec-WebSocket-Accept": "wrong" if mode == "bad_websocket" else accept})
                return
            if self.headers.get("Range"):
                if mode == "ignored_range":
                    # A client that attempts to consume this 200 body stalls.
                    self.send_response(200)
                    self.send_header("Content-Length", "999999999999")
                    self.end_headers()
                    self.wfile.flush()
                    self.connection.settimeout(2)
                    try:
                        if self.rfile.read(1) == b"":
                            closed_ignored_range.set()
                    except (OSError, ConnectionResetError):
                        pass
                    return
                body = b"x" * 16 if mode == "wrong_image" else PREFIX
                if mode == "oversized_range":
                    body += b"PRIVATE_EXTRA_BODY_CANARY"
                content_range = "bytes 1-16/100" if mode == "wrong_range" else "bytes 0-15/100"
                self.send(206, body, {"Content-Range": content_range})
                return
            if mode == "redirect" and self.path == "/app/":
                self.send(302, b"PRIVATE_REDIRECT_BODY", {"Location": "https://outside.invalid/PRIVATE_LOCATION_CANARY"})
                return
            if self.path.endswith(".xml"):
                mime = "text/html" if mode == "html_sidecar" else "application/xml"
                self.send(200, b"<xml>PRIVATE_METADATA_CANARY</xml>", {"Content-Type": mime})
            else:
                self.send(200, b"<!doctype html><html>PRIVATE_HTML_CANARY</html>", {"Content-Type": "text/html; charset=utf-8"})

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    server.daemon_threads = True
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    target = verifier.validate_url(f"https://127.0.0.1:{server.server_port}/app/")
    try:
        yield target, requests, closed_ignored_range
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def run_synthetic(target, **kwargs):
    return verifier.verify_access(target, AUTH, IMAGE, PREFIX_SHA,
                                  timeout=0.5, transport=LocalHttpTransport(target), **kwargs)


class ReviewerAccessTests(unittest.TestCase):
    def test_whole_site_auth_range_sidecar_and_websocket_pass(self):
        with synthetic_server() as (target, requests, _):
            report = run_synthetic(target, sidecar_paths=[SIDECAR])
        self.assertEqual(report["status"], "passed")
        self.assertEqual(report["counts"], {"checks": 13, "passed": 13, "failed": 0})
        self.assertEqual(len(requests), 13)
        self.assertTrue(all(auth is None for _, auth in requests[:8]))
        self.assertTrue(all(auth == AUTH.header() for _, auth in requests[8:]))
        self.assertNotIn("PRIVATE_", json.dumps(report))
        self.assertNotIn(AUTH.header(), json.dumps(report))

    def test_one_public_image_alias_fails_without_reading_protected_body(self):
        with synthetic_server("public_alias") as (target, _, _):
            report = run_synthetic(target)
        result = next(check for check in report["checks"] if check["check"] == "anonymous_img_data")
        self.assertEqual(result, {"check": "anonymous_img_data", "status": "failed", "http_status": 200, "body_bytes_read": 0})
        self.assertEqual(report["status"], "failed")
        self.assertNotIn("PRIVATE_", json.dumps(report))

    def test_ignored_range_closes_without_attempting_to_read_huge_tiff(self):
        with synthetic_server("ignored_range") as (target, _, closed):
            report = run_synthetic(target)
            self.assertTrue(closed.wait(1))
        result = next(check for check in report["checks"] if check["check"] == "authenticated_image_range")
        self.assertEqual(result["http_status"], 200)
        self.assertEqual(result["body_bytes_read"], 0)
        self.assertEqual(result["status"], "failed")

    def test_wrong_bytes_bad_range_and_excess_body_fail_bounded_checks(self):
        for mode, bytes_read in (("wrong_image", 16), ("wrong_range", 0), ("oversized_range", 17)):
            with self.subTest(mode=mode), synthetic_server(mode) as (target, _, _):
                report = run_synthetic(target)
            result = next(check for check in report["checks"] if check["check"] == "authenticated_image_range")
            self.assertEqual(result["status"], "failed")
            self.assertEqual(result["body_bytes_read"], bytes_read)

    def test_websocket_status_alone_does_not_accept_a_wrong_handshake(self):
        with synthetic_server("bad_websocket") as (target, _, _):
            report = run_synthetic(target)
        result = report["checks"][-1]
        self.assertEqual(result["http_status"], 101)
        self.assertEqual(result["status"], "failed")

    def test_proxy_html_fallback_does_not_pass_as_a_sidecar(self):
        with synthetic_server("html_sidecar") as (target, _, _):
            report = run_synthetic(target, sidecar_paths=[SIDECAR])
        result = next(check for check in report["checks"] if check["check"] == "authenticated_sidecar_1")
        self.assertEqual(result["http_status"], 200)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["body_bytes_read"], 0)

    def test_redirect_is_not_followed_or_printed(self):
        with synthetic_server("redirect") as (target, requests, _):
            report = run_synthetic(target)
        self.assertEqual(len(requests), 11)
        self.assertEqual(next(check for check in report["checks"] if check["check"] == "authenticated_app")["http_status"], 302)
        self.assertNotIn("PRIVATE_", json.dumps(report))
        self.assertNotIn("outside.invalid", json.dumps(report))

    def test_production_urls_and_relative_paths_reject_unsafe_inputs(self):
        for url in ("http://host/app/", "https://user:secret@host/app/", "https://host/app/?secret=x",
                    "https://host/app/#x", "https://host/app/?", "https://host/app", "https://host/a/../app/",
                    "https://host/%2e%2e/app/", "https://user%3Asecret%40host/app/", "https://host/app/\n"):
            with self.subTest(url=url), self.assertRaises(verifier.VerificationInputError):
                verifier.validate_url(url)
        for path in ("../image.ome.tiff", "/images/a.ome.tiff", "https://outside.invalid/a.ome.tiff",
                     "images/%2e%2e/a.ome.tiff", "images/a%2fb.ome.tiff", "images/a%252fb.ome.tiff",
                     "images/a.ome.tiff?x=y", "images/a.ome.tiff#x", "images/a.xml"):
            with self.subTest(path=path), self.assertRaises(verifier.VerificationInputError):
                verifier.validate_relative_path(path, image=True)

    def test_credentials_require_owned_0600_regular_json_and_hide_values(self):
        with tempfile.TemporaryDirectory() as directory:
            auth = Path(directory) / "auth.json"
            auth.write_text(json.dumps({"username": AUTH.username, "password": AUTH.password}))
            auth.chmod(0o600)
            self.assertEqual(verifier.read_credentials(auth).header(), AUTH.header())
            self.assertNotIn("PRIVATE_", repr(verifier.read_credentials(auth)))
            auth.chmod(0o640)
            with self.assertRaises(verifier.VerificationInputError):
                verifier.read_credentials(auth)
            auth.chmod(0o600)
            link = Path(directory) / "link.json"
            link.symlink_to(auth)
            with self.assertRaises(verifier.VerificationInputError):
                verifier.read_credentials(link)
            auth.write_text('{"password":"PRIVATE_ERROR_CANARY"}')
            with self.assertRaises(verifier.VerificationInputError):
                verifier.read_credentials(auth)
            fifo = Path(directory) / "fifo"
            os.mkfifo(fifo, 0o600)
            with self.assertRaises(verifier.VerificationInputError):
                verifier.read_credentials(fifo)

    def test_argument_and_transport_errors_never_serialize_raw_values(self):
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            self.assertEqual(verifier.main(["--PRIVATE_ARGUMENT_CANARY"]), 2)
        self.assertNotIn("PRIVATE_", stderr.getvalue())
        class BrokenTransport:
            def connect(self, target, timeout):
                raise OSError("PRIVATE_EXCEPTION_CANARY")
        report = verifier.verify_access(verifier.validate_url("https://example.invalid/app/"), AUTH,
                                        IMAGE, PREFIX_SHA, transport=BrokenTransport())
        self.assertEqual(report["counts"]["failed"], 11)
        self.assertNotIn("PRIVATE_", json.dumps(report))


if __name__ == "__main__":
    unittest.main(verbosity=2)
