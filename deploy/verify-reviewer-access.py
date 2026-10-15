#!/usr/bin/env python3
"""Check the external HTTPS/BasicAuth boundary without emitting protected values."""
from __future__ import annotations

import argparse
import base64
from dataclasses import dataclass, field
import hashlib
import http.client
import json
import os
from pathlib import Path
import re
import secrets
import ssl
import stat
import sys
from urllib.parse import quote, unquote, urlsplit


class VerificationInputError(Exception):
    """An intentionally value-free configuration error."""


class SafeArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise VerificationInputError("invalid_arguments")


@dataclass(frozen=True)
class Credentials:
    username: str = field(repr=False)
    password: str = field(repr=False)

    def header(self) -> str:
        encoded = base64.b64encode(f"{self.username}:{self.password}".encode("utf-8")).decode("ascii")
        return "Basic " + encoded


@dataclass(frozen=True)
class Target:
    hostname: str
    port: int
    origin: str
    prefix: str


def validate_url(value: str) -> Target:
    try:
        if re.search(r"[\s\\?#\x00-\x1f\x7f]", value):
            raise ValueError
        parsed = urlsplit(value)
        if (parsed.scheme != "https" or not parsed.hostname or parsed.username is not None
                or parsed.password is not None or parsed.query or parsed.fragment
                or not parsed.path.startswith("/") or not parsed.path.endswith("/")
                or "%" in parsed.netloc or not parsed.hostname.isascii()):
            raise ValueError
        port = parsed.port or 443
        if not 1 <= port <= 65535:
            raise ValueError
        segments = parsed.path.split("/")[1:-1]
        if any(not re.fullmatch(r"[A-Za-z0-9_~-]+(?:\.[A-Za-z0-9_~-]+)*", part) for part in segments):
            raise ValueError
        return Target(parsed.hostname, port, f"https://{parsed.netloc}", parsed.path)
    except (ValueError, TypeError):
        raise VerificationInputError("invalid_https_app_url") from None


def validate_relative_path(value: str, *, image=False, sidecar=False) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise VerificationInputError("invalid_relative_path")
    parsed = urlsplit(value)
    if parsed.scheme or parsed.netloc or parsed.query or parsed.fragment or value.startswith("/"):
        raise VerificationInputError("invalid_relative_path")
    decoded = unquote(value, errors="strict")
    # Inspect before URL normalization or a proxy's second decoding pass.
    if re.search(r"[\\%?#\x00-\x1f\x7f]", decoded) or re.search(r"%2f", value, re.I):
        raise VerificationInputError("invalid_relative_path")
    parts = decoded.split("/")
    if parts[-1] == "":
        parts.pop()
    if not parts or any(part in ("", ".", "..") for part in parts):
        raise VerificationInputError("invalid_relative_path")
    if image or sidecar:
        if len(parts) != 2 or parts[0] not in ("local_images", "images", "img_data"):
            raise VerificationInputError("invalid_static_path")
        extension = r"\.ome\.tiff?$" if image else r"\.(?:json|xml)$"
        if not re.search(extension, parts[1], re.I) or decoded.endswith("/"):
            raise VerificationInputError("invalid_static_path")
    return quote(decoded, safe="/-._~")


def read_credentials(path: Path) -> Credentials:
    descriptor = None
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        info = os.fstat(descriptor)
        if (not stat.S_ISREG(info.st_mode) or stat.S_IMODE(info.st_mode) != 0o600
                or info.st_uid != os.geteuid() or info.st_size > 8192):
            raise ValueError
        with os.fdopen(descriptor, "rb") as stream:
            descriptor = None
            content = json.loads(stream.read(8193))
        if not isinstance(content, dict) or set(content) != {"username", "password"}:
            raise ValueError
        username, password = content["username"], content["password"]
        if (not isinstance(username, str) or not isinstance(password, str)
                or not username or not password or ":" in username
                or any(char in username + password for char in "\r\n\x00")):
            raise ValueError
        return Credentials(username, password)
    except (OSError, ValueError, TypeError, UnicodeError):
        raise VerificationInputError("credential_file_invalid") from None
    finally:
        if descriptor is not None:
            os.close(descriptor)


class HttpsTransport:
    """Standard certificate/hostname validation; no redirects or proxy inheritance."""
    def connect(self, target: Target, timeout: float):
        return http.client.HTTPSConnection(target.hostname, target.port, timeout=timeout,
                                           context=ssl.create_default_context())


def websocket_headers(origin: str) -> tuple[dict, str]:
    key = base64.b64encode(secrets.token_bytes(16)).decode("ascii")
    # SHA-1 is required by the WebSocket handshake, not used as a password hash.
    accepted = base64.b64encode(hashlib.sha1(
        (key + "258EAFA5-E914-47DA-95CA-C5AB0DC85B11").encode("ascii"),
        usedforsecurity=False).digest()).decode("ascii")
    return {"Upgrade": "websocket", "Connection": "Upgrade", "Sec-WebSocket-Key": key,
            "Sec-WebSocket-Version": "13", "Origin": origin}, accepted


def verify_access(target: Target, credentials: Credentials, image_path: str,
                  prefix_sha256: str, *, sidecar_paths=(), websocket_path="websocket/",
                  download_path="session/auth-boundary-probe/download/auth-boundary-probe",
                  timeout=10.0, transport=None) -> dict:
    """Transport injection is for synthetic tests only; the CLI always uses TLS."""
    transport = transport or HttpsTransport()
    image_path = validate_relative_path(image_path, image=True)
    sidecar_paths = [validate_relative_path(path, sidecar=True) for path in sidecar_paths]
    websocket_path = validate_relative_path(websocket_path)
    download_path = validate_relative_path(download_path)
    if not re.fullmatch(r"[0-9a-fA-F]{64}", prefix_sha256):
        raise VerificationInputError("invalid_prefix_sha256")
    checks = []

    def request(label, relative, *, authenticated=False, kind="anonymous"):
        result = {"check": label, "status": "failed", "http_status": None, "body_bytes_read": 0}
        connection = None
        try:
            headers = {"Accept-Encoding": "identity", "User-Agent": "Islet-Explorer-Acceptance/1"}
            expected_accept = None
            if "websocket" in kind:
                extra, expected_accept = websocket_headers(target.origin)
                headers.update(extra)
            if authenticated:
                headers["Authorization"] = credentials.header()
            if kind == "range":
                headers["Range"] = "bytes=0-15"
            connection = transport.connect(target, timeout)
            connection.request("GET", target.prefix + relative, headers=headers)
            response = connection.getresponse()
            result["http_status"] = response.status
            # Never follow redirects, including same-origin redirects. Never
            # read an error/anonymous response body or a Range-ignored HTTP 200.
            passed = False
            if kind.startswith("anonymous"):
                challenge = response.getheader("WWW-Authenticate", "")
                passed = response.status == 401 and bool(re.search(r"\bBasic\b", challenge, re.I))
            elif kind == "html" and response.status == 200:
                content_type = response.getheader("Content-Type", "").split(";", 1)[0].strip().lower()
                if content_type == "text/html":
                    body = response.read(8192)
                    result["body_bytes_read"] = len(body)
                    passed = bool(re.search(br"<html\b|<!doctype\s+html\b", body, re.I))
            elif kind == "range" and response.status == 206:
                content_range = re.fullmatch(r"bytes 0-15/([0-9]+)", response.getheader("Content-Range", ""))
                if content_range and int(content_range[1]) >= 16:
                    body = response.read(17) # Bounded even if the server lies.
                    result["body_bytes_read"] = len(body)
                    passed = len(body) == 16 and hashlib.sha256(body).hexdigest() == prefix_sha256.lower()
            elif kind == "sidecar":
                content_type = response.getheader("Content-Type", "").split(";", 1)[0].strip().lower()
                allowed = ("application/json",) if relative.lower().endswith(".json") else ("application/xml", "text/xml")
                passed = response.status == 200 and content_type in allowed
            elif kind == "websocket":
                tokens = {part.strip().lower() for part in response.getheader("Connection", "").split(",")}
                passed = (response.status == 101 and "upgrade" in tokens
                          and response.getheader("Upgrade", "").lower() == "websocket"
                          and response.getheader("Sec-WebSocket-Accept", "") == expected_accept)
            if passed:
                result["status"] = "passed"
        except (OSError, ValueError, http.client.HTTPException):
            # Exception strings can contain URLs, credentials or response data.
            pass
        finally:
            if connection is not None:
                connection.close()
        checks.append(result)

    basename = image_path.split("/", 1)[1]
    anonymous = [("anonymous_app", ""), ("anonymous_avivator", "avivator/index.html")]
    anonymous += [("anonymous_" + alias, alias + "/" + basename) for alias in ("local_images", "images", "img_data")]
    anonymous.append(("anonymous_download", download_path))
    anonymous += [(f"anonymous_sidecar_{i}", path) for i, path in enumerate(sidecar_paths, 1)]
    for label, path in anonymous:
        request(label, path)
    request("anonymous_websocket", websocket_path, kind="anonymous_websocket")
    request("authenticated_app", "", authenticated=True, kind="html")
    request("authenticated_avivator", "avivator/index.html", authenticated=True, kind="html")
    request("authenticated_image_range", image_path, authenticated=True, kind="range")
    for i, path in enumerate(sidecar_paths, 1):
        request(f"authenticated_sidecar_{i}", path, authenticated=True, kind="sidecar")
    request("authenticated_websocket", websocket_path, authenticated=True, kind="websocket")
    passed = sum(item["status"] == "passed" for item in checks)
    return {"schema_version": 1, "status": "passed" if passed == len(checks) else "failed",
            "counts": {"checks": len(checks), "passed": passed, "failed": len(checks) - passed}, "checks": checks}


def main(argv=None) -> int:
    parser = SafeArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--image-path", required=True)
    parser.add_argument("--image-prefix-sha256", required=True)
    parser.add_argument("--auth-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sidecar-path", action="append", default=[])
    parser.add_argument("--websocket-path", default="websocket/")
    parser.add_argument("--download-path", default="session/auth-boundary-probe/download/auth-boundary-probe")
    parser.add_argument("--timeout", type=float, default=10.0)
    try:
        args = parser.parse_args(argv)
        target = validate_url(args.url)
        if not 0 < args.timeout <= 30 or args.output.exists():
            raise VerificationInputError("invalid_output_or_timeout")
        credentials = read_credentials(args.auth_file)
        report = verify_access(target, credentials, args.image_path, args.image_prefix_sha256,
                               sidecar_paths=args.sidecar_path, websocket_path=args.websocket_path,
                               download_path=args.download_path, timeout=args.timeout)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        descriptor = os.open(args.output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2)
            stream.write("\n")
        print("External reviewer acceptance: " + report["status"])
        return 0 if report["status"] == "passed" else 1
    except (VerificationInputError, OSError, ValueError, UnicodeError):
        print("External reviewer acceptance could not run: check HTTPS paths, protected credential file and output settings.", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
