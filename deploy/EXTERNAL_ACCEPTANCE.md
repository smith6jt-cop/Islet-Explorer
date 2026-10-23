# External reviewer access acceptance

Run this after UFRC configures the canonical HTTPS endpoint and whole-site Basic
Authentication. UFRC's vhosts went live on 2026-10-23; from the campus network,
with UFRC's temporary test account, all 15 checks passed on both hostnames that
day (`deploy/README.md` § 4). No external pass is claimed yet: the run from a
clean non-UF network with the real reviewer credential file is still required.
Python 3 on Linux is sufficient.

Create an owner-only **0600** regular JSON file outside the checkout containing
`{"username":"<reviewer username>","password":"<reviewer password>"}`. Edit it
locally; never put either value in command arguments, environment variables,
version control or evidence. Symlinks and files with other permissions are rejected.

Compute the SHA-256 of exactly the first 16 bytes of the approved local TIFF:

```sh
head -c 16 /restricted/approved-image.ome.tiff | sha256sum
```

Run against the canonical app prefix, including its trailing slash:

```sh
python3 deploy/verify-reviewer-access.py \
  --url https://reviewer.example.ufl.edu/islet-explorer/ \
  --image-path local_images/approved-image.ome.tiff \
  --image-prefix-sha256 '<64-character digest from the preceding command>' \
  --sidecar-path local_images/approved-image.ome.xml \
  --auth-file /restricted/reviewer-auth.json \
  --output /restricted/evidence/external-acceptance.json
```

Repeat `--sidecar-path` for known XML/JSON sidecars. The default WebSocket route is
Shiny Server's SockJS endpoint (`__sockjs__/n=<token>/<server>/<session>/websocket`,
with placeholder values; the bare `websocket/` route of a plain Shiny worker hangs
behind Shiny Server until the proxy times out); use `--websocket-path` when the
deployed route differs. The anonymous download probe defaults to
`session/auth-boundary-probe/download/auth-boundary-probe`; it checks that proxy
authentication runs before Shiny session lookup. `--download-path` can instead
name a real app-relative download route observed in an active reviewer session.
This probe does not claim to validate an authenticated scientific download.

Acceptance requires anonymous HTTP **401 with a Basic challenge** for the app,
Avivator entry point, all three `local_images`/`images`/`img_data` aliases, supplied
sidecars, the download route and WebSocket upgrade. Authenticated checks require
app/viewer HTML, known sidecar HTTP 200 with XML/JSON content type, exactly the approved TIFF prefix with
HTTP 206 and `Content-Range: bytes 0-15/TOTAL`, and HTTP 101 with the correct
WebSocket accept value. Authentication is sent only to the original HTTPS origin;
certificate validation is mandatory and redirects are never followed. The client
closes immediately if a Range request returns HTTP 200, without reading its body.
The authenticated Range check uses the supplied image path; all three aliases
are checked anonymously for whole-site authentication coverage.

Exit codes are 0 for all checks passing, 1 for a failed acceptance check, and 2
for invalid inputs or an evidence-write failure. The new JSON file contains only
fixed check labels, pass/fail status, HTTP statuses and byte/check counts. No URL,
filename, credential, response body, scientific row or raw exception is emitted.
Keep the operator's mapping from this evidence to the exact deployed image/data
release in the restricted deployment record. Local synthetic HTTP tests exercise
the protocol logic through an injected transport; the production CLI exposes no
HTTP or certificate-validation bypass.
