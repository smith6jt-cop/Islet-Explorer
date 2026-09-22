#!/usr/bin/env python3
"""Exercise the staged server after npm and development artifacts are omitted."""
import json
import os
from pathlib import Path
import pwd
import socket
import subprocess
import sys
import tempfile
import time
from urllib.error import URLError
from urllib.request import Request, urlopen


def main():
    server = Path(sys.argv[1]).resolve()
    report = Path(sys.argv[2]).resolve()
    config_template = Path(sys.argv[3]).read_text()
    wrapper_source = Path(sys.argv[4]).read_text()
    with tempfile.TemporaryDirectory(prefix="islet-server-smoke-") as temporary:
        root = Path(temporary)
        (root / "app").mkdir()
        (root / "logs").mkdir()
        (root / "bookmarks").mkdir()
        (root / "images").mkdir()
        (root / "r-temp").mkdir()
        (root / "images/range.bin").write_bytes(bytes(range(256)))
        profile = root / "worker-environment.json"
        wrapper = root / "r-worker"
        assert 'Path("/tmp/islet-explorer/worker-environment.json")' in wrapper_source
        wrapper.write_text(wrapper_source.replace(
            'Path("/tmp/islet-explorer/worker-environment.json")', f"Path({str(profile)!r})"
        ))
        wrapper.chmod(0o555)
        environment = os.environ.copy()
        expected = {
            "ISLET_RUNTIME_MODE": "production", "ISLET_DATA_ROOT": str(root / "data with spaces"),
            "ISLET_DATA_SCHEMA": "1", "LOCAL_IMAGE_ROOT": str(root / "images"),
            "RETICULATE_PYTHON": "/configured/python", "ISLET_ENABLE_AI": "FALSE",
            "ISLET_RELEASE": 'review "quoted" \\ ${HOME}; # literal\nvalue',
            "TMPDIR": str(root / "r-temp"), "OMP_NUM_THREADS": "1",
            "R_ENVIRON_USER": "/dev/null", "R_PROFILE_USER": "/dev/null",
        }
        environment.update(expected)
        environment.update(KEY="synthetic-not-a-credential", OPENAI_API_KEY="synthetic-not-a-credential",
                           ISLET_UNAPPROVED="must-not-reach-R")
        subprocess.run([str(wrapper), "--write-profile"], env=environment, check=True)
        assert profile.stat().st_mode & 0o777 == 0o600
        assert not {"KEY", "OPENAI_API_KEY", "ISLET_UNAPPROVED"}.intersection(json.loads(profile.read_text()))
        (root / "app/expected-environment.json").write_text(json.dumps(expected))
        environment["R"] = str(wrapper)
        (root / "app/app.R").write_text(
            'expected <- unlist(jsonlite::fromJSON("expected-environment.json"))\n'
            'stopifnot(identical(Sys.getenv(names(expected)), expected), '
            '!nzchar(Sys.getenv("KEY")), !nzchar(Sys.getenv("OPENAI_API_KEY")), '
            '!nzchar(Sys.getenv("ISLET_UNAPPROVED")))\n'
            'shiny::shinyApp(shiny::fluidPage("islet-server-smoke"), '
            'function(input, output, session) {})\n'
        )
        with socket.socket() as listener:
            listener.bind(("127.0.0.1", 0))
            port = listener.getsockname()[1]
        config = root / "shiny-server.conf"
        replacements = {
            "run_as shiny;": f"run_as {pwd.getpwuid(os.getuid()).pw_name};",
            "listen 3838;": f"listen {port} 127.0.0.1;",
            "/srv/islet-explorer/app/shiny_app/www/local_images": f"{root}/images",
            "/srv/islet-explorer/app/shiny_app": f"{root}/app",
            "/var/log/shiny-server": f"{root}/logs",
        }
        for original, replacement in replacements.items():
            assert original in config_template, f"production config changed: {original}"
            config_template = config_template.replace(original, replacement)
        config.write_text(f"bookmark_state_dir {root}/bookmarks;\n" + config_template)
        with (root / "server.log").open("w+") as log:
            process = subprocess.Popen([str(server / "bin/shiny-server"), str(config)],
                                       stdout=log, stderr=log, env=environment)
            try:
                base_url = f"http://127.0.0.1:{port}/islet-explorer/"
                deadline = time.monotonic() + 90
                while True:
                    if process.poll() is not None:
                        raise RuntimeError("staged Shiny Server exited during startup")
                    try:
                        with urlopen(base_url, timeout=65) as response:
                            page = response.read()
                        break
                    except URLError:
                        if time.monotonic() >= deadline:
                            raise
                        time.sleep(0.2)
                assert b"islet-server-smoke" in page, page[:1000]
                for asset in ["shared/shiny.js", "__assets__/shiny-server.css"]:
                    with urlopen(base_url + asset, timeout=10) as response:
                        assert response.status == 200 and len(response.read()) > 100
                for alias in ["local_images", "images", "img_data"]:
                    request = Request(base_url + alias + "/range.bin", headers={"Range": "bytes=0-15"})
                    with urlopen(request, timeout=10) as response:
                        assert response.status == 206
                        assert response.headers["Content-Range"] == "bytes 0-15/256"
                        assert response.read() == bytes(range(16))
                report.write_text(json.dumps({
                    "status": "passed", "checks": ["real R worker receives exact allowlisted production configuration including literal special characters", "private mode-0600 profile omits credentials and unknown variables", "R user startup files disabled before initialization", "real R application page", "Shiny JavaScript through proxy", "Shiny Server CSS under application prefix", "all three production image routes return HTTP 206 and exactly 16 requested bytes"],
                    "runtime": str(server), "npm_cli_present": (server / "ext/node/lib/node_modules/npm").exists()
                }, indent=2) + "\n")
            except BaseException:
                log.flush()
                log.seek(0)
                print(log.read(), file=sys.stderr)
                raise
            finally:
                process.terminate()
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()


if __name__ == "__main__":
    main()
