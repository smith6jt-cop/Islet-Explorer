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
    with tempfile.TemporaryDirectory(prefix="islet-server-smoke-") as temporary:
        root = Path(temporary)
        (root / "app").mkdir()
        (root / "logs").mkdir()
        (root / "bookmarks").mkdir()
        (root / "images").mkdir()
        (root / "images/range.bin").write_bytes(bytes(range(256)))
        (root / "app/app.R").write_text(
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
            process = subprocess.Popen([str(server / "bin/shiny-server"), str(config)], stdout=log, stderr=log)
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
                    "status": "passed", "checks": ["real R application page", "Shiny JavaScript through proxy", "Shiny Server CSS under application prefix", "all three production image routes return HTTP 206 and exactly 16 requested bytes"],
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
