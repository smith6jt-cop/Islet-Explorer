#!/usr/bin/python3
"""Pass approved container configuration through Shiny Server's R environment reset."""
import json
import os
from pathlib import Path
import shutil
import stat
import sys
import tempfile

# Private deployment runtime directory; mkstemp creates mode 0600 before replacement.
# Worker reads also verify the file owner and mode.
PROFILE = Path("/tmp/islet-explorer/worker-environment.json")  # nosec B108
ALLOWED = frozenset({
    "ISLET_RELEASE", "ISLET_RUNTIME_MODE", "ISLET_DATA_ROOT", "ISLET_DATA_SCHEMA",
    "LOCAL_IMAGE_ROOT", "RETICULATE_PYTHON", "PYTHONNOUSERSITE", "ISLET_ENABLE_AI",
    "SHOW_STATISTICS_TAB", "SHOW_INFLAMMATION_TAB", "DEBUG_CREDENTIALS", "VIEWER_DEBUG",
    "TMPDIR", "R_USER", "XDG_CACHE_HOME", "MPLCONFIGDIR", "NUMBA_CACHE_DIR",
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    "ISLET_DUCKDB_THREADS", "ISLET_DUCKDB_MEMORY_LIMIT", "ISLET_DUCKDB_TEMP_DIR",
    "ISLET_VERIFY_RELEASE_HASHES",
})


def main():
    if sys.argv[1:] == ["--write-profile"]:
        profile = {key: os.environ[key] for key in sorted(ALLOWED) if key in os.environ}
        PROFILE.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
        descriptor, temporary = tempfile.mkstemp(prefix="worker-environment-", dir=PROFILE.parent)
        try:
            with os.fdopen(descriptor, "w") as stream:
                json.dump(profile, stream)
                stream.write("\n")
            os.replace(temporary, PROFILE)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        return

    metadata = PROFILE.stat()
    if metadata.st_uid != os.getuid() or stat.S_IMODE(metadata.st_mode) != 0o600:
        raise ValueError("Invalid profile permissions")
    profile = json.loads(PROFILE.read_text())
    if not isinstance(profile, dict) or not set(profile).issubset(ALLOWED) or any(
        not isinstance(value, str) or "\0" in value for value in profile.values()
    ):
        raise ValueError("Invalid profile content")
    environment = {key: os.environ[key] for key in ("HOME", "LANG", "PATH") if key in os.environ}
    environment.update(profile)
    # Disable user startup files before R initializes. Credentials and arbitrary
    # environment variables are never copied from the server or profile.
    environment.update(R_ENVIRON_USER="/dev/null", R_PROFILE_USER="/dev/null")
    r_executable = shutil.which("R", path=environment.get("PATH"))
    if not r_executable:
        raise ValueError("R executable unavailable")
    # Fixed R program, argument array, and allowlisted environment; no shell.
    os.execve(r_executable, [r_executable, *sys.argv[1:]], environment)  # nosec B606


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, TypeError):
        print("ERROR: R worker environment profile is unavailable or invalid.", file=sys.stderr)
        sys.exit(78)
