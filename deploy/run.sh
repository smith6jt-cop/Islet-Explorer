#!/usr/bin/env bash
# Launch the Islet Explorer container on a UF Research Computing PubApps VM.
#
# RC's RSE team typically wires the *persistent* service (systemd/quadlet unit) for you —
# this script documents the exact bind mounts + secret injection that unit must perform,
# and is directly runnable for manual testing on the VM.
#
# Prereqs on the VM:
#   - deploy/islet-explorer.sif copied to /podman/$SERVICEUSER (or PODMAN alt: an image tag)
#   - Code + data + images staged under /pubapps/$PROJECT/islet-explorer (see deploy/README.md)
#   - A secrets env-file (chmod 600, NOT in git) holding the rotated UF Navigator key:
#         KEY=sk-...rotated...
#     Do NOT put KEY in the .Renviron on the shared volume; inject it here so it stays out
#     of the world-readable /pubapps tree.
set -euo pipefail

PROJECT="${PROJECT:?set PROJECT to your /pubapps project name}"
SERVICEUSER="${SERVICEUSER:-$USER}"
PUBAPPS_ROOT="/pubapps/${PROJECT}/islet-explorer"
SIF="${SIF:-/podman/${SERVICEUSER}/islet-explorer.sif}"
SECRET_ENV_FILE="${SECRET_ENV_FILE:-/pubapps/${PROJECT}/secrets/islet.env}"

# Persistent log directory. --writable-tmpfs makes the container filesystem ephemeral, so
# without this bind every shiny-server worker log is discarded when the instance stops —
# which is exactly what you need when a worker dies at startup (the classic failure here is
# a missing R package or a bad RETICULATE_PYTHON). README verification step 6 reads these.
LOG_DIR="${LOG_DIR:-/pubapps/${PROJECT}/logs}"
mkdir -p "${LOG_DIR}"

# --- Apptainer (default) ----------------------------------------------------------------
# Apptainer uses host networking, so shiny-server on :3838 is reachable on the VM directly —
# no port publishing needed. RC's reverse proxy then targets 127.0.0.1:3838.
#
# `instance start` runs the image's %startscript (NOT %runscript) — see deploy/islet-explorer.def.
apptainer instance start \
  --bind "${PUBAPPS_ROOT}:/srv/islet-explorer" \
  --bind "${LOG_DIR}:/var/log/shiny-server" \
  --env-file "${SECRET_ENV_FILE}" \
  --writable-tmpfs \
  "${SIF}" islet-explorer

echo "Started Apptainer instance 'islet-explorer'."
echo "  App:  http://127.0.0.1:3838/islet-explorer/   (RC proxy -> https://islet-explorer.rc.ufl.edu/...)"
echo "  Stop: apptainer instance stop islet-explorer"
echo "  Logs: tail -f ${LOG_DIR}/*.log   (bind-mounted from /var/log/shiny-server in the instance)"
echo "  List: apptainer instance list"

# --- Podman alternative (rootless) ------------------------------------------------------
# If RC standardizes on Podman + quadlet instead of Apptainer, the equivalent is:
#
#   podman run -d --name islet-explorer \
#     --publish 127.0.0.1:3838:3838 \
#     --volume "${PUBAPPS_ROOT}:/srv/islet-explorer:ro,z" \
#     --volume "${LOG_DIR}:/var/log/shiny-server:z" \
#     --env-file "${SECRET_ENV_FILE}" \
#     --userns=keep-id \
#
# The app volume is :ro deliberately — the documented update flow (`git pull` + `touch
# restart.txt`) runs on the HOST side of the bind, so read-only inside the container is
# sufficient and keeps a public endpoint from writing to the shared volume. The log volume
# must be read-write, hence the separate mount.
#     localhost/islet-explorer:latest
#
# Note --userns=keep-id (or a --uidmap) so the container's shiny worker can read the
# $serviceuser-owned files on /pubapps. Prefer a quadlet .container unit for auto-restart.
