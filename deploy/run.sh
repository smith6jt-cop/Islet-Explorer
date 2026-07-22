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

# --- Apptainer (default) ----------------------------------------------------------------
# Apptainer uses host networking, so shiny-server on :3838 is reachable on the VM directly —
# no port publishing needed. RC's reverse proxy then targets 127.0.0.1:3838.
apptainer instance start \
  --bind "${PUBAPPS_ROOT}:/srv/islet-explorer" \
  --env-file "${SECRET_ENV_FILE}" \
  --writable-tmpfs \
  "${SIF}" islet-explorer

echo "Started Apptainer instance 'islet-explorer'."
echo "  App:  http://127.0.0.1:3838/islet-explorer/   (RC proxy -> https://islet-explorer.rc.ufl.edu/...)"
echo "  Stop: apptainer instance stop islet-explorer"
echo "  Logs: apptainer instance list; tail -f ${PUBAPPS_ROOT}/app/shiny_app/var/log/... (or /var/log/shiny-server inside the instance)"

# --- Podman alternative (rootless) ------------------------------------------------------
# If RC standardizes on Podman + quadlet instead of Apptainer, the equivalent is:
#
#   podman run -d --name islet-explorer \
#     --publish 127.0.0.1:3838:3838 \
#     --volume "${PUBAPPS_ROOT}:/srv/islet-explorer:ro,z" \
#     --env-file "${SECRET_ENV_FILE}" \
#     --userns=keep-id \
#     localhost/islet-explorer:latest
#
# Note --userns=keep-id (or a --uidmap) so the container's shiny worker can read the
# $serviceuser-owned files on /pubapps. Prefer a quadlet .container unit for auto-restart.
