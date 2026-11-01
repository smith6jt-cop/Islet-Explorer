#!/usr/bin/env bash
# run_in_site.sh [--no-seal] CMD [ARGS...] -- run CMD on HiPerGator inside the site container, at the workstation's
# paths (hpg_env.sh), from /data/panc_CODEX/Islet-Explorer. By default the site identity hook is active
# (CODEX_SITE_SEAL with the pinned sha256); --no-seal runs without it, for site_seal.py, the smoke test and negative
# controls. The environment is cleaned; only the variables below and the SLURM ids are passed in.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=hpg_env.sh
source "$HERE/hpg_env.sh"

seal=1
if [[ ${1:-} == --no-seal ]]; then seal=0; shift; fi
(($#)) || { echo "usage: run_in_site.sh [--no-seal] CMD [ARGS...]" >&2; exit 2; }

env=(--env "TMPDIR=${TMPDIR:-$SITE/tmp}" --env "LD_LIBRARY_PATH=$SITE_LIBSTDCXX"
     --env "PYTHONPATH=$WS_ROOT/scripts/hpg" --env "PYTHONDONTWRITEBYTECODE=1"
     --env "PATH=$(dirname "$SITE_PYTHON"):/usr/local/bin:/usr/bin:/bin")
for v in SLURM_JOB_ID SLURM_JOB_NAME SLURM_CPUS_PER_TASK SLURM_MEM_PER_NODE SLURM_JOB_NODELIST JAVA_TOOL_OPTIONS; do
    if [[ -n ${!v:-} ]]; then env+=(--env "$v=${!v}"); fi
done
if ((seal)); then
    [[ -n $CODEX_SITE_SEAL_SHA256 ]] || { echo "run_in_site.sh: CODEX_SITE_SEAL_SHA256 is not pinned; seal first" >&2; exit 78; }
    env+=(--env "CODEX_SITE_SEAL=$CODEX_SITE_SEAL" --env "CODEX_SITE_SEAL_SHA256=$CODEX_SITE_SEAL_SHA256")
fi

exec apptainer exec --cleanenv "${env[@]}" \
    --bind "$SITE_ROOT:$WS_ROOT" --bind "$SITE_QPI:$WS_QPI:ro" --bind /usr:/usr:ro \
    --pwd "$WS_ROOT" "$SITE_CONTAINER" "$@"
