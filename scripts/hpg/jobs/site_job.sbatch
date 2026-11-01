#!/bin/bash
# One step of the HiPerGator typing copy: runs its arguments through run_in_site.sh, timed. Resources are given at
# submission, for example (QOS maigan or clive; maigan-b allowed for this work, owner 2026-10-31):
#   sbatch --job-name=gates --cpus-per-task=8 --mem=64G --time=02:00:00 site_job.sbatch \
#       python scripts/hpg/gate_dump.py --out output/hpg_transfer_20261031/gate_dump_hpg.json
#   sbatch --job-name=seal --cpus-per-task=32 --mem=64G --time=03:00:00 site_job.sbatch --no-seal \
#       python scripts/hpg/site_seal.py verify-and-seal ...
# Submit from a login node: ssh hpg 'cd /blue/maigan/smith6jt/islet-explorer-hpg/root && /opt/slurm/bin/sbatch ...'.
#SBATCH --account=maigan
#SBATCH --qos=maigan
#SBATCH --partition=hpg-turin
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/blue/maigan/smith6jt/islet-explorer-hpg/logs/%x-%j.out
set -eo pipefail
source /etc/profile.d/modules.sh
SCRIPTS=/blue/maigan/smith6jt/islet-explorer-hpg/root/scripts/hpg
mkdir -p /blue/maigan/smith6jt/islet-explorer-hpg/tmp
export TMPDIR=/blue/maigan/smith6jt/islet-explorer-hpg/tmp
set -u
echo "host=$(hostname) job=$SLURM_JOB_ID cpus=${SLURM_CPUS_PER_TASK:-} start=$(date -Is) args=$*"
/usr/bin/time -v "$SCRIPTS/run_in_site.sh" "$@"
echo "end=$(date -Is)"
