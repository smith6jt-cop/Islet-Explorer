#!/bin/bash
#SBATCH --job-name=codex-scvi-21
#SBATCH --account=maigan
#SBATCH --partition=hpg-rtx6000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --gres=gpu:rtx_pro_6000:1
#SBATCH --time=12:00:00
#SBATCH --output=/blue/maigan/smith6jt/joint-scvi-20260929/logs/slurm-%j.out
#SBATCH --error=/blue/maigan/smith6jt/joint-scvi-20260929/logs/slurm-%j.err

set -euo pipefail

JOB_ROOT=/blue/maigan/smith6jt/joint-scvi-20260929
PYTHON=/blue/maigan/smith6jt/miniforge3/envs/scvi-env/bin/python
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPLCONFIGDIR="$JOB_ROOT/mpl-cache"

printf 'Job %s on %s at %s\n' "$SLURM_JOB_ID" "$(hostname)" "$(date -Is)"
nvidia-smi
srun --ntasks=1 "$PYTHON" "$JOB_ROOT/retrain_joint_scvi.py" train \
  --output "$JOB_ROOT" --max-epochs 400 --batch-size 2048 --workers 8
printf 'Finished at %s\n' "$(date -Is)"
