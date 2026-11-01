# hpg_env.sh -- settings for the HiPerGator copy of the typing pipeline (owner, 2026-10-31). Source it; do not run it.
#
# The copy lives under $SITE and is presented at the workstation's paths inside an Apptainer container
# (run_in_site.sh), so records, QuPath image URIs and absolute symlinks need no rewriting. The file identities the
# guards compare come from the seal (site_seal.py, site_identity.py). See docs/hpg_processing_20261031.md.
# shellcheck shell=bash

export SITE=${SITE:-/blue/maigan/smith6jt/islet-explorer-hpg}
export SITE_ROOT=$SITE/root                    # bound as $WS_ROOT, one mount (activation renames across $RUN and output/)
export SITE_QPI=$SITE/qpi                      # bound read-only as $WS_QPI (source TIFFs, pixel classifiers)
export WS_ROOT=/data/panc_CODEX/Islet-Explorer
export WS_QPI=/data/QuPath_Projects/QuPath_project_IsletExplorer
export SITE_TRANSFER=$WS_ROOT/output/hpg_transfer_20261031       # capture and seal, as seen inside the container
export CODEX_SITE_SEAL=$SITE_TRANSFER/hpg_site_seal.json
# The seal of 2026-11-01 (job 48766607: 21,421 files and 57 symlinks verified against capture/workstation_identities.json
# 434295bf...). run_in_site.sh refuses to run with the hook if this is empty; a reseal pins its own sha256 here.
export CODEX_SITE_SEAL_SHA256=${CODEX_SITE_SEAL_SHA256:-54cad8c3c38e80eaf2405c3df34e4fac3168901b9fe892959691e56955f52ef0}
# QuPath 0.7's OpenCV natives need GLIBCXX_3.4.30; RHEL 9 stops at 3.4.29. One symlink, as in IO60's hpg_env.sh:
# $SITE/lib_stdcxx/libstdc++.so.6 -> /blue/maigan/smith6jt/miniforge3/lib/libstdc++.so.6.0.34
export SITE_LIBSTDCXX=$SITE/lib_stdcxx
# /home/smith6jt/miniconda3/envs/scvi-env on HPG is a symlink to the HPG scvi-env (same package versions), so the
# interpreter path hard-coded in run_remaining_qc_pipeline.py works inside and outside the container.
export SITE_PYTHON=/home/smith6jt/miniconda3/envs/scvi-env/bin/python
# Apptainer refuses the host root as a sandbox ("/ as sandbox is not authorized", job 48752194), so the container is a
# Rocky Linux 9 image (apptainer pull docker://rockylinux/rockylinux:9 on 2026-11-01: sif sha256 5d8efaca..., base digest
# sha256:8101994123cf...) with the host's /usr bound read-only (run_in_site.sh): git and the system libraries are the
# host's, as for IO60's QuPath runs on these nodes.
export SITE_CONTAINER=${SITE_CONTAINER:-$SITE/site/rocky9.sif}
