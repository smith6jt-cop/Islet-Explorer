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
# Pinned after the seal job (step E1); run_in_site.sh refuses to run with the hook while it is empty.
export CODEX_SITE_SEAL_SHA256=${CODEX_SITE_SEAL_SHA256:-}
# QuPath 0.7's OpenCV natives need GLIBCXX_3.4.30; RHEL 9 stops at 3.4.29. One symlink, as in IO60's hpg_env.sh:
# $SITE/lib_stdcxx/libstdc++.so.6 -> /blue/maigan/smith6jt/miniforge3/lib/libstdc++.so.6.0.34
export SITE_LIBSTDCXX=$SITE/lib_stdcxx
# /home/smith6jt/miniconda3/envs/scvi-env on HPG is a symlink to the HPG scvi-env (same package versions), so the
# interpreter path hard-coded in run_remaining_qc_pipeline.py works inside and outside the container.
export SITE_PYTHON=/home/smith6jt/miniconda3/envs/scvi-env/bin/python
export SITE_CONTAINER=${SITE_CONTAINER:-/}     # the host root as a directory container; a Rocky Linux 9 image if refused
