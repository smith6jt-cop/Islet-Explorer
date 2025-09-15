#!/usr/bin/env bash
set -euo pipefail

export STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

streamlit run app_v2.py \
  --server.address 0.0.0.0 \
  --server.port "${PORT:-8501}"

