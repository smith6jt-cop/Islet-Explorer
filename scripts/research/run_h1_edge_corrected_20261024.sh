#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
H1_PYTHON=${H1_PYTHON:-/home/smith6jt/miniconda3/envs/scvi-env/bin/python}
H1_OUT=${H1_OUT:-output/h1_edge_corrected_20261024}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$H1_OUT"
"$H1_PYTHON" -m pytest -q scripts/tests/test_h1_edge_corrections.py > "$H1_OUT/tests_run.log" 2>&1
"$H1_PYTHON" scripts/research/h1_validate_edge_20261024.py --output "$H1_OUT" --snapshot-only
"$H1_PYTHON" scripts/research/h1_edge_eligibility_20261024.py \
  --source-root data/development/uniform-20donor-20261015 \
  --app-root data/development/uniform-20donor-20261022-v10 \
  --output "$H1_OUT/eligibility" --workers 4 > "$H1_OUT/edge_run.log" 2>&1
"$H1_PYTHON" scripts/research/h1_features_20261022.py \
  --data-root data/development/uniform-20donor-20261022-v10 \
  --output-root "$H1_OUT" --eligibility-root "$H1_OUT/eligibility" \
  --workers 4 > "$H1_OUT/features_run.log" 2>&1
"$H1_PYTHON" scripts/research/h1_validate_edge_20261024.py \
  --output "$H1_OUT" > "$H1_OUT/validation_run.log" 2>&1
"$H1_PYTHON" scripts/research/h1_analysis_20261022.py \
  --features "$H1_OUT/h1_object_features.csv" --output "$H1_OUT/models" \
  --bootstrap 999 --workers 4 > "$H1_OUT/analysis_run.log" 2>&1
"$H1_PYTHON" scripts/research/h1_matched_20261022.py \
  --features "$H1_OUT/h1_object_features.csv" --output "$H1_OUT/independent" \
  > "$H1_OUT/matched_run.log" 2>&1
"$H1_PYTHON" scripts/research/h1_heterogeneity_20261022.py \
  --features "$H1_OUT/h1_object_features.csv" --output "$H1_OUT/heterogeneity" \
  > "$H1_OUT/heterogeneity_run.log" 2>&1
"$H1_PYTHON" scripts/research/h1_edge_report_20261024.py --output "$H1_OUT"
