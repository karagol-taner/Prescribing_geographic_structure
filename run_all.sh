#!/usr/bin/env bash
# Runs steps 1 to 24 in order, from the aggregated data in data/ to every number, table and figure.
# About 2 to 2.5 hours on one core; most of the time is spent in need_similarity.py (about 30 minutes) and
# placebo_covariates.py (75 to 100 minutes). All permutation tests use fixed seeds, so the outputs match
# the files in data/.
#
# Quicker run (about 50 minutes): PLACEBO_ARGS="10 5 5" bash run_all.sh
# This uses 10, 5 and 5 placebo sets instead of 200, 100 and 100 and so overwrites the placebo results
# (data/placebo_summary.json, data/placebo_runs.csv), the placebo numbers in data/need_summary_numbers.txt and
# Figure S8; in a git clone, restore them with
#   git checkout -- data/placebo_summary.json data/placebo_runs.csv data/need_summary_numbers.txt figures/figS8_placebo.png
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}

step() {
    local start=$(date +%s)
    echo "== $*"
    "$@"
    echo "   finished in $(( $(date +%s) - start )) s"
}

step $PY code/analyze_af.py
step $PY code/temporal_af.py
step $PY code/fdr_refine.py
step $PY code/robustness_loo.py
step $PY code/mantel_84.py
step $PY code/spatial_weights.py
step $PY code/spatial_analysis.py
step $PY code/spatial_sensitivity.py
step $PY code/trends.py
step $PY code/make_geography_figures.py
step $PY code/make_figures.py
echo "== code/summary_numbers.py > data/summary_numbers.txt"
$PY code/summary_numbers.py > data/summary_numbers.txt
step $PY code/phylogeography_tables.py
step $PY code/need_adjustment.py
step $PY code/need_similarity.py
step $PY code/placebo_covariates.py ${PLACEBO_ARGS:-}
step $PY code/signflip_check.py
step $PY code/class_need.py
step $PY code/class_robustness.py
step $PY code/section_choice.py
step $PY code/lisa.py
echo "== code/need_summary_numbers.py > data/need_summary_numbers.txt"
$PY code/need_summary_numbers.py > data/need_summary_numbers.txt
step $PY code/supplementary_tables.py
step $PY code/make_need_figures.py
echo "All steps finished."
