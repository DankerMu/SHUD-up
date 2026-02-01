#!/usr/bin/env bash
set -euo pipefail

sd="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${sd}/_common.sh"

need_cmd "${sd}/../../shud"

root="$(repo_root)"

# Experiments: keep END/DT_QE_ET from input/ccw by default; only vary MAX_SOLVER_STEP.
# Use distinct output dirs to avoid overwriting your main ccw.tsr results.
for max_step in 5 2; do
  export SHUD_VALIDATION_MAX_SOLVER_STEP_MIN="${max_step}"
  run_ccw_case 1 "output/ccw.tsr.ms${max_step}" "validation/tsr/tmp/ccw.tsr.ms${max_step}" "tsr-ms${max_step}"
done

echo "OK: generated TSR max-step experiment outputs under:"
echo "- ${root}/output/ccw.tsr.ms5"
echo "- ${root}/output/ccw.tsr.ms2"

