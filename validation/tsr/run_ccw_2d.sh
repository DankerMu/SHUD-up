#!/usr/bin/env bash
set -euo pipefail

sd="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${sd}/_common.sh"

# 2-day short run used by docs/water_balance_verification.md (hourly DT_QE_ET).
export SHUD_VALIDATION_END_DAYS="${SHUD_VALIDATION_END_DAYS:-2}"
export SHUD_VALIDATION_DT_QE_ET_MIN="${SHUD_VALIDATION_DT_QE_ET_MIN:-60}"

run_ccw_case 0 "output/ccw.base.2d" "validation/tsr/tmp/ccw.base.2d" "baseline-2d"
run_ccw_case 1 "output/ccw.tsr.2d" "validation/tsr/tmp/ccw.tsr.2d" "tsr-2d"

