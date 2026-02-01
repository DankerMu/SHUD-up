#!/usr/bin/env bash
set -euo pipefail

sd="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${sd}/_common.sh"

need_cmd "${sd}/../../shud"

echo "==> Full-run water-balance diagnostics with different integrators (TSR=ON)"
echo "    - Backward-Euler (default sampling integration)"
echo "    - Trapezoidal integration (SHUD_WB_DIAG_TRAPZ=1)"
echo "    - CVODE internal-step integration for basin fluxes (SHUD_WB_DIAG_QUAD=1; slow)"

export SHUD_WB_DIAG=1

export SHUD_WB_DIAG_TRAPZ=0
export SHUD_WB_DIAG_QUAD=0
run_ccw_case 1 "output/ccw.tsr.wb.be" "validation/tsr/tmp/ccw.tsr.wb.be" "tsr-wb-be"

export SHUD_WB_DIAG_TRAPZ=1
export SHUD_WB_DIAG_QUAD=0
run_ccw_case 1 "output/ccw.tsr.wb.trapz" "validation/tsr/tmp/ccw.tsr.wb.trapz" "tsr-wb-trapz"

export SHUD_WB_DIAG_TRAPZ=1
export SHUD_WB_DIAG_QUAD=1
run_ccw_case 1 "output/ccw.tsr.wb.trapzquad" "validation/tsr/tmp/ccw.tsr.wb.trapzquad" "tsr-wb-trapzquad"

echo "OK: outputs:"
echo "- output/ccw.tsr.wb.be"
echo "- output/ccw.tsr.wb.trapz"
echo "- output/ccw.tsr.wb.trapzquad (includes *.basinwbfull_quad.dat)"

