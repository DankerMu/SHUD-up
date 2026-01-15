#!/usr/bin/env bash
set -euo pipefail

sd="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${sd}/_common.sh"

run_ccw_case 1 "output/ccw.tsr" "validation/tsr/tmp/ccw.tsr" "tsr"
