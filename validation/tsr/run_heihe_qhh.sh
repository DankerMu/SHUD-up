#!/usr/bin/env bash
set -euo pipefail

sd="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${sd}/_common.sh"

write_project_file() {
  local project="$1"
  local project_file="$2"
  local in_rel="$3"
  local out_rel="$4"

  cat >"$project_file" <<EOF
PRJ	 ${project}
INPATH	 ${in_rel}
OUTPATH	 ${out_rel}
MESH	 ${in_rel}/${project}.sp.mesh
ATT	 ${in_rel}/${project}.sp.att
RIV	 ${in_rel}/${project}.sp.riv
RIVSEG	 ${in_rel}/${project}.sp.rivseg
CALIB	 ${in_rel}/${project}.cfg.calib
PARA	 ${in_rel}/${project}.cfg.para
INIT	 ${in_rel}/${project}.cfg.ic
LC	 ${in_rel}/${project}.para.lc
SOIL	 ${in_rel}/${project}.para.soil
GEOL	 ${in_rel}/${project}.para.geol
FORC	 ${in_rel}/${project}.tsd.forc
LAI	 ${in_rel}/${project}.tsd.lai
MF	 ${in_rel}/${project}.tsd.mf
RL	 ${in_rel}/${project}.tsd.rl
EOF

  # qhh uses the lake module and its lake spatial file is named <project>.lake.sp
  if [[ "$project" == "qhh" ]]; then
    # Insert LAKE after ATT line to match IO.cpp's project-file grammar.
    local tmp
    tmp="$(mktemp)"
    awk -v lake_path="${in_rel}/${project}.lake.sp" '
      {print}
      /^[[:space:]]*ATT[[:space:]]/ {print "LAKE\t " lake_path}
    ' "$project_file" >"$tmp"
    mv "$tmp" "$project_file"
  fi
}

run_case() {
  local project="$1"          # heihe / qhh
  local tsr_flag="$2"         # 0/1
  local out_rel="$3"          # e.g. output/heihe.base
  local tmp_rel="$4"          # e.g. validation/tsr/tmp/heihe.base
  local label="$5"            # baseline/tsr
  local end_override_days="${6:-}" # optional

  local root
  root="$(repo_root)"

  local shud_bin="${root}/shud"
  ensure_executable_file "$shud_bin"

  local input_src="${root}/input/${project}"
  ensure_dir_exists "$input_src"

  local out_dir="${root}/${out_rel}"
  ensure_dir_writable "${root}/output"
  backup_dir_if_exists "$out_dir"
  ensure_dir_writable "$out_dir"

  local tmp_dir="${root}/${tmp_rel}"
  mkdir -p "${root}/validation/tsr/tmp"
  copy_tree "$input_src" "$tmp_dir"

  local cfg="${tmp_dir}/${project}.cfg.para"
  local forc="${tmp_dir}/${project}.tsd.forc"
  [[ -f "$cfg" ]] || die "missing cfg.para in tmp copy: $cfg"
  [[ -f "$forc" ]] || die "missing tsd.forc in tmp copy: $forc"

  if [[ -n "$end_override_days" ]]; then
    upsert_cfg_kv "$cfg" "END" "$end_override_days"
  fi
  upsert_cfg_kv "$cfg" "TERRAIN_RADIATION" "${tsr_flag}"

  # Key fix: compute TSR factor as forcing-interval equivalent.
  if [[ "$tsr_flag" == "1" ]]; then
    upsert_cfg_kv "$cfg" "TSR_FACTOR_MODE" "FORCING_INTERVAL"
    upsert_cfg_kv "$cfg" "TSR_INTEGRATION_STEP_MIN" "60"
  fi

  set_forcing_csv_basepath "$forc" "./${tmp_rel}"

  local project_file="${tmp_dir}/${project}.SHUD"
  write_project_file "$project" "$project_file" "./${tmp_rel}" "./${out_rel}"

  local end_days
  end_days="$(get_cfg_kv "$cfg" "END")"
  info "running ${project} ${label}: TSR=${tsr_flag}, END=${end_days:-?} day(s)"
  (
    cd "$root"
    "${shud_bin}" -p "${project_file}" >"${out_dir}/run.log" 2>&1
  )

  for f in "${project}.rn_h.dat" "${project}.rn_t.dat" "${project}.rn_factor.dat"; do
    [[ -f "${out_dir}/${f}" ]] || die "missing expected output file: ${out_dir}/${f}"
  done
  if [[ -n "${SHUD_WB_DIAG:-}" ]]; then
    [[ -f "${out_dir}/${project}.basinwbfull.dat" ]] || die "missing expected wb output: ${out_dir}/${project}.basinwbfull.dat"
  fi

  info "done: ${project} ${label} -> ${out_rel}"
}

usage() {
  cat <<EOF
Usage: bash validation/tsr/run_heihe_qhh.sh [options]

Options:
  --heihe         Run heihe only
  --qhh           Run qhh only
  --baseline-only Run baseline only (TSR=OFF)
  --tsr-only      Run TSR only (TSR=ON)
  -h, --help      Show this help
EOF
}

run_heihe=1
run_qhh=1
run_baseline=1
run_tsr=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --heihe)
      run_qhh=0
      ;;
    --qhh)
      run_heihe=0
      ;;
    --baseline-only)
      run_tsr=0
      ;;
    --tsr-only)
      run_baseline=0
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
  shift
done

if [[ "$run_heihe" == "1" && "$run_baseline" == "1" ]]; then
  run_case "heihe" 0 "output/heihe.base" "validation/tsr/tmp/heihe.base" "baseline" "9495"
fi
if [[ "$run_heihe" == "1" && "$run_tsr" == "1" ]]; then
  run_case "heihe" 1 "output/heihe.tsr" "validation/tsr/tmp/heihe.tsr" "tsr" "9495"
fi

if [[ "$run_qhh" == "1" && "$run_baseline" == "1" ]]; then
  run_case "qhh" 0 "output/qhh.base" "validation/tsr/tmp/qhh.base" "baseline"
fi
if [[ "$run_qhh" == "1" && "$run_tsr" == "1" ]]; then
  run_case "qhh" 1 "output/qhh.tsr" "validation/tsr/tmp/qhh.tsr" "tsr"
fi
