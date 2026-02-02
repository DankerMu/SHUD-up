#!/usr/bin/env bash
set -euo pipefail

sd="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${sd}/_common.sh"

root="$(repo_root)"

shud_bin="${root}/shud"
ensure_executable_file "$shud_bin"

input_src="${root}/input/ccw"
ensure_dir_exists "$input_src"

out_rel="output/ccw.tsr.prcp5"
tmp_rel="validation/tsr/tmp/ccw.tsr.prcp5"

out_dir="${root}/${out_rel}"
ensure_dir_writable "${root}/output"
backup_dir_if_exists "$out_dir"
ensure_dir_writable "$out_dir"

tmp_dir="${root}/${tmp_rel}"
mkdir -p "${root}/validation/tsr/tmp"
copy_tree "$input_src" "$tmp_dir"

cfg="${tmp_dir}/ccw.cfg.para"
calib="${tmp_dir}/ccw.cfg.calib"
forc="${tmp_dir}/ccw.tsd.forc"
[[ -f "$cfg" ]] || die "missing cfg.para in tmp copy: $cfg"
[[ -f "$calib" ]] || die "missing cfg.calib in tmp copy: $calib"
[[ -f "$forc" ]] || die "missing tsd.forc in tmp copy: $forc"

upsert_cfg_kv "$cfg" "TERRAIN_RADIATION" "1"

# Pressure test: scale precipitation forcing by 5.
upsert_cfg_kv "$calib" "TS_PRCP" "5"

set_forcing_csv_basepath "$forc" "./${tmp_rel}"

project_file="${tmp_dir}/ccw.SHUD"
write_ccw_project_file "$project_file" "./${tmp_rel}" "./${out_rel}"

info "running ccw prcp5: TSR=1, TS_PRCP=5, END=$(get_cfg_kv "$cfg" "END") day(s), DT_QE_ET=$(get_cfg_kv "$cfg" "DT_QE_ET") min"
(
  cd "$root"
  "${shud_bin}" -p "${project_file}" >"${out_dir}/run.log" 2>&1
)

for f in "ccw.rn_h.dat" "ccw.rn_t.dat" "ccw.rn_factor.dat"; do
  [[ -f "${out_dir}/${f}" ]] || die "missing expected output file: ${out_dir}/${f}"
done

info "done: ccw prcp5 -> ${out_rel}"

