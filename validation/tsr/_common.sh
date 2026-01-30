#!/usr/bin/env bash
set -euo pipefail

die() {
  echo "ERROR: $*" >&2
  exit 1
}

info() {
  echo "==> $*" >&2
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "missing required command: $1"
}

script_dir() {
  cd "$(dirname "${BASH_SOURCE[0]}")" && pwd
}

repo_root() {
  local sd
  sd="$(script_dir)"
  cd "${sd}/../.." && pwd
}

ensure_executable_file() {
  local path="$1"
  [[ -f "$path" ]] || die "required file missing: $path"
  [[ -x "$path" ]] || die "required file is not executable: $path"
}

ensure_dir_exists() {
  local path="$1"
  [[ -d "$path" ]] || die "required directory missing: $path"
}

ensure_dir_writable() {
  local path="$1"
  mkdir -p "$path" || die "failed to create directory: $path"
  local probe="${path}/.writetest.$$"
  : >"$probe" 2>/dev/null || die "directory not writable: $path"
  rm -f "$probe"
}

backup_dir_if_exists() {
  local path="$1"
  if [[ -d "$path" ]]; then
    if find "$path" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
      local ts
      ts="$(date +"%Y%m%d-%H%M%S")"
      local backup="${path}.bak.${ts}"
      info "output dir exists; moving to ${backup}"
      mv "$path" "$backup"
    else
      rmdir "$path" 2>/dev/null || true
    fi
  fi
}

copy_tree() {
  local src="$1"
  local dst="$2"
  ensure_dir_exists "$src"
  rm -rf "$dst"
  mkdir -p "$dst"
  if command -v rsync >/dev/null 2>&1; then
    rsync -a "${src}/" "${dst}/"
  else
    cp -R "${src}/." "${dst}/"
  fi
}

upsert_cfg_kv() {
  local file="$1"
  local key="$2"
  local value="$3"
  local tmp
  tmp="$(mktemp)"
  awk -v key="$key" -v value="$value" '
    BEGIN{found=0}
    /^[[:space:]]*#/ {print; next}
    /^[[:space:]]*$/ {print; next}
    {
      k=$1
      if (toupper(k)==toupper(key)) {
        print key "\t" value
        found=1
      } else {
        print
      }
      next
    }
    END{
      if(found==0){
        print key "\t" value
      }
    }
  ' "$file" >"$tmp"
  mv "$tmp" "$file"
}

get_cfg_kv() {
  local file="$1"
  local key="$2"
  awk -v key="$key" '
    BEGIN{IGNORECASE=1}
    /^[[:space:]]*#/ {next}
    /^[[:space:]]*$/ {next}
    {
      k=$1
      v=$2
      if (toupper(k)==toupper(key)) {
        print v
        exit
      }
    }
  ' "$file"
}

set_forcing_csv_basepath() {
  local forc_list="$1"
  local new_path="$2"
  [[ -f "$forc_list" ]] || die "missing forcing list file: $forc_list"
  local tmp
  tmp="$(mktemp)"
  awk -v new_path="$new_path" 'NR==2{print new_path; next} {print}' "$forc_list" >"$tmp"
  mv "$tmp" "$forc_list"
}

write_ccw_project_file() {
  local project_file="$1"
  local in_rel="$2"
  local out_rel="$3"
  cat >"$project_file" <<EOF
PRJ	 ccw
INPATH	 ${in_rel}
OUTPATH	 ${out_rel}
MESH	 ${in_rel}/ccw.sp.mesh
ATT	 ${in_rel}/ccw.sp.att
RIV	 ${in_rel}/ccw.sp.riv
RIVSEG	 ${in_rel}/ccw.sp.rivseg
CALIB	 ${in_rel}/ccw.cfg.calib
PARA	 ${in_rel}/ccw.cfg.para
INIT	 ${in_rel}/ccw.cfg.ic
LC	 ${in_rel}/ccw.para.lc
SOIL	 ${in_rel}/ccw.para.soil
GEOL	 ${in_rel}/ccw.para.geol
FORC	 ${in_rel}/ccw.tsd.forc
LAI	 ${in_rel}/ccw.tsd.lai
MF	 ${in_rel}/ccw.tsd.mf
RL	 ${in_rel}/ccw.tsd.rl
EOF
}

run_ccw_case() {
  local tsr_flag="$1"          # 0/1
  local out_rel="$2"           # e.g. output/ccw.base
  local tmp_rel="$3"           # e.g. validation/tsr/tmp/ccw.base
  local label="$4"             # e.g. baseline/tsr

  local root
  root="$(repo_root)"

  local shud_bin="${root}/shud"
  ensure_executable_file "$shud_bin"

  local input_src="${root}/input/ccw"
  ensure_dir_exists "$input_src"

  local out_dir="${root}/${out_rel}"
  ensure_dir_writable "${root}/output"
  backup_dir_if_exists "$out_dir"
  ensure_dir_writable "$out_dir"

  local tmp_dir="${root}/${tmp_rel}"
  mkdir -p "${root}/validation/tsr/tmp"
  copy_tree "$input_src" "$tmp_dir"

  local cfg="${tmp_dir}/ccw.cfg.para"
  local forc="${tmp_dir}/ccw.tsd.forc"
  [[ -f "$cfg" ]] || die "missing cfg.para in tmp copy: $cfg"
  [[ -f "$forc" ]] || die "missing tsd.forc in tmp copy: $forc"

  # Default: keep values from input/ccw/ccw.cfg.para (copied into tmp_dir).
  # Optional overrides for fast validation:
  if [[ -n "${SHUD_VALIDATION_END_DAYS:-}" ]]; then
    upsert_cfg_kv "$cfg" "END" "${SHUD_VALIDATION_END_DAYS}"
  fi
  if [[ -n "${SHUD_VALIDATION_DT_QE_ET_MIN:-}" ]]; then
    upsert_cfg_kv "$cfg" "DT_QE_ET" "${SHUD_VALIDATION_DT_QE_ET_MIN}"
  fi
  upsert_cfg_kv "$cfg" "TERRAIN_RADIATION" "${tsr_flag}"

  set_forcing_csv_basepath "$forc" "./${tmp_rel}"

  local project_file="${tmp_dir}/ccw.SHUD"
  write_ccw_project_file "$project_file" "./${tmp_rel}" "./${out_rel}"

  local end_days
  local dt_qe_et_min
  end_days="$(get_cfg_kv "$cfg" "END")"
  dt_qe_et_min="$(get_cfg_kv "$cfg" "DT_QE_ET")"
  info "running ${label}: TSR=${tsr_flag}, END=${end_days:-?} day(s), DT_QE_ET=${dt_qe_et_min:-?} min"
  (
    cd "$root"
    "${shud_bin}" -p "${project_file}" >"${out_dir}/run.log" 2>&1
  )

  for f in "ccw.rn_h.dat" "ccw.rn_t.dat" "ccw.rn_factor.dat"; do
    [[ -f "${out_dir}/${f}" ]] || die "missing expected output file: ${out_dir}/${f}"
  done

  info "done: ${label} -> ${out_rel}"
}
