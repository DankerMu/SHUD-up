#!/usr/bin/env bash
set -euo pipefail

sd="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${sd}/_common.sh"

need_cmd python3

root="$(cd "${sd}/../.." && pwd)"
report_path="${sd}/test_report.txt"

# Keep integration test fast by default. Users can override these env vars if needed.
export SHUD_VALIDATION_END_DAYS="${SHUD_VALIDATION_END_DAYS:-2}"
export SHUD_VALIDATION_DT_QE_ET_MIN="${SHUD_VALIDATION_DT_QE_ET_MIN:-60}"

{
  echo "TSR Integration Test Report"
  echo "Timestamp: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  echo "Repo root: ${root}"
  echo "END days: ${SHUD_VALIDATION_END_DAYS}"
  echo "DT_QE_ET min: ${SHUD_VALIDATION_DT_QE_ET_MIN}"
  echo
} >"${report_path}"

if ! baseline_out="$(
  (run_ccw_case 0 "output/ccw.base.itest" "validation/tsr/tmp/ccw.base.itest" "baseline-itest") 2>&1
)"; then
  {
    echo "Step: baseline-itest: FAIL"
    echo "${baseline_out}"
    echo
    echo "RESULT: FAIL"
  } >>"${report_path}"
  echo "${baseline_out}" >&2
  exit 1
fi
{
  echo "Step: baseline-itest: OK"
  echo
} >>"${report_path}"

if ! tsr_out="$(
  (run_ccw_case 1 "output/ccw.tsr.itest" "validation/tsr/tmp/ccw.tsr.itest" "tsr-itest") 2>&1
)"; then
  {
    echo "Step: tsr-itest: FAIL"
    echo "${tsr_out}"
    echo
    echo "RESULT: FAIL"
  } >>"${report_path}"
  echo "${tsr_out}" >&2
  exit 1
fi
{
  echo "Step: tsr-itest: OK"
  echo
} >>"${report_path}"

SHUD_REPO_ROOT="${root}" SHUD_REPORT_PATH="${report_path}" python3 - <<'PY'
from __future__ import annotations

import math
import os
import struct
from pathlib import Path
from dataclasses import dataclass
import traceback


@dataclass
class DatMeta:
    path: Path
    header_text: str
    start_time: float
    num_var: int
    num_records: int
    dt_min: float
    data_offset: int
    record_size: int


class ValidationError(RuntimeError):
    pass


def report_append(report_path: Path, lines: list[str]) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("a", encoding="utf-8") as f:
        for line in lines:
            f.write(f"{line}\n")


def header_lines(text: str) -> list[str]:
    cleaned = text.replace("\x00", "").strip()
    if not cleaned:
        return []
    return [ln.rstrip() for ln in cleaned.splitlines() if ln.strip()]


def read_meta(path: Path) -> DatMeta:
    if not path.exists():
        raise ValidationError(f"missing .dat file: {path}")
    size = path.stat().st_size
    if size <= 0:
        raise ValidationError(f"empty .dat file: {path}")

    with path.open("rb") as f:
        header_raw = f.read(1024)
        if len(header_raw) != 1024:
            raise ValidationError(f"short header (<1024B) in {path}")
        header_text = header_raw.decode("utf-8", errors="ignore")

        start_time_raw = f.read(8)
        num_var_raw = f.read(8)
        if len(start_time_raw) != 8 or len(num_var_raw) != 8:
            raise ValidationError(f"truncated metadata in {path}")
        start_time = struct.unpack("d", start_time_raw)[0]
        num_var = int(struct.unpack("d", num_var_raw)[0])
        if num_var <= 0:
            raise ValidationError(f"invalid num_var={num_var} in {path}")

    data_offset = 1024 + 16 + 8 * num_var
    record_size = 8 * (1 + num_var)
    if size < data_offset:
        raise ValidationError(f"file too small for header+cols in {path} (size={size}, need>={data_offset})")
    rem = size - data_offset
    if rem % record_size != 0:
        raise ValidationError(f"data section not aligned in {path} (rem={rem}, rec_size={record_size})")
    num_records = rem // record_size
    if num_records <= 0:
        raise ValidationError(f"no records in {path}")

    # Derive dt from first two time stamps (minutes) if possible.
    dt_min = 0.0
    if num_records >= 2:
        with path.open("rb") as f:
            f.seek(data_offset, os.SEEK_SET)
            first = f.read(record_size)
            second = f.read(record_size)
            if len(first) != record_size or len(second) != record_size:
                raise ValidationError(f"truncated records in {path}")
            t0 = struct.unpack_from("d", first, 0)[0]
            t1 = struct.unpack_from("d", second, 0)[0]
            dt_min = t1 - t0
    else:
        dt_min = 0.0

    return DatMeta(
        path=path,
        header_text=header_text,
        start_time=start_time,
        num_var=num_var,
        num_records=num_records,
        dt_min=dt_min,
        data_offset=data_offset,
        record_size=record_size,
    )


def iter_records(meta: DatMeta):
    with meta.path.open("rb") as f:
        f.seek(meta.data_offset, os.SEEK_SET)
        for _ in range(meta.num_records):
            chunk = f.read(meta.record_size)
            if len(chunk) != meta.record_size:
                raise ValidationError(f"truncated record in {meta.path}")
            t = struct.unpack_from("d", chunk, 0)[0]
            values = struct.unpack_from(f"{meta.num_var}d", chunk, 8)
            yield t, values


def scan_stats(meta: DatMeta, abs_cap: float = 1e7) -> dict[str, float]:
    mn = math.inf
    mx = -math.inf
    non_finite = 0
    non_zero = 0
    total = 0
    eps = 1e-12
    for _, values in iter_records(meta):
        for x in values:
            total += 1
            if not math.isfinite(x):
                non_finite += 1
                continue
            if abs(x) > abs_cap:
                raise ValidationError(f"value magnitude too large in {meta.path} (|x|>{abs_cap:g})")
            if abs(x) > eps:
                non_zero += 1
            if x < mn:
                mn = x
            if x > mx:
                mx = x
    if total == 0:
        raise ValidationError(f"no numeric payload in {meta.path}")
    if non_finite != 0:
        raise ValidationError(f"non-finite values detected in {meta.path} (count={non_finite})")
    if non_zero == 0:
        raise ValidationError(f"all-zero payload detected in {meta.path}")
    return {"min": mn, "max": mx, "total": float(total)}


def max_abs_from_constant(meta: DatMeta, c: float) -> float:
    m = 0.0
    for _, values in iter_records(meta):
        for x in values:
            d = abs(x - c)
            if d > m:
                m = d
    return m


def max_abs_diff(a: DatMeta, b: DatMeta) -> float:
    if a.num_var != b.num_var or a.num_records != b.num_records:
        raise ValidationError("shape mismatch for diff")
    m = 0.0
    for (_, va), (_, vb) in zip(iter_records(a), iter_records(b)):
        for xa, xb in zip(va, vb):
            d = abs(xa - xb)
            if d > m:
                m = d
    return m


def max_abs_err_relation(rn_h: DatMeta, rn_t: DatMeta, factor: DatMeta) -> float:
    if rn_h.num_var != rn_t.num_var or rn_h.num_var != factor.num_var:
        raise ValidationError("num_var mismatch in relation")
    if rn_h.num_records != rn_t.num_records or rn_h.num_records != factor.num_records:
        raise ValidationError("record count mismatch in relation")

    m = 0.0
    for (th, vh), (tt, vt), (tf, vf) in zip(iter_records(rn_h), iter_records(rn_t), iter_records(factor)):
        if not (th == tt == tf):
            raise ValidationError("time index mismatch across rn_h/rn_t/factor")
        for h, t, f in zip(vh, vt, vf):
            d = abs(t - (h * f))
            if d > m:
                m = d
    return m


root = Path(os.environ.get("SHUD_REPO_ROOT", os.getcwd()))
report_path = Path(os.environ.get("SHUD_REPORT_PATH", root / "validation" / "tsr" / "test_report.txt"))
expected_end_days = int(os.environ.get("SHUD_VALIDATION_END_DAYS", "2"))
expected_dt_min = float(os.environ.get("SHUD_VALIDATION_DT_QE_ET_MIN", "60"))

base_dir = root / "output" / "ccw.base.itest"
tsr_dir = root / "output" / "ccw.tsr.itest"


def ensure_dir_nonempty(path: Path) -> None:
    if not path.is_dir():
        raise ValidationError(f"missing output directory: {path}")
    if not any(path.iterdir()):
        raise ValidationError(f"output directory is empty: {path}")


def case_metas(out_dir: Path) -> tuple[DatMeta, DatMeta, DatMeta]:
    rn_h = read_meta(out_dir / "ccw.rn_h.dat")
    rn_t = read_meta(out_dir / "ccw.rn_t.dat")
    factor = read_meta(out_dir / "ccw.rn_factor.dat")
    return rn_h, rn_t, factor


def format_meta(meta: DatMeta) -> str:
    dt_part = f"{meta.dt_min:g} min" if meta.dt_min > 0 else "n/a"
    return f"{meta.path.name}: rows={meta.num_records}, cols={meta.num_var}, dt={dt_part}, start={meta.start_time:g}"


lines: list[str] = []

try:
    ensure_dir_nonempty(base_dir)
    ensure_dir_nonempty(tsr_dir)

    base_rn_h, base_rn_t, base_factor = case_metas(base_dir)
    tsr_rn_h, tsr_rn_t, tsr_factor = case_metas(tsr_dir)

    # Header parsing: rows/cols/dt + basic consistency checks.
    for m in (base_rn_h, base_rn_t, base_factor, tsr_rn_h, tsr_rn_t, tsr_factor):
        if m.dt_min <= 0 and m.num_records >= 2:
            raise ValidationError(f"non-positive dt derived from {m.path}")
        if m.num_records != expected_end_days * (24 * 60) // int(expected_dt_min):
            raise ValidationError(
                f"unexpected record count in {m.path.name}: got {m.num_records}, "
                f"expected {expected_end_days * (24 * 60) // int(expected_dt_min)}"
            )
        if abs(m.dt_min - expected_dt_min) > 1e-9 and m.num_records >= 2:
            raise ValidationError(f"unexpected dt in {m.path.name}: got {m.dt_min:g}, expected {expected_dt_min:g}")

    # Baseline and TSR shapes must match.
    if base_factor.num_var != tsr_factor.num_var or base_factor.num_records != tsr_factor.num_records:
        raise ValidationError("baseline vs tsr shape mismatch in rn_factor")

    base_hdr = "\n".join(header_lines(base_factor.header_text))
    tsr_hdr = "\n".join(header_lines(tsr_factor.header_text))
    if "Terrain radiation (TSR): OFF" not in base_hdr:
        raise ValidationError("baseline header check failed: expected 'Terrain radiation (TSR): OFF'")
    if "Terrain radiation (TSR): ON" not in tsr_hdr:
        raise ValidationError("tsr header check failed: expected 'Terrain radiation (TSR): ON'")

    # Numeric sanity: finite + non-zero payload + broad magnitude cap.
    base_rn_h_stats = scan_stats(base_rn_h)
    base_rn_t_stats = scan_stats(base_rn_t)
    base_factor_stats = scan_stats(base_factor)
    tsr_rn_h_stats = scan_stats(tsr_rn_h)
    tsr_rn_t_stats = scan_stats(tsr_rn_t)
    tsr_factor_stats = scan_stats(tsr_factor)

    # Physics relation: rn_t ≈ rn_h * factor.
    err_base = max_abs_err_relation(base_rn_h, base_rn_t, base_factor)
    err_tsr = max_abs_err_relation(tsr_rn_h, tsr_rn_t, tsr_factor)
    if err_base > 1e-9:
        raise ValidationError(f"baseline rn_t != rn_h * factor (max abs err={err_base:.3e})")
    if err_tsr > 1e-9:
        raise ValidationError(f"tsr rn_t != rn_h * factor (max abs err={err_tsr:.3e})")

    # Factor expectations.
    base_factor_maxdev = max_abs_from_constant(base_factor, 1.0)
    if base_factor_maxdev > 1e-9:
        raise ValidationError(f"baseline factor not ~1 (max |f-1|={base_factor_maxdev:.3e})")

    tsr_mn = tsr_factor_stats["min"]
    tsr_mx = tsr_factor_stats["max"]
    if tsr_mx <= 1.001:
        raise ValidationError(f"tsr factor shows no variation (max={tsr_mx:.6f})")
    if tsr_mn < -1e-12:
        raise ValidationError(f"tsr factor has negative values (min={tsr_mn:.6f})")
    if tsr_mx > 5.0 + 1e-6:
        raise ValidationError(f"tsr factor exceeds RAD_FACTOR_CAP=5 (max={tsr_mx:.6f})")

    # TSR=ON vs OFF must differ.
    diff_factor = max_abs_diff(base_factor, tsr_factor)
    diff_rn_t = max_abs_diff(base_rn_t, tsr_rn_t)
    if diff_factor <= 1e-6 and diff_rn_t <= 1e-6:
        raise ValidationError("baseline vs tsr outputs are unexpectedly identical")

    lines.append("Validation: .dat header/shape")
    lines.append(f"- expected END={expected_end_days} day(s), DT_QE_ET={expected_dt_min:g} min")
    lines.append(f"- baseline {format_meta(base_factor)}")
    lines.append(f"- tsr      {format_meta(tsr_factor)}")
    lines.append("")
    lines.append("Validation: numeric sanity (min/max)")
    lines.append(f"- baseline rn_h  range=[{base_rn_h_stats['min']:.6g}, {base_rn_h_stats['max']:.6g}]")
    lines.append(f"- baseline rn_t  range=[{base_rn_t_stats['min']:.6g}, {base_rn_t_stats['max']:.6g}]")
    lines.append(f"- baseline factor range=[{base_factor_stats['min']:.6g}, {base_factor_stats['max']:.6g}]")
    lines.append(f"- tsr rn_h       range=[{tsr_rn_h_stats['min']:.6g}, {tsr_rn_h_stats['max']:.6g}]")
    lines.append(f"- tsr rn_t       range=[{tsr_rn_t_stats['min']:.6g}, {tsr_rn_t_stats['max']:.6g}]")
    lines.append(f"- tsr factor     range=[{tsr_factor_stats['min']:.6g}, {tsr_factor_stats['max']:.6g}]")
    lines.append("")
    lines.append("Validation: physical relation & diff")
    lines.append(f"- baseline: max |factor-1|={base_factor_maxdev:.3e}, rn_t=rn_h*factor err={err_base:.3e}")
    lines.append(f"- tsr: factor range=[{tsr_mn:.6f}, {tsr_mx:.6f}], rn_t=rn_h*factor err={err_tsr:.3e}")
    lines.append(f"- baseline vs tsr: max |factor_base-factor_tsr|={diff_factor:.3e}")
    lines.append(f"- baseline vs tsr: max |rn_t_base-rn_t_tsr|={diff_rn_t:.3e}")
    lines.append("")
    lines.append("RESULT: PASS")
    report_append(report_path, lines)
    print("PASS: TSR integration validation")
except ValidationError as e:
    lines.append("RESULT: FAIL")
    lines.append(f"Reason: {e}")
    report_append(report_path, lines)
    raise SystemExit(1)
except Exception:
    lines.append("RESULT: FAIL")
    lines.append("Reason: unexpected exception while validating outputs")
    lines.append(traceback.format_exc().rstrip())
    report_append(report_path, lines)
    raise SystemExit(1)
PY
