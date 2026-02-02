#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "validation" / "tsr" / "py"))
from shud_reader import DatReader  # noqa: E402


class SummarizeError(RuntimeError):
    pass


def _read_matrix_fast(path: Path) -> tuple[np.ndarray, np.ndarray, DatReader]:
    r = DatReader(path)
    m = r.meta
    dtype = np.dtype(m.endianness + "f8")
    total = m.num_records * (1 + m.num_var)
    with m.path.open("rb") as f:
        f.seek(m.data_offset, 0)
        data = np.fromfile(f, dtype=dtype, count=total)
    if data.size != total:
        raise SummarizeError(f"truncated read: {path} (got {data.size}, expected {total})")
    data = data.reshape((m.num_records, 1 + m.num_var))
    t_min = data[:, 0].astype(float)
    values = data[:, 1:].astype(float)
    return t_min, values, r


def _read_mesh_areas(mesh_path: Path) -> np.ndarray:
    lines = mesh_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise SummarizeError(f"empty mesh file: {mesh_path}")
    num_ele = int(lines[0].split()[0])

    elements: list[tuple[int, int, int, int]] = []
    for k in range(num_ele):
        cols = lines[2 + k].split()
        if len(cols) < 4:
            raise SummarizeError(f"invalid element row in {mesh_path}: {lines[2+k]!r}")
        elements.append((int(cols[0]), int(cols[1]), int(cols[2]), int(cols[3])))

    node_hdr_idx = 2 + num_ele
    num_node = int(lines[node_hdr_idx].split()[0])
    node_start = node_hdr_idx + 2

    x = np.zeros(num_node + 1)
    y = np.zeros(num_node + 1)
    for k in range(num_node):
        cols = lines[node_start + k].split()
        if len(cols) < 3:
            raise SummarizeError(f"invalid node row in {mesh_path}: {lines[node_start+k]!r}")
        nid = int(cols[0])
        x[nid] = float(cols[1])
        y[nid] = float(cols[2])

    areas = np.zeros(num_ele)
    for eid, n1, n2, n3 in elements:
        x1, y1 = x[n1], y[n1]
        x2, y2 = x[n2], y[n2]
        x3, y3 = x[n3], y[n3]
        areas[eid - 1] = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    return areas


def _rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(math.sqrt(float(np.nanmean(x * x))))


def _pct_abs(x: np.ndarray, p: float) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    return float(np.nanpercentile(np.abs(x), p))


def _parse_runtime_seconds(run_log: Path) -> Optional[float]:
    if not run_log.exists():
        return None
    text = run_log.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"Time used by model:\s*([0-9.+-eE]+)\s*seconds", text)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


@dataclass(frozen=True)
class BasinResidStats:
    n: int
    max_abs_mm: float
    rms_mm: float
    p50_abs_mm: float
    p95_abs_mm: float
    p99_abs_mm: float
    mean_mm: float
    sum_mm: float
    sum_abs_mm: float
    peak_t_min: float
    peak_day: float
    peak_resid_mm: float


def basin_resid_stats_mm(*, basinwb_path: Path, mesh_area_m2: float) -> BasinResidStats:
    t_min, vals, _ = _read_matrix_fast(basinwb_path)
    if vals.shape[1] < 9:
        raise SummarizeError(f"unexpected basinwb var count: {basinwb_path} (num_var={vals.shape[1]})")
    resid_m3 = vals[:, 8]
    resid_mm = resid_m3 / mesh_area_m2 * 1000.0
    abs_mm = np.abs(resid_mm)
    finite = np.isfinite(resid_mm)
    if not np.any(finite):
        raise SummarizeError(f"no finite residuals in {basinwb_path}")

    i_peak = int(np.nanargmax(abs_mm))
    peak_t = float(t_min[i_peak])
    peak_day = peak_t / 1440.0
    peak_resid = float(resid_mm[i_peak])

    return BasinResidStats(
        n=int(resid_mm.size),
        max_abs_mm=float(np.nanmax(abs_mm)),
        rms_mm=_rms(resid_mm),
        p50_abs_mm=_pct_abs(resid_mm, 50.0),
        p95_abs_mm=_pct_abs(resid_mm, 95.0),
        p99_abs_mm=_pct_abs(resid_mm, 99.0),
        mean_mm=float(np.nanmean(resid_mm)),
        sum_mm=float(np.nansum(resid_mm)),
        sum_abs_mm=float(np.nansum(abs_mm)),
        peak_t_min=peak_t,
        peak_day=float(peak_day),
        peak_resid_mm=peak_resid,
    )


@dataclass(frozen=True)
class BasinTotalsMm:
    dS_mm: float
    P_mm: float
    ET_mm: float
    Qout_mm: float
    resid_mm: float


def basin_totals_mm(*, basinwb_path: Path, mesh_area_m2: float) -> BasinTotalsMm:
    _, vals, _ = _read_matrix_fast(basinwb_path)
    total_m3 = vals.sum(axis=0)
    mm = total_m3 / mesh_area_m2 * 1000.0
    return BasinTotalsMm(
        dS_mm=float(mm[0]),
        P_mm=float(mm[1]),
        ET_mm=float(mm[2]),
        Qout_mm=float(mm[3]),
        resid_mm=float(mm[8]),
    )


@dataclass(frozen=True)
class AbsDistStats:
    max_abs: float
    p95_abs: float
    p99_abs: float


def abs_dist_stats(values: np.ndarray) -> AbsDistStats:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return AbsDistStats(max_abs=float("nan"), p95_abs=float("nan"), p99_abs=float("nan"))
    av = np.abs(values)
    return AbsDistStats(
        max_abs=float(np.max(av)),
        p95_abs=float(np.percentile(av, 95.0)),
        p99_abs=float(np.percentile(av, 99.0)),
    )


def flatten_dat_values(path: Path) -> np.ndarray:
    _, vals, _ = _read_matrix_fast(path)
    return vals.reshape(-1)


def corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 2:
        return float("nan")
    aa = a[m]
    bb = b[m]
    aa = aa - float(np.mean(aa))
    bb = bb - float(np.mean(bb))
    denom = float(np.sqrt(np.sum(aa * aa) * np.sum(bb * bb)))
    if denom <= 0.0:
        return float("nan")
    return float(np.sum(aa * bb) / denom)


def main() -> int:
    ap = argparse.ArgumentParser(description="Summarize docs/water_balance_verification.md numbers from outputs.")
    ap.add_argument("--mesh", type=Path, default=Path("input/ccw/ccw.sp.mesh"))
    ap.add_argument("--base", type=Path, default=Path("output/ccw.base"))
    ap.add_argument("--tsr", type=Path, default=Path("output/ccw.tsr"))
    ap.add_argument("--base-2d", type=Path, default=Path("output/ccw.base.2d"))
    ap.add_argument("--tsr-2d", type=Path, default=Path("output/ccw.tsr.2d"))
    ap.add_argument("--tsr-ms5", type=Path, default=Path("output/ccw.tsr.ms5"))
    ap.add_argument("--tsr-ms2", type=Path, default=Path("output/ccw.tsr.ms2"))
    ap.add_argument("--wb-be", type=Path, default=Path("output/ccw.tsr.wb.be"))
    ap.add_argument("--wb-trapz", type=Path, default=Path("output/ccw.tsr.wb.trapz"))
    ap.add_argument("--wb-trapzquad", type=Path, default=Path("output/ccw.tsr.wb.trapzquad"))
    ap.add_argument("--tsr-prcp5", type=Path, default=Path("output/ccw.tsr.prcp5"))
    args = ap.parse_args()

    areas = _read_mesh_areas(args.mesh)
    A = float(areas.sum())
    if not (A > 0.0):
        raise SummarizeError(f"invalid mesh area: {A}")

    def s_short(out_dir: Path) -> BasinResidStats:
        return basin_resid_stats_mm(basinwb_path=out_dir / "ccw.basinwbfull.dat", mesh_area_m2=A)

    def s_full(out_dir: Path) -> BasinResidStats:
        return basin_resid_stats_mm(basinwb_path=out_dir / "ccw.basinwbfull.dat", mesh_area_m2=A)

    def runtime(out_dir: Path) -> Optional[float]:
        return _parse_runtime_seconds(out_dir / "run.log")

    # 2-day short run
    base2 = s_short(args.base_2d)
    tsr2 = s_short(args.tsr_2d)
    ele_a_base2 = abs_dist_stats(flatten_dat_values(args.base_2d / "ccw.elewb3_resid.dat") * 1000.0)
    ele_a_tsr2 = abs_dist_stats(flatten_dat_values(args.tsr_2d / "ccw.elewb3_resid.dat") * 1000.0)
    ele_b_base2 = abs_dist_stats(flatten_dat_values(args.base_2d / "ccw.elewb3_budget_resid.dat") * 1000.0)
    ele_b_tsr2 = abs_dist_stats(flatten_dat_values(args.tsr_2d / "ccw.elewb3_budget_resid.dat") * 1000.0)

    # Full run (TSR=ON)
    tsr_full = s_full(args.tsr)
    tsr_tot = basin_totals_mm(basinwb_path=args.tsr / "ccw.basinwbfull.dat", mesh_area_m2=A)

    # Full run per-element residual stats (meters)
    ele_full_a = abs_dist_stats(flatten_dat_values(args.tsr / "ccw.elewb3_resid.dat"))
    ele_full_b = abs_dist_stats(flatten_dat_values(args.tsr / "ccw.elewb3_budget_resid.dat"))

    # Max-step experiments (TSR=ON)
    ms10 = tsr_full
    ms5 = s_full(args.tsr_ms5)
    ms2 = s_full(args.tsr_ms2)

    # Integrators experiments
    wb_be = s_full(args.wb_be)
    wb_tr = s_full(args.wb_trapz)
    wb_quad = basin_resid_stats_mm(basinwb_path=args.wb_trapzquad / "ccw.basinwbfull_quad.dat", mesh_area_m2=A)

    # PRCP scaling
    prcp1 = basin_totals_mm(basinwb_path=args.tsr / "ccw.basinwbfull.dat", mesh_area_m2=A)
    prcp5 = basin_totals_mm(basinwb_path=args.tsr_prcp5 / "ccw.basinwbfull.dat", mesh_area_m2=A)
    prcp1_stats = s_full(args.tsr)
    prcp5_stats = s_full(args.tsr_prcp5)

    ele_prcp1_a = abs_dist_stats(flatten_dat_values(args.tsr / "ccw.elewb3_resid.dat") * 1000.0)
    ele_prcp1_b = abs_dist_stats(flatten_dat_values(args.tsr / "ccw.elewb3_budget_resid.dat") * 1000.0)
    ele_prcp5_a = abs_dist_stats(flatten_dat_values(args.tsr_prcp5 / "ccw.elewb3_resid.dat") * 1000.0)
    ele_prcp5_b = abs_dist_stats(flatten_dat_values(args.tsr_prcp5 / "ccw.elewb3_budget_resid.dat") * 1000.0)

    # Correlations (TSR full run, daily intervals)
    _, wb_vals, _ = _read_matrix_fast(args.tsr / "ccw.basinwbfull.dat")
    resid_mm = wb_vals[:, 8] / A * 1000.0
    abs_resid_mm = np.abs(resid_mm)
    P_mm = wb_vals[:, 1] / A * 1000.0
    dS_mm = wb_vals[:, 0] / A * 1000.0
    Qout_mm = wb_vals[:, 3] / A * 1000.0
    corr_P = corr(abs_resid_mm, P_mm)
    corr_dS = corr(abs_resid_mm, dS_mm)
    corr_Q = corr(abs_resid_mm, Qout_mm)

    # Contribution concentration: Top-K by |resid|
    abs_sorted = np.sort(abs_resid_mm[np.isfinite(abs_resid_mm)])[::-1]
    total_abs = float(np.sum(abs_sorted))
    top10_pct = float(np.sum(abs_sorted[:10]) / total_abs * 100.0) if total_abs > 0.0 and abs_sorted.size >= 10 else float("nan")
    top100_pct = float(np.sum(abs_sorted[:100]) / total_abs * 100.0) if total_abs > 0.0 and abs_sorted.size >= 100 else float("nan")

    # Outlier in PRCP5: max |elewb3_resid| (mm) + location
    t_min_p5, ele_p5_vals, _ = _read_matrix_fast(args.tsr_prcp5 / "ccw.elewb3_resid.dat")
    ele_p5_mm = ele_p5_vals * 1000.0
    abs_ele_p5 = np.abs(ele_p5_mm)
    idx_flat = int(np.nanargmax(abs_ele_p5))
    k_time = int(idx_flat // abs_ele_p5.shape[1])
    k_ele = int(idx_flat % abs_ele_p5.shape[1])
    outlier_mm = float(ele_p5_mm[k_time, k_ele])
    outlier_day = float(t_min_p5[k_time] / 1440.0)
    outlier_ele_id_1based = k_ele + 1

    # Basin residual at that day for prcp5
    t_min_b_p5, wb_p5, _ = _read_matrix_fast(args.tsr_prcp5 / "ccw.basinwbfull.dat")
    resid_p5_mm = wb_p5[:, 8] / A * 1000.0
    P_p5_mm = wb_p5[:, 1] / A * 1000.0
    # match by exact t_min
    basin_i = int(np.where(t_min_b_p5 == t_min_p5[k_time])[0][0]) if np.any(t_min_b_p5 == t_min_p5[k_time]) else -1
    basin_resid_on_outlier_day = float(resid_p5_mm[basin_i]) if basin_i >= 0 else float("nan")
    basin_P_on_outlier_day = float(P_p5_mm[basin_i]) if basin_i >= 0 else float("nan")

    # Print a machine-greppable summary (also easy to paste into docs).
    print(f"Mesh area A = {A:.3f} m^2")

    print("\n[ccw 2-day short run] basin residual stats (mm / 60min interval)")
    for name, st in [("BASE", base2), ("TSR", tsr2)]:
        print(
            f"{name}: max|resid|={st.max_abs_mm:.7g} mm @ t={st.peak_t_min:.0f} min, "
            f"mean={st.mean_mm:+.7g} mm, rms={st.rms_mm:.7g} mm, "
            f"p95|·|={st.p95_abs_mm:.7g} mm, p99|·|={st.p99_abs_mm:.7g} mm, Σresid={st.sum_mm:+.7g} mm"
        )

    print("\n[ccw 2-day short run] per-element |resid| stats (mm / 60min interval)")
    print(f"A (elewb3_resid): BASE max={ele_a_base2.max_abs:.7g} p95={ele_a_base2.p95_abs:.7g} p99={ele_a_base2.p99_abs:.7g}")
    print(f"A (elewb3_resid): TSR  max={ele_a_tsr2.max_abs:.7g} p95={ele_a_tsr2.p95_abs:.7g} p99={ele_a_tsr2.p99_abs:.7g}")
    print(f"B (elewb3_budget_resid): BASE max={ele_b_base2.max_abs:.7g} p95={ele_b_base2.p95_abs:.7g} p99={ele_b_base2.p99_abs:.7g}")
    print(f"B (elewb3_budget_resid): TSR  max={ele_b_tsr2.max_abs:.7g} p95={ele_b_tsr2.p95_abs:.7g} p99={ele_b_tsr2.p99_abs:.7g}")

    print("\n[ccw full TSR=ON] basinwbfull resid stats (mm / 1440min interval)")
    print(
        f"n={tsr_full.n}, max|resid|={tsr_full.max_abs_mm:.7g} mm @ Day {tsr_full.peak_day:.0f} "
        f"(t={tsr_full.peak_t_min:.0f} min, resid={tsr_full.peak_resid_mm:+.7g} mm), "
        f"mean={tsr_full.mean_mm:+.7g} mm, rms={tsr_full.rms_mm:.7g} mm, "
        f"p50|·|={tsr_full.p50_abs_mm:.7g} mm, p95|·|={tsr_full.p95_abs_mm:.7g} mm, p99|·|={tsr_full.p99_abs_mm:.7g} mm, "
        f"Σresid={tsr_full.sum_mm:+.7g} mm"
    )
    print(
        f"Totals (mm): ΣP={tsr_tot.P_mm:.5g}, ΣET={tsr_tot.ET_mm:.5g}, ΣQout={tsr_tot.Qout_mm:.5g}, "
        f"ΔS={tsr_tot.dS_mm:.5g}, Σresid={tsr_tot.resid_mm:+.5g}"
    )

    print("\n[ccw full TSR=ON] per-element |resid| stats (meters / 1440min interval)")
    print(f"A (elewb3_resid): max={ele_full_a.max_abs:.7g} p95={ele_full_a.p95_abs:.7g} p99={ele_full_a.p99_abs:.7g}")
    print(f"B (elewb3_budget_resid): max={ele_full_b.max_abs:.7g} p95={ele_full_b.p95_abs:.7g} p99={ele_full_b.p99_abs:.7g}")

    print("\n[ccw TSR=ON] MAX_SOLVER_STEP convergence (mm / day interval)")
    for label, st, rt in [
        ("10", ms10, runtime(args.tsr)),
        ("5", ms5, runtime(args.tsr_ms5)),
        ("2", ms2, runtime(args.tsr_ms2)),
    ]:
        print(
            f"MS{label}: max|resid|={st.max_abs_mm:.7g} mm, rms={st.rms_mm:.7g} mm, "
            f"p95={st.p95_abs_mm:.7g} mm, p99={st.p99_abs_mm:.7g} mm, Σresid={st.sum_mm:+.7g} mm, "
            f"Σ|resid|={st.sum_abs_mm:.7g} mm, peak day={st.peak_day:.0f}, runtime_s={rt if rt is not None else float('nan'):.3f}"
        )

    print("\n[ccw TSR=ON] water-balance integrators (mm / day interval)")
    for label, st, rt in [
        ("BE(sampled)", wb_be, runtime(args.wb_be)),
        ("TRAPZ(sampled)", wb_tr, runtime(args.wb_trapz)),
        ("QUAD(internal-step)", wb_quad, runtime(args.wb_trapzquad)),
    ]:
        print(
            f"{label}: max|resid|={st.max_abs_mm:.7g} mm, rms={st.rms_mm:.7g} mm, p95={st.p95_abs_mm:.7g} mm, "
            f"p99={st.p99_abs_mm:.7g} mm, Σresid={st.sum_mm:+.7g} mm, Σ|resid|={st.sum_abs_mm:.7g} mm, "
            f"peak day={st.peak_day:.0f}, runtime_s={rt if rt is not None else float('nan'):.3f}"
        )

    print("\n[ccw TSR=ON] PRCP scaling (TS_PRCP)")
    print(
        f"TS_PRCP=1: ΣP={prcp1.P_mm:.5g} mm, max|resid|={prcp1_stats.max_abs_mm:.7g} mm, rms={prcp1_stats.rms_mm:.7g} mm, "
        f"p95={prcp1_stats.p95_abs_mm:.7g} mm, p99={prcp1_stats.p99_abs_mm:.7g} mm, Σresid={prcp1_stats.sum_mm:+.7g} mm, Σ|resid|={prcp1_stats.sum_abs_mm:.7g} mm"
    )
    print(
        f"TS_PRCP=5: ΣP={prcp5.P_mm:.5g} mm, max|resid|={prcp5_stats.max_abs_mm:.7g} mm, rms={prcp5_stats.rms_mm:.7g} mm, "
        f"p95={prcp5_stats.p95_abs_mm:.7g} mm, p99={prcp5_stats.p99_abs_mm:.7g} mm, Σresid={prcp5_stats.sum_mm:+.7g} mm, Σ|resid|={prcp5_stats.sum_abs_mm:.7g} mm"
    )

    print("\n[ccw TSR full run] correlation + concentration")
    print(f"corr(|resid|, P) = {corr_P:.6g}")
    print(f"corr(|resid|, dS_total) = {corr_dS:.6g}")
    print(f"corr(|resid|, Qout) = {corr_Q:.6g}")
    print(f"Top10 contribution by Σ|resid| = {top10_pct:.4g}%")
    print(f"Top100 contribution by Σ|resid| = {top100_pct:.4g}%")

    print("\n[ccw TS_PRCP=5] worst per-element A residual (elewb3_resid)")
    print(f"max |elewb3_resid| = {abs(float(outlier_mm)):.6g} mm at day={outlier_day:.0f}, element={outlier_ele_id_1based}, value={outlier_mm:+.6g} mm")
    print(f"same day basin resid ≈ {basin_resid_on_outlier_day:+.6g} mm, P ≈ {basin_P_on_outlier_day:.6g} mm/day")

    print("\n[ccw TS_PRCP] per-element |resid| distribution (mm / day interval)")
    print(
        f"A (elewb3_resid): TS_PRCP=1 p95={ele_prcp1_a.p95_abs:.6g} p99={ele_prcp1_a.p99_abs:.6g} max={ele_prcp1_a.max_abs:.6g}"
    )
    print(
        f"A (elewb3_resid): TS_PRCP=5 p95={ele_prcp5_a.p95_abs:.6g} p99={ele_prcp5_a.p99_abs:.6g} max={ele_prcp5_a.max_abs:.6g}"
    )
    print(
        f"B (elewb3_budget_resid): TS_PRCP=1 p95={ele_prcp1_b.p95_abs:.6g} p99={ele_prcp1_b.p99_abs:.6g} max={ele_prcp1_b.max_abs:.6g}"
    )
    print(
        f"B (elewb3_budget_resid): TS_PRCP=5 p95={ele_prcp5_b.p95_abs:.6g} p99={ele_prcp5_b.p99_abs:.6g} max={ele_prcp5_b.max_abs:.6g}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
