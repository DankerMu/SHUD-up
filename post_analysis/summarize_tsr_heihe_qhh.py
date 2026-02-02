#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "validation" / "tsr" / "py"))
from shud_reader import DatReader  # noqa: E402
from tsr_core import read_mesh_normals  # noqa: E402


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


def _compute_aspect_slope(normals: Dict[int, Tuple[float, float, float]]) -> tuple[np.ndarray, np.ndarray]:
    n = len(normals)
    nx = np.zeros(n)
    ny = np.zeros(n)
    nz = np.zeros(n)
    for eid, (a, b, c) in normals.items():
        nx[eid - 1] = float(a)
        ny[eid - 1] = float(b)
        nz[eid - 1] = float(c)

    nz_clamped = np.clip(nz, 0.0, 1.0)
    slope = np.arctan2(np.hypot(nx, ny), nz_clamped)

    # Match src/classes/Element.cpp (North=0 at +y, East=pi/2 at +x).
    aspect = np.arctan2(nx, ny)
    aspect = np.where(aspect < 0.0, aspect + 2.0 * math.pi, aspect)
    aspect = np.where(aspect >= 2.0 * math.pi, aspect - 2.0 * math.pi, aspect)
    aspect = np.where(slope < 1e-6, 0.0, aspect)
    return aspect, slope


def _circular_diff(a: np.ndarray, center: float) -> np.ndarray:
    return np.abs(((a - center + math.pi) % (2.0 * math.pi)) - math.pi)


@dataclass(frozen=True)
class AspectGroups:
    mask_south: np.ndarray
    mask_north: np.ndarray
    mask_flat: np.ndarray
    mask_other: np.ndarray
    areas: np.ndarray


def build_aspect_groups(*, mesh_path: Path, slope_deg_min: float = 5.0, aspect_halfwidth_deg: float = 45.0) -> AspectGroups:
    areas = _read_mesh_areas(mesh_path)
    aspect, slope = _compute_aspect_slope(read_mesh_normals(mesh_path))

    slope_thr = math.radians(float(slope_deg_min))
    half = math.radians(float(aspect_halfwidth_deg))

    valid = slope >= slope_thr
    flat = ~valid
    south = valid & (_circular_diff(aspect, math.pi) <= half)
    north = valid & (_circular_diff(aspect, 0.0) <= half)
    other = valid & ~(south | north)

    return AspectGroups(mask_south=south, mask_north=north, mask_flat=flat, mask_other=other, areas=areas)


def _aw_mean_series(values: np.ndarray, *, areas: np.ndarray, mask: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    w = np.asarray(areas, dtype=float) * mask.astype(float)

    finite = np.isfinite(values)
    if values.ndim != 2 or w.ndim != 1 or values.shape[1] != w.size:
        raise SummarizeError(f"shape mismatch: values={values.shape}, weights={w.shape}")

    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        num = (np.where(finite, values, 0.0) @ w).astype(float)
        den = (finite.astype(float) @ w).astype(float)
        return np.where(den > 0.0, num / den, np.nan)


def _flatten_stats(x: np.ndarray) -> dict[str, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {"min": float("nan"), "p50": float("nan"), "p95": float("nan"), "max": float("nan")}
    return {
        "min": float(np.min(x)),
        "p50": float(np.percentile(x, 50.0)),
        "p95": float(np.percentile(x, 95.0)),
        "max": float(np.max(x)),
    }


def basin_totals_mm(*, basinwb_path: Path, mesh_area_m2: float) -> dict[str, float]:
    _, vals, _ = _read_matrix_fast(basinwb_path)
    total_m3 = vals.sum(axis=0)
    mm = total_m3 / mesh_area_m2 * 1000.0
    keys = ["dS", "P", "ET", "Qout", "Qedge", "Qbc", "Qss", "NonConsEdge", "Resid"]
    return {k: float(v) for k, v in zip(keys, mm)}


def basin_resid_stats_mm(*, basinwb_path: Path, mesh_area_m2: float) -> dict[str, float]:
    t_min, vals, _ = _read_matrix_fast(basinwb_path)
    resid_mm = vals[:, 8] / mesh_area_m2 * 1000.0
    abs_mm = np.abs(resid_mm)
    finite = np.isfinite(resid_mm)
    if not np.any(finite):
        raise SummarizeError(f"no finite residuals in {basinwb_path}")
    return {
        "max_abs": float(np.max(abs_mm[finite])),
        "rms": float(np.sqrt(np.mean(resid_mm[finite] ** 2))),
        "p95_abs": float(np.percentile(abs_mm[finite], 95.0)),
        "p99_abs": float(np.percentile(abs_mm[finite], 99.0)),
        "sum": float(np.sum(resid_mm[finite])),
    }


def _read_lake_init_stage(cfg_ic: Path) -> float:
    lines = cfg_ic.read_text(encoding="utf-8", errors="ignore").splitlines()
    for i, ln in enumerate(lines):
        if "LakeStage" in ln:
            for j in range(i + 1, len(lines)):
                parts = lines[j].split()
                if len(parts) >= 2:
                    return float(parts[1])
    raise SummarizeError(f"failed to parse initial LakeStage from {cfg_ic}")


def lake_corrected_cumsum_resid_mm(*, out_dir: Path, basin: str, mesh_area_m2: float, cfg_ic: Path) -> float:
    _, wb, _ = _read_matrix_fast(out_dir / f"{basin}.basinwbfull.dat")
    resid_m3 = wb[:, 8]

    _, stage, _ = _read_matrix_fast(out_dir / f"{basin}.lakystage.dat")
    _, area, _ = _read_matrix_fast(out_dir / f"{basin}.lakatop.dat")
    _, evap, _ = _read_matrix_fast(out_dir / f"{basin}.lakvevap.dat")

    stage = stage[:, 0]
    area = area[:, 0]
    evap = evap[:, 0]

    init_stage = _read_lake_init_stage(cfg_ic)
    init_area = float(area[0])

    ds_lake = np.empty_like(resid_m3)
    prev_stage = float(init_stage)
    prev_area = float(init_area)
    for k in range(ds_lake.size):
        h = float(stage[k])
        a = float(area[k])
        ds_lake[k] = 0.5 * (prev_area + a) * (h - prev_stage)
        prev_stage = h
        prev_area = a

    evap_vol = evap * area
    resid_full_m3 = resid_m3 + ds_lake + evap_vol
    resid_full_mm = resid_full_m3 / mesh_area_m2 * 1000.0
    return float(np.sum(resid_full_mm))


def main() -> int:
    ap = argparse.ArgumentParser(description="Summarize docs/tsr_heihe_qhh_analysis.md numbers from outputs.")
    ap.add_argument("--slope-deg-min", type=float, default=5.0)
    ap.add_argument("--aspect-halfwidth-deg", type=float, default=45.0)

    ap.add_argument("--heihe-mesh", type=Path, default=Path("input/heihe/heihe.sp.mesh"))
    ap.add_argument("--heihe-base", type=Path, default=Path("output/heihe.base"))
    ap.add_argument("--heihe-tsr", type=Path, default=Path("output/heihe.tsr"))

    ap.add_argument("--qhh-mesh", type=Path, default=Path("input/qhh/qhh.sp.mesh"))
    ap.add_argument("--qhh-base", type=Path, default=Path("output/qhh.base"))
    ap.add_argument("--qhh-tsr", type=Path, default=Path("output/qhh.tsr"))
    ap.add_argument("--qhh-ic", type=Path, default=Path("input/qhh/qhh.cfg.ic"))
    args = ap.parse_args()

    def summarize_basin(basin: str, mesh: Path, out_base: Path, out_tsr: Path) -> None:
        groups = build_aspect_groups(mesh_path=mesh, slope_deg_min=args.slope_deg_min, aspect_halfwidth_deg=args.aspect_halfwidth_deg)
        mask_other_all = ~(groups.mask_south | groups.mask_north)
        A = float(groups.areas.sum())
        if not (A > 0.0):
            raise SummarizeError(f"invalid mesh area for {basin}: {A}")

        area_s = float(np.sum(groups.areas[groups.mask_south]))
        area_n = float(np.sum(groups.areas[groups.mask_north]))
        area_f = float(np.sum(groups.areas[groups.mask_flat]))
        area_o = float(np.sum(groups.areas[groups.mask_other]))
        n_ele = int(groups.areas.size)

        print(f"\n== {basin} ==")
        print(f"A_total = {A:.3f} m^2, N_ele = {n_ele}")
        print(
            "area fractions: "
            f"South={area_s/A*100.0:.2f}%, North={area_n/A*100.0:.2f}%, Flat={area_f/A*100.0:.2f}%, Other={area_o/A*100.0:.2f}%"
        )

        # rn_factor distribution (TSR=ON; across time x elements)
        _, rn_factor, _ = _read_matrix_fast(out_tsr / f"{basin}.rn_factor.dat")
        st = _flatten_stats(rn_factor.reshape(-1))
        print(f"rn_factor (TSR=ON): min={st['min']:.3f}, p50={st['p50']:.3f}, p95={st['p95']:.3f}, max={st['max']:.3f}")

        # Group-wise TSR effect on rn_t / ETa (area-weighted time mean)
        _, rn_base, _ = _read_matrix_fast(out_base / f"{basin}.rn_t.dat")
        _, rn_tsr, _ = _read_matrix_fast(out_tsr / f"{basin}.rn_t.dat")

        rn_b_s = _aw_mean_series(rn_base, areas=groups.areas, mask=groups.mask_south)
        rn_b_n = _aw_mean_series(rn_base, areas=groups.areas, mask=groups.mask_north)
        rn_b_o = _aw_mean_series(rn_base, areas=groups.areas, mask=mask_other_all)
        rn_t_s = _aw_mean_series(rn_tsr, areas=groups.areas, mask=groups.mask_south)
        rn_t_n = _aw_mean_series(rn_tsr, areas=groups.areas, mask=groups.mask_north)
        rn_t_o = _aw_mean_series(rn_tsr, areas=groups.areas, mask=mask_other_all)

        print("TSR effect on rn_t (TSR−BASE), area-weighted time mean:")
        print(f"  South: {float(np.nanmean(rn_t_s - rn_b_s)):+.4g} W/m^2")
        print(f"  North: {float(np.nanmean(rn_t_n - rn_b_n)):+.4g} W/m^2")
        print(f"  Other: {float(np.nanmean(rn_t_o - rn_b_o)):+.4g} W/m^2")

        _, et_base, _ = _read_matrix_fast(out_base / f"{basin}.eleveta.dat")
        _, et_tsr, _ = _read_matrix_fast(out_tsr / f"{basin}.eleveta.dat")
        et_b_s = _aw_mean_series(et_base, areas=groups.areas, mask=groups.mask_south) * 1000.0
        et_b_n = _aw_mean_series(et_base, areas=groups.areas, mask=groups.mask_north) * 1000.0
        et_b_o = _aw_mean_series(et_base, areas=groups.areas, mask=mask_other_all) * 1000.0
        et_t_s = _aw_mean_series(et_tsr, areas=groups.areas, mask=groups.mask_south) * 1000.0
        et_t_n = _aw_mean_series(et_tsr, areas=groups.areas, mask=groups.mask_north) * 1000.0
        et_t_o = _aw_mean_series(et_tsr, areas=groups.areas, mask=mask_other_all) * 1000.0
        print("TSR effect on ETa (TSR−BASE), area-weighted time mean:")
        print(f"  South: {float(np.nanmean(et_t_s - et_b_s)):+.4g} mm/day")
        print(f"  North: {float(np.nanmean(et_t_n - et_b_n)):+.4g} mm/day")
        print(f"  Other: {float(np.nanmean(et_t_o - et_b_o)):+.4g} mm/day")

        # Basin totals
        tb = basin_totals_mm(basinwb_path=out_base / f"{basin}.basinwbfull.dat", mesh_area_m2=A)
        tt = basin_totals_mm(basinwb_path=out_tsr / f"{basin}.basinwbfull.dat", mesh_area_m2=A)
        print("basin totals (mm over mesh area):")  # key subset
        print(f"  BASE: P={tb['P']:.3f}, ET={tb['ET']:.3f}, Qout={tb['Qout']:.3f}, dS={tb['dS']:.3f}, Resid={tb['Resid']:.3f}")
        print(f"  TSR : P={tt['P']:.3f}, ET={tt['ET']:.3f}, Qout={tt['Qout']:.3f}, dS={tt['dS']:.3f}, Resid={tt['Resid']:.3f}")
        print(f"  TSR-BASE: ΔET={tt['ET']-tb['ET']:+.3f}, ΔQout={tt['Qout']-tb['Qout']:+.3f}, ΔdS={tt['dS']-tb['dS']:+.3f}")

        # Storage aspect deltas (yUnsat / yGW)
        def aspect_delta_time_mean(out_dir: Path, var: str) -> float:
            _, vv, _ = _read_matrix_fast(out_dir / f"{basin}.{var}.dat")
            ss = _aw_mean_series(vv, areas=groups.areas, mask=groups.mask_south)
            nn = _aw_mean_series(vv, areas=groups.areas, mask=groups.mask_north)
            return float(np.nanmean(ss - nn))

        du_b = aspect_delta_time_mean(out_base, "eleyunsat")
        du_t = aspect_delta_time_mean(out_tsr, "eleyunsat")
        dg_b = aspect_delta_time_mean(out_base, "eleygw")
        dg_t = aspect_delta_time_mean(out_tsr, "eleygw")
        print("Δ storage (South-North), area-weighted time mean:")
        print(f"  yUnsat: BASE={du_b:+.5g} m, TSR={du_t:+.5g} m, effect={du_t-du_b:+.5g} m")
        print(f"  yGW   : BASE={dg_b:+.5g} m, TSR={dg_t:+.5g} m, effect={dg_t-dg_b:+.5g} m")

        # Residual stats
        rb = basin_resid_stats_mm(basinwb_path=out_base / f"{basin}.basinwbfull.dat", mesh_area_m2=A)
        rt = basin_resid_stats_mm(basinwb_path=out_tsr / f"{basin}.basinwbfull.dat", mesh_area_m2=A)
        print("basin residual stats (mm / interval):")
        print(f"  BASE: max|resid|={rb['max_abs']:.3f}, rms={rb['rms']:.4g}, p95={rb['p95_abs']:.4g}, p99={rb['p99_abs']:.4g}, Σresid={rb['sum']:+.3f}")
        print(f"  TSR : max|resid|={rt['max_abs']:.3f}, rms={rt['rms']:.4g}, p95={rt['p95_abs']:.4g}, p99={rt['p99_abs']:.4g}, Σresid={rt['sum']:+.3f}")

        if basin == "qhh":
            csum_b = lake_corrected_cumsum_resid_mm(out_dir=out_base, basin=basin, mesh_area_m2=A, cfg_ic=args.qhh_ic)
            csum_t = lake_corrected_cumsum_resid_mm(out_dir=out_tsr, basin=basin, mesh_area_m2=A, cfg_ic=args.qhh_ic)
            print(f"qhh corrected Σresid (raw + ΔS_lake + E_lake): BASE={csum_b:+.3f} mm, TSR={csum_t:+.3f} mm")

    summarize_basin("heihe", args.heihe_mesh, args.heihe_base, args.heihe_tsr)
    summarize_basin("qhh", args.qhh_mesh, args.qhh_base, args.qhh_tsr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
