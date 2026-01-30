#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

# Reuse the repo's SHUD .dat reader + mesh normal computation (kept C++-consistent).
_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "validation" / "tsr" / "py"))
from shud_reader import DatReader  # noqa: E402
from tsr_core import read_mesh_normals  # noqa: E402


class PlotError(RuntimeError):
    pass


def _setup_fonts() -> None:
    preferred = [
        "PingFang SC",
        "Arial Unicode MS",
        "Noto Sans CJK SC",
        "SimHei",
        "Microsoft YaHei",
    ]
    available = {f.name for f in fm.fontManager.ttflist}
    for name in preferred:
        if name in available:
            plt.rcParams["font.sans-serif"] = [name]
            break
    plt.rcParams["axes.unicode_minus"] = False


def _read_matrix_fast(path: Path) -> tuple[np.ndarray, np.ndarray, DatReader]:
    r = DatReader(path)
    m = r.meta
    dtype = np.dtype(m.endianness + "f8")
    total = m.num_records * (1 + m.num_var)
    with m.path.open("rb") as f:
        f.seek(m.data_offset, 0)
        data = np.fromfile(f, dtype=dtype, count=total)
    if data.size != total:
        raise PlotError(f"truncated read: {path} (got {data.size}, expected {total})")
    data = data.reshape((m.num_records, 1 + m.num_var))
    t_min = data[:, 0].astype(float)
    values = data[:, 1:].astype(float)
    return t_min, values, r


def _read_mesh_areas(mesh_path: Path) -> np.ndarray:
    lines = mesh_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise PlotError(f"empty mesh file: {mesh_path}")
    num_ele = int(lines[0].split()[0])

    elements: list[tuple[int, int, int, int]] = []
    for k in range(num_ele):
        cols = lines[2 + k].split()
        if len(cols) < 4:
            raise PlotError(f"invalid element row in {mesh_path}: {lines[2+k]!r}")
        elements.append((int(cols[0]), int(cols[1]), int(cols[2]), int(cols[3])))

    node_hdr_idx = 2 + num_ele
    num_node = int(lines[node_hdr_idx].split()[0])
    node_start = node_hdr_idx + 2

    x = np.zeros(num_node + 1)
    y = np.zeros(num_node + 1)
    for k in range(num_node):
        cols = lines[node_start + k].split()
        if len(cols) < 3:
            raise PlotError(f"invalid node row in {mesh_path}: {lines[node_start+k]!r}")
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


def _read_mesh_geom(mesh_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Return (areas, aq_depth_ele).
    - areas: element triangle area (m^2)
    - aq_depth_ele: element-mean aquifer depth (m), averaged from node AqDepth
    """
    lines = mesh_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise PlotError(f"empty mesh file: {mesh_path}")
    num_ele = int(lines[0].split()[0])

    elements: list[tuple[int, int, int, int]] = []
    for k in range(num_ele):
        cols = lines[2 + k].split()
        if len(cols) < 4:
            raise PlotError(f"invalid element row in {mesh_path}: {lines[2+k]!r}")
        elements.append((int(cols[0]), int(cols[1]), int(cols[2]), int(cols[3])))

    node_hdr_idx = 2 + num_ele
    num_node = int(lines[node_hdr_idx].split()[0])
    node_start = node_hdr_idx + 2

    x = np.zeros(num_node + 1)
    y = np.zeros(num_node + 1)
    aqd = np.zeros(num_node + 1)
    for k in range(num_node):
        cols = lines[node_start + k].split()
        if len(cols) < 4:
            raise PlotError(f"invalid node row in {mesh_path}: {lines[node_start+k]!r}")
        nid = int(cols[0])
        x[nid] = float(cols[1])
        y[nid] = float(cols[2])
        aqd[nid] = float(cols[3])

    areas = np.zeros(num_ele)
    aqd_ele = np.zeros(num_ele)
    for eid, n1, n2, n3 in elements:
        x1, y1 = x[n1], y[n1]
        x2, y2 = x[n2], y[n2]
        x3, y3 = x[n3], y[n3]
        areas[eid - 1] = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
        aqd_ele[eid - 1] = (aqd[n1] + aqd[n2] + aqd[n3]) / 3.0

    return areas, aqd_ele


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


def build_aspect_groups(
    *, mesh_path: Path, slope_deg_min: float = 5.0, aspect_halfwidth_deg: float = 45.0
) -> AspectGroups:
    areas = _read_mesh_areas(mesh_path)
    aspect, slope = _compute_aspect_slope(read_mesh_normals(mesh_path))

    slope_thr = math.radians(float(slope_deg_min))
    half = math.radians(float(aspect_halfwidth_deg))

    valid = slope >= slope_thr
    flat = ~valid
    south = valid & (_circular_diff(aspect, math.pi) <= half)
    north = valid & (_circular_diff(aspect, 0.0) <= half)
    other = valid & ~(south | north)

    return AspectGroups(
        mask_south=south,
        mask_north=north,
        mask_flat=flat,
        mask_other=other,
        areas=areas,
    )


def _aw_mean_series(values: np.ndarray, *, areas: np.ndarray, mask: np.ndarray) -> np.ndarray:
    w = areas * mask.astype(float)
    s = float(w.sum())
    if not (s > 0.0):
        raise PlotError("empty group (area sum=0)")
    return values.dot(w) / s


def plot_basin_residual_timeseries(
    *,
    out_base: Path,
    out_tsr: Path,
    mesh_path: Path,
    out_path: Path,
    title: str,
) -> None:
    areas = _read_mesh_areas(mesh_path)
    A = float(areas.sum())

    def read_resid_mm(out_dir: Path) -> tuple[np.ndarray, np.ndarray]:
        t_min, vals, _ = _read_matrix_fast(out_dir / "ccw.basinwbfull.dat")
        resid_m3 = vals[:, 8]
        resid_mm = resid_m3 / A * 1000.0
        day = t_min / 1440.0
        return day, resid_mm

    day_b, resid_b = read_resid_mm(out_base)
    day_t, resid_t = read_resid_mm(out_tsr)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.5, 7.2), sharex=True)

    ax1.plot(day_b, resid_b, lw=0.8, label="Baseline (TSR=OFF)")
    ax1.plot(day_t, resid_t, lw=0.8, label="TSR=ON")
    ax1.axhline(0.0, color="k", lw=0.7, alpha=0.5)
    ax1.set_ylabel("Residual (mm / interval)")
    ax1.set_title(title)
    ax1.legend(loc="upper right", frameon=True)

    ax2.plot(day_b, np.cumsum(resid_b), lw=0.9, label="Baseline cumulative")
    ax2.plot(day_t, np.cumsum(resid_t), lw=0.9, label="TSR cumulative")
    ax2.axhline(0.0, color="k", lw=0.7, alpha=0.5)
    ax2.set_xlabel("Day index (left endpoint)")
    ax2.set_ylabel("Cumulative residual (mm)")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_basin_abs_cdf(
    *,
    out_base: Path,
    out_tsr: Path,
    mesh_path: Path,
    out_path: Path,
    title: str,
) -> None:
    areas = _read_mesh_areas(mesh_path)
    A = float(areas.sum())

    def read_abs_mm(out_dir: Path) -> np.ndarray:
        _, vals, _ = _read_matrix_fast(out_dir / "ccw.basinwbfull.dat")
        resid_m3 = vals[:, 8]
        resid_mm = resid_m3 / A * 1000.0
        return np.abs(resid_mm)

    a_b = read_abs_mm(out_base)
    a_t = read_abs_mm(out_tsr)

    def cdf(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        xs = np.sort(x)
        ys = (np.arange(xs.size) + 1) / xs.size
        return xs, ys

    xb, yb = cdf(a_b)
    xt, yt = cdf(a_t)

    fig, ax = plt.subplots(1, 1, figsize=(10.5, 4.6))
    ax.plot(xb, yb, lw=1.2, label="Baseline (TSR=OFF)")
    ax.plot(xt, yt, lw=1.2, label="TSR=ON")
    ax.set_xscale("log")
    ax.set_xlabel("|Residual| (mm / interval)")
    ax.set_ylabel("CDF")
    ax.set_title(title)
    ax.grid(True, which="both", ls="--", lw=0.5, alpha=0.5)
    ax.legend(loc="lower right", frameon=True)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_element_residual_cdf(
    *,
    out_base: Path,
    out_tsr: Path,
    out_path: Path,
    title: str,
) -> None:
    def read_abs_mm(out_dir: Path, filename: str) -> np.ndarray:
        _, vals, _ = _read_matrix_fast(out_dir / filename)
        return np.abs(vals).ravel() * 1000.0  # m -> mm

    base_a = read_abs_mm(out_base, "ccw.elewb3_resid.dat")
    base_b = read_abs_mm(out_base, "ccw.elewb3_budget_resid.dat")
    tsr_a = read_abs_mm(out_tsr, "ccw.elewb3_resid.dat")
    tsr_b = read_abs_mm(out_tsr, "ccw.elewb3_budget_resid.dat")

    def cdf(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        xs = np.sort(x)
        ys = (np.arange(xs.size) + 1) / xs.size
        return xs, ys

    fig, ax = plt.subplots(1, 1, figsize=(10.5, 4.6))
    for x, label, ls in [
        (base_a, "Baseline A: ΔS−∫dS/dt", "-"),
        (base_b, "Baseline B: ΔS−∫RHS_budget", "--"),
        (tsr_a, "TSR A: ΔS−∫dS/dt", "-"),
        (tsr_b, "TSR B: ΔS−∫RHS_budget", "--"),
    ]:
        xs, ys = cdf(x)
        ax.plot(xs, ys, lw=1.0, ls=ls, label=label)

    ax.set_xscale("log")
    ax.set_xlabel("|Residual| (mm / interval)")
    ax.set_ylabel("CDF (over time×element samples)")
    ax.set_title(title)
    ax.grid(True, which="both", ls="--", lw=0.5, alpha=0.5)
    ax.legend(loc="lower right", fontsize=8, frameon=True)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_aspect_group_fractions(
    *,
    mesh_path: Path,
    slope_deg_min: float,
    aspect_halfwidth_deg: float,
    out_path: Path,
    title: str,
) -> None:
    g = build_aspect_groups(mesh_path=mesh_path, slope_deg_min=slope_deg_min, aspect_halfwidth_deg=aspect_halfwidth_deg)
    A = float(g.areas.sum())

    labels = ["South", "North", "Flat(<thr)", "Other(aspect)"]
    masks = [g.mask_south, g.mask_north, g.mask_flat, g.mask_other]
    fracs = [float(g.areas[m].sum() / A) for m in masks]

    fig, ax = plt.subplots(1, 1, figsize=(10.5, 4.6))
    bars = ax.bar(labels, fracs)
    ax.set_ylim(0.0, max(fracs) * 1.25)
    ax.set_ylabel("Area fraction (-)")
    ax.set_title(title)
    for b, f in zip(bars, fracs):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(), f"{f*100:.1f}%", ha="center", va="bottom", fontsize=10)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_aspect_delta_timeseries(
    *,
    out_base: Path,
    out_tsr: Path,
    mesh_path: Path,
    slope_deg_min: float,
    aspect_halfwidth_deg: float,
    out_path: Path,
    title: str,
) -> None:
    g = build_aspect_groups(mesh_path=mesh_path, slope_deg_min=slope_deg_min, aspect_halfwidth_deg=aspect_halfwidth_deg)

    def delta_series(out_dir: Path, dat_name: str) -> tuple[np.ndarray, np.ndarray]:
        t_min, vals, _ = _read_matrix_fast(out_dir / dat_name)
        south = _aw_mean_series(vals, areas=g.areas, mask=g.mask_south)
        north = _aw_mean_series(vals, areas=g.areas, mask=g.mask_north)
        day = t_min / 1440.0
        return day, (south - north)

    day_b, rn_b = delta_series(out_base, "ccw.rn_t.dat")
    day_t, rn_t = delta_series(out_tsr, "ccw.rn_t.dat")
    _, et_b_mday = delta_series(out_base, "ccw.eleveta.dat")
    _, et_t_mday = delta_series(out_tsr, "ccw.eleveta.dat")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.5, 7.2), sharex=True)

    ax1.plot(day_b, rn_b, lw=0.9, label="Baseline (TSR=OFF)")
    ax1.plot(day_t, rn_t, lw=0.9, label="TSR=ON")
    ax1.axhline(0.0, color="k", lw=0.7, alpha=0.5)
    ax1.set_ylabel("Δ rn_t (W/m²)  South−North")
    ax1.set_title(title)
    ax1.legend(loc="upper right", frameon=True)

    ax2.plot(day_b, et_b_mday * 1000.0, lw=0.9, label="Baseline (TSR=OFF)")
    ax2.plot(day_t, et_t_mday * 1000.0, lw=0.9, label="TSR=ON")
    ax2.axhline(0.0, color="k", lw=0.7, alpha=0.5)
    ax2.set_xlabel("Day index (left endpoint)")
    ax2.set_ylabel("Δ ETa (mm/day)  South−North")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_aspect_storage_delta_timeseries(
    *,
    out_base: Path,
    out_tsr: Path,
    mesh_path: Path,
    slope_deg_min: float,
    aspect_halfwidth_deg: float,
    out_path: Path,
    title: str,
) -> None:
    areas, aqd_ele = _read_mesh_geom(mesh_path)
    g = build_aspect_groups(
        mesh_path=mesh_path, slope_deg_min=slope_deg_min, aspect_halfwidth_deg=aspect_halfwidth_deg
    )

    def delta_series(out_dir: Path, dat_name: str) -> tuple[np.ndarray, np.ndarray]:
        t_min, vals, _ = _read_matrix_fast(out_dir / dat_name)
        south = _aw_mean_series(vals, areas=g.areas, mask=g.mask_south)
        north = _aw_mean_series(vals, areas=g.areas, mask=g.mask_north)
        day = t_min / 1440.0
        return day, (south - north)

    def delta_unsat_ratio(out_dir: Path) -> tuple[np.ndarray, np.ndarray]:
        t_min, uns, _ = _read_matrix_fast(out_dir / "ccw.eleyunsat.dat")
        _, gw, _ = _read_matrix_fast(out_dir / "ccw.eleygw.dat")
        deficit = aqd_ele.reshape((1, -1)) - gw
        ratio = np.where(deficit > 1e-12, uns / deficit, np.nan)
        # area-weighted mean ignoring NaNs
        def wmean(mask: np.ndarray) -> np.ndarray:
            m = mask.reshape((1, -1)) & np.isfinite(ratio)
            w = areas.reshape((1, -1)) * m.astype(float)
            ws = w.sum(axis=1)
            return (ratio * w).sum(axis=1) / ws

        south = wmean(g.mask_south)
        north = wmean(g.mask_north)
        day = t_min / 1440.0
        return day, (south - north)

    day_b, du_b = delta_series(out_base, "ccw.eleyunsat.dat")
    day_t, du_t = delta_series(out_tsr, "ccw.eleyunsat.dat")
    _, dg_b = delta_series(out_base, "ccw.eleygw.dat")
    _, dg_t = delta_series(out_tsr, "ccw.eleygw.dat")
    _, dr_b = delta_unsat_ratio(out_base)
    _, dr_t = delta_unsat_ratio(out_tsr)

    # TSR effect on contrast (difference-of-differences)
    du_eff = du_t - du_b
    dr_eff = dr_t - dr_b
    dg_eff = dg_t - dg_b

    fig, axes = plt.subplots(3, 1, figsize=(10.5, 8.6), sharex=True)

    ax = axes[0]
    ax.plot(day_t, du_t, lw=0.9, label="TSR=ON")
    ax.plot(day_b, du_b, lw=1.0, ls="--", label="Baseline (TSR=OFF)")
    ax.plot(day_t, du_eff, lw=0.8, color="k", alpha=0.6, label="TSR effect (ΔTSR−ΔBASE)")
    ax.axhline(0.0, color="k", lw=0.7, alpha=0.4)
    ax.set_ylabel("Δ yUnsat (m)\nSouth−North")
    ax.set_title(title)
    ax.legend(loc="upper right", frameon=True, fontsize=9)

    ax = axes[1]
    ax.plot(day_t, dr_t, lw=0.9, label="TSR=ON")
    ax.plot(day_b, dr_b, lw=1.0, ls="--", label="Baseline (TSR=OFF)")
    ax.plot(day_t, dr_eff, lw=0.8, color="k", alpha=0.6, label="TSR effect (ΔTSR−ΔBASE)")
    ax.axhline(0.0, color="k", lw=0.7, alpha=0.4)
    ax.set_ylabel("Δ UnsatRatio (-)\nSouth−North")

    ax = axes[2]
    # ΔyGW is O(0.1m) but TSR effect (difference-of-differences) is O(1e-4m).
    # Plot ΔyGW on the left axis and TSR effect on a right axis (mm) so neither becomes visually invisible.
    ax.plot(day_t, dg_t, lw=0.9, label="TSR=ON")
    ax.plot(day_b, dg_b, lw=1.0, ls="--", label="Baseline (TSR=OFF)")
    ax.set_xlabel("Day index (left endpoint)")
    ax.set_ylabel("Δ yGW (m)\nSouth−North")
    axr = ax.twinx()
    axr.plot(day_t, dg_eff * 1000.0, lw=0.9, color="k", alpha=0.6, label="TSR effect (ΔTSR−ΔBASE)")
    axr.axhline(0.0, color="k", lw=0.7, alpha=0.4)
    axr.set_ylabel("TSR effect on ΔyGW (mm)\n(ΔTSR−ΔBASE)")
    # Reasonable visibility window; auto-scale but keep a minimum range.
    eff_mm = dg_eff * 1000.0
    eff_max = float(np.nanmax(np.abs(eff_mm))) if eff_mm.size else 0.0
    y = max(1.0, eff_max * 1.25)
    axr.set_ylim(-y, y)
    # Combine legends from left and right axes (so baseline is always visible in legend).
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = axr.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper right", frameon=True, fontsize=9)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate plots for docs/water_balance_verification.md")
    ap.add_argument("--mesh", type=Path, default=Path("input/ccw/ccw.sp.mesh"))
    ap.add_argument("--base", type=Path, default=Path("output/ccw.base"))
    ap.add_argument("--tsr", type=Path, default=Path("output/ccw.tsr"))
    ap.add_argument("--base-2d", type=Path, default=Path("output/ccw.base.2d"))
    ap.add_argument("--tsr-2d", type=Path, default=Path("output/ccw.tsr.2d"))
    ap.add_argument("--outdir", type=Path, default=Path("docs/figures/water_balance"))
    ap.add_argument("--slope-deg-min", type=float, default=5.0)
    ap.add_argument("--aspect-halfwidth-deg", type=float, default=45.0)
    args = ap.parse_args()

    _setup_fonts()
    plt.style.use("seaborn-v0_8-whitegrid")

    # Full run basin plots
    plot_basin_residual_timeseries(
        out_base=args.base,
        out_tsr=args.tsr,
        mesh_path=args.mesh,
        out_path=args.outdir / "basin_residual_full.png",
        title="Basin water-balance residual (full run)",
    )
    plot_basin_abs_cdf(
        out_base=args.base,
        out_tsr=args.tsr,
        mesh_path=args.mesh,
        out_path=args.outdir / "basin_residual_abs_cdf_full.png",
        title="Basin |residual| CDF (full run)",
    )

    # Short run basin plot (2 days)
    if (args.base_2d / "ccw.basinwbfull.dat").exists() and (args.tsr_2d / "ccw.basinwbfull.dat").exists():
        plot_basin_residual_timeseries(
            out_base=args.base_2d,
            out_tsr=args.tsr_2d,
            mesh_path=args.mesh,
            out_path=args.outdir / "basin_residual_2d.png",
            title="Basin water-balance residual (2-day short run)",
        )

    # Per-element residual distribution (full run)
    plot_element_residual_cdf(
        out_base=args.base,
        out_tsr=args.tsr,
        out_path=args.outdir / "element_residual_abs_cdf_full.png",
        title="Per-element |residual| CDF (full run)",
    )

    # Aspect groups and delta time series
    plot_aspect_group_fractions(
        mesh_path=args.mesh,
        slope_deg_min=args.slope_deg_min,
        aspect_halfwidth_deg=args.aspect_halfwidth_deg,
        out_path=args.outdir / "aspect_group_area_fractions.png",
        title=f"Aspect groups (slope≥{args.slope_deg_min:g}°, halfwidth=±{args.aspect_halfwidth_deg:g}°)",
    )
    plot_aspect_delta_timeseries(
        out_base=args.base,
        out_tsr=args.tsr,
        mesh_path=args.mesh,
        slope_deg_min=args.slope_deg_min,
        aspect_halfwidth_deg=args.aspect_halfwidth_deg,
        out_path=args.outdir / "aspect_delta_timeseries_full.png",
        title="Aspect response (South−North) time series (full run)",
    )
    plot_aspect_storage_delta_timeseries(
        out_base=args.base,
        out_tsr=args.tsr,
        mesh_path=args.mesh,
        slope_deg_min=args.slope_deg_min,
        aspect_halfwidth_deg=args.aspect_halfwidth_deg,
        out_path=args.outdir / "aspect_storage_delta_timeseries_full.png",
        title="Aspect storage response (South−North) time series (full run)",
    )

    print(f"OK: wrote figures under {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
