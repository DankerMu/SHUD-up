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
import matplotlib.tri as mtri
from matplotlib.colors import TwoSlopeNorm
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


# Color-blind friendly palette (Okabe–Ito)
_C_BASE = "#0072B2"  # blue
_C_TSR = "#D55E00"  # vermillion
_C_GREEN = "#009E73"
_C_PURPLE = "#CC79A7"
_C_GREY = "#666666"


def _pub_rc() -> dict:
    # A lightweight "paper figure" style.
    return {
        "font.size": 12,
        "axes.titlesize": 12,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.8,
        "lines.linewidth": 1.2,
        "savefig.transparent": False,
    }


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


def _read_mesh_triangulation(mesh_path: Path) -> tuple[mtri.Triangulation, np.ndarray]:
    """
    Return (triangulation, ele_ids_in_order).

    - triangulation uses node XY coordinates converted to km.
    - ele_ids_in_order matches the triangle order (1-based element ids).
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

    node_ids: list[int] = []
    xs: list[float] = []
    ys: list[float] = []
    for k in range(num_node):
        cols = lines[node_start + k].split()
        if len(cols) < 3:
            raise PlotError(f"invalid node row in {mesh_path}: {lines[node_start+k]!r}")
        node_ids.append(int(cols[0]))
        xs.append(float(cols[1]) / 1000.0)
        ys.append(float(cols[2]) / 1000.0)

    id_to_idx = {nid: i for i, nid in enumerate(node_ids)}
    triangles = np.zeros((num_ele, 3), dtype=int)
    ele_ids = np.zeros(num_ele, dtype=int)
    for i, (eid, n1, n2, n3) in enumerate(elements):
        try:
            triangles[i, :] = (id_to_idx[n1], id_to_idx[n2], id_to_idx[n3])
        except KeyError as e:
            raise PlotError(f"mesh element references missing node: {e}") from e
        ele_ids[i] = eid

    tri = mtri.Triangulation(np.asarray(xs), np.asarray(ys), triangles=triangles)
    return tri, ele_ids


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


def _nanmean_time(values: np.ndarray) -> np.ndarray:
    v = values.astype(float, copy=False)
    v = np.where(v <= -9000.0, np.nan, v)
    return np.nanmean(v, axis=0)


def _nanintegrate_time(values: np.ndarray, t_min: np.ndarray) -> np.ndarray:
    """
    Integrate time-averaged series to cumulative energy/amount per column.

    - `values`: shape (nt, ncol), interpreted as interval mean.
    - `t_min`: left-endpoint minutes, shape (nt,)
    Returns: cumulative sum over time of `values * dt`, where dt is derived from `t_min` (seconds).
    """
    if values.ndim != 2 or t_min.ndim != 1 or values.shape[0] != t_min.size:
        raise PlotError(f"invalid shapes for integration: values={values.shape}, t_min={t_min.shape}")

    v = values.astype(float, copy=False)
    v = np.where(v <= -9000.0, np.nan, v)

    if t_min.size < 2:
        return np.nansum(v, axis=0) * 0.0

    dt_min = np.diff(t_min.astype(float))
    dt_last = float(np.nanmedian(dt_min))
    if not (dt_last > 0.0):
        dt_last = float(dt_min[-1]) if (dt_min.size and dt_min[-1] > 0.0) else 0.0
    dt_s = np.concatenate([dt_min, np.asarray([dt_last])]) * 60.0

    return np.nansum(v * dt_s.reshape((-1, 1)), axis=0)


def _robust_sym_lim(values: np.ndarray, *, q: float = 99.0) -> float:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 1.0
    lim = float(np.nanpercentile(np.abs(v), q))
    if not (lim > 0.0):
        lim = float(np.nanmax(np.abs(v)))
    return lim if lim > 0.0 else 1.0


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

    cum_b = np.cumsum(resid_b)
    cum_t = np.cumsum(resid_t)

    with plt.rc_context(_pub_rc()):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12.0, 7.2), sharex=True)

        ax1.plot(day_b, resid_b, color=_C_BASE, lw=1.0, alpha=0.9, label="Baseline (TSR=OFF)")
        ax1.plot(day_t, resid_t, color=_C_TSR, lw=1.0, alpha=0.9, ls="--", label="TSR=ON")
        ax1.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax1.set_ylabel("Residual (mm / interval)")
        ax1.set_title("Interval residual")

        ax2.plot(day_b, cum_b, color=_C_BASE, lw=1.1, label="Baseline (TSR=OFF)")
        ax2.plot(day_t, cum_t, color=_C_TSR, lw=1.1, ls="--", label="TSR=ON")
        ax2.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax2.set_xlabel("Simulation day (left endpoint)")
        ax2.set_ylabel("Cumulative residual (mm)")
        ax2.set_title("Cumulative residual")

        ax2.text(
            0.02,
            0.98,
            f"Σresid(base) = {float(cum_b[-1]):+.2f} mm\nΣresid(tsr)  = {float(cum_t[-1]):+.2f} mm",
            transform=ax2.transAxes,
            va="top",
            ha="left",
            bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.9, "pad": 5},
        )

        handles, labels = ax1.get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncols=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
        fig.suptitle(title, y=1.08)

        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
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

    with plt.rc_context(_pub_rc()):
        fig, ax = plt.subplots(1, 1, figsize=(10.8, 4.8))
        ax.plot(xb, yb, color=_C_BASE, lw=1.4, label="Baseline (TSR=OFF)")
        ax.plot(xt, yt, color=_C_TSR, lw=1.4, ls="--", label="TSR=ON")
        ax.set_xscale("log")
        ax.set_xlabel(r"$|residual|$ (mm / interval)")
        ax.set_ylabel("CDF")
        ax.set_title(title)
        ax.grid(True, which="major")
        ax.grid(True, which="minor", alpha=0.12)
        ax.legend(loc="lower right", frameon=True)
        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
        plt.close(fig)


def plot_basin_abs_cdf_multi(
    *,
    cases: list[tuple[str, Path]],
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

    def cdf(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        xs = np.sort(x)
        ys = (np.arange(xs.size) + 1) / xs.size
        return xs, ys

    colors = [_C_BASE, _C_TSR, _C_GREEN, _C_PURPLE, "#E69F00", _C_GREY]
    with plt.rc_context(_pub_rc()):
        fig, ax = plt.subplots(1, 1, figsize=(10.8, 4.8))
        for i, (label, out_dir) in enumerate(cases):
            x = read_abs_mm(out_dir)
            xs, ys = cdf(x)
            ax.plot(xs, ys, color=colors[i % len(colors)], lw=1.4, label=label)

        ax.set_xscale("log")
        ax.set_xlabel(r"$|residual|$ (mm / interval)")
        ax.set_ylabel("CDF")
        ax.set_title(title)
        ax.grid(True, which="major")
        ax.grid(True, which="minor", alpha=0.12)
        ax.legend(loc="lower right", frameon=True)
        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
        plt.close(fig)


def plot_basin_abs_cdf_multi_dat(
    *,
    cases: list[tuple[str, Path]],
    mesh_path: Path,
    out_path: Path,
    title: str,
) -> None:
    areas = _read_mesh_areas(mesh_path)
    A = float(areas.sum())

    def read_abs_mm(dat_path: Path) -> np.ndarray:
        _, vals, _ = _read_matrix_fast(dat_path)
        resid_m3 = vals[:, 8]
        resid_mm = resid_m3 / A * 1000.0
        return np.abs(resid_mm)

    def cdf(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        xs = np.sort(x)
        ys = (np.arange(xs.size) + 1) / xs.size
        return xs, ys

    colors = [_C_BASE, _C_TSR, _C_GREEN, _C_PURPLE, "#E69F00", _C_GREY]
    with plt.rc_context(_pub_rc()):
        fig, ax = plt.subplots(1, 1, figsize=(10.8, 4.8))
        for i, (label, dat_path) in enumerate(cases):
            x = read_abs_mm(dat_path)
            xs, ys = cdf(x)
            ax.plot(xs, ys, color=colors[i % len(colors)], lw=1.4, label=label)

        ax.set_xscale("log")
        ax.set_xlabel(r"$|residual|$ (mm / interval)")
        ax.set_ylabel("CDF")
        ax.set_title(title)
        ax.grid(True, which="major")
        ax.grid(True, which="minor", alpha=0.12)
        ax.legend(loc="lower right", frameon=True)
        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
        plt.close(fig)


def plot_basin_cumsum_multi_dat(
    *,
    cases: list[tuple[str, Path]],
    mesh_path: Path,
    out_path: Path,
    title: str,
) -> None:
    areas = _read_mesh_areas(mesh_path)
    A = float(areas.sum())

    colors = [_C_BASE, _C_TSR, _C_GREEN, _C_PURPLE, "#E69F00", _C_GREY]
    with plt.rc_context(_pub_rc()):
        fig, ax = plt.subplots(1, 1, figsize=(10.8, 4.8))
        for i, (label, dat_path) in enumerate(cases):
            t_min, vals, _ = _read_matrix_fast(dat_path)
            resid_mm = vals[:, 8] / A * 1000.0
            day = t_min / 1440.0
            cum = np.cumsum(resid_mm)
            ax.plot(day, cum, color=colors[i % len(colors)], lw=1.4, label=f"{label} (Σ={float(cum[-1]):+.2f} mm)")

        ax.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax.set_xlabel("Simulation day (left endpoint)")
        ax.set_ylabel("Cumulative residual (mm)")
        ax.set_title(title)
        ax.legend(loc="best", frameon=True)
        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
        plt.close(fig)


def plot_basin_residual_zoom_multi_dat(
    *,
    cases: list[tuple[str, Path]],
    mesh_path: Path,
    out_path: Path,
    title: str,
    center_day: int,
    window_days: int = 12,
) -> None:
    areas = _read_mesh_areas(mesh_path)
    A = float(areas.sum())

    colors = [_C_BASE, _C_TSR, _C_GREEN, _C_PURPLE, "#E69F00", _C_GREY]
    with plt.rc_context(_pub_rc()):
        fig, ax = plt.subplots(1, 1, figsize=(10.8, 4.8))
        for i, (label, dat_path) in enumerate(cases):
            t_min, vals, _ = _read_matrix_fast(dat_path)
            resid_mm = vals[:, 8] / A * 1000.0
            day = (t_min / 1440.0).astype(int)

            lo = center_day - window_days
            hi = center_day + window_days
            m = (day >= lo) & (day <= hi)
            ax.plot(day[m], resid_mm[m], color=colors[i % len(colors)], lw=1.4, label=label)

        ax.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax.set_xlabel("Simulation day (left endpoint)")
        ax.set_ylabel("Residual (mm / interval)")
        ax.set_title(title + f" (zoom: Day {center_day-window_days}–{center_day+window_days})")
        ax.legend(loc="upper right", frameon=True)
        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
        plt.close(fig)


def plot_basin_cumsum_multi(
    *,
    cases: list[tuple[str, Path]],
    mesh_path: Path,
    out_path: Path,
    title: str,
) -> None:
    areas = _read_mesh_areas(mesh_path)
    A = float(areas.sum())

    colors = [_C_BASE, _C_TSR, _C_GREEN, _C_PURPLE, "#E69F00", _C_GREY]
    with plt.rc_context(_pub_rc()):
        fig, ax = plt.subplots(1, 1, figsize=(10.8, 4.8))
        for i, (label, out_dir) in enumerate(cases):
            t_min, vals, _ = _read_matrix_fast(out_dir / "ccw.basinwbfull.dat")
            resid_mm = vals[:, 8] / A * 1000.0
            day = t_min / 1440.0
            cum = np.cumsum(resid_mm)
            ax.plot(day, cum, color=colors[i % len(colors)], lw=1.4, label=f"{label} (Σ={float(cum[-1]):+.2f} mm)")

        ax.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax.set_xlabel("Simulation day (left endpoint)")
        ax.set_ylabel("Cumulative residual (mm)")
        ax.set_title(title)
        ax.legend(loc="best", frameon=True)
        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
        plt.close(fig)


def plot_basin_residual_zoom_multi(
    *,
    cases: list[tuple[str, Path]],
    mesh_path: Path,
    out_path: Path,
    title: str,
    center_day: int,
    window_days: int = 12,
) -> None:
    areas = _read_mesh_areas(mesh_path)
    A = float(areas.sum())

    colors = [_C_BASE, _C_TSR, _C_GREEN, _C_PURPLE, "#E69F00", _C_GREY]
    with plt.rc_context(_pub_rc()):
        fig, ax = plt.subplots(1, 1, figsize=(10.8, 4.8))
        for i, (label, out_dir) in enumerate(cases):
            t_min, vals, _ = _read_matrix_fast(out_dir / "ccw.basinwbfull.dat")
            resid_mm = vals[:, 8] / A * 1000.0
            day = (t_min / 1440.0).astype(int)

            lo = center_day - window_days
            hi = center_day + window_days
            m = (day >= lo) & (day <= hi)
            ax.plot(day[m], resid_mm[m], color=colors[i % len(colors)], lw=1.4, label=label)

        ax.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax.set_xlabel("Simulation day (left endpoint)")
        ax.set_ylabel("Residual (mm / interval)")
        ax.set_title(title + f" (zoom: Day {center_day-window_days}–{center_day+window_days})")
        ax.legend(loc="upper right", frameon=True)
        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
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

    with plt.rc_context(_pub_rc()):
        fig, ax = plt.subplots(1, 1, figsize=(10.8, 4.8))
        for x, label, color, ls in [
            (base_a, "Baseline A (solver dS/dt)", _C_BASE, "-"),
            (base_b, "Baseline B (budget RHS)", _C_BASE, "--"),
            (tsr_a, "TSR A (solver dS/dt)", _C_TSR, "-"),
            (tsr_b, "TSR B (budget RHS)", _C_TSR, "--"),
        ]:
            xs, ys = cdf(x)
            ax.plot(xs, ys, color=color, lw=1.4, ls=ls, label=label)

        ax.set_xscale("log")
        ax.set_xlabel(r"$|residual|$ (mm / interval)")
        ax.set_ylabel("CDF (time×element samples)")
        ax.set_title(title)
        ax.grid(True, which="major")
        ax.grid(True, which="minor", alpha=0.12)
        ax.legend(loc="lower right", frameon=True)
        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
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

    with plt.rc_context(_pub_rc()):
        fig, ax = plt.subplots(1, 1, figsize=(10.8, 4.6))
        fracs_pct = [f * 100.0 for f in fracs]
        bars = ax.bar(labels, fracs_pct, color=_C_BASE, alpha=0.9)
        ax.set_ylim(0.0, max(fracs_pct) * 1.25)
        ax.set_ylabel("Area fraction (%)")
        ax.set_title(title)
        ax.grid(True, axis="y")
        ax.grid(False, axis="x")
        for b, f in zip(bars, fracs_pct):
            ax.text(
                b.get_x() + b.get_width() / 2,
                b.get_height(),
                f"{f:.1f}%",
                ha="center",
                va="bottom",
                fontsize=10,
            )
        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
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
    mask_other_all = ~(g.mask_south | g.mask_north)

    def group_means(out_dir: Path, dat_name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        t_min, vals, _ = _read_matrix_fast(out_dir / dat_name)
        south = _aw_mean_series(vals, areas=g.areas, mask=g.mask_south)
        north = _aw_mean_series(vals, areas=g.areas, mask=g.mask_north)
        other = _aw_mean_series(vals, areas=g.areas, mask=mask_other_all)
        day = t_min / 1440.0
        return day, south, north, other

    day_b, rn_b_s, rn_b_n, rn_b_o = group_means(out_base, "ccw.rn_t.dat")
    day_t, rn_t_s, rn_t_n, rn_t_o = group_means(out_tsr, "ccw.rn_t.dat")
    _, et_b_s, et_b_n, et_b_o = group_means(out_base, "ccw.eleveta.dat")
    _, et_t_s, et_t_n, et_t_o = group_means(out_tsr, "ccw.eleveta.dat")

    if day_b.shape != day_t.shape or not np.allclose(day_b, day_t, rtol=0.0, atol=0.0):
        raise PlotError("Baseline and TSR time axes differ; cannot plot aspect TSR effect.")

    rn_eff_s = rn_t_s - rn_b_s
    rn_eff_n = rn_t_n - rn_b_n
    rn_eff_o = rn_t_o - rn_b_o
    et_eff_s = (et_t_s - et_b_s) * 1000.0
    et_eff_n = (et_t_n - et_b_n) * 1000.0
    et_eff_o = (et_t_o - et_b_o) * 1000.0

    with plt.rc_context(_pub_rc()):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12.0, 7.2), sharex=True)

        ax1.plot(day_b, rn_eff_s, color=_C_TSR, lw=1.0, alpha=0.9, label="South (TSR−BASE)")
        ax1.plot(day_b, rn_eff_n, color=_C_BASE, lw=1.0, alpha=0.9, label="North (TSR−BASE)")
        ax1.plot(day_b, rn_eff_o, color=_C_GREY, lw=1.0, alpha=0.9, label="Other (TSR−BASE)")
        ax1.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax1.set_ylabel(r"$\Delta rn_t$ (W m$^{-2}$)  TSR−BASE")
        ax1.set_title("TSR effect on radiation (group means)")

        ax2.plot(day_b, et_eff_s, color=_C_TSR, lw=1.0, alpha=0.9, label="South (TSR−BASE)")
        ax2.plot(day_b, et_eff_n, color=_C_BASE, lw=1.0, alpha=0.9, label="North (TSR−BASE)")
        ax2.plot(day_b, et_eff_o, color=_C_GREY, lw=1.0, alpha=0.9, label="Other (TSR−BASE)")
        ax2.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax2.set_xlabel("Simulation day (left endpoint)")
        ax2.set_ylabel(r"$\Delta ET_a$ (mm / day)  TSR−BASE")
        ax2.set_title("TSR effect on ET (group means)")

        handles, labels = ax1.get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncols=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
        fig.suptitle(title, y=1.08)

        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
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

    with plt.rc_context(_pub_rc()):
        fig, axes = plt.subplots(4, 1, figsize=(12.0, 9.6), sharex=True, height_ratios=[1.0, 1.0, 1.0, 0.9])

        ax = axes[0]
        ax.plot(day_b, du_b, color=_C_BASE, lw=1.0, alpha=0.9, label="Baseline (TSR=OFF)")
        ax.plot(day_t, du_t, color=_C_TSR, lw=1.0, alpha=0.9, ls="--", label="TSR=ON")
        ax.plot(day_t, du_eff, color=_C_GREY, lw=1.2, label="TSR effect (ΔTSR−ΔBASE)")
        ax.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax.set_ylabel(r"$\Delta y_{Unsat}$ (m)  South−North")
        ax.set_title(r"Unsat storage contrast ($South-North$)")

        ax = axes[1]
        ax.plot(day_b, dr_b, color=_C_BASE, lw=1.0, alpha=0.9, label="Baseline (TSR=OFF)")
        ax.plot(day_t, dr_t, color=_C_TSR, lw=1.0, alpha=0.9, ls="--", label="TSR=ON")
        ax.plot(day_t, dr_eff, color=_C_GREY, lw=1.2, label="TSR effect (ΔTSR−ΔBASE)")
        ax.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax.set_ylabel(r"$\Delta (Unsat/Deficit)$ (–)  South−North")
        ax.set_title(r"Normalized unsat contrast ($South-North$)")

        ax = axes[2]
        ax.plot(day_b, dg_b, color=_C_BASE, lw=1.0, alpha=0.9, label="Baseline (TSR=OFF)")
        ax.plot(day_t, dg_t, color=_C_TSR, lw=1.0, alpha=0.9, ls="--", label="TSR=ON")
        ax.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax.set_ylabel(r"$\Delta y_{GW}$ (m)  South−North")
        ax.set_title(r"GW storage contrast ($South-North$)")

        ax = axes[3]
        eff_mm = dg_eff * 1000.0
        ax.plot(day_t, eff_mm, color=_C_GREY, lw=1.2, label="TSR effect (ΔTSR−ΔBASE)")
        ax.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax.set_xlabel("Simulation day (left endpoint)")
        ax.set_ylabel(r"$\Delta\Delta y_{GW}$ (mm)")
        ax.set_title(r"TSR effect on GW contrast ($\Delta_{TSR}-\Delta_{BASE}$)")

        handles, labels = axes[0].get_legend_handles_labels()
        handles2, labels2 = axes[3].get_legend_handles_labels()
        fig.legend(
            handles + handles2,
            labels + labels2,
            loc="upper center",
            ncols=3,
            frameon=False,
            bbox_to_anchor=(0.5, 1.02),
        )
        fig.suptitle(title, y=1.07)

        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
        plt.close(fig)


def plot_storage_delta_spatial(
    *,
    out_base: Path,
    out_tsr: Path,
    mesh_path: Path,
    out_path: Path,
    title: str,
) -> None:
    tri, ele_ids = _read_mesh_triangulation(mesh_path)

    def delta_mm(dat_name: str) -> np.ndarray:
        _, vb, _ = _read_matrix_fast(out_base / dat_name)
        _, vt, _ = _read_matrix_fast(out_tsr / dat_name)
        if vb.shape != vt.shape:
            raise PlotError(f"shape mismatch for {dat_name}: {vb.shape} vs {vt.shape}")
        d = (_nanmean_time(vt) - _nanmean_time(vb)) * 1000.0
        return d[ele_ids - 1]

    du = delta_mm("ccw.eleyunsat.dat")
    dg = delta_mm("ccw.eleygw.dat")

    panels = [
        (du, r"$\Delta y_{Unsat}$ (TSR−BASE)", r"$\Delta y_{Unsat}$ (mm)"),
        (dg, r"$\Delta y_{GW}$ (TSR−BASE)", r"$\Delta y_{GW}$ (mm)"),
    ]

    with plt.rc_context(_pub_rc()):
        fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0), constrained_layout=True)

        for ax, (data, subtitle, cblab) in zip(axes, panels):
            lim = _robust_sym_lim(data, q=99.0)
            norm = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)

            pc = ax.tripcolor(
                tri,
                facecolors=data,
                shading="flat",
                cmap="RdBu_r",
                norm=norm,
                edgecolors="none",
                rasterized=True,
            )
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("X (km)")
            ax.grid(False)
            ax.set_title(subtitle)

            if ax is axes[0]:
                ax.set_ylabel("Y (km)")
            else:
                ax.set_ylabel("")
                ax.tick_params(labelleft=False)

            cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.02)
            cb.set_label(cblab)

            mean = float(np.nanmean(data))
            p99 = float(np.nanpercentile(np.abs(data[np.isfinite(data)]), 99.0))
            ax.text(
                0.02,
                0.02,
                f"mean={mean:+.2f} mm\np99(|·|)={p99:.2f} mm",
                transform=ax.transAxes,
                va="bottom",
                ha="left",
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85),
            )

        fig.suptitle(title, y=1.02)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


def plot_radiation_cumsum_delta_spatial(
    *,
    out_base: Path,
    out_tsr: Path,
    mesh_path: Path,
    out_path: Path,
    title: str,
) -> None:
    """
    Per-element cumulative radiation change due to terrain correction:
    - baseline is the horizontal-plane radiation `rn_h`
    - TSR is `rn_t`

    We integrate interval-mean W/m^2 over the whole run to MJ/m^2, then map:
    - ΔE = E_t - E_h (MJ/m^2)
    - Δ% = (E_t/E_h - 1) * 100 (%)
    """
    tri, ele_ids = _read_mesh_triangulation(mesh_path)

    t_h, rn_h, _ = _read_matrix_fast(out_base / "ccw.rn_h.dat")
    t_t, rn_t, _ = _read_matrix_fast(out_tsr / "ccw.rn_t.dat")
    if rn_h.shape != rn_t.shape:
        raise PlotError(f"shape mismatch for rn_h vs rn_t: {rn_h.shape} vs {rn_t.shape}")
    if t_h.shape != t_t.shape or not np.allclose(t_h, t_t):
        raise PlotError("time axis mismatch for rn_h vs rn_t")

    E_h = _nanintegrate_time(rn_h, t_h)  # J/m^2
    E_t = _nanintegrate_time(rn_t, t_t)  # J/m^2

    dE = (E_t - E_h) / 1e6  # MJ/m^2
    dPct = np.where(E_h > 0.0, (E_t / E_h - 1.0) * 100.0, np.nan)  # %

    dE = dE[ele_ids - 1]
    dPct = dPct[ele_ids - 1]

    panels = [
        (dE, r"$\Delta E = \int (rn_t - rn_h)\,dt$", r"$\Delta E$ (MJ m$^{-2}$)"),
        (dPct, r"$\Delta\% = (E_t/E_h - 1)\times 100$", r"$\Delta\%$ (%)"),
    ]

    with plt.rc_context(_pub_rc()):
        fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0), constrained_layout=True)

        for ax, (data, subtitle, cblab) in zip(axes, panels):
            lim = _robust_sym_lim(data, q=99.0)
            if not (lim > 0.0):
                lim = 1.0
            norm = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)

            pc = ax.tripcolor(
                tri,
                facecolors=data,
                shading="flat",
                cmap="RdBu_r",
                norm=norm,
                edgecolors="none",
                rasterized=True,
            )
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("X (km)")
            ax.grid(False)
            ax.set_title(subtitle)

            if ax is axes[0]:
                ax.set_ylabel("Y (km)")
            else:
                ax.set_ylabel("")
                ax.tick_params(labelleft=False)

            cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.02)
            cb.set_label(cblab)

            vv = data[np.isfinite(data)]
            if vv.size:
                mean = float(np.nanmean(vv))
                p99 = float(np.nanpercentile(np.abs(vv), 99.0))
                ax.text(
                    0.02,
                    0.02,
                    f"mean={mean:+.2f}\np99(|·|)={p99:.2f}",
                    transform=ax.transAxes,
                    va="bottom",
                    ha="left",
                    fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85),
                )

        fig.suptitle(title, y=1.02)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate plots for docs/water_balance_verification.md")
    ap.add_argument("--mesh", type=Path, default=Path("input/ccw/ccw.sp.mesh"))
    ap.add_argument("--base", type=Path, default=Path("output/ccw.base"))
    ap.add_argument("--tsr", type=Path, default=Path("output/ccw.tsr"))
    ap.add_argument("--tsr-ms5", type=Path, default=Path("output/ccw.tsr.ms5"))
    ap.add_argument("--tsr-ms2", type=Path, default=Path("output/ccw.tsr.ms2"))
    ap.add_argument("--base-2d", type=Path, default=Path("output/ccw.base.2d"))
    ap.add_argument("--tsr-2d", type=Path, default=Path("output/ccw.tsr.2d"))
    ap.add_argument("--wb-be", type=Path, default=Path("output/ccw.tsr.wb.be"))
    ap.add_argument("--wb-trapz", type=Path, default=Path("output/ccw.tsr.wb.trapz"))
    ap.add_argument("--wb-trapzquad", type=Path, default=Path("output/ccw.tsr.wb.trapzquad"))
    ap.add_argument("--outdir", type=Path, default=Path("docs/figures/water_balance"))
    ap.add_argument("--slope-deg-min", type=float, default=5.0)
    ap.add_argument("--aspect-halfwidth-deg", type=float, default=45.0)
    args = ap.parse_args()

    _setup_fonts()

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

    # MAX_SOLVER_STEP convergence (TSR=ON) if experiment outputs exist.
    conv_cases: list[tuple[str, Path]] = []
    if (args.tsr / "ccw.basinwbfull.dat").exists():
        conv_cases.append(("MAX_SOLVER_STEP=10 min", args.tsr))
    if (args.tsr_ms5 / "ccw.basinwbfull.dat").exists():
        conv_cases.append(("MAX_SOLVER_STEP=5 min", args.tsr_ms5))
    if (args.tsr_ms2 / "ccw.basinwbfull.dat").exists():
        conv_cases.append(("MAX_SOLVER_STEP=2 min", args.tsr_ms2))

    if len(conv_cases) >= 2:
        plot_basin_abs_cdf_multi(
            cases=conv_cases,
            mesh_path=args.mesh,
            out_path=args.outdir / "basin_residual_convergence_abs_cdf.png",
            title="Basin |residual| CDF vs MAX_SOLVER_STEP (TSR=ON)",
        )

        plot_basin_cumsum_multi(
            cases=conv_cases,
            mesh_path=args.mesh,
            out_path=args.outdir / "basin_residual_convergence_cumsum.png",
            title="Basin cumulative residual vs MAX_SOLVER_STEP (TSR=ON)",
        )

        # Center around the largest residual day in the default (10-min) run.
        areas = _read_mesh_areas(args.mesh)
        A = float(areas.sum())
        t_min, vals, _ = _read_matrix_fast(args.tsr / "ccw.basinwbfull.dat")
        resid_mm = vals[:, 8] / A * 1000.0
        center_day = int(t_min[int(np.argmax(np.abs(resid_mm)))] / 1440.0)
        plot_basin_residual_zoom_multi(
            cases=conv_cases,
            mesh_path=args.mesh,
            out_path=args.outdir / "basin_residual_convergence_zoom.png",
            title="Basin residual vs MAX_SOLVER_STEP (TSR=ON)",
            center_day=center_day,
            window_days=12,
        )

    # Diagnostic integrator comparison (BE vs TRAPZ vs internal-step QUAD) if outputs exist.
    wb_be = args.wb_be / "ccw.basinwbfull.dat"
    wb_trapz = args.wb_trapz / "ccw.basinwbfull.dat"
    wb_quad_sampling = args.wb_trapzquad / "ccw.basinwbfull.dat"
    wb_quad_internal = args.wb_trapzquad / "ccw.basinwbfull_quad.dat"

    if wb_be.exists() and wb_trapz.exists() and wb_quad_internal.exists():
        cases_int = [
            ("Backward Euler (sampled)", wb_be),
            ("Trapezoidal (sampled)", wb_trapz),
            ("Internal-step (quad)", wb_quad_internal),
        ]
        plot_basin_abs_cdf_multi_dat(
            cases=cases_int,
            mesh_path=args.mesh,
            out_path=args.outdir / "basin_residual_integrators_abs_cdf.png",
            title="Basin |residual| CDF vs diagnostic integrator",
        )
        plot_basin_cumsum_multi_dat(
            cases=cases_int,
            mesh_path=args.mesh,
            out_path=args.outdir / "basin_residual_integrators_cumsum.png",
            title="Basin cumulative residual vs diagnostic integrator",
        )

        areas = _read_mesh_areas(args.mesh)
        A = float(areas.sum())
        t_min, vals, _ = _read_matrix_fast(wb_be)
        resid_mm = vals[:, 8] / A * 1000.0
        center_day = int(t_min[int(np.argmax(np.abs(resid_mm)))] / 1440.0)
        plot_basin_residual_zoom_multi_dat(
            cases=cases_int,
            mesh_path=args.mesh,
            out_path=args.outdir / "basin_residual_integrators_zoom.png",
            title="Basin residual vs diagnostic integrator",
            center_day=center_day,
            window_days=12,
        )

    if wb_quad_sampling.exists() and wb_quad_internal.exists():
        cases_quad = [
            ("Trapezoidal (same run)", wb_quad_sampling),
            ("Internal-step (same run)", wb_quad_internal),
        ]
        plot_basin_abs_cdf_multi_dat(
            cases=cases_quad,
            mesh_path=args.mesh,
            out_path=args.outdir / "basin_residual_quadrun_abs_cdf.png",
            title="Basin |residual| CDF (quad run: sampled vs internal-step)",
        )
        plot_basin_cumsum_multi_dat(
            cases=cases_quad,
            mesh_path=args.mesh,
            out_path=args.outdir / "basin_residual_quadrun_cumsum.png",
            title="Basin cumulative residual (quad run: sampled vs internal-step)",
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
        title="Aspect-group TSR effect time series (TSR−BASE; full run)",
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
    plot_storage_delta_spatial(
        out_base=args.base,
        out_tsr=args.tsr,
        mesh_path=args.mesh,
        out_path=args.outdir / "storage_delta_spatial.png",
        title="Spatial mean storage response (TSR−BASE)",
    )
    plot_radiation_cumsum_delta_spatial(
        out_base=args.base,
        out_tsr=args.tsr,
        mesh_path=args.mesh,
        out_path=args.outdir / "radiation_cumsum_delta_spatial.png",
        title="Cumulative terrain radiation change (TSR vs plane)",
    )

    print(f"OK: wrote figures under {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
