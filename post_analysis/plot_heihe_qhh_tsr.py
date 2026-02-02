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


def _read_mesh_aq_depth(mesh_path: Path) -> np.ndarray:
    """
    Return element-mean aquifer depth (m), averaged from node AqDepth column.
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

    aqd = np.zeros(num_node + 1)
    for k in range(num_node):
        cols = lines[node_start + k].split()
        if len(cols) < 4:
            raise PlotError(f"invalid node row in {mesh_path}: {lines[node_start+k]!r}")
        nid = int(cols[0])
        aqd[nid] = float(cols[3])

    aqd_ele = np.zeros(num_ele)
    for eid, n1, n2, n3 in elements:
        aqd_ele[eid - 1] = (aqd[n1] + aqd[n2] + aqd[n3]) / 3.0
    return aqd_ele


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

    return AspectGroups(mask_south=south, mask_north=north, mask_flat=flat, mask_other=other, areas=areas)


def _aw_mean_series(values: np.ndarray, *, areas: np.ndarray, mask: np.ndarray) -> np.ndarray:
    w = areas * mask.astype(float)
    s = float(w.sum())
    if not (s > 0.0):
        raise PlotError("empty group (area sum=0)")
    return values.dot(w) / s


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

    labels = ["South", "North", f"Flat(<{slope_deg_min:g}°)", "Other(aspect)"]
    masks = [g.mask_south, g.mask_north, g.mask_flat, g.mask_other]
    fracs = [float(g.areas[m].sum() / A) for m in masks]
    fracs_pct = [f * 100.0 for f in fracs]

    with plt.rc_context(_pub_rc()):
        fig, ax = plt.subplots(1, 1, figsize=(10.8, 4.6))
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
    basin: str,
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

    day_b, rn_b = delta_series(out_base, f"{basin}.rn_t.dat")
    day_t, rn_t = delta_series(out_tsr, f"{basin}.rn_t.dat")
    _, et_b_mday = delta_series(out_base, f"{basin}.eleveta.dat")
    _, et_t_mday = delta_series(out_tsr, f"{basin}.eleveta.dat")
    with plt.rc_context(_pub_rc()):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12.0, 7.2), sharex=True)

        ax1.plot(day_b, rn_b, color=_C_BASE, lw=1.0, alpha=0.9, label="Baseline (TSR=OFF)")
        ax1.plot(day_t, rn_t, color=_C_TSR, lw=1.0, alpha=0.9, ls="--", label="TSR=ON")
        ax1.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax1.set_ylabel(r"$\Delta rn_t$ (W m$^{-2}$)  South−North")
        ax1.set_title(r"Radiation contrast ($South-North$)")

        ax2.plot(day_b, et_b_mday * 1000.0, color=_C_BASE, lw=1.0, alpha=0.9, label="Baseline (TSR=OFF)")
        ax2.plot(day_t, et_t_mday * 1000.0, color=_C_TSR, lw=1.0, alpha=0.9, ls="--", label="TSR=ON")
        ax2.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax2.set_xlabel("Simulation day (left endpoint)")
        ax2.set_ylabel(r"$\Delta ET_a$ (mm / day)  South−North")
        ax2.set_title(r"ET contrast ($South-North$)")

        handles, labels = ax1.get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncols=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
        fig.suptitle(title, y=1.08)

        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
        plt.close(fig)


def plot_aspect_storage_delta_timeseries(
    *,
    basin: str,
    out_base: Path,
    out_tsr: Path,
    mesh_path: Path,
    slope_deg_min: float,
    aspect_halfwidth_deg: float,
    out_path: Path,
    title: str,
) -> None:
    aqd_ele = _read_mesh_aq_depth(mesh_path)
    g = build_aspect_groups(mesh_path=mesh_path, slope_deg_min=slope_deg_min, aspect_halfwidth_deg=aspect_halfwidth_deg)

    def delta_series(out_dir: Path, dat_name: str) -> tuple[np.ndarray, np.ndarray]:
        t_min, vals, _ = _read_matrix_fast(out_dir / dat_name)
        south = _aw_mean_series(vals, areas=g.areas, mask=g.mask_south)
        north = _aw_mean_series(vals, areas=g.areas, mask=g.mask_north)
        day = t_min / 1440.0
        return day, (south - north)

    def delta_unsat_ratio(out_dir: Path) -> tuple[np.ndarray, np.ndarray]:
        t_min, uns, _ = _read_matrix_fast(out_dir / f"{basin}.eleyunsat.dat")
        _, gw, _ = _read_matrix_fast(out_dir / f"{basin}.eleygw.dat")
        deficit = aqd_ele.reshape((1, -1)) - gw
        ratio = np.where(deficit > 1e-12, uns / deficit, np.nan)

        def wmean(mask: np.ndarray) -> np.ndarray:
            m = mask.reshape((1, -1)) & np.isfinite(ratio)
            w = g.areas.reshape((1, -1)) * m.astype(float)
            ws = w.sum(axis=1)
            return (ratio * w).sum(axis=1) / ws

        south = wmean(g.mask_south)
        north = wmean(g.mask_north)
        day = t_min / 1440.0
        return day, (south - north)

    day_b, du_b = delta_series(out_base, f"{basin}.eleyunsat.dat")
    day_t, du_t = delta_series(out_tsr, f"{basin}.eleyunsat.dat")
    _, dg_b = delta_series(out_base, f"{basin}.eleygw.dat")
    _, dg_t = delta_series(out_tsr, f"{basin}.eleygw.dat")
    _, dr_b = delta_unsat_ratio(out_base)
    _, dr_t = delta_unsat_ratio(out_tsr)

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
    basin: str,
    out_base: Path,
    out_tsr: Path,
    mesh_path: Path,
    out_path: Path,
    title: str,
) -> None:
    tri, ele_ids = _read_mesh_triangulation(mesh_path)

    def delta_mm(var: str) -> np.ndarray:
        _, vb, _ = _read_matrix_fast(out_base / f"{basin}.{var}.dat")
        _, vt, _ = _read_matrix_fast(out_tsr / f"{basin}.{var}.dat")
        if vb.shape != vt.shape:
            raise PlotError(f"shape mismatch for {basin}.{var}: {vb.shape} vs {vt.shape}")
        d = (_nanmean_time(vt) - _nanmean_time(vb)) * 1000.0
        return d[ele_ids - 1]

    du = delta_mm("eleyunsat")
    dg = delta_mm("eleygw")

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
    basin: str,
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

    t_h, rn_h, _ = _read_matrix_fast(out_base / f"{basin}.rn_h.dat")
    t_t, rn_t, _ = _read_matrix_fast(out_tsr / f"{basin}.rn_t.dat")
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


def plot_basin_residual_timeseries(
    *,
    basin: str,
    out_base: Path,
    out_tsr: Path,
    mesh_path: Path,
    out_path: Path,
    title: str,
) -> None:
    areas = _read_mesh_areas(mesh_path)
    A = float(areas.sum())

    def read_resid_mm(out_dir: Path) -> tuple[np.ndarray, np.ndarray]:
        t_min, vals, _ = _read_matrix_fast(out_dir / f"{basin}.basinwbfull.dat")
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


def _read_lake_init_stage(cfg_ic: Path) -> float:
    lines = cfg_ic.read_text(encoding="utf-8", errors="ignore").splitlines()
    for i, ln in enumerate(lines):
        if "LakeStage" in ln:
            # next non-empty line should be: <idx> <stage>
            for j in range(i + 1, len(lines)):
                parts = lines[j].split()
                if len(parts) >= 2:
                    return float(parts[1])
    raise PlotError(f"failed to parse initial LakeStage from {cfg_ic}")


def lake_corrected_residual_mm(
    *, out_dir: Path, basin: str, mesh_area_m2: float, cfg_ic: Path
) -> tuple[np.ndarray, np.ndarray]:
    """
    Augment basinwbfull residual with lake terms (qhh has lake module enabled):
      resid_full ~= resid_raw + dS_lake + E_lake

    Why no "-P_lake"?
    - basinwbfull.P is computed as sum(qElePrep * area) over ALL elements, including lake elements,
      so lake precipitation is already included in basin P.
    - basinwbfull.ET is computed from (qEleE_IC_raw + qEs + qEu + qEg + qTu + qTg); for lake elements
      those terms are 0, so open-water evaporation is missing and must be added back (E_lake).
    - basinwbfull storage only includes element Sfull + river storage; lake storage change must be added (dS_lake).

    Notes:
    - Lake outputs are in depth (m) for P/E, and stage (m) + top area (m^2).
    - dS_lake is approximated via trapezoid in (area, stage).
    """
    t_min, wb, _ = _read_matrix_fast(out_dir / f"{basin}.basinwbfull.dat")
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
    day = t_min / 1440.0
    return day, resid_full_mm


def plot_qhh_residual_with_lake_correction(
    *,
    out_base: Path,
    out_tsr: Path,
    mesh_path: Path,
    cfg_ic: Path,
    out_path: Path,
    title: str,
) -> None:
    areas = _read_mesh_areas(mesh_path)
    A = float(areas.sum())

    def read_raw(out_dir: Path) -> tuple[np.ndarray, np.ndarray]:
        t_min, vals, _ = _read_matrix_fast(out_dir / "qhh.basinwbfull.dat")
        resid_m3 = vals[:, 8]
        return t_min / 1440.0, resid_m3 / A * 1000.0

    day_b, raw_b = read_raw(out_base)
    day_t, raw_t = read_raw(out_tsr)
    day_bc, corr_b = lake_corrected_residual_mm(out_dir=out_base, basin="qhh", mesh_area_m2=A, cfg_ic=cfg_ic)
    day_tc, corr_t = lake_corrected_residual_mm(out_dir=out_tsr, basin="qhh", mesh_area_m2=A, cfg_ic=cfg_ic)

    if not np.array_equal(day_b, day_bc) or not np.array_equal(day_t, day_tc):
        raise PlotError("time axis mismatch between raw and corrected residuals")

    # Publication-style figure: separate raw vs corrected panels to avoid over-plotting,
    # and avoid mixing vastly different cumulative scales in one axis.
    base_color = _C_BASE
    tsr_color = _C_TSR

    raw_cum_b = np.cumsum(raw_b)
    raw_cum_t = np.cumsum(raw_t)
    corr_cum_b = np.cumsum(corr_b)
    corr_cum_t = np.cumsum(corr_t)

    sum_raw_b = float(raw_cum_b[-1])
    sum_raw_t = float(raw_cum_t[-1])
    sum_corr_b = float(corr_cum_b[-1])
    sum_corr_t = float(corr_cum_t[-1])

    with plt.rc_context(_pub_rc()):
        fig = plt.figure(figsize=(13.0, 7.5))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], wspace=0.16, hspace=0.28)

        ax_raw = fig.add_subplot(gs[0, 0])
        ax_corr = fig.add_subplot(gs[0, 1], sharex=ax_raw, sharey=ax_raw)
        ax_raw_c = fig.add_subplot(gs[1, 0], sharex=ax_raw)
        ax_corr_c = fig.add_subplot(gs[1, 1], sharex=ax_raw)

        # Interval residuals
        ax_raw.plot(day_b, raw_b, color=base_color, lw=0.9, alpha=0.85, label="Baseline (TSR=OFF)")
        ax_raw.plot(day_t, raw_t, color=tsr_color, lw=0.9, alpha=0.85, ls="--", label="TSR=ON")
        ax_raw.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax_raw.set_title("Raw diagnostic residual (basinwbfull)")
        ax_raw.set_ylabel("Residual (mm / day)")
        ax_raw.grid(True, alpha=0.25)

        ax_corr.plot(day_b, corr_b, color=base_color, lw=0.9, alpha=0.9, label="Baseline (TSR=OFF)")
        ax_corr.plot(day_t, corr_t, color=tsr_color, lw=0.9, alpha=0.9, ls="--", label="TSR=ON")
        ax_corr.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax_corr.set_title(r"Corrected residual = raw + $\Delta S_{lake}$ + $E_{lake}$")
        ax_corr.grid(True, alpha=0.25)

        # Cumulative residuals
        ax_raw_c.plot(day_b, raw_cum_b, color=base_color, lw=1.1, label="Baseline (TSR=OFF)")
        ax_raw_c.plot(day_t, raw_cum_t, color=tsr_color, lw=1.1, ls="--", label="TSR=ON")
        ax_raw_c.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax_raw_c.set_title("Cumulative raw residual")
        ax_raw_c.set_xlabel("Simulation day (left endpoint)")
        ax_raw_c.set_ylabel("Cumulative residual (mm)")
        ax_raw_c.grid(True, alpha=0.25)
        ax_raw_c.text(
            0.02,
            0.98,
            f"Σresid(base) = {sum_raw_b:+.1f} mm\nΣresid(tsr)  = {sum_raw_t:+.1f} mm",
            transform=ax_raw_c.transAxes,
            va="top",
            ha="left",
            bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.9, "pad": 5},
        )

        ax_corr_c.plot(day_b, corr_cum_b, color=base_color, lw=1.1, label="Baseline (TSR=OFF)")
        ax_corr_c.plot(day_t, corr_cum_t, color=tsr_color, lw=1.1, ls="--", label="TSR=ON")
        ax_corr_c.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax_corr_c.set_title(r"Cumulative corrected residual (+$\Delta S_{lake}$ + $E_{lake}$)")
        ax_corr_c.set_xlabel("Simulation day (left endpoint)")
        ax_corr_c.grid(True, alpha=0.25)
        ax_corr_c.text(
            0.02,
            0.98,
            f"Σresid(base) = {sum_corr_b:+.1f} mm\nΣresid(tsr)  = {sum_corr_t:+.1f} mm",
            transform=ax_corr_c.transAxes,
            va="top",
            ha="left",
            bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.9, "pad": 5},
        )

        # Common legend
        handles, labels = ax_raw.get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncols=2, frameon=False, bbox_to_anchor=(0.5, 0.99))
        fig.suptitle(title, y=1.03)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
        plt.close(fig)


def plot_basin_component_totals(
    *,
    basin: str,
    out_base: Path,
    out_tsr: Path,
    mesh_path: Path,
    out_path: Path,
    title: str,
) -> None:
    """
    Plot whole-run totals (mm) for: P, ET, Qout, ΔS, Resid.
    """
    areas = _read_mesh_areas(mesh_path)
    A = float(areas.sum())

    def totals_mm(out_dir: Path) -> dict[str, float]:
        _, vals, _ = _read_matrix_fast(out_dir / f"{basin}.basinwbfull.dat")
        total = vals.sum(axis=0)  # m3
        mm = total / A * 1000.0
        keys = ["dS", "P", "ET", "Qout", "Qedge", "Qbc", "Qss", "NonConsEdge", "Resid"]
        return {k: float(v) for k, v in zip(keys, mm)}

    b = totals_mm(out_base)
    t = totals_mm(out_tsr)

    big = ["P", "ET", "Qout", "dS"]
    small = ["Resid"]

    with plt.rc_context(_pub_rc()):
        fig = plt.figure(figsize=(12.0, 4.8))
        gs = fig.add_gridspec(1, 2, width_ratios=[4.0, 1.2], wspace=0.25)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])

        w = 0.38
        x1 = np.arange(len(big))
        ax1.bar(x1 - w / 2, [b[c] for c in big], width=w, color=_C_BASE, alpha=0.85, label="Baseline (TSR=OFF)")
        ax1.bar(x1 + w / 2, [t[c] for c in big], width=w, color=_C_TSR, alpha=0.85, label="TSR=ON")
        ax1.set_xticks(x1)
        ax1.set_xticklabels(big)
        ax1.set_ylabel("Total (mm over mesh area)")
        ax1.set_title("Major terms")
        ax1.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax1.grid(True, axis="y")
        ax1.grid(False, axis="x")

        for xi, c in enumerate(big):
            for dx, v in [(-w / 2, b[c]), (w / 2, t[c])]:
                ax1.text(
                    xi + dx,
                    v,
                    f"{v:.0f}",
                    ha="center",
                    va="bottom" if v >= 0 else "top",
                    fontsize=9,
                )

        x2 = np.arange(len(small))
        ax2.bar(x2 - w / 2, [b[c] for c in small], width=w, color=_C_BASE, alpha=0.85)
        ax2.bar(x2 + w / 2, [t[c] for c in small], width=w, color=_C_TSR, alpha=0.85)
        ax2.set_xticks(x2)
        ax2.set_xticklabels(small)
        ax2.set_title("Residual")
        ax2.axhline(0.0, color="k", lw=0.8, alpha=0.5)
        ax2.grid(True, axis="y")
        ax2.grid(False, axis="x")

        resid_min = min(b["Resid"], t["Resid"])
        resid_max = max(b["Resid"], t["Resid"])
        resid_span = max(abs(resid_min), abs(resid_max), 1e-6)
        ax2.set_ylim(-1.25 * resid_span, 1.25 * resid_span)

        for dx, v in [(-w / 2, b["Resid"]), (w / 2, t["Resid"])]:
            ax2.text(
                0 + dx,
                v,
                f"{v:+.2f}",
                ha="center",
                va="bottom" if v >= 0 else "top",
                fontsize=9,
            )

        handles, labels = ax1.get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncols=2, frameon=False, bbox_to_anchor=(0.5, 1.05))
        fig.suptitle(title, y=1.12)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=220, bbox_inches="tight")
        plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate TSR + water-balance figures for heihe and qhh basins.")
    ap.add_argument("--outdir", type=Path, default=Path("docs/figures/tsr_heihe_qhh"))
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
    _setup_fonts()

    # heihe
    heihe_dir = args.outdir / "heihe"
    plot_aspect_group_fractions(
        mesh_path=args.heihe_mesh,
        slope_deg_min=args.slope_deg_min,
        aspect_halfwidth_deg=args.aspect_halfwidth_deg,
        out_path=heihe_dir / "aspect_group_area_fractions.png",
        title="heihe: aspect group area fractions",
    )
    plot_aspect_delta_timeseries(
        basin="heihe",
        out_base=args.heihe_base,
        out_tsr=args.heihe_tsr,
        mesh_path=args.heihe_mesh,
        slope_deg_min=args.slope_deg_min,
        aspect_halfwidth_deg=args.aspect_halfwidth_deg,
        out_path=heihe_dir / "aspect_delta_timeseries.png",
        title="heihe: aspect response (South−North) time series",
    )
    plot_aspect_storage_delta_timeseries(
        basin="heihe",
        out_base=args.heihe_base,
        out_tsr=args.heihe_tsr,
        mesh_path=args.heihe_mesh,
        slope_deg_min=args.slope_deg_min,
        aspect_halfwidth_deg=args.aspect_halfwidth_deg,
        out_path=heihe_dir / "aspect_storage_delta_timeseries.png",
        title="heihe: storage aspect response (South−North)",
    )
    plot_storage_delta_spatial(
        basin="heihe",
        out_base=args.heihe_base,
        out_tsr=args.heihe_tsr,
        mesh_path=args.heihe_mesh,
        out_path=heihe_dir / "storage_delta_spatial.png",
        title="heihe: spatial mean storage response (TSR−BASE)",
    )
    plot_radiation_cumsum_delta_spatial(
        basin="heihe",
        out_base=args.heihe_base,
        out_tsr=args.heihe_tsr,
        mesh_path=args.heihe_mesh,
        out_path=heihe_dir / "radiation_cumsum_delta_spatial.png",
        title="heihe: cumulative terrain radiation change (TSR vs plane)",
    )
    plot_basin_component_totals(
        basin="heihe",
        out_base=args.heihe_base,
        out_tsr=args.heihe_tsr,
        mesh_path=args.heihe_mesh,
        out_path=heihe_dir / "basin_component_totals.png",
        title="heihe: whole-run water balance totals",
    )
    plot_basin_residual_timeseries(
        basin="heihe",
        out_base=args.heihe_base,
        out_tsr=args.heihe_tsr,
        mesh_path=args.heihe_mesh,
        out_path=heihe_dir / "basin_residual_timeseries.png",
        title="heihe: basin water-balance residual",
    )

    # qhh
    qhh_dir = args.outdir / "qhh"
    plot_aspect_group_fractions(
        mesh_path=args.qhh_mesh,
        slope_deg_min=args.slope_deg_min,
        aspect_halfwidth_deg=args.aspect_halfwidth_deg,
        out_path=qhh_dir / "aspect_group_area_fractions.png",
        title="qhh: aspect group area fractions",
    )
    plot_aspect_delta_timeseries(
        basin="qhh",
        out_base=args.qhh_base,
        out_tsr=args.qhh_tsr,
        mesh_path=args.qhh_mesh,
        slope_deg_min=args.slope_deg_min,
        aspect_halfwidth_deg=args.aspect_halfwidth_deg,
        out_path=qhh_dir / "aspect_delta_timeseries.png",
        title="qhh: aspect response (South−North) time series",
    )
    plot_aspect_storage_delta_timeseries(
        basin="qhh",
        out_base=args.qhh_base,
        out_tsr=args.qhh_tsr,
        mesh_path=args.qhh_mesh,
        slope_deg_min=args.slope_deg_min,
        aspect_halfwidth_deg=args.aspect_halfwidth_deg,
        out_path=qhh_dir / "aspect_storage_delta_timeseries.png",
        title="qhh: storage aspect response (South−North)",
    )
    plot_storage_delta_spatial(
        basin="qhh",
        out_base=args.qhh_base,
        out_tsr=args.qhh_tsr,
        mesh_path=args.qhh_mesh,
        out_path=qhh_dir / "storage_delta_spatial.png",
        title="qhh: spatial mean storage response (TSR−BASE)",
    )
    plot_radiation_cumsum_delta_spatial(
        basin="qhh",
        out_base=args.qhh_base,
        out_tsr=args.qhh_tsr,
        mesh_path=args.qhh_mesh,
        out_path=qhh_dir / "radiation_cumsum_delta_spatial.png",
        title="qhh: cumulative terrain radiation change (TSR vs plane)",
    )
    plot_basin_component_totals(
        basin="qhh",
        out_base=args.qhh_base,
        out_tsr=args.qhh_tsr,
        mesh_path=args.qhh_mesh,
        out_path=qhh_dir / "basin_component_totals.png",
        title="qhh: whole-run water balance totals",
    )
    plot_basin_residual_timeseries(
        basin="qhh",
        out_base=args.qhh_base,
        out_tsr=args.qhh_tsr,
        mesh_path=args.qhh_mesh,
        out_path=qhh_dir / "basin_residual_timeseries.png",
        title="qhh: basin water-balance residual (raw diagnostic)",
    )
    plot_qhh_residual_with_lake_correction(
        out_base=args.qhh_base,
        out_tsr=args.qhh_tsr,
        mesh_path=args.qhh_mesh,
        cfg_ic=args.qhh_ic,
        out_path=qhh_dir / "basin_residual_timeseries_lake_corrected.png",
        title="qhh: basin residual (raw) vs +lake terms",
    )

    print(f"OK: wrote figures under {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
