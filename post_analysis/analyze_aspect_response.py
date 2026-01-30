#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

# Reuse the repo's SHUD .dat reader + mesh normal computation (kept C++-consistent).
# validation/tsr/py is not a Python package, so add it to sys.path explicitly.
_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "validation" / "tsr" / "py"))
from shud_reader import DatReader  # noqa: E402
from tsr_core import read_mesh_normals  # noqa: E402


class AspectAnalysisError(RuntimeError):
    pass


@dataclass(frozen=True)
class AspectGroups:
    mask_south: np.ndarray
    mask_north: np.ndarray
    areas: np.ndarray
    aq_depth_ele: np.ndarray
    z_surf_ele: np.ndarray


def _circular_diff(a: np.ndarray, center: float) -> np.ndarray:
    return np.abs(((a - center + math.pi) % (2.0 * math.pi)) - math.pi)


def _read_mesh_geom(mesh_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = mesh_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise AspectAnalysisError(f"empty mesh file: {mesh_path}")
    num_ele = int(lines[0].split()[0])

    elements: list[tuple[int, int, int, int]] = []
    for k in range(num_ele):
        cols = lines[2 + k].split()
        if len(cols) < 4:
            raise AspectAnalysisError(f"invalid element row in {mesh_path}: {lines[2+k]!r}")
        elements.append((int(cols[0]), int(cols[1]), int(cols[2]), int(cols[3])))

    node_hdr_idx = 2 + num_ele
    num_node = int(lines[node_hdr_idx].split()[0])
    node_start = node_hdr_idx + 2

    x = np.zeros(num_node + 1)
    y = np.zeros(num_node + 1)
    aqd = np.zeros(num_node + 1)
    zmax = np.zeros(num_node + 1)
    for k in range(num_node):
        cols = lines[node_start + k].split()
        if len(cols) < 5:
            raise AspectAnalysisError(f"invalid node row in {mesh_path}: {lines[node_start+k]!r}")
        nid = int(cols[0])
        x[nid] = float(cols[1])
        y[nid] = float(cols[2])
        aqd[nid] = float(cols[3])
        zmax[nid] = float(cols[4])

    areas = np.zeros(num_ele)
    aqd_ele = np.zeros(num_ele)
    zsurf_ele = np.zeros(num_ele)
    for eid, n1, n2, n3 in elements:
        x1, y1 = x[n1], y[n1]
        x2, y2 = x[n2], y[n2]
        x3, y3 = x[n3], y[n3]
        areas[eid - 1] = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
        aqd_ele[eid - 1] = (aqd[n1] + aqd[n2] + aqd[n3]) / 3.0
        zsurf_ele[eid - 1] = (zmax[n1] + zmax[n2] + zmax[n3]) / 3.0

    return areas, aqd_ele, zsurf_ele


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


def build_aspect_groups(
    *,
    mesh_path: Path,
    slope_deg_min: float = 5.0,
    aspect_halfwidth_deg: float = 45.0,
) -> AspectGroups:
    areas, aq_depth_ele, z_surf_ele = _read_mesh_geom(mesh_path)
    aspect, slope = _compute_aspect_slope(read_mesh_normals(mesh_path))

    slope_thr = math.radians(float(slope_deg_min))
    half = math.radians(float(aspect_halfwidth_deg))
    valid = slope >= slope_thr

    mask_south = valid & (_circular_diff(aspect, math.pi) <= half)
    mask_north = valid & (_circular_diff(aspect, 0.0) <= half)

    if not bool(np.any(mask_south)):
        raise AspectAnalysisError("empty south-facing group (adjust thresholds)")
    if not bool(np.any(mask_north)):
        raise AspectAnalysisError("empty north-facing group (adjust thresholds)")

    return AspectGroups(
        mask_south=mask_south,
        mask_north=mask_north,
        areas=areas,
        aq_depth_ele=aq_depth_ele,
        z_surf_ele=z_surf_ele,
    )


def _aw_mean(values: np.ndarray, *, areas: np.ndarray, mask: np.ndarray) -> float:
    a = areas[mask]
    return float((values[mask] * a).sum() / a.sum())


def mean_series_by_group(dat_path: Path, *, groups: AspectGroups) -> tuple[float, float]:
    r = DatReader(dat_path)
    if r.meta.num_var != groups.areas.size:
        raise AspectAnalysisError(
            f"shape mismatch: {dat_path} cols={r.meta.num_var} vs mesh elements={groups.areas.size}"
        )

    south = 0.0
    north = 0.0
    nrec = r.meta.num_records
    for _, vals in r.iter_records():
        v = np.asarray(vals, dtype=float)
        south += _aw_mean(v, areas=groups.areas, mask=groups.mask_south)
        north += _aw_mean(v, areas=groups.areas, mask=groups.mask_north)

    return south / nrec, north / nrec


def mean_unsat_ratio_by_group(out_dir: Path, *, groups: AspectGroups) -> tuple[float, float]:
    r_uns = DatReader(out_dir / "ccw.eleyunsat.dat")
    r_gw = DatReader(out_dir / "ccw.eleygw.dat")
    if r_uns.meta.num_records != r_gw.meta.num_records:
        raise AspectAnalysisError("Unsat/GW record mismatch")

    south = 0.0
    north = 0.0
    nrec = r_uns.meta.num_records

    for (_, v_uns), (_, v_gw) in zip(r_uns.iter_records(), r_gw.iter_records()):
        uns = np.asarray(v_uns, dtype=float)
        gw = np.asarray(v_gw, dtype=float)
        deficit = groups.aq_depth_ele - gw
        ratio = np.where(deficit > 1e-12, uns / deficit, np.nan)

        def wmean(mask: np.ndarray) -> float:
            m = mask & np.isfinite(ratio)
            a = groups.areas[m]
            return float((ratio[m] * a).sum() / a.sum())

        south += wmean(groups.mask_south)
        north += wmean(groups.mask_north)

    return south / nrec, north / nrec


def seasonal_delta_stats(out_dir: Path, *, dat_name: str, groups: AspectGroups) -> dict[str, dict[str, float]]:
    base_date = dt.date(2000, 1, 1)
    path = out_dir / dat_name
    r = DatReader(path)

    ts = []
    delta = []
    for t, vals in r.iter_records():
        v = np.asarray(vals, dtype=float)
        ts.append(float(t))
        delta.append(
            _aw_mean(v, areas=groups.areas, mask=groups.mask_south)
            - _aw_mean(v, areas=groups.areas, mask=groups.mask_north)
        )
    t_min = np.asarray(ts)
    delta = np.asarray(delta)

    days = (t_min / 1440.0).astype(int)
    months = np.asarray([(base_date + dt.timedelta(days=int(d))).month for d in days], dtype=int)
    summer = (months >= 4) & (months <= 9)
    winter = ~summer

    def stats(x: np.ndarray) -> dict[str, float]:
        return {
            "mean": float(np.mean(x)),
            "p50": float(np.quantile(x, 0.5)),
            "p95": float(np.quantile(x, 0.95)),
        }

    return {"all": stats(delta), "summer": stats(delta[summer]), "winter": stats(delta[winter])}


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Analyze aspect (south vs north slopes) response in SHUD outputs (ccw case)."
    )
    ap.add_argument("--mesh", type=Path, default=Path("input/ccw/ccw.sp.mesh"))
    ap.add_argument("--base", type=Path, default=Path("output/ccw.base"))
    ap.add_argument("--tsr", type=Path, default=Path("output/ccw.tsr"))
    ap.add_argument("--slope-deg-min", type=float, default=5.0)
    ap.add_argument("--aspect-halfwidth-deg", type=float, default=45.0)
    args = ap.parse_args()

    groups = build_aspect_groups(
        mesh_path=args.mesh,
        slope_deg_min=args.slope_deg_min,
        aspect_halfwidth_deg=args.aspect_halfwidth_deg,
    )

    area_total = float(groups.areas.sum())
    area_south = float(groups.areas[groups.mask_south].sum())
    area_north = float(groups.areas[groups.mask_north].sum())
    z_south = float((groups.z_surf_ele[groups.mask_south] * groups.areas[groups.mask_south]).sum() / area_south)
    z_north = float((groups.z_surf_ele[groups.mask_north] * groups.areas[groups.mask_north]).sum() / area_north)

    print("Aspect groups")
    print(f"- slope>= {args.slope_deg_min:g} deg, halfwidth=±{args.aspect_halfwidth_deg:g} deg")
    print(f"- south elems={int(groups.mask_south.sum())}, area_frac={area_south/area_total:.6f}")
    print(f"- north elems={int(groups.mask_north.sum())}, area_frac={area_north/area_total:.6f}")
    print(f"- z_surf mean (area-weighted): south={z_south:.3f} m, north={z_north:.3f} m")

    for label, outdir in [("BASE", args.base), ("TSR", args.tsr)]:
        rn_t_s, rn_t_n = mean_series_by_group(outdir / "ccw.rn_t.dat", groups=groups)
        et_s, et_n = mean_series_by_group(outdir / "ccw.eleveta.dat", groups=groups)
        uns_s, uns_n = mean_series_by_group(outdir / "ccw.eleyunsat.dat", groups=groups)
        gw_s, gw_n = mean_series_by_group(outdir / "ccw.eleygw.dat", groups=groups)
        ur_s, ur_n = mean_unsat_ratio_by_group(outdir, groups=groups)

        print(f"\n{label} (area-weighted time mean)")
        print(f"- rn_t (W/m2): south={rn_t_s:.6g}, north={rn_t_n:.6g}, delta={rn_t_s-rn_t_n:.6g}")
        print(f"- ETa (m/day): south={et_s:.6g}, north={et_n:.6g}, delta={et_s-et_n:.6g}")
        print(f"- yUnsat (m):  south={uns_s:.6g}, north={uns_n:.6g}, delta={uns_s-uns_n:.6g}")
        print(f"- yGW (m):     south={gw_s:.6g}, north={gw_n:.6g}, delta={gw_s-gw_n:.6g}")
        print(f"- yUnsat/(AqDepth-yGW): south={ur_s:.6g}, north={ur_n:.6g}, delta={ur_s-ur_n:.6g}")

        if label == "TSR":
            rn_t_stats = seasonal_delta_stats(outdir, dat_name="ccw.rn_t.dat", groups=groups)
            et_stats = seasonal_delta_stats(outdir, dat_name="ccw.eleveta.dat", groups=groups)
            print("\nTSR seasonal delta (south-north)")
            for name, st in [("rn_t (W/m2)", rn_t_stats), ("ETa (m/day)", et_stats)]:
                print(f"- {name}: all_mean={st['all']['mean']:.6g}, summer_mean={st['summer']['mean']:.6g}, winter_mean={st['winter']['mean']:.6g}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
