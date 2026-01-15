#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from shud_reader import DatFormatError, DatReader
from tsr_core import (
    NA_VALUE,
    TimeContext,
    read_mesh_normals,
    solar_position,
    solar_update_bucket,
    terrain_factor,
)


DEFAULT_SOLAR_UPDATE_INTERVAL_MIN = 60
DEFAULT_RAD_FACTOR_CAP = 5.0
DEFAULT_RAD_COSZ_MIN = 0.05
DEFAULT_SOLAR_LONLAT_MODE = "FORCING_FIRST"


class CompareError(RuntimeError):
    pass


def _repo_root() -> Path:
    # validation/tsr/py/compare_tsr.py -> repo root is 3 parents up from py/.
    return Path(__file__).resolve().parents[3]


def _resolve_repo_path(repo_root: Path, raw: str) -> Path:
    p = Path(raw)
    if p.is_absolute():
        return p
    return (repo_root / p).resolve()


def _find_unique(out_dir: Path, suffix: str, *, case: Optional[str]) -> Path:
    if case is not None:
        p = out_dir / f"{case}.{suffix}.dat"
        if p.exists():
            return p
        raise CompareError(f"missing expected file: {p}")

    matches = sorted(out_dir.glob(f"*.{suffix}.dat"))
    if not matches:
        raise CompareError(f"missing '*.{suffix}.dat' in {out_dir}")
    if len(matches) != 1:
        names = ", ".join(p.name for p in matches[:10])
        extra = "..." if len(matches) > 10 else ""
        raise CompareError(
            f"multiple '*.{suffix}.dat' in {out_dir}; pass --case to disambiguate (found: {names}{extra})"
        )
    return matches[0]


def _parse_shud_file(path: Path) -> dict[str, str]:
    if not path.exists():
        raise CompareError(f"missing SHUD project file: {path}")
    m: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        key, value = parts[0].strip(), parts[1].strip()
        if key:
            m[key] = value
    return m


def _read_key_value_text(path: Path) -> dict[str, str]:
    m: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        key = parts[0]
        val = parts[1]
        m[key] = val
    return m


def _parse_solar_lonlat_mode(mode: str) -> str:
    s = mode.strip().upper()
    if s in ("FORCING_FIRST", "FORCING_MEAN", "FIXED"):
        return s
    # Allow numeric modes for parity with C++ parser.
    if s in ("0", "1", "2"):
        return {"0": "FORCING_FIRST", "1": "FORCING_MEAN", "2": "FIXED"}[s]
    return DEFAULT_SOLAR_LONLAT_MODE


def _parse_int(m: dict[str, str], key: str, default: int) -> int:
    raw = m.get(key)
    if raw is None:
        return default
    try:
        return int(float(raw))
    except Exception:
        return default


def _parse_float(m: dict[str, str], key: str, default: float) -> float:
    raw = m.get(key)
    if raw is None:
        return default
    try:
        return float(raw)
    except Exception:
        return default


def _read_forcing_list(path: Path) -> tuple[int, int, list[tuple[float, float]]]:
    if not path.exists():
        raise CompareError(f"missing forcing list file: {path}")

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise CompareError(f"forcing list file too short: {path}")

    head = lines[0].split()
    if len(head) < 2:
        raise CompareError(f"invalid forcing list header line: {lines[0]!r}")
    num_forc = int(head[0])
    forc_start = int(float(head[1]))

    # Skip lines[1] (path) and lines[2] (table header); parse station records.
    stations: list[tuple[float, float]] = []
    for raw in lines[3:]:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split()
        if len(cols) < 3:
            continue
        try:
            lon = float(cols[1])
            lat = float(cols[2])
        except Exception:
            continue
        stations.append((lon, lat))
        if len(stations) >= num_forc:
            break

    if len(stations) != num_forc:
        raise CompareError(
            f"forcing list station count mismatch in {path}: expected {num_forc}, got {len(stations)}"
        )
    return num_forc, forc_start, stations


def _select_solar_lonlat(
    *,
    mode: str,
    stations: list[tuple[float, float]],
    lon_fixed: float,
    lat_fixed: float,
) -> tuple[float, float]:
    # Mirrors MD_readin.cpp selection logic.
    if mode == "FIXED":
        if lon_fixed == NA_VALUE or lat_fixed == NA_VALUE:
            raise CompareError("SOLAR_LONLAT_MODE=FIXED but SOLAR_LON_DEG/SOLAR_LAT_DEG is missing")
        return lon_fixed, lat_fixed

    if mode == "FORCING_MEAN":
        sum_lon = 0.0
        sum_lat = 0.0
        n = 0
        for lo, la in stations:
            if lo == NA_VALUE or la == NA_VALUE:
                continue
            if lo < -180.0 or lo > 180.0 or la < -90.0 or la > 90.0:
                continue
            sum_lon += lo
            sum_lat += la
            n += 1
        if n > 0:
            return sum_lon / n, sum_lat / n
        raise CompareError("SOLAR_LONLAT_MODE=FORCING_MEAN but no valid Lon/Lat found in forcing list")

    # FORCING_FIRST (default)
    if not stations:
        raise CompareError("forcing list is empty; cannot select solar lon/lat")
    return stations[0]


@dataclass(frozen=True)
class ErrorStats:
    n: int
    rmse: float
    mae: float
    max_abs: float
    max_time_min: float
    max_col: int  # 0-based column index
    max_id: int  # element id (1-based)
    max_cpp: float
    max_py: float


def _stats(
    *,
    times: list[float],
    col_ids: list[int],
    cpp: list[list[float]],
    py: list[list[float]],
) -> ErrorStats:
    if len(times) != len(cpp) or len(times) != len(py):
        raise CompareError("row count mismatch in stats()")
    if not times:
        raise CompareError("no time records in stats()")

    n_total = 0
    sum_abs = 0.0
    sum_sq = 0.0
    max_abs = -1.0
    max_i = 0
    max_j = 0
    max_cpp = 0.0
    max_py = 0.0

    ncol = len(col_ids)
    for i, t in enumerate(times):
        row_cpp = cpp[i]
        row_py = py[i]
        if len(row_cpp) != ncol or len(row_py) != ncol:
            raise CompareError("column count mismatch in stats()")
        for j in range(ncol):
            d = row_py[j] - row_cpp[j]
            ad = abs(d)
            n_total += 1
            sum_abs += ad
            sum_sq += d * d
            if ad > max_abs:
                max_abs = ad
                max_i = i
                max_j = j
                max_cpp = row_cpp[j]
                max_py = row_py[j]

    rmse = math.sqrt(sum_sq / n_total) if n_total else float("nan")
    mae = sum_abs / n_total if n_total else float("nan")
    return ErrorStats(
        n=n_total,
        rmse=rmse,
        mae=mae,
        max_abs=max_abs,
        max_time_min=times[max_i],
        max_col=max_j,
        max_id=col_ids[max_j],
        max_cpp=max_cpp,
        max_py=max_py,
    )


def _read_matrix(path: Path) -> tuple[list[float], list[list[float]], list[int]]:
    r = DatReader(path)
    meta = r.meta
    times, rows = r.read_matrix()
    col_ids: list[int] = []
    for x in meta.col_ids:
        xi = int(round(x))
        if abs(x - xi) > 1e-9:
            raise CompareError(f"non-integer column id in {path}: {x!r}")
        col_ids.append(int(xi))
    return times, rows, col_ids


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Recompute TSR factor and rn_t in Python (C++-consistent) and compare pointwise to SHUD outputs."
    )
    parser.add_argument("output_dir", type=Path, help="SHUD output directory (e.g., output/ccw.tsr)")
    parser.add_argument("--case", type=str, default=None, help="Case prefix (e.g., ccw) if output_dir has multiple cases")
    parser.add_argument("--tol", type=float, default=1e-10, help="Max-abs tolerance for pass/fail (default: 1e-10)")
    parser.add_argument("--shud", type=Path, default=None, help="Path to <case>.SHUD (defaults to the one in output_dir)")
    parser.add_argument("--mesh", type=Path, default=None, help="Override mesh path (ccw.sp.mesh)")
    parser.add_argument("--para", type=Path, default=None, help="Override para cfg path (ccw.cfg.para)")
    parser.add_argument("--forc", type=Path, default=None, help="Override forcing list path (ccw.tsd.forc)")
    args = parser.parse_args(argv)

    out_dir: Path = args.output_dir
    if not out_dir.is_dir():
        print(f"ERROR: not a directory: {out_dir}", file=sys.stderr)
        return 2

    try:
        rn_factor_path = _find_unique(out_dir, "rn_factor", case=args.case)
        rn_h_path = _find_unique(out_dir, "rn_h", case=args.case)
        rn_t_path = _find_unique(out_dir, "rn_t", case=args.case)

        # Locate and parse .SHUD to find mesh/para/forc, unless explicitly provided.
        shud_path: Optional[Path] = args.shud
        if shud_path is None:
            shuds = sorted(out_dir.glob("*.SHUD"))
            if len(shuds) == 1:
                shud_path = shuds[0]
            elif len(shuds) == 0:
                shud_path = None
            else:
                raise CompareError("multiple *.SHUD files in output_dir; pass --shud or --case")

        repo_root = _repo_root()
        proj: dict[str, str] = {}
        if shud_path is not None:
            proj = _parse_shud_file(shud_path)

        mesh_path = args.mesh or (_resolve_repo_path(repo_root, proj["MESH"]) if "MESH" in proj else None)
        para_path = args.para or (_resolve_repo_path(repo_root, proj["PARA"]) if "PARA" in proj else None)
        forc_path = args.forc or (_resolve_repo_path(repo_root, proj["FORC"]) if "FORC" in proj else None)
        if mesh_path is None or para_path is None or forc_path is None:
            raise CompareError("missing mesh/para/forcing paths; pass --mesh/--para/--forc or provide a valid .SHUD")

        para_kv = _read_key_value_text(para_path)
        mode = _parse_solar_lonlat_mode(para_kv.get("SOLAR_LONLAT_MODE", DEFAULT_SOLAR_LONLAT_MODE))
        lon_fixed = _parse_float(para_kv, "SOLAR_LON_DEG", float(NA_VALUE))
        lat_fixed = _parse_float(para_kv, "SOLAR_LAT_DEG", float(NA_VALUE))

        solar_update_interval = _parse_int(
            para_kv, "SOLAR_UPDATE_INTERVAL", DEFAULT_SOLAR_UPDATE_INTERVAL_MIN
        )
        rad_factor_cap = _parse_float(para_kv, "RAD_FACTOR_CAP", DEFAULT_RAD_FACTOR_CAP)
        rad_cosz_min = _parse_float(para_kv, "RAD_COSZ_MIN", DEFAULT_RAD_COSZ_MIN)

        _, forc_start, stations = _read_forcing_list(forc_path)
        solar_lon_deg, solar_lat_deg = _select_solar_lonlat(
            mode=mode, stations=stations, lon_fixed=lon_fixed, lat_fixed=lat_fixed
        )

        # Read C++ outputs.
        times_f, factor_cpp, col_ids = _read_matrix(rn_factor_path)
        times_h, rn_h_cpp, col_ids_h = _read_matrix(rn_h_path)
        times_t, rn_t_cpp, col_ids_t = _read_matrix(rn_t_path)
        if times_f != times_h or times_f != times_t:
            raise CompareError("time index mismatch across rn_factor/rn_h/rn_t")
        if col_ids != col_ids_h or col_ids != col_ids_t:
            raise CompareError("column id mismatch across rn_factor/rn_h/rn_t")

        # Time base date is stored in .dat StartTime. Validate vs forcing list.
        start_time_dat = int(float(DatReader(rn_factor_path).meta.start_time))
        if start_time_dat != forc_start:
            raise CompareError(
                f"StartTime mismatch: rn_factor.dat has {start_time_dat}, forcing list has {forc_start}"
            )

        normals = read_mesh_normals(mesh_path)
        try:
            normals_selected = [normals[eid] for eid in col_ids]
        except KeyError as e:
            raise CompareError(f"mesh missing element id referenced by output: {e}") from e

        tc = TimeContext(start_time_dat)

        # Recompute factor and rn_t.
        factor_py: list[list[float]] = []
        rn_t_py: list[list[float]] = []
        for i, t in enumerate(times_f):
            _, t_aligned = solar_update_bucket(t, solar_update_interval)
            sp = solar_position(t_aligned, solar_lat_deg, solar_lon_deg, tc, timezone_hours=0.0)
            row_factor = [
                terrain_factor(nx, ny, nz, sp, rad_factor_cap, rad_cosz_min)
                for (nx, ny, nz) in normals_selected
            ]
            factor_py.append(row_factor)

            row_rn_t = [h * f for h, f in zip(rn_h_cpp[i], row_factor)]
            rn_t_py.append(row_rn_t)

        stats_factor = _stats(times=times_f, col_ids=col_ids, cpp=factor_cpp, py=factor_py)
        stats_rn_t = _stats(times=times_f, col_ids=col_ids, cpp=rn_t_cpp, py=rn_t_py)

        print("TSR pointwise comparison")
        print(f"- output_dir: {out_dir}")
        print(f"- mesh: {mesh_path}")
        print(f"- para: {para_path}")
        print(f"- forc: {forc_path}")
        print(f"- solar_lon/lat_mode: {mode}")
        print(f"- solar_lon_deg: {solar_lon_deg:.10f}")
        print(f"- solar_lat_deg: {solar_lat_deg:.10f}")
        print(f"- solar_update_interval_min: {solar_update_interval}")
        print(f"- rad_factor_cap: {rad_factor_cap:g}")
        print(f"- rad_cosz_min: {rad_cosz_min:g}")
        print()

        def _print_stats(name: str, s: ErrorStats) -> None:
            print(f"{name}: n={s.n} rmse={s.rmse:.3e} mae={s.mae:.3e} max_abs={s.max_abs:.3e}")
            print(
                f"  max@time_min={s.max_time_min:g} ele_id={s.max_id} "
                f"cpp={s.max_cpp:.17g} py={s.max_py:.17g} abs_err={s.max_abs:.3e}"
            )

        _print_stats("factor", stats_factor)
        _print_stats("rn_t", stats_rn_t)

        tol = float(args.tol)
        ok = stats_factor.max_abs <= tol and stats_rn_t.max_abs <= tol
        if not ok:
            print(f"\nRESULT: FAIL (tol={tol:g})", file=sys.stderr)
            return 1
        print(f"\nRESULT: PASS (tol={tol:g})")
        return 0

    except (CompareError, DatFormatError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main(sys.argv[1:]))

