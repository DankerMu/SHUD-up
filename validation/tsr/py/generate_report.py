#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as _dt
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import compare_tsr
from shud_reader import DatFormatError, DatReader

try:
    import numpy as np
except ImportError as e:  # pragma: no cover
    raise SystemExit("ERROR: missing dependency 'numpy' (try: python3 -m pip install numpy)") from e

try:
    import pandas as pd
except ImportError as e:  # pragma: no cover
    raise SystemExit("ERROR: missing dependency 'pandas' (try: python3 -m pip install pandas)") from e

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as e:  # pragma: no cover
    raise SystemExit("ERROR: missing dependency 'matplotlib' (try: python3 -m pip install matplotlib)") from e


DEFAULT_REPORT_DIRNAME = "report"
DEFAULT_TOLERANCE = 1e-10
DEFAULT_MAX_ANOMALIES = 200


class ReportError(RuntimeError):
    pass


@dataclass(frozen=True)
class MetricRow:
    name: str
    n: int
    rmse: float
    mae: float
    q50: float
    q90: float
    q99: float
    max_abs: float


@dataclass(frozen=True)
class ReportArtifacts:
    report_dir: Path
    markdown_path: Path
    plot_paths: list[Path]
    metrics: list[MetricRow]
    anomaly_counts: dict[str, int]


def _fmt_float(x: float) -> str:
    if not np.isfinite(x):
        return "nan"
    ax = abs(float(x))
    if ax != 0.0 and (ax < 1e-3 or ax >= 1e4):
        return f"{x:.3e}"
    return f"{x:.6g}"


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _downsample(x: np.ndarray, y: np.ndarray, *, max_points: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    if x.size <= max_points:
        return x, y
    idx = np.linspace(0, x.size - 1, num=max_points, dtype=int)
    return x[idx], y[idx]


def _metric_row(name: str, cpp: np.ndarray, py: np.ndarray) -> MetricRow:
    if cpp.shape != py.shape:
        raise ReportError(f"shape mismatch for {name}: cpp={cpp.shape} py={py.shape}")
    diff = py - cpp
    rmse = float(np.sqrt(np.mean(diff * diff)))
    mae = float(np.mean(np.abs(diff)))
    abs_err = np.abs(diff).ravel()
    qs = pd.Series(abs_err).quantile([0.5, 0.9, 0.99])
    return MetricRow(
        name=name,
        n=int(abs_err.size),
        rmse=rmse,
        mae=mae,
        q50=float(qs.loc[0.5]),
        q90=float(qs.loc[0.9]),
        q99=float(qs.loc[0.99]),
        max_abs=float(abs_err.max()) if abs_err.size else 0.0,
    )


def _write_markdown(path: Path, lines: Sequence[str]) -> None:
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def _find_first(out_dir: Path, suffix: str, *, case: Optional[str]) -> Optional[Path]:
    try:
        return compare_tsr._find_unique(out_dir, suffix, case=case)  # noqa: SLF001
    except compare_tsr.CompareError:
        return None


def _read_matrix(path: Path) -> tuple[list[float], np.ndarray, list[int]]:
    times, rows, col_ids = compare_tsr._read_matrix(path)  # noqa: SLF001
    return times, np.asarray(rows, dtype=float), col_ids


def _pick_sample_elements(col_ids: list[int], requested: Optional[list[int]]) -> list[int]:
    if requested:
        return list(dict.fromkeys(int(x) for x in requested))
    if not col_ids:
        return []
    picks = [col_ids[0], col_ids[len(col_ids) // 2], col_ids[-1]]
    return list(dict.fromkeys(int(x) for x in picks))


def _count_and_collect_nonfinite(
    *,
    name: str,
    times: list[float],
    col_ids: list[int],
    data: np.ndarray,
    max_items: int,
) -> tuple[int, list[dict[str, str]]]:
    if data.size == 0:
        return 0, []
    bad = ~np.isfinite(data)
    n_bad = int(bad.sum())
    if n_bad == 0:
        return 0, []
    idx = np.argwhere(bad)
    items: list[dict[str, str]] = []
    for i, j in idx[:max_items]:
        t = float(times[int(i)])
        eid = int(col_ids[int(j)])
        v = data[int(i), int(j)]
        items.append(
            {
                "type": "NaN/Inf",
                "var": name,
                "time_min": f"{t:g}",
                "element": str(eid),
                "value": repr(float(v)),
                "detail": "non-finite",
            }
        )
    return n_bad, items


def _count_and_collect_no_daylight_nonzero(
    *,
    times: list[float],
    col_ids: list[int],
    den: np.ndarray,
    factor_cpp: np.ndarray,
    rn_t_cpp: np.ndarray,
    tol: float,
    max_items: int,
) -> tuple[int, list[dict[str, str]]]:
    no_daylight = den <= 0.0
    if not bool(np.any(no_daylight)):
        return 0, []
    f_bad = np.abs(factor_cpp[no_daylight, :]) > float(tol)
    r_bad = np.abs(rn_t_cpp[no_daylight, :]) > float(tol)
    bad = f_bad | r_bad
    n_bad = int(bad.sum())
    if n_bad == 0:
        return 0, []
    time_idx = np.nonzero(no_daylight)[0]
    idx = np.argwhere(bad)
    items: list[dict[str, str]] = []
    for ii, jj in idx[:max_items]:
        i = int(time_idx[int(ii)])
        j = int(jj)
        t = float(times[i])
        eid = int(col_ids[j])
        items.append(
            {
                "type": "NoDaylight!=0",
                "var": "rn_factor/rn_t",
                "time_min": f"{t:g}",
                "element": str(eid),
                "value": f"factor={factor_cpp[i, j]:.6g}, rn_t={rn_t_cpp[i, j]:.6g}, den={den[i]:.6g}",
                "detail": f"no daylight (den<=0) but abs(value)>{tol:g}",
            }
        )
    return n_bad, items


def _count_and_collect_horizontal_factor(
    *,
    times: list[float],
    col_ids: list[int],
    min_sz: np.ndarray,
    cosz_threshold: float,
    normals: list[tuple[float, float, float]],
    factor_cpp: np.ndarray,
    tol: float,
    max_items: int,
) -> tuple[int, list[dict[str, str]]]:
    if not col_ids:
        return 0, []
    nx = np.asarray([n[0] for n in normals], dtype=float)
    ny = np.asarray([n[1] for n in normals], dtype=float)
    nz = np.asarray([n[2] for n in normals], dtype=float)
    norm2 = nx * nx + ny * ny + nz * nz
    # Only treat *nearly perfectly* horizontal surfaces as "horizontal".
    # (Avoid flagging gentle slopes where factor != 1 is expected.)
    horiz = (
        (np.abs(norm2 - 1.0) <= 1e-6)
        & (nz >= 0.999999)
        & (np.abs(nx) <= 1e-6)
        & (np.abs(ny) <= 1e-6)
    )
    if not bool(np.any(horiz)):
        return 0, []

    # Only check intervals where *all* included solar samples satisfy cosZ >= cosz_threshold.
    day = min_sz >= float(cosz_threshold)
    if not bool(np.any(day)):
        return 0, []

    cols = np.nonzero(horiz)[0]
    day_idx = np.nonzero(day)[0]
    sub = factor_cpp[np.ix_(day_idx, cols)]
    bad = np.abs(sub - 1.0) > float(tol)
    n_bad = int(bad.sum())
    if n_bad == 0:
        return 0, []

    idx = np.argwhere(bad)
    items: list[dict[str, str]] = []
    for ii, jj in idx[:max_items]:
        i = int(day_idx[int(ii)])
        col = int(cols[int(jj)])
        t = float(times[i])
        eid = int(col_ids[col])
        items.append(
            {
                "type": "Horizontal!=1",
                "var": "rn_factor",
                "time_min": f"{t:g}",
                "element": str(eid),
                "value": f"factor={factor_cpp[i, col]:.6g}, min_sz={min_sz[i]:.6g}, nz={nz[col]:.6g}",
                "detail": f"horizontal (nz≈1) but abs(factor-1)>{tol:g}",
            }
        )
    return n_bad, items


def _plot_element(
    *,
    out_png: Path,
    t_min: np.ndarray,
    element_id: int,
    factor_cpp: np.ndarray,
    factor_py: np.ndarray,
    rn_h_cpp: np.ndarray,
    rn_t_cpp: np.ndarray,
    rn_t_py: np.ndarray,
    etp_cpp: Optional[np.ndarray],
) -> None:
    t_hr_full = t_min / 60.0
    t_hr, factor_cpp_s = _downsample(t_hr_full, factor_cpp)
    _, factor_py_s = _downsample(t_hr_full, factor_py)
    _, rn_h_s = _downsample(t_hr_full, rn_h_cpp)
    _, rn_t_cpp_s = _downsample(t_hr_full, rn_t_cpp)
    _, rn_t_py_s = _downsample(t_hr_full, rn_t_py)

    nrows = 4 if etp_cpp is not None else 3
    fig, axes = plt.subplots(nrows=nrows, ncols=1, figsize=(12, 2.5 * nrows), sharex=True)
    if nrows == 1:
        axes = [axes]

    axes[0].plot(t_hr, factor_cpp_s, label="factor (C++)", linewidth=1.2)
    axes[0].plot(t_hr, factor_py_s, label="factor (py)", linewidth=1.0, linestyle="--")
    axes[0].set_ylabel("factor [-]")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="best")

    axes[1].plot(t_hr, rn_h_s, label="rn_h (C++)", linewidth=1.2)
    axes[1].set_ylabel("rn_h")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(loc="best")

    axes[2].plot(t_hr, rn_t_cpp_s, label="rn_t (C++)", linewidth=1.2)
    axes[2].plot(t_hr, rn_t_py_s, label="rn_t (py)", linewidth=1.0, linestyle="--")
    axes[2].set_ylabel("rn_t")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(loc="best")

    if etp_cpp is not None:
        _, etp_s = _downsample(t_hr_full, etp_cpp)
        axes[3].plot(t_hr, etp_s, label="ETP (C++)", linewidth=1.2)
        axes[3].set_ylabel("ETP")
        axes[3].grid(True, alpha=0.3)
        axes[3].legend(loc="best")

    axes[-1].set_xlabel("t (hour)")
    fig.suptitle(f"Element {element_id} time series", y=0.995)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def generate_report(
    output_dir: Path,
    *,
    report_dir: Optional[Path] = None,
    tolerance: float = DEFAULT_TOLERANCE,
    sample_elements: Optional[list[int]] = None,
    case: Optional[str] = None,
    max_anomalies: int = DEFAULT_MAX_ANOMALIES,
) -> ReportArtifacts:
    if not output_dir.is_dir():
        raise ReportError(f"not a directory: {output_dir}")

    report_dir = report_dir if report_dir is not None else (output_dir / DEFAULT_REPORT_DIRNAME)
    plots_dir = report_dir / "plots"
    _ensure_dir(plots_dir)
    md_path = report_dir / "report.md"

    rn_factor_path = compare_tsr._find_unique(output_dir, "rn_factor", case=case)  # noqa: SLF001
    rn_h_path = compare_tsr._find_unique(output_dir, "rn_h", case=case)  # noqa: SLF001
    rn_t_path = compare_tsr._find_unique(output_dir, "rn_t", case=case)  # noqa: SLF001
    elevetp_path = _find_first(output_dir, "elevetp", case=case)

    shud_path: Optional[Path] = None
    shuds = sorted(output_dir.glob("*.SHUD"))
    if len(shuds) == 1:
        shud_path = shuds[0]
    elif len(shuds) > 1:
        raise ReportError("multiple *.SHUD files in output_dir; pass --case and/or keep one")

    repo_root = compare_tsr._repo_root()  # noqa: SLF001
    proj: dict[str, str] = compare_tsr._parse_shud_file(shud_path) if shud_path is not None else {}  # noqa: SLF001
    mesh_path = compare_tsr._resolve_repo_path(repo_root, proj["MESH"]) if "MESH" in proj else None  # noqa: SLF001
    para_path = compare_tsr._resolve_repo_path(repo_root, proj["PARA"]) if "PARA" in proj else None  # noqa: SLF001
    forc_path = compare_tsr._resolve_repo_path(repo_root, proj["FORC"]) if "FORC" in proj else None  # noqa: SLF001
    if mesh_path is None or para_path is None or forc_path is None:
        raise ReportError("missing mesh/para/forcing paths; ensure output_dir has a valid *.SHUD file")

    para_kv = compare_tsr._read_key_value_text(para_path)  # noqa: SLF001
    mode = compare_tsr._parse_solar_lonlat_mode(  # noqa: SLF001
        para_kv.get("SOLAR_LONLAT_MODE", compare_tsr.DEFAULT_SOLAR_LONLAT_MODE)
    )
    lon_fixed = compare_tsr._parse_float(para_kv, "SOLAR_LON_DEG", float(compare_tsr.NA_VALUE))  # noqa: SLF001
    lat_fixed = compare_tsr._parse_float(para_kv, "SOLAR_LAT_DEG", float(compare_tsr.NA_VALUE))  # noqa: SLF001
    tsr_integration_step_min = compare_tsr._parse_int(  # noqa: SLF001
        para_kv, "TSR_INTEGRATION_STEP_MIN", compare_tsr.DEFAULT_TSR_INTEGRATION_STEP_MIN
    )
    if "TSR_INTEGRATION_STEP_MIN" not in para_kv and "SOLAR_UPDATE_INTERVAL" in para_kv:
        # Deprecated alias kept for backward compatibility.
        tsr_integration_step_min = compare_tsr._parse_int(  # noqa: SLF001
            para_kv, "SOLAR_UPDATE_INTERVAL", compare_tsr.DEFAULT_TSR_INTEGRATION_STEP_MIN
        )
    rad_factor_cap = compare_tsr._parse_float(  # noqa: SLF001
        para_kv, "RAD_FACTOR_CAP", compare_tsr.DEFAULT_RAD_FACTOR_CAP
    )
    rad_cosz_min = compare_tsr._parse_float(  # noqa: SLF001
        para_kv, "RAD_COSZ_MIN", compare_tsr.DEFAULT_RAD_COSZ_MIN
    )

    _, forc_start, stations = compare_tsr._read_forcing_list(forc_path)  # noqa: SLF001
    solar_lon_deg, solar_lat_deg = compare_tsr._select_solar_lonlat(  # noqa: SLF001
        mode=mode, stations=stations, lon_fixed=lon_fixed, lat_fixed=lat_fixed
    )

    times_f, factor_cpp, col_ids = _read_matrix(rn_factor_path)
    times_h, rn_h_cpp, col_ids_h = _read_matrix(rn_h_path)
    times_t, rn_t_cpp, col_ids_t = _read_matrix(rn_t_path)
    if times_f != times_h or times_f != times_t:
        raise ReportError("time index mismatch across rn_factor/rn_h/rn_t")
    if col_ids != col_ids_h or col_ids != col_ids_t:
        raise ReportError("column id mismatch across rn_factor/rn_h/rn_t")

    start_time_dat = int(float(DatReader(rn_factor_path).meta.start_time))
    if start_time_dat != forc_start:
        raise ReportError(f"StartTime mismatch: rn_factor.dat has {start_time_dat}, forcing list has {forc_start}")

    normals_map = compare_tsr.read_mesh_normals(mesh_path)
    try:
        normals_selected = [normals_map[eid] for eid in col_ids]
    except KeyError as e:
        raise ReportError(f"mesh missing element id referenced by output: {e}") from e

    tc = compare_tsr.TimeContext(start_time_dat)
    no_daylight_den: list[float] = []
    min_sz_in_interval: list[float] = []
    factor_py_rows: list[list[float]] = []
    rn_t_py_rows: list[list[float]] = []

    station_files = compare_tsr._read_forcing_station_files(forc_list=forc_path, repo_root=repo_root)  # noqa: SLF001
    if not station_files:
        raise ReportError(f"no forcing station files listed in {forc_path}")
    station_times_min = compare_tsr._read_station_times_min(station_files[0])  # noqa: SLF001

    cached_interval: Optional[tuple[float, float, int]] = None
    cached_samples: list[tuple[float, float, float, float]] = []
    cached_den = 0.0
    cached_min_sz = -1.0

    for i, t in enumerate(times_f):
        t0, t1 = compare_tsr._forcing_interval(station_times_min=station_times_min, t_min=t)  # noqa: SLF001
        key = (t0, t1, int(tsr_integration_step_min))
        if cached_interval != key:
            cached_interval = key
            cached_samples, cached_den = compare_tsr._solar_samples_for_interval(  # noqa: SLF001
                t0_min=t0,
                t1_min=t1,
                dt_int_min=int(tsr_integration_step_min),
                lat_deg=solar_lat_deg,
                lon_deg=solar_lon_deg,
                tc=tc,
            )
            cached_min_sz = min((float(sz) for _, _, sz, _ in cached_samples), default=-1.0)

        no_daylight_den.append(float(cached_den))
        min_sz_in_interval.append(float(cached_min_sz))

        row_factor = [
            compare_tsr._factor_from_samples(  # noqa: SLF001
                nx=nx,
                ny=ny,
                nz=nz,
                samples=cached_samples,
                den=cached_den,
                cap=rad_factor_cap,
                cosz_min=rad_cosz_min,
            )
            for (nx, ny, nz) in normals_selected
        ]
        factor_py_rows.append(row_factor)
        rn_t_py_rows.append([h * f for h, f in zip(rn_h_cpp[i, :].tolist(), row_factor)])

    factor_py = np.asarray(factor_py_rows, dtype=float)
    rn_t_py = np.asarray(rn_t_py_rows, dtype=float)
    den_np = np.asarray(no_daylight_den, dtype=float)
    min_sz_np = np.asarray(min_sz_in_interval, dtype=float)

    metrics = [
        _metric_row("rn_factor", factor_cpp, factor_py),
        _metric_row("rn_t", rn_t_cpp, rn_t_py),
    ]

    anomaly_items: list[dict[str, str]] = []
    anomaly_counts: dict[str, int] = {}

    for var_name, data in (("rn_factor", factor_cpp), ("rn_h", rn_h_cpp), ("rn_t", rn_t_cpp), ("rn_t_py", rn_t_py)):
        c, items = _count_and_collect_nonfinite(
            name=var_name, times=times_f, col_ids=col_ids, data=data, max_items=max_anomalies
        )
        anomaly_counts[f"NaN/Inf:{var_name}"] = c
        anomaly_items.extend(items)

    nodl_c, nodl_items = _count_and_collect_no_daylight_nonzero(
        times=times_f,
        col_ids=col_ids,
        den=den_np,
        factor_cpp=factor_cpp,
        rn_t_cpp=rn_t_cpp,
        tol=tolerance,
        max_items=max_anomalies,
    )
    anomaly_counts["NoDaylight!=0"] = nodl_c
    anomaly_items.extend(nodl_items)

    horiz_c, horiz_items = _count_and_collect_horizontal_factor(
        times=times_f,
        col_ids=col_ids,
        min_sz=min_sz_np,
        cosz_threshold=rad_cosz_min,
        normals=normals_selected,
        factor_cpp=factor_cpp,
        tol=tolerance,
        max_items=max_anomalies,
    )
    anomaly_counts["Horizontal!=1"] = horiz_c
    anomaly_items.extend(horiz_items)

    # Optional: ETP, only used for plotting.
    etp_cpp: Optional[np.ndarray] = None
    if elevetp_path is not None:
        times_e, etp_m, col_ids_e = _read_matrix(elevetp_path)
        if times_e == times_f and col_ids_e == col_ids:
            etp_cpp = etp_m

    col_index = {eid: i for i, eid in enumerate(col_ids)}
    chosen = [eid for eid in _pick_sample_elements(col_ids, sample_elements) if eid in col_index]
    if not chosen and col_ids:
        chosen = [col_ids[0]]

    plot_paths: list[Path] = []
    times_np = np.asarray(times_f, dtype=float)
    for eid in chosen:
        j = col_index[eid]
        out_png = plots_dir / f"element_{eid}.png"
        _plot_element(
            out_png=out_png,
            t_min=times_np,
            element_id=eid,
            factor_cpp=factor_cpp[:, j],
            factor_py=factor_py[:, j],
            rn_h_cpp=rn_h_cpp[:, j],
            rn_t_cpp=rn_t_cpp[:, j],
            rn_t_py=rn_t_py[:, j],
            etp_cpp=(etp_cpp[:, j] if etp_cpp is not None else None),
        )
        plot_paths.append(out_png)

    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    t0 = tc.format_iso(times_f[0]) if times_f else "n/a"
    t1 = tc.format_iso(times_f[-1]) if times_f else "n/a"

    lines: list[str] = []
    lines.append("# TSR 验证报告")
    lines.append("")
    lines.append(f"- 生成时间: {now}")
    lines.append(f"- 输出目录: `{output_dir}`")
    lines.append(f"- 时间范围: {t0} ~ {t1} (n={len(times_f)})")
    lines.append(f"- Element 数: {len(col_ids)}")
    lines.append("")
    lines.append("## 关键统计（C++ 输出 vs Python 复算）")
    lines.append("")
    lines.append("| 变量 | N | RMSE | MAE | Q50(|err|) | Q90(|err|) | Q99(|err|) | Max(|err|) |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for m in metrics:
        lines.append(
            "| "
            + " | ".join(
                [
                    m.name,
                    str(m.n),
                    _fmt_float(m.rmse),
                    _fmt_float(m.mae),
                    _fmt_float(m.q50),
                    _fmt_float(m.q90),
                    _fmt_float(m.q99),
                    _fmt_float(m.max_abs),
                ]
            )
            + " |"
        )
    lines.append("")
    lines.append("## 典型 Element 时间序列")
    lines.append("")
    lines.append(f"- 样本 elements: {', '.join(str(x) for x in chosen) if chosen else 'n/a'}")
    lines.append("")
    for p in plot_paths:
        lines.append(f"### {p.stem}")
        lines.append("")
        lines.append(f"![{p.name}](plots/{p.name})")
        lines.append("")

    lines.append("## 异常列表（自动检测）")
    lines.append("")
    lines.append(f"- 容差: `{tolerance:g}`")
    lines.append(f"- 夜间阈值: `cosZ < {rad_cosz_min:g}`")
    lines.append("")
    lines.append("| 类型 | 计数 |")
    lines.append("|---|---:|")
    for k in sorted(anomaly_counts.keys()):
        lines.append(f"| {k} | {anomaly_counts[k]} |")
    lines.append("")
    if anomaly_items:
        lines.append(f"仅展示前 {max_anomalies} 条：")
        lines.append("")
        lines.append("| type | var | time_min | element | value | detail |")
        lines.append("|---|---|---:|---:|---|---|")
        for it in anomaly_items[:max_anomalies]:
            lines.append(
                "| "
                + " | ".join(
                    [
                        it.get("type", ""),
                        it.get("var", ""),
                        it.get("time_min", ""),
                        it.get("element", ""),
                        it.get("value", "").replace("\n", " "),
                        it.get("detail", "").replace("\n", " "),
                    ]
                )
                + " |"
            )
        lines.append("")
    else:
        lines.append("未检测到异常。")
        lines.append("")

    lines.append("## 结论与下一步建议")
    lines.append("")
    ok = all(m.max_abs <= float(tolerance) for m in metrics)
    if ok and all(v == 0 for v in anomaly_counts.values()):
        lines.append(f"- 结论: `PASS`（关键变量误差 <= {tolerance:g}，且未发现物理约束异常）")
    else:
        lines.append("- 结论: `CHECK`（存在误差或异常项，建议按下方步骤排查）")
    lines.append("- 建议:")
    lines.append("  - 先运行 `python3 validation/tsr/py/compare_tsr.py <output_dir> --tol <tolerance>` 复核点位误差")
    lines.append("  - 若出现 NoDaylight!=0：检查 solar 几何/经纬度与 forcing 时间序列的区间边界是否正确")
    lines.append("  - 若出现 Horizontal!=1：检查 mesh 法向量计算/归一化与 `RAD_FACTOR_CAP`")
    lines.append("  - 若出现 NaN/Inf：优先查看对应输出文件与 `run.log`")
    lines.append("")
    lines.append("## 运行配置（摘录）")
    lines.append("")
    lines.append("```text")
    lines.append(f"mesh={mesh_path}")
    lines.append(f"para={para_path}")
    lines.append(f"forc={forc_path}")
    lines.append(f"SOLAR_LONLAT_MODE={mode}")
    lines.append(f"solar_lon_deg={solar_lon_deg:.10f}")
    lines.append(f"solar_lat_deg={solar_lat_deg:.10f}")
    lines.append(f"TSR_INTEGRATION_STEP_MIN={int(tsr_integration_step_min)}")
    lines.append(f"RAD_FACTOR_CAP={rad_factor_cap:g}")
    lines.append(f"RAD_COSZ_MIN={rad_cosz_min:g}")
    lines.append("```")

    _write_markdown(md_path, lines)

    return ReportArtifacts(
        report_dir=report_dir,
        markdown_path=md_path,
        plot_paths=plot_paths,
        metrics=metrics,
        anomaly_counts=anomaly_counts,
    )


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Generate TSR validation report (Markdown + plots).")
    parser.add_argument("output_dir", type=Path, help="SHUD output directory (e.g., output/ccw.tsr)")
    parser.add_argument(
        "--output-dir",
        dest="report_dir",
        type=Path,
        default=None,
        help=f"Report output directory (default: <output_dir>/{DEFAULT_REPORT_DIRNAME}/)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help=f"Tolerance for anomaly checks and pass/fail (default: {DEFAULT_TOLERANCE:g})",
    )
    parser.add_argument(
        "--sample-elements",
        type=int,
        nargs="*",
        default=None,
        help="Element IDs to plot (default: auto pick)",
    )
    parser.add_argument(
        "--case",
        type=str,
        default=None,
        help="Case prefix (e.g., ccw) if output_dir has multiple cases",
    )
    parser.add_argument(
        "--max-anomalies",
        type=int,
        default=DEFAULT_MAX_ANOMALIES,
        help=f"Max anomaly rows to include in markdown (default: {DEFAULT_MAX_ANOMALIES})",
    )
    args = parser.parse_args(argv)

    out_dir: Path = args.output_dir
    if not out_dir.is_dir():
        print(f"ERROR: not a directory: {out_dir}", file=sys.stderr)
        return 2

    try:
        artifacts = generate_report(
            out_dir,
            report_dir=args.report_dir,
            tolerance=float(args.tolerance),
            sample_elements=args.sample_elements,
            case=args.case,
            max_anomalies=int(args.max_anomalies),
        )
    except (ReportError, compare_tsr.CompareError, DatFormatError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    print("TSR validation report")
    print(f"- output_dir: {out_dir}")
    print(f"- report_dir: {artifacts.report_dir}")
    print(f"- report_md: {artifacts.markdown_path}")
    for m in artifacts.metrics:
        print(f"- {m.name}: rmse={m.rmse:.3e} mae={m.mae:.3e} max_abs={m.max_abs:.3e} (n={m.n})")
    print("- anomalies:")
    for k in sorted(artifacts.anomaly_counts.keys()):
        print(f"  - {k}: {artifacts.anomaly_counts[k]}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main(sys.argv[1:]))
