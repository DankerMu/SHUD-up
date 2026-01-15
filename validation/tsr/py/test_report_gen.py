from __future__ import annotations

import contextlib
import io
import struct
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Optional

# Allow running as: python3 -m unittest validation/tsr/py/test_report_gen.py
_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

import compare_tsr  # noqa: E402
from generate_report import generate_report, main as report_main  # noqa: E402
from tsr_core import TimeContext, read_mesh_normals, solar_position, solar_update_bucket, terrain_factor  # noqa: E402


def _write_dat(
    path: Path,
    *,
    header_lines: list[str],
    start_time: float,
    col_ids: list[float],
    records: list[tuple[float, list[float]]],
    endianness: str = "<",
) -> None:
    header_text = "\n".join(header_lines).rstrip() + "\n"
    header_raw = header_text.encode("utf-8")
    if len(header_raw) > 1024:
        raise ValueError("header too long for test helper")
    header_raw = header_raw.ljust(1024, b"\x00")

    num_var = len(col_ids)
    with path.open("wb") as f:
        f.write(header_raw)
        f.write(struct.pack(endianness + "d", float(start_time)))
        f.write(struct.pack(endianness + "d", float(num_var)))
        for x in col_ids:
            f.write(struct.pack(endianness + "d", float(x)))
        for t, vals in records:
            if len(vals) != num_var:
                raise ValueError("record width mismatch")
            f.write(struct.pack(endianness + "d", float(t)))
            for v in vals:
                f.write(struct.pack(endianness + "d", float(v)))


def _write_mesh(path: Path) -> None:
    # 2 elements, 4 nodes.
    # - ele 1: horizontal (z=0) => normal (0,0,1)
    # - ele 2: sloped (z=y) => normal (0,-1,1)/sqrt(2)
    path.write_text(
        "\n".join(
            [
                "2\t8",
                "ID\tNode1\tNode2\tNode3\tNabr1\tNabr2\tNabr3\tZmax",
                "1\t1\t2\t3\t0\t0\t0\t0",
                "2\t1\t2\t4\t0\t0\t0\t0",
                "4\t5",
                "ID\tX\tY\tAqDepth\tElevation",
                "1\t0\t0\t30\t0",
                "2\t1\t0\t30\t0",
                "3\t0\t1\t30\t0",
                "4\t0\t1\t30\t1",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _write_para(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "SOLAR_LONLAT_MODE FIXED",
                "SOLAR_LON_DEG 0.0",
                "SOLAR_LAT_DEG 0.0",
                "SOLAR_UPDATE_INTERVAL 60",
                "RAD_FACTOR_CAP 5.0",
                "RAD_COSZ_MIN 0.05",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _write_forc(path: Path, *, forc_start: int) -> None:
    path.write_text(
        "\n".join(
            [
                f"1 {float(forc_start):.0f}",
                "/dev/null",
                "ID Lon Lat",
                "1 0.0 0.0",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _write_shud(path: Path, *, mesh: Path, para: Path, forc: Path) -> None:
    path.write_text(
        "\n".join(
            [
                f"MESH {mesh}",
                f"PARA {para}",
                f"FORC {forc}",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _build_case(
    base_dir: Path,
    *,
    case: str = "case",
    start_yyyymmdd: int = 20000101,
    with_anomalies: bool = False,
    with_mismatch: bool = False,
) -> Path:
    out_dir = base_dir / "out"
    out_dir.mkdir(parents=True, exist_ok=True)

    mesh = base_dir / f"{case}.sp.mesh"
    para = base_dir / f"{case}.cfg.para"
    forc = base_dir / f"{case}.tsd.forc"
    shud = out_dir / f"{case}.SHUD"
    _write_mesh(mesh)
    _write_para(para)
    _write_forc(forc, forc_start=start_yyyymmdd)
    _write_shud(shud, mesh=mesh, para=para, forc=forc)

    normals = read_mesh_normals(mesh)
    normals_selected = [normals[1], normals[2]]
    tc = TimeContext(start_yyyymmdd)
    solar_update_interval = 60
    lat = 0.0
    lon = 0.0
    cap = 5.0
    cosz_min = 0.05

    times = [0.0, 60.0, 720.0]
    factor_rows: list[list[float]] = []
    rn_h_rows: list[list[float]] = []
    rn_t_rows: list[list[float]] = []
    etp_rows: list[list[float]] = []
    for t in times:
        _, t_aligned = solar_update_bucket(t, solar_update_interval)
        sp = solar_position(t_aligned, lat, lon, tc, timezone_hours=0.0)
        row_factor = [terrain_factor(nx, ny, nz, sp, cap, cosz_min) for (nx, ny, nz) in normals_selected]
        rn_h = max(0.0, float(sp.cosZ)) * 100.0
        row_rn_h = [rn_h, rn_h]
        row_rn_t = [h * f for h, f in zip(row_rn_h, row_factor)]
        row_etp = [0.01 * v for v in row_rn_t]

        factor_rows.append(row_factor)
        rn_h_rows.append(row_rn_h)
        rn_t_rows.append(row_rn_t)
        etp_rows.append(row_etp)

    if with_anomalies:
        # 1) NaN/Inf: inject NaN in rn_factor at t=0, ele=1.
        factor_rows[0][0] = float("nan")
        rn_t_rows[0][0] = float("nan")
        # 2) Night!=0: inject non-zero at night for ele=2.
        factor_rows[0][1] = 0.1
        rn_t_rows[0][1] = rn_h_rows[0][1] * factor_rows[0][1]
        # 3) Horizontal!=1: inject deviation for horizontal ele=1 at midday.
        factor_rows[2][0] = 1.2
        rn_t_rows[2][0] = rn_h_rows[2][0] * factor_rows[2][0]

    if with_mismatch:
        # Nudge one value slightly so compare_tsr fails strict tol checks.
        factor_rows[2][1] += 1.0e-6
        rn_t_rows[2][1] = rn_h_rows[2][1] * factor_rows[2][1]

    start_time = float(start_yyyymmdd)
    col_ids = [1.0, 2.0]
    header = ["# SHUD output", "# Terrain radiation (TSR): ON"]
    _write_dat(
        out_dir / f"{case}.rn_factor.dat",
        header_lines=header,
        start_time=start_time,
        col_ids=col_ids,
        records=list(zip(times, factor_rows)),
    )
    _write_dat(
        out_dir / f"{case}.rn_h.dat",
        header_lines=header,
        start_time=start_time,
        col_ids=col_ids,
        records=list(zip(times, rn_h_rows)),
    )
    _write_dat(
        out_dir / f"{case}.rn_t.dat",
        header_lines=header,
        start_time=start_time,
        col_ids=col_ids,
        records=list(zip(times, rn_t_rows)),
    )
    _write_dat(
        out_dir / f"{case}.elevetp.dat",
        header_lines=header,
        start_time=start_time,
        col_ids=col_ids,
        records=list(zip(times, etp_rows)),
    )
    return out_dir


class TestReportGeneration(unittest.TestCase):
    def test_generate_report_happy_path(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            out_dir = _build_case(Path(td))
            artifacts = generate_report(out_dir, sample_elements=[1], tolerance=1e-12, max_anomalies=10)

            self.assertTrue(artifacts.markdown_path.exists())
            self.assertTrue(any(p.exists() for p in artifacts.plot_paths))
            self.assertEqual([m.name for m in artifacts.metrics], ["rn_factor", "rn_t"])
            self.assertTrue(all(m.max_abs <= 1e-12 for m in artifacts.metrics))
            self.assertTrue(all(v == 0 for v in artifacts.anomaly_counts.values()))

            text = artifacts.markdown_path.read_text(encoding="utf-8")
            self.assertIn("# TSR 验证报告", text)
            self.assertIn("## 关键统计", text)
            self.assertIn("## 典型 Element 时间序列", text)
            self.assertIn("plots/element_1.png", text)

    def test_report_main_success_writes_to_custom_dir(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            out_dir = _build_case(Path(td))
            report_dir = Path(td) / "custom_report"
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                rc = report_main([str(out_dir), "--output-dir", str(report_dir), "--tolerance", "1e-12"])
            self.assertEqual(rc, 0)
            self.assertTrue((report_dir / "report.md").exists())

    def test_generate_report_detects_anomalies(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            out_dir = _build_case(Path(td), with_anomalies=True)
            artifacts = generate_report(out_dir, sample_elements=[1, 2], tolerance=1e-6, max_anomalies=50)

            self.assertGreater(artifacts.anomaly_counts.get("Night!=0", 0), 0)
            self.assertGreater(artifacts.anomaly_counts.get("Horizontal!=1", 0), 0)
            self.assertGreater(artifacts.anomaly_counts.get("NaN/Inf:rn_factor", 0), 0)

            text = artifacts.markdown_path.read_text(encoding="utf-8")
            self.assertIn("## 异常列表", text)
            self.assertIn("Night!=0", text)
            self.assertIn("Horizontal!=1", text)
            self.assertIn("NaN/Inf:rn_factor", text)

    def test_report_main_bad_dir(self) -> None:
        buf = io.StringIO()
        with contextlib.redirect_stderr(buf):
            rc = report_main(["/path/does/not/exist"])
        self.assertEqual(rc, 2)
        self.assertIn("not a directory", buf.getvalue())


class TestCompareTsrCLI(unittest.TestCase):
    def test_compare_tsr_main_pass(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            out_dir = _build_case(Path(td))
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                rc = compare_tsr.main([str(out_dir), "--tol", "1e-12"])
            self.assertEqual(rc, 0)

    def test_compare_tsr_main_fail(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            out_dir = _build_case(Path(td), with_mismatch=True)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                rc = compare_tsr.main([str(out_dir), "--tol", "1e-12"])
            self.assertEqual(rc, 1)


def load_tests(
    loader: unittest.TestLoader, tests: unittest.TestSuite, pattern: Optional[str]
) -> unittest.TestSuite:
    # When invoked directly (pattern is None), include the existing test suites so that:
    # - coverage stays representative for the whole validation/tsr/py code
    # - `python3 -m unittest validation/tsr/py/test_report_gen.py` is sufficient
    #
    # Under discovery (`python3 -m unittest discover ...`), avoid duplicating tests.
    if pattern is not None:
        return tests

    suite = unittest.TestSuite()
    suite.addTests(tests)
    for mod_name in ("test_tsr_core", "test_shud_reader"):
        mod = __import__(mod_name)
        suite.addTests(loader.loadTestsFromModule(mod))
    return suite
