from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path
import struct
import contextlib
import io

from tsr_core import (
    K_PI,
    NA_VALUE,
    SolarPosition,
    TimeContext,
    read_mesh_normals,
    solar_position,
    solar_update_bucket,
    terrain_factor,
)


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


class TestTimeContext(unittest.TestCase):
    def test_julian_day_basic(self) -> None:
        tc = TimeContext(20000101)
        self.assertEqual(tc.julian_day(0.0), 1)
        self.assertEqual(tc.julian_day(1440.0), 2)

    def test_julian_day_leap_year(self) -> None:
        # 2000 is a leap year.
        tc = TimeContext(20000228)
        self.assertEqual(tc.julian_day(0.0), 59)  # Feb 28
        self.assertEqual(tc.julian_day(1440.0), 60)  # Feb 29
        self.assertEqual(tc.julian_day(2 * 1440.0), 61)  # Mar 1

    def test_negative_minutes_truncation_semantics(self) -> None:
        tc = TimeContext(20000102)
        # -1 min should be previous day 23:59, matching C++ trunc-toward-zero + fixup.
        self.assertEqual(tc.format_iso(-1.0), "2000-01-01 23:59")
        self.assertEqual(tc.julian_day(-1.0), 1)

    def test_invalid_base_date(self) -> None:
        tc = TimeContext()
        tc.set_base_date(0)
        self.assertEqual(tc.julian_day(0.0), 0)
        self.assertEqual(tc.format_date(0.0), "0000-00-00")


class TestSolarBucket(unittest.TestCase):
    def test_bucket_alignment_epsilon(self) -> None:
        # Mirrors MD_ET.cpp bucket_eps logic: t=59.999999 + 1e-6 => 60.0 => bucket 1.
        b, t_aligned = solar_update_bucket(59.999999, 60)
        self.assertEqual(b, 1)
        self.assertEqual(t_aligned, 60.0)

    def test_bucket_nonpositive_interval_defaults(self) -> None:
        b, t_aligned = solar_update_bucket(0.0, 0)
        self.assertEqual(b, 0)
        self.assertEqual(t_aligned, 0.0)


class TestSolarPosition(unittest.TestCase):
    def test_solar_position_ranges(self) -> None:
        tc = TimeContext(20000101)
        sp = solar_position(0.0, 39.195, -122.71, tc, timezone_hours=0.0)
        self.assertTrue(math.isfinite(sp.cosZ))
        self.assertTrue(-1.0 <= sp.cosZ <= 1.0)
        self.assertTrue(0.0 <= sp.zenith <= K_PI)
        self.assertTrue(0.0 <= sp.azimuth < 2.0 * K_PI)
        self.assertTrue(math.isfinite(sp.declination))
        self.assertTrue(math.isfinite(sp.hourAngle))

    def test_timezone_fallback_when_missing(self) -> None:
        # Exercise the "timezone not provided" branch (timezone inferred from lon).
        tc = TimeContext(20000101)
        sp_inferred = solar_position(0.0, 0.0, 7.5, tc, timezone_hours=None)  # lon/15 = 0.5
        # std::round(0.5) => 1.0 (half away from zero), so inferred tz should match 1.0.
        sp_tz1 = solar_position(0.0, 0.0, 7.5, tc, timezone_hours=1.0)
        self.assertAlmostEqual(sp_inferred.hourAngle, sp_tz1.hourAngle, places=15)

        # But it should differ from tz=0.0.
        sp_tz0 = solar_position(0.0, 0.0, 7.5, tc, timezone_hours=0.0)
        self.assertNotAlmostEqual(sp_inferred.hourAngle, sp_tz0.hourAngle, places=12)

    def test_nonfinite_inputs_are_sanitized(self) -> None:
        tc = TimeContext(20000101)
        sp = solar_position(float("nan"), float("inf"), float("nan"), tc, timezone_hours=float("nan"))
        # Inputs should not propagate NaNs; all fields must remain finite.
        self.assertTrue(all(math.isfinite(x) for x in (sp.cosZ, sp.zenith, sp.azimuth, sp.declination, sp.hourAngle)))


class TestTerrainFactor(unittest.TestCase):
    def test_horizontal_surface_gives_factor_one(self) -> None:
        sp = SolarPosition(cosZ=0.8, zenith=0.0, azimuth=0.0, declination=0.0, hourAngle=0.0)
        f = terrain_factor(0.0, 0.0, 1.0, sp, cap=5.0, cosz_min=0.05)
        self.assertAlmostEqual(f, 1.0, places=15)

    def test_cosz_nonpositive_gives_zero(self) -> None:
        sp = SolarPosition(cosZ=0.0, zenith=0.0, azimuth=0.0, declination=0.0, hourAngle=0.0)
        self.assertEqual(terrain_factor(0.0, 0.0, 1.0, sp, cap=5.0, cosz_min=0.05), 0.0)

    def test_negative_cosi_gives_zero(self) -> None:
        # Sun direction: azimuth=pi/2 => sx=sinz, sy=0. Choose normal opposite sx to make cosi < 0.
        sp = SolarPosition(cosZ=0.5, zenith=0.0, azimuth=K_PI / 2.0, declination=0.0, hourAngle=0.0)
        self.assertEqual(terrain_factor(-1.0, 0.0, 0.0, sp, cap=5.0, cosz_min=0.05), 0.0)

    def test_cap_and_nonfinite_cap(self) -> None:
        # Normal aligned with sun direction => cosi=1, denom=cosz.
        sp = SolarPosition(cosZ=0.08, zenith=0.0, azimuth=0.0, declination=0.0, hourAngle=0.0)
        # For azimuth=0: sun vector = (sx=0, sy=sinz, sz=cosz).
        sinz = math.sqrt(max(0.0, 1.0 - sp.cosZ * sp.cosZ))
        nx, ny, nz = 0.0, sinz, sp.cosZ

        f1 = terrain_factor(nx, ny, nz, sp, cap=1.5, cosz_min=0.05)
        self.assertAlmostEqual(f1, 1.5, places=15)

        f2 = terrain_factor(nx, ny, nz, sp, cap=float("nan"), cosz_min=0.05)
        self.assertAlmostEqual(f2, 10.0, places=15)  # non-finite cap => 10.0

    def test_degenerate_normal_is_handled(self) -> None:
        sp = SolarPosition(cosZ=0.8, zenith=0.0, azimuth=0.0, declination=0.0, hourAngle=0.0)
        # Zero-length normal => fallback to (0,0,1) => factor=1.
        self.assertAlmostEqual(terrain_factor(0.0, 0.0, 0.0, sp, cap=5.0, cosz_min=0.05), 1.0, places=15)

    def test_invalid_inputs(self) -> None:
        sp = SolarPosition(cosZ=0.8, zenith=0.0, azimuth=0.0, declination=0.0, hourAngle=0.0)
        self.assertEqual(terrain_factor(float("nan"), 0.0, 1.0, sp, cap=5.0, cosz_min=0.05), 0.0)
        self.assertEqual(terrain_factor(0.0, 0.0, 1.0, sp, cap=5.0, cosz_min=float("inf")), 1.0)
        self.assertEqual(terrain_factor(0.0, 0.0, 1.0, sp, cap=-1.0, cosz_min=0.05), 0.0)  # cap->0 => 0


class TestMeshNormals(unittest.TestCase):
    def test_read_mesh_normals(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "case.sp.mesh"
            p.write_text(
                "\n".join(
                    [
                        "1\t8",
                        "ID\tNode1\tNode2\tNode3\tNabr1\tNabr2\tNabr3\tZmax",
                        "1\t1\t2\t3\t0\t0\t0\t0",
                        "3\t5",
                        "ID\tX\tY\tAqDepth\tElevation",
                        "1\t0\t0\t30\t0",
                        "2\t1\t0\t30\t0",
                        "3\t0\t1\t30\t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            normals = read_mesh_normals(p)
            nx, ny, nz = normals[1]
            # Expected normal for plane z=y: (0, -1, 1) normalized.
            self.assertAlmostEqual(nx, 0.0, places=15)
            self.assertAlmostEqual(ny, -1.0 / math.sqrt(2.0), places=15)
            self.assertAlmostEqual(nz, 1.0 / math.sqrt(2.0), places=15)

    def test_read_mesh_normals_degenerate_triangle(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "degenerate.sp.mesh"
            p.write_text(
                "\n".join(
                    [
                        "1 8",
                        "ID Node1 Node2 Node3 Nabr1 Nabr2 Nabr3 Zmax",
                        "1 1 2 3 0 0 0 0",
                        "3 5",
                        "ID X Y AqDepth Elevation",
                        "1 0 0 30 0",
                        "2 1 0 30 0",
                        "3 2 0 30 0",
                    ]
                ),
                encoding="utf-8",
            )
            normals = read_mesh_normals(p)
            self.assertEqual(normals[1], (0.0, 0.0, 1.0))


class TestCompareTSR(unittest.TestCase):
    def test_compare_script_synthetic(self) -> None:
        from compare_tsr import main as compare_main

        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td) / "out"
            out_dir.mkdir()

            mesh = Path(td) / "case.sp.mesh"
            mesh.write_text(
                "\n".join(
                    [
                        "1 8",
                        "ID Node1 Node2 Node3 Nabr1 Nabr2 Nabr3 Zmax",
                        "1 1 2 3 0 0 0 0",
                        "3 5",
                        "ID X Y AqDepth Elevation",
                        "1 0 0 30 0",
                        "2 1 0 30 0",
                        "3 0 1 30 0",
                    ]
                ),
                encoding="utf-8",
            )

            para = Path(td) / "case.cfg.para"
            para.write_text(
                "\n".join(
                    [
                        "SOLAR_LONLAT_MODE FIXED",
                        "SOLAR_LON_DEG 0",
                        "SOLAR_LAT_DEG 0",
                        "SOLAR_UPDATE_INTERVAL 60",
                        "RAD_FACTOR_CAP 5.0",
                        "RAD_COSZ_MIN 0.05",
                    ]
                ),
                encoding="utf-8",
            )

            forc = Path(td) / "case.tsd.forc"
            forc.write_text(
                "\n".join(
                    [
                        "1 20000101",
                        ".",
                        "ID Lon Lat X Y Z Filename",
                        "1 0 0 0 0 -9999 forcing.csv",
                    ]
                ),
                encoding="utf-8",
            )

            tc = TimeContext(20000101)
            t = 720.0
            sp = solar_position(t, 0.0, 0.0, tc, timezone_hours=0.0)
            f = terrain_factor(0.0, 0.0, 1.0, sp, cap=5.0, cosz_min=0.05)
            rn_h = 123.456

            header = ["# SHUD output", "# Terrain radiation (TSR): ON"]
            _write_dat(
                out_dir / "case.rn_factor.dat",
                header_lines=header,
                start_time=20000101.0,
                col_ids=[1.0],
                records=[(t, [f])],
            )
            _write_dat(
                out_dir / "case.rn_h.dat",
                header_lines=header,
                start_time=20000101.0,
                col_ids=[1.0],
                records=[(t, [rn_h])],
            )
            _write_dat(
                out_dir / "case.rn_t.dat",
                header_lines=header,
                start_time=20000101.0,
                col_ids=[1.0],
                records=[(t, [rn_h * f])],
            )

            buf = io.StringIO()
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                rc = compare_main(
                    [
                        str(out_dir),
                        "--case",
                        "case",
                        "--mesh",
                        str(mesh),
                        "--para",
                        str(para),
                        "--forc",
                        str(forc),
                        "--tol",
                        "1e-12",
                    ]
                )
            self.assertEqual(rc, 0)


class TestCompareTSRCoverage(unittest.TestCase):
    def test_helper_branches(self) -> None:
        import compare_tsr as ct

        self.assertEqual(ct._parse_solar_lonlat_mode("1"), "FORCING_MEAN")
        self.assertEqual(ct._parse_solar_lonlat_mode("2"), "FIXED")
        self.assertEqual(ct._parse_solar_lonlat_mode("bogus"), ct.DEFAULT_SOLAR_LONLAT_MODE)

        self.assertEqual(ct._parse_int({"X": "10"}, "X", 3), 10)
        self.assertEqual(ct._parse_int({"X": "bad"}, "X", 3), 3)
        self.assertEqual(ct._parse_int({}, "X", 3), 3)
        self.assertAlmostEqual(ct._parse_float({"X": "1.5"}, "X", 0.0), 1.5)
        self.assertAlmostEqual(ct._parse_float({"X": "bad"}, "X", 0.0), 0.0)
        self.assertAlmostEqual(ct._parse_float({}, "X", 0.0), 0.0)

        # _resolve_repo_path: absolute paths are returned as-is.
        repo = ct._repo_root()
        abs_p = Path("/tmp/abs_path")
        self.assertEqual(ct._resolve_repo_path(repo, str(abs_p)), abs_p)

    def test_find_unique_errors(self) -> None:
        import compare_tsr as ct

        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            # Case specified but file missing.
            with self.assertRaises(ct.CompareError):
                _ = ct._find_unique(d, "rn_factor", case="ccw")

            # No matches.
            with self.assertRaises(ct.CompareError):
                _ = ct._find_unique(d, "rn_factor", case=None)

            # Multiple matches.
            (d / "a.rn_factor.dat").write_bytes(b"")
            (d / "b.rn_factor.dat").write_bytes(b"")
            with self.assertRaises(ct.CompareError):
                _ = ct._find_unique(d, "rn_factor", case=None)

    def test_parse_shud_file_missing_and_present(self) -> None:
        import compare_tsr as ct

        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "missing.SHUD"
            with self.assertRaises(ct.CompareError):
                _ = ct._parse_shud_file(p)

            p.write_text(
                "\n".join(
                    [
                        "# comment",
                        "MESH ./mesh.sp.mesh",
                        "PARA ./case.cfg.para",
                        "FORC ./case.tsd.forc",
                        "BADLINE",
                    ]
                ),
                encoding="utf-8",
            )
            m = ct._parse_shud_file(p)
            self.assertEqual(m["MESH"], "./mesh.sp.mesh")
            self.assertIn("PARA", m)

    def test_forcing_list_and_solar_selection(self) -> None:
        import compare_tsr as ct

        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "case.tsd.forc"
            p.write_text(
                "\n".join(
                    [
                        "2 20000101",
                        ".",
                        "ID Lon Lat X Y Z Filename",
                        "1 0 0 0 0 -9999 forcing.csv",
                    ]
                ),
                encoding="utf-8",
            )
            with self.assertRaises(ct.CompareError):
                _ = ct._read_forcing_list(p)

        lon, lat = ct._select_solar_lonlat(
            mode="FORCING_MEAN",
            stations=[(0.0, 0.0), (10.0, 20.0)],
            lon_fixed=float(NA_VALUE),
            lat_fixed=float(NA_VALUE),
        )
        self.assertAlmostEqual(lon, 5.0)
        self.assertAlmostEqual(lat, 10.0)

        with self.assertRaises(ct.CompareError):
            _ = ct._select_solar_lonlat(
                mode="FIXED",
                stations=[(0.0, 0.0)],
                lon_fixed=float(NA_VALUE),
                lat_fixed=float(NA_VALUE),
            )

    def test_read_matrix_non_integer_cols(self) -> None:
        import compare_tsr as ct

        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "case.var.dat"
            _write_dat(
                p,
                header_lines=["# SHUD output"],
                start_time=20000101.0,
                col_ids=[1.25],
                records=[(0.0, [1.0])],
            )
            with self.assertRaises(ct.CompareError):
                _ = ct._read_matrix(p)

    def test_main_error_paths(self) -> None:
        from compare_tsr import main as compare_main

        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            rc = compare_main(["/path/does/not/exist"])
        self.assertEqual(rc, 2)

        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td) / "out"
            out_dir.mkdir()
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                rc = compare_main([str(out_dir)])
            self.assertEqual(rc, 2)

    def test_main_mismatched_times_and_starttime(self) -> None:
        from compare_tsr import main as compare_main

        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td) / "out"
            out_dir.mkdir()

            mesh = Path(td) / "case.sp.mesh"
            mesh.write_text(
                "\n".join(
                    [
                        "1 8",
                        "ID Node1 Node2 Node3 Nabr1 Nabr2 Nabr3 Zmax",
                        "1 1 2 3 0 0 0 0",
                        "3 5",
                        "ID X Y AqDepth Elevation",
                        "1 0 0 30 0",
                        "2 1 0 30 0",
                        "3 0 1 30 0",
                    ]
                ),
                encoding="utf-8",
            )

            para = Path(td) / "case.cfg.para"
            para.write_text(
                "\n".join(
                    [
                        "SOLAR_LONLAT_MODE FIXED",
                        "SOLAR_LON_DEG 0",
                        "SOLAR_LAT_DEG 0",
                    ]
                ),
                encoding="utf-8",
            )

            forc = Path(td) / "case.tsd.forc"
            forc.write_text(
                "\n".join(
                    [
                        "1 20000101",
                        ".",
                        "ID Lon Lat X Y Z Filename",
                        "1 0 0 0 0 -9999 forcing.csv",
                    ]
                ),
                encoding="utf-8",
            )

            header = ["# SHUD output"]
            _write_dat(
                out_dir / "case.rn_factor.dat",
                header_lines=header,
                start_time=20000101.0,
                col_ids=[1.0],
                records=[(0.0, [1.0])],
            )
            _write_dat(
                out_dir / "case.rn_h.dat",
                header_lines=header,
                start_time=20000101.0,
                col_ids=[1.0],
                records=[(60.0, [1.0])],  # mismatched time
            )
            _write_dat(
                out_dir / "case.rn_t.dat",
                header_lines=header,
                start_time=20000101.0,
                col_ids=[1.0],
                records=[(0.0, [1.0])],
            )

            buf = io.StringIO()
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                rc = compare_main(
                    [
                        str(out_dir),
                        "--case",
                        "case",
                        "--mesh",
                        str(mesh),
                        "--para",
                        str(para),
                        "--forc",
                        str(forc),
                    ]
                )
            self.assertEqual(rc, 2)

        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td) / "out"
            out_dir.mkdir()

            mesh = Path(td) / "case.sp.mesh"
            mesh.write_text(
                "\n".join(
                    [
                        "1 8",
                        "ID Node1 Node2 Node3 Nabr1 Nabr2 Nabr3 Zmax",
                        "1 1 2 3 0 0 0 0",
                        "3 5",
                        "ID X Y AqDepth Elevation",
                        "1 0 0 30 0",
                        "2 1 0 30 0",
                        "3 0 1 30 0",
                    ]
                ),
                encoding="utf-8",
            )

            para = Path(td) / "case.cfg.para"
            para.write_text(
                "\n".join(
                    [
                        "SOLAR_LONLAT_MODE FIXED",
                        "SOLAR_LON_DEG 0",
                        "SOLAR_LAT_DEG 0",
                    ]
                ),
                encoding="utf-8",
            )

            forc = Path(td) / "case.tsd.forc"
            forc.write_text(
                "\n".join(
                    [
                        "1 20000101",
                        ".",
                        "ID Lon Lat X Y Z Filename",
                        "1 0 0 0 0 -9999 forcing.csv",
                    ]
                ),
                encoding="utf-8",
            )

            header = ["# SHUD output"]
            _write_dat(
                out_dir / "case.rn_factor.dat",
                header_lines=header,
                start_time=19991231.0,  # mismatch
                col_ids=[1.0],
                records=[(0.0, [1.0])],
            )
            _write_dat(
                out_dir / "case.rn_h.dat",
                header_lines=header,
                start_time=19991231.0,
                col_ids=[1.0],
                records=[(0.0, [1.0])],
            )
            _write_dat(
                out_dir / "case.rn_t.dat",
                header_lines=header,
                start_time=19991231.0,
                col_ids=[1.0],
                records=[(0.0, [1.0])],
            )

            buf = io.StringIO()
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                rc = compare_main(
                    [
                        str(out_dir),
                        "--case",
                        "case",
                        "--mesh",
                        str(mesh),
                        "--para",
                        str(para),
                        "--forc",
                        str(forc),
                    ]
                )
            self.assertEqual(rc, 2)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
