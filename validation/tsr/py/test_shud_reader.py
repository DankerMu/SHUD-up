from __future__ import annotations

import contextlib
import io
import struct
import tempfile
import unittest
from pathlib import Path

from shud_reader import DatFormatError, DatReader
from inspect_output import main as inspect_main


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


class TestDatReader(unittest.TestCase):
    def test_meta_parsing_and_shape(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "ccw.rn_factor.dat"
            _write_dat(
                p,
                header_lines=[
                    "# SHUD output",
                    "# Terrain radiation (TSR): ON",
                ],
                start_time=20000101.0,
                col_ids=[1.0, 2.0, 5.0],
                records=[
                    (0.0, [1.0, 2.0, 3.0]),
                    (60.0, [4.0, 5.0, 6.0]),
                    (120.0, [7.0, 8.0, 9.0]),
                ],
            )

            r = DatReader(p)
            m = r.meta
            self.assertEqual(m.variable_name, "rn_factor")
            self.assertEqual(m.num_var, 3)
            self.assertEqual(m.num_records, 3)
            self.assertEqual(m.col_ids, (1.0, 2.0, 5.0))
            self.assertAlmostEqual(m.start_time, 20000101.0)
            self.assertAlmostEqual(m.dt_min or 0.0, 60.0)
            self.assertIs(m.terrain_radiation, True)

    def test_iter_records_and_read_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "case.var.dat"
            _write_dat(
                p,
                header_lines=["# SHUD output"],
                start_time=1.0,
                col_ids=[1.0, 2.0],
                records=[
                    (0.0, [10.0, 20.0]),
                    (30.0, [11.0, 21.0]),
                ],
            )
            r = DatReader(p)
            recs = list(r.iter_records())
            self.assertEqual(recs[0][0], 0.0)
            self.assertEqual(recs[0][1], (10.0, 20.0))
            self.assertEqual(recs[1][0], 30.0)
            self.assertEqual(recs[1][1], (11.0, 21.0))

            t, rows = r.read_matrix()
            self.assertEqual(t, [0.0, 30.0])
            self.assertEqual(rows, [[10.0, 20.0], [11.0, 21.0]])

    def test_to_dataframe_time_index_columns_auto(self) -> None:
        import pandas as pd  # noqa: F401

        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "ccw.rn_factor.dat"
            _write_dat(
                p,
                header_lines=["# SHUD output"],
                start_time=20000101.0,
                col_ids=[1.0, 2.0, 5.0],
                records=[
                    (0.0, [1.0, 2.0, 3.0]),
                    (60.0, [4.0, 5.0, 6.0]),
                ],
            )

            df = DatReader(p).to_dataframe(time_mode="index", column_naming="auto")
            self.assertEqual(df.shape, (2, 3))
            self.assertEqual(df.index.name, "Time_min")
            self.assertEqual(list(df.index.values), [0.0, 60.0])
            self.assertEqual(list(df.columns), ["X1", "X2", "X5"])
            self.assertEqual(df.loc[0.0, "X1"], 1.0)
            self.assertEqual(df.loc[60.0, "X5"], 6.0)

    def test_to_dataframe_time_column_sequential(self) -> None:
        import pandas as pd  # noqa: F401

        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "case.var.dat"
            _write_dat(
                p,
                header_lines=["# SHUD output"],
                start_time=1.0,
                col_ids=[10.0, 20.0, 30.0],
                records=[
                    (0.0, [1.0, 2.0, 3.0]),
                    (15.0, [4.0, 5.0, 6.0]),
                ],
            )
            df = DatReader(p).to_dataframe(time_mode="column", column_naming="sequential")
            self.assertEqual(df.shape, (2, 4))
            self.assertEqual(list(df.columns), ["Time_min", "X1", "X2", "X3"])
            self.assertEqual(list(df["Time_min"].values), [0.0, 15.0])
            self.assertEqual(df.loc[1, "X2"], 5.0)

    def test_big_endian_auto_detect(self) -> None:
        import pandas as pd  # noqa: F401

        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "case.var.dat"
            _write_dat(
                p,
                header_lines=["# SHUD output"],
                start_time=123.0,
                col_ids=[1.0, 2.0],
                records=[
                    (0.0, [1.0, 2.0]),
                    (1.0, [3.0, 4.0]),
                ],
                endianness=">",
            )
            r = DatReader(p, endianness="auto")
            self.assertEqual(r.meta.endianness, ">")
            df = r.to_dataframe(time_mode="index")
            self.assertEqual(list(df.index.values), [0.0, 1.0])
            self.assertEqual(df.loc[1.0, "X2"], 4.0)

    def test_empty_file_error(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "empty.dat"
            p.write_bytes(b"")
            with self.assertRaises(DatFormatError) as ctx:
                _ = DatReader(p).meta
            self.assertIn("empty .dat file", str(ctx.exception))

    def test_short_header_error(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "short.dat"
            p.write_bytes(b"# SHUD output\n")
            with self.assertRaises(DatFormatError) as ctx:
                _ = DatReader(p).meta
            self.assertIn("short header", str(ctx.exception))

    def test_numvar_not_integer_error(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "bad_numvar.dat"
            header = b"# SHUD output\n".ljust(1024, b"\x00")
            with p.open("wb") as f:
                f.write(header)
                f.write(struct.pack("<d", 0.0))  # StartTime
                f.write(struct.pack("<d", 3.5))  # NumVar (invalid)
            with self.assertRaises(DatFormatError) as ctx:
                _ = DatReader(p, endianness="<").meta
            self.assertIn("NumVar is not an integer", str(ctx.exception))

    def test_alignment_error(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "misaligned.dat"
            _write_dat(
                p,
                header_lines=["# SHUD output"],
                start_time=0.0,
                col_ids=[1.0],
                records=[(0.0, [1.0])],
            )
            with p.open("ab") as f:
                f.write(b"\x00")  # break alignment
            with self.assertRaises(DatFormatError) as ctx:
                _ = DatReader(p).meta
            self.assertIn("data section not aligned", str(ctx.exception))

    def test_non_integer_col_ids_column_naming(self) -> None:
        import pandas as pd  # noqa: F401

        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "case.var.dat"
            _write_dat(
                p,
                header_lines=["# SHUD output"],
                start_time=0.0,
                col_ids=[1.25, 2.5],
                records=[(0.0, [10.0, 20.0])],
            )
            df = DatReader(p).to_dataframe(time_mode="index", column_naming="icol")
            self.assertEqual(list(df.columns), ["X1.25", "X2.5"])


class TestInspectOutput(unittest.TestCase):
    def test_list_and_show_header(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td)
            _write_dat(
                out_dir / "case.var.dat",
                header_lines=["# SHUD output", "# Terrain radiation (TSR): OFF"],
                start_time=0.0,
                col_ids=[1.0],
                records=[(0.0, [1.0])],
            )
            # Include one unreadable .dat to cover error path.
            (out_dir / "empty.dat").write_bytes(b"")

            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = inspect_main([str(out_dir), "--show-header"])
            self.assertEqual(rc, 0)
            out = buf.getvalue()
            self.assertIn("Found 2 .dat files", out)
            self.assertIn("case.var.dat: rows=1, cols=1", out)
            self.assertIn("empty.dat: ERROR:", out)
            self.assertIn("# Terrain radiation (TSR): OFF", out)

    def test_export_to_csv(self) -> None:
        import pandas as pd  # noqa: F401

        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td)
            dat_name = "case.var.dat"
            _write_dat(
                out_dir / dat_name,
                header_lines=["# SHUD output"],
                start_time=0.0,
                col_ids=[1.0, 2.0],
                records=[(0.0, [10.0, 20.0])],
            )
            csv_path = out_dir / "export.csv"

            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = inspect_main([str(out_dir), "--export", dat_name, "--csv", str(csv_path)])
            self.assertEqual(rc, 0)
            self.assertTrue(csv_path.exists())
            self.assertIn("Exported:", buf.getvalue())
            self.assertEqual(csv_path.read_text(encoding="utf-8").splitlines()[0], "Time_min,X1,X2")

    def test_errors(self) -> None:
        # Not a directory.
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            rc = inspect_main(["/path/does/not/exist"])
        self.assertEqual(rc, 2)
        self.assertIn("not a directory", stderr.getvalue())

        # Missing exported file within directory.
        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td)
            stdout = io.StringIO()
            stderr = io.StringIO()
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                rc = inspect_main([str(out_dir), "--export", "missing.dat"])
            self.assertEqual(rc, 2)
            self.assertIn("missing .dat file", stderr.getvalue())


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
