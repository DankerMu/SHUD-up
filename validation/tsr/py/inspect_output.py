#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from shud_reader import DatFormatError, DatReader


def _list_dat_files(out_dir: Path) -> list[Path]:
    return sorted([p for p in out_dir.iterdir() if p.is_file() and p.suffix.lower() == ".dat"])


def _fmt_bool(v):
    if v is None:
        return "n/a"
    return "ON" if v else "OFF"


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Inspect SHUD output .dat files (list metadata / export to CSV)."
    )
    parser.add_argument("output_dir", type=Path, help="SHUD output directory (e.g., output/ccw.base)")
    parser.add_argument(
        "--export",
        type=str,
        default=None,
        help="Export the given .dat filename (within output_dir) to CSV",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="CSV output path (default: <dat>.csv next to the .dat file)",
    )
    parser.add_argument(
        "--show-header",
        action="store_true",
        help="Print the 1024B ASCII header (cleaned) for each .dat file",
    )
    args = parser.parse_args(argv)

    out_dir: Path = args.output_dir
    if not out_dir.is_dir():
        print(f"ERROR: not a directory: {out_dir}", file=sys.stderr)
        return 2

    dat_files = _list_dat_files(out_dir)
    if dat_files:
        print(f"Found {len(dat_files)} .dat files in {out_dir}:\n")
        for p in dat_files:
            try:
                r = DatReader(p)
                m = r.meta
            except DatFormatError as e:
                print(f"- {p.name}: ERROR: {e}")
                continue

            dt_part = f"{m.dt_min:g} min" if m.dt_min is not None else "n/a"
            print(
                f"- {p.name}: rows={m.num_records}, cols={m.num_var}, dt={dt_part}, TSR={_fmt_bool(m.terrain_radiation)}"
            )
            if args.show_header and m.header_text:
                print("  header:")
                for ln in m.header_lines:
                    print(f"    {ln}")
    else:
        print(f"No .dat files found in {out_dir}")

    if args.export is None:
        return 0

    dat_path = out_dir / args.export
    if not dat_path.exists():
        print(f"ERROR: missing .dat file: {dat_path}", file=sys.stderr)
        return 2

    csv_path: Path = args.csv if args.csv is not None else dat_path.with_suffix(".csv")
    reader = DatReader(dat_path)
    df = reader.to_dataframe(time_mode="column")
    df.to_csv(csv_path, index=False)
    print(f"\nExported: {dat_path.name} -> {csv_path}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main(sys.argv[1:]))
