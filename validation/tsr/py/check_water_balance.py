#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from shud_reader import DatReader


class CheckError(RuntimeError):
    pass


def _scan_max_abs(path: Path) -> tuple[float, float, int]:
    r = DatReader(path)
    meta = r.meta
    if meta.num_var <= 0:
        raise CheckError(f"{path}: NumVar={meta.num_var} (expected > 0)")
    if meta.num_records <= 0:
        raise CheckError(f"{path}: NumRecords={meta.num_records} (expected > 0)")

    max_abs = 0.0
    max_t = float("nan")
    max_i = -1

    for t, values in r.iter_records():
        for i, x in enumerate(values):
            if not math.isfinite(x):
                raise CheckError(f"{path}: non-finite residual at t={t}, i={i+1}: {x!r}")
            ax = abs(x)
            if ax > max_abs:
                max_abs = ax
                max_t = t
                max_i = i + 1  # 1-based element index

    return max_abs, max_t, max_i


def _find_single(output_dir: Path, suffix: str) -> Path:
    matches = sorted(output_dir.glob(f"*.{suffix}.dat"))
    if not matches:
        raise CheckError(f"{output_dir}: missing *.{suffix}.dat (did you rebuild+rerun?)")
    if len(matches) > 1:
        raise CheckError(f"{output_dir}: multiple *.{suffix}.dat found; please keep only one: {matches}")
    return matches[0]


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Check per-element water-balance residual outputs (*.elewb*_resid.dat) are within tolerance."
    )
    ap.add_argument(
        "output_dirs",
        nargs="+",
        type=Path,
        help="One or more SHUD output directories (e.g., output/ccw.base output/ccw.tsr)",
    )
    ap.add_argument("--tol", type=float, default=1e-6, help="Absolute tolerance in meters (default: 1e-6)")
    args = ap.parse_args()

    tol = float(args.tol)
    if not math.isfinite(tol) or tol <= 0:
        raise SystemExit("--tol must be finite and > 0")

    failures: list[str] = []

    for outdir in args.output_dirs:
        if not outdir.is_dir():
            raise SystemExit(f"not a directory: {outdir}")

        wb3 = _find_single(outdir, "elewb3_resid")
        wbf = _find_single(outdir, "elewbfull_resid")

        for label, path in [("wb3", wb3), ("wbfull", wbf)]:
            mx, t, i = _scan_max_abs(path)
            if mx > tol:
                failures.append(f"{outdir} ({label}): max_abs={mx:.6g} m at t={t:.1f} min, ele={i}")
            else:
                print(f"OK: {outdir} ({label}): max_abs={mx:.6g} m (<= {tol:g})")

    if failures:
        print("\nFAIL:")
        for line in failures:
            print(f"- {line}")
        return 2

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CheckError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise SystemExit(2)

