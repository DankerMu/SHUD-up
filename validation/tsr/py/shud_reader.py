from __future__ import annotations

import math
import os
import re
import struct
import sys
from array import array
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Literal, Optional, Sequence, Union

HEADER_BYTES = 1024


class DatFormatError(RuntimeError):
    def __init__(self, message: str, *, path: Optional[Path] = None) -> None:
        prefix = f"{path}: " if path is not None else ""
        super().__init__(f"{prefix}{message}")
        self.path = path


@dataclass(frozen=True)
class DatMeta:
    path: Path
    header_text: str
    header_lines: tuple[str, ...]
    start_time: float
    num_var: int
    col_ids: tuple[float, ...]
    num_records: int
    data_offset: int
    record_size: int
    dt_min: Optional[float]
    endianness: Literal["<", ">"]
    time_col_index: int = 0
    time_col_name: str = "Time_min"
    variable_name: Optional[str] = None
    dims: Optional[str] = None
    terrain_radiation: Optional[bool] = None


def _clean_header_text(raw: bytes) -> str:
    text = raw.decode("utf-8", errors="ignore")
    return text.replace("\x00", "").strip()


def _header_lines(text: str) -> tuple[str, ...]:
    cleaned = text.replace("\x00", "").strip()
    if not cleaned:
        return ()
    return tuple(ln.rstrip() for ln in cleaned.splitlines() if ln.strip())


_TSR_RE = re.compile(r"^#\s*Terrain radiation \(TSR\):\s*(ON|OFF)\s*$")


def _parse_tsr(header_lines: Sequence[str]) -> Optional[bool]:
    for line in header_lines:
        m = _TSR_RE.match(line)
        if m:
            return m.group(1) == "ON"
    return None


def _infer_variable_name(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".dat"):
        name = name[:-4]
    # Typical naming: <case>.<var> (e.g., ccw.rn_factor)
    if "." in name:
        return name.split(".", 1)[1]
    return name


def _is_int_like(x: float, *, tol: float = 1e-9) -> bool:
    if not math.isfinite(x):
        return False
    return abs(x - round(x)) <= tol


def _validate_and_shape(
    *, path: Path, size: int, num_var: int
) -> tuple[int, int, int]:
    if num_var < 0:
        raise DatFormatError(f"invalid NumVar={num_var} (must be >= 0)", path=path)

    data_offset = HEADER_BYTES + 16 + 8 * num_var
    if size < data_offset:
        raise DatFormatError(
            f"file too small for header+metadata (size={size}, need>={data_offset})",
            path=path,
        )

    record_size = 8 * (1 + num_var)
    rem = size - data_offset
    if rem % record_size != 0:
        raise DatFormatError(
            "data section not aligned: "
            f"(size={size}, data_offset={data_offset}, record_size={record_size}, remainder={rem})",
            path=path,
        )
    num_records = rem // record_size
    return data_offset, record_size, num_records


def _read_exact(f, n: int, *, path: Path, what: str) -> bytes:
    b = f.read(n)
    if len(b) != n:
        raise DatFormatError(f"truncated {what}: expected {n} bytes, got {len(b)}", path=path)
    return b


def _read_meta_for_endianness(path: Path, endianness: Literal["<", ">"]) -> DatMeta:
    if not path.exists():
        raise DatFormatError("missing .dat file", path=path)
    size = path.stat().st_size
    if size <= 0:
        raise DatFormatError("empty .dat file", path=path)
    if size < HEADER_BYTES:
        raise DatFormatError(f"short header (<{HEADER_BYTES}B)", path=path)

    with path.open("rb") as f:
        header_raw = _read_exact(f, HEADER_BYTES, path=path, what="header")
        header_text = _clean_header_text(header_raw)
        lines = _header_lines(header_text)

        start_time = struct.unpack(endianness + "d", _read_exact(f, 8, path=path, what="StartTime"))[0]
        num_var_raw = struct.unpack(endianness + "d", _read_exact(f, 8, path=path, what="NumVar"))[0]
        if not math.isfinite(num_var_raw):
            raise DatFormatError(f"invalid NumVar={num_var_raw!r}", path=path)

        num_var_int = int(num_var_raw)
        if abs(num_var_raw - num_var_int) > 1e-9:
            raise DatFormatError(f"NumVar is not an integer: {num_var_raw!r}", path=path)

        # Validate shape against file size before reading potentially-large col_ids.
        data_offset, record_size, num_records = _validate_and_shape(path=path, size=size, num_var=num_var_int)

        col_ids_raw = _read_exact(f, 8 * num_var_int, path=path, what="column id array (icol)")
        col_ids = struct.unpack(endianness + f"{num_var_int}d", col_ids_raw) if num_var_int else ()

        dt_min: Optional[float] = None
        if num_records >= 2:
            f.seek(data_offset, os.SEEK_SET)
            rec1 = _read_exact(f, record_size, path=path, what="record[0]")
            rec2 = _read_exact(f, record_size, path=path, what="record[1]")
            t0 = struct.unpack_from(endianness + "d", rec1, 0)[0]
            t1 = struct.unpack_from(endianness + "d", rec2, 0)[0]
            dt_min = t1 - t0
            if not math.isfinite(dt_min) or dt_min <= 0:
                raise DatFormatError(
                    f"non-positive dt derived from first two records: dt_min={dt_min!r}",
                    path=path,
                )

    return DatMeta(
        path=path,
        header_text=header_text,
        header_lines=lines,
        start_time=start_time,
        num_var=num_var_int,
        col_ids=tuple(float(x) for x in col_ids),
        num_records=int(num_records),
        data_offset=int(data_offset),
        record_size=int(record_size),
        dt_min=dt_min,
        endianness=endianness,
        variable_name=_infer_variable_name(path),
        dims=f"time={num_records}, var={num_var_int}",
        terrain_radiation=_parse_tsr(lines),
    )


def _detect_endianness(path: Path) -> Literal["<", ">"]:
    # SHUD writes using native endianness (typically little-endian). We try
    # little-endian first, then fall back to big-endian.
    le_err: Optional[DatFormatError] = None
    try:
        _read_meta_for_endianness(path, "<")
        return "<"
    except DatFormatError as e:
        le_err = e

    try:
        _read_meta_for_endianness(path, ">")
        return ">"
    except DatFormatError:
        # Prefer the original (little-endian) error, which is usually the most
        # actionable (e.g., empty file, short header, misalignment).
        raise le_err


class DatReader:
    def __init__(self, path: Union[str, Path], *, endianness: Literal["auto", "<", ">"] = "auto") -> None:
        self.path = Path(path)
        self._endianness: Literal["auto", "<", ">"] = endianness
        self._meta: Optional[DatMeta] = None

    @property
    def meta(self) -> DatMeta:
        if self._meta is None:
            end = _detect_endianness(self.path) if self._endianness == "auto" else self._endianness
            self._meta = _read_meta_for_endianness(self.path, end)
        return self._meta

    def iter_records(self) -> Iterable[tuple[float, tuple[float, ...]]]:
        meta = self.meta
        if meta.num_records <= 0:
            return

        with meta.path.open("rb") as f:
            f.seek(meta.data_offset, os.SEEK_SET)
            for i in range(meta.num_records):
                chunk = f.read(meta.record_size)
                if len(chunk) != meta.record_size:
                    raise DatFormatError(f"truncated record[{i}]", path=meta.path)
                t = struct.unpack_from(meta.endianness + "d", chunk, 0)[0]
                if meta.num_var:
                    values = struct.unpack_from(meta.endianness + f"{meta.num_var}d", chunk, 8)
                else:
                    values = ()
                yield t, tuple(float(x) for x in values)

    def read_matrix(self) -> tuple[list[float], list[list[float]]]:
        meta = self.meta
        times: list[float] = []
        rows: list[list[float]] = []
        for t, values in self.iter_records():
            times.append(float(t))
            rows.append([float(x) for x in values])
        return times, rows

    def to_dataframe(
        self,
        *,
        time_mode: Literal["index", "column"] = "index",
        column_naming: Literal["auto", "icol", "sequential"] = "auto",
        column_prefix: str = "X",
    ) -> Any:
        meta = self.meta
        try:
            import pandas as pd  # type: ignore
        except Exception as e:  # pragma: no cover
            raise RuntimeError(
                "pandas is required for DatReader.to_dataframe(); install it via `pip install pandas`"
            ) from e

        # Fast path: read the full payload as doubles.
        total_doubles = meta.num_records * (1 + meta.num_var)
        payload_offset = meta.data_offset

        a = array("d")
        with meta.path.open("rb") as f:
            f.seek(payload_offset, os.SEEK_SET)
            try:
                a.fromfile(f, total_doubles)
            except EOFError as e:
                raise DatFormatError("truncated data section while reading payload", path=meta.path) from e

        if len(a) != total_doubles:
            raise DatFormatError(
                f"truncated data section: expected {total_doubles} doubles, got {len(a)}",
                path=meta.path,
            )

        file_little = meta.endianness == "<"
        host_little = sys.byteorder == "little"
        if file_little != host_little:
            a.byteswap()

        import numpy as np  # type: ignore

        data = np.frombuffer(a, dtype=np.float64).reshape(meta.num_records, 1 + meta.num_var)
        time_min = data[:, 0]
        values = data[:, 1:]

        columns = self._column_names(column_naming=column_naming, column_prefix=column_prefix)
        df = pd.DataFrame(values, columns=columns)

        if time_mode == "index":
            df.index = pd.Index(time_min, name=meta.time_col_name)
        elif time_mode == "column":
            df.insert(0, meta.time_col_name, time_min)
        else:  # pragma: no cover
            raise ValueError(f"unsupported time_mode={time_mode!r}")

        return df

    def _column_names(
        self,
        *,
        column_naming: Literal["auto", "icol", "sequential"],
        column_prefix: str,
    ) -> list[str]:
        meta = self.meta

        if meta.num_var <= 0:
            return []

        if column_naming == "sequential":
            ids = list(range(1, meta.num_var + 1))
            return [f"{column_prefix}{i}" for i in ids]

        col_ids = list(meta.col_ids)
        if column_naming == "auto":
            if all(_is_int_like(x) for x in col_ids):
                return [f"{column_prefix}{int(round(x))}" for x in col_ids]
            ids = list(range(1, meta.num_var + 1))
            return [f"{column_prefix}{i}" for i in ids]

        # column_naming == "icol"
        if all(_is_int_like(x) for x in col_ids):
            return [f"{column_prefix}{int(round(x))}" for x in col_ids]
        return [f"{column_prefix}{x:g}" for x in col_ids]


def read_dat(
    path: Union[str, Path],
    *,
    time_mode: Literal["index", "column"] = "index",
    endianness: Literal["auto", "<", ">"] = "auto",
    column_naming: Literal["auto", "icol", "sequential"] = "auto",
    column_prefix: str = "X",
) -> Any:
    return DatReader(path, endianness=endianness).to_dataframe(
        time_mode=time_mode, column_naming=column_naming, column_prefix=column_prefix
    )
