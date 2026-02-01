from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple, Union

# -----------------------------------------------------------------------------
# Constants (match C++ where relevant)
# -----------------------------------------------------------------------------

# SolarRadiation.cpp uses its own high-precision pi constant.
K_PI = 3.141592653589793238462643383279502884
K_TWO_PI = 2.0 * K_PI
K_DEG2RAD = K_PI / 180.0

# src/Model/Macros.hpp
ZERO = 1.0e-10
NA_VALUE = -9999


# -----------------------------------------------------------------------------
# TimeContext (Python reimplementation of src/classes/TimeContext.*)
# -----------------------------------------------------------------------------


def _is_finite(x: float) -> bool:
    return math.isfinite(x)


def _trunc_div_int(a: int, b: int) -> int:
    # C++ integer division truncates toward 0; Python // floors for negatives.
    if b == 0:
        raise ZeroDivisionError("division by zero")
    if a >= 0:
        return a // b
    return -((-a) // b)


@dataclass
class TimeContext:
    base_yyyymmdd: int = 0
    base_year: int = 0
    base_month: int = 0
    base_day: int = 0
    base_days_since_epoch: int = 0

    def __post_init__(self) -> None:
        if self.base_yyyymmdd:
            self.set_base_date(self.base_yyyymmdd)

    def set_base_date(self, yyyymmdd: int) -> None:
        y = 0
        m = 0
        d = 0
        ok = self._parse_yyyymmdd(yyyymmdd)
        if ok is None:
            self.base_yyyymmdd = 0
            self.base_year = 0
            self.base_month = 0
            self.base_day = 0
            self.base_days_since_epoch = 0
            return
        y, m, d = ok

        self.base_yyyymmdd = int(yyyymmdd)
        self.base_year = int(y)
        self.base_month = int(m)
        self.base_day = int(d)
        self.base_days_since_epoch = int(self._days_from_civil(self.base_year, self.base_month, self.base_day))

    def base_date(self) -> int:
        return int(self.base_yyyymmdd)

    def julian_day(self, t_min: float) -> int:
        y, m, d, _, _ = self._to_civil(t_min)
        if y == 0 or m == 0 or d == 0:
            return 0
        return int(self._day_of_year(y, m, d))

    def format_iso(self, t_min: float) -> str:
        y, m, d, hour, minute = self._to_civil(t_min)
        return f"{y:04d}-{m:02d}-{d:02d} {hour:02d}:{minute:02d}"

    def format_date(self, t_min: float) -> str:
        y, m, d, _, _ = self._to_civil(t_min)
        return f"{y:04d}-{m:02d}-{d:02d}"

    @staticmethod
    def _is_leap_year(year: int) -> bool:
        if (year % 4) != 0:
            return False
        if (year % 100) != 0:
            return True
        return (year % 400) == 0

    @classmethod
    def _days_in_month(cls, year: int, month: int) -> int:
        dim = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        if month < 1 or month > 12:
            return 0
        if month == 2 and cls._is_leap_year(year):
            return 29
        return int(dim[month - 1])

    @classmethod
    def _day_of_year(cls, year: int, month: int, day: int) -> int:
        cum = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
        if month < 1 or month > 12:
            return 0
        doy = int(cum[month - 1] + day)
        if month > 2 and cls._is_leap_year(year):
            doy += 1
        return int(doy)

    @staticmethod
    def _days_from_civil(y: int, m: int, d: int) -> int:
        # Matches src/classes/TimeContext.cpp exactly.
        y -= 1 if m <= 2 else 0
        era = (y if y >= 0 else y - 399) // 400
        yoe = y - era * 400
        doy = (153 * (m + (-3 if m > 2 else 9)) + 2) // 5 + d - 1
        doe = yoe * 365 + yoe // 4 - yoe // 100 + doy
        return int(era * 146097 + doe - 719468)

    @staticmethod
    def _civil_from_days(z: int) -> tuple[int, int, int]:
        # Matches src/classes/TimeContext.cpp exactly.
        z += 719468
        era = (z if z >= 0 else z - 146096) // 146097
        doe = z - era * 146097
        yoe = (doe - doe // 1460 + doe // 36524 - doe // 146096) // 365
        y = int(yoe + era * 400)
        doy = doe - (365 * yoe + yoe // 4 - yoe // 100)
        mp = (5 * doy + 2) // 153
        d = int(doy - (153 * mp + 2) // 5 + 1)
        m = int(mp + (3 if mp < 10 else -9))
        y += 1 if m <= 2 else 0
        return y, m, d

    @classmethod
    def _parse_yyyymmdd(cls, yyyymmdd: int) -> Optional[tuple[int, int, int]]:
        if yyyymmdd <= 0:
            return None
        y = int(yyyymmdd // 10000)
        md = int(yyyymmdd % 10000)
        m = int(md // 100)
        d = int(md % 100)
        if m < 1 or m > 12:
            return None
        dim = cls._days_in_month(y, m)
        if dim <= 0 or d < 1 or d > dim:
            return None
        return y, m, d

    def _to_civil(self, t_min: float) -> tuple[int, int, int, int, int]:
        # Matches src/classes/TimeContext.cpp exactly (incl. truncation rules).
        if self.base_yyyymmdd <= 0:
            return 0, 0, 0, 0, 0

        total_minutes = 0
        if _is_finite(t_min):
            total_minutes = int(t_min)  # trunc toward 0, like C++ (long long)t_min

        day_offset = _trunc_div_int(total_minutes, 1440)
        minute_of_day = total_minutes - day_offset * 1440  # C++ remainder semantics
        if minute_of_day < 0:
            minute_of_day += 1440
            day_offset -= 1

        days = int(self.base_days_since_epoch + day_offset)
        y, m, d = self._civil_from_days(days)
        hour = int(minute_of_day // 60)
        minute = int(minute_of_day % 60)
        return int(y), int(m), int(d), hour, minute


# -----------------------------------------------------------------------------
# Solar geometry + TSR factor (Python reimplementation of src/Equations/SolarRadiation.*)
# -----------------------------------------------------------------------------


def _clamp(x: float, lo: float, hi: float) -> float:
    if x < lo:
        return lo
    if x > hi:
        return hi
    return x


def _wrap_degrees_180(lon_deg: float) -> float:
    if not _is_finite(lon_deg):
        return 0.0
    x = math.fmod(lon_deg, 360.0)
    if x > 180.0:
        x -= 360.0
    elif x < -180.0:
        x += 360.0
    return x


def _clamp_latitude(lat_deg: float) -> float:
    if not _is_finite(lat_deg):
        return 0.0
    return _clamp(lat_deg, -90.0, 90.0)


def _minute_of_day(t_min: float) -> float:
    if not _is_finite(t_min):
        return 0.0
    m = math.fmod(t_min, 1440.0)
    if m < 0.0:
        m += 1440.0
    return m


def _wrap_minutes_1440(minutes: float) -> float:
    if not _is_finite(minutes):
        return 0.0
    m = math.fmod(minutes, 1440.0)
    if m < 0.0:
        m += 1440.0
    return m


def _wrap_radians_2pi(rad: float) -> float:
    if not _is_finite(rad):
        return 0.0
    x = math.fmod(rad, K_TWO_PI)
    if x < 0.0:
        x += K_TWO_PI
    return x


def _round_half_away_from_zero(x: float) -> float:
    # std::round: halfway cases away from zero.
    if not _is_finite(x):
        return 0.0
    ax = abs(x)
    y = math.floor(ax + 0.5)
    return y if x >= 0.0 else -y


def _approximate_timezone_hours(lon_deg: float) -> float:
    if not _is_finite(lon_deg):
        return 0.0
    return _round_half_away_from_zero(lon_deg / 15.0)


def _safe_acos(x: float) -> float:
    return math.acos(_clamp(x, -1.0, 1.0))


@dataclass(frozen=True)
class SolarPosition:
    cosZ: float
    zenith: float
    azimuth: float
    declination: float
    hourAngle: float


def solar_position(
    t_min: float,
    lat_deg: float,
    lon_deg: float,
    tc: TimeContext,
    timezone_hours: Optional[float] = None,
) -> SolarPosition:
    # Mirrors solarPositionImpl() in SolarRadiation.cpp.
    cosz = 0.0
    zenith = K_PI / 2.0
    azimuth = 0.0
    declination = 0.0
    hour_angle = 0.0

    lat = _clamp_latitude(lat_deg)
    lon = _wrap_degrees_180(lon_deg)

    tz = float(timezone_hours) if timezone_hours is not None else 0.0
    timezone_provided = timezone_hours is not None
    if (not timezone_provided) or (not _is_finite(tz)):
        tz = _approximate_timezone_hours(lon)

    doy = int(tc.julian_day(t_min))
    if doy < 1 or doy > 366:
        doy = 1

    mod_min = _minute_of_day(t_min)
    hour = mod_min / 60.0

    gamma = (K_TWO_PI / 365.0) * (float(doy - 1) + (hour - 12.0) / 24.0)
    sin_g = math.sin(gamma)
    cos_g = math.cos(gamma)
    sin_2g = math.sin(2.0 * gamma)
    cos_2g = math.cos(2.0 * gamma)
    sin_3g = math.sin(3.0 * gamma)
    cos_3g = math.cos(3.0 * gamma)

    eq_time_min = 229.18 * (
        0.000075
        + 0.001868 * cos_g
        - 0.032077 * sin_g
        - 0.014615 * cos_2g
        - 0.040849 * sin_2g
    )

    decl = (
        0.006918
        - 0.399912 * cos_g
        + 0.070257 * sin_g
        - 0.006758 * cos_2g
        + 0.000907 * sin_2g
        - 0.002697 * cos_3g
        + 0.00148 * sin_3g
    )
    declination = decl

    time_offset_min = eq_time_min + 4.0 * lon - 60.0 * tz
    true_solar_time_min = _wrap_minutes_1440(mod_min + time_offset_min)

    ha_deg = true_solar_time_min / 4.0 - 180.0
    ha = ha_deg * K_DEG2RAD
    hour_angle = ha

    lat_rad = lat * K_DEG2RAD
    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)

    sin_decl = math.sin(decl)
    cos_decl = math.cos(decl)
    sin_ha = math.sin(ha)
    cos_ha = math.cos(ha)

    cosz_raw = sin_lat * sin_decl + cos_lat * cos_decl * cos_ha
    cosz = _clamp(cosz_raw, -1.0, 1.0)
    zenith = _safe_acos(cosz)

    east = -cos_decl * sin_ha
    north = cos_lat * sin_decl - sin_lat * cos_decl * cos_ha
    azimuth = _wrap_radians_2pi(math.atan2(east, north))

    sp = SolarPosition(
        cosZ=cosz,
        zenith=zenith,
        azimuth=azimuth,
        declination=declination,
        hourAngle=hour_angle,
    )
    if not (
        _is_finite(sp.cosZ)
        and _is_finite(sp.zenith)
        and _is_finite(sp.azimuth)
        and _is_finite(sp.declination)
        and _is_finite(sp.hourAngle)
    ):
        return SolarPosition(
            cosZ=0.0,
            zenith=K_PI / 2.0,
            azimuth=0.0,
            declination=0.0,
            hourAngle=0.0,
        )
    return sp


def terrain_factor(
    nx: float,
    ny: float,
    nz: float,
    sp: SolarPosition,
    cap: float,
    cosz_min: float,
) -> float:
    # Mirrors terrainFactor() in SolarRadiation.cpp.
    if not (
        _is_finite(nx)
        and _is_finite(ny)
        and _is_finite(nz)
        and _is_finite(sp.cosZ)
        and _is_finite(sp.azimuth)
    ):
        return 0.0

    if not _is_finite(cap):
        cap = 10.0
    elif cap < 0.0:
        cap = 0.0

    if (not _is_finite(cosz_min)) or cosz_min < 0.0:
        cosz_min = 0.0
    cosz_min = _clamp(cosz_min, 0.0, 1.0)

    cosz = sp.cosZ
    if not (cosz > 0.0):
        return 0.0

    nlen = math.sqrt(nx * nx + ny * ny + nz * nz)
    if not (nlen > 0.0) or (not _is_finite(nlen)):
        nx = 0.0
        ny = 0.0
        nz = 1.0
    else:
        nx /= nlen
        ny /= nlen
        nz /= nlen
        if nz < 0.0:
            nx = -nx
            ny = -ny
            nz = -nz

    cosz_clamped = _clamp(cosz, -1.0, 1.0)
    sinz = math.sqrt(max(0.0, 1.0 - cosz_clamped * cosz_clamped))

    sin_az = math.sin(sp.azimuth)
    cos_az = math.cos(sp.azimuth)
    sx = sinz * sin_az
    sy = sinz * cos_az
    sz = cosz_clamped

    cosi = nx * sx + ny * sy + nz * sz
    if not (cosi > 0.0) or (not _is_finite(cosi)):
        return 0.0

    denom = cosz_clamped
    if denom < cosz_min:
        denom = cosz_min
    if not (denom > 0.0) or (not _is_finite(denom)):
        return 0.0

    factor = cosi / denom
    if (not _is_finite(factor)) or factor <= 0.0:
        return 0.0

    if factor > cap:
        factor = cap
    if (not _is_finite(factor)) or factor <= 0.0:
        return 0.0
    return factor


def solar_update_bucket(t_min: float, interval_min: int) -> tuple[int, float]:
    # Matches the bucket alignment in src/ModelData/MD_ET.cpp.
    if interval_min <= 0:
        interval_min = 60
    if not _is_finite(t_min):
        return 0, 0.0
    bucket_eps = 1.0e-6
    bucket = int(math.floor((t_min + bucket_eps) / float(interval_min)))
    t_aligned = float(bucket) * float(interval_min)
    return bucket, t_aligned


def forcing_interval_factor(
    nx: float,
    ny: float,
    nz: float,
    *,
    t0_min: float,
    t1_min: float,
    lat_deg: float,
    lon_deg: float,
    tc: TimeContext,
    cap: float,
    cosz_min: float,
    dt_int_min: int = 60,
    timezone_hours: float = 0.0,
) -> float:
    """
    Effective TSR factor over forcing interval [t0, t1), assuming forcing shortwave is an interval-mean flux.

    This mirrors the intended C++ behavior for TSR_FACTOR_MODE=FORCING_INTERVAL:
      F_eff ≈ (∫ w(t) f(t) dt) / (∫ w(t) dt), where w(t)=max(cosZ,0).
    """
    if not (_is_finite(t0_min) and _is_finite(t1_min)):
        return 0.0
    if not (t1_min > t0_min):
        return 0.0
    if dt_int_min <= 0:
        dt_int_min = 60

    dt_forc = float(t1_min - t0_min)
    dt_int = min(float(dt_int_min), dt_forc)
    n = int(math.ceil(dt_forc / dt_int)) if dt_int > 0.0 else 1
    if n < 1:
        n = 1
    dt_seg = dt_forc / float(n)

    num = 0.0
    den = 0.0
    for k in range(n):
        tk = float(t0_min + (k + 0.5) * dt_seg)
        sp = solar_position(tk, lat_deg, lon_deg, tc, timezone_hours=timezone_hours)
        w = max(0.0, float(sp.cosZ))
        if not (w > 0.0):
            continue
        fk = terrain_factor(nx, ny, nz, sp, cap=cap, cosz_min=cosz_min)
        num += w * fk * dt_seg
        den += w * dt_seg

    if not (den > 0.0):
        return 0.0
    feff = num / den
    if (not _is_finite(feff)) or (not (feff > 0.0)):
        return 0.0
    return float(feff)


# -----------------------------------------------------------------------------
# Mesh reader: compute terrain normal vectors like src/classes/Element.cpp
# -----------------------------------------------------------------------------


class MeshFormatError(RuntimeError):
    pass


def read_mesh_normals(mesh_path: Union[str, Path]) -> dict[int, tuple[float, float, float]]:
    path = Path(mesh_path)
    if not path.exists():
        raise MeshFormatError(f"missing mesh file: {path}")

    with path.open("r", encoding="utf-8", errors="replace") as f:
        first = f.readline()
        if not first:
            raise MeshFormatError(f"empty mesh file: {path}")
        parts = first.split()
        if len(parts) < 1:
            raise MeshFormatError(f"invalid mesh header line: {first!r}")
        num_ele = int(parts[0])

        # Element table header.
        hdr = f.readline()
        if not hdr:
            raise MeshFormatError(f"missing element header row in: {path}")

        elements: dict[int, tuple[int, int, int]] = {}
        for _ in range(num_ele):
            line = f.readline()
            if not line:
                raise MeshFormatError(f"unexpected EOF in element table: {path}")
            cols = line.split()
            if len(cols) < 4:
                raise MeshFormatError(f"invalid element row: {line!r}")
            eid = int(cols[0])
            n1 = int(cols[1])
            n2 = int(cols[2])
            n3 = int(cols[3])
            elements[eid] = (n1, n2, n3)

        # Node table header: "<NumNode> 5"
        node_hdr = f.readline()
        if not node_hdr:
            raise MeshFormatError(f"missing node table header in: {path}")
        node_hdr_parts = node_hdr.split()
        if len(node_hdr_parts) < 1:
            raise MeshFormatError(f"invalid node table header line: {node_hdr!r}")
        num_node = int(node_hdr_parts[0])

        # Node table column header.
        node_cols = f.readline()
        if not node_cols:
            raise MeshFormatError(f"missing node header row in: {path}")

        # Store x,y,zmax by node id.
        x: dict[int, float] = {}
        y: dict[int, float] = {}
        zmax: dict[int, float] = {}
        for _ in range(num_node):
            line = f.readline()
            if not line:
                raise MeshFormatError(f"unexpected EOF in node table: {path}")
            cols = line.split()
            if len(cols) < 5:
                raise MeshFormatError(f"invalid node row: {line!r}")
            nid = int(cols[0])
            x[nid] = float(cols[1])
            y[nid] = float(cols[2])
            zmax[nid] = float(cols[4])

    normals: dict[int, tuple[float, float, float]] = {}
    for eid, (n1, n2, n3) in elements.items():
        # Element.cpp uses zmax (surface elevation) for the terrain plane.
        p1x = x[n1]
        p1y = y[n1]
        p1z = zmax[n1]
        p2x = x[n2]
        p2y = y[n2]
        p2z = zmax[n2]
        p3x = x[n3]
        p3y = y[n3]
        p3z = zmax[n3]

        v1x = p2x - p1x
        v1y = p2y - p1y
        v1z = p2z - p1z
        v2x = p3x - p1x
        v2y = p3y - p1y
        v2z = p3z - p1z

        nx_raw = v1y * v2z - v1z * v2y
        ny_raw = v1z * v2x - v1x * v2z
        nz_raw = v1x * v2y - v1y * v2x

        nlen = math.sqrt(nx_raw * nx_raw + ny_raw * ny_raw + nz_raw * nz_raw)
        if nlen <= ZERO:
            nxv, nyv, nzv = 0.0, 0.0, 1.0
        else:
            nxv = nx_raw / nlen
            nyv = ny_raw / nlen
            nzv = nz_raw / nlen
            if nzv < 0.0:
                nxv = -nxv
                nyv = -nyv
                nzv = -nzv
        normals[int(eid)] = (float(nxv), float(nyv), float(nzv))

    return normals
