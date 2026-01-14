#include "SolarRadiation.hpp"

#include "TimeContext.hpp"

#include <algorithm>
#include <cmath>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kDeg2Rad = kPi / 180.0;

inline bool isFinite(double x) {
    return std::isfinite(x) != 0;
}

inline double clamp(double x, double lo, double hi) {
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

inline double wrapDegrees180(double lon_deg) {
    if (!isFinite(lon_deg)) {
        return 0.0;
    }
    double x = std::fmod(lon_deg, 360.0);
    if (x > 180.0) {
        x -= 360.0;
    } else if (x < -180.0) {
        x += 360.0;
    }
    return x;
}

inline double clampLatitude(double lat_deg) {
    if (!isFinite(lat_deg)) {
        return 0.0;
    }
    return clamp(lat_deg, -90.0, 90.0);
}

inline double minuteOfDay(double t_min) {
    if (!isFinite(t_min)) {
        return 0.0;
    }
    double m = std::fmod(t_min, 1440.0);
    if (m < 0.0) {
        m += 1440.0;
    }
    return m;
}

inline double wrapMinutes1440(double minutes) {
    if (!isFinite(minutes)) {
        return 0.0;
    }
    double m = std::fmod(minutes, 1440.0);
    if (m < 0.0) {
        m += 1440.0;
    }
    return m;
}

inline double wrapRadians2Pi(double rad) {
    if (!isFinite(rad)) {
        return 0.0;
    }
    double x = std::fmod(rad, kTwoPi);
    if (x < 0.0) {
        x += kTwoPi;
    }
    return x;
}

inline double approximateTimezoneHours(double lon_deg) {
    if (!isFinite(lon_deg)) {
        return 0.0;
    }
    return std::round(lon_deg / 15.0);
}

inline double safeAcos(double x) {
    return std::acos(clamp(x, -1.0, 1.0));
}

SolarPosition solarPositionImpl(double t_min,
                                double lat_deg,
                                double lon_deg,
                                const TimeContext& tc,
                                double timezone_hours,
                                bool timezone_provided) {
    SolarPosition sp{};
    sp.cosZ = 0.0;
    sp.zenith = kPi / 2.0;
    sp.azimuth = 0.0;
    sp.declination = 0.0;
    sp.hourAngle = 0.0;

    const double lat = clampLatitude(lat_deg);
    const double lon = wrapDegrees180(lon_deg);

    double tz = timezone_hours;
    if (!timezone_provided || !isFinite(tz)) {
        tz = approximateTimezoneHours(lon);
    }

    int doy = tc.julianDay(t_min);
    if (doy < 1 || doy > 366) {
        doy = 1;
    }

    const double mod_min = minuteOfDay(t_min);
    const double hour = mod_min / 60.0;

    const double gamma = (kTwoPi / 365.0) * (static_cast<double>(doy - 1) + (hour - 12.0) / 24.0);
    const double sin_g = std::sin(gamma);
    const double cos_g = std::cos(gamma);
    const double sin_2g = std::sin(2.0 * gamma);
    const double cos_2g = std::cos(2.0 * gamma);
    const double sin_3g = std::sin(3.0 * gamma);
    const double cos_3g = std::cos(3.0 * gamma);

    const double eq_time_min =
        229.18 * (0.000075 + 0.001868 * cos_g - 0.032077 * sin_g - 0.014615 * cos_2g - 0.040849 * sin_2g);

    const double decl =
        0.006918 - 0.399912 * cos_g + 0.070257 * sin_g - 0.006758 * cos_2g + 0.000907 * sin_2g - 0.002697 * cos_3g +
        0.00148 * sin_3g;
    sp.declination = decl;

    const double time_offset_min = eq_time_min + 4.0 * lon - 60.0 * tz;
    const double true_solar_time_min = wrapMinutes1440(mod_min + time_offset_min);

    const double ha_deg = true_solar_time_min / 4.0 - 180.0;
    const double ha = ha_deg * kDeg2Rad;
    sp.hourAngle = ha;

    const double lat_rad = lat * kDeg2Rad;
    const double sin_lat = std::sin(lat_rad);
    const double cos_lat = std::cos(lat_rad);

    const double sin_decl = std::sin(decl);
    const double cos_decl = std::cos(decl);
    const double sin_ha = std::sin(ha);
    const double cos_ha = std::cos(ha);

    const double cosz_raw = sin_lat * sin_decl + cos_lat * cos_decl * cos_ha;
    const double cosz = clamp(cosz_raw, -1.0, 1.0);
    sp.cosZ = cosz;
    sp.zenith = safeAcos(cosz);

    const double east = -cos_decl * sin_ha;
    const double north = cos_lat * sin_decl - sin_lat * cos_decl * cos_ha;
    const double az = std::atan2(east, north);
    sp.azimuth = wrapRadians2Pi(az);

    if (!isFinite(sp.cosZ) || !isFinite(sp.zenith) || !isFinite(sp.azimuth) || !isFinite(sp.declination) || !isFinite(sp.hourAngle)) {
        SolarPosition safe{};
        safe.cosZ = 0.0;
        safe.zenith = kPi / 2.0;
        safe.azimuth = 0.0;
        safe.declination = 0.0;
        safe.hourAngle = 0.0;
        return safe;
    }

    return sp;
}

} // namespace

SolarPosition solarPosition(double t_min, double lat_deg, double lon_deg, const TimeContext& tc) {
    return solarPositionImpl(t_min, lat_deg, lon_deg, tc, 0.0, false);
}

SolarPosition solarPosition(double t_min, double lat_deg, double lon_deg, const TimeContext& tc, double timezone_hours) {
    return solarPositionImpl(t_min, lat_deg, lon_deg, tc, timezone_hours, true);
}

double terrainFactor(double nx, double ny, double nz, const SolarPosition& sp, double cap, double cosz_min) {
    if (!isFinite(nx) || !isFinite(ny) || !isFinite(nz) || !isFinite(sp.cosZ) || !isFinite(sp.azimuth)) {
        return 0.0;
    }

    if (!isFinite(cap)) {
        // Treat missing/invalid cap as "effectively uncapped" while still keeping a finite bound.
        cap = 10.0;
    } else if (cap < 0.0) {
        cap = 0.0;
    }
    if (!isFinite(cosz_min) || cosz_min < 0.0) {
        cosz_min = 0.0;
    }
    cosz_min = clamp(cosz_min, 0.0, 1.0);

    const double cosz = sp.cosZ;
    if (!(cosz > 0.0)) {
        return 0.0;
    }

    const double nlen = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (!(nlen > 0.0) || !isFinite(nlen)) {
        nx = 0.0;
        ny = 0.0;
        nz = 1.0;
    } else {
        nx /= nlen;
        ny /= nlen;
        nz /= nlen;
        if (nz < 0.0) {
            nx = -nx;
            ny = -ny;
            nz = -nz;
        }
    }

    const double cosz_clamped = clamp(cosz, -1.0, 1.0);
    const double sinz = std::sqrt(std::max(0.0, 1.0 - cosz_clamped * cosz_clamped));

    const double sin_az = std::sin(sp.azimuth);
    const double cos_az = std::cos(sp.azimuth);
    const double sx = sinz * sin_az;
    const double sy = sinz * cos_az;
    const double sz = cosz_clamped;

    const double cosi = nx * sx + ny * sy + nz * sz;
    if (!(cosi > 0.0) || !isFinite(cosi)) {
        return 0.0;
    }

    double denom = cosz_clamped;
    if (denom < cosz_min) {
        denom = cosz_min;
    }
    if (!(denom > 0.0) || !isFinite(denom)) {
        return 0.0;
    }

    double factor = cosi / denom;
    if (!isFinite(factor) || factor <= 0.0) {
        return 0.0;
    }

    if (factor > cap) {
        factor = cap;
    }
    if (!isFinite(factor) || factor <= 0.0) {
        return 0.0;
    }

    return factor;
}
