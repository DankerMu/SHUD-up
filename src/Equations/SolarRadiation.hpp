//  SolarRadiation.hpp
//  Created by SHUD-up contributors.
//
#ifndef SolarRadiation_hpp
#define SolarRadiation_hpp

class TimeContext;

struct SolarPosition {
    double cosZ;        /* cos(zenith) [-] */
    double zenith;      /* zenith angle [rad], [0, pi] */
    double azimuth;     /* azimuth [rad], [0, 2*pi), North=0, East=pi/2 */
    double declination; /* solar declination [rad] */
    double hourAngle;   /* hour angle [rad], 0 at solar noon, positive afternoon */
};

/*
 * Compute solar position at simulation time.
 *
 * t_min is simulation time in minutes (not local clock time).
 *
 * Note: this overload estimates timezone from longitude using round(lon_deg / 15).
 * This geometric estimate can be ~1 hour off in political timezones.
 */
SolarPosition solarPosition(double t_min,
                            double lat_deg,
                            double lon_deg,
                            const TimeContext& tc);

/*
 * Compute solar position at simulation time.
 *
 * t_min is simulation time in minutes (not local clock time).
 *
 * timezone_hours is the UTC offset in hours (e.g., UTC+8 = 8.0).
 */
SolarPosition solarPosition(double t_min,
                            double lat_deg,
                            double lon_deg,
                            const TimeContext& tc,
                            double timezone_hours);

double terrainFactor(double nx,
                     double ny,
                     double nz,
                     const SolarPosition& sp,
                     double cap,
                     double cosz_min);

#endif /* SolarRadiation_hpp */
