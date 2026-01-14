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

SolarPosition solarPosition(double t_min,
                            double lat_deg,
                            double lon_deg,
                            const TimeContext& tc);

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

