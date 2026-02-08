//  NetcdfForcingProvider.hpp
//  SHUD
//
//  NetCDF forcing providers (Phase A)
//
#ifndef NetcdfForcingProvider_hpp
#define NetcdfForcingProvider_hpp

#include "ForcingProvider.hpp"

#include <string>
#include <vector>

#ifdef _NETCDF_ON

struct ForcingStationMeta {
    double lon_deg = NA_VALUE;
    double lat_deg = NA_VALUE;
    double z_m = NA_VALUE;
};

class NetcdfForcingProvider final : public ForcingProvider {
public:
    NetcdfForcingProvider(const char *forcing_cfg_path,
                          const std::vector<ForcingStationMeta> &stations,
                          long forc_start_yyyymmdd,
                          double sim_start_min,
                          double sim_end_min);

    ~NetcdfForcingProvider() override;

    int numStations() const override;
    void movePointer(double t_min) override;
    double get(int station_idx, int column) const override;

    double currentTimeMin(int station_idx) const override;
    double nextTimeMin(int station_idx) const override;

    double lon(int station_idx) const override;
    double lat(int station_idx) const override;
    double z(int station_idx) const override;

    // Coverage (minutes relative to ForcStartTime), used for validation.
    double minTimeMin() const;
    double maxTimeMin() const;

private:
    struct Impl;
    Impl *impl_;
};

#endif /* _NETCDF_ON */

#endif /* NetcdfForcingProvider_hpp */

