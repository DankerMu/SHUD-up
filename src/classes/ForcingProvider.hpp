//  ForcingProvider.hpp
//  SHUD
//
//  Forcing provider abstraction for supporting multiple forcing sources
//  (baseline CSV today; NetCDF in Phase A).
//
//  The forcing contract is the existing SHUD 5-variable forcing set:
//    Precip(mm/day), Temp(C), RH(0-1), Wind(m/s), RN(W/m2)
//
#ifndef ForcingProvider_hpp
#define ForcingProvider_hpp

#include "TimeSeriesData.hpp"

class ForcingProvider {
public:
    virtual ~ForcingProvider() = default;

    virtual int numStations() const = 0;
    virtual void movePointer(double t_min) = 0;

    // Step-function semantics: return the current-interval value after movePointer(t_min).
    // Column indices follow SHUD forcing macros: i_prcp..i_rn (1..5).
    virtual double get(int station_idx, int column) const = 0;

    // Forcing interval bounds (minutes) for the active pointer (used by TSR).
    virtual double currentTimeMin(int station_idx) const = 0;
    virtual double nextTimeMin(int station_idx) const = 0;

    // Station metadata (lon/lat in degrees; z in meters)
    virtual double lon(int station_idx) const = 0;
    virtual double lat(int station_idx) const = 0;
    virtual double z(int station_idx) const = 0;
};

class CsvForcingProvider final : public ForcingProvider {
public:
    CsvForcingProvider(_TimeSeriesData *tsd_weather, int num_forc)
        : tsd_weather_(tsd_weather), num_forc_(num_forc) {}

    int numStations() const override { return num_forc_; }

    void movePointer(double t_min) override
    {
        if (tsd_weather_ == nullptr || num_forc_ <= 0) {
            return;
        }
        for (int i = 0; i < num_forc_; i++) {
            tsd_weather_[i].movePointer(t_min);
        }
    }

    double get(int station_idx, int column) const override
    {
        return tsd_weather_[station_idx].getX(0.0, column);
    }

    double currentTimeMin(int station_idx) const override { return tsd_weather_[station_idx].currentTimeMin(); }
    double nextTimeMin(int station_idx) const override { return tsd_weather_[station_idx].nextTimeMin(); }

    double lon(int station_idx) const override { return tsd_weather_[station_idx].lon(); }
    double lat(int station_idx) const override { return tsd_weather_[station_idx].lat(); }
    double z(int station_idx) const override { return tsd_weather_[station_idx].xyz[2]; }

private:
    _TimeSeriesData *tsd_weather_ = nullptr;
    int num_forc_ = 0;
};

#endif /* ForcingProvider_hpp */

