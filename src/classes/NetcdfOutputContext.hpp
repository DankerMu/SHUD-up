//  NetcdfOutputContext.hpp
//  SHUD
//
//  NetCDF output helpers (Phase B)
//
#ifndef NetcdfOutputContext_hpp
#define NetcdfOutputContext_hpp

#ifdef _NETCDF_ON

#include <memory>

class IPrintSink;

class NetcdfOutputContext {
public:
    explicit NetcdfOutputContext(const char *ncoutput_cfg_path);
    ~NetcdfOutputContext();

    IPrintSink *createElementSink();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

#endif /* _NETCDF_ON */

#endif /* NetcdfOutputContext_hpp */

