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
class _Element;
class _Node;

class NetcdfOutputContext {
public:
    NetcdfOutputContext(const char *ncoutput_cfg_path,
                        const _Node *nodes,
                        int num_nodes,
                        const _Element *elements,
                        int num_elements);
    ~NetcdfOutputContext();

    IPrintSink *createElementSink();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

#endif /* _NETCDF_ON */

#endif /* NetcdfOutputContext_hpp */
