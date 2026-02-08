#include "NetcdfOutputContext.hpp"

#ifdef _NETCDF_ON

#include "Model_Control.hpp"

#include <netcdf.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "funPlatform.hpp"

namespace {

static std::string trim(const std::string &s)
{
    const size_t first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const size_t last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

static std::string toUpper(std::string s)
{
    for (char &c : s) {
        c = (char)std::toupper((unsigned char)c);
    }
    return s;
}

static bool isAbsPath(const std::string &p)
{
    if (p.empty()) {
        return false;
    }
    if (p[0] == '/' || p[0] == '\\') {
        return true;
    }
    if (p.size() >= 2 && std::isalpha((unsigned char)p[0]) && p[1] == ':') {
        return true;
    }
    return false;
}

static std::string dirnameOf(const std::string &p)
{
    const size_t pos = p.find_last_of("/\\");
    if (pos == std::string::npos) {
        return ".";
    }
    if (pos == 0) {
        return "/";
    }
    return p.substr(0, pos);
}

static std::string basenameOf(const std::string &p)
{
    const size_t pos = p.find_last_of("/\\");
    if (pos == std::string::npos) {
        return p;
    }
    return p.substr(pos + 1);
}

static std::string joinPath(const std::string &a, const std::string &b)
{
    if (a.empty()) {
        return b;
    }
    if (b.empty()) {
        return a;
    }
    if (a.back() == '/') {
        return a + b;
    }
    return a + "/" + b;
}

static void ncCheck(int status, const char *what, const char *path)
{
    if (status == NC_NOERR) {
        return;
    }
    fprintf(stderr, "\n  Fatal Error: NetCDF output failure.\n");
    fprintf(stderr, "  Op: %s\n", what);
    if (path != nullptr && path[0] != '\0') {
        fprintf(stderr, "  File: %s\n", path);
    }
    fprintf(stderr, "  NetCDF: %s\n", nc_strerror(status));
    myexit(ERRFileIO);
}

static std::map<std::string, std::string> readKvCfg(const char *path)
{
    std::ifstream f(path);
    if (!f.is_open()) {
        fprintf(stderr, "\n  Fatal Error: Failed to open NetCDF output cfg.\n");
        fprintf(stderr, "  File: %s\n", path);
        myexit(ERRFileIO);
    }
    std::map<std::string, std::string> kv;
    std::string line;
    long lineNo = 0;
    while (std::getline(f, line)) {
        lineNo++;
        const std::string t = trim(line);
        if (t.empty() || t[0] == '#') {
            continue;
        }
        std::istringstream iss(t);
        std::string k;
        std::string v;
        if (!(iss >> k)) {
            continue;
        }
        if (!(iss >> v)) {
            fprintf(stderr, "\n  Fatal Error: Invalid KEY VALUE ncoutput cfg line (missing value).\n");
            fprintf(stderr, "  File: %s\n", path);
            fprintf(stderr, "  Line: %ld\n", lineNo);
            fprintf(stderr, "  Content: %s\n", t.c_str());
            myexit(ERRFileIO);
        }
        kv[toUpper(k)] = v;
    }
    return kv;
}

static bool parseYYYYMMDD(long yyyymmdd, int &y, int &m, int &d)
{
    if (yyyymmdd <= 0) {
        return false;
    }
    y = (int)(yyyymmdd / 10000);
    m = (int)((yyyymmdd / 100) % 100);
    d = (int)(yyyymmdd % 100);
    if (m < 1 || m > 12 || d < 1 || d > 31) {
        return false;
    }
    return true;
}

static std::string isoDateFromYYYYMMDD(long yyyymmdd)
{
    int y = 0;
    int m = 0;
    int d = 0;
    if (!parseYYYYMMDD(yyyymmdd, y, m, d)) {
        return "0000-00-00";
    }
    char buf[32];
    snprintf(buf, sizeof(buf), "%04d-%02d-%02d", y, m, d);
    return std::string(buf);
}

static std::string deriveVarName(const char *legacy_basename)
{
    if (legacy_basename == nullptr || legacy_basename[0] == '\0') {
        fprintf(stderr, "\n  Fatal Error: Empty legacy_basename for NetCDF output.\n");
        myexit(ERRFileIO);
    }
    const std::string leaf = basenameOf(std::string(legacy_basename));
    const size_t dot = leaf.find_last_of('.');
    if (dot == std::string::npos || dot + 1 >= leaf.size()) {
        fprintf(stderr, "\n  Fatal Error: Cannot derive NetCDF variable name from legacy basename.\n");
        fprintf(stderr, "  legacy_basename: %s\n", legacy_basename);
        fprintf(stderr, "  Expected: <prefix>.<varname> (e.g., qhh.eleygw)\n");
        myexit(ERRFileIO);
    }
    return leaf.substr(dot + 1);
}

static std::string derivePrefixLeaf(const char *legacy_basename)
{
    const std::string leaf = basenameOf(std::string(legacy_basename));
    const size_t dot = leaf.find_last_of('.');
    if (dot == std::string::npos || dot == 0) {
        fprintf(stderr, "\n  Fatal Error: Cannot derive NetCDF output prefix from legacy basename.\n");
        fprintf(stderr, "  legacy_basename: %s\n", legacy_basename);
        fprintf(stderr, "  Expected: <prefix>.<varname> (e.g., qhh.eleygw)\n");
        myexit(ERRFileIO);
    }
    return leaf.substr(0, dot);
}

class NetcdfElementFile {
public:
    explicit NetcdfElementFile(std::string out_dir_abs, float fill_value)
        : out_dir_abs_(std::move(out_dir_abs)), fill_value_(fill_value)
    {
    }

    ~NetcdfElementFile() { close(); }

    void openIfNeeded(const char *legacy_basename, long forc_start_yyyymmdd, int n_all)
    {
        if (ncid_ != -1) {
            if (n_all_ != n_all || forc_start_yyyymmdd_ != forc_start_yyyymmdd) {
                fprintf(stderr, "\n  Fatal Error: Inconsistent NetCDF element output context.\n");
                fprintf(stderr, "  Existing: n_all=%d ForcStartTime=%ld\n", n_all_, forc_start_yyyymmdd_);
                fprintf(stderr, "  New:      n_all=%d ForcStartTime=%ld\n", n_all, forc_start_yyyymmdd);
                myexit(ERRFileIO);
            }
            return;
        }

        n_all_ = n_all;
        forc_start_yyyymmdd_ = forc_start_yyyymmdd;

        const std::string prefix_leaf = derivePrefixLeaf(legacy_basename);

        std::string out_dir = out_dir_abs_.empty() ? dirnameOf(std::string(legacy_basename)) : out_dir_abs_;
        if (out_dir.empty()) {
            out_dir = ".";
        }

        std::vector<char> out_dir_mut(out_dir.begin(), out_dir.end());
        out_dir_mut.push_back('\0');
        mkdir_p(out_dir_mut.data(), 0777);

        file_path_ = joinPath(out_dir, prefix_leaf + ".ele.nc");

        const int cmode = NC_NETCDF4 | NC_CLASSIC_MODEL | NC_CLOBBER;
        ncCheck(nc_create(file_path_.c_str(), cmode, &ncid_), "nc_create", file_path_.c_str());

        ncCheck(nc_def_dim(ncid_, "time", NC_UNLIMITED, &dim_time_), "nc_def_dim(time)", file_path_.c_str());
        ncCheck(nc_def_dim(ncid_, "mesh_face", (size_t)n_all_, &dim_face_), "nc_def_dim(mesh_face)", file_path_.c_str());

        const int time_dims[1] = {dim_time_};
        ncCheck(nc_def_var(ncid_, "time", NC_DOUBLE, 1, time_dims, &var_time_), "nc_def_var(time)", file_path_.c_str());

        const std::string iso = isoDateFromYYYYMMDD(forc_start_yyyymmdd_);
        const std::string units = "minutes since " + iso + " 00:00:00 UTC";
        ncCheck(nc_put_att_text(ncid_, var_time_, "units", units.size(), units.c_str()),
                "nc_put_att_text(time.units)",
                file_path_.c_str());
        const char *calendar = "standard";
        ncCheck(nc_put_att_text(ncid_, var_time_, "calendar", strlen(calendar), calendar),
                "nc_put_att_text(time.calendar)",
                file_path_.c_str());

        ncCheck(nc_enddef(ncid_), "nc_enddef", file_path_.c_str());
    }

    int ensureVar(const std::string &name)
    {
        const auto it = varids_.find(name);
        if (it != varids_.end()) {
            return it->second;
        }
        if (ncid_ == -1) {
            fprintf(stderr, "\n  Fatal Error: NetCDF element file is not initialized before defining variables.\n");
            myexit(ERRFileIO);
        }

        ncCheck(nc_redef(ncid_), "nc_redef", file_path_.c_str());
        const int dims[2] = {dim_time_, dim_face_};
        int varid = -1;
        ncCheck(nc_def_var(ncid_, name.c_str(), NC_FLOAT, 2, dims, &varid),
                "nc_def_var(data)",
                file_path_.c_str());
        ncCheck(nc_put_att_float(ncid_, varid, "_FillValue", NC_FLOAT, 1, &fill_value_),
                "nc_put_att_float(_FillValue)",
                file_path_.c_str());
        ncCheck(nc_enddef(ncid_), "nc_enddef", file_path_.c_str());

        varids_[name] = varid;
        return varid;
    }

    size_t ensureTimeIndex(long long t_min)
    {
        const auto it = time_to_idx_.find(t_min);
        if (it != time_to_idx_.end()) {
            return it->second;
        }
        if (ncid_ == -1) {
            fprintf(stderr, "\n  Fatal Error: NetCDF element file is not initialized before writing time.\n");
            myexit(ERRFileIO);
        }

        const size_t idx = times_.size();
        times_.push_back((double)t_min);
        time_to_idx_[t_min] = idx;

        const size_t start[1] = {idx};
        const double v = (double)t_min;
        ncCheck(nc_put_var1_double(ncid_, var_time_, start, &v), "nc_put_var1_double(time)", file_path_.c_str());
        return idx;
    }

    void writeVar(int varid, size_t time_idx, const std::vector<float> &data)
    {
        if (ncid_ == -1) {
            fprintf(stderr, "\n  Fatal Error: NetCDF element file is not initialized before writing data.\n");
            myexit(ERRFileIO);
        }
        if ((int)data.size() != n_all_) {
            fprintf(stderr, "\n  Fatal Error: NetCDF write buffer size mismatch.\n");
            fprintf(stderr, "  Expected n_all=%d, got=%zu\n", n_all_, data.size());
            myexit(ERRFileIO);
        }

        const size_t start[2] = {time_idx, 0};
        const size_t count[2] = {1, (size_t)n_all_};
        ncCheck(nc_put_vara_float(ncid_, varid, start, count, data.data()), "nc_put_vara_float(data)", file_path_.c_str());
    }

    void close()
    {
        if (ncid_ != -1) {
            nc_close(ncid_);
            ncid_ = -1;
        }
    }

private:
    std::string out_dir_abs_;
    float fill_value_ = 9.96921e36f;

    std::string file_path_;
    int ncid_ = -1;
    int dim_time_ = -1;
    int dim_face_ = -1;
    int var_time_ = -1;
    int n_all_ = 0;
    long forc_start_yyyymmdd_ = 0;

    std::map<std::string, int> varids_;
    std::map<long long, size_t> time_to_idx_;
    std::vector<double> times_;
};

class NetcdfElementVarSink final : public IPrintSink {
public:
    explicit NetcdfElementVarSink(NetcdfElementFile *file) : file_(file) {}

    void onInit(const char *legacy_basename,
                long start_yyyymmdd,
                int interval_min,
                int n_all,
                int num_selected,
                const double *icol_1based,
                int /* radiation_input_mode */,
                int /* terrain_radiation */,
                int /* solar_lonlat_mode */,
                double /* solar_lon_deg */,
                double /* solar_lat_deg */) override
    {
        (void)interval_min;
        var_name_ = deriveVarName(legacy_basename);
        n_all_ = n_all;
        icol_.clear();
        icol_.reserve((size_t)std::max(0, num_selected));
        for (int i = 0; i < num_selected; i++) {
            icol_.push_back((int)llround(icol_1based[i]));
        }

        file_->openIfNeeded(legacy_basename, start_yyyymmdd, n_all_);
        varid_ = file_->ensureVar(var_name_);
        initialized_ = true;
    }

    void onWrite(double t_quantized_min, int num_selected, const double *buffer) override
    {
        if (!initialized_ || file_ == nullptr || buffer == nullptr) {
            return;
        }

        const long long t_min = (long long)llround(t_quantized_min);
        const size_t time_idx = file_->ensureTimeIndex(t_min);

        std::vector<float> out((size_t)n_all_, fill_value_);
        const int nsel = std::min<int>(num_selected, (int)icol_.size());
        for (int j = 0; j < nsel; j++) {
            const int idx0 = icol_[(size_t)j] - 1;
            if (idx0 < 0 || idx0 >= n_all_) {
                continue;
            }
            out[(size_t)idx0] = (float)buffer[j];
        }

        file_->writeVar(varid_, time_idx, out);
    }

    void onClose() override {}

private:
    NetcdfElementFile *file_ = nullptr;
    bool initialized_ = false;

    std::string var_name_;
    int varid_ = -1;
    int n_all_ = 0;
    float fill_value_ = 9.96921e36f;

    std::vector<int> icol_;
};

struct ParsedNcOutputCfg {
    std::string schema_abs;
    std::string out_dir_abs;
};

static ParsedNcOutputCfg parseNcOutputCfg(const std::string &cfg_path)
{
    ParsedNcOutputCfg cfg;
    const std::map<std::string, std::string> kv = readKvCfg(cfg_path.c_str());

    auto get = [&](const char *k) -> std::string {
        auto it = kv.find(toUpper(k));
        if (it == kv.end()) {
            return "";
        }
        return it->second;
    };

    const std::string cfg_dir = dirnameOf(cfg_path);
    const std::string input_dir = dirnameOf(cfg_dir);
    const std::string run_dir = dirnameOf(input_dir);

    const std::string schema_cfg = get("SCHEMA");
    if (!schema_cfg.empty()) {
        cfg.schema_abs = isAbsPath(schema_cfg) ? schema_cfg : joinPath(run_dir, schema_cfg);
    }

    const std::string out_dir_cfg = get("OUT_DIR");
    if (!out_dir_cfg.empty()) {
        cfg.out_dir_abs = isAbsPath(out_dir_cfg) ? out_dir_cfg : joinPath(run_dir, out_dir_cfg);
    }

    return cfg;
}

} // namespace

struct NetcdfOutputContext::Impl {
    std::string ncoutput_cfg_path;
    ParsedNcOutputCfg cfg;
    NetcdfElementFile ele_file;
    std::vector<std::unique_ptr<IPrintSink>> sinks;

    explicit Impl(const char *path)
        : ncoutput_cfg_path(path ? path : ""),
          cfg(parseNcOutputCfg(ncoutput_cfg_path)),
          ele_file(cfg.out_dir_abs, 9.96921e36f)
    {
    }

    IPrintSink *createElementSink()
    {
        auto s = std::make_unique<NetcdfElementVarSink>(&ele_file);
        IPrintSink *ptr = s.get();
        sinks.push_back(std::move(s));
        return ptr;
    }
};

NetcdfOutputContext::NetcdfOutputContext(const char *ncoutput_cfg_path)
    : impl_(std::make_unique<Impl>(ncoutput_cfg_path))
{
}

NetcdfOutputContext::~NetcdfOutputContext() = default;

IPrintSink *NetcdfOutputContext::createElementSink()
{
    return impl_->createElementSink();
}

#endif /* _NETCDF_ON */

