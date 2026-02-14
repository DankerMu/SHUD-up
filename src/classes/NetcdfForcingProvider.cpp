#include "NetcdfForcingProvider.hpp"

#ifdef _NETCDF_ON

#include <netcdf.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <glob.h>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

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

static std::string toLower(std::string s)
{
    for (char &c : s) {
        c = (char)std::tolower((unsigned char)c);
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

static void replaceAll(std::string &s, const std::string &from, const std::string &to)
{
    if (from.empty()) {
        return;
    }
    size_t start = 0;
    while ((start = s.find(from, start)) != std::string::npos) {
        s.replace(start, from.length(), to);
        start += to.length();
    }
}

static std::string zfillInt(int v, int width)
{
    std::ostringstream oss;
    oss << v;
    std::string s = oss.str();
    if ((int)s.size() < width) {
        s = std::string((size_t)(width - (int)s.size()), '0') + s;
    }
    return s;
}

static std::string resolveSingleGlob(const std::string &pattern)
{
    glob_t g;
    memset(&g, 0, sizeof(g));
    const int rc = glob(pattern.c_str(), 0, nullptr, &g);
    if (rc != 0 || g.gl_pathc == 0) {
        globfree(&g);
        fprintf(stderr, "\n  Fatal Error: NetCDF forcing file not found.\n");
        fprintf(stderr, "  Glob: %s\n", pattern.c_str());
        myexit(ERRFileIO);
    }
    std::vector<std::string> matches;
    matches.reserve(g.gl_pathc);
    for (size_t i = 0; i < g.gl_pathc; i++) {
        matches.emplace_back(g.gl_pathv[i]);
    }
    globfree(&g);
    std::sort(matches.begin(), matches.end());
    if (matches.size() != 1) {
        fprintf(stderr, "\n  Fatal Error: NetCDF forcing file glob is ambiguous (%zu matches).\n", matches.size());
        fprintf(stderr, "  Glob: %s\n", pattern.c_str());
        const size_t nshow = std::min<size_t>(matches.size(), 5);
        for (size_t i = 0; i < nshow; i++) {
            fprintf(stderr, "  Match[%zu]: %s\n", i, matches[i].c_str());
        }
        fprintf(stderr, "  Fix: narrow the pattern to match exactly one file.\n");
        myexit(ERRFileIO);
    }
    return matches[0];
}

// Howard Hinnant's civil calendar algorithms (public domain)
static long long daysFromCivil(int y, unsigned m, unsigned d)
{
    y -= m <= 2;
    const int era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = (unsigned)(y - era * 400);                // [0, 399]
    const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1; // [0, 365]
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;    // [0, 146096]
    return (long long)era * 146097LL + (long long)doe - 719468LL;  // days since 1970-01-01
}

static void civilFromDays(long long z, int &y, unsigned &m, unsigned &d)
{
    z += 719468LL;
    const long long era = (z >= 0 ? z : z - 146096) / 146097;
    const unsigned doe = (unsigned)(z - era * 146097);             // [0, 146096]
    const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365; // [0, 399]
    y = (int)yoe + (int)era * 400;
    const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);  // [0, 365]
    const unsigned mp = (5 * doy + 2) / 153;                       // [0, 11]
    d = doy - (153 * mp + 2) / 5 + 1;                              // [1, 31]
    m = mp + (mp < 10 ? 3 : -9);                                   // [1, 12]
    y += (m <= 2);
}

static bool parseYYYYMMDD(long yyyymmdd, int &y, unsigned &m, unsigned &d)
{
    if (yyyymmdd <= 0) {
        return false;
    }
    y = (int)(yyyymmdd / 10000);
    m = (unsigned)((yyyymmdd / 100) % 100);
    d = (unsigned)(yyyymmdd % 100);
    if (m < 1 || m > 12 || d < 1 || d > 31) {
        return false;
    }
    return true;
}

struct UnitsSince {
    double factor_to_minutes = 0.0; // multiply unit value by this to get minutes
    long long base_minutes_since_epoch = 0; // epoch: 1970-01-01 00:00 UTC
};

static bool parseIsoDateTime(const std::string &s, int &y, unsigned &m, unsigned &d, int &hh, int &mm, double &ss)
{
    // Accept:
    //  YYYY-MM-DD
    //  YYYY-MM-DD HH:MM
    //  YYYY-MM-DD HH:MM:SS(.fff)
    //  YYYY-MM-DDTHH:MM:SS
    const std::string t = trim(s);
    if (t.size() < 10) {
        return false;
    }
    if (sscanf(t.c_str(), "%d-%u-%u", &y, &m, &d) < 3) {
        return false;
    }
    hh = 0;
    mm = 0;
    ss = 0.0;
    if (t.size() == 10) {
        return true;
    }

    size_t pos = 10;
    while (pos < t.size() && (t[pos] == ' ' || t[pos] == 'T')) {
        pos++;
    }
    if (pos >= t.size()) {
        return true;
    }

    int h = 0;
    int mi = 0;
    double se = 0.0;
    const char *p = t.c_str() + pos;
    const int n = sscanf(p, "%d:%d:%lf", &h, &mi, &se);
    if (n < 2) {
        return false;
    }
    hh = h;
    mm = mi;
    ss = (n >= 3) ? se : 0.0;
    return true;
}

static UnitsSince parseUnitsSince(const std::string &units_attr)
{
    const std::string u = toLower(trim(units_attr));
    const size_t since_pos = u.find("since");
    if (since_pos == std::string::npos) {
        fprintf(stderr, "\n  Fatal Error: NetCDF time.units missing 'since'.\n");
        fprintf(stderr, "  units: %s\n", units_attr.c_str());
        myexit(ERRFileIO);
    }
    const std::string unit_part = trim(u.substr(0, since_pos));
    const std::string base_part = trim(u.substr(since_pos + 5));

    double factor_to_min = 0.0;
    if (unit_part.find("second") == 0) {
        factor_to_min = 1.0 / 60.0;
    } else if (unit_part.find("minute") == 0) {
        factor_to_min = 1.0;
    } else if (unit_part.find("hour") == 0) {
        factor_to_min = 60.0;
    } else if (unit_part.find("day") == 0) {
        factor_to_min = 1440.0;
    } else {
        fprintf(stderr, "\n  Fatal Error: Unsupported NetCDF time unit.\n");
        fprintf(stderr, "  units: %s\n", units_attr.c_str());
        myexit(ERRFileIO);
    }

    int y = 0;
    unsigned m = 0;
    unsigned d = 0;
    int hh = 0;
    int mm = 0;
    double ss = 0.0;
    if (!parseIsoDateTime(base_part, y, m, d, hh, mm, ss)) {
        fprintf(stderr, "\n  Fatal Error: Failed to parse NetCDF time base date.\n");
        fprintf(stderr, "  units: %s\n", units_attr.c_str());
        myexit(ERRFileIO);
    }

    const long long base_days = daysFromCivil(y, m, d);
    const long long base_minutes = base_days * 1440LL + (long long)hh * 60LL + (long long)mm;
    UnitsSince out;
    out.factor_to_minutes = factor_to_min;
    out.base_minutes_since_epoch = base_minutes;
    // Note: ignore seconds; they are rarely used in the current forcing products.
    (void)ss;
    return out;
}

static std::string ncStrError(int status)
{
    const char *msg = nc_strerror(status);
    return msg ? std::string(msg) : std::string("unknown netcdf error");
}

static std::string ncGetAttText(int ncid, int varid, const char *name)
{
    size_t len = 0;
    const int rc_len = nc_inq_attlen(ncid, varid, name, &len);
    if (rc_len != NC_NOERR) {
        return "";
    }
    std::string out(len, '\0');
    const int rc = nc_get_att_text(ncid, varid, name, &out[0]);
    if (rc != NC_NOERR) {
        return "";
    }
    // NetCDF text attrs are not guaranteed to be NUL-terminated.
    return trim(out);
}

static bool ncGetAttDouble(int ncid, int varid, const char *name, double &out_val)
{
    const int rc = nc_get_att_double(ncid, varid, name, &out_val);
    return rc == NC_NOERR;
}

struct NcVarPointReader {
    std::string file;
    int ncid = -1;
    int varid = -1;
    int ndims = 0;
    int time_dim_pos = -1;
    int lat_dim_pos = -1;
    int lon_dim_pos = -1;

    bool has_scale = false;
    bool has_offset = false;
    double scale = 1.0;
    double offset = 0.0;

    bool has_fill = false;
    bool has_missing = false;
    double fill = 0.0;
    double missing = 0.0;

    std::string var_name;

    void close()
    {
        if (ncid >= 0) {
            nc_close(ncid);
            ncid = -1;
        }
    }
};

static NcVarPointReader openPointReader(const std::string &file,
                                        const std::string &var_name,
                                        const std::string &dim_time,
                                        const std::string &dim_lat,
                                        const std::string &dim_lon)
{
    NcVarPointReader r;
    r.file = file;
    r.var_name = var_name;

    int ncid = -1;
    int rc = nc_open(file.c_str(), NC_NOWRITE, &ncid);
    if (rc != NC_NOERR) {
        fprintf(stderr, "\n  Fatal Error: Failed to open NetCDF file.\n");
        fprintf(stderr, "  File: %s\n", file.c_str());
        fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
        myexit(ERRFileIO);
    }
    r.ncid = ncid;

    int varid = -1;
    rc = nc_inq_varid(ncid, var_name.c_str(), &varid);
    if (rc != NC_NOERR) {
        fprintf(stderr, "\n  Fatal Error: NetCDF variable not found.\n");
        fprintf(stderr, "  File: %s\n", file.c_str());
        fprintf(stderr, "  Var: %s\n", var_name.c_str());
        fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
        myexit(ERRFileIO);
    }
    r.varid = varid;

    char name_buf[NC_MAX_NAME + 1];
    nc_type xtype;
    int ndims = 0;
    int dimids[NC_MAX_VAR_DIMS];
    int natts = 0;
    rc = nc_inq_var(ncid, varid, name_buf, &xtype, &ndims, dimids, &natts);
    if (rc != NC_NOERR) {
        fprintf(stderr, "\n  Fatal Error: nc_inq_var failed.\n");
        fprintf(stderr, "  File: %s\n", file.c_str());
        fprintf(stderr, "  Var: %s\n", var_name.c_str());
        fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
        myexit(ERRFileIO);
    }
    r.ndims = ndims;
    if (ndims != 3) {
        fprintf(stderr, "\n  Fatal Error: NetCDF forcing variable must be 3-D (time/lat/lon).\n");
        fprintf(stderr, "  File: %s\n", file.c_str());
        fprintf(stderr, "  Var: %s\n", var_name.c_str());
        fprintf(stderr, "  ndims: %d\n", ndims);
        myexit(ERRFileIO);
    }

    for (int i = 0; i < ndims; i++) {
        char dimname[NC_MAX_NAME + 1];
        rc = nc_inq_dimname(ncid, dimids[i], dimname);
        if (rc != NC_NOERR) {
            fprintf(stderr, "\n  Fatal Error: nc_inq_dimname failed.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Var: %s\n", var_name.c_str());
            fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
            myexit(ERRFileIO);
        }
        const std::string dn(dimname);
        if (dn == dim_time) {
            r.time_dim_pos = i;
        } else if (dn == dim_lat) {
            r.lat_dim_pos = i;
        } else if (dn == dim_lon) {
            r.lon_dim_pos = i;
        }
    }
    if (r.time_dim_pos < 0 || r.lat_dim_pos < 0 || r.lon_dim_pos < 0) {
        fprintf(stderr, "\n  Fatal Error: NetCDF forcing variable dims must include time/lat/lon by name.\n");
        fprintf(stderr, "  File: %s\n", file.c_str());
        fprintf(stderr, "  Var: %s\n", var_name.c_str());
        fprintf(stderr, "  Expect dims: time='%s', lat='%s', lon='%s'\n",
                dim_time.c_str(), dim_lat.c_str(), dim_lon.c_str());
        fprintf(stderr, "  Got dim positions: time=%d lat=%d lon=%d\n",
                r.time_dim_pos, r.lat_dim_pos, r.lon_dim_pos);
        myexit(ERRFileIO);
    }

    double scale = 1.0;
    if (ncGetAttDouble(ncid, varid, "scale_factor", scale)) {
        r.has_scale = true;
        r.scale = scale;
    }
    double offset = 0.0;
    if (ncGetAttDouble(ncid, varid, "add_offset", offset)) {
        r.has_offset = true;
        r.offset = offset;
    }
    double fill = 0.0;
    if (ncGetAttDouble(ncid, varid, "_FillValue", fill)) {
        r.has_fill = true;
        r.fill = fill;
    }
    double missing = 0.0;
    if (ncGetAttDouble(ncid, varid, "missing_value", missing)) {
        r.has_missing = true;
        r.missing = missing;
    }

    return r;
}

static double readPoint(const NcVarPointReader &r,
                        size_t time_idx,
                        size_t lat_idx,
                        size_t lon_idx,
                        int station_idx,
                        double station_lon,
                        double station_lat,
                        double grid_lon,
                        double grid_lat)
{
    std::vector<size_t> start((size_t)r.ndims, 0);
    std::vector<size_t> count((size_t)r.ndims, 1);
    start[(size_t)r.time_dim_pos] = time_idx;
    start[(size_t)r.lat_dim_pos] = lat_idx;
    start[(size_t)r.lon_dim_pos] = lon_idx;

    double raw = std::numeric_limits<double>::quiet_NaN();
    const int rc = nc_get_vara_double(r.ncid, r.varid, start.data(), count.data(), &raw);
    if (rc != NC_NOERR) {
        fprintf(stderr, "\n  Fatal Error: Failed to read NetCDF forcing value.\n");
        fprintf(stderr, "  File: %s\n", r.file.c_str());
        fprintf(stderr, "  Var: %s\n", r.var_name.c_str());
        fprintf(stderr, "  Index: time=%zu lat=%zu lon=%zu\n", time_idx, lat_idx, lon_idx);
        fprintf(stderr, "  Station[%d]: lon=%.6f lat=%.6f (grid_lon=%.6f grid_lat=%.6f)\n",
                station_idx, station_lon, station_lat, grid_lon, grid_lat);
        fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
        myexit(ERRFileIO);
    }

    if (!std::isfinite(raw)) {
        fprintf(stderr, "\n  Fatal Error: NetCDF forcing value is not finite.\n");
        fprintf(stderr, "  File: %s\n", r.file.c_str());
        fprintf(stderr, "  Var: %s\n", r.var_name.c_str());
        fprintf(stderr, "  Index: time=%zu lat=%zu lon=%zu\n", time_idx, lat_idx, lon_idx);
        fprintf(stderr, "  Station[%d]: lon=%.6f lat=%.6f (grid_lon=%.6f grid_lat=%.6f)\n",
                station_idx, station_lon, station_lat, grid_lon, grid_lat);
        myexit(ERRDATAIN);
    }

    if (r.has_fill && raw == r.fill) {
        fprintf(stderr, "\n  Fatal Error: NetCDF forcing value is _FillValue.\n");
        fprintf(stderr, "  File: %s\n", r.file.c_str());
        fprintf(stderr, "  Var: %s\n", r.var_name.c_str());
        fprintf(stderr, "  Index: time=%zu lat=%zu lon=%zu\n", time_idx, lat_idx, lon_idx);
        fprintf(stderr, "  Station[%d]: lon=%.6f lat=%.6f (grid_lon=%.6f grid_lat=%.6f)\n",
                station_idx, station_lon, station_lat, grid_lon, grid_lat);
        myexit(ERRDATAIN);
    }
    if (r.has_missing && raw == r.missing) {
        fprintf(stderr, "\n  Fatal Error: NetCDF forcing value is missing_value.\n");
        fprintf(stderr, "  File: %s\n", r.file.c_str());
        fprintf(stderr, "  Var: %s\n", r.var_name.c_str());
        fprintf(stderr, "  Index: time=%zu lat=%zu lon=%zu\n", time_idx, lat_idx, lon_idx);
        fprintf(stderr, "  Station[%d]: lon=%.6f lat=%.6f (grid_lon=%.6f grid_lat=%.6f)\n",
                station_idx, station_lon, station_lat, grid_lon, grid_lat);
        myexit(ERRDATAIN);
    }

    double v = raw;
    if (r.has_scale) {
        v *= r.scale;
    }
    if (r.has_offset) {
        v += r.offset;
    }
    return v;
}

static std::map<std::string, std::string> readKvCfg(const char *path)
{
    std::ifstream f(path);
    if (!f.is_open()) {
        fprintf(stderr, "\n  Fatal Error: Failed to open forcing cfg.\n");
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
            fprintf(stderr, "\n  Fatal Error: Invalid KEY VALUE forcing cfg line (missing value).\n");
            fprintf(stderr, "  File: %s\n", path);
            fprintf(stderr, "  Line: %ld\n", lineNo);
            fprintf(stderr, "  Content: %s\n", t.c_str());
            myexit(ERRFileIO);
        }
        kv[toUpper(k)] = v;
    }
    return kv;
}

} // namespace

struct NetcdfForcingProvider::Impl {
    std::vector<ForcingStationMeta> stations;
    long forc_start_yyyymmdd = 0;
    double sim_start_min = 0.0;
    double sim_end_min = 0.0;

    std::string cfg_path;

    // Parsed cfg
    std::string product;
    std::string data_root_abs;

    std::string layout_file_pattern;
    std::map<std::string, std::string> layout_var_dir; // key: PREC/TEMP/...

    std::string dim_time;
    std::string dim_lat;
    std::string dim_lon;

    std::string time_var;
    std::string lat_var;
    std::string lon_var;

    std::map<std::string, std::string> nc_var; // key: PREC/TEMP/...

    std::string radiation_kind;
    std::string cmfd_precip_units; // AUTO/KG_M2_S/MM_HR/MM_DAY
    bool layout_year_subdir = false; // ERA5: DATA_ROOT/<yyyy>/... when true (fallback supported)

    // Coordinates (global, derived from first month file)
    std::vector<double> grid_lat;
    std::vector<double> grid_lon;
    bool grid_lon_is_0360 = false;

    // Station -> grid index
    std::vector<size_t> station_lat_idx;
    std::vector<size_t> station_lon_idx;
    std::vector<double> station_grid_lat;
    std::vector<double> station_grid_lon;

    // Global time axis (minutes since forc_start_yyyymmdd)
    struct TimeMapItem {
        int file_idx = -1;         // CMFD2: month index; ERA5: day index
        size_t local_time_idx = 0; // time index within that file
    };
    std::vector<double> time_min;
    std::vector<TimeMapItem> time_map;

    int now_time_idx = 0;
    int loaded_time_idx = -1;

    // Cache arrays (SHUD forcing contract units)
    std::vector<double> prcp_mm_day;
    std::vector<double> temp_c;
    std::vector<double> rh_1;
    std::vector<double> wind_ms;
    std::vector<double> rn_wm2;

    // CMFD2 file set per month
    struct CmfdMonthFiles {
        std::string yyyymm;
        std::string f_prec;
        std::string f_temp;
        std::string f_shum;
        std::string f_srad;
        std::string f_wind;
        std::string f_pres;
    };
    std::vector<CmfdMonthFiles> cmfd_months;

    // ERA5 file set per day (single NetCDF per day)
    struct Era5DayFile {
        std::string yyyymmdd;
        std::string file;
    };
    std::vector<Era5DayFile> era5_days;

    int open_month_idx = -1;
    NcVarPointReader r_prec;
    NcVarPointReader r_temp;
    NcVarPointReader r_shum;
    NcVarPointReader r_srad;
    NcVarPointReader r_wind;
    NcVarPointReader r_pres;

    ~Impl()
    {
        closeMonthReaders();
        closeEra5Readers();
    }

    void closeMonthReaders()
    {
        r_prec.close();
        r_temp.close();
        r_shum.close();
        r_srad.close();
        r_wind.close();
        r_pres.close();
        open_month_idx = -1;
    }

    // ERA5 readers (single file per day)
    int open_day_idx = -1;
    NcVarPointReader r_tp;
    NcVarPointReader r_t2m;
    NcVarPointReader r_d2m;
    NcVarPointReader r_u10;
    NcVarPointReader r_v10;
    NcVarPointReader r_ssr;
    NcVarPointReader r_sp;

    void closeEra5Readers()
    {
        r_tp.close();
        r_t2m.close();
        r_d2m.close();
        r_u10.close();
        r_v10.close();
        r_ssr.close();
        r_sp.close();
        open_day_idx = -1;
    }

    void parseCfg()
    {
        const std::map<std::string, std::string> kv = readKvCfg(cfg_path.c_str());

        auto must = [&](const char *k) -> std::string {
            auto it = kv.find(toUpper(k));
            if (it == kv.end() || it->second.empty()) {
                fprintf(stderr, "\n  Fatal Error: Missing required forcing cfg key.\n");
                fprintf(stderr, "  File: %s\n", cfg_path.c_str());
                fprintf(stderr, "  Key: %s\n", k);
                myexit(ERRFileIO);
            }
            return it->second;
        };
        auto get = [&](const char *k, const char *defval) -> std::string {
            auto it = kv.find(toUpper(k));
            if (it == kv.end() || it->second.empty()) {
                return defval ? std::string(defval) : std::string();
            }
            return it->second;
        };

        product = toUpper(must("PRODUCT"));
        radiation_kind = toUpper(get("RADIATION_KIND", ""));

        // Resolve DATA_ROOT relative to run_dir = <forcing_cfg_dir>/../..
        const std::string data_root_cfg = must("DATA_ROOT");
        if (isAbsPath(data_root_cfg)) {
            data_root_abs = data_root_cfg;
        } else {
            const std::string cfg_dir = dirnameOf(cfg_path);
            const std::string input_dir = dirnameOf(cfg_dir);    // .../input
            const std::string run_dir = dirnameOf(input_dir);    // run root
            data_root_abs = joinPath(run_dir, data_root_cfg);
        }

        // Common names (rendered by meta-runner today)
        layout_file_pattern = get("LAYOUT_FILE_PATTERN", "");
        dim_time = get("NC_DIM_TIME", "time");
        dim_lat = get("NC_DIM_LAT", "lat");
        dim_lon = get("NC_DIM_LON", "lon");

        time_var = get("TIME_VAR", dim_time.c_str());
        lat_var = get("LAT_VAR", dim_lat.c_str());
        lon_var = get("LON_VAR", dim_lon.c_str());

        cmfd_precip_units = toUpper(get("CMFD_PRECIP_UNITS", "AUTO"));
        {
            std::string ys = get("LAYOUT_YEAR_SUBDIR", "");
            if (ys.empty()) {
                ys = get("ERA5_YEAR_SUBDIR", "");
            }
            ys = toUpper(trim(ys));
            layout_year_subdir = (ys == "1" || ys == "TRUE" || ys == "YES");
        }

        // Layout var dirs
        for (const auto &it : kv) {
            const std::string &k = it.first;
            const std::string &v = it.second;
            const std::string prefix = "LAYOUT_VAR_DIR_";
            if (k.size() > prefix.size() && k.compare(0, prefix.size(), prefix) == 0) {
                layout_var_dir[k.substr(prefix.size())] = v;
            }
        }

        // NetCDF variable names
        for (const auto &it : kv) {
            const std::string &k = it.first;
            const std::string &v = it.second;
            const std::string prefix = "NC_VAR_";
            if (k.size() > prefix.size() && k.compare(0, prefix.size(), prefix) == 0) {
                nc_var[k.substr(prefix.size())] = v;
            }
        }

        // Backward-compatible aliases (from docs/SPEC_阶段A_...; allow overriding)
        if (layout_file_pattern.empty()) {
            if (product == "CMFD2") {
                layout_file_pattern = get("CMFD_FILE_PATTERN", "");
            } else if (product == "ERA5") {
                layout_file_pattern = get("ERA5_FILE_PATTERN", "");
            }
        }
    }

    void requireCmfdKey(const char *k, const std::string &v) const
    {
        if (v.empty()) {
            fprintf(stderr, "\n  Fatal Error: Missing required CMFD2 forcing cfg key.\n");
            fprintf(stderr, "  File: %s\n", cfg_path.c_str());
            fprintf(stderr, "  Key: %s\n", k);
            myexit(ERRFileIO);
        }
    }

    void initCmfdMonths()
    {
        // Determine required months from simulation interval.
        int base_y = 0;
        unsigned base_m = 0;
        unsigned base_d = 0;
        if (!parseYYYYMMDD(forc_start_yyyymmdd, base_y, base_m, base_d)) {
            fprintf(stderr, "\n  Fatal Error: Invalid ForcStartTime for NetCDF forcing: %ld\n", forc_start_yyyymmdd);
            myexit(ERRDATAIN);
        }
        const long long base_days = daysFromCivil(base_y, base_m, base_d);

        const long long start_days = base_days + (long long)std::floor(sim_start_min / 1440.0);
        // Simulation END is an exclusive bound for forcing file discovery: the model never
        // queries forcing at t == sim_end_min, only for times in [sim_start_min, sim_end_min).
        // Avoid requesting the month of the END boundary when it lands exactly on a month/day boundary.
        double end_min_excl = sim_end_min;
        if (sim_end_min > sim_start_min + 1e-12) {
            end_min_excl = std::nextafter(sim_end_min, -std::numeric_limits<double>::infinity());
        }
        const long long end_days = base_days + (long long)std::floor(end_min_excl / 1440.0);
        int y0 = 0, y1 = 0;
        unsigned m0 = 0, m1 = 0;
        unsigned dtmp = 0;
        civilFromDays(start_days, y0, m0, dtmp);
        civilFromDays(end_days, y1, m1, dtmp);

        std::vector<std::string> months;
        int y = y0;
        unsigned m = m0;
        while (y < y1 || (y == y1 && m <= m1)) {
            months.push_back(zfillInt(y, 4) + zfillInt((int)m, 2));
            m++;
            if (m > 12) {
                m = 1;
                y++;
            }
        }

        // Required keys
        requireCmfdKey("LAYOUT_FILE_PATTERN", layout_file_pattern);
        const auto get_dir = [&](const char *key, const char *fallback) -> std::string {
            const std::string k_up = toUpper(std::string(key));
            auto it = layout_var_dir.find(k_up);
            if (it != layout_var_dir.end()) {
                return it->second;
            }
            // backward-compat keys
            // (not currently rendered by runner, but allowed for manual override)
            // fallback should be like "Prec"
            (void)fallback;
            return "";
        };
        const auto get_var = [&](const char *key) -> std::string {
            const std::string k_up = toUpper(std::string(key));
            auto it = nc_var.find(k_up);
            if (it != nc_var.end()) {
                return it->second;
            }
            return "";
        };

        const std::string dir_prec = get_dir("PREC", "Prec");
        const std::string dir_temp = get_dir("TEMP", "Temp");
        const std::string dir_shum = get_dir("SHUM", "SHum");
        const std::string dir_srad = get_dir("SRAD", "SRad");
        const std::string dir_wind = get_dir("WIND", "Wind");
        const std::string dir_pres = get_dir("PRES", "Pres");

        requireCmfdKey("LAYOUT_VAR_DIR_PREC", dir_prec);
        requireCmfdKey("LAYOUT_VAR_DIR_TEMP", dir_temp);
        requireCmfdKey("LAYOUT_VAR_DIR_SHUM", dir_shum);
        requireCmfdKey("LAYOUT_VAR_DIR_SRAD", dir_srad);
        requireCmfdKey("LAYOUT_VAR_DIR_WIND", dir_wind);
        requireCmfdKey("LAYOUT_VAR_DIR_PRES", dir_pres);

        const std::string v_prec = get_var("PREC");
        const std::string v_temp = get_var("TEMP");
        const std::string v_shum = get_var("SHUM");
        const std::string v_srad = get_var("SRAD");
        const std::string v_wind = get_var("WIND");
        const std::string v_pres = get_var("PRES");
        requireCmfdKey("NC_VAR_PREC", v_prec);
        requireCmfdKey("NC_VAR_TEMP", v_temp);
        requireCmfdKey("NC_VAR_SHUM", v_shum);
        requireCmfdKey("NC_VAR_SRAD", v_srad);
        requireCmfdKey("NC_VAR_WIND", v_wind);
        requireCmfdKey("NC_VAR_PRES", v_pres);

        cmfd_months.clear();
        cmfd_months.reserve(months.size());
        for (const std::string &yyyymm : months) {
            auto resolve = [&](const std::string &dir, const std::string &var_lower) -> std::string {
                std::string pat = layout_file_pattern;
                replaceAll(pat, "{var_lower}", toLower(var_lower));
                replaceAll(pat, "{yyyymm}", yyyymm);
                const std::string globpat = joinPath(joinPath(data_root_abs, dir), pat);
                return resolveSingleGlob(globpat);
            };
            CmfdMonthFiles mfiles;
            mfiles.yyyymm = yyyymm;
            mfiles.f_prec = resolve(dir_prec, v_prec);
            mfiles.f_temp = resolve(dir_temp, v_temp);
            mfiles.f_shum = resolve(dir_shum, v_shum);
            mfiles.f_srad = resolve(dir_srad, v_srad);
            mfiles.f_wind = resolve(dir_wind, v_wind);
            mfiles.f_pres = resolve(dir_pres, v_pres);
            cmfd_months.push_back(mfiles);
        }
    }

    std::vector<double> readCoord1d(const std::string &file, const std::string &var_name)
    {
        int ncid = -1;
        int rc = nc_open(file.c_str(), NC_NOWRITE, &ncid);
        if (rc != NC_NOERR) {
            fprintf(stderr, "\n  Fatal Error: Failed to open NetCDF file for coordinate read.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
            myexit(ERRFileIO);
        }
        int varid = -1;
        rc = nc_inq_varid(ncid, var_name.c_str(), &varid);
        if (rc != NC_NOERR) {
            fprintf(stderr, "\n  Fatal Error: Coordinate variable not found.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Var: %s\n", var_name.c_str());
            fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
            nc_close(ncid);
            myexit(ERRFileIO);
        }
        char name_buf[NC_MAX_NAME + 1];
        nc_type xtype;
        int ndims = 0;
        int dimids[NC_MAX_VAR_DIMS];
        int natts = 0;
        rc = nc_inq_var(ncid, varid, name_buf, &xtype, &ndims, dimids, &natts);
        if (rc != NC_NOERR) {
            fprintf(stderr, "\n  Fatal Error: nc_inq_var failed for coordinate.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Var: %s\n", var_name.c_str());
            fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
            nc_close(ncid);
            myexit(ERRFileIO);
        }
        if (ndims != 1) {
            fprintf(stderr, "\n  Fatal Error: Coordinate variable must be 1-D.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Var: %s\n", var_name.c_str());
            fprintf(stderr, "  ndims: %d\n", ndims);
            nc_close(ncid);
            myexit(ERRFileIO);
        }
        size_t dimlen = 0;
        rc = nc_inq_dimlen(ncid, dimids[0], &dimlen);
        if (rc != NC_NOERR) {
            fprintf(stderr, "\n  Fatal Error: nc_inq_dimlen failed for coordinate.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Var: %s\n", var_name.c_str());
            fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
            nc_close(ncid);
            myexit(ERRFileIO);
        }
        std::vector<double> out(dimlen, 0.0);
        rc = nc_get_var_double(ncid, varid, out.data());
        if (rc != NC_NOERR) {
            fprintf(stderr, "\n  Fatal Error: Failed to read coordinate variable.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Var: %s\n", var_name.c_str());
            fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
            nc_close(ncid);
            myexit(ERRFileIO);
        }
        nc_close(ncid);
        return out;
    }

    std::vector<double> readTimeAxisMin(const std::string &file)
    {
        int ncid = -1;
        int rc = nc_open(file.c_str(), NC_NOWRITE, &ncid);
        if (rc != NC_NOERR) {
            fprintf(stderr, "\n  Fatal Error: Failed to open NetCDF file for time read.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
            myexit(ERRFileIO);
        }
        int varid = -1;
        rc = nc_inq_varid(ncid, time_var.c_str(), &varid);
        if (rc != NC_NOERR) {
            fprintf(stderr, "\n  Fatal Error: Time variable not found.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  TIME_VAR: %s\n", time_var.c_str());
            fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
            nc_close(ncid);
            myexit(ERRFileIO);
        }

        std::string units = ncGetAttText(ncid, varid, "units");
        if (units.empty()) {
            fprintf(stderr, "\n  Fatal Error: NetCDF time variable missing units attribute.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Var: %s\n", time_var.c_str());
            nc_close(ncid);
            myexit(ERRFileIO);
        }

        // Determine time length from its dim (assume 1-D)
        char name_buf[NC_MAX_NAME + 1];
        nc_type xtype;
        int ndims = 0;
        int dimids[NC_MAX_VAR_DIMS];
        int natts = 0;
        rc = nc_inq_var(ncid, varid, name_buf, &xtype, &ndims, dimids, &natts);
        if (rc != NC_NOERR || ndims != 1) {
            fprintf(stderr, "\n  Fatal Error: NetCDF time variable must be 1-D.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Var: %s\n", time_var.c_str());
            fprintf(stderr, "  ndims: %d\n", ndims);
            fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
            nc_close(ncid);
            myexit(ERRFileIO);
        }
        size_t len = 0;
        rc = nc_inq_dimlen(ncid, dimids[0], &len);
        if (rc != NC_NOERR || len == 0) {
            fprintf(stderr, "\n  Fatal Error: NetCDF time dimension length invalid.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Var: %s\n", time_var.c_str());
            fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
            nc_close(ncid);
            myexit(ERRFileIO);
        }

        std::vector<double> tv(len, 0.0);
        rc = nc_get_var_double(ncid, varid, tv.data());
        if (rc != NC_NOERR) {
            fprintf(stderr, "\n  Fatal Error: Failed to read NetCDF time variable values.\n");
            fprintf(stderr, "  File: %s\n", file.c_str());
            fprintf(stderr, "  Var: %s\n", time_var.c_str());
            fprintf(stderr, "  Error: %s\n", ncStrError(rc).c_str());
            nc_close(ncid);
            myexit(ERRFileIO);
        }
        nc_close(ncid);

        // Convert to minutes relative to ForcStartTime
        int y0 = 0;
        unsigned m0 = 0;
        unsigned d0 = 0;
        if (!parseYYYYMMDD(forc_start_yyyymmdd, y0, m0, d0)) {
            fprintf(stderr, "\n  Fatal Error: Invalid ForcStartTime for time axis conversion: %ld\n", forc_start_yyyymmdd);
            myexit(ERRDATAIN);
        }
        const long long forc_base_min = daysFromCivil(y0, m0, d0) * 1440LL;
        const UnitsSince us = parseUnitsSince(units);

        std::vector<double> out;
        out.reserve(tv.size());
        double prev = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < tv.size(); i++) {
            const double abs_min = (double)us.base_minutes_since_epoch + tv[i] * us.factor_to_minutes;
            const double t_min = abs_min - (double)forc_base_min;
            if (i > 0 && t_min + 1e-9 < prev) {
                fprintf(stderr, "\n  Fatal Error: NetCDF time axis is not monotonic non-decreasing.\n");
                fprintf(stderr, "  File: %s\n", file.c_str());
                fprintf(stderr, "  Var: %s\n", time_var.c_str());
                fprintf(stderr, "  Index[%zu]=%.6f min, Index[%zu]=%.6f min\n", i - 1, prev, i, t_min);
                myexit(ERRDATAIN);
            }
            out.push_back(t_min);
            prev = t_min;
        }
        return out;
    }

    void discoverGridAndStationsCmfd()
    {
        if (cmfd_months.empty()) {
            fprintf(stderr, "\n  Fatal Error: CMFD2 month file list is empty.\n");
            myexit(ERRFileIO);
        }
        const std::string &grid_file = cmfd_months[0].f_prec;
        grid_lat = readCoord1d(grid_file, lat_var);
        grid_lon = readCoord1d(grid_file, lon_var);
        if (grid_lat.empty() || grid_lon.empty()) {
            fprintf(stderr, "\n  Fatal Error: Failed to read grid coordinate arrays.\n");
            fprintf(stderr, "  File: %s\n", grid_file.c_str());
            fprintf(stderr, "  LAT_VAR=%s (n=%zu)\n", lat_var.c_str(), grid_lat.size());
            fprintf(stderr, "  LON_VAR=%s (n=%zu)\n", lon_var.c_str(), grid_lon.size());
            myexit(ERRFileIO);
        }

        double lon_min = grid_lon[0];
        double lon_max = grid_lon[0];
        for (double v : grid_lon) {
            lon_min = std::min(lon_min, v);
            lon_max = std::max(lon_max, v);
        }
        grid_lon_is_0360 = (lon_min >= 0.0 && lon_max > 180.0);

        const size_t nst = stations.size();
        station_lat_idx.assign(nst, 0);
        station_lon_idx.assign(nst, 0);
        station_grid_lat.assign(nst, NA_VALUE);
        station_grid_lon.assign(nst, NA_VALUE);

        for (size_t i = 0; i < nst; i++) {
            double slon = stations[i].lon_deg;
            double slat = stations[i].lat_deg;
            if (grid_lon_is_0360) {
                if (slon < 0.0) {
                    slon += 360.0;
                }
                while (slon >= 360.0) {
                    slon -= 360.0;
                }
            }

            size_t best_j = 0;
            double best_d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < grid_lon.size(); j++) {
                const double d = std::fabs(grid_lon[j] - slon);
                if (d < best_d) {
                    best_d = d;
                    best_j = j;
                }
            }
            station_lon_idx[i] = best_j;
            station_grid_lon[i] = grid_lon[best_j];

            size_t best_k = 0;
            double best_dlat = std::numeric_limits<double>::infinity();
            for (size_t k = 0; k < grid_lat.size(); k++) {
                const double d = std::fabs(grid_lat[k] - slat);
                if (d < best_dlat) {
                    best_dlat = d;
                    best_k = k;
                }
            }
            station_lat_idx[i] = best_k;
            station_grid_lat[i] = grid_lat[best_k];
        }

        const size_t nlog = std::min<size_t>(stations.size(), 3);
        fprintf(stdout, "\tNetCDF forcing: PRODUCT=%s, stations=%zu, grid=(lat=%zu, lon=%zu)\n",
                product.c_str(), stations.size(), grid_lat.size(), grid_lon.size());
        for (size_t i = 0; i < nlog; i++) {
            fprintf(stdout,
                    "\tNetCDF forcing map[%zu]: station(lon=%.6f, lat=%.6f) -> grid(idx_lat=%zu idx_lon=%zu; lon=%.6f lat=%.6f)\n",
                    i,
                    stations[i].lon_deg,
                    stations[i].lat_deg,
                    station_lat_idx[i],
                    station_lon_idx[i],
                    station_grid_lon[i],
                    station_grid_lat[i]);
        }
    }

    void buildGlobalTimeAxisCmfd()
    {
        time_min.clear();
        time_map.clear();
        for (size_t mi = 0; mi < cmfd_months.size(); mi++) {
            const std::vector<double> t = readTimeAxisMin(cmfd_months[mi].f_prec);
            for (size_t k = 0; k < t.size(); k++) {
                if (!time_min.empty() && t[k] + 1e-9 < time_min.back()) {
                    fprintf(stderr, "\n  Fatal Error: NetCDF time axis across files is not monotonic.\n");
                    fprintf(stderr, "  Previous max t_min: %.6f\n", time_min.back());
                    fprintf(stderr, "  Current  t_min: %.6f\n", t[k]);
                    fprintf(stderr, "  File: %s\n", cmfd_months[mi].f_prec.c_str());
                    myexit(ERRDATAIN);
                }
                time_min.push_back(t[k]);
                TimeMapItem item;
                item.file_idx = (int)mi;
                item.local_time_idx = k;
                time_map.push_back(item);
            }
        }
        if (time_min.empty()) {
            fprintf(stderr, "\n  Fatal Error: NetCDF time axis is empty.\n");
            myexit(ERRFileIO);
        }
        fprintf(stdout, "\tNetCDF forcing time coverage: [%.3f, %.3f] min ([%.6f, %.6f] day)\n",
                time_min.front(), time_min.back(), time_min.front() / 1440.0, time_min.back() / 1440.0);
    }

    void openCmfdMonth(int month_idx)
    {
        if (month_idx < 0 || month_idx >= (int)cmfd_months.size()) {
            fprintf(stderr, "\n  Fatal Error: Invalid month index for NetCDF forcing: %d\n", month_idx);
            myexit(ERRDATAIN);
        }
        if (open_month_idx == month_idx) {
            return;
        }
        closeMonthReaders();

        const CmfdMonthFiles &m = cmfd_months[(size_t)month_idx];
        r_prec = openPointReader(m.f_prec, nc_var["PREC"], dim_time, dim_lat, dim_lon);
        r_temp = openPointReader(m.f_temp, nc_var["TEMP"], dim_time, dim_lat, dim_lon);
        r_shum = openPointReader(m.f_shum, nc_var["SHUM"], dim_time, dim_lat, dim_lon);
        r_srad = openPointReader(m.f_srad, nc_var["SRAD"], dim_time, dim_lat, dim_lon);
        r_wind = openPointReader(m.f_wind, nc_var["WIND"], dim_time, dim_lat, dim_lon);
        r_pres = openPointReader(m.f_pres, nc_var["PRES"], dim_time, dim_lat, dim_lon);
        open_month_idx = month_idx;
    }

    enum CmfdPrecipUnits { CMFD_AUTO, CMFD_KG_M2_S, CMFD_MM_HR, CMFD_MM_DAY };

    CmfdPrecipUnits determineCmfdPrecipUnits() const
    {
        if (cmfd_precip_units == "KG_M2_S") return CMFD_KG_M2_S;
        if (cmfd_precip_units == "MM_HR" || cmfd_precip_units == "MM/HR" || cmfd_precip_units == "MM_H-1") return CMFD_MM_HR;
        if (cmfd_precip_units == "MM_DAY" || cmfd_precip_units == "MM/DAY" || cmfd_precip_units == "MM_D-1") return CMFD_MM_DAY;
        if (cmfd_precip_units != "AUTO") {
            fprintf(stderr, "\n  Fatal Error: Invalid CMFD_PRECIP_UNITS value: %s\n", cmfd_precip_units.c_str());
            fprintf(stderr, "  Valid: AUTO | KG_M2_S | MM_HR | MM_DAY\n");
            myexit(ERRFileIO);
        }

        // Auto-detect from NetCDF metadata in the current open month.
        if (open_month_idx < 0) {
            return CMFD_AUTO;
        }
        std::string units = ncGetAttText(r_prec.ncid, r_prec.varid, "units");
        units = toLower(trim(units));
        if (units.find("kg") != std::string::npos &&
            (units.find("m-2") != std::string::npos || units.find("m**-2") != std::string::npos) &&
            (units.find("s-1") != std::string::npos || units.find("s**-1") != std::string::npos)) {
            return CMFD_KG_M2_S;
        }
        if (units.find("mm") != std::string::npos &&
            (units.find("hr") != std::string::npos || units.find("h-1") != std::string::npos || units.find("h**-1") != std::string::npos)) {
            return CMFD_MM_HR;
        }
        if (units.find("mm") != std::string::npos &&
            (units.find("day") != std::string::npos || units.find("d-1") != std::string::npos || units.find("d**-1") != std::string::npos)) {
            return CMFD_MM_DAY;
        }

        fprintf(stderr, "\n  Fatal Error: Failed to auto-detect CMFD2 precip units from NetCDF metadata.\n");
        fprintf(stderr, "  File: %s\n", r_prec.file.c_str());
        fprintf(stderr, "  Var: %s\n", r_prec.var_name.c_str());
        fprintf(stderr, "  units: %s\n", units.c_str());
        fprintf(stderr, "  Fix: set CMFD_PRECIP_UNITS in %s (AUTO|KG_M2_S|MM_HR|MM_DAY).\n", cfg_path.c_str());
        myexit(ERRFileIO);
        return CMFD_AUTO;
    }

    void loadCacheForTimeIndexCmfd(int t_idx)
    {
        if (t_idx < 0 || t_idx >= (int)time_map.size()) {
            fprintf(stderr, "\n  Fatal Error: Invalid time index for NetCDF forcing cache load: %d\n", t_idx);
            myexit(ERRDATAIN);
        }
        const TimeMapItem &tm = time_map[(size_t)t_idx];
        openCmfdMonth(tm.file_idx);
        const CmfdPrecipUnits pu = determineCmfdPrecipUnits();

        const size_t nst = stations.size();
        const size_t t_local = tm.local_time_idx;
        for (size_t i = 0; i < nst; i++) {
            const size_t ilat = station_lat_idx[i];
            const size_t ilon = station_lon_idx[i];

            const double slon = stations[i].lon_deg;
            const double slat = stations[i].lat_deg;
            const double glon = station_grid_lon[i];
            const double glat = station_grid_lat[i];

            const double prec_raw = readPoint(r_prec, t_local, ilat, ilon, (int)i, slon, slat, glon, glat);
            const double temp_k = readPoint(r_temp, t_local, ilat, ilon, (int)i, slon, slat, glon, glat);
            const double shum = readPoint(r_shum, t_local, ilat, ilon, (int)i, slon, slat, glon, glat);
            const double srad = readPoint(r_srad, t_local, ilat, ilon, (int)i, slon, slat, glon, glat);
            const double wind = readPoint(r_wind, t_local, ilat, ilon, (int)i, slon, slat, glon, glat);
            const double pres = readPoint(r_pres, t_local, ilat, ilon, (int)i, slon, slat, glon, glat);

            double prec_mm_day = 0.0;
            if (pu == CMFD_KG_M2_S) {
                prec_mm_day = prec_raw * 86400.0;
            } else if (pu == CMFD_MM_HR) {
                prec_mm_day = prec_raw * 24.0;
            } else if (pu == CMFD_MM_DAY) {
                prec_mm_day = prec_raw;
            } else {
                // Should never happen; CMFD_AUTO must resolve.
                fprintf(stderr, "\n  Fatal Error: CMFD precip units unresolved.\n");
                myexit(ERRFileIO);
            }
            if (!std::isfinite(prec_mm_day)) {
                prec_mm_day = 0.0;
            }
            if (prec_mm_day < 0.0) {
                prec_mm_day = 0.0;
            }
            // Match AutoSHUD baseline forcing CSV quantization (4 decimals) before thresholding.
            prec_mm_day = std::nearbyint(prec_mm_day * 10000.0) / 10000.0;
            const double min_mm_day = 0.0001; // match AutoSHUD behavior
            if (prec_mm_day < min_mm_day) {
                prec_mm_day = 0.0;
            }

            double temp_c_val = temp_k - 273.15;
            if (!std::isfinite(temp_c_val)) {
                temp_c_val = 0.0;
            }
            // Match AutoSHUD baseline forcing CSV quantization (2 decimals).
            temp_c_val = std::nearbyint(temp_c_val * 100.0) / 100.0;

            double rh_percent = 0.263 * pres * shum /
                                std::exp(17.67 * (temp_k - 273.15) / (temp_k - 29.65));
            if (!std::isfinite(rh_percent)) {
                rh_percent = 0.0;
            }
            if (rh_percent > 100.0) {
                rh_percent = 100.0;
            }
            if (rh_percent < 0.0) {
                rh_percent = 0.0;
            }
            double rh_val = rh_percent / 100.0;
            // Match AutoSHUD baseline forcing CSV quantization (4 decimals).
            rh_val = std::nearbyint(rh_val * 10000.0) / 10000.0;
            if (rh_val > 1.0) {
                rh_val = 1.0;
            }
            if (rh_val < 0.0) {
                rh_val = 0.0;
            }

            double wind_val = std::fabs(wind);
            if (!std::isfinite(wind_val)) {
                wind_val = 0.0;
            }
            // Match AutoSHUD baseline forcing CSV quantization (2 decimals) + min clamp.
            wind_val = std::nearbyint(wind_val * 100.0) / 100.0;
            const double min_wind_ms = 0.05;
            if (wind_val < min_wind_ms) {
                wind_val = min_wind_ms;
            }

            double rn_val = srad;
            if (!std::isfinite(rn_val)) {
                rn_val = 0.0;
            }
            if (rn_val < 0.0) {
                rn_val = 0.0;
            }
            // Match AutoSHUD baseline forcing CSV quantization (integer).
            rn_val = std::nearbyint(rn_val);

            prcp_mm_day[i] = prec_mm_day;
            temp_c[i] = temp_c_val;
            rh_1[i] = rh_val;
            wind_ms[i] = wind_val;
            rn_wm2[i] = rn_val;
        }
        loaded_time_idx = t_idx;
    }

    void requireEra5Key(const char *k, const std::string &v) const
    {
        if (v.empty()) {
            fprintf(stderr, "\n  Fatal Error: Missing required ERA5 forcing cfg key.\n");
            fprintf(stderr, "  File: %s\n", cfg_path.c_str());
            fprintf(stderr, "  Key: %s\n", k);
            myexit(ERRFileIO);
        }
    }

    void initEra5Days()
    {
        int base_y = 0;
        unsigned base_m = 0;
        unsigned base_d = 0;
        if (!parseYYYYMMDD(forc_start_yyyymmdd, base_y, base_m, base_d)) {
            fprintf(stderr, "\n  Fatal Error: Invalid ForcStartTime for NetCDF forcing: %ld\n", forc_start_yyyymmdd);
            myexit(ERRDATAIN);
        }
        const long long base_days = daysFromCivil(base_y, base_m, base_d);
        const long long start_days = base_days + (long long)std::floor(sim_start_min / 1440.0);
        const long long end_days = base_days + (long long)std::floor(sim_end_min / 1440.0);

        requireEra5Key("LAYOUT_FILE_PATTERN", layout_file_pattern);

        auto tryResolveSingle = [&](const std::string &pattern, std::string &out) -> bool {
            glob_t g;
            memset(&g, 0, sizeof(g));
            const int rc = glob(pattern.c_str(), 0, nullptr, &g);
            if (rc != 0 || g.gl_pathc == 0) {
                globfree(&g);
                return false;
            }
            std::vector<std::string> matches;
            matches.reserve(g.gl_pathc);
            for (size_t i = 0; i < g.gl_pathc; i++) {
                matches.emplace_back(g.gl_pathv[i]);
            }
            globfree(&g);
            std::sort(matches.begin(), matches.end());
            if (matches.size() != 1) {
                fprintf(stderr, "\n  Fatal Error: NetCDF forcing file glob is ambiguous (%zu matches).\n", matches.size());
                fprintf(stderr, "  Glob: %s\n", pattern.c_str());
                myexit(ERRFileIO);
            }
            out = matches[0];
            return true;
        };

        era5_days.clear();
        era5_days.reserve((size_t)std::max<long long>(0LL, end_days - start_days + 1));
        for (long long z = start_days; z <= end_days; z++) {
            int y = 0;
            unsigned m = 0;
            unsigned d = 0;
            civilFromDays(z, y, m, d);
            const std::string yyyy = zfillInt(y, 4);
            const std::string yyyymmdd = yyyy + zfillInt((int)m, 2) + zfillInt((int)d, 2);

            std::string pat = layout_file_pattern;
            replaceAll(pat, "{yyyymmdd}", yyyymmdd);

            std::string resolved;
            if (layout_year_subdir) {
                const std::string p1 = joinPath(joinPath(data_root_abs, yyyy), pat);
                if (!tryResolveSingle(p1, resolved)) {
                    // Fallback: some datasets do not use the year subdir in practice.
                    const std::string p2 = joinPath(data_root_abs, pat);
                    if (!tryResolveSingle(p2, resolved)) {
                        fprintf(stderr, "\n  Fatal Error: ERA5 NetCDF file not found.\n");
                        fprintf(stderr, "  Tried: %s\n", p1.c_str());
                        fprintf(stderr, "  Tried: %s\n", p2.c_str());
                        myexit(ERRFileIO);
                    }
                }
            } else {
                const std::string p = joinPath(data_root_abs, pat);
                resolved = resolveSingleGlob(p);
            }

            Era5DayFile df;
            df.yyyymmdd = yyyymmdd;
            df.file = resolved;
            era5_days.push_back(df);
        }

        // Ensure required variable mappings exist (sp optional)
        requireEra5Key("NC_VAR_TP", nc_var["TP"]);
        requireEra5Key("NC_VAR_T2M", nc_var["T2M"]);
        requireEra5Key("NC_VAR_D2M", nc_var["D2M"]);
        requireEra5Key("NC_VAR_U10", nc_var["U10"]);
        requireEra5Key("NC_VAR_V10", nc_var["V10"]);
        requireEra5Key("NC_VAR_SSR", nc_var["SSR"]);
    }

    void discoverGridAndStationsEra5()
    {
        if (era5_days.empty()) {
            fprintf(stderr, "\n  Fatal Error: ERA5 day file list is empty.\n");
            myexit(ERRFileIO);
        }
        const std::string &grid_file = era5_days[0].file;
        grid_lat = readCoord1d(grid_file, lat_var);
        grid_lon = readCoord1d(grid_file, lon_var);
        if (grid_lat.empty() || grid_lon.empty()) {
            fprintf(stderr, "\n  Fatal Error: Failed to read ERA5 grid coordinate arrays.\n");
            fprintf(stderr, "  File: %s\n", grid_file.c_str());
            fprintf(stderr, "  LAT_VAR=%s (n=%zu)\n", lat_var.c_str(), grid_lat.size());
            fprintf(stderr, "  LON_VAR=%s (n=%zu)\n", lon_var.c_str(), grid_lon.size());
            myexit(ERRFileIO);
        }

        double lon_min = grid_lon[0];
        double lon_max = grid_lon[0];
        for (double v : grid_lon) {
            lon_min = std::min(lon_min, v);
            lon_max = std::max(lon_max, v);
        }
        grid_lon_is_0360 = (lon_min >= 0.0 && lon_max > 180.0);

        const size_t nst = stations.size();
        station_lat_idx.assign(nst, 0);
        station_lon_idx.assign(nst, 0);
        station_grid_lat.assign(nst, NA_VALUE);
        station_grid_lon.assign(nst, NA_VALUE);

        for (size_t i = 0; i < nst; i++) {
            double slon = stations[i].lon_deg;
            double slat = stations[i].lat_deg;
            if (grid_lon_is_0360) {
                if (slon < 0.0) {
                    slon += 360.0;
                }
                while (slon >= 360.0) {
                    slon -= 360.0;
                }
            }

            size_t best_j = 0;
            double best_d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < grid_lon.size(); j++) {
                const double d = std::fabs(grid_lon[j] - slon);
                if (d < best_d) {
                    best_d = d;
                    best_j = j;
                }
            }
            station_lon_idx[i] = best_j;
            station_grid_lon[i] = grid_lon[best_j];

            size_t best_k = 0;
            double best_dlat = std::numeric_limits<double>::infinity();
            for (size_t k = 0; k < grid_lat.size(); k++) {
                const double d = std::fabs(grid_lat[k] - slat);
                if (d < best_dlat) {
                    best_dlat = d;
                    best_k = k;
                }
            }
            station_lat_idx[i] = best_k;
            station_grid_lat[i] = grid_lat[best_k];
        }

        const size_t nlog = std::min<size_t>(stations.size(), 3);
        fprintf(stdout, "\tNetCDF forcing: PRODUCT=%s, stations=%zu, grid=(lat=%zu, lon=%zu)\n",
                product.c_str(), stations.size(), grid_lat.size(), grid_lon.size());
        for (size_t i = 0; i < nlog; i++) {
            fprintf(stdout,
                    "\tNetCDF forcing map[%zu]: station(lon=%.6f, lat=%.6f) -> grid(idx_lat=%zu idx_lon=%zu; lon=%.6f lat=%.6f)\n",
                    i,
                    stations[i].lon_deg,
                    stations[i].lat_deg,
                    station_lat_idx[i],
                    station_lon_idx[i],
                    station_grid_lon[i],
                    station_grid_lat[i]);
        }
    }

    void buildGlobalTimeAxisEra5()
    {
        time_min.clear();
        time_map.clear();
        for (size_t di = 0; di < era5_days.size(); di++) {
            const std::vector<double> t = readTimeAxisMin(era5_days[di].file);
            for (size_t k = 0; k < t.size(); k++) {
                if (!time_min.empty() && t[k] + 1e-9 < time_min.back()) {
                    fprintf(stderr, "\n  Fatal Error: ERA5 time axis across files is not monotonic.\n");
                    fprintf(stderr, "  Previous max t_min: %.6f\n", time_min.back());
                    fprintf(stderr, "  Current  t_min: %.6f\n", t[k]);
                    fprintf(stderr, "  File: %s\n", era5_days[di].file.c_str());
                    myexit(ERRDATAIN);
                }
                time_min.push_back(t[k]);
                TimeMapItem item;
                item.file_idx = (int)di;
                item.local_time_idx = k;
                time_map.push_back(item);
            }
        }
        if (time_min.empty()) {
            fprintf(stderr, "\n  Fatal Error: ERA5 time axis is empty.\n");
            myexit(ERRFileIO);
        }
        fprintf(stdout, "\tNetCDF forcing time coverage: [%.3f, %.3f] min ([%.6f, %.6f] day)\n",
                time_min.front(), time_min.back(), time_min.front() / 1440.0, time_min.back() / 1440.0);
    }

    void openEra5Day(int day_idx)
    {
        if (day_idx < 0 || day_idx >= (int)era5_days.size()) {
            fprintf(stderr, "\n  Fatal Error: Invalid ERA5 day index for NetCDF forcing: %d\n", day_idx);
            myexit(ERRDATAIN);
        }
        if (open_day_idx == day_idx) {
            return;
        }
        closeEra5Readers();

        const std::string &f = era5_days[(size_t)day_idx].file;
        r_tp = openPointReader(f, nc_var["TP"], dim_time, dim_lat, dim_lon);
        r_t2m = openPointReader(f, nc_var["T2M"], dim_time, dim_lat, dim_lon);
        r_d2m = openPointReader(f, nc_var["D2M"], dim_time, dim_lat, dim_lon);
        r_u10 = openPointReader(f, nc_var["U10"], dim_time, dim_lat, dim_lon);
        r_v10 = openPointReader(f, nc_var["V10"], dim_time, dim_lat, dim_lon);
        r_ssr = openPointReader(f, nc_var["SSR"], dim_time, dim_lat, dim_lon);
        if (!nc_var["SP"].empty()) {
            r_sp = openPointReader(f, nc_var["SP"], dim_time, dim_lat, dim_lon);
        }
        open_day_idx = day_idx;
    }

    void loadCacheForTimeIndexEra5(int t_idx)
    {
        if (t_idx < 0 || t_idx >= (int)time_map.size()) {
            fprintf(stderr, "\n  Fatal Error: Invalid time index for ERA5 forcing cache load: %d\n", t_idx);
            myexit(ERRDATAIN);
        }
        const TimeMapItem &tm0 = time_map[(size_t)t_idx];
        openEra5Day(tm0.file_idx);

        const bool has_next = ((size_t)t_idx + 1 < time_map.size());
        TimeMapItem tm1;
        tm1.file_idx = -1;
        tm1.local_time_idx = 0;
        double dt_sec = 3600.0;
        if (has_next) {
            tm1 = time_map[(size_t)t_idx + 1];
            const double dt_min = time_min[(size_t)t_idx + 1] - time_min[(size_t)t_idx];
            dt_sec = dt_min * 60.0;
            if (!(dt_sec > 0.0)) {
                fprintf(stderr, "\n  Fatal Error: ERA5 forcing dt_sec is invalid (<=0).\n");
                fprintf(stderr, "  t_idx=%d, t0=%.6f, t1=%.6f\n", t_idx, time_min[(size_t)t_idx], time_min[(size_t)t_idx + 1]);
                myexit(ERRDATAIN);
            }
        }

        NcVarPointReader tp_next;
        NcVarPointReader ssr_next;
        const bool next_in_other_file = (has_next && tm1.file_idx != tm0.file_idx);
        if (next_in_other_file) {
            const std::string &f1 = era5_days[(size_t)tm1.file_idx].file;
            tp_next = openPointReader(f1, nc_var["TP"], dim_time, dim_lat, dim_lon);
            ssr_next = openPointReader(f1, nc_var["SSR"], dim_time, dim_lat, dim_lon);
        }

        const size_t nst = stations.size();
        const size_t t0_local = tm0.local_time_idx;
        const size_t t1_local = has_next ? tm1.local_time_idx : tm0.local_time_idx;

        for (size_t i = 0; i < nst; i++) {
            const size_t ilat = station_lat_idx[i];
            const size_t ilon = station_lon_idx[i];

            const double slon = stations[i].lon_deg;
            const double slat = stations[i].lat_deg;
            const double glon = station_grid_lon[i];
            const double glat = station_grid_lat[i];

            const double t2m_k = readPoint(r_t2m, t0_local, ilat, ilon, (int)i, slon, slat, glon, glat);
            const double d2m_k = readPoint(r_d2m, t0_local, ilat, ilon, (int)i, slon, slat, glon, glat);
            const double u10 = readPoint(r_u10, t0_local, ilat, ilon, (int)i, slon, slat, glon, glat);
            const double v10 = readPoint(r_v10, t0_local, ilat, ilon, (int)i, slon, slat, glon, glat);

            const double tp0 = readPoint(r_tp, t0_local, ilat, ilon, (int)i, slon, slat, glon, glat);
            const double ssr0 = readPoint(r_ssr, t0_local, ilat, ilon, (int)i, slon, slat, glon, glat);

            double tp1 = tp0;
            double ssr1 = ssr0;
            if (has_next) {
                if (next_in_other_file) {
                    tp1 = readPoint(tp_next, t1_local, ilat, ilon, (int)i, slon, slat, glon, glat);
                    ssr1 = readPoint(ssr_next, t1_local, ilat, ilon, (int)i, slon, slat, glon, glat);
                } else {
                    tp1 = readPoint(r_tp, t1_local, ilat, ilon, (int)i, slon, slat, glon, glat);
                    ssr1 = readPoint(r_ssr, t1_local, ilat, ilon, (int)i, slon, slat, glon, glat);
                }
            }

            // Accumulated -> interval increment
            double tp_inc_m = 0.0;
            double ssr_inc_jm2 = 0.0;
            if (has_next) {
                const double d_tp = tp1 - tp0;
                tp_inc_m = (d_tp >= 0.0) ? d_tp : tp1;
                const double d_ssr = ssr1 - ssr0;
                ssr_inc_jm2 = (d_ssr >= 0.0) ? d_ssr : ssr1;
            }

            double prcp_mm_day_val = 0.0;
            double rn_wm2_val = 0.0;
            if (has_next) {
                prcp_mm_day_val = tp_inc_m * 1000.0 * (86400.0 / dt_sec);
                rn_wm2_val = ssr_inc_jm2 / dt_sec;
            }
            if (!std::isfinite(prcp_mm_day_val) || prcp_mm_day_val < 0.0) {
                prcp_mm_day_val = 0.0;
            }
            if (!std::isfinite(rn_wm2_val)) {
                rn_wm2_val = 0.0;
            }

            const double temp_c_val = t2m_k - 273.15;
            const double td_c = d2m_k - 273.15;

            const double es = 6.112 * std::exp(17.67 * temp_c_val / (temp_c_val + 243.5));
            const double ea = 6.112 * std::exp(17.67 * td_c / (td_c + 243.5));
            double rh_val = 0.0;
            if (std::isfinite(es) && es > 0.0 && std::isfinite(ea)) {
                rh_val = ea / es;
            }
            if (!std::isfinite(rh_val)) {
                rh_val = 0.0;
            }
            rh_val = std::min(1.0, std::max(0.0, rh_val));

            const double wind_val = std::sqrt(u10 * u10 + v10 * v10);

            prcp_mm_day[i] = prcp_mm_day_val;
            temp_c[i] = temp_c_val;
            rh_1[i] = rh_val;
            wind_ms[i] = std::fabs(wind_val);
            rn_wm2[i] = rn_wm2_val;
        }

        if (next_in_other_file) {
            tp_next.close();
            ssr_next.close();
        }

        loaded_time_idx = t_idx;
    }

    void loadCacheForTimeIndex(int t_idx)
    {
        if (product == "CMFD2") {
            loadCacheForTimeIndexCmfd(t_idx);
            return;
        }
        if (product == "ERA5") {
            loadCacheForTimeIndexEra5(t_idx);
            return;
        }
        fprintf(stderr, "\n  Fatal Error: Unsupported NetCDF forcing PRODUCT: %s\n", product.c_str());
        myexit(ERRFileIO);
    }

    void init()
    {
        parseCfg();
        if (product == "CMFD2") {
            initCmfdMonths();
            discoverGridAndStationsCmfd();
            buildGlobalTimeAxisCmfd();
        } else if (product == "ERA5") {
            initEra5Days();
            discoverGridAndStationsEra5();
            buildGlobalTimeAxisEra5();
        } else {
            fprintf(stderr, "\n  Fatal Error: NetCDF forcing PRODUCT is not supported: %s\n", product.c_str());
            fprintf(stderr, "  Implemented: CMFD2, ERA5\n");
            myexit(ERRFileIO);
        }

        prcp_mm_day.assign(stations.size(), 0.0);
        temp_c.assign(stations.size(), 0.0);
        rh_1.assign(stations.size(), 0.0);
        wind_ms.assign(stations.size(), 0.0);
        rn_wm2.assign(stations.size(), 0.0);

        now_time_idx = 0;
        loaded_time_idx = -1;
        open_month_idx = -1;
        open_day_idx = -1;
    }

    void movePointer(double t_min_in)
    {
        if (time_min.empty()) {
            return;
        }
        // Advance pointer (monotonic calls expected)
        while (now_time_idx + 1 < (int)time_min.size() && t_min_in + 1e-12 >= time_min[(size_t)now_time_idx + 1]) {
            now_time_idx++;
        }
        if (loaded_time_idx != now_time_idx) {
            loadCacheForTimeIndex(now_time_idx);
        }
    }

    double currentTimeMin() const
    {
        if (time_min.empty()) {
            return NA_VALUE;
        }
        return time_min[(size_t)now_time_idx];
    }

    double nextTimeMin() const
    {
        if (time_min.empty()) {
            return NA_VALUE;
        }
        const size_t next = (size_t)now_time_idx + 1;
        if (next >= time_min.size()) {
            return NA_VALUE;
        }
        return time_min[next];
    }

    double minTimeMin() const { return time_min.empty() ? NA_VALUE : time_min.front(); }
    double maxTimeMin() const
    {
        if (time_min.empty()) {
            return NA_VALUE;
        }
        const double t_max = time_min.back();
        // Step-function forcing: the last record covers one more interval beyond its timestamp.
        // Estimate that final interval width from the last positive dt in the time axis.
        double dt_last = 0.0;
        for (size_t i = time_min.size(); i > 1; i--) {
            const double dt = time_min[i - 1] - time_min[i - 2];
            if (dt > 1e-9) {
                dt_last = dt;
                break;
            }
        }
        return t_max + dt_last;
    }
};

NetcdfForcingProvider::NetcdfForcingProvider(const char *forcing_cfg_path,
                                             const std::vector<ForcingStationMeta> &stations,
                                             long forc_start_yyyymmdd,
                                             double sim_start_min,
                                             double sim_end_min)
{
    impl_ = new Impl();
    impl_->cfg_path = forcing_cfg_path ? std::string(forcing_cfg_path) : std::string();
    impl_->stations = stations;
    impl_->forc_start_yyyymmdd = forc_start_yyyymmdd;
    impl_->sim_start_min = sim_start_min;
    impl_->sim_end_min = sim_end_min;

    impl_->init();
}

NetcdfForcingProvider::~NetcdfForcingProvider()
{
    delete impl_;
    impl_ = nullptr;
}

int NetcdfForcingProvider::numStations() const
{
    return (int)impl_->stations.size();
}

void NetcdfForcingProvider::movePointer(double t_min)
{
    impl_->movePointer(t_min);
}

double NetcdfForcingProvider::get(int station_idx, int column) const
{
    if (station_idx < 0 || station_idx >= (int)impl_->stations.size()) {
        return NA_VALUE;
    }
    const size_t i = (size_t)station_idx;
    switch (column) {
        case i_prcp:
            return impl_->prcp_mm_day[i];
        case i_temp:
            return impl_->temp_c[i];
        case i_rh:
            return impl_->rh_1[i];
        case i_wind:
            return impl_->wind_ms[i];
        case i_rn:
            return impl_->rn_wm2[i];
        default:
            return NA_VALUE;
    }
}

double NetcdfForcingProvider::currentTimeMin(int station_idx) const
{
    (void)station_idx;
    return impl_->currentTimeMin();
}

double NetcdfForcingProvider::nextTimeMin(int station_idx) const
{
    (void)station_idx;
    return impl_->nextTimeMin();
}

double NetcdfForcingProvider::lon(int station_idx) const
{
    if (station_idx < 0 || station_idx >= (int)impl_->stations.size()) {
        return NA_VALUE;
    }
    return impl_->stations[(size_t)station_idx].lon_deg;
}

double NetcdfForcingProvider::lat(int station_idx) const
{
    if (station_idx < 0 || station_idx >= (int)impl_->stations.size()) {
        return NA_VALUE;
    }
    return impl_->stations[(size_t)station_idx].lat_deg;
}

double NetcdfForcingProvider::z(int station_idx) const
{
    if (station_idx < 0 || station_idx >= (int)impl_->stations.size()) {
        return NA_VALUE;
    }
    return impl_->stations[(size_t)station_idx].z_m;
}

double NetcdfForcingProvider::minTimeMin() const
{
    return impl_->minTimeMin();
}

double NetcdfForcingProvider::maxTimeMin() const
{
    return impl_->maxTimeMin();
}

#endif /* _NETCDF_ON */
