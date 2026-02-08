#include "NetcdfOutputContext.hpp"

#ifdef _NETCDF_ON

#include "Element.hpp"
#include "Model_Control.hpp"
#include "Node.hpp"

#include <netcdf.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <ctime>
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

static std::string readTextFileOrDie(const std::string &path, const char *what)
{
    std::ifstream f(path);
    if (!f.is_open()) {
        fprintf(stderr, "\n  Fatal Error: Failed to open %s.\n", what ? what : "text file");
        fprintf(stderr, "  File: %s\n", path.c_str());
        myexit(ERRFileIO);
    }
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
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

static std::string nowUtcIso8601()
{
    std::time_t t = std::time(nullptr);
    std::tm tm{};
    gmtime_r(&t, &tm);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return std::string(buf);
}

struct VarMeta {
    const char *units;
    const char *long_name;
};

static bool lookupVarMeta(const std::string &name, VarMeta &out)
{
    static const std::map<std::string, VarMeta> kMeta = {
        // Element storage (m)
        {"eleyic", {"m", "interception storage"}},
        {"eleysnow", {"m", "snow depth"}},
        {"eleysurf", {"m", "surface water depth"}},
        {"eleyunsat", {"m", "unsaturated zone storage depth"}},
        {"eleygw", {"m", "groundwater head"}},

        // Element flux (m/day)
        {"elevprcp", {"m/day", "precipitation to land"}},
        {"elevnetprcp", {"m/day", "net precipitation"}},
        {"elevetp", {"m/day", "potential evapotranspiration"}},
        {"eleveta", {"m/day", "actual evapotranspiration"}},
        {"elevrech", {"m/day", "recharge"}},
        {"elevinfil", {"m/day", "infiltration"}},
        {"elevexfil", {"m/day", "exfiltration"}},
        {"elevetic", {"m/day", "evapotranspiration: interception"}},
        {"elevettr", {"m/day", "evapotranspiration: transpiration"}},
        {"elevetev", {"m/day", "evapotranspiration: evaporation"}},

        // Element/river exchange (m3/day)
        {"eleqrsurf", {"m3/day", "element to river surface flow"}},
        {"eleqrsub", {"m3/day", "element to river subsurface flow"}},

        // Element inter-cell flux (m3/day)
        {"eleqsub", {"m3/day", "subsurface flow: total"}},
        {"eleqsub1", {"m3/day", "subsurface flow: direction 1"}},
        {"eleqsub2", {"m3/day", "subsurface flow: direction 2"}},
        {"eleqsub3", {"m3/day", "subsurface flow: direction 3"}},
        {"eleqsurf", {"m3/day", "surface flow: total"}},
        {"eleqsurf1", {"m3/day", "surface flow: direction 1"}},
        {"eleqsurf2", {"m3/day", "surface flow: direction 2"}},
        {"eleqsurf3", {"m3/day", "surface flow: direction 3"}},

        // River (m3/day or m)
        {"rivqdown", {"m3/day", "river downstream discharge"}},
        {"rivqup", {"m3/day", "river upstream discharge"}},
        {"rivqsurf", {"m3/day", "river surface discharge"}},
        {"rivqsub", {"m3/day", "river subsurface discharge"}},
        {"rivystage", {"m", "river stage"}},

        // Lake (m / m2 / m3/day / m/day)
        {"lakystage", {"m", "lake stage"}},
        {"lakatop", {"m2", "lake top area"}},
        {"lakvevap", {"m/day", "lake evaporation"}},
        {"lakvprcp", {"m/day", "lake precipitation"}},
        {"lakqrivin", {"m3/day", "lake river inflow"}},
        {"lakqrivout", {"m3/day", "lake river outflow"}},
        {"lakqsurf", {"m3/day", "lake surface discharge"}},
        {"lakqsub", {"m3/day", "lake subsurface discharge"}},

        // TSR debug outputs
        {"rn_h", {"W m-2", "shortwave radiation on horizontal surface"}},
        {"rn_t", {"W m-2", "terrain-corrected shortwave radiation"}},
        {"rn_factor", {"1", "terrain radiation correction factor"}},
    };

    const auto it = kMeta.find(name);
    if (it == kMeta.end()) {
        return false;
    }
    out = it->second;
    return true;
}

static void applyVarMetadata(int ncid, int varid, const std::string &name, const char *path)
{
    VarMeta meta{};
    const bool found = lookupVarMeta(name, meta);

    const std::string long_name = found ? std::string(meta.long_name) : ("SHUD output variable: " + name);
    ncCheck(nc_put_att_text(ncid, varid, "long_name", long_name.size(), long_name.c_str()),
            "nc_put_att_text(var.long_name)",
            path);

    if (found && meta.units != nullptr && meta.units[0] != '\0') {
        ncCheck(nc_put_att_text(ncid, varid, "units", strlen(meta.units), meta.units),
                "nc_put_att_text(var.units)",
                path);
    }
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
    NetcdfElementFile(std::string out_dir_abs,
                      float fill_value,
                      std::string crs_wkt,
                      const _Node *nodes,
                      int num_nodes,
                      const _Element *elements,
                      int num_elements)
        : out_dir_abs_(std::move(out_dir_abs)),
          fill_value_(fill_value),
          crs_wkt_(std::move(crs_wkt)),
          nodes_(nodes),
          num_nodes_(num_nodes),
          elements_(elements),
          num_elements_(num_elements)
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
        if (num_elements_ > 0 && n_all_ != num_elements_) {
            fprintf(stderr, "\n  Fatal Error: mesh_face dimension mismatch for NetCDF element output.\n");
            fprintf(stderr, "  Expected NumEle=%d, got n_all=%d\n", num_elements_, n_all_);
            myexit(ERRFileIO);
        }

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
        ncCheck(nc_def_dim(ncid_, "mesh_node", (size_t)num_nodes_, &dim_node_),
                "nc_def_dim(mesh_node)",
                file_path_.c_str());
        ncCheck(nc_def_dim(ncid_, "max_face_nodes", 3, &dim_max_face_nodes_),
                "nc_def_dim(max_face_nodes)",
                file_path_.c_str());

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

        // UGRID mesh topology and coordinates (Phase B)
        const int node_dims[1] = {dim_node_};
        ncCheck(nc_def_var(ncid_, "mesh_node_x", NC_DOUBLE, 1, node_dims, &var_node_x_),
                "nc_def_var(mesh_node_x)",
                file_path_.c_str());
        ncCheck(nc_def_var(ncid_, "mesh_node_y", NC_DOUBLE, 1, node_dims, &var_node_y_),
                "nc_def_var(mesh_node_y)",
                file_path_.c_str());
        const char *xname = "projection_x_coordinate";
        const char *yname = "projection_y_coordinate";
        ncCheck(nc_put_att_text(ncid_, var_node_x_, "standard_name", strlen(xname), xname),
                "nc_put_att_text(mesh_node_x.standard_name)",
                file_path_.c_str());
        ncCheck(nc_put_att_text(ncid_, var_node_y_, "standard_name", strlen(yname), yname),
                "nc_put_att_text(mesh_node_y.standard_name)",
                file_path_.c_str());
        const char *xy_units = "m";
        ncCheck(nc_put_att_text(ncid_, var_node_x_, "units", strlen(xy_units), xy_units),
                "nc_put_att_text(mesh_node_x.units)",
                file_path_.c_str());
        ncCheck(nc_put_att_text(ncid_, var_node_y_, "units", strlen(xy_units), xy_units),
                "nc_put_att_text(mesh_node_y.units)",
                file_path_.c_str());

        const int face_node_dims[2] = {dim_face_, dim_max_face_nodes_};
        ncCheck(nc_def_var(ncid_, "mesh_face_nodes", NC_INT, 2, face_node_dims, &var_face_nodes_),
                "nc_def_var(mesh_face_nodes)",
                file_path_.c_str());
        const int start_index = 1;
        ncCheck(nc_put_att_int(ncid_, var_face_nodes_, "start_index", NC_INT, 1, &start_index),
                "nc_put_att_int(mesh_face_nodes.start_index)",
                file_path_.c_str());

        const int face_dims[1] = {dim_face_};
        ncCheck(nc_def_var(ncid_, "mesh_face_x", NC_DOUBLE, 1, face_dims, &var_face_x_),
                "nc_def_var(mesh_face_x)",
                file_path_.c_str());
        ncCheck(nc_def_var(ncid_, "mesh_face_y", NC_DOUBLE, 1, face_dims, &var_face_y_),
                "nc_def_var(mesh_face_y)",
                file_path_.c_str());
        ncCheck(nc_put_att_text(ncid_, var_face_x_, "standard_name", strlen(xname), xname),
                "nc_put_att_text(mesh_face_x.standard_name)",
                file_path_.c_str());
        ncCheck(nc_put_att_text(ncid_, var_face_y_, "standard_name", strlen(yname), yname),
                "nc_put_att_text(mesh_face_y.standard_name)",
                file_path_.c_str());
        ncCheck(nc_put_att_text(ncid_, var_face_x_, "units", strlen(xy_units), xy_units),
                "nc_put_att_text(mesh_face_x.units)",
                file_path_.c_str());
        ncCheck(nc_put_att_text(ncid_, var_face_y_, "units", strlen(xy_units), xy_units),
                "nc_put_att_text(mesh_face_y.units)",
                file_path_.c_str());

        if (!crs_wkt_.empty()) {
            ncCheck(nc_def_var(ncid_, "crs", NC_INT, 0, nullptr, &var_crs_),
                    "nc_def_var(crs)",
                    file_path_.c_str());
            const char *crs_long_name = "coordinate reference system";
            ncCheck(nc_put_att_text(ncid_, var_crs_, "long_name", strlen(crs_long_name), crs_long_name),
                    "nc_put_att_text(crs.long_name)",
                    file_path_.c_str());
            ncCheck(nc_put_att_text(ncid_, var_crs_, "spatial_ref", crs_wkt_.size(), crs_wkt_.c_str()),
                    "nc_put_att_text(crs.spatial_ref)",
                    file_path_.c_str());
            ncCheck(nc_put_att_text(ncid_, var_crs_, "crs_wkt", crs_wkt_.size(), crs_wkt_.c_str()),
                    "nc_put_att_text(crs.crs_wkt)",
                    file_path_.c_str());
        }

        int var_mesh_ = -1;
        ncCheck(nc_def_var(ncid_, "mesh", NC_INT, 0, nullptr, &var_mesh_), "nc_def_var(mesh)", file_path_.c_str());
        const char *cf_role = "mesh_topology";
        ncCheck(nc_put_att_text(ncid_, var_mesh_, "cf_role", strlen(cf_role), cf_role),
                "nc_put_att_text(mesh.cf_role)",
                file_path_.c_str());
        const int topo_dim = 2;
        ncCheck(nc_put_att_int(ncid_, var_mesh_, "topology_dimension", NC_INT, 1, &topo_dim),
                "nc_put_att_int(mesh.topology_dimension)",
                file_path_.c_str());
        const char *node_coords = "mesh_node_x mesh_node_y";
        ncCheck(nc_put_att_text(ncid_, var_mesh_, "node_coordinates", strlen(node_coords), node_coords),
                "nc_put_att_text(mesh.node_coordinates)",
                file_path_.c_str());
        const char *face_nodes = "mesh_face_nodes";
        ncCheck(nc_put_att_text(ncid_, var_mesh_, "face_node_connectivity", strlen(face_nodes), face_nodes),
                "nc_put_att_text(mesh.face_node_connectivity)",
                file_path_.c_str());
        const char *face_coords = "mesh_face_x mesh_face_y";
        ncCheck(nc_put_att_text(ncid_, var_mesh_, "face_coordinates", strlen(face_coords), face_coords),
                "nc_put_att_text(mesh.face_coordinates)",
                file_path_.c_str());

        // Element index and mask (Phase B: full-dimension outputs)
        ncCheck(nc_def_var(ncid_, "element_id", NC_INT, 1, face_dims, &var_element_id_),
                "nc_def_var(element_id)",
                file_path_.c_str());
        const int element_start_index = 1;
        ncCheck(nc_put_att_int(ncid_, var_element_id_, "start_index", NC_INT, 1, &element_start_index),
                "nc_put_att_int(element_id.start_index)",
                file_path_.c_str());

        ncCheck(nc_def_var(ncid_, "element_output_mask", NC_INT, 1, face_dims, &var_output_mask_),
                "nc_def_var(element_output_mask)",
                file_path_.c_str());

        const char *conventions = "CF-1.10";
        ncCheck(nc_put_att_text(ncid_, NC_GLOBAL, "Conventions", strlen(conventions), conventions),
                "nc_put_att_text(global.Conventions)",
                file_path_.c_str());

        const std::string title = "SHUD output: " + prefix_leaf;
        ncCheck(nc_put_att_text(ncid_, NC_GLOBAL, "title", title.size(), title.c_str()),
                "nc_put_att_text(global.title)",
                file_path_.c_str());
        const char *source = "SHUD";
        ncCheck(nc_put_att_text(ncid_, NC_GLOBAL, "source", strlen(source), source),
                "nc_put_att_text(global.source)",
                file_path_.c_str());
        const std::string history = nowUtcIso8601() + ": created by SHUD";
        ncCheck(nc_put_att_text(ncid_, NC_GLOBAL, "history", history.size(), history.c_str()),
                "nc_put_att_text(global.history)",
                file_path_.c_str());

        ncCheck(nc_enddef(ncid_), "nc_enddef", file_path_.c_str());

        // Write mesh coordinate/connectivity data.
        if (num_nodes_ > 0 && nodes_ == nullptr) {
            fprintf(stderr, "\n  Fatal Error: mesh nodes are not available for NetCDF output.\n");
            myexit(ERRFileIO);
        }
        if (num_elements_ > 0 && elements_ == nullptr) {
            fprintf(stderr, "\n  Fatal Error: mesh elements are not available for NetCDF output.\n");
            myexit(ERRFileIO);
        }

        std::vector<double> node_x((size_t)num_nodes_, NA_VALUE);
        std::vector<double> node_y((size_t)num_nodes_, NA_VALUE);
        for (int i = 0; i < num_nodes_; i++) {
            const int idx1 = nodes_[i].index;
            if (idx1 < 1 || idx1 > num_nodes_) {
                fprintf(stderr, "\n  Fatal Error: invalid node index %d (NumNode=%d).\n", idx1, num_nodes_);
                myexit(ERRFileIO);
            }
            node_x[(size_t)(idx1 - 1)] = nodes_[i].x;
            node_y[(size_t)(idx1 - 1)] = nodes_[i].y;
        }
        ncCheck(nc_put_var_double(ncid_, var_node_x_, node_x.data()), "nc_put_var_double(mesh_node_x)", file_path_.c_str());
        ncCheck(nc_put_var_double(ncid_, var_node_y_, node_y.data()), "nc_put_var_double(mesh_node_y)", file_path_.c_str());

        std::vector<int> face_nodes_buf((size_t)num_elements_ * 3u, 0);
        std::vector<double> face_x((size_t)num_elements_, NA_VALUE);
        std::vector<double> face_y((size_t)num_elements_, NA_VALUE);
        for (int e = 0; e < num_elements_; e++) {
            for (int k = 0; k < 3; k++) {
                const int n1 = elements_[e].node[k];
                if (n1 < 1 || n1 > num_nodes_) {
                    fprintf(stderr, "\n  Fatal Error: invalid element node index %d (NumNode=%d).\n", n1, num_nodes_);
                    myexit(ERRFileIO);
                }
                face_nodes_buf[(size_t)e * 3u + (size_t)k] = n1;
            }
            face_x[(size_t)e] = elements_[e].x;
            face_y[(size_t)e] = elements_[e].y;
        }
        ncCheck(nc_put_var_int(ncid_, var_face_nodes_, face_nodes_buf.data()),
                "nc_put_var_int(mesh_face_nodes)",
                file_path_.c_str());
        ncCheck(nc_put_var_double(ncid_, var_face_x_, face_x.data()), "nc_put_var_double(mesh_face_x)", file_path_.c_str());
        ncCheck(nc_put_var_double(ncid_, var_face_y_, face_y.data()), "nc_put_var_double(mesh_face_y)", file_path_.c_str());

        std::vector<int> element_ids((size_t)n_all_, 0);
        for (int i = 0; i < n_all_; i++) {
            element_ids[(size_t)i] = i + 1;
        }
        ncCheck(nc_put_var_int(ncid_, var_element_id_, element_ids.data()), "nc_put_var_int(element_id)", file_path_.c_str());
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
        const char *mesh = "mesh";
        ncCheck(nc_put_att_text(ncid_, varid, "mesh", strlen(mesh), mesh),
                "nc_put_att_text(data.mesh)",
                file_path_.c_str());
        const char *loc = "face";
        ncCheck(nc_put_att_text(ncid_, varid, "location", strlen(loc), loc),
                "nc_put_att_text(data.location)",
                file_path_.c_str());
        const char *coords = "mesh_face_x mesh_face_y";
        ncCheck(nc_put_att_text(ncid_, varid, "coordinates", strlen(coords), coords),
                "nc_put_att_text(data.coordinates)",
                file_path_.c_str());
        if (var_crs_ != -1) {
            const char *grid_mapping = "crs";
            ncCheck(nc_put_att_text(ncid_, varid, "grid_mapping", strlen(grid_mapping), grid_mapping),
                    "nc_put_att_text(data.grid_mapping)",
                    file_path_.c_str());
        }
        applyVarMetadata(ncid_, varid, name, file_path_.c_str());
        ncCheck(nc_enddef(ncid_), "nc_enddef", file_path_.c_str());

        varids_[name] = varid;
        return varid;
    }

    void writeOutputMask(const std::vector<int> &icol_1based)
    {
        if (mask_written_) {
            return;
        }
        if (ncid_ == -1 || var_output_mask_ == -1) {
            fprintf(stderr, "\n  Fatal Error: NetCDF element file is not initialized before writing mask.\n");
            myexit(ERRFileIO);
        }

        std::vector<int> mask((size_t)n_all_, 0);
        for (const int idx1 : icol_1based) {
            const int idx0 = idx1 - 1;
            if (idx0 < 0 || idx0 >= n_all_) {
                continue;
            }
            mask[(size_t)idx0] = 1;
        }

        ncCheck(nc_put_var_int(ncid_, var_output_mask_, mask.data()),
                "nc_put_var_int(element_output_mask)",
                file_path_.c_str());
        mask_written_ = true;
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
    std::string crs_wkt_;

    std::string file_path_;
    int ncid_ = -1;
    int dim_time_ = -1;
    int dim_face_ = -1;
    int dim_node_ = -1;
    int dim_max_face_nodes_ = -1;
    int var_time_ = -1;
    int var_crs_ = -1;
    int var_node_x_ = -1;
    int var_node_y_ = -1;
    int var_face_nodes_ = -1;
    int var_face_x_ = -1;
    int var_face_y_ = -1;
    int var_element_id_ = -1;
    int var_output_mask_ = -1;
    int n_all_ = 0;
    long forc_start_yyyymmdd_ = 0;
    bool mask_written_ = false;

    const _Node *nodes_ = nullptr;
    int num_nodes_ = 0;
    const _Element *elements_ = nullptr;
    int num_elements_ = 0;

    std::map<std::string, int> varids_;
    std::map<long long, size_t> time_to_idx_;
    std::vector<double> times_;
};

class NetcdfSimpleFile {
public:
    NetcdfSimpleFile(std::string out_dir_abs,
                     float fill_value,
                     std::string obj_dim_name,
                     std::string obj_id_var_name,
                     std::string file_suffix)
        : out_dir_abs_(std::move(out_dir_abs)),
          fill_value_(fill_value),
          obj_dim_name_(std::move(obj_dim_name)),
          obj_id_var_name_(std::move(obj_id_var_name)),
          file_suffix_(std::move(file_suffix))
    {
    }

    ~NetcdfSimpleFile() { close(); }

    void openIfNeeded(const char *legacy_basename, long forc_start_yyyymmdd, int n_all)
    {
        if (ncid_ != -1) {
            if (n_all_ != n_all || forc_start_yyyymmdd_ != forc_start_yyyymmdd) {
                fprintf(stderr, "\n  Fatal Error: Inconsistent NetCDF output context.\n");
                fprintf(stderr, "  Existing: n_all=%d ForcStartTime=%ld\n", n_all_, forc_start_yyyymmdd_);
                fprintf(stderr, "  New:      n_all=%d ForcStartTime=%ld\n", n_all, forc_start_yyyymmdd);
                myexit(ERRFileIO);
            }
            return;
        }

        if (n_all <= 0) {
            fprintf(stderr, "\n  Fatal Error: Invalid output dimension for NetCDF output (%s).\n", obj_dim_name_.c_str());
            fprintf(stderr, "  n_all: %d\n", n_all);
            myexit(ERRFileIO);
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

        file_path_ = joinPath(out_dir, prefix_leaf + file_suffix_);

        const int cmode = NC_NETCDF4 | NC_CLASSIC_MODEL | NC_CLOBBER;
        ncCheck(nc_create(file_path_.c_str(), cmode, &ncid_), "nc_create", file_path_.c_str());

        ncCheck(nc_def_dim(ncid_, "time", NC_UNLIMITED, &dim_time_), "nc_def_dim(time)", file_path_.c_str());
        ncCheck(nc_def_dim(ncid_, obj_dim_name_.c_str(), (size_t)n_all_, &dim_obj_),
                "nc_def_dim(object)",
                file_path_.c_str());

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

        const int obj_dims[1] = {dim_obj_};
        ncCheck(nc_def_var(ncid_, obj_id_var_name_.c_str(), NC_INT, 1, obj_dims, &var_obj_id_),
                "nc_def_var(object_id)",
                file_path_.c_str());
        const int start_index = 1;
        ncCheck(nc_put_att_int(ncid_, var_obj_id_, "start_index", NC_INT, 1, &start_index),
                "nc_put_att_int(object_id.start_index)",
                file_path_.c_str());

        const std::string mask_name = obj_dim_name_ + "_output_mask";
        ncCheck(nc_def_var(ncid_, mask_name.c_str(), NC_INT, 1, obj_dims, &var_output_mask_),
                "nc_def_var(output_mask)",
                file_path_.c_str());

        const char *conventions = "CF-1.10 UGRID-1.0";
        ncCheck(nc_put_att_text(ncid_, NC_GLOBAL, "Conventions", strlen(conventions), conventions),
                "nc_put_att_text(global.Conventions)",
                file_path_.c_str());

        const std::string title = "SHUD output: " + prefix_leaf;
        ncCheck(nc_put_att_text(ncid_, NC_GLOBAL, "title", title.size(), title.c_str()),
                "nc_put_att_text(global.title)",
                file_path_.c_str());
        const char *source = "SHUD";
        ncCheck(nc_put_att_text(ncid_, NC_GLOBAL, "source", strlen(source), source),
                "nc_put_att_text(global.source)",
                file_path_.c_str());
        const std::string history = nowUtcIso8601() + ": created by SHUD";
        ncCheck(nc_put_att_text(ncid_, NC_GLOBAL, "history", history.size(), history.c_str()),
                "nc_put_att_text(global.history)",
                file_path_.c_str());

        ncCheck(nc_enddef(ncid_), "nc_enddef", file_path_.c_str());

        std::vector<int> ids((size_t)n_all_, 0);
        for (int i = 0; i < n_all_; i++) {
            ids[(size_t)i] = i + 1;
        }
        ncCheck(nc_put_var_int(ncid_, var_obj_id_, ids.data()), "nc_put_var_int(object_id)", file_path_.c_str());
    }

    int ensureVar(const std::string &name)
    {
        const auto it = varids_.find(name);
        if (it != varids_.end()) {
            return it->second;
        }
        if (ncid_ == -1) {
            fprintf(stderr, "\n  Fatal Error: NetCDF output file is not initialized before defining variables.\n");
            myexit(ERRFileIO);
        }

        ncCheck(nc_redef(ncid_), "nc_redef", file_path_.c_str());
        const int dims[2] = {dim_time_, dim_obj_};
        int varid = -1;
        ncCheck(nc_def_var(ncid_, name.c_str(), NC_FLOAT, 2, dims, &varid),
                "nc_def_var(data)",
                file_path_.c_str());
        ncCheck(nc_put_att_float(ncid_, varid, "_FillValue", NC_FLOAT, 1, &fill_value_),
                "nc_put_att_float(_FillValue)",
                file_path_.c_str());
        applyVarMetadata(ncid_, varid, name, file_path_.c_str());
        ncCheck(nc_enddef(ncid_), "nc_enddef", file_path_.c_str());

        varids_[name] = varid;
        return varid;
    }

    void writeOutputMask(const std::vector<int> &icol_1based)
    {
        if (mask_written_) {
            return;
        }
        if (ncid_ == -1 || var_output_mask_ == -1) {
            fprintf(stderr, "\n  Fatal Error: NetCDF output file is not initialized before writing mask.\n");
            myexit(ERRFileIO);
        }

        std::vector<int> mask((size_t)n_all_, 0);
        for (const int idx1 : icol_1based) {
            const int idx0 = idx1 - 1;
            if (idx0 < 0 || idx0 >= n_all_) {
                continue;
            }
            mask[(size_t)idx0] = 1;
        }

        const std::string mask_name = obj_dim_name_ + "_output_mask";
        ncCheck(nc_put_var_int(ncid_, var_output_mask_, mask.data()), ("nc_put_var_int(" + mask_name + ")").c_str(), file_path_.c_str());
        mask_written_ = true;
    }

    size_t ensureTimeIndex(long long t_min)
    {
        const auto it = time_to_idx_.find(t_min);
        if (it != time_to_idx_.end()) {
            return it->second;
        }
        if (ncid_ == -1) {
            fprintf(stderr, "\n  Fatal Error: NetCDF output file is not initialized before writing time.\n");
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
            fprintf(stderr, "\n  Fatal Error: NetCDF output file is not initialized before writing data.\n");
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
    std::string obj_dim_name_;
    std::string obj_id_var_name_;
    std::string file_suffix_;

    std::string file_path_;
    int ncid_ = -1;
    int dim_time_ = -1;
    int dim_obj_ = -1;
    int var_time_ = -1;
    int var_obj_id_ = -1;
    int var_output_mask_ = -1;
    int n_all_ = 0;
    long forc_start_yyyymmdd_ = 0;
    bool mask_written_ = false;

    std::map<std::string, int> varids_;
    std::map<long long, size_t> time_to_idx_;
    std::vector<double> times_;
};

class NetcdfElementVarSink final : public IPrintSink {
public:
    explicit NetcdfElementVarSink(std::shared_ptr<NetcdfElementFile> file) : file_(std::move(file)) {}

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
        file_->writeOutputMask(icol_);
        varid_ = file_->ensureVar(var_name_);
        initialized_ = true;
    }

    void onWrite(double t_quantized_min, int num_selected, const double *buffer) override
    {
        if (!initialized_ || !file_ || buffer == nullptr) {
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
    std::shared_ptr<NetcdfElementFile> file_;
    bool initialized_ = false;

    std::string var_name_;
    int varid_ = -1;
    int n_all_ = 0;
    float fill_value_ = 9.96921e36f;

    std::vector<int> icol_;
};

class NetcdfSimpleVarSink final : public IPrintSink {
public:
    explicit NetcdfSimpleVarSink(std::shared_ptr<NetcdfSimpleFile> file) : file_(std::move(file)) {}

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
        file_->writeOutputMask(icol_);
        varid_ = file_->ensureVar(var_name_);
        initialized_ = true;
    }

    void onWrite(double t_quantized_min, int num_selected, const double *buffer) override
    {
        if (!initialized_ || !file_ || buffer == nullptr) {
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
    std::shared_ptr<NetcdfSimpleFile> file_;
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
    std::string crs_wkt_abs;
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

    const std::string crs_wkt_cfg = get("CRS_WKT");
    if (!crs_wkt_cfg.empty()) {
        cfg.crs_wkt_abs = isAbsPath(crs_wkt_cfg) ? crs_wkt_cfg : joinPath(run_dir, crs_wkt_cfg);
    }

    return cfg;
}

} // namespace

struct NetcdfOutputContext::Impl {
    std::string ncoutput_cfg_path;
    ParsedNcOutputCfg cfg;
    std::string crs_wkt;
    std::shared_ptr<NetcdfElementFile> ele_file;
    std::shared_ptr<NetcdfSimpleFile> riv_file;
    std::shared_ptr<NetcdfSimpleFile> lake_file;

    Impl(const char *path, const _Node *nodes, int num_nodes, const _Element *elements, int num_elements)
        : ncoutput_cfg_path(path ? path : ""),
          cfg(parseNcOutputCfg(ncoutput_cfg_path))
    {
        if (!cfg.schema_abs.empty()) {
            fprintf(stderr,
                    "WARNING: NCOUTPUT_CFG SCHEMA is currently ignored (not implemented): %s\n",
                    cfg.schema_abs.c_str());
        }
        if (!cfg.crs_wkt_abs.empty()) {
            crs_wkt = readTextFileOrDie(cfg.crs_wkt_abs, "CRS_WKT");
        } else {
            fprintf(stderr,
                    "WARNING: NetCDF element outputs do not embed CRS information. "
                    "Set NCOUTPUT_CFG CRS_WKT <path> to include a WKT CRS.\n");
        }

        ele_file = std::make_shared<NetcdfElementFile>(cfg.out_dir_abs,
                                                       9.96921e36f,
                                                       crs_wkt,
                                                       nodes,
                                                       num_nodes,
                                                       elements,
                                                       num_elements);
        riv_file = std::make_shared<NetcdfSimpleFile>(cfg.out_dir_abs, 9.96921e36f, "river", "river_id", ".riv.nc");
        lake_file = std::make_shared<NetcdfSimpleFile>(cfg.out_dir_abs, 9.96921e36f, "lake", "lake_id", ".lake.nc");
    }

    std::shared_ptr<IPrintSink> createElementSink()
    {
        return std::make_shared<NetcdfElementVarSink>(ele_file);
    }

    std::shared_ptr<IPrintSink> createRiverSink()
    {
        return std::make_shared<NetcdfSimpleVarSink>(riv_file);
    }

    std::shared_ptr<IPrintSink> createLakeSink()
    {
        return std::make_shared<NetcdfSimpleVarSink>(lake_file);
    }
};

NetcdfOutputContext::NetcdfOutputContext(const char *ncoutput_cfg_path,
                                         const _Node *nodes,
                                         int num_nodes,
                                         const _Element *elements,
                                         int num_elements)
    : impl_(std::make_unique<Impl>(ncoutput_cfg_path, nodes, num_nodes, elements, num_elements))
{
}

NetcdfOutputContext::~NetcdfOutputContext() = default;

std::shared_ptr<IPrintSink> NetcdfOutputContext::createElementSink()
{
    return impl_->createElementSink();
}

std::shared_ptr<IPrintSink> NetcdfOutputContext::createRiverSink()
{
    return impl_->createRiverSink();
}

std::shared_ptr<IPrintSink> NetcdfOutputContext::createLakeSink()
{
    return impl_->createLakeSink();
}

#endif /* _NETCDF_ON */
