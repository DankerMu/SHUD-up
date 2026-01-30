#include "WaterBalanceDiag.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <strings.h>

#include "Model_Data.hpp"

static int choose_interval_min(const Control_Data &cs)
{
    if (cs.dt_qe_et > 0) {
        return cs.dt_qe_et;
    }
    if (cs.dt_ye_surf > 0) {
        return cs.dt_ye_surf;
    }
    return 1440;
}

static bool env_truthy(const char *s)
{
    if (s == nullptr || s[0] == '\0') {
        return false;
    }
    if (strcmp(s, "0") == 0) {
        return false;
    }
    if (strcasecmp(s, "false") == 0) {
        return false;
    }
    if (strcasecmp(s, "off") == 0) {
        return false;
    }
    return true;
}

WaterBalanceDiag::WaterBalanceDiag(Model_Data *md_)
    : md(md_)
{
}

WaterBalanceDiag::~WaterBalanceDiag()
{
    closeFiles();
}

void WaterBalanceDiag::enable()
{
    enabled = 1;
    openFiles();
}

int WaterBalanceDiag::isEnabled() const
{
    return enabled;
}

std::string WaterBalanceDiag::outputPrefix() const
{
    if (md == nullptr || md->pf_out == nullptr) {
        return "";
    }
    return std::string(md->pf_out->outpath) + "/" + md->pf_out->projectname + md->pf_out->suffix;
}

std::vector<double> WaterBalanceDiag::makeIcol(int n) const
{
    std::vector<double> v;
    v.reserve((size_t)n);
    for (int i = 0; i < n; i++) {
        v.push_back((double)(i + 1));
    }
    return v;
}

void WaterBalanceDiag::writeDatHeader(FILE *fp,
                                     const char *label,
                                     int num_var,
                                     const std::vector<double> &col_ids) const
{
    if (fp == nullptr || md == nullptr) {
        return;
    }
    if (num_var < 0) {
        return;
    }
    if ((int)col_ids.size() != num_var) {
        return;
    }

    char header[1024];
    memset(header, 0, sizeof(header));
    snprintf(header,
             sizeof(header),
             "# SHUD output\n"
             "# Water-balance diagnostic: %s\n"
             "# Radiation input mode: %s\n"
             "# Terrain radiation (TSR): %s\n"
             "# Solar lon/lat mode: %s\n"
             "# Solar lon/lat (deg): lon=%.6f, lat=%.6f\n",
             (label != nullptr ? label : "unknown"),
             md->CS.radiation_input_mode == SWNET ? "SWNET" : "SWDOWN",
             md->CS.terrain_radiation ? "ON" : "OFF",
             SolarLonLatModeName(md->CS.solar_lonlat_mode),
             md->CS.solar_lon_deg,
             md->CS.solar_lat_deg);

    fwrite(header, sizeof(char), 1024, fp);
    const double start_time = (double)md->ForcStartTime;
    const double num_var_d = (double)num_var;
    fwrite(&start_time, sizeof(double), 1, fp);
    fwrite(&num_var_d, sizeof(double), 1, fp);
    if (num_var > 0) {
        fwrite(col_ids.data(), sizeof(double), (size_t)num_var, fp);
    }
}

void WaterBalanceDiag::writeDatRecord(FILE *fp, double t_min, int num_var, const std::vector<double> &values) const
{
    if (fp == nullptr || md == nullptr) {
        return;
    }
    if (num_var < 0) {
        return;
    }
    if ((int)values.size() != num_var) {
        return;
    }
    fwrite(&t_min, sizeof(double), 1, fp);
    if (num_var > 0) {
        fwrite(values.data(), sizeof(double), (size_t)num_var, fp);
    }
    if (global_fflush_mode) {
        fflush(fp);
    }
}

void WaterBalanceDiag::initBasinOutlets()
{
    outlet_riv_idx.clear();
    if (md == nullptr) {
        return;
    }
    for (int i = 0; i < md->NumRiv; i++) {
        if (md->Riv[i].toLake >= 0) {
            continue;
        }
        if (md->Riv[i].down < 0) {
            outlet_riv_idx.push_back(i);
        }
    }
}

double WaterBalanceDiag::basinElementStorageFull_m3() const
{
    if (md == nullptr) {
        return 0.0;
    }
    double s = 0.0;
    for (int i = 0; i < md->NumEle; i++) {
        const double sy = md->Ele[i].Sy;
        const double depth =
            md->yEleSurf[i] + sy * md->yEleUnsat[i] + sy * md->yEleGW[i] + md->yEleSnow[i] + md->yEleIS[i];
        s += depth * md->Ele[i].area;
    }
    return s;
}

double WaterBalanceDiag::basinRiverStorage_m3()
{
    if (md == nullptr) {
        return 0.0;
    }
    double s = 0.0;
    for (int i = 0; i < md->NumRiv; i++) {
        const double stage = md->yRivStg[i];
        md->Riv[i].updateRiver(stage);
        s += md->Riv[i].u_CSarea * md->Riv[i].Length;
    }
    return s;
}

double WaterBalanceDiag::basinOutletDischarge_m3min()
{
    if (md == nullptr) {
        return 0.0;
    }
    if (outlet_riv_idx.empty()) {
        return 0.0;
    }

    double q = 0.0;
    for (int idx : outlet_riv_idx) {
        if (idx < 0 || idx >= md->NumRiv) {
            continue;
        }

        const double stage = md->yRivStg[idx];
        md->Riv[idx].updateRiver(stage);

        const double perem = md->Riv[idx].u_CSperem;
        const double csarea = md->Riv[idx].u_CSarea;
        const double r = (perem <= 0.0) ? 0.0 : (csarea / perem);

        switch (md->Riv[idx].down) {
            case -4:
                q += csarea * sqrt(GRAV * stage) * 60.0;
                break;
            default: {
                const double s = md->Riv[idx].BedSlope + stage * 2.0 / md->Riv[idx].Length;
                q += ManningEquation(csarea, md->Riv[idx].avgRough, r, s);
            } break;
        }
    }
    return q;
}

double WaterBalanceDiag::basinBoundaryEdgeOutflow_m3min() const
{
    if (md == nullptr) {
        return 0.0;
    }
    if (md->CS.CloseBoundary) {
        return 0.0;
    }
    if (md->QeleSurf == nullptr || md->QeleSub == nullptr || md->Ele == nullptr) {
        return 0.0;
    }

    double q = 0.0;
    for (int i = 0; i < md->NumEle; i++) {
        for (int j = 0; j < 3; j++) {
            const int inabr = md->Ele[i].nabr[j] - 1;
            const int ilake = md->Ele[i].lakenabr[j] - 1;
            if (inabr >= 0 || ilake >= 0) {
                continue;
            }
            q += md->QeleSurf[i][j] + md->QeleSub[i][j];
        }
    }
    return q;
}

void WaterBalanceDiag::openFiles()
{
    if (!enabled || md == nullptr) {
        return;
    }

    /* Allow disabling via env var without code changes */
    if (!env_truthy(getenv("SHUD_WB_DIAG"))) {
        enabled = 0;
        return;
    }

    const int n = md->NumEle;
    interval_min = choose_interval_min(md->CS);

    ic_raw.assign((size_t)n, 0.0);
    flux3_acc.assign((size_t)n, 0.0);
    fluxfull_acc.assign((size_t)n, 0.0);
    budget3_acc.assign((size_t)n, 0.0);
    budgetfull_acc.assign((size_t)n, 0.0);
    s3_prev.assign((size_t)n, 0.0);
    sfull_prev.assign((size_t)n, 0.0);
    resid3.assign((size_t)n, 0.0);
    residfull.assign((size_t)n, 0.0);
    resid3_budget.assign((size_t)n, 0.0);
    residfull_budget.assign((size_t)n, 0.0);

    icol = makeIcol(n);
    basin_icol = makeIcol(9);
    basin_values.assign(9, 0.0);

    /* Initial storage snapshot (after IC is loaded, before any forcing updates) */
    for (int i = 0; i < n; i++) {
        const double sy = md->Ele[i].Sy;
        const double s3 = md->yEleSurf[i] + sy * md->yEleUnsat[i] + sy * md->yEleGW[i];
        s3_prev[(size_t)i] = s3;
        sfull_prev[(size_t)i] = s3 + md->yEleSnow[i] + md->yEleIS[i];
    }

    const std::string prefix = outputPrefix();
    if (prefix.empty()) {
        enabled = 0;
        return;
    }

    const std::string f3 = prefix + ".elewb3_resid.dat";
    const std::string ffull = prefix + ".elewbfull_resid.dat";
    const std::string f3b = prefix + ".elewb3_budget_resid.dat";
    const std::string ffullb = prefix + ".elewbfull_budget_resid.dat";
    const std::string fbasin = prefix + ".basinwbfull.dat";
    fp3 = fopen(f3.c_str(), "wb");
    CheckFile(fp3, f3.c_str());
    fpfull = fopen(ffull.c_str(), "wb");
    CheckFile(fpfull, ffull.c_str());
    fp3_budget = fopen(f3b.c_str(), "wb");
    CheckFile(fp3_budget, f3b.c_str());
    fpfull_budget = fopen(ffullb.c_str(), "wb");
    CheckFile(fpfull_budget, ffullb.c_str());
    fpbasin = fopen(fbasin.c_str(), "wb");
    CheckFile(fpbasin, fbasin.c_str());

    writeDatHeader(fp3, "elewb3_resid", md->NumEle, icol);
    writeDatHeader(fpfull, "elewbfull_resid", md->NumEle, icol);
    writeDatHeader(fp3_budget, "elewb3_budget_resid", md->NumEle, icol);
    writeDatHeader(fpfull_budget, "elewbfull_budget_resid", md->NumEle, icol);
    writeDatHeader(fpbasin, "basinwbfull (m3)", 9, basin_icol);

    initBasinOutlets();
    basin_s_ele_prev_m3 = basinElementStorageFull_m3();
    basin_s_riv_prev_m3 = basinRiverStorage_m3();
    basin_p_acc_m3 = 0.0;
    basin_et_acc_m3 = 0.0;
    basin_qout_acc_m3 = 0.0;
    basin_qedge_acc_m3 = 0.0;
    basin_qbc_acc_m3 = 0.0;
    basin_qss_acc_m3 = 0.0;
    basin_noncons_edge_acc_m3 = 0.0;

    last_t = md->CS.StartTime;
    has_prev = 1;
    last_written_floor = -1;
}

void WaterBalanceDiag::closeFiles()
{
    if (fp3 != nullptr) {
        fclose(fp3);
        fp3 = nullptr;
    }
    if (fpfull != nullptr) {
        fclose(fpfull);
        fpfull = nullptr;
    }
    if (fp3_budget != nullptr) {
        fclose(fp3_budget);
        fp3_budget = nullptr;
    }
    if (fpfull_budget != nullptr) {
        fclose(fpfull_budget);
        fpfull_budget = nullptr;
    }
    if (fpbasin != nullptr) {
        fclose(fpbasin);
        fpbasin = nullptr;
    }
}

void WaterBalanceDiag::onETUpdate()
{
    if (!enabled || md == nullptr) {
        return;
    }
    const int n = md->NumEle;
    for (int i = 0; i < n; i++) {
        ic_raw[(size_t)i] = md->qEleE_IC[i];
    }
}

void WaterBalanceDiag::sample(double t_min, const double *DY)
{
    if (!enabled || md == nullptr || DY == nullptr) {
        return;
    }

    const int n = md->NumEle;

    const double dt_min = t_min - last_t;
    if (!(dt_min > 0.0) || !has_prev) {
        last_t = t_min;
        has_prev = 1;
        return;
    }

    for (int i = 0; i < n; i++) {
        const int isf = i;
        const int ius = i + n;
        const int igw = i + 2 * n;
        const double sy = md->Ele[i].Sy;
        const double area = md->Ele[i].area;

        const double net3 = DY[isf] + sy * DY[ius] + sy * DY[igw];
        const double netfull = net3 + (md->qElePrep[i] - md->qEleNetPrep[i]) - ic_raw[(size_t)i];

        /* Backward Euler integration on accepted solution samples */
        flux3_acc[(size_t)i] += net3 * dt_min;
        fluxfull_acc[(size_t)i] += netfull * dt_min;

        const double et3 = md->qEs[i] + md->qEu[i] + md->qEg[i] + md->qTu[i] + md->qTg[i];
        const double qlat3 = (md->QeleSurfTot[i] + md->QeleSubTot[i]) / area;
        const double qbc = (md->Ele[i].iBC < 0) ? (md->Ele[i].QBC / area) : 0.0;
        const double qss = (md->Ele[i].iSS != 0) ? (md->Ele[i].QSS / area) : 0.0;

        const double net3_budget = md->qEleNetPrep[i] - et3 - qlat3 + qbc + qss;
        const double netfull_budget = md->qElePrep[i] - (ic_raw[(size_t)i] + et3) - qlat3 + qbc + qss;
        budget3_acc[(size_t)i] += net3_budget * dt_min;
        budgetfull_acc[(size_t)i] += netfull_budget * dt_min;
    }

    {
        double p_rate_m3min = 0.0;
        double et_rate_m3min = 0.0;
        double qbc_rate_m3min = 0.0;
        double qss_rate_m3min = 0.0;
        double noncons_edge_rate_m3min = 0.0;

        for (int i = 0; i < n; i++) {
            const double area = md->Ele[i].area;
            p_rate_m3min += md->qElePrep[i] * area;
            et_rate_m3min +=
                (ic_raw[(size_t)i] + md->qEs[i] + md->qEu[i] + md->qEg[i] + md->qTu[i] + md->qTg[i]) * area;
            if (md->Ele[i].iBC < 0) {
                qbc_rate_m3min += md->Ele[i].QBC;
            }
            if (md->Ele[i].iSS != 0) {
                qss_rate_m3min += md->Ele[i].QSS;
            }
            for (int j = 0; j < 3; j++) {
                const int inabr = md->Ele[i].nabr[j] - 1;
                if (inabr >= 0) {
                    noncons_edge_rate_m3min += md->QeleSurf[i][j] + md->QeleSub[i][j];
                }
            }
        }
        for (int i = 0; i < md->NumRiv; i++) {
            if (md->Riv[i].BC < 0) {
                qbc_rate_m3min += md->Riv[i].qBC;
            }
        }

        const double qout_rate_m3min = basinOutletDischarge_m3min();
        const double qedge_rate_m3min = basinBoundaryEdgeOutflow_m3min();

        basin_p_acc_m3 += p_rate_m3min * dt_min;
        basin_et_acc_m3 += et_rate_m3min * dt_min;
        basin_qout_acc_m3 += qout_rate_m3min * dt_min;
        basin_qedge_acc_m3 += qedge_rate_m3min * dt_min;
        basin_qbc_acc_m3 += qbc_rate_m3min * dt_min;
        basin_qss_acc_m3 += qss_rate_m3min * dt_min;
        basin_noncons_edge_acc_m3 += noncons_edge_rate_m3min * dt_min;
    }

    last_t = t_min;
}

void WaterBalanceDiag::maybeWrite(double t_min)
{
    if (!enabled || md == nullptr) {
        return;
    }

    if (interval_min <= 0) {
        return;
    }

    const long long t_floor = (long long)floor(t_min);
    if (t_floor < interval_min) {
        return;
    }
    if ((t_floor % interval_min) != 0) {
        return;
    }
    if (t_floor == last_written_floor) {
        return;
    }

    const int n = md->NumEle;
    for (int i = 0; i < n; i++) {
        const double sy = md->Ele[i].Sy;
        const double s3 = md->yEleSurf[i] + sy * md->yEleUnsat[i] + sy * md->yEleGW[i];
        const double sfull = s3 + md->yEleSnow[i] + md->yEleIS[i];

        resid3[(size_t)i] = s3 - s3_prev[(size_t)i] - flux3_acc[(size_t)i];
        residfull[(size_t)i] = sfull - sfull_prev[(size_t)i] - fluxfull_acc[(size_t)i];
        resid3_budget[(size_t)i] = s3 - s3_prev[(size_t)i] - budget3_acc[(size_t)i];
        residfull_budget[(size_t)i] = sfull - sfull_prev[(size_t)i] - budgetfull_acc[(size_t)i];

        s3_prev[(size_t)i] = s3;
        sfull_prev[(size_t)i] = sfull;
        flux3_acc[(size_t)i] = 0.0;
        fluxfull_acc[(size_t)i] = 0.0;
        budget3_acc[(size_t)i] = 0.0;
        budgetfull_acc[(size_t)i] = 0.0;
    }

    const double t_quantized = (double)(t_floor - (long long)interval_min);
    writeDatRecord(fp3, t_quantized, md->NumEle, resid3);
    writeDatRecord(fpfull, t_quantized, md->NumEle, residfull);
    writeDatRecord(fp3_budget, t_quantized, md->NumEle, resid3_budget);
    writeDatRecord(fpfull_budget, t_quantized, md->NumEle, residfull_budget);

    const double s_ele_m3 = basinElementStorageFull_m3();
    const double s_riv_m3 = basinRiverStorage_m3();
    const double ds_ele_m3 = s_ele_m3 - basin_s_ele_prev_m3;
    const double ds_riv_m3 = s_riv_m3 - basin_s_riv_prev_m3;
    const double ds_total_m3 = ds_ele_m3 + ds_riv_m3;
    const double net_in_minus_out_m3 =
        basin_p_acc_m3 + basin_qbc_acc_m3 + basin_qss_acc_m3 - basin_et_acc_m3 - basin_qout_acc_m3 - basin_qedge_acc_m3;
    const double resid_m3 = ds_total_m3 - net_in_minus_out_m3;

    basin_values[0] = ds_total_m3;
    basin_values[1] = basin_p_acc_m3;
    basin_values[2] = basin_et_acc_m3;
    basin_values[3] = basin_qout_acc_m3;
    basin_values[4] = basin_qedge_acc_m3;
    basin_values[5] = basin_qbc_acc_m3;
    basin_values[6] = basin_qss_acc_m3;
    basin_values[7] = basin_noncons_edge_acc_m3;
    basin_values[8] = resid_m3;
    writeDatRecord(fpbasin, t_quantized, 9, basin_values);

    basin_s_ele_prev_m3 = s_ele_m3;
    basin_s_riv_prev_m3 = s_riv_m3;
    basin_p_acc_m3 = 0.0;
    basin_et_acc_m3 = 0.0;
    basin_qout_acc_m3 = 0.0;
    basin_qedge_acc_m3 = 0.0;
    basin_qbc_acc_m3 = 0.0;
    basin_qss_acc_m3 = 0.0;
    basin_noncons_edge_acc_m3 = 0.0;

    last_written_floor = t_floor;
}
