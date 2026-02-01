//
//  MD_ET.cpp
//  SHUD
//
//  Created by Lele Shu on 10/26/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Model_Data.hpp"

#include <algorithm>
#include <cmath>

void Model_Data::updateforcing(double t){
    int i;
    for(i = 0; i < NumEle; i++){
        Ele[i].updateElement(uYsf[i], uYus[i], uYgw[i]);
        tReadForcing(t,i);
    }
}
void Model_Data::tReadForcing(double t, int i){
    int idx = Ele[i].iForc - 1;
    double etp, ra, rs, t0, hc, U2, Uz, Zmeasure, lai;
    double GroundHeatFlux, RG;
    t_prcp[i] = tsd_weather[idx].getX(t, i_prcp) * gc.cPrep;
    t0= tsd_weather[idx].getX(t, i_temp);
    t_temp[i] = TemperatureOnElevation(t0, Ele[i].z_surf, tsd_weather[idx].xyz[2]) +  gc.cTemp;
    t_lai[i] = tsd_LAI.getX(t, Ele[i].iLC) * gc.cLAItsd ;
    lai = t_lai[i];
    t_mf[i] = tsd_MF.getX(t, Ele[i].iMF) * gc.cMF / 1440.;  /*  [m/day/C] to [m/min/C].
                                                            1.6 ~ 6.0 mm/day/C is typical value in USDA book
                                                            Input is 1.4 ~ 3.0 mm/d/c */
    /* Shortwave radiation forcing (W/m^2) */
    const double dswrf_h = tsd_weather[idx].getX(t, i_rn);
    double dswrf_t = dswrf_h;
    double factor = 1.0;
    if (CS.terrain_radiation) {
        if (CS.tsr_factor_mode == TSR_FORCING_INTERVAL) {
            const double t0 = tsd_weather[idx].currentTimeMin();
            double t1 = tsd_weather[idx].nextTimeMin();
            if (!std::isfinite(t0)) {
                factor = 0.0;
            } else {
                if (!std::isfinite(t1) || !(t1 > t0)) {
                    // Fall back to a single-step interval if forcing time is not available.
                    t1 = t0 + CS.SolverStep;
                }

                /* Bucket key: forcing interval left endpoint (minute) */
                const long long bucket = (long long)llround(t0);

                /* Precompute solar samples for this forcing interval (shared across elements). */
                int dt_int_min = CS.tsr_integration_step_min;
                if (dt_int_min <= 0) {
                    dt_int_min = 60;
                }

                if (bucket != tsr_forcing_bucket || t0 != tsr_forcing_t0 || t1 != tsr_forcing_t1 ||
                    dt_int_min != tsr_forcing_dt_int_min) {
                    tsr_forcing_bucket = bucket;
                    tsr_forcing_t0 = t0;
                    tsr_forcing_t1 = t1;
                    tsr_forcing_dt_int_min = dt_int_min;

                    const double dt_forc = t1 - t0;
                    double dt_int = (double)dt_int_min;
                    if (dt_int > dt_forc) {
                        dt_int = dt_forc;
                    }

                    int n = (int)std::ceil(dt_forc / dt_int);
                    if (n < 1) {
                        n = 1;
                    }
                    tsr_forcing_n = n;
                    const double dt_seg = dt_forc / (double)n;

                    tsr_forcing_sx.assign((size_t)n, 0.0);
                    tsr_forcing_sy.assign((size_t)n, 0.0);
                    tsr_forcing_sz.assign((size_t)n, 0.0);
                    tsr_forcing_wdt.assign((size_t)n, 0.0);
                    tsr_forcing_den = 0.0;

                    for (int k = 0; k < n; k++) {
                        const double tk = t0 + (k + 0.5) * dt_seg;

                        // forcing 时间为 UTC；必须显式传入 timezone_hours=0.0
                        // 避免 solarPosition() 根据 lon 推算本地时区（round(lon/15)）导致相位偏移
                        const SolarPosition sp = solarPosition(tk, CS.solar_lat_deg, CS.solar_lon_deg, Time, 0.0);
                        const double cosz = sp.cosZ;
                        if (!(cosz > 0.0) || !std::isfinite(cosz) || !std::isfinite(sp.azimuth)) {
                            continue;
                        }

                        const double cosz_clamped = std::min(1.0, std::max(-1.0, cosz));
                        const double sinz = std::sqrt(std::max(0.0, 1.0 - cosz_clamped * cosz_clamped));
                        const double sin_az = std::sin(sp.azimuth);
                        const double cos_az = std::cos(sp.azimuth);

                        const double sx = sinz * sin_az;
                        const double sy = sinz * cos_az;
                        const double sz = cosz_clamped;

                        const double wdt = std::max(0.0, cosz_clamped) * dt_seg; // weight = max(cosZ,0)
                        if (!(wdt > 0.0) || !std::isfinite(wdt)) {
                            continue;
                        }

                        tsr_forcing_sx[(size_t)k] = sx;
                        tsr_forcing_sy[(size_t)k] = sy;
                        tsr_forcing_sz[(size_t)k] = sz;
                        tsr_forcing_wdt[(size_t)k] = wdt;
                        tsr_forcing_den += wdt;
                    }
                }

                if (tsr_factor_bucket != nullptr && tsr_factor != nullptr) {
                    if (tsr_factor_bucket[i] != bucket) {
                        double num = 0.0;
                        const double cap = CS.rad_factor_cap;
                        const double cosz_min = CS.rad_cosz_min;

                        if (tsr_forcing_den > 0.0 && tsr_forcing_n > 0) {
                            const int n = tsr_forcing_n;
                            const double nx = Ele[i].nx;
                            const double ny = Ele[i].ny;
                            const double nz = Ele[i].nz;
                            for (int k = 0; k < n; k++) {
                                const double wdt = tsr_forcing_wdt[(size_t)k];
                                if (!(wdt > 0.0)) {
                                    continue;
                                }
                                const double sx = tsr_forcing_sx[(size_t)k];
                                const double sy = tsr_forcing_sy[(size_t)k];
                                const double sz = tsr_forcing_sz[(size_t)k];

                                const double cosi = nx * sx + ny * sy + nz * sz;
                                if (!(cosi > 0.0) || !std::isfinite(cosi)) {
                                    continue;
                                }

                                double denom = sz;
                                if (denom < cosz_min) {
                                    denom = cosz_min;
                                }
                                if (!(denom > 0.0) || !std::isfinite(denom)) {
                                    continue;
                                }

                                double fk = cosi / denom;
                                if (!std::isfinite(fk) || !(fk > 0.0)) {
                                    continue;
                                }
                                if (fk > cap) {
                                    fk = cap;
                                }
                                num += wdt * fk;
                            }
                        }

                        double feff = 0.0;
                        if (tsr_forcing_den > 0.0) {
                            feff = num / tsr_forcing_den;
                            if (!std::isfinite(feff) || !(feff > 0.0)) {
                                feff = 0.0;
                            }
                            if (feff > CS.rad_factor_cap) {
                                feff = CS.rad_factor_cap;
                            }
                        }

                        tsr_factor[i] = feff;
                        tsr_factor_bucket[i] = bucket;
                    }
                    factor = tsr_factor[i];
                }
            }
        } else { /* TSR_INSTANT (default) */
            int interval_min = CS.solar_update_interval;
            if (interval_min <= 0) {
                interval_min = 60;
            }

            /* Timestamp alignment: left endpoint of SOLAR_UPDATE_INTERVAL bucket */
            const double bucket_eps = 1.0e-6; /* [min] */
            const long long bucket = (long long)floor((t + bucket_eps) / (double)interval_min);
            const double t_aligned = (double)bucket * (double)interval_min;
            if (bucket != tsr_solar_bucket) {
                tsr_solar_bucket = bucket;
                tsr_solar_t_aligned = t_aligned;
                /*
                 * TSR 太阳位置近似（Terrain Shortwave Radiation）
                 * - 使用单个全局太阳位置：来自 B2a 的 solar_lon_deg/solar_lat_deg（CS.solar_lon_deg/lat_deg）
                 * - 物理假设：流域内太阳位置（高度角/方位角）的空间差异可忽略
                 * - 适用范围：推荐用于特征长度 <200km 的流域
                 * - 空间误差估计：<50km 优秀(<1%)；50-200km 可接受(<5%)；>200km 需改进(10-20%)
                 */
                // forcing 时间为 UTC；必须显式传入 timezone_hours=0.0
                // 避免 solarPosition() 根据 lon 推算本地时区（round(lon/15)）导致相位偏移
                tsr_solar_pos = solarPosition(t_aligned, CS.solar_lat_deg, CS.solar_lon_deg, Time, 0.0);
            }

            if (tsr_factor_bucket != nullptr && tsr_factor != nullptr) {
                if (tsr_factor_bucket[i] != bucket) {
                    tsr_factor[i] = terrainFactor(Ele[i].nx,
                                                  Ele[i].ny,
                                                  Ele[i].nz,
                                                  tsr_solar_pos,
                                                  CS.rad_factor_cap,
                                                  CS.rad_cosz_min);
                    tsr_factor_bucket[i] = bucket;
                }
                factor = tsr_factor[i];
            }
        }

        dswrf_t = dswrf_h * factor;
    }
    if (ele_rn_h_wm2 != nullptr) ele_rn_h_wm2[i] = dswrf_h;
    if (ele_rn_t_wm2 != nullptr) ele_rn_t_wm2[i] = dswrf_t;
    if (ele_rn_factor != nullptr) ele_rn_factor[i] = factor;
    if (CS.radiation_input_mode == SWNET) {
        // SWNET mode: forcing 第6列已是净短波，不再乘 (1-Albedo)
        t_rn[i] = dswrf_t;
    } else {
        // SWDOWN mode (default): forcing 第6列是下行短波，需乘 (1-Albedo) 净化
        t_rn[i] = dswrf_t * (1 - Ele[i].Albedo);
    }
    Uz = t_wind[i] = (fabs(tsd_weather[idx].getX(t, i_wind) ) + 0.001); // +.001 voids ZERO.
    t_rh[i] = tsd_weather[idx].getX(t, i_rh);
//    t_hc[i] = tsd_RL.getX(t, Ele[i].iLC);
//    t_hc[i] = max(t_hc[i], CONSt_hc);
    /* Precipitation  */
    t_prcp[i]   = t_prcp[i] * 0.001 / 1440. ; // [mm d-1] to [m min-1]
    /* Potential ET */
    t_rn[i]     = t_rn[i] * 1.0e-6;  // [W m-2] to [MJ m-2 s-1]
    /*
     t_wind [m s-1] ;
     t_temp  [C] ;
     t_rh  [0-1] ;
     */
    t_rh[i]     = min(max(t_rh[i], CONST_RH), 1.0); // [value is b/w 0~1 ]
    
    qElePrep[i] = t_prcp[i];
    double lambda = LatentHeat(t_temp[i]);                      // eq 4.2.1  [MJ/kg]
    double Gamma = PsychrometricConstant(Ele[i].FixPressure, lambda); // eq 4.2.28  [kPa C-1]
    double es = VaporPressure_Sat(t_temp[i]);                   // eq 4.2.2 [kpa]
    double ea = es * t_rh[i];   // [kPa]
    double ed = es - ea ;  // [kPa]
    double Delta = SlopeSatVaporPressure(t_temp[i], es);        // eq 4.2.3 [kPa C-1]
    double rho = AirDensity(Ele[i].FixPressure, t_temp[i]);;    // eq 4.2.4 [kg m-3]
    /* R - G in the PM equation.*/
    if(Ele[i].iLake > 0 ){
        GroundHeatFlux = 0.;
        RG = t_rn[i];
    }else{
        if(lai > 0){
            GroundHeatFlux = 0.4 * exp(-0.5 * lai) * t_rn[i];
        }else{
            GroundHeatFlux = 0.1 * t_rn[i];
        }
    }
    RG = t_rn[i] - GroundHeatFlux;
    U2 = WindProfile(2.0, t_wind[i], Ele[i].windH, 0., ROUGHNESS_WATER); // [m s-1]
    qPotEvap[i] = gc.cETP * PET_PM_openwater(Delta, Gamma, lambda, RG, ed, U2) * 60.; // eq 4.2.30
    if(Ele[i].iLake > 0){        /* Open-water */
        qPotTran[i] = gc.cETP * 0.;
        etp = qPotEvap[i];
    }else if(lai <= 0.){        /* Bare soiln */
        qPotTran[i] = gc.cETP * 0.;
        etp = qPotEvap[i];
    }else{
//        hc = lai2hc(lai);
        hc = lai * 0.5;
        Zmeasure = hc*1.3333; /* When hc > Zm, Zm = hc + 5.0m */
        ra = AerodynamicResistance(Uz, hc, Zmeasure, Zmeasure); // eq 4.2.25  [s m-1]
//        if( Zmeasure > HeightWindMeasure){ /* Veg Height > Wind Measure Height */
//            ra = AerodynamicResistance(Uz, hc, Zmeasure, Zmeasure); // eq 4.2.25  [s m-1]
//        }else{
//            ra = AerodynamicResistance(Uz, hc, HeightWindMeasure, 2.0); // eq 4.2.25  [s m-1]
//        }
//        ra = min(300., ra);
//        if(ra < 0){
//            ra = AerodynamicResistance(Uz, hc, Zmeasure, Zmeasure); // eq 4.2.25  [s m-1]
//            ra = AerodynamicResistance(Uz, hc, HeightWindMeasure, 2.0); // eq 4.2.25  [s m-1]
//        }
        CheckNonZero(ra, i, "Aerodynamic Resistance");
//        CheckNANi(ra, i, "Aerodynamic Resistance");
        rs = BulkSurfaceResistance(lai);  // eq 4.2.22 & 4.2.25  [s m-1]
        qPotTran[i] = gc.cETP * PET_Penman_Monteith(RG, rho, ed, Delta, ra, rs, Gamma, lambda) * 60.;// eq 4.2.27
        etp = qPotTran[i] * Ele[i].VegFrac + qPotEvap[i] * (1. - Ele[i].VegFrac);
        CheckNANi(qPotTran[i], i, "qPotTran[i]");
    }
    qEleETP[i] = etp;
}
void Model_Data::ET(double t, double tnext){
    double  T=NA_VALUE,  LAI=NA_VALUE, MF =NA_VALUE, prcp = NA_VALUE;
    double  snFrac, snAcc, snMelt, snStg;
    double  icAcc, icEvap, icStg, icMax, vgFrac;
    double  DT_min = tnext - t;
    double  ta_surf, ta_sub;
    int i;
#ifdef _OPENMP_ON
#pragma omp for
#endif
    for(i = 0; i < NumEle; i++) {
        T = t_temp[i];
        prcp = t_prcp[i];
        /* Snow Accumulation */
        MF = t_mf[i];
        snStg = yEleSnow[i];
        /* Snow Accumulation/Melt Calculation*/
        snFrac  = FrozenFraction(T, Train, Tsnow);
        
        if(CS.cryosphere){
            AccT_surf[i].push(T, t);
            AccT_sub[i].push(T, t);
            ta_surf = AccT_surf[i].getACC();
            ta_sub  = AccT_sub[i].getACC();
            fu_Sub[i] = 1. - FrozenFraction(ta_sub, AccT_sub_max, AccT_sub_min);
            fu_Surf[i] = 1. - FrozenFraction(ta_surf, AccT_surf_max, AccT_surf_min);
        }else{
            fu_Sub[i] = 1.;
            fu_Surf[i] = 1.;
        }
        
        snAcc = snFrac * prcp;
        snMelt = (T > To ? (T - To) * MF : 0.);    /* eq. 7.3.14 in Maidment */
        snMelt = min(max(0., snStg / DT_min), max(0., snMelt));
//        CheckNonNegative(snMelt, i, "Snow Melting");
        snStg += (snAcc - snMelt) * DT_min;
        
        /* Interception */
        LAI = t_lai[i];
        vgFrac = Ele[i].VegFrac;
        icStg = (vgFrac > ZERO) ? (yEleIS[i] / vgFrac) : 0.0;
        if(LAI > ZERO){
            icMax = gc.cISmax * IC_MAX * LAI;
            icAcc = min(prcp - snAcc, max(0., (icMax - icStg) / DT_min) );
            icEvap = min(max(0., icStg / DT_min), qPotEvap[i]);
        }else{
            icAcc = 0.;
            icEvap = 0.;
        }
        icStg += (icAcc - icEvap) * DT_min;
        
        /* Update the storage value and net precipitaion */
        yEleIS[i] = icStg * vgFrac;
        yEleSnow[i] = snStg;
        qEleE_IC[i] = icEvap * vgFrac;
        qEleNetPrep[i] = (1. - snFrac) * prcp + snMelt - icAcc * vgFrac ;

//        CheckNonNegative(qEleNetPrep[i], i, "Net Precipitation");
//        CheckNonNegative(qEleE_IC[i], i, "qEleE_IC");
    }
}
void Model_Data::f_etFlux(int i, double t){
    double Es = 0., Eu = 0., Tu = 0., Eg = 0., Tg = 0.;
    double va = Ele[i].VegFrac, vb = 1. - Ele[i].VegFrac;
    double pj = 1. - Ele[i].ImpAF;
    iBeta[i] = SoilMoistureStress(Soil[(Ele[i].iSoil - 1)].ThetaS, Soil[(Ele[i].iSoil - 1)].ThetaR, Ele[i].u_satn);
    /* Evaporation from SURFACE ponding water */
    Es = min(max(0., uYsf[i]), qPotEvap[i]) * vb;
    if(Es < qPotEvap[i]){
        /* Some PET is extracted by surface Evaporation, so PET - Es is the effective PET now. */
        if(uYgw[i] > Ele[i].WetlandLevel){
            /* Evporation from GroundWater, ONLY when gw above wetland level*/
            Eg = min(max(0., uYgw[i]), qPotEvap[i] - Es) * pj * vb;
            Eu = 0.;
        }else{
            Eg = 0.;
            /* Evaporation from Unsaturated Zone. */
            Eu = min(max(0., uYus[i]), iBeta[i] * (qPotEvap[i] - Es)) * pj * vb;
        }
    }else{
        /* All evporation is from land surface ONLY */
        Eg = 0.;
        Eu = 0.;
    }
    /* Vegetation Transpiration */
    if(t_lai[i] > ZERO){
        if(qEleE_IC[i] >= qPotTran[i]){
            Tg = Tu = 0.;
            qEleE_IC[i] = qPotTran[i] * pj * va;
        }else{
            if(uYgw[i] > Ele[i].RootReachLevel){
                Tg = min(max(0., uYgw[i]), (qPotTran[i] - qEleE_IC[i]) ) * pj * va;
                Tu = 0.;
            }else{
                Tg = 0.;
                Tu = min(max(0., uYus[i]), iBeta[i] * (qPotTran[i] - qEleE_IC[i]) ) * pj * va;
            }
        }
    }else{
        Tg = Tu = qEleE_IC[i] = 0.;
    }
    qEs[i] = Es;
    qEu[i] = Eu;
    qEg[i] = Eg;
    qTu[i] = Tu;
    qTg[i] = Tg;
    qEleTrans[i] = Tg + Tu;
    qEleEvapo[i] = Eu + Eg + Es;  
    qEleETA[i] = qEleE_IC[i] + qEleEvapo[i] + qEleTrans[i];
    if(qEleETA[i] > qEleETP[i] * 2.){
        printf("Warning: More AET(%.3E) than PET(%.3E) on Element (%d).", qEleETA[i], qEleETP[i], i+1);
    }
    CheckNonNegative(Es, i, "Es"); // Debug Only
    CheckNonNegative(Eu, i, "Eu");
    CheckNonNegative(Eg, i, "Eg");
    CheckNonNegative(Tu, i, "Tu");
    CheckNonNegative(Tg, i, "Tg");
    CheckNANi(qEleETA[i], i, "Potential ET (Model_Data::EvapoTranspiration)");
    CheckNANi(qEleEvapo[i], i, "Transpiration (Model_Data::EvapoTranspiration)");
    CheckNANi(qEleTrans[i], i, "Soil Evaporation (Model_Data::EvapoTranspiration)");
#ifdef DEBUG
#endif
}
