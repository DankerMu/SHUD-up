//  Model_Control.cpp
//
//  Created by Lele Shu on 7/17/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Model_Control.hpp"

#include <cctype>
#include <cmath>

static bool _read_kv_cfg_value(const char *path, const char *key, char *out, size_t out_sz)
{
    if (out == NULL || out_sz == 0) {
        return false;
    }
    out[0] = '\0';

    if (path == NULL || path[0] == '\0' || key == NULL || key[0] == '\0') {
        return false;
    }

    FILE *fp = fopen(path, "r");
    if (fp == NULL) {
        return false;
    }

    char line[MAXLEN];
    while (fgets(line, (int)sizeof(line), fp) != NULL) {
        char *p = line;
        while (*p != '\0' && std::isspace((unsigned char)*p)) {
            p++;
        }
        if (*p == '\0' || *p == '\n' || *p == '#') {
            continue;
        }

        char k[MAXLEN] = "";
        char v[MAXLEN] = "";
        if (sscanf(p, "%s %s", k, v) != 2) {
            continue;
        }
        if (strcasecmp(k, key) == 0) {
            strncpy(out, v, out_sz - 1);
            out[out_sz - 1] = '\0';
            fclose(fp);
            return true;
        }
    }

    fclose(fp);
    return false;
}
void PrintOutDt::defaultmode(){
    int dt = 1440;
    /* Element storage */
    dt_ye_gw = dt;
    dt_ye_surf = dt;
    dt_ye_snow = dt;
    dt_ye_ic = 0;
    dt_ye_unsat = dt;
    
    /* Element Fluxes */
    dt_qe_prcp = dt; // default output PRCP.
    dt_qe_infil = dt;
    dt_qe_et = dt;
    
    dt_qe_rech = dt;
    dt_qe_etp = dt;
    dt_qe_eta = dt;
    
    /* Element volume Fluxes */
    dt_Qe_sub = 0;
    dt_Qe_surf = 0;
    
    /* River Stage */
    dt_yr_stage = dt;
    /* River volume Fluxes */
    dt_Qr_up = dt;
    dt_Qr_down = dt;
    dt_Qr_sub = dt;
    dt_Qr_surf = dt;
    
    /* Lake  stage */
    dt_lake     = dt;
}
void PrintOutDt::calibmode(int dt ){
    /* Element storage */
    dt_ye_gw = dt;
    dt_ye_surf = dt;
    dt_ye_snow = 0;
    dt_ye_ic = 0;
    dt_ye_unsat = 0;
    
    /* Element Fluxes */
    dt_qe_prcp = dt; // default output PRCP.
    dt_qe_infil = 0;
    dt_qe_et = dt;
    
    dt_qe_rech = 0;
    dt_qe_etp = dt;
    dt_qe_eta = dt;
    
    /* Element volume Fluxes */
    dt_Qe_sub = 0;
    dt_Qe_surf = 0;
    
    /* River Stage */
    dt_yr_stage = dt;
    /* River volume Fluxes */
    dt_Qr_down = dt;
    dt_Qr_sub = dt;
    dt_Qr_surf = dt;
    
    /* Lake  stage */
    dt_lake     = dt;
}
Control_Data::Control_Data(){
}
Control_Data::~Control_Data(){
//    delete Tout;
}
void Control_Data::ExportResults(double t){
    for (int i = 0; i < NumPrint; i++){
        PCtrl[i].PrintData(dt, t);
    }
}

void Control_Data::updateSimPeriod(){
    updateSimPeriod(DayStart, DayEnd);
}
void Control_Data::updateSimPeriod(double day0, double day1){
    DayStart = day0;
    DayEnd = day1;
    StartTime = DayStart  * 1440;
    EndTime = (DayEnd) * 1440;
    NumSteps = (unsigned long) (EndTime - StartTime) / SolverStep;
}


void Control_Data::read(const char *fn){
    char    str[MAXLEN];
    char    optstr[MAXLEN];
    bool radiation_input_mode_user_set = false;
    bool radiation_input_mode_from_forcing_cfg = false;
    char forcing_radiation_kind_eff[MAXLEN] = "";
#ifdef _OPENMP_ON
    num_threads = omp_get_max_threads();; /*Default number of threads for OpenMP*/
#else
    num_threads = 0;
#endif
    
    FILE *fp;
    fp = fopen(fn, "r");
    CheckFile(fp, fn);
    /* Read through parameter file to find parameters */
    double val;
    while (fgets(str, MAXLEN, fp)) {
        if (str[0] == '#' || str[0] == '\n' || str[0] == '\0' || str[0] == ' ')
        {
            continue;
        }
        sscanf (str, "%s %lf", optstr, &val);
        /* Get Model Parameters */
        if (strcasecmp ("ASCII_OUTPUT", optstr) == 0)
            Ascii=  (int) val;
        else if (strcasecmp ("BINARY_OUTPUT", optstr) == 0)
            Binary = (int)  val;
        else if (strcasecmp ("SpinupDay", optstr) == 0)
            Spinup =  (int) val;
        else if (strcasecmp ("SCR_INTV", optstr) == 0)
            screenIntv = (int)  val;
        else if (strcasecmp ("VERBOSE", optstr) == 0)
            Verbose = (int)  val;
        else if (strcasecmp ("CloseBoundary", optstr) == 0)
            CloseBoundary = (int)  val;
        else if (strcasecmp ("INIT_MODE", optstr) == 0)
            init_type = (int)  val;
        else if (strcasecmp ("NUM_OPENMP", optstr) == 0)
            num_threads = (int)  val;
        else if (strcasecmp ("ABSTOL", optstr) == 0)
            abstol =  val;
        else if (strcasecmp ("RELTOL", optstr) == 0)
            reltol =  val;
        else if (strcasecmp ("INIT_SOLVER_STEP", optstr) == 0)
            InitStep =  val;
        else if (strcasecmp ("MAX_SOLVER_STEP", optstr) == 0)
            MaxStep =  val;
        else if (strcasecmp ("Update_IC_STEP", optstr) == 0)
            UpdateICStep =  val;
        else if (strcasecmp ("ET_Mode", optstr) == 0)
            ET_Mode =  val;
        else if (strcasecmp("FORCING_MODE", optstr) == 0) {
            const ForcingMode default_mode = FORCING_CSV;
            char mode_str[MAXLEN] = "";
            forcing_mode = default_mode;
            if (sscanf(str, "%s %s", optstr, mode_str) != 2) {
                fprintf(stderr,
                        "WARNING: FORCING_MODE missing value in %s; using default CSV (%d).\n",
                        fn,
                        default_mode);
            } else if (strcasecmp(mode_str, "CSV") == 0) {
                forcing_mode = FORCING_CSV;
            } else if (strcasecmp(mode_str, "NETCDF") == 0) {
                forcing_mode = FORCING_NETCDF;
            } else {
                char *endptr = NULL;
                const double mode_val = strtod(mode_str, &endptr);
                if (endptr != NULL && *endptr == '\0' && (mode_val == 0.0 || mode_val == 1.0)) {
                    forcing_mode = (mode_val == 1.0) ? FORCING_NETCDF : FORCING_CSV;
                } else {
                    fprintf(stderr,
                            "ERROR: invalid FORCING_MODE value '%s' in %s. Valid values: CSV/NETCDF or 0/1.\n",
                            mode_str,
                            fn);
                    myexit(ERRFileIO);
                }
            }
        }
        else if (strcasecmp("FORCING_CFG", optstr) == 0) {
            char cfg_str[MAXLEN] = "";
            if (sscanf(str, "%s %s", optstr, cfg_str) != 2) {
                fprintf(stderr, "ERROR: FORCING_CFG missing value in %s.\n", fn);
                myexit(ERRFileIO);
            }
            strncpy(forcing_cfg, cfg_str, MAXLEN - 1);
            forcing_cfg[MAXLEN - 1] = '\0';
        }
        else if (strcasecmp("OUTPUT_MODE", optstr) == 0) {
            const OutputMode default_mode = OUTPUT_LEGACY;
            char mode_str[MAXLEN] = "";
            output_mode = default_mode;
            if (sscanf(str, "%s %s", optstr, mode_str) != 2) {
                fprintf(stderr,
                        "WARNING: OUTPUT_MODE missing value in %s; using default LEGACY (%d).\n",
                        fn,
                        default_mode);
            } else if (strcasecmp(mode_str, "LEGACY") == 0) {
                output_mode = OUTPUT_LEGACY;
            } else if (strcasecmp(mode_str, "NETCDF") == 0) {
                output_mode = OUTPUT_NETCDF;
            } else if (strcasecmp(mode_str, "BOTH") == 0) {
                output_mode = OUTPUT_BOTH;
            } else {
                char *endptr = NULL;
                const double mode_val = strtod(mode_str, &endptr);
                if (endptr != NULL && *endptr == '\0' && (mode_val == 0.0 || mode_val == 1.0 || mode_val == 2.0)) {
                    output_mode = (mode_val == 1.0) ? OUTPUT_NETCDF : (mode_val == 2.0 ? OUTPUT_BOTH : OUTPUT_LEGACY);
                } else {
                    fprintf(stderr,
                            "ERROR: invalid OUTPUT_MODE value '%s' in %s. Valid values: LEGACY/NETCDF/BOTH or 0/1/2.\n",
                            mode_str,
                            fn);
                    myexit(ERRFileIO);
                }
            }
        }
        else if (strcasecmp("NCOUTPUT_CFG", optstr) == 0) {
            char cfg_str[MAXLEN] = "";
            if (sscanf(str, "%s %s", optstr, cfg_str) != 2) {
                fprintf(stderr, "ERROR: NCOUTPUT_CFG missing value in %s.\n", fn);
                myexit(ERRFileIO);
            }
            strncpy(ncoutput_cfg, cfg_str, MAXLEN - 1);
            ncoutput_cfg[MAXLEN - 1] = '\0';
        }
        else if (strcasecmp("RADIATION_INPUT_MODE", optstr) == 0) {
            const int default_mode = SWDOWN;
            char mode_str[MAXLEN] = "";
            radiation_input_mode = default_mode;
            radiation_input_mode_user_set = false;
            if (sscanf(str, "%s %s", optstr, mode_str) != 2) {
                fprintf(stderr,
                        "WARNING: RADIATION_INPUT_MODE missing value in %s; using default %s (%d).\n",
                        fn,
                        default_mode == SWNET ? "SWNET" : "SWDOWN",
                        default_mode);
            } else if (strcasecmp(mode_str, "SWDOWN") == 0) {
                radiation_input_mode = SWDOWN;
                radiation_input_mode_user_set = true;
            } else if (strcasecmp(mode_str, "SWNET") == 0) {
                radiation_input_mode = SWNET;
                radiation_input_mode_user_set = true;
            } else {
                char *endptr = NULL;
                const double mode_val = strtod(mode_str, &endptr);
                if (endptr != NULL && *endptr == '\0' && (mode_val == 0.0 || mode_val == 1.0)) {
                    radiation_input_mode = (mode_val == 1.0) ? SWNET : SWDOWN;
                    radiation_input_mode_user_set = true;
                } else {
                    fprintf(stderr,
                            "WARNING: invalid RADIATION_INPUT_MODE value '%s' in %s; using default %s (%d). "
                            "Valid values: SWDOWN/SWNET or 0/1.\n",
                            mode_str,
                            fn,
                            default_mode == SWNET ? "SWNET" : "SWDOWN",
                            default_mode);
                }
            }
        }
        else if (strcasecmp("SOLAR_LONLAT_MODE", optstr) == 0) {
            const SolarLonLatMode default_mode = FORCING_FIRST;
            char mode_str[MAXLEN] = "";
            solar_lonlat_mode = default_mode;
            if (sscanf(str, "%s %s", optstr, mode_str) != 2) {
                fprintf(stderr,
                        "WARNING: SOLAR_LONLAT_MODE missing value in %s; using default %s (%d).\n",
                        fn,
                        SolarLonLatModeName(default_mode),
                        default_mode);
            } else if (strcasecmp(mode_str, "FORCING_FIRST") == 0) {
                solar_lonlat_mode = FORCING_FIRST;
            } else if (strcasecmp(mode_str, "FORCING_MEAN") == 0) {
                solar_lonlat_mode = FORCING_MEAN;
            } else if (strcasecmp(mode_str, "FIXED") == 0) {
                solar_lonlat_mode = FIXED;
            } else {
                char *endptr = NULL;
                const double mode_val = strtod(mode_str, &endptr);
                if (endptr != NULL && *endptr == '\0' &&
                    (mode_val == 0.0 || mode_val == 1.0 || mode_val == 2.0)) {
                    solar_lonlat_mode = static_cast<SolarLonLatMode>(static_cast<int>(mode_val));
                } else {
                    fprintf(stderr,
                            "WARNING: invalid SOLAR_LONLAT_MODE value '%s' in %s; using default %s (%d). "
                            "Valid values: FORCING_FIRST/FORCING_MEAN/FIXED or 0/1/2.\n",
                            mode_str,
                            fn,
                            SolarLonLatModeName(default_mode),
                            default_mode);
                    solar_lonlat_mode = default_mode;
                }
            }
        }
        else if (strcasecmp("SOLAR_LON_DEG", optstr) == 0)
            solar_lon_deg_fixed = val;
        else if (strcasecmp("SOLAR_LAT_DEG", optstr) == 0)
            solar_lat_deg_fixed = val;
        else if (strcasecmp("TERRAIN_RADIATION", optstr) == 0) {
            const int flag = (int)val;
            if (flag == 0 || flag == 1) {
                terrain_radiation = flag;
            } else {
                fprintf(stderr,
                        "WARNING: invalid TERRAIN_RADIATION value %.3f in %s; using %d. Valid values: 0/1.\n",
                        val,
                        fn,
                        terrain_radiation);
            }
        }
        else if (strcasecmp("SOLAR_UPDATE_INTERVAL", optstr) == 0) {
            const int interval = (int)val;
            if (interval > 0) {
                /* Deprecated: SOLAR_UPDATE_INTERVAL used to control an "instant" TSR bucket update.
                 * TSR now always uses forcing-interval equivalent factor; reuse this value as the
                 * integration step for backward compatibility.
                 */
                tsr_integration_step_min = interval;
                fprintf(stderr,
                        "WARNING: SOLAR_UPDATE_INTERVAL is deprecated; treating it as TSR_INTEGRATION_STEP_MIN=%d (min).\n",
                        tsr_integration_step_min);
            } else {
                fprintf(stderr,
                        "WARNING: invalid SOLAR_UPDATE_INTERVAL value %.3f in %s; ignoring. Must be > 0.\n",
                        val,
                        fn);
            }
        }
        else if (strcasecmp("RAD_FACTOR_CAP", optstr) == 0) {
            if (std::isfinite(val) && val > 0.0) {
                rad_factor_cap = val;
            } else {
                fprintf(stderr,
                        "WARNING: invalid RAD_FACTOR_CAP value %.3f in %s; using %.3f. Must be finite and > 0.\n",
                        val,
                        fn,
                        rad_factor_cap);
            }
        }
        else if (strcasecmp("RAD_COSZ_MIN", optstr) == 0) {
            if (std::isfinite(val) && val >= 0.0) {
                rad_cosz_min = val;
                if (rad_cosz_min > 1.0) {
                    rad_cosz_min = 1.0;
                }
            } else {
                fprintf(stderr,
                        "WARNING: invalid RAD_COSZ_MIN value %.3f in %s; using %.3f. Must be finite and >= 0.\n",
                        val,
                        fn,
                        rad_cosz_min);
            }
        }
        else if (strcasecmp("TSR_FACTOR_MODE", optstr) == 0) {
            char mode_str[MAXLEN] = "";
            if (sscanf(str, "%s %s", optstr, mode_str) == 2) {
                if (strcasecmp(mode_str, "INSTANT") == 0 || strcmp(mode_str, "0") == 0) {
                    fprintf(stderr,
                            "WARNING: TSR_FACTOR_MODE=INSTANT is deprecated and no longer supported; using forcing-interval factor.\n");
                }
            } else {
                fprintf(stderr,
                        "WARNING: TSR_FACTOR_MODE is deprecated; ignoring (TSR uses forcing-interval factor).\n");
            }
        }
        else if (strcasecmp("TSR_INTEGRATION_STEP_MIN", optstr) == 0) {
            const int dt_int = (int)val;
            if (dt_int > 0) {
                tsr_integration_step_min = dt_int;
            } else {
                fprintf(stderr,
                        "WARNING: invalid TSR_INTEGRATION_STEP_MIN value %.3f in %s; using %d (min). Must be > 0.\n",
                        val,
                        fn,
                        tsr_integration_step_min);
            }
        }
        else if (strcasecmp ("ET_STEP", optstr) == 0 || strcasecmp ("LSM_STEP", optstr) == 0)
            ETStep =  val;
        else if (strcasecmp ("START", optstr) == 0)
            DayStart =  val;/* Convert days to Minutes */
        else if (strcasecmp ("END", optstr) == 0)
            DayEnd =  val;    /* Convert days to Minutes */
        else if (strcasecmp ("Exfiltration", optstr) == 0)
            exfiltration =  val;
        else if (strcasecmp ("cryosphere", optstr) == 0)
            cryosphere =  val;
//        else if (strcasecmp ("STEPSIZE_FACTOR", optstr) == 0)
//            a =  val;
//        else if (strcasecmp ("MODEL_STEPSIZE", optstr) == 0)
//            b =  val;
        
        /* Element print out control */
        // Y ele
        else if (strcasecmp ("dt_ye_ic", optstr) == 0)
            dt_ye_ic =  val;
        else if (strcasecmp ("dt_ye_SNOW", optstr) == 0)
            dt_ye_snow =  val;
        else if (strcasecmp ("dt_ye_SURF", optstr) == 0)
            dt_ye_surf =  val;
        else if (strcasecmp ("dt_ye_UNSAT", optstr) == 0)
            dt_ye_unsat =  val;
        else if (strcasecmp ("dt_ye_GW", optstr) == 0)
            dt_ye_gw =  val;
        // q Ele
        else if (strcasecmp ("dt_qe_PRCP", optstr) == 0)
            dt_qe_prcp =  val;
        else if (strcasecmp ("dt_qe_ET", optstr) == 0){
            dt_qe_et  =  val;
            dt_qe_etp =  val;
            dt_qe_eta =  val;
        }
        else if (strcasecmp ("dt_qe_rech", optstr) == 0)
            dt_qe_rech =  val;
        else if (strcasecmp ("dt_qe_infil", optstr) == 0)
            dt_qe_infil =  val;
        // Q Ele
        else if (strcasecmp ("dt_Qe_sub", optstr) == 0)
            dt_Qe_sub =  val;
        else if (strcasecmp ("dt_Qe_subx", optstr) == 0)
            dt_Qe_subx =  val;
        else if (strcasecmp ("dt_Qe_surf", optstr) == 0)
            dt_Qe_surf =  val;
        else if (strcasecmp ("dt_Qe_surfx", optstr) == 0)
            dt_Qe_surfx =  val;
        else if (strcasecmp ("dt_Qe_rsub", optstr) == 0)
            dt_Qe_rsub =  val;
        else if (strcasecmp ("dt_Qe_rsurf", optstr) == 0)
            dt_Qe_rsurf =  val;
        
        /* River print out control */
        // y Riv
        else if (strcasecmp ("dt_yr_stage", optstr) == 0)
            dt_yr_stage =  val;
        // Q riv
        else if (strcasecmp ("dt_Qr_Surf", optstr) == 0)
            dt_Qr_surf=  val;
        else if (strcasecmp ("dt_Qr_Sub", optstr) == 0)
            dt_Qr_sub =  val;
        else if (strcasecmp ("dt_Qr_down", optstr) == 0)
            dt_Qr_down =  val;
        else if (strcasecmp ("dt_Qr_up", optstr) == 0)
            dt_Qr_up =  val;
        /* Lake print out control */
        // Y lake
        else if (strcasecmp ("dt_lake", optstr) == 0)
            dt_lake =  val;
        /* Unrecognized Parameter Flag */
        else{
            printf
            ("\n  Parameter:%s cannot be recognized. Please see User's Manual for more details!\n",
             optstr);
//            printf("Continue? (y/N)\n");
//            char cc = getchar();
//            if(cc =='Y' || cc=='y' ){
//                /* Void */
//            }else{
//                myexit(ERRFileIO);
//            }
        }
    }
    SolverStep = MaxStep; // Improve it in future for overcasting.
    fclose (fp);

    if (forcing_mode == FORCING_NETCDF) {
        if (forcing_cfg[0] == '\0') {
            fprintf(stderr,
                    "ERROR: FORCING_MODE=NETCDF requires FORCING_CFG <path> in %s (relative paths are relative to the input directory).\n",
                    fn);
            myexit(ERRFileIO);
        }

        char base_dir[MAXLEN];
        strncpy(base_dir, fn, MAXLEN - 1);
        base_dir[MAXLEN - 1] = '\0';
        char *last_slash = strrchr(base_dir, '/');
        char *last_bslash = strrchr(base_dir, '\\');
        if (last_bslash != NULL && (last_slash == NULL || last_bslash > last_slash)) {
            last_slash = last_bslash;
        }
        if (last_slash != NULL) {
            *last_slash = '\0';
        } else {
            strcpy(base_dir, ".");
        }

        const bool is_abs = (forcing_cfg[0] == '/' || forcing_cfg[0] == '\\' ||
                             (isalpha((unsigned char)forcing_cfg[0]) && forcing_cfg[1] == ':'));
        char resolved[MAXLEN];
        if (is_abs) {
            strncpy(resolved, forcing_cfg, MAXLEN - 1);
            resolved[MAXLEN - 1] = '\0';
        } else {
            snprintf(resolved, MAXLEN, "%s/%s", base_dir, forcing_cfg);
        }

        FILE *fp_cfg = fopen(resolved, "r");
        if (fp_cfg == NULL) {
            fprintf(stderr,
                    "ERROR: FORCING_CFG is not readable: %s (from %s).\n",
                    resolved,
                    fn);
            myexit(ERRFileIO);
        }
        fclose(fp_cfg);
        strncpy(forcing_cfg, resolved, MAXLEN - 1);
        forcing_cfg[MAXLEN - 1] = '\0';

        // Derive radiation semantics from FORCING_CFG.RADIATION_KIND when user did not specify
        // RADIATION_INPUT_MODE explicitly. Error on conflict to avoid silent double-albedo.
        char product[MAXLEN] = "";
        char radiation_kind[MAXLEN] = "";
        (void)_read_kv_cfg_value(forcing_cfg, "PRODUCT", product, sizeof(product));
        const bool has_kind = _read_kv_cfg_value(forcing_cfg, "RADIATION_KIND", radiation_kind, sizeof(radiation_kind));

        const bool is_era5 = (product[0] != '\0' && strcasecmp(product, "ERA5") == 0);
        const char *default_kind = is_era5 ? "SWNET" : "SWDOWN";
        const char *effective_kind = (has_kind && radiation_kind[0] != '\0') ? radiation_kind : default_kind;

        // Store for startup logs.
        strncpy(forcing_radiation_kind_eff, effective_kind, sizeof(forcing_radiation_kind_eff) - 1);
        forcing_radiation_kind_eff[sizeof(forcing_radiation_kind_eff) - 1] = '\0';

        int desired_mode = NA_VALUE;
        if (strcasecmp(effective_kind, "SWDOWN") == 0) {
            desired_mode = SWDOWN;
        } else if (strcasecmp(effective_kind, "SWNET") == 0) {
            desired_mode = SWNET;
        } else {
            fprintf(stderr, "\n  Fatal Error: invalid RADIATION_KIND in FORCING_CFG.\n");
            fprintf(stderr, "  FORCING_CFG: %s\n", forcing_cfg);
            fprintf(stderr, "  RADIATION_KIND: %s\n", effective_kind);
            fprintf(stderr, "  Valid: SWDOWN | SWNET\n");
            myexit(ERRFileIO);
        }

        if (radiation_input_mode_user_set) {
            if (radiation_input_mode != desired_mode) {
                fprintf(stderr, "\n  Fatal Error: RADIATION_INPUT_MODE conflicts with FORCING_CFG RADIATION_KIND.\n");
                fprintf(stderr, "  cfg.para RADIATION_INPUT_MODE: %s\n", radiation_input_mode == SWNET ? "SWNET" : "SWDOWN");
                fprintf(stderr, "  FORCING_CFG RADIATION_KIND:    %s\n", effective_kind);
                fprintf(stderr, "  Fix: make them consistent. For NetCDF forcing, ET expects they match.\n");
                myexit(ERRFileIO);
            }
        } else {
            radiation_input_mode = desired_mode;
            radiation_input_mode_from_forcing_cfg = true;
        }
    }

    if (output_mode == OUTPUT_NETCDF || output_mode == OUTPUT_BOTH) {
        if (ncoutput_cfg[0] == '\0') {
            fprintf(stderr,
                    "ERROR: OUTPUT_MODE includes NETCDF and requires NCOUTPUT_CFG <path> in %s (relative paths are relative to the input directory).\n",
                    fn);
            myexit(ERRFileIO);
        }

        char base_dir[MAXLEN];
        strncpy(base_dir, fn, MAXLEN - 1);
        base_dir[MAXLEN - 1] = '\0';
        char *last_slash = strrchr(base_dir, '/');
        char *last_bslash = strrchr(base_dir, '\\');
        if (last_bslash != NULL && (last_slash == NULL || last_bslash > last_slash)) {
            last_slash = last_bslash;
        }
        if (last_slash != NULL) {
            *last_slash = '\0';
        } else {
            strcpy(base_dir, ".");
        }

        const bool is_abs = (ncoutput_cfg[0] == '/' || ncoutput_cfg[0] == '\\' ||
                             (isalpha((unsigned char)ncoutput_cfg[0]) && ncoutput_cfg[1] == ':'));
        char resolved[MAXLEN];
        if (is_abs) {
            strncpy(resolved, ncoutput_cfg, MAXLEN - 1);
            resolved[MAXLEN - 1] = '\0';
        } else {
            snprintf(resolved, MAXLEN, "%s/%s", base_dir, ncoutput_cfg);
        }

        FILE *fp_cfg = fopen(resolved, "r");
        if (fp_cfg == NULL) {
            fprintf(stderr,
                    "ERROR: NCOUTPUT_CFG is not readable: %s (from %s).\n",
                    resolved,
                    fn);
            myexit(ERRFileIO);
        }
        fclose(fp_cfg);
        strncpy(ncoutput_cfg, resolved, MAXLEN - 1);
        ncoutput_cfg[MAXLEN - 1] = '\0';
    }

    updateSimPeriod();
    if (Verbose || global_verbose_mode){
        fprintf(stdout, "* \t ETStep: %.2f min\n", ETStep);
    }
    fprintf(stdout, "* \t FORCING_MODE: %s\n", forcing_mode == FORCING_NETCDF ? "NETCDF" : "CSV");
    const char *out_mode_str =
        output_mode == OUTPUT_NETCDF ? "NETCDF" : (output_mode == OUTPUT_BOTH ? "BOTH" : "LEGACY");
    fprintf(stdout, "* \t OUTPUT_MODE: %s\n", out_mode_str);
    {
        char note[MAXLEN] = "";
        if (forcing_mode == FORCING_NETCDF && forcing_radiation_kind_eff[0] != '\0') {
            if (radiation_input_mode_from_forcing_cfg) {
                snprintf(note, sizeof(note), " (derived from FORCING_CFG RADIATION_KIND=%s)", forcing_radiation_kind_eff);
            } else if (radiation_input_mode_user_set) {
                snprintf(note, sizeof(note), " (cfg.para; FORCING_CFG RADIATION_KIND=%s)", forcing_radiation_kind_eff);
            } else {
                snprintf(note, sizeof(note), " (FORCING_CFG RADIATION_KIND=%s)", forcing_radiation_kind_eff);
            }
        }
        fprintf(stdout, "* \t RADIATION_INPUT_MODE: %s%s\n",
                radiation_input_mode == SWNET ? "SWNET" : "SWDOWN",
                note);
    }
    fprintf(stdout, "* \t SOLAR_LONLAT_MODE: %s\n", SolarLonLatModeName(solar_lonlat_mode));
    if (solar_lonlat_mode == FIXED) {
        fprintf(stdout, "* \t SOLAR_LON_DEG: %.6f\n", solar_lon_deg_fixed);
        fprintf(stdout, "* \t SOLAR_LAT_DEG: %.6f\n", solar_lat_deg_fixed);
    }
    fprintf(stdout, "* \t TERRAIN_RADIATION: %d\n", terrain_radiation);
    fprintf(stdout, "* \t RAD_FACTOR_CAP: %.6f\n", rad_factor_cap);
    fprintf(stdout, "* \t RAD_COSZ_MIN: %.6f\n", rad_cosz_min);
    if (terrain_radiation) {
        fprintf(stdout, "* \t TSR_INTEGRATION_STEP_MIN: %d min\n", tsr_integration_step_min);
    }
}
void Control_Data::write(const char *fn){
    
}
void Control_Data::getValue(const char *varname){
    
}
Print_Ctrl::Print_Ctrl(){}
void Print_Ctrl::setHeader(const char *s){
    strcpy(header, s);
}
void Print_Ctrl::open_file(int a,
                           int b,
                           int radiation_input_mode,
                           int terrain_radiation,
                           SolarLonLatMode solar_lonlat_mode,
                           double solar_lon_deg,
                           double solar_lat_deg){
    Ascii = a;
    Binary = b;
    double tmp;
    if(strlen(filename) < 1){
        fprintf(stderr, "WARNING: filename (%s)is empty.\n;", filename);
    }
    sprintf(filea, "%s.csv", filename);
    sprintf(fileb, "%s.dat", filename);

    /* Binary header is fixed-size (1024 bytes) and should be deterministic. */
    memset(header, 0, sizeof(header));
    snprintf(header, sizeof(header),
             "# SHUD output\n"
             "# Radiation input mode: %s\n"
             "# Terrain radiation (TSR): %s\n"
             "# Solar lon/lat mode: %s\n"
             "# Solar lon/lat (deg): lon=%.6f, lat=%.6f\n",
             radiation_input_mode == SWNET ? "SWNET" : "SWDOWN",
             terrain_radiation ? "ON" : "OFF",
             SolarLonLatModeName(solar_lonlat_mode),
             solar_lon_deg,
             solar_lat_deg);

    if (sink) {
        sink->onInit(filename,
                     StartTime,
                     Interval,
                     NumAll,
                     NumVar,
                     icol,
                     radiation_input_mode,
                     terrain_radiation,
                     static_cast<int>(solar_lonlat_mode),
                     solar_lon_deg,
                     solar_lat_deg);
    }

    if (Binary){
        fid_bin = fopen (fileb, "wb");
        CheckFile(fid_bin, fileb);
        fwrite(header, sizeof(char), 1024, fid_bin);
        tmp = (double) StartTime;
        fwrite( &tmp, sizeof(tmp), 1, fid_bin);
        tmp = (double) NumVar;
        fwrite( &tmp, sizeof(tmp), 1, fid_bin);
        fwrite( icol, sizeof(double), NumVar, fid_bin);
        if(global_fflush_mode){
            fflush(fid_bin);
        }
    }
    if (Ascii){
        fid_asc = fopen (filea, "w");
        CheckFile(fid_asc, filea);
        fprintf(fid_asc, "# Timestamp semantics: left endpoint (t-Interval)\n");
        fprintf(fid_asc, "%d\t %d\t %ld\n", 0, NumVar, StartTime);
        fprintf(fid_asc, "# Radiation input mode: %s\n",
                radiation_input_mode == SWNET ? "SWNET" : "SWDOWN");
        fprintf(fid_asc, "# Terrain radiation (TSR): %s\n",
                terrain_radiation ? "ON" : "OFF");
        fprintf(fid_asc, "# Solar lon/lat mode: %s\n", SolarLonLatModeName(solar_lonlat_mode));
        fprintf(fid_asc, "# Solar lon/lat (deg): lon=%.6f, lat=%.6f\n", solar_lon_deg, solar_lat_deg);
        fprintf(fid_asc, "%s", "Time_min");
        for(int i = 0; i < NumVar; i++){
            fprintf(fid_asc, " \tX%d", i + 1);
        }
        fprintf(fid_asc, "\n");
        if(global_fflush_mode){
            fflush(fid_asc);
        }
    }
}
void Print_Ctrl::Init(long st, int n, const char *s, int dt, double *x, int iFlux){
    StartTime = st;
    NumVar  = n;
    NumAll  = n;
    PrintVar = new double*[NumVar];
    buffer  = new double[NumVar];
    icol    = new double[NumVar];
    strcpy(filename, s);
    if(strlen(filename) < 1){
        fprintf(stderr, "WARNING: filename (%s)is empty.\n;", filename);
    }
    if(dt == 0 ){
        myexit(ERRCONSIS);
    }
    Interval = dt;
    for(int i=0; i<NumVar; i++){
        icol[i] = (double) (i + 1);
        PrintVar[i] = &x[i];
//        *(PrintVar[i]) = 0.0;
        buffer[i] = 0.0;
    }
    if(iFlux){
        tau = 1440.;
    }else{
        tau = 1.;
    }
}
void Print_Ctrl::InitIJ(long st, int n, const char *s, int dt, double **x, int j, int iFlux){
    StartTime = st;
    NumVar  = n;
    NumAll  = n;
    PrintVar = new double*[NumVar];
    buffer  = new double[NumVar];
    icol    = new double[NumVar];
    strcpy(filename, s);
    if(dt == 0 ){
        myexit(ERRCONSIS);
    }
    Interval = dt;
    for(int i=0; i<NumVar; i++){
        icol[i] = (double) (i + 1);
        PrintVar[i] = &(x[i][j]);
//        *(PrintVar[i]) = 0.0;
        buffer[i] = 0.0;
    }
    if(iFlux){
        tau = 1440.;
    }else{
        tau = 1;
    }
}

void Print_Ctrl::Init(long st, int n, const char *s, int dt, double *x, int iFlux, int *flag_IO){
    StartTime = st;
    NumAll = n;
    strcpy(filename, s);
    if(strlen(filename) < 1){
        fprintf(stderr, "WARNING: filename (%s)is empty.\n;", filename);
    }
    if(dt == 0 ){
        myexit(ERRCONSIS);
    }
    
    Interval = dt;
    NumVar = 0;
    for(int i = 0; i < n; i++){
        if(flag_IO[i]){ /* IO is TRUE*/
            NumVar++;
        }
    }
    if(NumVar <= 0){
        fprintf(stderr, "WARNING: Empty columns in %s.\n;", filename);
    }
    buffer = new double[NumVar];
    PrintVar = new double*[NumVar];
    icol    = new double[NumVar];
    int k = 0;
    for(int i = 0; i < n; i++){
        if(flag_IO[i]){ /* IO is TRUE*/
            PrintVar[k] = &x[i];
            icol[k] = (double) (i + 1);
            buffer[k] = 0.0;
            k++;
        }
    }
    if(iFlux){
        tau = 1440.;
    }else{
        tau = 1.;
    }
}

void Print_Ctrl::InitIJ(long st, int n, const char *s, int dt, double **x, int j, int iFlux, int *flag_IO){
    StartTime = st;
    NumAll = n;
    strcpy(filename, s);
    if(dt == 0 ){
        myexit(ERRCONSIS);
    }
    Interval = dt;
    NumVar = 0;
    for(int i = 0; i < n; i++){
        if(flag_IO[i]){ /* IO is TRUE*/
            NumVar++;
        }
    }
    if(NumVar <= 0){
        fprintf(stderr, "WARNING: Empty columns in %s.\n;", filename);
    }
    buffer = new double[NumVar];
    PrintVar = new double*[NumVar];
    icol    = new double[NumVar];
    int k = 0;
    for(int i = 0; i < n; i++){
        if(flag_IO[i]){ /* IO is TRUE*/
            PrintVar[k] = &(x[i][j]);
            icol[k] = (double) (i + 1);
            buffer[k] = 0.0;
            k++;
        }
    }
    
    if(iFlux){
        tau = 1440.;
    }else{
        tau = 1;
    }
}
Print_Ctrl::~Print_Ctrl(){
    close_file();
    if (PrintVar != NULL) delete[] PrintVar;
    if (buffer != NULL) delete[] buffer;
    if (icol != NULL) delete[] icol;
}
void Print_Ctrl::fun_printBINARY(double t, double dt){
    fwrite (&t, sizeof (double), 1, fid_bin);
    fwrite (buffer, sizeof (double), NumVar, fid_bin);
    if(global_fflush_mode){
        fflush (fid_bin); // DEBUG
    }
}
void Print_Ctrl::fun_printASCII(double t, double dt){
    fprintf(fid_asc, "%.1f\t", t);
    for (int i = 0; i < NumVar; i++){
        fprintf (fid_asc, "%e\t", buffer[i]);
    }
    fprintf (fid_asc, "\n");
    if(global_fflush_mode){
//    fflush (fid_asc); // DEBUG
    }
}
void Print_Ctrl::close_file(){
    if (Binary){
        if(fid_bin != NULL) {
            fclose(fid_bin);
            fid_bin = NULL;
        }
    }
    if (Ascii){
        if(fid_asc != NULL) {
            fclose(fid_asc);
            fid_asc = NULL;
        }
    }

    if (sink && !sink_closed) {
        sink->onClose();
        sink_closed = true;
    }
    
}
void Print_Ctrl::PrintData(double dt, double t){
    long long t_floor;
    NumUpdate++; /* Number of times to push data into the buffer*/
    for (int i = 0; i < NumVar; i++){
        buffer[i] += *(PrintVar[i]);
    }
    /* Small epsilon (minutes) for numerical stability when detecting output boundaries.
       Prevents missing output points due to floating-point errors (e.g., 59.9999999 -> 60)
       while avoiding early triggers. Value 0.001 min (0.06 s) is far below output
       resolution (Time_min printed to 0.1 min) and typical solver steps. */
    constexpr double OUTPUT_TRIGGER_EPSILON = 0.001;
    t_floor = static_cast<long long>(floor(t + OUTPUT_TRIGGER_EPSILON));
    if ((t_floor % Interval) == 0){
        for (int i = 0; i < NumVar; i++){
            buffer[i] *= tau / NumUpdate ; /* Get the mean in the time-interval*/
        }
        NumUpdate = 0;
        const double t_quantized =
            static_cast<double>(t_floor - static_cast<long long>(Interval));
        if(Ascii){
            fun_printASCII(t_quantized, dt);
        }
        if(Binary){
            fun_printBINARY(t_quantized, dt);
        }
        if (sink) {
            sink->onWrite(t_quantized, NumVar, buffer);
        }
        for (int i = 0; i < NumVar; i++){
            buffer[i] = 0.;  /* Reset the buffer */
        }
    }
}
