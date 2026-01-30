#ifndef WaterBalanceDiag_hpp
#define WaterBalanceDiag_hpp

#include <cstdio>
#include <string>
#include <vector>

class Model_Data;

class WaterBalanceDiag {
public:
    explicit WaterBalanceDiag(Model_Data *md);
    ~WaterBalanceDiag();

    WaterBalanceDiag(const WaterBalanceDiag &) = delete;
    WaterBalanceDiag &operator=(const WaterBalanceDiag &) = delete;

    void enable();
    int isEnabled() const;

    void onETUpdate();
    void sample(double t_min, const double *DY);
    void maybeWrite(double t_min);

private:
    Model_Data *md = nullptr;
    int enabled = 0;

    int interval_min = 0;
    double last_t = 0.0;
    int has_prev = 0;
    long long last_written_floor = -1;

    std::vector<double> ic_raw; /* qEleE_IC right after ET() */

    std::vector<double> flux3_acc;
    std::vector<double> fluxfull_acc;

    std::vector<double> budget3_acc;
    std::vector<double> budgetfull_acc;

    std::vector<double> s3_prev;
    std::vector<double> sfull_prev;

    std::vector<double> resid3;
    std::vector<double> residfull;
    std::vector<double> resid3_budget;
    std::vector<double> residfull_budget;
    std::vector<double> icol;

    FILE *fp3 = nullptr;
    FILE *fpfull = nullptr;
    FILE *fp3_budget = nullptr;
    FILE *fpfull_budget = nullptr;
    FILE *fpbasin = nullptr;

    std::vector<int> outlet_riv_idx;

    std::vector<double> basin_values;
    std::vector<double> basin_icol;

    double basin_s_ele_prev_m3 = 0.0;
    double basin_s_riv_prev_m3 = 0.0;

    double basin_p_acc_m3 = 0.0;
    double basin_et_acc_m3 = 0.0;
    double basin_qout_acc_m3 = 0.0;
    double basin_qedge_acc_m3 = 0.0;
    double basin_qbc_acc_m3 = 0.0;
    double basin_qss_acc_m3 = 0.0;
    double basin_noncons_edge_acc_m3 = 0.0;

    void openFiles();
    void closeFiles();

    void writeDatHeader(FILE *fp,
                        const char *label,
                        int num_var,
                        const std::vector<double> &col_ids) const;
    void writeDatRecord(FILE *fp, double t_min, int num_var, const std::vector<double> &values) const;

    std::vector<double> makeIcol(int n) const;
    std::string outputPrefix() const;

    void initBasinOutlets();
    double basinElementStorageFull_m3() const;
    double basinRiverStorage_m3();
    double basinOutletDischarge_m3min();
    double basinBoundaryEdgeOutflow_m3min() const;
};

#endif /* WaterBalanceDiag_hpp */
