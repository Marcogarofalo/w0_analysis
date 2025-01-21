#define functions_w0_C
#include "functions_w0.hpp"

double lhs_function_w0_eg(int j, double ****in, int t, struct fit_type fit_info)
{
    int id = fit_info.corr_id[0];
    double r = in[j][id][t][0];
    return r;
}

double lhs_function_Wt_p_dmcorr(int j, double ****in, int t, struct fit_type fit_info)
{
    int id = fit_info.corr_id[0];
    int id_cor = fit_info.corr_id[1];
    int reim = fit_info.myen[0];
    double dmu = fit_info.ext_P[0][j];

    double Wt = in[j][id][t][0];
    double loop = in[j][id_cor][0][reim];
    double loop_Wt = in[j][id_cor][t][reim];

    double r = Wt + dmu * (loop_Wt - loop * Wt);
    return r;
}
