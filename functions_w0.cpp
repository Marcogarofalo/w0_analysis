#define functions_w0_C
#include "functions_w0.hpp"

double lhs_function_w0_eg(int j, double ****in, int t, struct fit_type fit_info)
{
    int id = fit_info.corr_id[0];
    double r = in[j][id][t][0];
    // if(j==fit_info.Njack-1) printf("%d   %g\n",t,r);
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
    // if(j==fit_info.Njack-1) printf("%d   %g\n",t,r);
    return r;
}

double lhs_function_Wt_p_dm_all_corr(int j, double ****in, int t, struct fit_type fit_info)
{
    int id = fit_info.corr_id[0];
    int id_cor_l = fit_info.corr_id[1];
    int id_cor_s = fit_info.corr_id[2];
    int id_cor_c = fit_info.corr_id[3];
    int reim = fit_info.myen[0];
    double dmul = fit_info.ext_P[0][j];
    double dmus = fit_info.ext_P[1][j];
    double dmuc = fit_info.ext_P[2][j];

    double Wt = in[j][id][t][0];

    double loop_l = in[j][id_cor_l][0][reim];
    double loop_Wt_l = in[j][id_cor_l][t][reim];

    double loop_s = in[j][id_cor_s][0][reim];
    double loop_Wt_s = in[j][id_cor_s][t][reim];

    double loop_c = in[j][id_cor_l][0][reim];
    double loop_Wt_c = in[j][id_cor_l][t][reim];

    double r = Wt + dmul * (loop_Wt_l - loop_l * Wt) +
               dmus * (loop_Wt_s - loop_s * Wt) +
               dmuc * (loop_Wt_c - loop_c * Wt);
    // if(j==fit_info.Njack-1) printf("%d   %g\n",t,r);
    return r;
}
