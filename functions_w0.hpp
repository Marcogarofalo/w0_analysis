#ifndef functions_w0_H
#define functions_w0_H
#include "non_linear_fit.hpp"

double ****bin_intoN_exp(double ****data, int ivar, int T, int Nconf_in, int Nb);
double ****bin_intoN_exp1(double ****data, double ****data_noexp, int ivar, int ivar_noexp, int T, int Nconf_in, int Nb);



double lhs_function_w0_eg(int j, double**** in, int t, struct fit_type fit_info);
double lhs_function_Wt_p_dmcorr(int j, double ****in, int t, struct fit_type fit_info);
double lhs_function_Wt_p_dm_all_corr(int j, double ****in, int t, struct fit_type fit_info);
double lhs_function_Wt_der_mu(int j, double ****in, int t, struct fit_type fit_info);
#endif
