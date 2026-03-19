#ifndef functions_w0_H
#define functions_w0_H
#include <array>
#include "non_linear_fit.hpp"
static constexpr std::array<double, 4> coeffs = { -1.0, -2.0,-5.0,-10.0 };
static constexpr std::array<double, 4> coeffs_mc = { -100, -200, -400 , -800};

static constexpr double hbarc = 197.326963; // MeV*fm
static constexpr double fpi_MeV = 130.5;
static constexpr double fpi_MeV_err = 0.04;

static constexpr double Mpi_MeV = 135;
static constexpr double Mpi_MeV_err = 0.2;

static constexpr double w0_fm = 0.17236;
static constexpr double w0_fm_err = 0.0000002; // 0.00070

static constexpr double w0_MeV = w0_fm / hbarc;
static constexpr double w0_MeV_err = w0_fm_err / hbarc; // 0.00070

static constexpr double MK_MeV = 494.6;
static constexpr double MK_MeV_err = 0.3;

static constexpr double MDs_MeV = 1967.0;
static constexpr double MDs_MeV_err = 0.4;

double ****bin_intoN_exp(double ****data, int ivar, int T, int Nconf_in, int Nb);
double ****bin_intoN_exp1(double ****data, double ****data_noexp, int ivar, int ivar_noexp, int T, int Nconf_in, int Nb);
void make_ratio_of_jacks(double ****final, int Njack, int i, int T, double ****num, int i1, double ****den, int i2);


double lhs_function_w0_eg(int j, double**** in, int t, struct fit_type fit_info);
double lhs_function_W_rew(int j, double ****in, int t, struct fit_type fit_info);
double lhs_function_Wt_p_dmcorr(int j, double ****in, int t, struct fit_type fit_info);
double lhs_function_Wt_p_dm_all_corr(int j, double ****in, int t, struct fit_type fit_info);
double lhs_function_Wt_der_mu(int j, double ****in, int t, struct fit_type fit_info);
double lhs_dM_sea(int j, double**** in, int t, struct fit_type fit_info);
double lhs_dfpi_sea(int j, double**** in, int t, struct fit_type fit_info);

double lhs_zero(int j, double ****in, int t, struct fit_type fit_info);
double lhs_mpcac(int j, double ****in, int t, struct fit_type fit_info);

double lhs_function_M_PS_p_dmcorr(int j, double ****in, int t, struct fit_type fit_info);
double lhs_function_f_PS_p_dmcorr(int j, double ****in, int t, struct fit_type fit_info);

double lhs_plateau_dM_dmu(int j, double ****in, int t, struct fit_type fit_info);
double lhs_plateau_df_dmu(int j, double ****in, int t, struct fit_type fit_info);

double lhs_plateau_ratio_dM_dmu(int j, double ****in, int t, struct fit_type fit_info);

double lhs_me(int j, double**** in, int t, struct fit_type fit_info);
double lhs_fpi_P5A0(int j, double**** in, int t, struct fit_type fit_info);
double lhs_plateau_fpi_P5A0(int j, double**** in, int t, struct fit_type fit_info);
#endif
