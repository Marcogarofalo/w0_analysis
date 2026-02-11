#define functions_w0_C
#include "functions_w0.hpp"
#include "correlators_analysis.hpp"
#include "read.hpp"

void make_ratio_of_jacks(double**** final, int Njack, int i, int T, double**** num, int i1, double**** den, int i2) {

    // double sum_r = 0;
    for (int j = 0; j < Njack; j++) {
        for (int tf = 0; tf < T; tf++) {
            for (int j1 = 0; j1 < Njack - 1; j1++) {
                if (j1 == j)
                    continue;
                double r = 0;
                for (int j2 = 0; j2 < Njack - 1; j2++) {
                    if (j2 == j)
                        continue;
                    else {
                        r += exp(den[j2][i2][0][0] - num[j1][i1][tf][0]);
                        // printf("r=%g   %g    %g \n", r, conf_jack_r[j2][0][0][0], conf_jack_rO[j1][0][tf][0]);
                    }
                }
                final[j][i][tf][0] += 1.0 / r;
            }
        }
    }
}

double**** bin_intoN_exp(double**** data, int ivar, int T, int Nconf_in, int Nb) {

    double clustSize = ((double)Nconf_in) / ((double)Nb);
    // double clustSize = ((double)confs.confs_after_binning) / ((double)Nb);
    double**** to_write = calloc_corr(Nb, 1, T);
    for (size_t iClust = 0; iClust < Nb; iClust++) {

        /// Initial time of the bin
        const double binBegin = iClust * clustSize;
        /// Final time of the bin
        const double binEnd = binBegin + clustSize;
        double binPos = binBegin;
        const size_t iConf0 = floor(binPos + 1e-10);
        do {
            /// Index of the configuration related to the time
            const size_t iConf = floor(binPos + 1e-10);

            /// Rectangle left point
            const double beg = binPos;

            /// Rectangle right point
            const double end = std::min(binEnd, iConf + 1.0);

            /// Rectangle horizontal size
            const double weight = end - beg;

            // Perform the operation passing the info
            for (int t = 0; t < T; t++) {
                // printf("iConf %ld  iConf0 %ld ivar %d  t %d\n", iConf, iConf0, ivar, t);
                to_write[iClust][ivar][t][0] += weight * exp(data[iConf][ivar][t][0] - data[iConf0][ivar][t][0]);
                // to_write[iClust][ivar][t][1] += weight * exp(data[iConf][ivar][t][1]- data[0][ivar][t][1]); // imag part not implemented
            }
            // printf("Cluster=%ld  iConf=%ld  weight=%g  size=%g  end=%g  beg=%g\n",iClust,iConf,weight,clustSize, end,beg);
            // Updates the position
            binPos = end;
        } while (binEnd - binPos > 1e-10);
        for (int t = 0; t < T; t++) {
            to_write[iClust][ivar][t][0] /= ((double)clustSize);
            // to_write[iClust][ivar][t][1] /= ((double)clustSize);

            to_write[iClust][ivar][t][0] = log(to_write[iClust][ivar][t][0]);
            to_write[iClust][ivar][t][0] += data[iConf0][ivar][t][0];
        }
    }
    return to_write;
}

double**** bin_intoN_exp1(double**** data, double**** data_noexp, int ivar, int ivar_noexp, int T, int Nconf_in, int Nb) {

    double clustSize = ((double)Nconf_in) / ((double)Nb);
    // double clustSize = ((double)confs.confs_after_binning) / ((double)Nb);
    double**** to_write = calloc_corr(Nb, 1, T);
    for (size_t iClust = 0; iClust < Nb; iClust++) {

        /// Initial time of the bin
        const double binBegin = iClust * clustSize;
        /// Final time of the bin
        const double binEnd = binBegin + clustSize;
        double binPos = binBegin;
        const size_t iConf0 = floor(binPos + 1e-10);
        do {
            /// Index of the configuration related to the time
            const size_t iConf = floor(binPos + 1e-10);

            /// Rectangle left point
            const double beg = binPos;

            /// Rectangle right point
            const double end = std::min(binEnd, iConf + 1.0);

            /// Rectangle horizontal size
            const double weight = end - beg;

            // Perform the operation passing the info
            for (int t = 0; t < T; t++) {
                to_write[iClust][ivar][t][0] += weight * data_noexp[iConf][ivar_noexp][t][0] * exp(data[iConf][ivar][0][0] - data[iConf0][ivar][0][0]);
                // to_write[iClust][ivar][t][1] += weight * exp(data[iConf][ivar][t][1]- data[0][ivar][t][1]); // imag part not implemented
            }
            // printf("Cluster=%ld  iConf=%ld  weight=%g  size=%g  end=%g  beg=%g\n",iClust,iConf,weight,clustSize, end,beg);
            // Updates the position
            binPos = end;
        } while (binEnd - binPos > 1e-10);
        for (int t = 0; t < T; t++) {
            to_write[iClust][ivar][t][0] /= ((double)clustSize);
            // to_write[iClust][ivar][t][1] /= ((double)clustSize);

            to_write[iClust][ivar][t][0] = log(to_write[iClust][ivar][t][0]);
            to_write[iClust][ivar][t][0] += data[iConf0][ivar][0][0];
        }
    }
    return to_write;
}

double lhs_function_w0_eg(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    double r = in[j][id][t][0];
    // if(j==fit_info.Njack-1) printf("%d   %g\n",t,r);
    return r;
}

double lhs_function_W_rew(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];  // <W(t)r>
    int id1 = fit_info.corr_id[1]; // <r>
    double r = in[j][id][t][0] / in[j][id1][0][0];
    return r;
}

double lhs_function_Wt_p_dmcorr(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    int id_cor = fit_info.corr_id[1];
    int reim = fit_info.myen[0];
    double dmu = fit_info.ext_P[0][j];

    double Wt = in[j][id][t][0];
    double loop = in[j][id_cor][0][reim];
    double loop_Wt = in[j][id_cor][t][reim];

    double r = Wt - dmu * (loop_Wt - loop * Wt);
    // if(j==fit_info.Njack-1) printf("%d   %g   Wt=%g   dm=%g  loop_Wt=%g   loop=%g    der=%g  corr=%g \n",t,r,Wt, dmu, loop_Wt, loop,
    // (loop_Wt - loop * Wt), (loop_Wt - loop * Wt)*dmu );
    return r;
}

double lhs_function_M_PS_p_dmcorr(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    int id_cor = fit_info.corr_id[1];
    int reim = fit_info.myen[0];
    double dmu = fit_info.ext_P[0][j];

    double Ct = in[j][id][t][0];
    double loop = in[j][id_cor][0][reim];
    double loop_Ct = in[j][id_cor][t][reim];

    double rt = Ct + dmu * (loop_Ct - loop * Ct);

    Ct = in[j][id][t + 1][0];
    loop = in[j][id_cor][0][reim];
    loop_Ct = in[j][id_cor][t + 1][reim];

    double rtp1 = Ct + dmu * (loop_Ct - loop * Ct);

    double mass = M_eff_T_ct_ctp1(t, fit_info.T, rt, rtp1);

    return mass;
}

double lhs_function_f_PS_p_dmcorr(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    int id_cor = fit_info.corr_id[1];
    int reim = fit_info.myen[0];
    double dmu = fit_info.ext_P[3][j];

    double Ct = in[j][id][t][0];
    double loop = in[j][id_cor][0][reim];
    double loop_Ct = in[j][id_cor][t][reim];

    double rt = Ct + dmu * (loop_Ct - loop * Ct);

    double M = fit_info.ext_P[0][j];
    double mu1 = fit_info.ext_P[1][j];
    double mu2 = fit_info.ext_P[2][j];

    double me = sqrt(rt * 2 * M / (exp(-t * M) + exp(-(fit_info.T - t) * M)));

    if (t == fit_info.T / 2 - 1 && j == fit_info.Njack - 1)
        printf("t=%d  fpi=%g   G=%g   A=%g\n", t, (mu1 + mu2) * me / (M * sinh(M)), me, (rt / (exp(-t * M) + exp(-(fit_info.T - t) * M))));
    return (mu1 + mu2) * me / (M * sinh(M));
}

double lhs_function_mpcac_p_dmcorr(int j, double**** in, int t, struct fit_type fit_info) {
    int idV = fit_info.corr_id[0];
    int idP = fit_info.corr_id[1];
    int id_corV = fit_info.corr_id[2];
    int id_corP = fit_info.corr_id[3];
    int reim = fit_info.myen[0];
    double dmu = fit_info.ext_P[0][j];

    double Ct = in[j][idV][t][0];
    double loop = in[j][id_corV][0][reim];
    double loop_Ct = in[j][id_corV][t][reim];

    double Vt = Ct + dmu * (loop_Ct - loop * Ct);

    Ct = in[j][idV][t + 1][0];
    loop = in[j][id_corV][0][reim];
    loop_Ct = in[j][id_corV][t + 1][reim];

    double Vtp1 = Ct + dmu * (loop_Ct - loop * Ct);

    Ct = in[j][idP][t + 1][0];
    loop = in[j][id_corP][0][reim];
    loop_Ct = in[j][id_corP][t + 1][reim];

    double PP = Ct + dmu * (loop_Ct - loop * Ct);

    double mpcac = -(Vtp1 - Vt) / (2. * PP);
    ;

    return mpcac;
}

double lhs_function_Wt_der_mu(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    int id_cor = fit_info.corr_id[1];
    int reim = fit_info.myen[0];
    // double dmu = fit_info.ext_P[0][j];

    double Wt = in[j][id][t][0];
    double loop = in[j][id_cor][0][reim];
    double loop_Wt = in[j][id_cor][t][reim];

    double r = -(loop_Wt - loop * Wt);
    // if(j==fit_info.Njack-1) printf("%d   %g\n",t,r);
    return r;
}

double lhs_function_Wt_p_dm_all_corr(int j, double**** in, int t, struct fit_type fit_info) {
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

    double r = Wt - dmul * (loop_Wt_l - loop_l * Wt) -
        dmus * (loop_Wt_s - loop_s * Wt) -
        dmuc * (loop_Wt_c - loop_c * Wt);
    // if(j==fit_info.Njack-1) printf("%d   %g\n",t,r);
    return r;
}
double lhs_zero(int j, double**** in, int t, struct fit_type fit_info) {
    return 0;
}

double lhs_mpcac(int j, double**** in, int t, struct fit_type fit_info) {
    int id_V = fit_info.corr_id[0];
    int id_P = fit_info.corr_id[1];
    double mu = fit_info.ave_P[0];
    // we should take -Im of V0P5 which is saved as real part in this data
    // double r = -(in[j][id_V][t + 1][0] - in[j][id_V][t][0]) / (2. * in[j][id_P][t][0]);
    //double r = -(in[j][id_V][t + 1][0] - in[j][id_V][t-1][0]) / (4. * in[j][id_P][t][0]);
    // printf("%d  %d   %d   %d\n", j, id_V, id_P, t);
    double r = -(in[j][id_V][t + 1][0] - in[j][id_V][(t - 1) % fit_info.T][0]) / (4. * mu * in[j][id_P][t][0]);
    return r;
}

double lhs_dM_sea(int j, double**** in, int t, struct fit_type fit_info) {
    // fit_info.myen = { 1, 1 }; // sign , reim   for mass correction of the pion

    int id = fit_info.corr_id[0];
    int id_cor = fit_info.corr_id[1];
    int reim = fit_info.myen[0];

    double Ct = in[j][id][t][0];
    double loop = in[j][id_cor][0][reim];
    double loop_Ct = in[j][id_cor][t][reim];

    double rt = (loop_Ct - loop * Ct) / Ct;

    Ct = in[j][id][t + 1][0];
    loop = in[j][id_cor][0][reim];
    loop_Ct = in[j][id_cor][t + 1][reim];

    double rtp1 = (loop_Ct - loop * Ct) / Ct;

    int T = fit_info.T;
    double M = fit_info.ext_P[0][j];

    double den = (T / 2 - (t + 1)) * tanh(M * (T / 2 - (t + 1)));
    den -= (T / 2 - (t)) * tanh(M * (T / 2 - (t)));
    double r = (rtp1)-(rt);
    r /= den;
    return r;
}
double lhs_dfpi_sea(int j, double**** in, int t, struct fit_type fit_info) {
    // fit_info.myen = { 1, 1 }; // sign , reim   for mass correction of the pion

    int id = fit_info.corr_id[0];
    int id_cor = fit_info.corr_id[1];
    int reim = fit_info.myen[0];
    int T = fit_info.T;
    double M = fit_info.ext_P[0][j];
    double fpi = fit_info.ext_P[1][j];
    double dM = fit_info.ext_P[2][j];
    double mu = fit_info.ext_P[3][j];

    double G = fpi * (M * sinh(M)) / (2 * mu);
    double A = G * G / (2 * M);

    double Ct = in[j][id][t][0];
    double loop = in[j][id_cor][0][reim];
    double loop_Ct = in[j][id_cor][t][reim];

    double dA = (loop_Ct - loop * Ct) / Ct;
    dA += dM * T / 2.0;
    dA -= dM * (T / 2.0 - t) * tanh(M * (T / 2.0 - t));
    dA *= A;

    double dG = M * dA / G + G * dM / (2 * M);

    // double df = (2 * mu) * dG / (M * M) + 2 * G / (M * M) - (4 * mu) * G * dM / (M * M * M);
    double df = (2 * mu) * dG / (M * sinh(M));
    // df += 2 * G / (M * sinh(M)); // do not know why this is not here, maybe it is a valence term
    df -= (2 * mu) * G * dM * (sinh(M) + M * cosh(M)) / ((M * sinh(M)) * (M * sinh(M)));

    if (t == T / 2 - 1 && j == fit_info.Njack - 1)
        printf("dM=%g  df=%g  fpi=%g  G=%g  A=%g\n", dM, df, fpi, G, A);
    return df;
}

double lhs_plateau_dM_dmu(int j, double**** in, int t, struct fit_type fit_info) {
    int idp = fit_info.corr_id[0];
    int id = fit_info.corr_id[1];

    double Mp = M_eff_T(t, fit_info.T, in[j][idp]);
    double M = M_eff_T(t, fit_info.T, in[j][id]);

    double dmu = fit_info.ave_P[0];

    // we should take -Im of V0P5 which is saved as real part in this data
    double r = (Mp - M) / dmu;

    return r;
}

double lhs_plateau_ratio_dM_dmu(int j, double**** in, int t, struct fit_type fit_info) {
    int idp = fit_info.corr_id[0];
    int idl = fit_info.corr_id[1];
    int id = fit_info.corr_id[2];

    double dmu2 = fit_info.ave_P[0];
    double dmul = fit_info.ave_P[1];

    double M2 = M_eff_T(t, fit_info.T, in[j][idp]);
    double Ml = M_eff_T(t, fit_info.T, in[j][idl]);
    double M = M_eff_T(t, fit_info.T, in[j][id]);

    // we should take -Im of V0P5 which is saved as real part in this data
    double r2 = (M2 - M) / dmu2;
    double rl = (Ml - M) / dmul;

    return r2 / rl;
}


double lhs_plateau_df_dmu(int j, double**** in, int t, struct fit_type fit_info) {
    std::vector<int> id(fit_info.corr_id);
    double* Mp = fit_info.ext_P[0];
    double* M = fit_info.ext_P[1];
    double* mu1 = fit_info.ext_P[2];
    double* mu2 = fit_info.ext_P[3];

    fit_info.ext_P[0] = Mp;
    fit_info.ext_P[1] = mu1;
    fit_info.ext_P[2] = mu2;
    double fp = lhs_function_f_PS(j, in, t, fit_info);

    fit_info.corr_id[0] = id[1];
    fit_info.ext_P[0] = M;
    fit_info.ext_P[1] = mu1;
    fit_info.ext_P[2] = mu2;
    double f = lhs_function_f_PS(j, in, t, fit_info);

    // rstore the original values
    fit_info.corr_id[0] = id[0];
    fit_info.ext_P[0] = Mp;
    fit_info.ext_P[1] = M;
    fit_info.ext_P[2] = mu1;
    fit_info.ext_P[3] = mu2;

    double dmu = fit_info.ave_P[0];

    return (fp - f) / dmu;
}
