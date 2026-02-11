#define CONTROL

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "global.hpp"
#include "read.hpp"
#include "resampling.hpp"
// #include "m_eff.hpp"
// #include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "mutils.hpp"
#include "various_fits.hpp"
// #include "correlators_analysis.hpp"
// #include "eigensystem.hpp"
#include "fit_all.hpp"
#include "global.hpp"
#include "non_linear_fit.hpp"
#include "resampling_new.hpp"
#include "tower.hpp"

#include "functions_w0.hpp"

#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

constexpr double Mpi_MeV = 135;
constexpr double Mpi_MeV_err = 0.2;
// #include "do_analysis_charm.hpp"

generic_header read_header(FILE* stream) {
    generic_header header;
    int ir = 0;
    ir += fread(&header.T, sizeof(int), 1, stream);
    ir += fread(&header.L, sizeof(int), 1, stream);
    int s;
    ir += fread(&s, sizeof(int), 1, stream);
    header.mus = std::vector<double>(s);
    for (int i = 0; i < s; i++) {
        ir += fread(&header.mus[i], sizeof(double), 1, stream);
    }
    ir += fread(&s, sizeof(int), 1, stream);
    header.thetas = std::vector<double>(s);
    for (int i = 0; i < s; i++) {
        ir += fread(&header.thetas[i], sizeof(double), 1, stream);
    }

    ir += fread(&header.Njack, sizeof(int), 1, stream);
    header.struct_size = ftell(stream);
    return header;
}

double read_single_Nobs(FILE* stream, int header_size, int Njack) {
    int Nobs;
    long int tmp;
    int s = header_size;

    // size_t i = fread(&Njack, sizeof(int), 1, stream);

    fseek(stream, 0, SEEK_END);
    tmp = ftell(stream);
    tmp -= header_size;

    s = Njack;

    Nobs = (tmp) / ((s) * sizeof(double));

    fseek(stream, header_size, SEEK_SET);

    return Nobs;
}

data_single read_single_dataj(FILE* stream) {

    int Njack;
    int Nobs;

    // read_single_Njack_Nobs(stream, header.header_size, Njack, Nobs);
    //  data_single dj(Nobs,Njack);
    data_single dj;
    dj.header.read_header_jack(stream);

    dj.Nobs = read_single_Nobs(stream, dj.header.struct_size, dj.header.Njack);
    dj.Njack = dj.header.Njack;
    printf("read_single_dataj: Njack=%d Nobs=%d\n", dj.Njack, dj.Nobs);
    dj.jack = double_malloc_2(dj.Nobs, dj.Njack);

    //
    size_t i = 0;
    for (int obs = 0; obs < dj.Nobs; obs++) {
        i += fread(dj.jack[obs], sizeof(double), dj.Njack, stream);
    }
    return dj;
}

data_all read_all_the_files(std::vector<std::string> files, const char* resampling) {
    data_all jackall;
    jackall.resampling = resampling;
    // jackall->en = (data_single*)malloc(sizeof(data_single) * files.size());
    jackall.en = new data_single[files.size()];
    jackall.ens = files.size();
    int count = 0;
    for (std::string s : files) {
        printf("reading %s\n", s.c_str());
        FILE* f = open_file(s.c_str(), "r");

        // read_single_dataj(f, params, &(jackall->en[count]));
        jackall.en[count] = read_single_dataj(f);
        jackall.en[count].resampling = resampling;
        count++;
        fclose(f);
    }
    return jackall;
}

double lhs_fun(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r = gjack.en[e].jack[fit_info.corr_id[0]][j]; // d(w0/a)/d(a*mu)
    return r;
}
double lhs_fun_ratio(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r = gjack.en[e].jack[fit_info.corr_id[0]][j];
    r /= gjack.en[e].jack[fit_info.corr_id[1]][j];
    return r;
}
double lhs_fun_diff(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r = gjack.en[e].jack[fit_info.corr_id[0]][j];
    r -= gjack.en[e].jack[fit_info.corr_id[1]][j];
    return r;
}

double lhs_fun_fpi(int n, int e, int j, data_all gjack, struct fit_type fit_info) {

    double r = gjack.en[e].jack[fit_info.corr_id[0]][j]; // fpi
    return r;
}
double rhs_1overamu(int n, int Nvar, double* x, int Npar, double* P) {
    double amu = x[0];
    double r = P[0] / amu;
    return r;
}



double rhs_1overamu_mu2_mu3(int n, int Nvar, double* x, int Npar, double* P) {
    double amu = x[0];
    double r = P[0] / amu + P[1] / (amu * amu) + P[2] / (amu * amu * amu);
    return r;
}
double rhs_1_mu(int n, int Nvar, double* x, int Npar, double* P) {
    double amu = x[0];
    double r = P[0] + P[1] * amu;
    return r;
}

double rhs_1overamu_1overamu2(int n, int Nvar, double* x, int Npar, double* P) {
    double amu = x[0];
    double r = P[0] / amu + P[1] / (amu * amu);
    return r;
}
double rhs_1overamu2(int n, int Nvar, double* x, int Npar, double* P) {
    double amu = x[0];
    double r = P[0] / (amu * amu);
    return r;
}

double rhs_1overamu_1overamu2_1overamu3(int n, int Nvar, double* x, int Npar, double* P) {
    double amu = x[0];
    double r = P[0] / amu + P[1] / (amu * amu) + P[2] / (amu * amu * amu);
    return r;
}


int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir \n argc=%d\n", argc);
    char namefile[NAMESIZE];

    std::vector<std::string> files;
    // std::vector<std::string> beta_names;

    files.emplace_back("../../data/jackknife/jack_scale_setting_system_B64");
    files.emplace_back("../../data/jackknife/jack_scale_setting_system_C80");
    files.emplace_back("../../data/jackknife/jack_scale_setting_system_D96");
    files.emplace_back("../../data/jackknife/jack_scale_setting_system_E112");

    std::vector<std::vector<int>> myen(1, std::vector<int>(files.size()));
    for (int i = 0; i < myen[0].size(); i++) {
        myen[0][i] = i;
    }

    data_all jackall = read_all_the_files(files, argv[1]);
    printf("we read all\n");
    jackall.create_generalised_resampling();

    data_all jackall_chi = read_all_the_files(files, argv[1]);
    jackall_chi.create_generalised_resampling();

    int Njack = jackall.en[0].Njack;
    if (strcmp(argv[1], "jack") == 0) {
        myres = new resampling_jack(Njack - 1);
    }
    else if (strcmp(argv[1], "boot") == 0) {
        myres = new resampling_boot(Njack - 1);
    }


    //////////////////////////////////////////////////////////////
    // make table results
    //////////////////////////////////////////////////////////////
    mysprintf(namefile, NAMESIZE, "%s/data_from_fpi.txt", argv[3]);
    FILE* summary_out = open_file(namefile, "w+");
    fprintf(summary_out, "ens   a[fm] da[fm]  w0/a dw0/a  amul damul amus damus  amuc  damuc\n");
    for (int i = 0;i < myen[0].size();i++) {
        size_t lastUnderscore = files[i].find_last_of('_');

        error(lastUnderscore == std::string::npos, 1, "error in file name", "cannot read ens from name %s", files[i].c_str());
        // 2. Estrae tutto ciò che segue l'ultimo underscore
        std::string result = files[i].substr(lastUnderscore + 1);
        fprintf(summary_out, "%s  ", result.c_str());
        double* tmp = jackall.en[i].jack[33];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[34];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[30];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[31];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[32];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        fprintf(summary_out, "\n");
    }
    fclose(summary_out);
    mysprintf(namefile, NAMESIZE, "%s/data_from_w0.txt", argv[3]);
    summary_out = open_file(namefile, "w+");
    fprintf(summary_out, "ens   a[fm] da[fm]  afpi dafpi  amul damul amus damus  amuc  damuc \n");
    for (int i = 0;i < myen[0].size();i++) {
        size_t lastUnderscore = files[i].find_last_of('_');

        error(lastUnderscore == std::string::npos, 1, "error in file name", "cannot read ens from name %s", files[i].c_str());
        // 2. Estrae tutto ciò che segue l'ultimo underscore
        std::string result = files[i].substr(lastUnderscore + 1);
        fprintf(summary_out, "%s  ", result.c_str());
        double* tmp = jackall.en[i].jack[38];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[39];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[35];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[36];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[37];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        fprintf(summary_out, "\n");
    }


    //////////////////////////////////////////////////////////////
    // fitting
    //////////////////////////////////////////////////////////////
    std::vector<std::string> obs = { "w0" ,"fpi" };
    std::vector<int> id_obs = { 34, 39 };
    std::vector<int> id_a = { 33, 38 };

    std::vector<std::string> obs2 = { "ms_over_ml" ,"mc_over_ms" };
    std::vector<std::vector<int>> id_obs_m = { {31,30}, {32,31}, {36,35}, {37,36} };


    int id_amuliso = 13;
    int id_a_fm = 16;
    int id_w0 = 1;
    int id_mu = 11;
    for (int i = 0; i < obs.size(); i++) {
        fit_type fit_info;

        fit_info.corr_id = { id_obs[i] };

        //////////////////////////////////////////////////////////////
        // charm small volume fit - constant
        //////////////////////////////////////////////////////////////

        fit_info.N = 1;
        fit_info.Nxen = std::vector<std::vector<int>>(fit_info.N);
        for (int n = 0; n < fit_info.N; n++) {
            fit_info.Nxen[n].resize(myen[n].size());
            for (int b = 0; b < myen[n].size(); b++) {
                fit_info.Nxen[n][b] = myen[n][b];  // the charm is in myen[beta][0]
            }
        }
        fit_info.init_N_etot_form_Nxen();

        fit_info.function = rhs_a2;
        fit_info.linear_fit = true;
        fit_info.Npar = 2;
        fit_info.Nvar = 1; // a2
        fit_info.Njack = jackall.en[0].Njack;
        fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);

        int count = 0;
        for (int n = 0; n < fit_info.N; n++) {
            for (int e : fit_info.Nxen[n]) {
                for (int j = 0; j < Njack; j++) {
                    fit_info.x[0][count][j] = pow(jackall.en[e].jack[id_a[i]][j], 2); // a^2 fm^2
                }
                count++;
            }
        }
        // fit_info.linear_fit = false;
        fit_info.verbosity = 0;
        // fit_info.covariancey = true;
        // fit_info.compute_cov_fit(argv, jackall, lhs_fun);
        // fit_info.make_covariance_block_diagonal_in_n();
        // fit_info.compute_cov1_fit();

        std::string namefit = "fit_" + obs[i] + "_cont_a2";

        fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});
        der_fpi_const_full.clear();
        for (size_t imr = 0; imr < obs2.size(); imr++) {
            fit_info.corr_id = id_obs_m[imr + obs2.size() * i];
            std::string namefit = "fit_" + obs2[imr] + "_from_" + obs[(i + 1) % 2] + "_cont_a2";
            fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun_ratio, fit_info, namefit.c_str());
            fit_info.band_range = { 0, 0.008145209846823482 };
            print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});

            der_fpi_const_full.clear();
        }
        fit_info.restore_default();

    }

    //////////////////////////////////////////////////////////////
    // lattice spacings
    //////////////////////////////////////////////////////////////
    {
        fit_type fit_info;

        fit_info.corr_id = { id_a[0], id_a[1] };
        fit_info.N = 1;
        fit_info.Nxen = std::vector<std::vector<int>>(fit_info.N);
        for (int n = 0; n < fit_info.N; n++) {
            fit_info.Nxen[n].resize(myen[n].size());
            for (int b = 0; b < myen[n].size(); b++) {
                fit_info.Nxen[n][b] = myen[n][b];  // the charm is in myen[beta][0]
            }
        }
        fit_info.init_N_etot_form_Nxen();

        fit_info.function = rhs_a2;
        fit_info.linear_fit = true;
        fit_info.Npar = 2;
        fit_info.Nvar = 1; // a2
        fit_info.Njack = jackall.en[0].Njack;
        fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);

        int count = 0;
        for (int n = 0; n < fit_info.N; n++) {
            for (int e : fit_info.Nxen[n]) {
                for (int j = 0; j < Njack; j++) {
                    fit_info.x[0][count][j] = pow(jackall.en[e].jack[id_a[0]][j], 2); // a^2 fm^2
                }
                count++;
            }
        }
        // fit_info.linear_fit = false;
        fit_info.verbosity = 0;
        // fit_info.covariancey = true;
        // fit_info.compute_cov_fit(argv, jackall, lhs_fun);
        // fit_info.make_covariance_block_diagonal_in_n();
        // fit_info.compute_cov1_fit();

        std::string namefit = "fit_afpi_over_aw0_cont_a2";

        fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun_ratio, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});

        der_fpi_const_full.clear();
        namefit = "fit_afpi_minus_aw0_cont_a2";

        der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun_diff, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});

        fit_info.restore_default();


    }
}
