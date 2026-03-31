#define CONTROL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #include <time.h>
// #include <complex.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>

#include "global.hpp"
#include "read.hpp"
#include "resampling_new.hpp"

#include "linear_fit.hpp"
#include "mutils.hpp"
#include "tower.hpp"
#include "various_fits.hpp"

#include "correlators_analysis.hpp"
#include "functions_w0.hpp"
#include <acb_calc.h>
#include <arb.h>
#include <mag.h>

#include "myarb.hpp"

struct kinematic kinematic_2pt;

generic_header read_head(FILE* stream) {
    generic_header header;
    return header;
}
void write_header_g2(FILE* jack_file, generic_header head) {
    int fi = 0;
    fi += fwrite(&head.T, sizeof(int), 1, jack_file);
    fi += fwrite(&head.L, sizeof(int), 1, jack_file);
    int nmu = head.mus.size();
    fi += fwrite(&nmu, sizeof(int), 1, jack_file);
    for (double mu : head.mus) {
        fi += fwrite(&mu, sizeof(double), 1, jack_file);
    }
}

char** argv_to_options(char** argv) {
    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, "read_plateaux"); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p");            // -p
    mysprintf(option[3], NAMESIZE, argv[2]);         // path
    mysprintf(option[4], NAMESIZE, argv[6]);         // resampling
    mysprintf(option[5], NAMESIZE, "no");            // pdf
    mysprintf(option[6], NAMESIZE, argv[3]);         // infile
    return option;
}

void init_global_head(generic_header head) {
    file_head.l1 = head.L;
    file_head.l0 = head.T;
    file_head.l2 = head.L;
    file_head.l3 = head.L;
    file_head.nk = 2;
    file_head.musea = head.mus[0];
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;
    file_head.k[1] = 0;
    file_head.k[2] = head.mus[0];
    file_head.k[3] = head.mus[0];

    file_head.nmoms = 1;
    file_head.mom = (double**)malloc(sizeof(double*) * file_head.nmoms);
    for (int i = 0; i < file_head.nmoms; i++) {
        file_head.mom[i] = (double*)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }
}

void read_twopt(FILE* stream, double*** to_write, generic_header head) {
    // write your function to read the data
    // int fi = 0;
    // for (int k = 0; k < head.ncorr; k++) {
    //     for (int t = 0; t < head.T;t++) {
    //         fi += fscanf(stream, "%lf  %lf", &to_write[k][t][0],
    //         &to_write[k][t][1]);
    //         // printf(" corr=%d t=%d %.12g   %.12g\n", k, t, to_write[k][t][0],
    //         to_write[k][t][1]);
    //     }
    // }
    //// binary
    int fi = 0;
    int id;
    int i = fread(&id, sizeof(int), 1, stream);
    for (int k = 0; k < head.ncorr; k++) {
        for (int t = 0; t < head.T; t++) {
            fi += fread(to_write[k][t], sizeof(double), 2, stream);
        }
    }
}

double int2flowt(double i) {
    return 0.010000 + i * 0.02;
}

double poly3(int n, int Nvar, double* x, int Npar, double* P) {
    double it = x[0];
    double tf = int2flowt(x[0]);
    return P[0] + P[1] * tf + P[2] * tf * tf + P[3] * tf * tf * tf;
}

double lhs_fun(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r = gjack.en[e].jack[fit_info.corr_id[n]][j]; // d(w0/a)/d(a*mu)
    return r;
}

double lin_fit(int n, int Nvar, double* x, int Npar, double* P) {

    return P[0] + P[1] * x[0];
}
double mu2_fit(int n, int Nvar, double* x, int Npar, double* P) {

    return P[0] + P[1] * x[0] + P[2] * x[0] * x[0];
}
double mu_shift(int n, int Nvar, double* x, int Npar, double* P) {

    return  P[0] * (x[0] - x[1]);
}


double mu2_shift(int n, int Nvar, double* x, int Npar, double* P) {

    double dx = x[0] - x[1];
    return  P[0] * dx + P[1] * dx * dx;
}
double mu3_shift(int n, int Nvar, double* x, int Npar, double* P) {

    double dx = x[0] - x[1];
    return  P[0] * dx + P[1] * dx * dx + P[2] * dx * dx * dx;
}
double mu4_shift(int n, int Nvar, double* x, int Npar, double* P) {

    double dx = x[0] - x[1];
    return  P[0] * dx + P[1] * dx * dx + P[2] * dx * dx * dx + P[3] * dx * dx * dx * dx;
}
double pade_2_2_shift(int n, int Nvar, double* x, int Npar, double* P) {

    double dx = x[0] - x[1];
    return  (P[0] * dx + P[1] * dx * dx) / (1 + P[2] * dx + P[3] * dx * dx);
}


double lhs_shift(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    r = gjack.en[e].jack[fit_info.corr_id[n]][j] - gjack.en[0].jack[0][j]; //

    return r;
}
int main(int argc, char** argv) {
    error(argc < 10, 1, "main ",
        "usage:././w0_rew -p path file -bin $bin  jack/boot   reweighting_factors  name_rew  rew_minun_m   [extra rew]\n separate "
        "path and file please");


    int Nrew = argc - 8; // number of reweighting files, here 2 for the w0 and w0_m
    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, argv[6]);
    printf("resampling: %s\n", resampling);

    char** option = argv_to_options(argv);

    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], option[6]);

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");

    FILE* infile = open_file(namefile, "r");

    //////////////////////////////////// read and setup header
    generic_header head;
    head.read_header_debug(infile);
    // head.print_header();
    init_global_head(head);

    //////////////////////////////////////////////////////////////
    // read the reweighting
    //////////////////////////////////////////////////////////////

    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[7]);
    printf("namefile rew_m: %s\n", namefile);
    FILE* infile_rew = open_file(namefile, "r");

    generic_header head_rew;
    head_rew.read_header(infile_rew);
    error(head.Njack != head_rew.Njack, 1, "main", "Najck w0 = %d   while Njack rew  = %d", head.Njack, head_rew.Njack);
    double**** data_rew = calloc_corr(head_rew.Njack, head_rew.ncorr, head_rew.T);
    for (int iconf = 0; iconf < head_rew.Njack; iconf++) {
        read_twopt(infile_rew, data_rew[iconf], head_rew);
        error(head.smearing[iconf].compare(head_rew.smearing[iconf]) != 0, 2, "main",
            "configuration order differ at %d\n flow file conf: %s\n rew conf: %s ", iconf,
            head.smearing[iconf].c_str(), head_rew.smearing[iconf].c_str());
    }


    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[9]);
    printf("namefile rew_m: %s\n", namefile);
    FILE* infile_rew_m = open_file(namefile, "r");

    generic_header head_rew_m;
    head_rew_m.read_header(infile_rew_m);
    error(head.Njack != head_rew_m.Njack, 1, "main", "Najck w0 = %d   while Njack rew  = %d", head.Njack, head_rew_m.Njack);
    double**** data_rew_m = calloc_corr(head_rew_m.Njack, head_rew_m.ncorr, head_rew_m.T);
    for (int iconf = 0; iconf < head_rew_m.Njack; iconf++) {
        read_twopt(infile_rew_m, data_rew_m[iconf], head_rew_m);
        error(head.smearing[iconf].compare(head_rew_m.smearing[iconf]) != 0, 2, "main",
            "configuration order differ at %d\n flow file conf: %s\n rew_m conf: %s ", iconf,
            head.smearing[iconf].c_str(), head_rew_m.smearing[iconf].c_str());
    }

    std::vector<generic_header> head_rew_m_vec(Nrew - 2);
    std::vector<double****> data_rew_m_vec(Nrew - 2);
    for (int i = 2; i < Nrew;i++) {
        mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[8 + i]);
        printf("namefile rew_m: %s\n", namefile);
        FILE* infile_rew_m = open_file(namefile, "r");
        head_rew_m_vec[i - 2].read_header(infile_rew_m);
        data_rew_m_vec[i - 2] = calloc_corr(head_rew_m_vec[i - 2].Njack, head_rew_m_vec[i - 2].ncorr, head_rew_m_vec[i - 2].T);
        for (int iconf = 0; iconf < head_rew_m_vec[i - 2].Njack; iconf++) {
            read_twopt(infile_rew_m, data_rew_m_vec[i - 2][iconf], head_rew_m_vec[i - 2]);
            error(head.smearing[iconf].compare(head_rew_m_vec[i - 2].smearing[iconf]) != 0, 2, "main",
                "configuration order differ at %d\n flow file conf: %s\n rew_m conf: %s ", iconf,
                head.smearing[iconf].c_str(), head_rew_m_vec[i - 2].smearing[iconf].c_str());
        }
    }
    //////////////////////////////////////////////////////////////
    // reading flow
    //////////////////////////////////////////////////////////////
    int ncorr_new = head.ncorr;    // current number of correlators
    int Max_corr = head.ncorr + 4; // max number of correlators
    for (int i = 2; i < Nrew;i++) Max_corr += 2;

    double**** data = calloc_corr(head.Njack, Max_corr, head.T);

    printf("confs=%d\n", head.Njack);
    printf("ncorr=%d\n", head.ncorr);
    printf("kappa=%g\n", head.kappa);
    for (int iconf = 0; iconf < head.Njack; iconf++) {
        read_twopt(infile, data[iconf], head);
    }


    //////////////////////////////////////////////////////////////
    // correcting w0
    //////////////////////////////////////////////////////////////
    // for (int i = 0; i < head_rew.ncorr; i++)
    // {
    double sum_r = 0;
    for (int j = 0; j < head.Njack; j++) {
        // for (int j1 = 0; j1 < head.Njack; j1++)
        // {
        //     r += exp(data_rew[j1][0][0][0] - data_rew[j][0][0][0]);
        // }
        // r = head.Njack / r;
        // sum_r += r;
        // printf("rew %d  %g  %g    -  %g %g\n",i, r,data_rew[j][0][0][1],data_rew[j][0][0][0],exp(data_rew[j][0][0][0]));
        // r = data_rew[j][0][0][1]; // the im part had the normalized reweighting factor
        // printf("rew %d %d %g\n", i, j, r);

        double r = data_rew[j][0][0][1]; // the im part is the exponentiated subtracted reweighting factor
        double rm = data_rew_m[j][0][0][1]; // the im part is the exponentiated subtracted reweighting factor
        for (int tf = 0; tf < head.T; tf++) {
            data[j][head.ncorr + 0][tf][0] = data[j][6][tf][0] * r;
            data[j][head.ncorr + 0][tf][1] = data[j][6][tf][1] * r;
            data[j][head.ncorr + 1][tf][0] = data_rew[j][0][0][1];
            data[j][head.ncorr + 1][tf][1] = 0;
            data[j][head.ncorr + 2][tf][0] = data[j][6][tf][0] * rm;
            data[j][head.ncorr + 2][tf][1] = data[j][6][tf][1] * rm;
            data[j][head.ncorr + 3][tf][0] = data_rew_m[j][0][0][1];
            data[j][head.ncorr + 3][tf][1] = 0;
            for (int i = 2; i < Nrew;i++) {
                double rv = data_rew_m_vec[i - 2][j][0][0][1]; // the im part is the exponentiated subtracted reweighting factor
                data[j][head.ncorr + i * 2][tf][0] = data[j][6][tf][0] * rv;
                data[j][head.ncorr + i * 2 + 1][tf][0] = rv;
            }
        }
    }
    // printf("sum_r = %g\n",sum_r/head.Njack);
    // }
    // free_corr(head_rew.Njack, head_rew.ncorr, head_rew.T, data_rew);

    //////////////////////////////////////////////////////////////
    // binning and resampling
    //////////////////////////////////////////////////////////////
    //////////////////////////////////// setup jackboot and binning
    int confs = head.Njack;
    int bin = atoi(argv[5]);
    // int Neff = confs / bin; // standard binning
    int Neff = bin; // bin2N
    int Njack;
    if (strcmp(argv[6], "jack") == 0) {
        Njack = Neff + 1;
        myres = new resampling_jack(Neff);
    }
    else if (strcmp(argv[6], "boot") == 0) {
        Njack = (Neff * 2 + 1);
        myres = new resampling_boot(Neff * 2);
    }
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[7]= %s is not jack or boot", argv[7]);
    }
    // now Njack need to be the number of jacks
    head.Nconf = head.Njack;
    head.Njack = Njack;
    //////////////////////////////////// setup output files
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_%s_output", option[3], option[6], argv[7], argv[9]);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s_%s_%s", option[3], option[4], option[6], argv[7], argv[9]);
    printf("writing output in :\n %s \n", namefile);
    FILE* jack_file = open_file(namefile, "w+");
    // write_header_g2(jack_file, head);
    head.write_header(jack_file);

    // double ****data_bin = binning(confs, Max_corr, head.T, data, bin);
    // double ****data_bin = binning_toNb(confs, Max_corr, head.T, data, bin); // binning into N=bin, cutting out the reminder
    // double ****data_bin = bin_intoN(data, Max_corr, head.T, confs, bin); // binning into N=bin with not integer
    // free_corr(confs, Max_corr, head.T, data);
    // double ****conf_jack = myres->create(Neff, Max_corr, head.T, data_bin);
    // free_corr(Neff, Max_corr, head.T, data_bin);

    //////////////////////////////////////////////////////////////
    // making custom jackknifes
    //////////////////////////////////////////////////////////////
    double**** data_bin = bin_intoN(data, Max_corr, head.T, confs, bin); // binning into N=bin with not integer
    double**** conf_jack = myres->create(Neff, Max_corr, head.T, data_bin);
    free_corr(Neff, Max_corr, head.T, data_bin);
    free_corr(head_rew.Njack, head_rew.ncorr, head_rew.T, data_rew);

    // <Or>/<r>
    for (int j = 0; j < head.Njack; j++) {

        for (int tf = 0; tf < head.T; tf++) {
            conf_jack[j][head.ncorr + 0][tf][0] /= conf_jack[j][head.ncorr + 1][tf][0];
            conf_jack[j][head.ncorr + 2][tf][0] /= conf_jack[j][head.ncorr + 3][tf][0];
            for (int i = 2; i < Nrew;i++) {
                conf_jack[j][head.ncorr + i * 2][tf][0] /= conf_jack[j][head.ncorr + i * 2 + 1][tf][0];
            }
        }
    }

    printf("Njack = %d  Nbins = %d  confs = %d\n", head.Njack, bin, confs);


    free_corr(head.Njack, 6, head.T, data);

    printf("corr at time 10 jacks\n");
    for (int j = 0; j < head.Njack; j++) {
        printf("%d %.12g  %.12g %.12g\n", j, conf_jack[j][head.ncorr + 0][10][0], conf_jack[j][head.ncorr + 1][0][0], conf_jack[j][6][10][0]);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print all the effective masses correlators
    // set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_%s_meff_correlators", option[3], option[6], argv[7], argv[9]);
    FILE* outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_%s_raw_correlators", option[3], option[6], argv[7], argv[9]);
    FILE* outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_%s_shifted_correlators", option[3], option[6], argv[7], argv[9]);
    FILE* outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_%s_log_meff_shifted", option[3], option[6], argv[7], argv[9]);
    FILE* outfile_log_meff_shifted = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_%s_gamma", option[3], option[6], argv[7], argv[9]);
    FILE* out_gamma = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_%s_HLT_kernel", option[3], option[6], argv[7], argv[9]);
    FILE* outfile_HLT_kernel = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_%s_HLT_AoverB", option[3], option[6], argv[7], argv[9]);
    FILE* outfile_HLT_AoverB = open_file(namefile, "w+");

    char save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");
    FILE* dev_null = open_file("/dev/null", "w");
    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;

    for (int icorr = 0; icorr < head.ncorr; icorr++) {
        // log effective mass
        double* tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_meff_corr, icorr, "log", M_eff_log, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        // raw correlator
        file_head.l0 = head.T * 2;
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity_im,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        file_head.l0 = head.T;
        // shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_shifted_corr, icorr, "shift_cor", shift_corr,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_log_meff_shifted, icorr, "log_shift",
            M_eff_log_shift, dev_null, fit_info_silent);
        free(tmp_meff_corr);
    }
    fit_info_silent.restore_default();
    sprintf(option[1], "%s", save_option); // restore option
    corr_counter = -1;

    ////////////////////////////////////////////////////////////
    // symmetrize
    ////////////////////////////////////////////////////////////
    // symmetrise_jackboot(Njack, 0, head.T, conf_jack);
    // symmetrise_jackboot(Njack, 1, head.T, conf_jack, -1);

    ////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////
    corr_counter = -1;
    file_head.l0 = head.T * 2;  // dou

    std::vector<int> id_w0 = { 6, head.ncorr + 0, head.ncorr + 2 };
    std::vector<std::string> name_w0 = { "w0", "w0_rew", "w0_rew_m" };
    std::vector<double> mus = { head_rew.mus[0], head_rew.oranges[0], head_rew_m.oranges[0] };

    for (int i = 2; i < Nrew;i++) {
        id_w0.push_back(head.ncorr + i * 2);
        name_w0.push_back("w0_rew_m_" + std::to_string(i - 2));
        mus.push_back(head_rew_m_vec[i - 2].oranges[0]);
    }
    double** w0 = malloc_2<double>(id_w0.size(), Njack);

    for (int i = 0;i < id_w0.size();i++) {

        // eg of fit to correlator
        struct fit_type fit_info;
        struct fit_result fit_out;

        fit_info.Nvar = 1;
        fit_info.Npar = 4;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.codeplateaux = true;

        fit_info.function = poly3;
        fit_info.linear_fit = true;
        fit_info.T = head.T * 2;
        fit_info.corr_id = { id_w0[i] };

        for (int t = 1; t < head.T; t++) {
            if (lhs_function_w0_eg(Njack - 1, conf_jack, t, fit_info) > 0.3) {
                fit_info.tmin = t - 2;
                fit_info.tmax = t + 1;
                break;
            }
        }

        struct fit_result fit_W = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_function_w0_eg, "W(t)", fit_info,
            jack_file);
        check_correlatro_counter(0 + i * 2);

        double** tif = swap_indices(fit_info.Npar, Njack, fit_W.P);
        std::vector<double> swapped_x(fit_info.Nvar);

        for (size_t j = 0; j < Njack; j++) {
            w0[i][j] = rtbis_func_eq_input(fit_info.function, 0 /*n*/, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[j], 0, 0.3, fit_info.tmin, fit_info.tmax, 1e-10, 2);
            w0[i][j] = int2flowt(w0[i][j]);
            w0[i][j] = std::sqrt(w0[i][j]);
        }
        printf("w0 = %.12g  %.12g\n", w0[i][Njack - 1], myres->comp_error(w0[i]));
        write_jack(w0[i], Njack, jack_file);
        check_correlatro_counter(1 + i * 2);




        print_result_in_file(outfile, w0[i], name_w0[i].c_str(), 0.0, fit_info.tmin, fit_info.tmax);
        free_2(Njack, tif);
    }
    // // i want to sort the w0 by the mu, so I need to sort the id_w0, name_w0, mus and w0 together
    // for (size_t i = 0; i < id_w0.size(); i++) {
    //     for (size_t j = i + 1; j < id_w0.size(); j++) {
    //         if (mus[j] < mus[i]) {
    //             std::swap(mus[i], mus[j]);
    //             std::swap(id_w0[i], id_w0[j]);
    //             std::swap(name_w0[i], name_w0[j]);
    //             std::swap(w0[i], w0[j]);
    //         }
    //     }
    // }
    // printf("sorted mus and w0\n");
    for (size_t i = 0; i < id_w0.size(); i++) {
        printf("%g   %.12g  %.12g\n", mus[i], w0[i][Njack - 1], myres->comp_error(w0[i]));
    }

    //////////////////////////////////////////////////////////////
    // reading onlinemeas parameters
    //////////////////////////////////////////////////////////////
    std::vector<double*> amuiso(3);
    double mean, err;
    int seed;
    std::string name_in(option[6]);
    std::string name_online(option[6]);
    name_online = std::regex_replace(name_online, std::regex("flow"), "onlinemeas");
    mysprintf(option[6], NAMESIZE, "%s", name_online.c_str());

    line_read_param(option, "muliso", mean, err, seed, namefile_plateaux);
    amuiso[0] = myres->create_fake(mean, err, seed);
    line_read_param(option, "musiso", mean, err, seed, namefile_plateaux);
    amuiso[1] = myres->create_fake(mean, err, seed);
    line_read_param(option, "muciso", mean, err, seed, namefile_plateaux);
    amuiso[2] = myres->create_fake(mean, err, seed);

    std::vector<double*> amusim(3);
    line_read_param(option, "mulsim", mean, err, seed, namefile_plateaux);
    amusim[0] = myres->create_fake(mean, err, seed);
    line_read_param(option, "mussim", mean, err, seed, namefile_plateaux);
    amusim[1] = myres->create_fake(mean, err, seed);
    line_read_param(option, "mucsim", mean, err, seed, namefile_plateaux);
    amusim[2] = myres->create_fake(mean, err, seed);

    line_read_param(option, "a", mean, err, seed, namefile_plateaux);
    double* a_fm = myres->create_fake(mean, err, seed);

    mysprintf(option[6], NAMESIZE, "%s", name_in.c_str());

    double* der_f = (double*)malloc(sizeof(double) * Njack);
    double* der_b = (double*)malloc(sizeof(double) * Njack);
    double* der_sym = (double*)malloc(sizeof(double) * Njack);
    double* der_2 = (double*)malloc(sizeof(double) * Njack);

    //// end reading onlinemeas parameters
    double dmu_f = mus[1] - mus[0];
    double dmu_b = mus[0] - mus[2];
    double dmu_sym = mus[1] - mus[2];
    error(std::fabs(dmu_f - dmu_b) > 1e-10, 1, "main", "unexpected mu spacing for 3 points, got dmu_f = %.12g   dmu_b = %.12g", dmu_f, dmu_b);
    printf("mus = %g   %g  %g\n", mus[0], mus[1], mus[2]);
    printf("dmu_f = %g  dmu_b = %g\n", dmu_f, dmu_b);
    for (size_t j = 0; j < Njack; j++) {
        der_f[j] = (w0[1][j] - w0[0][j]) / dmu_f;
        der_b[j] = (w0[0][j] - w0[2][j]) / dmu_b;
        der_sym[j] = (w0[1][j] - w0[2][j]) / dmu_sym;
    }
    printf("der_f = %g  %g\n", der_f[Njack - 1], myres->comp_error(der_f));
    printf("der_b = %g  %g\n", der_b[Njack - 1], myres->comp_error(der_b));
    printf("der_sym = %g  %g\n", der_sym[Njack - 1], myres->comp_error(der_sym));

    double* der_h4 = (double*)malloc(sizeof(double) * Njack);
    if (id_w0.size() > 4) {
        for (size_t j = 0; j < Njack; j++) {
            der_b[j] = (w0[2][j] - w0[4][j]) / dmu_b;
        }
        printf("der_bb = %g  %g\n", der_b[Njack - 1], myres->comp_error(der_b));
        for (size_t j = 0; j < Njack; j++) {
            der_b[j] = (w0[3][j] - w0[1][j]) / dmu_b;
        }
        printf("der_ff = %g  %g\n", der_b[Njack - 1], myres->comp_error(der_b));

        for (size_t j = 0; j < Njack; j++) {
            der_h4[j] = (-w0[3][j] + 8.0 * w0[1][j] - 8.0 * w0[2][j] + w0[4][j]) / (12 * dmu_f);
        }
        printf("der_h4 = %g  %g\n", der_h4[Njack - 1], myres->comp_error(der_h4));

    }

    for (size_t j = 0; j < Njack; j++) {
        der_2[j] = (w0[1][j] - 2.0 * w0[0][j] + w0[2][j]) / (dmu_f * dmu_b);
    }

    printf("der_2 = %g  %g\n", der_2[Njack - 1], myres->comp_error(der_2));

    double* der_2b = (double*)malloc(sizeof(double) * Njack);
    double* der_2f = (double*)malloc(sizeof(double) * Njack);
    double* der_2h4 = (double*)malloc(sizeof(double) * Njack);
    double* der_3b = (double*)malloc(sizeof(double) * Njack);
    double* der_4b = (double*)malloc(sizeof(double) * Njack);
    if (id_w0.size() > 4) {
        for (size_t j = 0; j < Njack; j++) {
            der_2b[j] = (w0[0][j] - 2.0 * w0[2][j] + w0[4][j]) / (dmu_f * dmu_b);
        }
        printf("der_2b = %g  %g\n", der_2b[Njack - 1], myres->comp_error(der_2b));
        for (size_t j = 0; j < Njack; j++) {
            der_2f[j] = (w0[3][j] - 2.0 * w0[1][j] + w0[0][j]) / (dmu_f * dmu_b);
        }
        printf("der_2f = %g  %g\n", der_2f[Njack - 1], myres->comp_error(der_2f));
        // do it with Arb 
        arb_t fph, fmh, f, dm, arb_der;
        arb_init(fph);
        arb_init(fmh);
        arb_init(f);
        arb_init(dm);
        arb_init(arb_der);

        int prec = 256;
        arb_set_d(dm, dmu_b);
        for (size_t j = 0; j < Njack; j++) {
            arb_set_d(fph, w0[3][j]);
            arb_set_d(fmh, w0[0][j]);
            arb_set_d(f, w0[1][j]);
            arb_sub(arb_der, fph, f, prec);
            arb_sub(arb_der, arb_der, f, prec);
            arb_add(arb_der, arb_der, fmh, prec);
            arb_div(arb_der, arb_der, dm, prec);
            arb_div(arb_der, arb_der, dm, prec);
            der_2f[j] = arbtod(arb_der);
        }
        printf("der_2f(ARB) = %g  %g\n", der_2f[Njack - 1], myres->comp_error(der_2f));

        for (size_t j = 0; j < Njack; j++) {
            der_2h4[j] = (-w0[3][j] + 16.0 * w0[1][j] - 30.0 * w0[0][j] + 16.0 * w0[2][j] - w0[4][j]) / (12.0 * dmu_f * dmu_b);
        }
        printf("der_2h4 = %g  %g\n", der_2h4[Njack - 1], myres->comp_error(der_2h4));

        double dmu_b2 = mus[2] - mus[4];
        double dmu_2 = mus[3] - mus[1];
        error(std::fabs(dmu_b2 - dmu_b) > 1e-10, 1, "main", "unexpected mu spacing for 5 points, got %g expected %g", dmu_b2, dmu_b);
        error(std::fabs(dmu_2 - dmu_f) > 1e-10, 1, "main", "unexpected mu spacing for 5 points, got %g expected %g", dmu_2, dmu_f);
        for (size_t j = 0; j < Njack; j++) {
            der_3b[j] = (w0[3][j] - 2.0 * w0[1][j] + 2.0 * w0[2][j] - w0[4][j]) / (dmu_b * dmu_b * dmu_b);
        }
        printf("der_3b = %g  %g\n", der_3b[Njack - 1], myres->comp_error(der_3b));

        for (size_t j = 0; j < Njack; j++) {
            der_4b[j] = (w0[3][j] - 4.0 * w0[1][j] + 6 * w0[0][j] - 4.0 * w0[2][j] - w0[4][j]) / (dmu_b * dmu_b * dmu_b * dmu_b);
        }
        printf("der_4b = %g  %g\n", der_4b[Njack - 1], myres->comp_error(der_4b));
    }



    printf("fitting  %ld ensambles\n", id_w0.size());

    data_all jackall;
    jackall.resampling = argv[6];
    jackall.ens = id_w0.size();
    jackall.en = new data_single[jackall.ens];
    for (int i = 0;i < jackall.ens;i++) {
        jackall.en[i].header = head;
        jackall.en[i].Nobs = 1;
        jackall.en[i].Njack = head.Njack;
        jackall.en[i].jack = (double**)malloc(sizeof(double*) * jackall.en[i].Nobs);
        jackall.en[i].jack[0] = w0[i];
    }


    struct fit_type fit_info;
    struct fit_result fit_out;

    fit_info.Npar = 2;
    fit_info.Nvar = 1;
    fit_info.Njack = head.Njack;
    fit_info.N = 1;
    fit_info.Nxen = std::vector<std::vector<int>>(fit_info.N);
    for (int n = 0; n < fit_info.N; n++) {
        fit_info.Nxen[n].resize(id_w0.size());
        for (int b = 0; b < id_w0.size(); b++) {
            fit_info.Nxen[n][b] = b;  // the charm is in myen[beta][0]
        }
    }
    fit_info.init_N_etot_form_Nxen();

    fit_info.malloc_x();
    int count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        // for (int e : fit_info.myen) {
        for (int e = 0; e < fit_info.Nxen[n].size(); e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][count][j] = mus[e];
            }
            count++;
        }
    }
    fit_info.corr_id = { 0 };
    fit_info.function = lin_fit;
    fit_info.linear_fit = true;
    fit_info.verbosity = 3;
    char** option_cov = malloc_2<char>(2, NAMESIZE);
    mysprintf(option_cov[1], NAMESIZE, resampling);
    // fit_info.compute_cov_fit(option_cov, jackall, lhs_identity);
    // fit_info.compute_cov1_fit();
    fit_info.verbosity = 0;
    char namefit[NAMESIZE];
    mysprintf(namefit, NAMESIZE, "w0_charm_lin");
    char** temp_argv = malloc_2<char>(5, NAMESIZE);
    mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
    mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// outdir

    fit_result fit_inter_alpha = fit_all_data(temp_argv, jackall, lhs_identity, fit_info, namefit);
    fit_info.band_range = { 0.22856, 0.23134 + 1e-4 };
    print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "mu", fit_inter_alpha, fit_inter_alpha, 0, fit_info.Nxen[0][0], 5e-6);

    fit_info.function = mu2_fit;
    fit_info.linear_fit = true;
    fit_info.Npar = 3;
    fit_info.precision_sum = 2;
    fit_info.acc = 1e-26;
    fit_info.h = 1e-12;
    fit_info.compute_cov_fit(option_cov, jackall, lhs_identity);
    fit_info.compute_cov1_fit();
    fit_info.covariancey = true;
    mysprintf(namefit, NAMESIZE, "w0_charm_mu2");
    fit_result fit_mu2 = fit_all_data(temp_argv, jackall, lhs_identity, fit_info, namefit);
    print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "mu", fit_mu2, fit_mu2, 0, fit_info.Nxen[0][0], 5e-6);


    //////////////////////////////////////////////////////////////
    // diff
    //////////////////////////////////////////////////////////////

    data_all jackall_m;
    jackall_m.resampling = argv[6];
    jackall_m.ens = id_w0.size() - 1;
    jackall_m.en = new data_single[jackall_m.ens];
    for (int i = 0;i < jackall_m.ens;i++) {
        jackall_m.en[i].header = head;
        jackall_m.en[i].Nobs = 2;
        jackall_m.en[i].Njack = head.Njack;
        jackall_m.en[i].jack = malloc_2<double>(jackall_m.en[i].Nobs, jackall_m.en[i].Njack);
    }
    for (int j = 0; j < Njack;j++) {
        for (int i = 0;i < jackall_m.ens;i++) {
            jackall_m.en[i].jack[0][j] = (w0[i + 1][j] - w0[0][j]);
            jackall_m.en[i].jack[1][j] = mus[i + 1];
        }
    }
    for (int i = 0;i < jackall_m.ens;i++) {
        printf("%.12g  %.12g   %.12g\n", jackall_m.en[i].jack[1][Njack - 1], jackall_m.en[i].jack[0][Njack - 1], myres->comp_error(jackall_m.en[i].jack[0]));
    }
    for (int n = 0; n < fit_info.N; n++) {
        fit_info.Nxen[n].resize(jackall_m.ens);
        for (int b = 0; b < jackall_m.ens; b++) {
            fit_info.Nxen[n][b] = b;  // the charm is in myen[beta][0]
        }
    }
    fit_info.init_N_etot_form_Nxen();

    fit_info.Nvar = 2;
    fit_info.malloc_x();
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        // for (int e : fit_info.myen) {
        for (int e = 0; e < fit_info.Nxen[n].size(); e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][count][j] = mus[e + 1];
                fit_info.x[1][count][j] = mus[0];
            }
            count++;
        }
    }

    fit_info.function = mu_shift;
    fit_info.linear_fit = true;
    fit_info.Npar = 1;
    fit_info.precision_sum = 2;
    mysprintf(namefit, NAMESIZE, "w0_charm_shift_mu");
    // fit_info.compute_cov_fit(option_cov, jackall, lhs_identity);
    // fit_info.compute_cov1_fit();

    fit_result fit_mu_shift = fit_all_data(temp_argv, jackall_m, lhs_identity, fit_info, namefit);
    print_fit_band(temp_argv, jackall_m, fit_info, fit_info, namefit, "mu", fit_mu_shift, fit_mu_shift, 0, fit_info.Nxen[0][0], 1e-6);

    fit_info.function = mu2_shift;
    fit_info.linear_fit = true;
    fit_info.Npar = 2;
    mysprintf(namefit, NAMESIZE, "w0_charm_shift_mu2");


    fit_result fit_mu2_shift = fit_all_data(temp_argv, jackall_m, lhs_identity, fit_info, namefit);
    print_fit_band(temp_argv, jackall_m, fit_info, fit_info, namefit, "mu", fit_mu2_shift, fit_mu2_shift, 0, fit_info.Nxen[0][0], 1e-6);


    fit_info.function = mu3_shift;
    fit_info.linear_fit = true;
    fit_info.Npar = 3;
    mysprintf(namefit, NAMESIZE, "w0_charm_shift_mu3");


    fit_result fit_mu3_shift = fit_all_data(temp_argv, jackall_m, lhs_identity, fit_info, namefit);
    print_fit_band(temp_argv, jackall_m, fit_info, fit_info, namefit, "mu", fit_mu3_shift, fit_mu3_shift, 0, fit_info.Nxen[0][0], 1e-6);


    fit_info.function = mu4_shift;
    fit_info.linear_fit = true;
    fit_info.Npar = 4;
    mysprintf(namefit, NAMESIZE, "w0_charm_shift_mu4");


    fit_result fit_mu4_shift = fit_all_data(temp_argv, jackall_m, lhs_identity, fit_info, namefit);
    print_fit_band(temp_argv, jackall_m, fit_info, fit_info, namefit, "mu", fit_mu4_shift, fit_mu4_shift, 0, fit_info.Nxen[0][0], 1e-6);

    fit_info.function = pade_2_2_shift;
    fit_info.linear_fit = false;
    fit_info.Npar = 4;
    mysprintf(namefit, NAMESIZE, "w0_charm_pade_2_2_shift_mu");


    fit_result fit_pade_2_2_shift = fit_all_data(temp_argv, jackall_m, lhs_identity, fit_info, namefit);
    print_fit_band(temp_argv, jackall_m, fit_info, fit_info, namefit, "mu", fit_pade_2_2_shift, fit_pade_2_2_shift, 0, fit_info.Nxen[0][0], 1e-6);

}