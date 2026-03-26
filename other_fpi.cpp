#define CONTROL
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
// #include <time.h>
// #include <complex.h>
#include <fstream>
#include <iostream>
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
#include "gamma_analysis.hpp"
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
    fi = fread(&id, sizeof(int), 1, stream);
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
    // double it = x[0];
    double tf = int2flowt(x[0]);
    return P[0] + P[1] * tf + P[2] * tf * tf + P[3] * tf * tf * tf;
}

int main(int argc, char** argv) {
    error(!(argc == 11), 1, "main ",
        "usage:././w0_rew -p path file -bin $bin  jack/boot   fileP5P5  fileA0P5_mu1  fileP5P5_mu2   basename_plateau\n separate "
        "path and file please   argc =%d", argc);

    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, argv[6]);
    printf("resampling: %s\n", resampling);

    char** option = argv_to_options(argv);


    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");




    //////////////////////////////////// read and setup header
    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], option[6]);
    FILE* infile_A0 = open_file(namefile, "r");
    generic_header head_A0;
    head_A0.read_header_debug(infile_A0);
    init_global_head(head_A0);

    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[7]);
    FILE* infile_P5 = open_file(namefile, "r");
    generic_header head_P5;
    head_P5.read_header_debug(infile_P5);
    // init_global_head(head_P5);

    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[8]);
    FILE* infile_A0_mu = open_file(namefile, "r");
    generic_header head_A0_mu;
    head_A0_mu.read_header_debug(infile_A0_mu);
    // init_global_head(head_A0_mu);

    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[9]);
    FILE* infile_P5_mu = open_file(namefile, "r");
    generic_header head_P5_mu;
    head_P5_mu.read_header_debug(infile_P5_mu);
    // init_global_head(head_P5);

    error(head_A0_mu.mus.size() != 2, 1, "main", "derivative file %s should have 2 mus, instead %d\n", argv[8], head_A0_mu.mus.size());
    error(head_P5_mu.mus.size() != 2, 1, "main", "derivative file %s should have 2 mus, instead %d\n", argv[9], head_P5_mu.mus.size());

    error(head_A0.T != head_P5.T, 1, "main", "A0 and P5 files should have the same T, instead %d and %d\n", head_A0.T, head_P5.T);
    error(head_A0_mu.T != head_P5_mu.T, 1, "main", "A0_mu and P5_mu files should have the same T, instead %d and %d\n", head_A0_mu.T, head_P5_mu.T);
    error(head_A0.T != head_P5_mu.T, 1, "main", "A0 and P5_mu files should have the same T, instead %d and %d\n", head_A0.T, head_P5_mu.T);
    //////////////////////////////////////////////////////////////
    // reading flow
    //////////////////////////////////////////////////////////////
    error(head_A0.Njack != head_P5.Njack, 1, "main", "A0 file does not have the same confs of P5 file\n");
    int ncorr_new = head_A0.ncorr + 3;        // current number of correlators
    int Max_corr = head_A0.ncorr + 3; // max number of correlators

    double**** data = calloc_corr(head_A0.Njack, Max_corr, head_A0.T);
    double*** d_tmp = malloc_3<double>(1, head_A0.T, 2);

    for (int iconf = 0; iconf < head_A0.Njack; iconf++) {
        read_twopt(infile_A0, data[iconf], head_A0);
        read_twopt(infile_P5, d_tmp, head_P5);
        for (int t = 0; t < head_A0.T; t++)
            data[iconf][1][t][0] = d_tmp[0][t][0];
    }

    double**** data_mu = calloc_corr(head_A0_mu.Njack, 2, head_A0_mu.T);
    double*** d_tmp_mu = malloc_3<double>(1, head_A0_mu.T, 2);
    for (int iconf = 0; iconf < head_A0_mu.Njack; iconf++) {
        read_twopt(infile_A0_mu, data_mu[iconf], head_A0_mu);
        read_twopt(infile_P5_mu, d_tmp_mu, head_P5_mu);
        for (int t = 0; t < head_A0_mu.T; t++)
            data_mu[iconf][1][t][0] = d_tmp_mu[0][t][0];
    }

    free_3(1, head_A0.T, d_tmp);
    free_3(1, head_A0.T, d_tmp_mu);

    //////////////////////////////////////////////////////////////
    // binning and resampling
    //////////////////////////////////////////////////////////////
    //////////////////////////////////// setup jackboot and binning
    int confs = head_A0.Njack;
    int bin = atoi(argv[5]);
    // int Neff = confs / bin; // standard binning
    int Neff = bin; // bin2N
    int Njack;
    if (strcmp(argv[6], "jack") == 0) {
        Njack = Neff + 1;
        myres = new resampling_jack(Neff);
    }
    else if (strcmp(argv[6], "boot") == 0) {
        Neff = 1000;
        Njack = Neff + 1;//(Neff * 2 + 1);
        myres = new resampling_boot(Neff);
    }
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[7]= %s is not jack or boot", argv[7]);
    }
    // now Njack need to be the number of jacks
    head_A0.Nconf = head_A0.Njack;
    head_A0.Njack = Njack;


    //////////////////////////////////// setup output files
    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", option[3], option[6], argv[7]);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s", option[3], option[4], option[6]);
    FILE* jack_file = open_file(namefile, "w+");
    // write_header_g2(jack_file, head);
    head_A0.write_header(jack_file);



    FILE* dev_null = open_file("/dev/null", "w");

    char** option_r = argv_to_options(argv);
    mysprintf(option_r[6], NAMESIZE, "/dev/null");





    //////////////////////////////////////////////////////////////
    // making custom jackknifes
    //////////////////////////////////////////////////////////////
    double**** data_bin = bin_intoN(data, Max_corr, head_A0.T, confs, bin); // binning into N=bin with not integer
    double**** conf_jack = myres->create(Neff, Max_corr, head_A0.T, data_bin);
    free_corr(Neff, Max_corr, head_A0.T, data_bin);
    free_corr(confs, Max_corr, head_A0.T, data);

    int confs_mu = head_A0_mu.Njack;
    head_A0_mu.Njack = Njack;
    double**** data_bin_mu = bin_intoN(data_mu, 2, head_A0_mu.T, confs_mu, bin); // binning into N=bin with not integer
    double**** conf_jack_mu = myres->create(Neff, 2, head_A0_mu.T, data_bin_mu);
    free_corr(Neff, 2, head_A0_mu.T, data_bin_mu);
    free_corr(confs_mu, 2, head_A0_mu.T, data_mu);


    /// normalize rew factor
    error(head_A0.mus[0] != head_P5.mus[0] != 0, 1, "error mu do not mathch ", "file  A0   mu =%g \nfile  P5   mu =%g ", head_A0.mus[0], head_P5.mus[0]);
    // head_A0.mus[0] = sim
    // head_A0_mu.mus[0] = extra
    double dmu = head_A0.mus[0] - head_A0_mu.mus[1];
    printf("mu 0 = %g\n", head_A0.mus[0]);
    printf("mu 1 = %g\n", head_A0_mu.mus[1]);
    printf("dmu = %g\n", dmu);
    for (int j = 0; j < head_A0.Njack; j++) {
        for (int tf = 0; tf < head_A0.T; tf++) {
            // conf_jack[j][2][tf][0] = conf_jack_mu[j][0][tf][0];
            // conf_jack[j][2][tf][1] = conf_jack_mu[j][0][tf][1];
            conf_jack[j][2][tf][0] = conf_jack[j][0][tf][0] - conf_jack_mu[j][0][tf][0] * dmu;
            conf_jack[j][2][tf][1] = conf_jack[j][0][tf][1] - conf_jack_mu[j][0][tf][1] * dmu;

            conf_jack[j][3][tf][0] = conf_jack[j][1][tf][0] - conf_jack_mu[j][1][tf][0] * dmu;
            conf_jack[j][3][tf][1] = conf_jack[j][1][tf][1] - conf_jack_mu[j][1][tf][1] * dmu;
        }
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print all the effective masses correlators
    // set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_meff_correlators", option[3], option[6], argv[7]);
    FILE* outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_raw_correlators", option[3], option[6], argv[7]);
    FILE* outfile_raw_corr = open_file(namefile, "w+");
    printf("writing in:\n %s \n", namefile);
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_shifted_correlators", option[3], option[6], argv[7]);
    FILE* outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_log_meff_shifted", option[3], option[6], argv[7]);
    FILE* outfile_log_meff_shifted = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_HLT_kernel", option[3], option[6], argv[7]);
    FILE* outfile_HLT_kernel = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_HLT_AoverB", option[3], option[6], argv[7]);
    FILE* outfile_HLT_AoverB = open_file(namefile, "w+");

    char save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");

    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;

    for (int icorr = 0; icorr < Max_corr; icorr++) {
        // log effective mass
        double* tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_meff_corr, icorr, "log", M_eff_log, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        // raw correlator
        file_head.l0 = head_A0.T * 2;
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
        file_head.l0 = head_A0.T;
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

    //////////////////////////////////////////////////////////////
    //  read m^iso
    //////////////////////////////////////////////////////////////
    std::vector<double*> amuiso(3);
    double mean, err;
    int seed;
    char** option_p = argv_to_options(argv);
    mysprintf(option_p[6], NAMESIZE, "%s", argv[10]);         // basename plateau
    line_read_param(option_p, "muliso", mean, err, seed, namefile_plateaux);
    amuiso[0] = myres->create_fake(mean, err, seed);
    line_read_param(option_p, "musiso", mean, err, seed, namefile_plateaux);
    amuiso[1] = myres->create_fake(mean, err, seed);
    line_read_param(option_p, "muciso", mean, err, seed, namefile_plateaux);
    amuiso[2] = myres->create_fake(mean, err, seed);

    std::vector<double*> amusim(3);
    line_read_param(option_p, "mulsim", mean, err, seed, namefile_plateaux);
    amusim[0] = myres->create_fake(mean, err, seed);
    line_read_param(option_p, "mussim", mean, err, seed, namefile_plateaux);
    amusim[1] = myres->create_fake(mean, err, seed);
    line_read_param(option_p, "mucsim", mean, err, seed, namefile_plateaux);
    amusim[2] = myres->create_fake(mean, err, seed);

    line_read_param(option_p, "a", mean, err, seed, namefile_plateaux);
    double* a_fm = myres->create_fake(mean, err, seed);

    std::string name_OStm = argv[3];
    std::string name_OStm1 = argv[7];

    size_t first = name_OStm.find('_');
    size_t second = name_OStm.find('_', first + 1);
    std::string label = name_OStm.substr(first + 1, second - first - 1);

    size_t first1 = name_OStm1.find('_');
    size_t second1 = name_OStm1.find('_', first1 + 1);
    std::string label1 = name_OStm1.substr(first1 + 1, second1 - first1 - 1);

    std::string Z_RIMON;
    if (label.compare("tm") == 0 && label1.compare("tm") == 0)
        Z_RIMON = "ZV_RIMOM";
    else if (label.compare("OS") == 0 && label1.compare("OS") == 0)
        Z_RIMON = "ZA_RIMOM";
    else
        error(1 == 1, 1, "main", "label of input files should be  both tm or OS");

    printf("Z_RIMON: %s\n", Z_RIMON.c_str());
    line_read_param(option_p, Z_RIMON.c_str(), mean, err, seed, namefile_plateaux);
    double* Z = myres->create_fake(mean, err, -1);


    ////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////
    corr_counter = -1;

    double* M_PS = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, 1, "M_{PS}", M_eff_T, jack_file);
    // free(M_PS);
    check_correlatro_counter(0);

    struct fit_type fit_info;
    fit_info.codeplateaux = true;
    line_read_param(option, "M_{PS}", fit_info.tmin, fit_info.tmax, seed, namefile_plateaux);

    double* M_PS_mu = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, 3, "M_{PS}_mu", M_eff_T, jack_file, fit_info);
    // free(M_PS);
    check_correlatro_counter(1);

    {
        double* deriv = myres->create_copy(M_PS);
        for (int j = 0; j < Njack;j++) {
            deriv[j] = (M_PS[j] - M_PS_mu[j]) / dmu;
        }
        printf("deriv: %.12g   %.12g\n", deriv[Njack - 1], myres->comp_error(deriv));
        free(deriv);
    }

    //////////////// me  and fpi

    struct fit_result fit_out;

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = head_A0.T;
    fit_info.n_ext_P = 1;
    // fit_info.malloc_ext_P();
    fit_info.ext_P = (double**)malloc(sizeof(double*) * 2);
    fit_info.ext_P[0] = M_PS;
    fit_info.corr_id = { 1 };

    struct fit_result me = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_me, "me_oPp", fit_info,
        jack_file);
    check_correlatro_counter(2);


    fit_info.ext_P[1] = me.P[0];
    fit_info.corr_id = { 0 ,1 };
    struct fit_result fpi = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_fpi_P5A0, "fpi_P5A0", fit_info,
        jack_file);
    check_correlatro_counter(3);


    fit_info.ext_P[0] = M_PS_mu;
    fit_info.corr_id = { 3 };

    struct fit_result me_mu = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_me, "me_oPp_mu", fit_info,
        jack_file);
    check_correlatro_counter(4);

    fit_info.ext_P[1] = me_mu.P[0];
    fit_info.corr_id = { 2 ,3 };
    struct fit_result fpi_mu = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_fpi_P5A0, "fpi_P5A0_mu", fit_info,
        jack_file);
    check_correlatro_counter(5);

    ///////// deriv

    double* deriv = myres->create_copy(M_PS);
    for (int j = 0; j < Njack;j++) {
        deriv[j] = Z[j] * (fpi.P[0][j] - fpi_mu.P[0][j]) / dmu;
    }
    double *Z_fpi = myres->create_zero();
    double *Z_fpi_mu = myres->create_zero();
    for (int j = 0; j < Njack;j++) {
        Z_fpi[j] = fpi.P[0][j]*Z[j];
        Z_fpi_mu[j] = fpi_mu.P[0][j]*Z[j];
    }
    printf("Z_fpi: %.12g   %.12g\n", Z_fpi[Njack - 1], myres->comp_error(Z_fpi));
    printf("Z_fpi_mu: %.12g   %.12g\n", Z_fpi_mu[Njack - 1], myres->comp_error(Z_fpi_mu));
    printf("deriv: %.12g   %.12g\n", deriv[Njack - 1], myres->comp_error(deriv));
    std::string name_jack_fpi = "deriv/deriv_val_fpi_P5A0_" + std::string(argv[3]) + ".jack_txt";
    myres->write_jack_in_file(deriv, name_jack_fpi.c_str());


    name_jack_fpi = "deriv/fpi_P5A0_" + std::string(argv[3]) + ".jack_txt";
    myres->write_jack_in_file(Z_fpi, name_jack_fpi.c_str());

    write_jack(deriv, Njack, jack_file);
    check_correlatro_counter(6);


    /// fpi
    // fit_info.restore_default();
    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 3;
    fit_info.ext_P = malloc_2<double>(fit_info.n_ext_P, fit_info.Njack);
    for (int j = 0; j < fit_info.Njack; j++) {
        fit_info.ext_P[0][j] = M_PS[j];
        fit_info.ext_P[1][j] = head_P5.mus[0];
        fit_info.ext_P[2][j] = head_P5.mus[0];
    }
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = head_P5.T;
    fit_info.corr_id = { 1 };

    struct fit_result f_PS = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_function_f_PS, "f_{PS}", fit_info, jack_file);
    check_correlatro_counter(7);
    // with reweighting
    fit_info.corr_id = { 3 };
    for (int j = 0; j < fit_info.Njack; j++) {
        fit_info.ext_P[0][j] = M_PS_mu[j];
        fit_info.ext_P[1][j] = head_P5_mu.mus[1];
        fit_info.ext_P[2][j] = head_P5_mu.mus[1];
    }
    mysprintf(namefile, NAMESIZE, "f_{PS}_mu");
    struct fit_result f_PS_rew = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_function_f_PS, namefile, fit_info, jack_file);
    check_correlatro_counter(8);

    // free_fit_result(fit_info, fit_out);

    double* deriv_WTI = myres->create_copy(M_PS);
    for (int j = 0; j < Njack;j++) {
        deriv_WTI[j] = (f_PS.P[0][j] - f_PS_rew.P[0][j]) / dmu;
    }
    printf("deriv_WTI: %.12g   %.12g\n", deriv_WTI[Njack - 1], myres->comp_error(deriv_WTI));
    write_jack(deriv_WTI, Njack, jack_file);
    check_correlatro_counter(9);

    name_jack_fpi = "deriv/deriv_val_fpi_WTI_" + std::string(argv[3]) + ".jack_txt";
    myres->write_jack_in_file(deriv_WTI, name_jack_fpi.c_str());



    // free_fit_result(fit_info, fit_out);


    double* M_PS_A0 = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, 0, "M_{PS}_A0", M_eff_sinh_T, jack_file, fit_info);
    // free(M_PS);
    check_correlatro_counter(10);

    fit_info.restore_default();
}