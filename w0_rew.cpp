#define CONTROL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #include <time.h>
// #include <complex.h>
#include <cstring>
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
struct kinematic kinematic_2pt;

generic_header read_head(FILE *stream)
{
    generic_header header;
    return header;
}
void write_header_g2(FILE *jack_file, generic_header head)
{
    int fi = 0;
    fi += fwrite(&head.T, sizeof(int), 1, jack_file);
    fi += fwrite(&head.L, sizeof(int), 1, jack_file);
    int nmu = head.mus.size();
    fi += fwrite(&nmu, sizeof(int), 1, jack_file);
    for (double mu : head.mus)
    {
        fi += fwrite(&mu, sizeof(double), 1, jack_file);
    }
}

char **argv_to_options(char **argv)
{
    char **option;
    option = (char **)malloc(sizeof(char *) * 7);
    option[0] = (char *)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char *)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char *)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char *)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char *)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char *)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char *)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, "read_plateaux"); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p");            // -p
    mysprintf(option[3], NAMESIZE, argv[2]);         // path
    mysprintf(option[4], NAMESIZE, argv[6]);         // resampling
    mysprintf(option[5], NAMESIZE, "no");            // pdf
    mysprintf(option[6], NAMESIZE, argv[3]);         // infile
    return option;
}

void init_global_head(generic_header head)
{
    file_head.l1 = head.L;
    file_head.l0 = head.T;
    file_head.l2 = head.L;
    file_head.l3 = head.L;
    file_head.nk = 2;
    file_head.musea = head.mus[0];
    file_head.k = (double *)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;
    file_head.k[1] = 0;
    file_head.k[2] = head.mus[0];
    file_head.k[3] = head.mus[0];

    file_head.nmoms = 1;
    file_head.mom = (double **)malloc(sizeof(double *) * file_head.nmoms);
    for (int i = 0; i < file_head.nmoms; i++)
    {
        file_head.mom[i] = (double *)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }
}

void read_twopt(FILE *stream, double ***to_write, generic_header head)
{
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
    for (int k = 0; k < head.ncorr; k++)
    {
        for (int t = 0; t < head.T; t++)
        {
            fi += fread(to_write[k][t], sizeof(double), 2, stream);
        }
    }
}

double int2flowt(double i)
{
    return 0.010000 + i * 0.02;
}

double poly3(int n, int Nvar, double *x, int Npar, double *P)
{
    double it = x[0];
    double tf = int2flowt(x[0]);
    return P[0] + P[1] * tf + P[2] * tf * tf + P[3] * tf * tf * tf;
}

int main(int argc, char **argv)
{
    error(argc != 9, 1, "main ",
          "usage:././w0_rew -p path file -bin $bin  jack/boot   reweighting_factors  name_rew\n separate "
          "path and file please");

    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, argv[6]);
    printf("resampling: %s\n", resampling);

    char **option = argv_to_options(argv);

    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], option[6]);

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");

    FILE *infile = open_file(namefile, "r");

    //////////////////////////////////// read and setup header
    generic_header head;
    head.read_header_debug(infile);
    // head.print_header();
    init_global_head(head);

    //////////////////////////////////////////////////////////////
    // read the reweighting
    //////////////////////////////////////////////////////////////

    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[7]);
    FILE *infile_rew = open_file(namefile, "r");

    generic_header head_rew;
    head_rew.read_header_debug(infile_rew);
    error(head.Njack != head_rew.Njack, 1, "main", "Najck w0 = %d   while Njack rew  = %d", head.Njack, head_rew.Njack);
    double ****data_rew = calloc_corr(head_rew.Njack, head_rew.ncorr, head_rew.T);
    for (int iconf = 0; iconf < head_rew.Njack; iconf++)
    {
        read_twopt(infile_rew, data_rew[iconf], head_rew);
        error(head.smearing[iconf].compare(head_rew.smearing[iconf]) != 0, 2, "main",
              "configuration order differ at %d\n flow file conf: %s\n loops conf: %s ", iconf,
              head.smearing[iconf].c_str(), head_rew.smearing[iconf].c_str());
    }

    //////////////////////////////////////////////////////////////
    // reading flow
    //////////////////////////////////////////////////////////////
    int ncorr_new = head.ncorr;    // current number of correlators
    int Max_corr = head.ncorr + 2; // max number of correlators

    double ****data = calloc_corr(head.Njack, Max_corr, head.T);

    printf("confs=%d\n", head.Njack);
    printf("ncorr=%d\n", head.ncorr);
    printf("kappa=%g\n", head.kappa);
    for (int iconf = 0; iconf < head.Njack; iconf++)
    {
        read_twopt(infile, data[iconf], head);
    }
    //////////////////////////////////////////////////////////////
    // removing confs
    //////////////////////////////////////////////////////////////
    // std::vector<std::string> to_remove = {"1776_r0", "1540_r0", "1088_r1"};
    // double ****data_removed = calloc_corr(head.Njack - to_remove.size(), Max_corr, head.T);
    // double ****rew_removed = calloc_corr(head.Njack - to_remove.size(), head_rew.ncorr, head_rew.T);
    // int countj = 0;
    // bool get_out = false;

    // for (int j = 0; j < head.Njack; j++)
    // {
    //     for (std::string conf : to_remove)
    //     {
    //         if (strcmp(head.smearing[j].c_str(), conf.c_str()) == 0)
    //         {
    //             get_out = true;
    //         }
    //     }
    //     if (get_out)
    //     {
    //         printf("removing %s \n", head.smearing[j].c_str());
    //         get_out = false;
    //         continue;
    //     }

    //     // if we arrive here it means that we keep the configuration
    //     for (int k = 0; k < Max_corr; k++)
    //     {
    //         for (int t = 0; t < head.T; t++)
    //         {
    //             data_removed[countj][k][t][0] = data[j][k][t][0];
    //             data_removed[countj][k][t][1] = data[j][k][t][1];
    //         }
    //     }
    //     for (int k = 0; k < head_rew.ncorr; k++)
    //     {
    //         for (int t = 0; t < head_rew.T; t++)
    //         {
    //             rew_removed[countj][k][t][0] = data_rew[j][k][t][0];
    //             rew_removed[countj][k][t][1] = data_rew[j][k][t][1];
    //         }
    //     }
    //     countj++;
    // }
    // free_corr(head.Njack, Max_corr, head.T, data);
    // free_corr(head.Njack, head_rew.ncorr, head_rew.T, data_rew);
    // head.Njack = head.Njack - to_remove.size();
    // head_rew.Njack = head_rew.Njack - to_remove.size(); // the reweighting has the same number of confs
    // data = data_removed;
    // data_rew = rew_removed;
    // printf("removed %ld confs\n", to_remove.size());

    //////////////////////////////////////////////////////////////
    // correcting w0
    //////////////////////////////////////////////////////////////
    // for (int i = 0; i < head_rew.ncorr; i++)
    // {
    double sum_r = 0;
    for (int j = 0; j < head.Njack; j++)
    {
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
        for (int tf = 0; tf < head.T; tf++)
        {
            data[j][head.ncorr + 0][tf][0] = data[j][6][tf][0] * r;
            data[j][head.ncorr + 0][tf][1] = data[j][6][tf][1] * r;
            data[j][head.ncorr + 1][tf][0] = data_rew[j][0][0][1];
            data[j][head.ncorr + 1][tf][1] = 0;
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
    if (strcmp(argv[6], "jack") == 0)
    {
        Njack = Neff + 1;
        myres = new resampling_jack(Neff);
    }
    else if (strcmp(argv[6], "boot") == 0)
    {
        Njack = (Neff * 2 + 1);
        myres = new resampling_boot(Neff * 2);
    }
    else
    {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[7]= %s is not jack or boot", argv[7]);
    }
    // now Njack need to be the number of jacks
    head.Nconf = head.Njack;
    head.Njack = Njack;
    //////////////////////////////////// setup output files
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_output", option[3], option[6], argv[7]);
    printf("writing output in :\n %s \n", namefile);
    FILE *outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s_%s", option[3], option[4], option[6], argv[7]);
    FILE *jack_file = open_file(namefile, "w+");
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
    double ****data_bin = bin_intoN(data, Max_corr, head.T, confs, bin); // binning into N=bin with not integer
    double ****conf_jack = myres->create(Neff, Max_corr, head.T, data_bin);
    free_corr(Neff, Max_corr, head.T, data_bin);

    double ****conf_jack_Orew = calloc_corr(head.Njack, head_rew.ncorr, head.T);

    double ****conf_jack_r = bin_intoN_exp(data_rew, 0, head_rew.T, confs, bin);
    double ****conf_jack_rO = bin_intoN_exp1(data_rew, data, 0, 6, head.T, confs, bin);

    free_corr(head.Njack, head_rew.ncorr, head_rew.T, data_rew);

    printf("Njack = %d  Nbins = %d  confs = %d\n", head.Njack, bin, confs);
    // for (int i = 0; i < head_rew.ncorr; i++)
    // {
    //     // double sum_r = 0;
    //     for (int j = 0; j < head.Njack; j++)
    //     {
    //         for (int tf = 0; tf < head.T; tf++)
    //         {
    //             for (int j1 = 0; j1 < bin; j1++)
    //             {
    //                 if (j1 == j)
    //                     continue;
    //                 double r = 0;
    //                 for (int j2 = 0; j2 < bin; j2++)
    //                 {
    //                     if (j2 == j)
    //                         continue;
    //                     else
    //                     {
    //                         r += exp(conf_jack_r[j2][0][0][0] - conf_jack_rO[j1][0][tf][0]);
    //                         // printf("r=%g   %g    %g \n", r, conf_jack_r[j2][0][0][0], conf_jack_rO[j1][0][tf][0]);
    //                     }
    //                 }
    //                 conf_jack_Orew[j][i][tf][0] += 1.0 / r;
    //             }
    //         }
    //     }
    // }
    make_ratio_of_jacks(conf_jack_Orew, head.Njack, 0, head.T, conf_jack_rO, 0, conf_jack_r, 0);
    free_corr(bin, head_rew.ncorr, head_rew.T, conf_jack_r);
    free_corr(bin, head_rew.ncorr, head.T, conf_jack_rO);

    // double ****conf_jack_Orew_s;
    // if (argc > 8)
    // {
    //     // opening file
    //     mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[8]);
    //     FILE *infile_rew = open_file(namefile, "r");

    //     // reading the reweighting
    //     generic_header head_rew_1;
    //     head_rew_1.read_header_debug(infile_rew);
    //     error(head.Njack != head_rew_1.Njack, 1, "main", "Najck w0 = %d   while Njack rew  = %d", head.Njack, head_rew_1.Njack);
    //     double ****data_rew = calloc_corr(head_rew_1.Njack, head_rew_1.ncorr, head_rew_1.T);
    //     for (int iconf = 0; iconf < head_rew_1.Njack; iconf++)
    //     {
    //         read_twopt(infile_rew, data_rew[iconf], head_rew_1);
    //         error(head.smearing[iconf].compare(head_rew_1.smearing[iconf]) != 0, 2, "main",
    //               "configuration order differ at %d\n flow file conf: %s\n loops conf: %s ", iconf,
    //               head.smearing[iconf].c_str(), head_rew_1.smearing[iconf].c_str());
    //     }

    //     double ****conf_jack_r = bin_intoN_exp(data_rew, 0, head_rew_1.T, confs, bin);
    //     double ****conf_jack_rO = bin_intoN_exp1(data_rew, data, 0, 6, head.T, confs, bin);

    //     conf_jack_Orew_s= calloc_corr(head.Njack, head_rew_1.ncorr, head.T);

    //     make_ratio_of_jacks(conf_jack_Orew, head.Njack, 0, head.T, conf_jack_rO, 0, conf_jack_r, 0);
    //     free_corr(bin, head_rew.ncorr, head_rew.T, conf_jack_r);
    //     free_corr(bin, head_rew.ncorr, head.T, conf_jack_rO);

    // }

    free_corr(head.Njack, 6, head.T, data);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print all the effective masses correlators
    // set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_meff_correlators", option[3], option[6], argv[7]);
    FILE *outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_raw_correlators", option[3], option[6], argv[7]);
    FILE *outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_shifted_correlators", option[3], option[6], argv[7]);
    FILE *outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_log_meff_shifted", option[3], option[6], argv[7]);
    FILE *outfile_log_meff_shifted = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_gamma", option[3], option[6], argv[7]);
    FILE *out_gamma = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_HLT_kernel", option[3], option[6], argv[7]);
    FILE *outfile_HLT_kernel = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_HLT_AoverB", option[3], option[6], argv[7]);
    FILE *outfile_HLT_AoverB = open_file(namefile, "w+");

    char save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");
    FILE *dev_null = open_file("/dev/null", "w");
    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;

    for (int icorr = 0; icorr < head.ncorr; icorr++)
    {
        // log effective mass
        double *tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_meff_corr, icorr, "log", M_eff_log, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        // raw correlator
        file_head.l0 = head.T * 2;
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity_im,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        file_head.l0 = head.T;
        // shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_shifted_corr, icorr, "shift_cor", shift_corr,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
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
    // std::vector<double> w0(Njack);
    // for (size_t j = 0; j < Njack; j++) {
    //     w0[j] = rtbis_func_eq_input(fit_info.function, l /*n*/, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[j], 0, Mpi2_fpi2_phys[j], 1e-6, 2, 1e-10, 2);
    // }

    // double* M_PS = plateau_correlator_function(
    //     option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
    //     namefile_plateaux, outfile, 0, "M_{PS}", M_eff_T, jack_file);
    // free(M_PS);
    // check_correlatro_counter(0);

    // eg of fit to correlator
    struct fit_type fit_info;
    struct fit_result fit_out;

    fit_info.Nvar = 1;
    fit_info.Npar = 4;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.codeplateaux = true;
    // fit_info.n_ext_P = head.T;
    // fit_info.ext_P = (double **)malloc(sizeof(double *) * fit_info.n_ext_P);
    // fit_info.ext_P[0] = something;
    // fit_info.ave_P = std::vector<double>(head.T);
    // for (int t = 0; t < head.T; t++)
    // {
    //     fit_info.ext_P[t] = (double *)malloc(sizeof(double) * Njack - 1);
    //     for (int j = 0; j < Njack; j++)
    //     {
    //         fit_info.ext_P[t][j] = conf_jack[Njack - 1][0][t][0];
    //     }
    // }
    fit_info.function = poly3;
    fit_info.linear_fit = true;
    fit_info.T = head.T;
    fit_info.corr_id = {6};

    for (int t = 1; t < head.T; t++)
    {
        if (lhs_function_w0_eg(Njack - 1, conf_jack, t, fit_info) > 0.3)
        {
            fit_info.tmin = t - 2;
            fit_info.tmax = t + 1;
            break;
        }
    }

    // // print for frezzotti
    // {
    //     double *tmp = (double *)malloc(sizeof(double) * Njack);
    //     for (int j = 0; j < Njack; j++)
    //     {
    //         tmp[j] = lhs_function_w0_eg(j, conf_jack, 150, fit_info);
    //     }
    //     char name[NAMESIZE];
    //     mysprintf(name, NAMESIZE, "%s/out/W_t150_jack.txt", option[3]);
    //     myres->write_jack_in_file(tmp, name);
    //     free(tmp);
    // }
    // c++ 0 || r 1
    struct fit_result fit_W = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_function_w0_eg, "W(t)", fit_info,
        jack_file);
    check_correlatro_counter(0);

    double **tif = swap_indices(fit_info.Npar, Njack, fit_W.P);
    std::vector<double> swapped_x(fit_info.Nvar);
    std::vector<double> w0(Njack);

    for (size_t j = 0; j < Njack; j++)
    {
        w0[j] = rtbis_func_eq_input(fit_info.function, 0 /*n*/, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[j], 0, 0.3, fit_info.tmin, fit_info.tmax, 1e-10, 2);
        w0[j] = int2flowt(w0[j]);
        w0[j] = std::sqrt(w0[j]);
    }
    printf("w0 = %g  %g\n", w0[Njack - 1], myres->comp_error(w0.data()));

    print_result_in_file(outfile, w0.data(), "w0", 0.0, fit_info.tmin, fit_info.tmax);
    free_2(Njack, tif);
    //////////////////////////////////////////////////////////////
    // charm reweighting
    //////////////////////////////////////////////////////////////
    {
        fit_info.corr_id = {head.ncorr + 0, head.ncorr + 1};
        // fit_info.linear_fit = false;

        for (int t = 1; t < head.T; t++)
        {
            if (lhs_function_W_rew(Njack - 1, conf_jack, t, fit_info) > 0.3)
            {
                fit_info.tmin = t - 5;
                fit_info.tmax = t + 5;
                break;
            }
            if (t == head.T - 1)
                printf("Warning: W(t) with charm reweighting never gets to 0.3\n");
        }

        // // print for frezzotti
        // {
        //     double *tmp = (double *)malloc(sizeof(double) * Njack);
        //     for (int j = 0; j < Njack; j++)
        //     {
        //         tmp[j] = lhs_function_W_rew(j, conf_jack, 150, fit_info);
        //     }
        //     char name[NAMESIZE];
        //     mysprintf(name, NAMESIZE, "%s/out/W_rew_charm_OS_t150_jack.txt", option[3]);
        //     myres->write_jack_in_file(tmp, name);
        //     free(tmp);
        // }
        char name_rew[NAMESIZE];
        mysprintf(name_rew, NAMESIZE, "W_%s(t)", argv[8]);
        struct fit_result fit_W = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_function_W_rew, name_rew, fit_info,
            jack_file);
        check_correlatro_counter(1);

        double **tif = swap_indices(fit_info.Npar, Njack, fit_W.P);
        std::vector<double> swapped_x(fit_info.Nvar);
        std::vector<double> w0(Njack);

        mysprintf(name_rew, NAMESIZE, "w0_%s", argv[8]);
        for (size_t j = 0; j < Njack; j++)
        {
            w0[j] = rtbis_func_eq_input(fit_info.function, 0 /*n*/, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[j], 0, 0.3, fit_info.tmin - 5, fit_info.tmax + 5, 1e-10, 2);
            w0[j] = int2flowt(w0[j]);
            w0[j] = std::sqrt(w0[j]);
        }
        printf("%s = %g  %g\n", name_rew, w0[Njack - 1], myres->comp_error(w0.data()));

        print_result_in_file(outfile, w0.data(), name_rew, 0.0, fit_info.tmin, fit_info.tmax);
        free_2(Njack, tif);
    }

    // //////////////////////////////////////////////////////////////
    // // charm reweighting
    // //////////////////////////////////////////////////////////////
    // {
    //     fit_info.corr_id = {0};
    //     fit_info.linear_fit = true;

    //     for (int t = 1; t < head.T; t++)
    //     {
    //         if (lhs_function_w0_eg(Njack - 1, conf_jack_Orew, t, fit_info) > 0.3)
    //         {
    //             fit_info.tmin = t - 5;
    //             fit_info.tmax = t + 5;
    //             break;
    //         }
    //         if (t == head.T - 1)
    //             printf("Warning: W(t) with charm reweighting never gets to 0.3\n");
    //     }

    //     // // print for frezzotti
    //     // {
    //     //     double *tmp = (double *)malloc(sizeof(double) * Njack);
    //     //     for (int j = 0; j < Njack; j++)
    //     //     {
    //     //         tmp[j] = lhs_function_w0_eg(j, conf_jack_Orew, 150, fit_info);
    //     //     }
    //     //     char name[NAMESIZE];
    //     //     mysprintf(name, NAMESIZE, "%s/out/W_rew_charm_OS_t150_jack.txt", option[3]);
    //     //     myres->write_jack_in_file(tmp, name);
    //     //     free(tmp);
    //     // }
    //     char name_rew[NAMESIZE];
    //     mysprintf(name_rew, NAMESIZE, "W_%s(t)", argv[8]);
    //     struct fit_result fit_W = fit_fun_to_fun_of_corr(
    //         option, kinematic_2pt, (char *)"P5P5", conf_jack_Orew, namefile_plateaux,
    //         outfile, lhs_function_w0_eg, name_rew, fit_info,
    //         jack_file);
    //     check_correlatro_counter(2);

    //     double **tif = swap_indices(fit_info.Npar, Njack, fit_W.P);
    //     std::vector<double> swapped_x(fit_info.Nvar);
    //     std::vector<double> w0(Njack);

    //     mysprintf(name_rew, NAMESIZE, "w0_%s", argv[8]);
    //     for (size_t j = 0; j < Njack; j++)
    //     {
    //         w0[j] = rtbis_func_eq_input(fit_info.function, 0 /*n*/, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[j], 0, 0.3, fit_info.tmin - 5, fit_info.tmax + 5, 1e-10, 2);
    //         w0[j] = int2flowt(w0[j]);
    //         w0[j] = std::sqrt(w0[j]);
    //     }
    //     printf("%s = %g  %g\n", name_rew, w0[Njack - 1], myres->comp_error(w0.data()));

    //     print_result_in_file(outfile, w0.data(), name_rew, 0.0, fit_info.tmin, fit_info.tmax);
    //     free_2(Njack, tif);
    // }

    fit_info.restore_default();
}