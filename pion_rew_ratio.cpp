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
    fi = fread(&id, sizeof(int), 1, stream);
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
    // double it = x[0];
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
    error(head.Njack != head_rew.Njack, 1, "main", "Najck w0 = %d   while Njack rew1  = %d", head.Njack, head_rew.Njack);
    double ****data_rew = calloc_corr(head_rew.Njack, head_rew.ncorr, head_rew.T);

    for (int iconf = 0; iconf < head_rew.Njack; iconf++)
    {
        read_twopt(infile_rew, data_rew[iconf], head_rew);
        error(head.smearing[iconf].compare(head_rew.smearing[iconf]) != 0, 2, "main",
              "configuration order differ at %d\n flow file conf: %s\n rew 1 conf: %s ", iconf,
              head.smearing[iconf].c_str(), head_rew.smearing[iconf].c_str());
    }

    //////////////////////////////////////////////////////////////
    // read the reweighting 2
    //////////////////////////////////////////////////////////////

    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[8]);
    FILE *infile_rew2 = open_file(namefile, "r");

    generic_header head_rew2;
    head_rew2.read_header_debug(infile_rew2);
    error(head.Njack != head_rew2.Njack, 1, "main", "Najck w0 = %d   while Njack rew2  = %d", head.Njack, head_rew2.Njack);
    double ****data_rew2 = calloc_corr(head_rew2.Njack, head_rew2.ncorr, head_rew2.T);

    for (int iconf = 0; iconf < head_rew2.Njack; iconf++)
    {
        read_twopt(infile_rew2, data_rew2[iconf], head_rew2);
        error(head.smearing[iconf].compare(head_rew2.smearing[iconf]) != 0, 2, "main",
              "configuration order differ at %d\n flow file conf: %s\n rew 2 conf: %s ", iconf,
              head.smearing[iconf].c_str(), head_rew2.smearing[iconf].c_str());
    }

    //////////////////////////////////////////////////////////////
    // reading flow
    //////////////////////////////////////////////////////////////
    int ncorr_new = head.ncorr;        // current number of correlators
    int Max_corr = head.ncorr * 3 + 2; // max number of correlators

    double ****data = calloc_corr(head.Njack, Max_corr, head.T);

    printf("confs=%d\n", head.Njack);
    printf("ncorr=%d\n", head.ncorr);
    printf("kappa=%g\n", head.kappa);
    for (int iconf = 0; iconf < head.Njack; iconf++)
    {
        read_twopt(infile, data[iconf], head);
    }
    symmetrise_corr(head.Njack, 0, head.T, data);
    //////////////////////////////////////////////////////////////
    // normalization of the correlator
    //////////////////////////////////////////////////////////////
    for (int i = 0; i < head.ncorr; i++)
    {
        for (int j = 0; j < head.Njack; j++)
        {

            for (int tf = 0; tf < head.T; tf++)
            {
                data[j][i][tf][0] = data[j][i][tf][0] * head.kappa * head.kappa * 2.0;
            }
        }
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
    for (int ir = 0; ir < 2; ir++)
    {
        for (int i = 0; i < head.ncorr; i++)
        {
            for (int j = 0; j < head.Njack; j++)
            {

                double r;
                if (ir == 0)
                    r = data_rew[j][0][0][1]; // the im part is the exponentiated subtracted reweighting factor
                else if (ir == 1)
                    r = data_rew2[j][0][0][1]; // the im part is the exponentiated subtracted reweighting factor
                for (int tf = 0; tf < head.T; tf++)
                {
                    data[j][head.ncorr + i + ir * head.ncorr][tf][0] = data[j][i][tf][0] * r;
                    data[j][head.ncorr + i + ir * head.ncorr][tf][1] = data[j][i][tf][1] * r;
                }
            }
        }
    }
    for (int j = 0; j < head.Njack; j++)
    {
        for (int tf = 0; tf < head.T; tf++)
        {
            data[j][head.ncorr * 3][tf][0] = data_rew[j][0][0][1];
            data[j][head.ncorr * 3][tf][1] = 0;
            data[j][head.ncorr * 3 + 1][tf][0] = data_rew2[j][0][0][1];
            data[j][head.ncorr * 3 + 1][tf][1] = 0;
        }
    }
    // printf("corr at time 10 jacks\n");
    // for (int j = 0; j < head.Njack; j++)
    // {
    //     printf("%d %.12g  %.12g %.12g\n", j, data[j][head.ncorr + 0][10][0], data[j][head.ncorr * 2][0][0], data[j][0][10][0]);
    // }
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
    mysprintf(namefile, NAMESIZE, "%s/out/ratio_%s_%s_output", option[3], option[6], argv[8]);
    printf("writing output in :\n %s \n", namefile);
    FILE *outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/ratio_%s_%s_%s", option[3], option[4], option[6], argv[8]);
    FILE *jack_file = open_file(namefile, "w+");
    // write_header_g2(jack_file, head);
    head.write_header(jack_file);

    //////////////////////////////////////////////////////////////
    // making custom jackknifes
    //////////////////////////////////////////////////////////////
    double ****data_bin = bin_intoN(data, Max_corr, head.T, confs, bin); // binning into N=bin with not integer
    double ****conf_jack = myres->create(Neff, Max_corr, head.T, data_bin);
    free_corr(Neff, Max_corr, head.T, data_bin);

    free_corr(confs, Max_corr, head.T, data);

    /// normalize rew factor
    for (int ir = 0; ir < 2; ir++)
    {
        for (int i = 0; i < head.ncorr; i++)
        {
            for (int j = 0; j < head.Njack; j++)
            {

                double r = conf_jack[j][head.ncorr * 3 + ir][0][0]; // the im part is the exponentiated subtracted reweighting factor
                int id = head.ncorr + i + ir * head.ncorr;
                for (int tf = 0; tf < head.T; tf++)
                {
                    conf_jack[j][id][tf][0] = conf_jack[j][id][tf][0] / r;
                    conf_jack[j][id][tf][1] = conf_jack[j][id][tf][1] / r;
                }
            }
        }
    }
    // printf("corr at time 10 jacks\n");
    // for (int j = 0; j < head.Njack; j++){
    //     printf("%d %.12g  %.12g %.12g\n",j, conf_jack[j][head.ncorr +0][10][0], conf_jack[j][head.ncorr *2][0][0],  conf_jack[j][0][10][0]);
    // }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print all the effective masses correlators
    // set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/ratio_%s_%s_meff_correlators", option[3], option[6], argv[7]);
    FILE *outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/ratio_%s_%s_raw_correlators", option[3], option[6], argv[7]);
    FILE *outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/ratio_%s_%s_shifted_correlators", option[3], option[6], argv[7]);
    FILE *outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/ratio_%s_%s_log_meff_shifted", option[3], option[6], argv[7]);
    FILE *outfile_log_meff_shifted = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/ratio_%s_%s_gamma", option[3], option[6], argv[7]);
    FILE *out_gamma = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/ratio_%s_%s_HLT_kernel", option[3], option[6], argv[7]);
    FILE *outfile_HLT_kernel = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/ratio_%s_%s_HLT_AoverB", option[3], option[6], argv[7]);
    FILE *outfile_HLT_AoverB = open_file(namefile, "w+");

    char save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");
    FILE *dev_null = open_file("/dev/null", "w");
    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;

    for (int icorr = 0; icorr < Max_corr; icorr++)
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

    double dmul = head_rew.oranges[0] - head_rew.mus[0];

    double dmu2 = head_rew2.oranges[0] - head_rew2.mus[0];

    //////////////////////////////////////////////////////////////
    //  read m^iso
    //////////////////////////////////////////////////////////////
    std::vector<double *> amuiso(3);
    double mean, err;
    int seed;
    line_read_param(option, "muliso", mean, err, seed, namefile_plateaux);
    amuiso[0] = myres->create_fake(mean, err, seed);
    line_read_param(option, "musiso", mean, err, seed, namefile_plateaux);
    amuiso[1] = myres->create_fake(mean, err, seed);
    line_read_param(option, "muciso", mean, err, seed, namefile_plateaux);
    amuiso[2] = myres->create_fake(mean, err, seed);

    std::vector<double *> amusim(3);
    line_read_param(option, "mulsim", mean, err, seed, namefile_plateaux);
    amusim[0] = myres->create_fake(mean, err, seed);
    line_read_param(option, "mussim", mean, err, seed, namefile_plateaux);
    amusim[1] = myres->create_fake(mean, err, seed);
    line_read_param(option, "mucsim", mean, err, seed, namefile_plateaux);
    amusim[2] = myres->create_fake(mean, err, seed);

    line_read_param(option, "a", mean, err, seed, namefile_plateaux);
    double *a_fm = myres->create_fake(mean, err, seed);

    ////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////
    corr_counter = -1;

    double *M_PS = plateau_correlator_function(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, 0, "M_{PS}", M_eff_T, jack_file);

    check_correlatro_counter(0);
    char name_rew[NAMESIZE];

    mysprintf(name_rew, NAMESIZE, "M_{PS}_%s", argv[7]);
    double *M_PS_rewl = plateau_correlator_function(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, head.ncorr + 0, name_rew, M_eff_T, jack_file);

    mysprintf(name_rew, NAMESIZE, "M_{PS}_%s", argv[8]);
    double *M_PS_rew2 = plateau_correlator_function(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, head.ncorr * 2 + 0, name_rew, M_eff_T, jack_file);

    double *zero = (double *)calloc(Njack, sizeof(double));
    write_jack(zero, Njack, jack_file);
    write_jack(zero, Njack, jack_file);
    write_jack(zero, Njack, jack_file);
    write_jack(zero, Njack, jack_file);
    write_jack(zero, Njack, jack_file);
    write_jack(zero, Njack, jack_file);
    write_jack(zero, Njack, jack_file);
    write_jack(zero, Njack, jack_file);

    double *tmp = myres->create_fake(head_rew.mus[0], 1e-20, 1);
    print_result_in_file(outfile, tmp, "mu_in", 0, 0, 0);
    write_jack(tmp, Njack, jack_file);
    check_correlatro_counter(11);
    free(tmp);
    tmp = myres->create_fake(head_rew.oranges[0], 1e-20, 1);
    print_result_in_file(outfile, tmp, "mu_out", 0, 0, 0);
    write_jack(tmp, Njack, jack_file);
    check_correlatro_counter(12);
    free(tmp);

    write_jack(amuiso[0], Njack, jack_file);
    check_correlatro_counter(13);
    write_jack(amuiso[1], Njack, jack_file);
    check_correlatro_counter(14);
    write_jack(amuiso[2], Njack, jack_file);
    check_correlatro_counter(15);
    write_jack(a_fm, Njack, jack_file);
    check_correlatro_counter(16);

    write_jack(amusim[0], Njack, jack_file);
    check_correlatro_counter(17);
    write_jack(amusim[1], Njack, jack_file);
    check_correlatro_counter(18);
    write_jack(amusim[2], Njack, jack_file);
    check_correlatro_counter(19);

    //////////////////////////////////////////////////////////////
    // plateau derivative ratio
    //////////////////////////////////////////////////////////////

    // fit_info.codeplateaux = true;
    // int tmin, tmax,sep;
    // mysprintf(name_rew, NAMESIZE, "M_{PS}_%s", argv[8]);
    // line_read_param(option, name_rew, fit_info.tmin, fit_info.tmax, fit_info.sep , namefile_plateaux);

    fit_type fit_info;
    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = head.T;
    fit_info.corr_id = {head.ncorr + 0, 0};
    fit_info.ave_P = {dmul};

    mysprintf(name_rew, NAMESIZE, "plateau_dM/d%s", argv[7]);
    struct fit_result fit_dM_dmul = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_plateau_dM_dmu, name_rew, fit_info,
        jack_file);
    check_correlatro_counter(20);

    fit_info.corr_id = {head.ncorr * 2 + 0, 0};
    fit_info.ave_P = {dmu2};

    mysprintf(name_rew, NAMESIZE, "plateau_dM/d%s", argv[8]);
    struct fit_result fit_dM_dmu2 = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_plateau_dM_dmu, name_rew, fit_info,
        jack_file);
    check_correlatro_counter(21);
    ///////

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = head.T;
    fit_info.corr_id = {head.ncorr * 2 + 0, head.ncorr + 0, 0};
    fit_info.ave_P = {dmu2, dmul};

    mysprintf(name_rew, NAMESIZE, "plateau_ratio_dM/d%s", argv[8]);
    struct fit_result fit_ratio_dM_dmu = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_plateau_ratio_dM_dmu, name_rew, fit_info,
        jack_file);
    check_correlatro_counter(22);
    // fit_info.restore_default();

    // fit_info.n_ext_P = 4;
    // fit_info.ext_P = malloc_2<double>(fit_info.n_ext_P, fit_info.Njack);
    // for (int j = 0; j < fit_info.Njack; j++)
    // {
    //     fit_info.ext_P[0][j] = M_PS_rew[j];
    //     fit_info.ext_P[1][j] = M_PS[j];
    //     fit_info.ext_P[2][j] = amusim[0][j];
    //     fit_info.ext_P[3][j] = amusim[0][j];
    // }

    // mysprintf(name_rew, NAMESIZE, "plateau_df/d%s", argv[8]);
    // struct fit_result fit_df_dmu = fit_fun_to_fun_of_corr(
    //     option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
    //     outfile, lhs_plateau_df_dmu, name_rew, fit_info,
    //     jack_file);
    // check_correlatro_counter(21);
}