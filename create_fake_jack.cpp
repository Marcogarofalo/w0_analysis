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

double lhs_mpcac(int j, double ****in, int t, struct fit_type fit_info)
{
    int id_V = fit_info.corr_id[0];
    int id_P = fit_info.corr_id[1];
    // we should take -Im of V0P5 which is saved as real part in this data
    double r = -(in[j][id_V][t + 1][0] - in[j][id_V][t][0]) / (2. * in[j][id_P][t][0]);

    return r;
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

    // FILE *infile = open_file(namefile, "r");

    //////////////////////////////////// read and setup header
    generic_header head;
    head.Njack = 100;
    char name[NAMESIZE];
    double mean, err;
    int seed;
    mysprintf(name, NAMESIZE, "mu_in_%s", argv[8]);
    line_read_param(option, name, mean, err, seed, namefile_plateaux);
    head.mus.push_back(mean);

    mysprintf(name, NAMESIZE, "mu_out_%s", argv[8]);
    line_read_param(option, name, mean, err, seed, namefile_plateaux);
    head.oranges.push_back(mean);

    // head.read_header_debug(infile);
    // head.print_header();
    init_global_head(head);

    //////////////////////////////////////////////////////////////
    // read the reweighting
    //////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////
    // reading flow
    //////////////////////////////////////////////////////////////
    int ncorr_new = head.ncorr;        // current number of correlators
    int Max_corr = head.ncorr * 2 + 1; // max number of correlators

    double ****data = calloc_corr(head.Njack, Max_corr, head.T);

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

    //////////////////////////////////////////////////////////////
    // making custom jackknifes
    //////////////////////////////////////////////////////////////
    double ****data_bin = bin_intoN(data, Max_corr, head.T, confs, bin); // binning into N=bin with not integer
    double ****conf_jack = myres->create(Neff, Max_corr, head.T, data_bin);
    free_corr(Neff, Max_corr, head.T, data_bin);

    free_corr(confs, Max_corr, head.T, data);

    // /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // // print all the effective masses correlators
    // // set the option to not read for a plateaux
    // mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_meff_correlators", option[3], option[6], argv[7]);
    // FILE *outfile_meff_corr = open_file(namefile, "w+");
    // mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_raw_correlators", option[3], option[6], argv[7]);
    // FILE *outfile_raw_corr = open_file(namefile, "w+");
    // mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_shifted_correlators", option[3], option[6], argv[7]);
    // FILE *outfile_shifted_corr = open_file(namefile, "w+");
    // mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_log_meff_shifted", option[3], option[6], argv[7]);
    // FILE *outfile_log_meff_shifted = open_file(namefile, "w+");
    // mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_gamma", option[3], option[6], argv[7]);
    // FILE *out_gamma = open_file(namefile, "w+");

    // mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_HLT_kernel", option[3], option[6], argv[7]);
    // FILE *outfile_HLT_kernel = open_file(namefile, "w+");
    // mysprintf(namefile, NAMESIZE, "%s/out/%s_%s_HLT_AoverB", option[3], option[6], argv[7]);
    // FILE *outfile_HLT_AoverB = open_file(namefile, "w+");

    // char save_option[NAMESIZE];
    // sprintf(save_option, "%s", option[1]);
    // sprintf(option[1], "blind");
    // FILE *dev_null = open_file("/dev/null", "w");
    // struct fit_type fit_info_silent;
    // fit_info_silent.verbosity = -1;
    // fit_info_silent.chi2_gap_jackboot = 1e+6;
    // fit_info_silent.guess_per_jack = 0;

    // for (int icorr = 0; icorr < Max_corr; icorr++)
    // {
    //     // log effective mass
    //     double *tmp_meff_corr = plateau_correlator_function(
    //         option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
    //         namefile_plateaux, outfile_meff_corr, icorr, "log", M_eff_log, dev_null,
    //         fit_info_silent);
    //     free(tmp_meff_corr);
    //     // raw correlator
    //     file_head.l0 = head.T * 2;
    //     tmp_meff_corr = plateau_correlator_function(
    //         option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
    //         namefile_plateaux, outfile_raw_corr, icorr, "cor", identity, dev_null,
    //         fit_info_silent);
    //     free(tmp_meff_corr);
    //     tmp_meff_corr = plateau_correlator_function(
    //         option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
    //         namefile_plateaux, outfile_raw_corr, icorr, "cor", identity_im,
    //         dev_null, fit_info_silent);
    //     free(tmp_meff_corr);
    //     file_head.l0 = head.T;
    //     // shifted correlator
    //     tmp_meff_corr = plateau_correlator_function(
    //         option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
    //         namefile_plateaux, outfile_shifted_corr, icorr, "shift_cor", shift_corr,
    //         dev_null, fit_info_silent);
    //     free(tmp_meff_corr);
    //     // log_meff shifted correlator
    //     tmp_meff_corr = plateau_correlator_function(
    //         option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
    //         namefile_plateaux, outfile_log_meff_shifted, icorr, "log_shift",
    //         M_eff_log_shift, dev_null, fit_info_silent);
    //     free(tmp_meff_corr);
    // }
    // fit_info_silent.restore_default();
    // sprintf(option[1], "%s", save_option); // restore option
    corr_counter = -1;

    ////////////////////////////////////////////////////////////
    // symmetrize
    ////////////////////////////////////////////////////////////
    // symmetrise_jackboot(Njack, 0, head.T, conf_jack);
    // symmetrise_jackboot(Njack, 1, head.T, conf_jack, -1);

    //////////////////////////////////////////////////////////////
    //  read m^iso
    //////////////////////////////////////////////////////////////
    std::vector<double *> amuiso(3);
    
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

    double *zeros = (double *)calloc(Njack, sizeof(double));
    mysprintf(name, NAMESIZE, "M_{PS}", argv[8]);
    line_read_param(option, name, mean, err, seed, namefile_plateaux);
    double *M_PS = myres->create_fake(mean, err, seed);
    write_jack(M_PS, Njack, jack_file);// M_PS

    // free(M_PS);
    check_correlatro_counter(0);

    //////////////// mpcac

    write_jack(zeros, Njack, jack_file);

    check_correlatro_counter(1);

    write_jack(zeros, Njack, jack_file);

    check_correlatro_counter(2);
    write_jack(zeros, Njack, jack_file);

    check_correlatro_counter(3);

    //////////////////////////////////////////////////////////////
    // detiv MPS
    //////////////////////////////////////////////////////////////
    mysprintf(name, NAMESIZE, "dM_{PS}/dmu%s", argv[8]);
    line_read_param(option, name, mean, err, seed, namefile_plateaux);
    double *derM = myres->create_fake(mean, err, seed);
    write_jack(derM, Njack, jack_file);

    check_correlatro_counter(4);

    /// correlator at time 20
    write_jack(zeros, Njack, jack_file);

    check_correlatro_counter(5);

    //////////////////////////////////////////////////////////////
    // fpi
    //////////////////////////////////////////////////////////////
    mysprintf(name, NAMESIZE, "f_{PS}", argv[8]);
    line_read_param(option, name, mean, err, seed, namefile_plateaux);
    double *f_PS = myres->create_fake(mean, err, seed);
    write_jack(f_PS, Njack, jack_file);

    check_correlatro_counter(6);
    // with reweighting
    write_jack(zeros, Njack, jack_file);

    check_correlatro_counter(7);
    // der mpi_fpi /dmu
    write_jack(zeros, Njack, jack_file);

    check_correlatro_counter(8);
    // der mpi_fpi /dmu
    mysprintf(name, NAMESIZE, "df_{PS}/dmu%s", argv[8]);
    line_read_param(option, name, mean, err, seed, namefile_plateaux);
    double *derf = myres->create_fake(mean, err, seed);
    write_jack(derf, Njack, jack_file);

    check_correlatro_counter(9);
    write_jack(zeros, Njack, jack_file);

    check_correlatro_counter(10);

    ///

    double *tmp = myres->create_fake(head.mus[0], 1e-20, 1);
    print_result_in_file(outfile, tmp, "mu_in", 0, 0, 0);
    write_jack(tmp, Njack, jack_file);
    check_correlatro_counter(11);
    free(tmp);
    tmp = myres->create_fake(head.oranges[0], 1e-20, 1);
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
    
}