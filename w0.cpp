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
    error(argc != 8, 1, "nissa_mpcac ",
          "usage:./nissa_mpcac -p path file -bin $bin  jack/boot   loops\n separate "
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
    // read the loops
    //////////////////////////////////////////////////////////////

    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[7]);
    FILE *infile_loop = open_file(namefile, "r");

    generic_header head_loops;
    head_loops.read_header_debug(infile_loop);
    error(head.Njack != head_loops.Njack, 1, "main", "Najck w0 = %d   while Njack loops  = %d", head.Njack, head_loops.Njack);
    double ****data_loop = calloc_corr(head_loops.Njack, head_loops.ncorr, head_loops.T);
    for (int iconf = 0; iconf < head_loops.Njack; iconf++)
    {
        read_twopt(infile_loop, data_loop[iconf], head_loops);
        error(head.smearing[iconf].compare(head_loops.smearing[iconf]) != 0, 2, "main",
              "configuration order differ at %d\n flow file conf: %s\n loops conf: %s ", iconf,
              head.smearing[iconf].c_str(), head_loops.smearing[iconf].c_str());
    }
    //////////////////////////////////////////////////////////////
    // reading flow
    //////////////////////////////////////////////////////////////
    int ncorr_new = head.ncorr;                   // current number of correlators
    int Max_corr = head.ncorr + head_loops.ncorr; // max number of correlators

    double ****data = calloc_corr(head.Njack, Max_corr, head.T);

    printf("confs=%d\n", head.Njack);
    printf("ncorr=%d\n", head.ncorr);
    printf("kappa=%g\n", head.kappa);
    for (int iconf = 0; iconf < head.Njack; iconf++)
    {
        read_twopt(infile, data[iconf], head);
    }
    //////////////////////////////////////////////////////////////
    // correcting w0
    //////////////////////////////////////////////////////////////
    for (int i = 0; i < head_loops.ncorr; i++)
    {
        for (int j = 0; j < head.Njack; j++)
        {
            std::complex<double> bubble(0, 0);
            for (int t = 0; t < head_loops.T; t++)
            {
                bubble += std::complex<double>(data_loop[j][i][t][0], data_loop[j][i][t][1]);
            }

            for (int tf = 0; tf < head.T; tf++)
            {
                std::complex<double> Wt(data[j][6][tf][0], data[j][6][tf][1]);
                Wt *= bubble;
                data[j][head.ncorr + i][tf][0] = Wt.real();
                data[j][head.ncorr + i][tf][1] = Wt.imag();
            }
            // at tf=0 we put <loop>
            data[j][head.ncorr + i][0][0] = bubble.real();
            data[j][head.ncorr + i][0][1] = bubble.imag();
            // if (i == 27 && j == head.Njack - 1)
            //     printf("bubble dmu = %g  %g       <W(0)*bubble>= %g %g\n", bubble.real(), bubble.imag(), data[j][head.ncorr + i][0][0], data[j][head.ncorr + i][0][1]);
        }
    }
    free_corr(head_loops.Njack, head_loops.ncorr, head_loops.T, data_loop);
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
    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", option[3], option[6]);
    printf("writing output in :\n %s \n", namefile);
    FILE *outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s", option[3], option[4],
              option[6]);
    FILE *jack_file = open_file(namefile, "w+");
    // write_header_g2(jack_file, head);
    head.write_header(jack_file);

    // double ****data_bin = binning(confs, Max_corr, head.T, data, bin);
    // double ****data_bin = binning_toNb(confs, Max_corr, head.T, data, bin); // binning into N=bin, cutting out the reminder
    double ****data_bin = bin_intoN(data, Max_corr, head.T, confs, bin); // binning into N=bin with not integer
    free_corr(confs, Max_corr, head.T, data);
    double ****conf_jack = myres->create(Neff, Max_corr, head.T, data_bin);
    free_corr(Neff, Max_corr, head.T, data_bin);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print all the effective masses correlators
    // set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_meff_correlators", option[3],
              option[6]);
    FILE *outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_raw_correlators", option[3],
              option[6]);
    FILE *outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_shifted_correlators", option[3],
              option[6]);
    FILE *outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_log_meff_shifted", option[3],
              option[6]);
    FILE *outfile_log_meff_shifted = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_gamma", option[3], option[6]);
    FILE *out_gamma = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_kernel", option[3],
              option[6]);
    FILE *outfile_HLT_kernel = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_AoverB", option[3],
              option[6]);
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

    // print for frezzotti
    {
        double *tmp = (double *)malloc(sizeof(double) * Njack);
        for (int j = 0; j < Njack; j++)
        {
            tmp[j] = lhs_function_w0_eg(j, conf_jack, 150, fit_info);
        }
        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "%s/out/W_t150_jack.txt", option[3]);
        myres->write_jack_in_file(tmp, name);
        free(tmp);
    }
    // c++ 0 || r 1
    struct fit_result fit_W = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_function_w0_eg, "W(t)", fit_info,
        jack_file);
    check_correlatro_counter(0);

    double **tif = swap_indices(fit_info.Npar, Njack, fit_W.P);
    std::vector<double> swapped_x(fit_info.Nvar);
    std::vector<double> w0(Njack);

    // printf("fit:\n");
    // for (int i = 0; i < fit_info.Npar; i++)
    // {
    //     printf("%g  %g\n", fit_W.P[i][Njack - 1], myres->comp_error(fit_W.P[0]));
    // }
    // printf("%g  \n", fit_info.function(0, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[Njack - 1]));
    for (size_t j = 0; j < Njack; j++)
    {
        w0[j] = rtbis_func_eq_input(fit_info.function, 0 /*n*/, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[j], 0, 0.3, fit_info.tmin, fit_info.tmax, 1e-10, 2);
        w0[j] = int2flowt(w0[j]);
        w0[j] = std::sqrt(w0[j]);
    }
    printf("w0 = %g  %g\n", w0[Njack - 1], myres->comp_error(w0.data()));
    // free_fit_result(fit_info, fit_out);
    // fit_info.restore_default();

    print_result_in_file(outfile, w0.data(), "w0", 0.0, fit_info.tmin, fit_info.tmax);

    free_2(Njack, tif);
    //////////////////////////////////////////////////////////////
    // light mass correction
    //////////////////////////////////////////////////////////////

    double mean, err;
    int seed;
    int Nquark = 3;
    std::vector<double *> delta_amul(3);
    std::vector<std::string> q_name = {"l", "s", "c"};
    line_read_param(option, "delta_amul", mean, err, seed, namefile_plateaux);
    delta_amul[0] = myres->create_fake(mean, err, seed);
    line_read_param(option, "delta_amus", mean, err, seed, namefile_plateaux);
    delta_amul[1] = myres->create_fake(mean, err, seed);
    line_read_param(option, "delta_amuc", mean, err, seed, namefile_plateaux);
    delta_amul[2] = myres->create_fake(mean, err, seed);
    fit_info.n_ext_P = 1;
    fit_info.ext_P = (double **)malloc(sizeof(double) * fit_info.n_ext_P);
    std::vector<double> w0pdmu(Njack);

    // for (int i = 0; i<head_loops.gammas.size();i++){
    //     if (strcmp(head_loops.gammas[i].c_str(),"std_g5")==0)
    //         printf("id =%d %s %d\n", i,head_loops.gammas[i].c_str(), strcmp(head_loops.gammas[i].c_str(),"std_g5"));
    // }

    for (int iq = 0; iq < Nquark; iq++)
    {

        fit_info.ext_P[0] = delta_amul[iq];

        printf("dmu%s= %g  %g\n", q_name[iq].c_str(), fit_info.ext_P[0][Njack - 1], myres->comp_error(fit_info.ext_P[0]));

        int Ng = head_loops.gammas.size();
        fit_info.corr_id = {6, head.ncorr + 27 + iq * Ng};
        fit_info.myen = {1}; // reim of corr_id[1] (the correction)

        for (int t = 1; t < head.T; t++)
        {
            if (lhs_function_Wt_p_dmcorr(Njack - 1, conf_jack, t, fit_info) > 0.3)
            {
                fit_info.tmin = t - 2;
                fit_info.tmax = t + 1;
                break;
            }
        }
        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "W+mu%s_correction(t)", q_name[iq].c_str());
        struct fit_result fit_Wpdmu = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_function_Wt_p_dmcorr, name, fit_info,
            jack_file);
        check_correlatro_counter(1 + iq);

        tif = swap_indices(fit_info.Npar, Njack, fit_Wpdmu.P);

        for (size_t j = 0; j < Njack; j++)
        {
            w0pdmu[j] = rtbis_func_eq_input(fit_info.function, 0 /*n*/, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[j], 0, 0.3, fit_info.tmin - 2, fit_info.tmax + 2, 1e-10, 2);
            w0pdmu[j] = int2flowt(w0pdmu[j]);
            w0pdmu[j] = std::sqrt(w0pdmu[j]);
        }
        mysprintf(name, NAMESIZE, "w0+mu%s_correction", q_name[iq].c_str());
        printf("%s = %g  %g\n", name, w0pdmu[Njack - 1], myres->comp_error(w0pdmu.data()));
        printf("unbiased %s = %g  %g\n", name, myres->comp_mean_unbias(w0pdmu.data()), myres->comp_error(w0pdmu.data()));
        print_result_in_file(outfile, w0pdmu.data(), name, 0.0, fit_info.tmin, fit_info.tmax);

        double *Delta_mul_w0 = myres->create_copy(w0.data());
        myres->sub(Delta_mul_w0, Delta_mul_w0, w0pdmu.data());
        mysprintf(name, NAMESIZE, "Delta_mu%s_w0", q_name[iq].c_str());
        print_result_in_file(outfile, Delta_mul_w0, name, 0.0, fit_info.tmin, fit_info.tmax);

        free_2(Njack, tif);
        free(Delta_mul_w0);
        fit_Wpdmu.clear();
    }
    fit_info.restore_default();

    //////////////////////////////////////////////////////////////
    // all mass corrections
    //////////////////////////////////////////////////////////////

    fit_info.Nvar = 1;
    fit_info.Npar = 4;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.codeplateaux = true;
    fit_info.function = poly3;
    fit_info.linear_fit = true;
    fit_info.T = head.T;
    int Ng = head_loops.gammas.size();
    fit_info.corr_id = {
        6,
        head.ncorr + 27 + 0 * Ng,
        head.ncorr + 27 + 1 * Ng,
        head.ncorr + 27 + 2 * Ng,
    };
    fit_info.n_ext_P = 3;
    fit_info.ext_P = (double **)malloc(sizeof(double) * fit_info.n_ext_P);
    fit_info.ext_P[0] = delta_amul[0];
    fit_info.ext_P[1] = delta_amul[1];
    fit_info.ext_P[2] = delta_amul[2];
    fit_info.myen = {1}; // reim of corr_id[1] (the correction)

    for (int t = 1; t < head.T; t++)
    {
        if (lhs_function_Wt_p_dm_all_corr(Njack - 1, conf_jack, t, fit_info) > 0.3)
        {
            fit_info.tmin = t - 2;
            fit_info.tmax = t + 1;
            break;
        }
    }
    char name[NAMESIZE];
    mysprintf(name, NAMESIZE, "W+all_mu_correction(t)");
    struct fit_result fit_Wpdmu = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_function_Wt_p_dm_all_corr, name, fit_info,
        jack_file);
    check_correlatro_counter(4);

    tif = swap_indices(fit_info.Npar, Njack, fit_Wpdmu.P);

    for (size_t j = 0; j < Njack; j++)
    {
        w0pdmu[j] = rtbis_func_eq_input(fit_info.function, 0 /*n*/, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[j], 0, 0.3, fit_info.tmin - 1, fit_info.tmax + 1, 1e-10, 2);
        w0pdmu[j] = int2flowt(w0pdmu[j]);
        w0pdmu[j] = std::sqrt(w0pdmu[j]);
    }
    mysprintf(name, NAMESIZE, "w0+all_mu_correction");
    printf("%s = %g  %g\n", name, w0pdmu[Njack - 1], myres->comp_error(w0pdmu.data()));
    print_result_in_file(outfile, w0pdmu.data(), name, 0.0, fit_info.tmin, fit_info.tmax);

    double *Delta_mul_w0 = myres->create_copy(w0.data());
    myres->sub(Delta_mul_w0, Delta_mul_w0, w0pdmu.data());
    mysprintf(name, NAMESIZE, "Delta_all_mu_w0");
    print_result_in_file(outfile, Delta_mul_w0, name, 0.0, fit_info.tmin, fit_info.tmax);

    free_2(Njack, tif);
    free(Delta_mul_w0);
    fit_Wpdmu.clear();

    fit_info.restore_default();

    //////////////////////////////////////////////////////////////
    //  mass derivate
    //////////////////////////////////////////////////////////////

    fit_info.Nvar = 1;
    fit_info.Npar = 4;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.codeplateaux = true;

    fit_info.function = poly3;
    fit_info.linear_fit = true;
    fit_info.T = head.T;
    fit_info.n_ext_P = 0;

    // for (int i = 0; i<head_loops.gammas.size();i++){
    //     if (strcmp(head_loops.gammas[i].c_str(),"std_g5")==0)
    //         printf("id =%d %s %d\n", i,head_loops.gammas[i].c_str(), strcmp(head_loops.gammas[i].c_str(),"std_g5"));
    // }
    fit_info.tmin = 1;
    fit_info.tmax = 2;
    double *tmp = (double *)malloc(sizeof(double) * Njack);
    for (int iq = 0; iq < Nquark; iq++)
    {

        int Ng = head_loops.gammas.size();
        fit_info.corr_id = {6, head.ncorr + 27 + iq * Ng};
        fit_info.myen = {1}; // reim of corr_id[1] (the correction)

        // print for frezzotti
        for (int j = 0; j < Njack; j++)
        {
            tmp[j] = lhs_function_Wt_der_mu(j, conf_jack, 150, fit_info);
        }
        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "%s/out/der_W_mu%s_correction_t150_jack", option[3], q_name[iq].c_str());
        myres->write_jack_in_file(tmp, name);

        // for (int t = 1; t < head.T; t++)
        // {
        //     if (lhs_function_Wt_p_dmcorr(Njack - 1, conf_jack, t, fit_info) > 0.3)
        //     {
        //         fit_info.tmin = t - 2;
        //         fit_info.tmax = t + 1;
        //         break;
        //     }
        // }
        mysprintf(name, NAMESIZE, "der_W_mu%s_correction(t)", q_name[iq].c_str());
        struct fit_result fit_Wpdmu = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_function_Wt_der_mu, name, fit_info,
            jack_file);
        check_correlatro_counter(5 + iq);

        // tif = swap_indices(fit_info.Npar, Njack, fit_Wpdmu.P);

        // for (size_t j = 0; j < Njack; j++)
        // {
        //     w0pdmu[j] = rtbis_func_eq_input(fit_info.function, 0 /*n*/, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[j], 0, 0.3, fit_info.tmin, fit_info.tmax, 1e-10, 2);
        //     w0pdmu[j] = int2flowt(w0pdmu[j]);
        //     w0pdmu[j] = std::sqrt(w0pdmu[j]);
        // }
        // mysprintf(name, NAMESIZE, "w0+mu%s_correction", q_name[iq].c_str());
        // printf("%s = %g  %g\n", name, w0pdmu[Njack - 1], myres->comp_error(w0pdmu.data()));
        // print_result_in_file(outfile, w0pdmu.data(), name, 0.0, fit_info.tmin, fit_info.tmax);

        // double *Delta_mul_w0 = myres->create_copy(w0.data());
        // myres->sub(Delta_mul_w0, Delta_mul_w0, w0pdmu.data());
        // mysprintf(name, NAMESIZE, "Delta_mu%s_w0", q_name[iq].c_str());
        // print_result_in_file(outfile, Delta_mul_w0, name, 0.0, fit_info.tmin, fit_info.tmax);

        // free_2(Njack, tif);
        // free(Delta_mul_w0);
        // fit_Wpdmu.clear();
    }
    free(tmp);
    fit_info.restore_default();
}