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
constexpr double hbarc = 1.97326979e-6; // GeV*fm

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
    int ncorr_new = head.ncorr;                 // current number of correlators
    int Max_corr = head.ncorr + head.ncorr * 3; // max number of correlators

    double ****data = calloc_corr(head.Njack, Max_corr, head.T);

    printf("confs=%d\n", head.Njack);
    printf("ncorr=%d\n", head.ncorr);
    printf("kappa=%g\n", head.kappa);
    for (int iconf = 0; iconf < head.Njack; iconf++)
    {
        read_twopt(infile, data[iconf], head);
    }
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
    // correcting P5P5 and V0P5
    //////////////////////////////////////////////////////////////
    for (int i = 0; i < head.ncorr; i++)
    {
        int Ng = head_loops.gammas.size();
        for (int iq = 0; iq < 3; iq++)
        {
            int idloop_q = 27 + iq * Ng;
            int id_final = head.ncorr + i + iq * head.ncorr;

            for (int j = 0; j < head.Njack; j++)
            {
                std::complex<double> bubble(0, 0);
                for (int t = 0; t < head_loops.T; t++)
                {
                    bubble += std::complex<double>(data_loop[j][idloop_q][t][0], data_loop[j][idloop_q][t][1]);
                }

                for (int tf = 0; tf < head.T; tf++)
                {
                    std::complex<double> Wt(data[j][i][tf][0], data[j][i][tf][1]);
                    Wt *= bubble;
                    data[j][id_final][tf][0] = Wt.real();
                    data[j][id_final][tf][1] = Wt.imag();
                }
                data[j][id_final][0][0] = bubble.real();
                data[j][id_final][0][1] = bubble.imag();
            }
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

    double *M_PS = plateau_correlator_function(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, 0, "M_{PS}", M_eff_T, jack_file);
    // free(M_PS);
    check_correlatro_counter(0);

    //////////////// mpcac

    struct fit_type fit_info;
    struct fit_result fit_out;

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 0;
    // fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
    // fit_info.ext_P[0] = something;
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = head.T;
    fit_info.corr_id = {1, 0};

    struct fit_result fit_mpcac = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_mpcac, "mpcac", fit_info,
        jack_file);
    check_correlatro_counter(1);
    // free_fit_result(fit_info, fit_out);
    // fit_info.restore_default(); // keep T and other stuff

    double mean, err;
    int seed;
    int Nquark = 3;
    std::vector<double *> delta_amul(3);
    std::vector<std::string> q_name = {"l", "s", "c"};
    std::string name_of_flow_run(option[6]);
    name_of_flow_run.replace(name_of_flow_run.begin(), name_of_flow_run.begin() + 10, "flow");
    std::cout << "name_of_flow_run: " << name_of_flow_run << "\n";

    std::string save_data_file_name(option[6]);
    mysprintf(option[6], NAMESIZE, "%s", name_of_flow_run.c_str());

    line_read_param(option, "delta_amul", mean, err, seed, namefile_plateaux);
    delta_amul[0] = myres->create_fake(mean, err, seed);
    line_read_param(option, "delta_amus", mean, err, seed, namefile_plateaux);
    delta_amul[1] = myres->create_fake(mean, err, seed);
    line_read_param(option, "delta_amuc", mean, err, seed, namefile_plateaux);
    delta_amul[2] = myres->create_fake(mean, err, seed);

    // restore option[6]: data file name
    mysprintf(option[6], NAMESIZE, "%s", save_data_file_name.c_str());

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

    fit_info.n_ext_P = 1;
    fit_info.ext_P = (double **)malloc(sizeof(double) * fit_info.n_ext_P);

    //////////////////////////////////////////////////////////////
    // corrections
    //////////////////////////////////////////////////////////////
    std::vector<fit_result> fit_M_PSpdmu(Nquark);
    for (int iq = 0; iq < Nquark; iq++)
    {

        fit_info.ext_P[0] = delta_amul[iq];

        printf("dmu%s= %g  %g\n", q_name[iq].c_str(), fit_info.ext_P[0][Njack - 1], myres->comp_error(fit_info.ext_P[0]));

        int Ng = head_loops.gammas.size();
        fit_info.corr_id = {0, head.ncorr + 0 + iq * head.ncorr};
        fit_info.myen = {1}; // reim of corr_id[1] (the correction)

        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "M_{PS}+mu%s_correction", q_name[iq].c_str());
        fit_M_PSpdmu[iq] = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_function_M_PS_p_dmcorr, name, fit_info,
            jack_file);
        check_correlatro_counter(2 + iq * 2);

        // print for frezzotti
        char name_f_jack[NAMESIZE];
        mysprintf(name_f_jack, NAMESIZE, "%s/out/%s_jack.txt", option[3], name);
        myres->write_jack_in_file(fit_M_PSpdmu[iq].P[0], name_f_jack);

        fit_info.corr_id = {1, 0, head.ncorr + 1 + iq * head.ncorr, head.ncorr + 0 + iq * head.ncorr};
        mysprintf(name, NAMESIZE, "mpcac+mu%s_correction", q_name[iq].c_str());
        struct fit_result fit_mpcacpdmu = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_function_Wt_p_dmcorr, name, fit_info,
            jack_file);

        // fit_M_PSpdmu.clear();
        fit_mpcacpdmu.clear();
    }
    fit_info.restore_default();

    double *derM = myres->create_copy(M_PS);
    for (int iq = 0; iq < Nquark; iq++)
    {

        for (int j = 0; j < Njack; j++)
        {
            double mr = fit_M_PSpdmu[iq].P[0][j];
            double m = M_PS[j];
            derM[j] = (mr - m) / delta_amul[iq][j];
            if (j == Njack - 1)
                printf("dM_{PS}/dmu%s = %g   %g   \n", q_name[iq].c_str(), derM[j], myres->comp_error(derM));
        }
        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "dM_{PS}/dmu%s_correction", q_name[iq].c_str());
        print_result_in_file(outfile, derM, name, 0, 0, 0);
    }

    //////////////////////////////////////////////////////////////
    // fpi
    //////////////////////////////////////////////////////////////

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.T = head.T;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 3;
    fit_info.ext_P = malloc_2<double>(fit_info.n_ext_P, fit_info.Njack);
    for (int j = 0; j < fit_info.Njack; j++)
    {
        fit_info.ext_P[0][j] = M_PS[j];
        fit_info.ext_P[1][j] = head.mus[0];
        fit_info.ext_P[2][j] = head.mus[0];
    }
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;

    fit_info.corr_id = {0};

    struct fit_result f_PS = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_function_f_PS, "f_{PS}", fit_info, jack_file);

    std::vector<fit_result> fit_f_PSpdmu(Nquark);
    free_2(fit_info.n_ext_P, fit_info.ext_P);
    fit_info.n_ext_P = 4;
    fit_info.ext_P = malloc_2<double>(fit_info.n_ext_P, fit_info.Njack);
    for (int j = 0; j < fit_info.Njack; j++)
    {
        fit_info.ext_P[1][j] = head.mus[0];
        fit_info.ext_P[2][j] = head.mus[0];
    }

    for (int iq = 0; iq < Nquark; iq++)
    {

        myres->copy(fit_info.ext_P[0], fit_M_PSpdmu[iq].P[0]);
        myres->copy(fit_info.ext_P[3], delta_amul[iq]);

        int Ng = head_loops.gammas.size();
        fit_info.corr_id = {0, head.ncorr + 0 + iq * head.ncorr};
        fit_info.myen = {1}; // reim of corr_id[1] (the correction)

        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "f_{PS}+mu%s_correction", q_name[iq].c_str());
        fit_f_PSpdmu[iq] = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_function_f_PS_p_dmcorr, name, fit_info,
            jack_file);
        // check_correlatro_counter(2 + iq * 2);
    }
    for (int iq = 0; iq < Nquark; iq++)
    {

        for (int j = 0; j < Njack; j++)
        {
            double mr = fit_M_PSpdmu[iq].P[0][j];
            double m = M_PS[j];
            double fr = fit_f_PSpdmu[iq].P[0][j];
            double f = f_PS.P[0][j];
            double dmu = delta_amul[iq][j];
            derM[j] = (mr - m) / dmu - (m / f) * ((fr - f) / dmu);
            // if (j == Njack - 1)
            //     printf(" %g   %g   %g   %g  %g  --> %g\n", mr, m, fr, f, dmu,derM[j]);
        }
        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "f_{PS}d(M_{PS}/f_{PS})/dmu%s_correction", q_name[iq].c_str());
        print_result_in_file(outfile, derM, name, 0, 0, 0);

        for (int j = 0; j < Njack; j++)
        {
            double fr = fit_f_PSpdmu[iq].P[0][j];
            double f = f_PS.P[0][j];
            double dmu = delta_amul[iq][j];
            derM[j] = (fr - f) / dmu;
        }
        mysprintf(name, NAMESIZE, "df_{PS}/dmu%s_correction", q_name[iq].c_str());
        print_result_in_file(outfile, derM, name, 0, 0, 0);
        printf("%s = %g +- %g\n", name, derM[Njack - 1], myres->comp_error(derM));

        // for (int j = 0; j < fit_info.Njack; j++)
        // {
        //     double fr = fit_f_PSpdmu[iq].P[0][j];
        //     double f = f_PS.P[0][j];
        //     double dmu = delta_amul[iq][j];
        //     derM[j] = f + derM[j] * dmu;
        // }
        // mysprintf(name, NAMESIZE, "f_{PS}(iso)_mu%s_correction", q_name[iq].c_str());
        // print_result_in_file(outfile, derM, name, 0, 0, 0);
        // printf("%s = %g +- %g\n", name, derM[Njack - 1], myres->comp_error(derM));
    }

    fit_info.restore_default();
    // compute derivative expanding
    ////////////////////////////////////////////////////////////
    //  mass correction
    ////////////////////////////////////////////////////////////
    std::vector<fit_result> Delta_M(Nquark);
    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 1;
    fit_info.ext_P = (double **)malloc(sizeof(double *) * 5);
    fit_info.ext_P[0] = M_PS;
    fit_info.ext_P[1] = f_PS.P[0];
    // fit_info.ext_P[2] = dM_PS;
    fit_info.ext_P[3] = (double *)malloc(sizeof(double) * Njack);
    fit_info.ext_P[4] = (double *)malloc(sizeof(double) * Njack);
    for (int j = 0; j < Njack; j++)
    {
        fit_info.ext_P[3][j] = head.mus[0];
        fit_info.ext_P[4][j] = head.mus[0];
    }

    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = head.T;

    for (int iq = 0; iq < Nquark; iq++)
    {

        fit_info.n_ext_P = 1;
        //////////////////////////////////////////  up
        fit_info.myen = {1}; // reim of corr_id[1] (the correction)
        
        int id_dmu_u_pi = head.ncorr + 0 + iq * head.ncorr;
        int id_PS = 0;
        fit_info.corr_id = {id_PS, id_dmu_u_pi };
        char namefit[NAMESIZE];
        mysprintf(namefit, NAMESIZE, "Delta_mu%s_M_{PS}", q_name[iq].c_str());

        Delta_M[iq] = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_dM_sea, namefit, fit_info,
            jack_file);

        fit_info.n_ext_P = 5;
        mysprintf(namefit, NAMESIZE, "Delta_mu%s_f_{PS}", q_name[iq].c_str());
        // correction to fpi
        fit_info.ext_P[2] = Delta_M[iq].P[0];
        Delta_M[iq] = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char *)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_dfpi_sea, namefit, fit_info,
            jack_file);
    }
}