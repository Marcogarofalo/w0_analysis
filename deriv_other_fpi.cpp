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
#include "match_confs.hpp"

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
    mysprintf(option[6], NAMESIZE, "%s_%s", argv[3], argv[8]);         // infile
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
    error(!(argc == 10), 1, "main ",
        "usage:././w0_rew -p path file -bin $bin  jack/boot   fileP5P5  rew  basename_palteau\n separate "
        "path and file please   argc =%d", argc);

    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, argv[6]);
    printf("resampling: %s\n", resampling);

    char** option = argv_to_options(argv);


    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");




    //////////////////////////////////// read and setup header
    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[3]);
    FILE* infile_A0 = open_file(namefile, "r");
    generic_header head_A0;
    head_A0.read_header(infile_A0);
    init_global_head(head_A0);

    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[7]);
    FILE* infile_P5 = open_file(namefile, "r");
    generic_header head_P5;
    head_P5.read_header(infile_P5);
    init_global_head(head_P5);

    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[8]);
    FILE* infile_rew = open_file(namefile, "r");
    generic_header head_rew;
    head_rew.read_header(infile_rew);
    // init_global_head(head_rew); // don't!! it will set the file_head.l0=1

    printf("mu = %g\n", head_A0.mus[0]);

    auto [idx_A0, idx_P5] = matching_indices(head_A0.smearing, head_P5.smearing);

    std::vector<std::string> confs_2pt(idx_A0.size());
    for (int i = 0;i < idx_A0.size();i++)
        confs_2pt[i] = head_A0.smearing[idx_A0[i]];

    error(hasDuplicates(head_A0.smearing), 1, "main", "duplicates in confs of A0");
    error(hasDuplicates(head_P5.smearing), 1, "main", "duplicates in confs of P5");
    error(hasDuplicates(head_rew.smearing), 1, "main", "duplicates in confs of rew");

    auto results = findAllIndices(head_A0.smearing, head_P5.smearing, head_rew.smearing);
    printf("confs A0P5 %ld\n", head_A0.smearing.size());
    printf("confs P5P5 %ld\n", head_P5.smearing.size());
    printf("confs rew  %ld\n", head_rew.smearing.size());
    printf("conf matching: %ld\n", results.size());
    // for (const auto& m : results) {
    //     std::cout << "Value: " << m.value 
    //               << " | Indices: v1[" << m.idx1 
    //               << "], v2[" << m.idx2 
    //               << "], v3[" << m.idx3 << "]\n";
    // }
    // error(head_A0.Njack != head_P5.Njack, 1, "main", "A0 file does not have the same confs of P5 file\n A0 %d\n P5 %d\n",head_A0.Njack ,head_P5.Njack );
    // for (int i = 0;i < head_A0.Njack;i++) {
    //     error(head_A0.smearing[i].compare(head_P5.smearing[i]) != 0, 2, "main",
    //         "configuration order differ at %d\n flow file conf: %s\n loops conf: %s ", i,
    //         head_A0.smearing[i].c_str(), head_P5.smearing[i].c_str());
    // }

    // auto [idx_a, idx_b] = matching_indices(head_A0.smearing, head_rew.smearing);
    // printf("matching confs: %ld\n", idx_a.size());


    //////////////////////////////////////////////////////////////
    // reading flow
    //////////////////////////////////////////////////////////////
    head_A0.Njack = results.size();
    head_P5.Njack = results.size();
    head_rew.Njack = results.size();

    int ncorr_new = head_A0.ncorr + 1;        // current number of correlators
    int Max_corr = head_A0.ncorr + 3 + 1; // max number of correlators

    double**** data = calloc_corr(head_A0.Njack, Max_corr, head_A0.T);
    double*** d_tmp = malloc_3<double>(1, head_A0.T, 2);

    for (int iconf = 0; iconf < head_A0.Njack; iconf++) {
        int seek_pos = head_A0.struct_size + results[iconf].idx1 * (sizeof(double) * (head_A0.T * head_A0.ncorr * 2) + sizeof(int));
        fseek(infile_A0, seek_pos, SEEK_SET);
        read_twopt(infile_A0, data[iconf], head_A0);

        seek_pos = head_P5.struct_size + results[iconf].idx2 * (sizeof(double) * (head_P5.T * head_P5.ncorr * 2) + sizeof(int));
        fseek(infile_P5, seek_pos, SEEK_SET);
        read_twopt(infile_P5, d_tmp, head_P5);

        for (int t = 0; t < head_A0.T; t++)
            data[iconf][1][t][0] = d_tmp[0][t][0];
    }

    double**** data_rew = calloc_corr(head_rew.Njack, 2, head_rew.T);
    for (int iconf = 0; iconf < head_rew.Njack; iconf++) {
        int seek_pos = head_rew.struct_size + results[iconf].idx3 * (sizeof(double) * (head_rew.T * head_rew.ncorr * 2) + sizeof(int));

        fseek(infile_rew, seek_pos, SEEK_SET);
        read_twopt(infile_rew, data_rew[iconf], head_rew);

    }
    free_3(1, head_A0.T, d_tmp);

    //////////////////////////////////////////////////////////////
    // rew
    //////////////////////////////////////////////////////////////
    for (int i = 0; i < ncorr_new; i++) {
        for (int j = 0; j < head_A0.Njack; j++) {

            double r = data_rew[j][0][0][1]; // the im part is the exponentiated subtracted reweighting factor
            for (int tf = 0; tf < head_A0.T; tf++) {
                data[j][ncorr_new + i][tf][0] = data[j][i][tf][0] * r;
                data[j][ncorr_new + i][tf][1] = data[j][i][tf][1] * r;
            }
        }
    }
    for (int j = 0; j < head_A0.Njack; j++) {
        for (int tf = 0; tf < head_A0.T; tf++) {
            data[j][ncorr_new * 2][tf][0] = data_rew[j][0][0][1];
            data[j][ncorr_new * 2][tf][1] = 0;
        }
    }

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
    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", option[3], option[6]);
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

    /// normalize rew factor
    for (int i = 0; i < ncorr_new; i++) {
        for (int j = 0; j < head_A0.Njack; j++) {

            double r = conf_jack[j][ncorr_new * 2][0][0]; // the im part is the exponentiated subtracted reweighting factor
            for (int tf = 0; tf < head_A0.T; tf++) {
                conf_jack[j][ncorr_new + i][tf][0] = conf_jack[j][ncorr_new + i][tf][0] / r;
                conf_jack[j][ncorr_new + i][tf][1] = conf_jack[j][ncorr_new + i][tf][1] / r;
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print all the effective masses correlators
    // set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_meff_correlators", option[3], option[6]);
    FILE* outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_raw_correlators", option[3], option[6]);
    FILE* outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_shifted_correlators", option[3], option[6]);
    FILE* outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_log_meff_shifted", option[3], option[6]);
    FILE* outfile_log_meff_shifted = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_kernel", option[3], option[6]);
    FILE* outfile_HLT_kernel = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_AoverB", option[3], option[6]);
    FILE* outfile_HLT_AoverB = open_file(namefile, "w+");

    char save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");

    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;

    for (int icorr = 0; icorr < ncorr_new; icorr++) {
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

    ////////////////////////////////////////////////////////////
    // symmetrize
    ////////////////////////////////////////////////////////////
    // symmetrise_jackboot(Njack, 0, head_A0T, conf_jack);
    // symmetrise_jackboot(Njack, 1, head_A0T, conf_jack, -1);

    //////////////////////////////////////////////////////////////
    //  read m^iso
    //////////////////////////////////////////////////////////////
    std::vector<double*> amuiso(3);
    double mean, err;
    int seed;
    char** option_p = argv_to_options(argv);
    mysprintf(option_p[6], NAMESIZE, "%s", argv[9]);         // basename plateau
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
    struct fit_type fit_info;
    fit_info.codeplateaux = true;

    line_read_param(option, "M_{PS}", fit_info.tmin, fit_info.tmax, seed, namefile_plateaux);

    double* M_PS = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, 1, "M_{PS}", M_eff_T, jack_file);
    // free(M_PS);
    check_correlatro_counter(0);


    double* M_PS_mu = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, 3, "M_{PS}_mu", M_eff_T, jack_file, fit_info);
    // free(M_PS);
    check_correlatro_counter(1);

    double dmu = head_rew.mus[0] - head_rew.oranges[0];
    printf("dmu: %.12g mu1 =  %g   mu2 = %g\n", dmu, head_rew.mus[0], head_rew.oranges[0]);
    double* deriv_M = myres->create_zero();
    for (int j = 0; j < Njack;j++) {
        deriv_M[j] = (M_PS[j] - M_PS_mu[j]) / dmu;
    }
    printf("deriv M: %.12g   %.12g\n", deriv_M[Njack - 1], myres->comp_error(deriv_M));

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
    fit_info.ext_P = (double**)malloc(sizeof(double*) * 3);
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

    double* Zfpi = myres->create_copy(fpi.P[0]);
    myres->mult(Zfpi, Z, fpi.P[0]);
    printf("Zfpi: %.12g   %.12g\n", Zfpi[Njack - 1], myres->comp_error(Zfpi));


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
    printf("deriv: %.12g   %.12g\n", deriv[Njack - 1], myres->comp_error(deriv));
    write_jack(deriv, Njack, jack_file);
    check_correlatro_counter(6);

    std::string name(option[6]);
    // 1. Find the first underscore
    first = name.find('_');

    // 2. Find the second underscore starting search from the position after the first
    second = name.find('_', first + 1);
    if (second != std::string::npos) {
        // 3. Extract from the start (0) up to the second underscore
        std::string result = name.substr(0, second);
        std::cout << result << std::endl; // Output: B64_TM
    }

    //////////////////////////////////////////////////////////////
    // fpi WTI
    //////////////////////////////////////////////////////////////
    std::string name_rew = argv[8];
    std::cout << name_rew << "\n";
    first = name_rew.find('_');
    second = name_rew.find('_', first + 1);
    label = name_rew.substr(first + 1, second - first - 1);
    std::cout << label << "\n";
    if (label.compare("charm") == 0)
        label = "rewcOS";
    if (label.compare("light") == 0)
        label = "rewlOS";
    if (label.compare("strange") == 0)
        label = "rewsOS";

    std::string name_jack_fpi = "deriv/deriv_fpi_P5A0_" + std::string(argv[3]) + "_" + std::string(argv[8]) + ".jack_txt";
    myres->write_jack_in_file(deriv, name_jack_fpi.c_str());

    name_jack_fpi = "deriv/fpi_P5A0_" + std::string(argv[3]) + ".jack_txt";
    double* Zf = myres->create_copy(fpi.P[0]);
    myres->mult(Zf, Z, fpi.P[0]);
    myres->write_jack_in_file(Zf, name_jack_fpi.c_str());



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
        fit_info.ext_P[1][j] = head_P5.mus[0];
        fit_info.ext_P[2][j] = head_P5.mus[0];
        // fit_info.ext_P[1][j] = head_rew.oranges[0];
        // fit_info.ext_P[2][j] = head_rew.oranges[0];
    }
    mysprintf(namefile, NAMESIZE, "f_{PS}_%s", label.c_str());
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

    write_jack(amuiso[0], Njack, jack_file); check_correlatro_counter(10);
    write_jack(amuiso[1], Njack, jack_file); check_correlatro_counter(11);
    write_jack(amuiso[2], Njack, jack_file); check_correlatro_counter(12);

    write_jack(amusim[0], Njack, jack_file); check_correlatro_counter(13);
    write_jack(amusim[1], Njack, jack_file); check_correlatro_counter(14);
    write_jack(amusim[2], Njack, jack_file); check_correlatro_counter(15);

    write_jack(a_fm, Njack, jack_file); check_correlatro_counter(16);

    double* mul_over_fpi_deriv = myres->create_copy(deriv);
    for (int j = 0; j < Njack;j++) {
        mul_over_fpi_deriv[j] = (amuiso[0][j] / fpi.P[0][j]) * deriv[j];
    }
    printf("mul_over_fpi_deriv: %.12g   %.12g\n", mul_over_fpi_deriv[Njack - 1], myres->comp_error(mul_over_fpi_deriv));
    double* mul_over_fpi_deriv_WTI = myres->create_copy(deriv_WTI);
    for (int j = 0; j < Njack;j++) {
        mul_over_fpi_deriv_WTI[j] = (amuiso[0][j] / fpi.P[0][j]) * deriv_WTI[j];
    }
    printf("mul_over_fpi_deriv_WTI: %.12g   %.12g\n", mul_over_fpi_deriv_WTI[Njack - 1], myres->comp_error(mul_over_fpi_deriv_WTI));
    write_jack(mul_over_fpi_deriv, Njack, jack_file); check_correlatro_counter(17);
    write_jack(mul_over_fpi_deriv_WTI, Njack, jack_file); check_correlatro_counter(18);

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 0;
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = head_A0.T;
    fit_info.n_ext_P = 5;
    // fit_info.malloc_ext_P();
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
    fit_info.ext_P[0] = M_PS;
    fit_info.ext_P[1] = me.P[0];
    fit_info.ext_P[2] = M_PS_mu;
    fit_info.ext_P[3] = me_mu.P[0];
    fit_info.ext_P[4] = Z;
    fit_info.ave_P = { dmu };


    fit_info.corr_id = { 0 ,2 }
    ;
    struct fit_result plateau_fpi = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_plateau_fpi_P5A0, "pateau_deriv_fpi_P5A0", fit_info,
        jack_file);
    check_correlatro_counter(19);


    double* mul = myres->create_fake(head_P5.mus[0], 1e-20, -1);
    write_jack(mul, Njack, jack_file); check_correlatro_counter(20);

    write_jack(Z, Njack, jack_file); check_correlatro_counter(21);

    fit_info.restore_default();

    write_jack(deriv_M, Njack, jack_file); check_correlatro_counter(22);
}