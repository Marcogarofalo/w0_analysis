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
#include <ranges>
#include <algorithm>
#include <iterator>

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

int read_single_Nobs(FILE* stream, int header_size, int Njack) {
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

double rhs_a2_a4(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    r = P[0] + a2 * P[1] + a2 * a2 * P[2];
    return r;
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

double lhs_fun_masses(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double m = gjack.en[e].jack[fit_info.corr_id[0 + n * 2]][j];
    double Z = gjack.en[e].jack[fit_info.corr_id[1 + n * 2]][j];

    int ix = get_count_x_in_lhs(n, e, fit_info);
    // if (j == 0)printf("n=%d  e=%d ix=%d\n", n, e, ix);
    double a = std::sqrt(fit_info.x[0][ix][j]) / hbarc;
    return m / (Z * a);
}
double rhs_masses(int n, int Nvar, double* x, int Npar, double* P) {
    double a2 = x[0];
    int Nz = 2;
    int Nferm = 3;
    int iz = n % 2;
    int ferm = n / 2;
    // P[0,1,2] = c(f)
    // c0(f)* (1 + a ^ 2 * c1(Z, f) + a ^ 2log(a ^ 2 * Lam ^ 2) * clog(Z))
    double cf = P[ferm];
    double c1 = P[Nferm + iz * Nferm + ferm];
    double clogZ = P[Nferm + Nz * Nferm + iz]; ;
    double lam = 295 / hbarc;
    double lam2 = lam * lam;

    double r = cf * (1 + a2 * c1);
    if (a2 > 1e-10)
        r += cf * (a2 * std::log(a2 * lam2) * clogZ);
    return r;
}
double  rhs_mass_light(int n, int Nvar, double* x, int Npar, double* P) {
    double a2 = x[0];
    int Nz = 2;
    int Nferm = 1;
    int iz = n % 2;
    int ferm = n / 2;
    // c0(f)* (1 + a ^ 2 * c1(Z, f) + a ^ 2log(a ^ 2 * Lam ^ 2) * clog(Z))
    double cf = P[ferm];
    double c1 = P[Nferm + iz * Nferm + ferm];
    double clogZ = P[Nferm + Nz * Nferm + iz]; ;
    double lam = 295 / hbarc;
    double lam2 = lam * lam;

    double r = cf * (1 + a2 * c1);
    if (a2 > 1e-10)
        r += cf * (a2 * std::log(a2 * lam2) * clogZ);
    return r;
}
double lhs_w0_sim_aiso(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    int id =(n == 0) ?fit_info.corr_id[0] : fit_info.corr_id[1];
    double w0_a = gjack.en[e].jack[id][j]; 

    double a = myres->mean(gjack.en[e].jack[fit_info.corr_id[2]]); 
    return w0_a*a;
}

double lhs_w0_sim(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double w0_a = gjack.en[e].jack[fit_info.corr_id[0]][j]; // d(w0/a)/d(a*mu)
    double afpi = (n == 0) ? gjack.en[e].jack[fit_info.corr_id[1]][j] : gjack.en[e].jack[fit_info.corr_id[2]][j];

    double a = afpi / (fpi_MeV / hbarc);
    return w0_a * afpi;
}
double lhs_w0_sim1(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double w0_a = (n == 0) ? gjack.en[e].jack[fit_info.corr_id[0]][j]: gjack.en[e].jack[fit_info.corr_id[3]][j] ; // d(w0/a)/d(a*mu)
    double afpi = (n == 0) ? gjack.en[e].jack[fit_info.corr_id[1]][j] : gjack.en[e].jack[fit_info.corr_id[2]][j];

    double a = afpi / (fpi_MeV / hbarc);
    return w0_a * afpi;
}
double lhs_w0_full(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double w0 =  gjack.en[e].jack[fit_info.corr_id[0]][j];
    double afpi_fpiphys =  gjack.en[e].jack[fit_info.corr_id[1]][j];

    return w0*(fpi_MeV/hbarc); // w0*fpi full corrected
}


double lhs_fun_N3(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r = gjack.en[e].jack[fit_info.corr_id[n]][j]; // d(w0/a)/d(a*mu)
    return r;
}

double rhs_a2_N3(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    r = P[0] + a2 * P[1 + n];
    return r;
}
double rhs_a2_N3_a4(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    r = P[0] + a2 * P[1 + n] + a2 * a2 * P[1 + 3 + n];
    return r;
}
double rhs_a2_N3_a4WTI(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    r = P[0] + a2 * P[1 + n];
    if (n == 0) r += a2 * a2 * P[4];
    return r;
}
double rhs_a2_N3_a4tm(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    r = P[0] + a2 * P[1 + n];
    if (n == 1) r += a2 * a2 * P[4];
    return r;
}
double rhs_a2_N3_a4OS(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    r = P[0] + a2 * P[1 + n];
    if (n == 2) r += a2 * a2 * P[4];
    return r;
}
double rhs_a2_N3_a4_noa4WTI(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    r = P[0] + a2 * P[1 + n];
    if (n == 1) r += a2 * a2 * P[4];
    if (n == 2) r += a2 * a2 * P[5];
    return r;
}
double rhs_a2_N3_a4_noa4tm(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    r = P[0] + a2 * P[1 + n];
    if (n == 0) r += a2 * a2 * P[4];
    if (n == 2) r += a2 * a2 * P[5];
    return r;
}
double rhs_a2_N3_a4_noa4OS(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    r = P[0] + a2 * P[1 + n];
    if (n == 0) r += a2 * a2 * P[4];
    if (n == 1) r += a2 * a2 * P[5];
    return r;
}


template <double C>
double rhs_a2_Husung(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    static constexpr double lam = 295 / hbarc;
    static constexpr double lam2 = lam * lam;


    r = P[0] + a2 * P[1] + a2 / std::pow(-std::log(a2 * lam2), C) * P[2];

    return r;
}

template <double C>
double rhs_a2_N3_Husung(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    static constexpr double lam = 295 / hbarc;
    static constexpr double lam2 = lam * lam;


    r = P[0] + a2 * P[1 + n] + a2 / std::pow(-std::log(a2 * lam2), C) * P[1 + 3 + n];

    return r;
}
template <double C>
double rhs_a2_N3_HusungWTI(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    static constexpr double lam = 295 / hbarc;
    static constexpr double lam2 = lam * lam;

    r = P[0] + a2 * P[1 + n];
    if (n == 0) r += a2 / std::pow(-std::log(a2 * lam2), C) * P[4];
    return r;
}
template <double C>
double rhs_a2_N3_Husungtm(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    static constexpr double lam = 295 / hbarc;
    static constexpr double lam2 = lam * lam;
    r = P[0] + a2 * P[1 + n];
    if (n == 1) r += a2 / std::pow(-std::log(a2 * lam2), C) * P[4];
    return r;
}
template <double C>
double rhs_a2_N3_HusungOS(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a2 = x[0];
    static constexpr double lam = 295 / hbarc;
    static constexpr double lam2 = lam * lam;
    r = P[0] + a2 * P[1 + n];
    if (n == 2) r += a2 / std::pow(-std::log(a2 * lam2), C) * P[4];
    return r;
}

void shuffle_jackknife(double* jack_samples, int Njack) {
    // 1. Seed the random number generator (typically once at program start)
    // srand(time(NULL)); 

    if (Njack <= 1) return;

    // start from Njack - 2 because the last one is the mean and we want to keep it at the end of the array
    for (int i = Njack - 2; i > 0; i--) {
        // 2. Generate a random index j such that 0 <= j <= i
        int j = rand() % (i + 1);

        // 3. Swap components at indices i and j
        double temp = jack_samples[i];
        jack_samples[i] = jack_samples[j];
        jack_samples[j] = temp;
    }
}

int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir \n argc=%d\n", argc);
    char namefile[NAMESIZE];

    std::vector<std::string> files;
    // std::vector<std::string> beta_names;

    files.emplace_back(std::string(argv[2]) + "jack_scale_setting_system_B64");
    files.emplace_back(std::string(argv[2]) + "jack_scale_setting_system_C80");
    files.emplace_back(std::string(argv[2]) + "jack_scale_setting_system_D96");
    files.emplace_back(std::string(argv[2]) + "jack_scale_setting_system_E112");

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
    // add Zp/Zs
    //////////////////////////////////////////////////////////////
    jackall.add_space_for_n_observables(2);
    double* Zps_B64 = myres->create_fake(0.79007090, 0.00029391, -1);
    double* Zps_C80 = myres->create_fake(0.82327801, 0.00019334, -1);
    double* Zps_D96 = myres->create_fake(0.85131058, 0.00011303, -1);
    double* Zps_E112 = myres->create_fake(0.87068160, 0.00012266, -1);

    double* Zs_B64 = myres->create_fake(0.6549, 0.0037, -1);
    double* Zs_C80 = myres->create_fake(0.6427, 0.0032, -1);
    double* Zs_D96 = myres->create_fake(0.6291, 0.0025, -1);
    double* Zs_E112 = myres->create_fake(0.6180, 0.0021, -1);

    // double* Zps2_B64 = myres->create_fake(0.75552769, 0.002385089, -1);
    // double* Zps2_C80 = myres->create_fake(0.79427359, 0.00185236, -1);
    // double* Zps2_D96 = myres->create_fake(0.83269926, 0.00145200, -1);
    // double* Zps2_E112 = myres->create_fake(0.85859601, 0.00088635, -1);

    double* Zp_B64 = myres->create_fake(0.4852, 0.0036, -1);
    double* Zp_C80 = myres->create_fake(0.4920, 0.0022, -1);
    double* Zp_D96 = myres->create_fake(0.5016, 0.0017, -1);
    double* Zp_E112 = myres->create_fake(0.5049, 0.0015, -1);


    int idZ1 = jackall.en[0].Nobs - 1;
    int idZ2 = jackall.en[0].Nobs - 2;
    for (int e = 0;e < jackall.ens;e++) {
        switch (e) {
        case 0:
            jackall.en[e].jack[idZ1] = Zp_B64;
            jackall.en[e].jack[idZ2] = Zps_B64;
            myres->mult(jackall.en[e].jack[idZ2], jackall.en[e].jack[idZ2], Zs_B64);
            /* code */
            break;

        case 1:
            jackall.en[e].jack[idZ1] = Zp_C80;
            jackall.en[e].jack[idZ2] = Zps_C80;
            myres->mult(jackall.en[e].jack[idZ2], jackall.en[e].jack[idZ2], Zs_C80);
            /* code */
            break;

        case 2:
            jackall.en[e].jack[idZ1] = Zp_D96;
            jackall.en[e].jack[idZ2] = Zps_D96;
            myres->mult(jackall.en[e].jack[idZ2], jackall.en[e].jack[idZ2], Zs_D96);
            /* code */
            break;

        case 3:
            jackall.en[e].jack[idZ1] = Zp_E112;
            jackall.en[e].jack[idZ2] = Zps_E112;
            myres->mult(jackall.en[e].jack[idZ2], jackall.en[e].jack[idZ2], Zs_E112);
            /* code */
            break;

        default:
            break;
        }

    }

    double convert_rip_to_21 = 1.16851551284137;
    double from_21_to_4 = 0.826399813665257; // divide the mass by this to convert from the 21 to the 4 definition
    double from_21_to_9 = 0.91708855448905;

    for (int e = 0;e < jackall.ens;e++) {
        for (int j = 0; j < Njack;j++) {
            jackall.en[e].jack[idZ1][j] *= convert_rip_to_21;
            jackall.en[e].jack[idZ2][j] *= convert_rip_to_21;
        }

    }

    data_all jack_interpol;
    jack_interpol.create(2, 2, Njack, argv[1]);
    for (int j = 0;j < Njack;j++) {
        jack_interpol.en[0].jack[0][j] = fpi_MeV;
        jack_interpol.en[1].jack[0][j] = fpi_MeV_interpol;
    }

    //////////////////////////////////////////////////////////////
    // make table results
    //////////////////////////////////////////////////////////////
    mysprintf(namefile, NAMESIZE, "%s/data_from_fpi.txt", argv[3]);
    FILE* summary_out = open_file(namefile, "w+");
    fprintf(summary_out, "ens   a[fm] da[fm]  w0/a dw0/a  amul damul amus damus  amuc  damuc   delta_amul  ddelta_amul  delta_amus  ddelta_amus  delta_amuc   ddelta_amuc\n");
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

        tmp = jackall.en[i].jack[40];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[41];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[42];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

        fprintf(summary_out, "\n");
    }
    fclose(summary_out);
    mysprintf(namefile, NAMESIZE, "%s/data_from_w0.txt", argv[3]);
    summary_out = open_file(namefile, "w+");
    fprintf(summary_out, "ens   a[fm] da[fm]  afpi dafpi  amul damul amus damus  amuc  damuc    delta_amul  ddelta_amul  delta_amus  ddelta_amus  delta_amuc   ddelta_amuc\n");
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

        tmp = jackall.en[i].jack[43];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[44];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[45];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

        fprintf(summary_out, "\n");
    }
    fclose(summary_out);


    mysprintf(namefile, NAMESIZE, "%s/data_from_w0_hybrid.txt", argv[3]);
    summary_out = open_file(namefile, "w+");
    fprintf(summary_out, "ens   a[fm] da[fm]  afpi dafpi  amul damul amus damus  amuc  damuc    delta_amul  ddelta_amul  delta_amus  ddelta_amus  delta_amuc   ddelta_amuc\n");
    for (int i = 0;i < myen[0].size();i++) {
        size_t lastUnderscore = files[i].find_last_of('_');

        error(lastUnderscore == std::string::npos, 1, "error in file name", "cannot read ens from name %s", files[i].c_str());
        // 2. Estrae tutto ciò che segue l'ultimo underscore
        std::string result = files[i].substr(lastUnderscore + 1);
        fprintf(summary_out, "%s  ", result.c_str());
        double* tmp = jackall.en[i].jack[51];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[52];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[48];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[49];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[50];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

        tmp = jackall.en[i].jack[53];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[54];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[55];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

        fprintf(summary_out, "\n");
    }

    mysprintf(namefile, NAMESIZE, "%s/data_from_wp25_deriv_mc_a2.txt", argv[3]);
    summary_out = open_file(namefile, "w+");
    fprintf(summary_out, "ens   a[fm] da[fm]  afpi dafpi  amul damul amus damus  amuc  damuc    delta_amul  ddelta_amul  delta_amus  ddelta_amus  delta_amuc   ddelta_amuc\n");
    for (int i = 0;i < myen[0].size();i++) {
        size_t lastUnderscore = files[i].find_last_of('_');

        error(lastUnderscore == std::string::npos, 1, "error in file name", "cannot read ens from name %s", files[i].c_str());
        // 2. Estrae tutto ciò che segue l'ultimo underscore
        std::string result = files[i].substr(lastUnderscore + 1);
        fprintf(summary_out, "%s  ", result.c_str());
        double* tmp = jackall.en[i].jack[61];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[62];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[58];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[59];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[60];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

        tmp = jackall.en[i].jack[63];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[64];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[65];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

        fprintf(summary_out, "\n");
    }

    mysprintf(namefile, NAMESIZE, "%s/data_from_wp25_deriv_mc_a2_la_RDs.txt", argv[3]);
    summary_out = open_file(namefile, "w+");
    fprintf(summary_out, "ens   a[fm] da[fm]  afpi dafpi  amul damul amus damus  amuc  damuc    delta_amul  ddelta_amul  delta_amus  ddelta_amus  delta_amuc   ddelta_amuc\n");
    for (int i = 0;i < myen[0].size();i++) {
        size_t lastUnderscore = files[i].find_last_of('_');

        error(lastUnderscore == std::string::npos, 1, "error in file name", "cannot read ens from name %s", files[i].c_str());
        // 2. Estrae tutto ciò che segue l'ultimo underscore
        std::string result = files[i].substr(lastUnderscore + 1);
        fprintf(summary_out, "%s  ", result.c_str());
        double* tmp = jackall.en[i].jack[69];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[70];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[66];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[67];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[68];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

        tmp = jackall.en[i].jack[71];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[72];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[73];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

        fprintf(summary_out, "\n");
    }

    mysprintf(namefile, NAMESIZE, "%s/data_from_wp25_deriv_mc_a2_la_mc.txt", argv[3]);
    summary_out = open_file(namefile, "w+");
    fprintf(summary_out, "ens   a[fm] da[fm]  afpi dafpi  amul damul amus damus  amuc  damuc    delta_amul  ddelta_amul  delta_amus  ddelta_amus  delta_amuc   ddelta_amuc\n");
    for (int i = 0;i < myen[0].size();i++) {
        size_t lastUnderscore = files[i].find_last_of('_');

        error(lastUnderscore == std::string::npos, 1, "error in file name", "cannot read ens from name %s", files[i].c_str());
        // 2. Estrae tutto ciò che segue l'ultimo underscore
        std::string result = files[i].substr(lastUnderscore + 1);
        fprintf(summary_out, "%s  ", result.c_str());
        double* tmp = jackall.en[i].jack[86];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[87];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[83];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[84];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[85];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

        tmp = jackall.en[i].jack[88];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[89];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
        tmp = jackall.en[i].jack[90];
        fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

        fprintf(summary_out, "\n");
    }


    for (std::size_t ic = 0; ic < coeffs.size(); ++ic) {

        mysprintf(namefile, NAMESIZE, "%s/data_from_wp25_lin_deriv_mc_la_MDs_C%g.txt", argv[3], coeffs[ic]);
        summary_out = open_file(namefile, "w+");
        fprintf(summary_out, "ens   a[fm] da[fm]  afpi dafpi  amul damul amus damus  amuc  damuc    delta_amul  ddelta_amul  delta_amus  ddelta_amus  delta_amuc   ddelta_amuc\n");
        for (int i = 0;i < myen[0].size();i++) {
            size_t lastUnderscore = files[i].find_last_of('_');

            error(lastUnderscore == std::string::npos, 1, "error in file name", "cannot read ens from name %s", files[i].c_str());
            // 2. Estrae tutto ciò che segue l'ultimo underscore
            std::string result = files[i].substr(lastUnderscore + 1);
            fprintf(summary_out, "%s  ", result.c_str());
            double* tmp = jackall.en[i].jack[95 + ic * 8];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[96 + ic * 8];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[92 + ic * 8];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[93 + ic * 8];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[94 + ic * 8];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

            tmp = jackall.en[i].jack[97 + ic * 8];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[98 + ic * 8];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[99 + ic * 8];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

            fprintf(summary_out, "\n");
        }
    }

    for (std::size_t ic = 0; ic < coeffs_mc.size(); ++ic) {

        mysprintf(namefile, NAMESIZE, "%s/data_from_wp25_deriv_mc_a2_lamcC%g.txt", argv[3], coeffs_mc[ic]);
        summary_out = open_file(namefile, "w+");
        fprintf(summary_out, "ens   a[fm] da[fm]  afpi dafpi  amul damul amus damus  amuc  damuc    delta_amul  ddelta_amul  delta_amus  ddelta_amus  delta_amuc   ddelta_amuc\n");
        for (int i = 0;i < myen[0].size();i++) {
            size_t lastUnderscore = files[i].find_last_of('_');

            error(lastUnderscore == std::string::npos, 1, "error in file name", "cannot read ens from name %s", files[i].c_str());
            // 2. Estrae tutto ciò che segue l'ultimo underscore
            std::string result = files[i].substr(lastUnderscore + 1);
            fprintf(summary_out, "%s  ", result.c_str());
            double* tmp = jackall.en[i].jack[95 + coeffs.size() * 8 + ic * 9];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[96 + coeffs.size() * 8 + ic * 9];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[92 + coeffs.size() * 8 + ic * 9];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[93 + coeffs.size() * 8 + ic * 9];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[94 + coeffs.size() * 8 + ic * 9];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

            tmp = jackall.en[i].jack[97 + coeffs.size() * 8 + ic * 9];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[98 + coeffs.size() * 8 + ic * 9];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));
            tmp = jackall.en[i].jack[99 + coeffs.size() * 8 + ic * 9];
            fprintf(summary_out, "%.12g    %.12g     ", myres->mean(tmp), myres->comp_error(tmp));

            fprintf(summary_out, "\n");
        }
    }

    //////////////////////////////////////////////////////////////
    // fitting
    //////////////////////////////////////////////////////////////
    std::vector<std::string> obs = { "w0" ,"fpi", "w0_ens", "w0_hybrid", "fpi_hybrid", "fpi_Mpi_wp25_hybrid", "fpi_wp25_lin", "fpi_wp25_lin_laRDs",  "fpi_wp25_lin_laRDs_exact",
         "MDs_wp25_lin_laRDs_exact", "fpi_wp25_lin_la_mc_exact", "MDs_wp25_lin_la_mc_exact" ,
        "w0_lin_deriv" ,"sqrtt0_from_fpi" ,"w0_fpi_wp25_Cm5" , "w0_fpi_interpol" };

    std::vector<int> id_obs = { 34, 39, 46 ,47, 52, 57, 62, 70, 78 ,
        82, 87, 91, id_w0_lin_deriv, id_sqrtt0_from_fpi, id_w0_fpi_wp25_Cm5, id_w0_fpi_interpol
    };
    std::vector<int> id_a = { 33, 38 , 33, 33 ,51, 56, 61, 69, 77 ,
         77, 86, 86,33, 33, id_a_fpi_wp25_Cm5, id_a_fpi_interpol
    };

    for (std::size_t ic = 0; ic < coeffs.size(); ++ic) {
        const double coeff = coeffs[ic];
        obs.push_back("fpi_wp25_C" + std::to_string(coeff));
        id_obs.push_back(96 + ic * 8);
        id_a.push_back(95 + ic * 8);
    }

    for (std::size_t ic = 0; ic < coeffs_mc.size(); ++ic) {
        const double coeff = coeffs_mc[ic];
        obs.push_back("fpi_wp25_mc_C" + std::to_string(coeff));
        id_obs.push_back(96 + coeffs.size() * 8 + ic * 9);
        id_a.push_back(95 + coeffs.size() * 8 + ic * 9);
    }
    for (std::size_t ic = 0; ic < coeffs_mc.size(); ++ic) {
        const double coeff = coeffs_mc[ic];
        obs.push_back("MDs_wp25_mc_C" + std::to_string(coeff));
        id_obs.push_back(100 + coeffs.size() * 8 + ic * 9);
        id_a.push_back(95 + coeffs.size() * 8 + ic * 9);
    }
    for (std::size_t ic = 0; ic < coeffs.size(); ++ic) {
        const double coeff = coeffs[ic];
        obs.push_back("sqrtt0_wp25_C" + std::to_string(coeff));
        id_obs.push_back(id_sqrtt0_from_fpi + 1 + ic);
        id_a.push_back(95 + ic * 8);
    }

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

        if (obs[i] == "w0_lin_deriv") {
            myres->copy(jack_interpol.en[0].jack[1], der_fpi_const_full.P[0]);
        }
        if (obs[i] == "w0_fpi_interpol") {
            myres->copy(jack_interpol.en[1].jack[1], der_fpi_const_full.P[0]);
        }


        der_fpi_const_full.clear();
        // if (i < 2)
        //     for (size_t imr = 0; imr < obs2.size(); imr++) {
        //         fit_info.corr_id = id_obs_m[imr + obs2.size() * i];
        //         std::string namefit = "fit_" + obs2[imr] + "_from_" + obs[(i + 1) % 2] + "_cont_a2";
        //         fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun_ratio, fit_info, namefit.c_str());
        //         fit_info.band_range = { 0, 0.008145209846823482 };
        //         print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});

        //         der_fpi_const_full.clear();
        //     }
        fit_info.restore_default();

    }


    //////////////////////////////////////////////////////////////
    // obs2
    //////////////////////////////////////////////////////////////
    std::vector<std::string> from = { "fpi" ,"w0", "w0_hybrid" };
    std::vector<std::string> obs2 = { "ms_over_ml" ,"mc_over_ms" };
    std::vector<std::vector<int>> id_obs_m = {
        {31,30}, {32,31},  // from fpi
        {36,35}, {37,36},  // from w0
        {49,48}, {50,49}   // from w0_hybrid
    };
    id_a = { 33, 38  ,51 };
    for (int ifrom = 0; ifrom < from.size(); ifrom++) {
        fit_type fit_info;

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
                    fit_info.x[0][count][j] = pow(jackall.en[e].jack[id_a[ifrom]][j], 2); // a^2 fm^2
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



        for (size_t iobs = 0; iobs < obs2.size(); iobs++) {
            fit_info.corr_id = id_obs_m[iobs + obs2.size() * ifrom];
            std::string namefit = "fit_" + obs2[iobs] + "_from_" + from[ifrom] + "_cont_a2";
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
    std::vector<std::string> diff_a = { "FLAG_wp25", "FLAG_wp25hybrid", "FLAG_wp25Mpihybrid" };
    std::vector<std::vector<int>> id_diff_a = {
        { 33, 38 }, // w0
        { 33, 51 },  // w0_hybrid
        { 33, 56 }  // w0_Mpi_hybrid
    };
    for (int idiff_a = 0; idiff_a < diff_a.size(); idiff_a++) {

        fit_type fit_info;

        fit_info.corr_id = { id_diff_a[idiff_a][0], id_diff_a[idiff_a][1] };
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
                    fit_info.x[0][count][j] = pow(jackall.en[e].jack[id_diff_a[idiff_a][0]][j], 2); // a^2 fm^2
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

        std::string namefit = "fit_a_ratio_" + diff_a[idiff_a] + "_cont_a2";

        fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun_ratio, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});

        der_fpi_const_full.clear();
        namefit = "fit_a_diff_" + diff_a[idiff_a] + "_cont_a2";

        der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun_diff, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});

        fit_info.restore_default();


    }

    //////////////////////////////////////////////////////////////
    // quark masses
    //////////////////////////////////////////////////////////////
    std::vector<std::vector<int>> id_masses = { {30,31,32}, {35,36,37}, {48,49,50} };

    for (int ifrom = 0; ifrom < from.size(); ifrom++) {
        fit_type fit_info;

        fit_info.N = 3 * 2;
        fit_info.Nxen = std::vector<std::vector<int>>(fit_info.N);
        for (int n = 0; n < fit_info.N; n++) {
            fit_info.Nxen[n].resize(myen[0].size());
            for (int b = 0; b < myen[0].size(); b++) {
                fit_info.Nxen[n][b] = myen[0][b];
            }
        }
        fit_info.init_N_etot_form_Nxen();

        fit_info.function = rhs_masses;
        fit_info.linear_fit = false;
        int Nz = 2;
        int Nferm = 3;
        fit_info.Npar = Nferm + Nz * Nferm + Nz;
        fit_info.Nvar = 1; // a2
        fit_info.Njack = jackall.en[0].Njack;
        fit_info.malloc_x();
        // double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);

        int count = 0;
        for (int n = 0; n < fit_info.N; n++) {
            for (int e : fit_info.Nxen[n]) {
                for (int j = 0; j < Njack; j++) {
                    fit_info.x[0][count][j] = pow(jackall.en[e].jack[id_a[ifrom]][j], 2); // a^2 fm^2
                }
                count++;
            }
        }
        fit_info.corr_id = {
            id_masses[ifrom][0], idZ1,
            id_masses[ifrom][0], idZ2,
            id_masses[ifrom][1], idZ1,
            id_masses[ifrom][1], idZ2,
            id_masses[ifrom][2], idZ1,
            id_masses[ifrom][2], idZ2
        };

        // fit_info.linear_fit = false;
        fit_info.verbosity = 0;
        // fit_info.covariancey = true;
        // fit_info.compute_cov_fit(argv, jackall, lhs_fun);
        // fit_info.make_covariance_block_diagonal_in_n();
        // fit_info.compute_cov1_fit();


        std::string namefit = "fit_masses_from_" + from[ifrom] + "_cont_a2";
        fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun_masses, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});

        der_fpi_const_full.clear();

        // fit_info.restore_default();

        Nferm = 1;
        fit_info.N = 1 * 2;
        fit_info.Npar = Nferm + Nz * Nferm + Nz;
        fit_info.Nxen = std::vector<std::vector<int>>(fit_info.N);
        for (int n = 0; n < fit_info.N; n++) {
            fit_info.Nxen[n].resize(myen[0].size());
            for (int b = 0; b < myen[0].size(); b++) {
                fit_info.Nxen[n][b] = myen[0][b];
            }
        }
        fit_info.init_N_etot_form_Nxen();

        fit_info.function = rhs_mass_light;

        namefit = "fit_mass_light_from_" + from[ifrom] + "_cont_a2";
        fit_result m_light_fit = fit_all_data(argv, jackall, lhs_fun_masses, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", m_light_fit, m_light_fit, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});


    }


    //////////////////////////////////////////////////////////////
    // AIC fpi
    //////////////////////////////////////////////////////////////

    obs = { "fpi_wp25_lin" };
    id_obs = { 62 };
    id_a = { 61 };
    std::vector<std::string> fits = { "a2", "a2_3pt", "a2_a4" , "a2_3pt_noC" };
    double (*funcArray[])(int, int, double*, int, double*) = { rhs_a2, rhs_a2, rhs_a2_a4 ,rhs_a2 };
    std::vector<int> Npars = { 2, 2, 3, 2 };
    std::vector<std::vector<int>> ens_to_fit = { {0,1,2,3}, {1,2,3}, {0,1,2,3}, {0,2,3} };

    for (auto [ifit, fit] : std::views::enumerate(fits)) {
        for (int i = 0; i < obs.size(); i++) {
            fit_type fit_info;

            fit_info.corr_id = { id_obs[i] };

            //////////////////////////////////////////////////////////////
            // charm small volume fit - constant
            //////////////////////////////////////////////////////////////

            fit_info.N = 1;
            fit_info.Nxen = std::vector<std::vector<int>>(fit_info.N);
            for (int n = 0; n < fit_info.N; n++) {
                // fit_info.Nxen[n].resize(ens_to_fit[ifit].size());
                // for (int b = 0; b < ens_to_fit[ifit].size(); b++) {
                //     fit_info.Nxen[n][b] = ens_to_fit[ifit][b];  // the charm is in myen[beta][0]
                // }
                fit_info.Nxen[n] = ens_to_fit[ifit];
            }
            fit_info.init_N_etot_form_Nxen();

            fit_info.function = funcArray[ifit];
            fit_info.linear_fit = true;
            fit_info.Npar = Npars[ifit];
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

            std::string namefit = "fit_" + obs[i] + "_" + fit;

            fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
            fit_info.band_range = { 0, 0.008145209846823482 };
            print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});
            der_fpi_const_full.clear();

            fit_info.restore_default();

        }
    }

    //////////////////////////////////////////////////////////////
    // multifit fpi
    //////////////////////////////////////////////////////////////
    std::vector<std::vector<std::vector<int>>> points = {
         {{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },
         {{1,2,3}, {1,2,3}, {1,2,3} },
         {{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },
         {{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },
         {{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },
         {{0,2,3}, {0,2,3}, {0,2,3} },
         {{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },
         {{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },
         {{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },
         {{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },{{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },{{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },
         {{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },{{0,1,2,3}, {0,1,2,3}, {0,1,2,3} },{{0,1,2,3}, {0,1,2,3}, {0,1,2,3} }
    };
    std::vector<std::string> fits_fpi = { "N3_a2", "N3_a2_noBOS_noBtm_noBWTI" ,
        "N3_a2_a4",
        "N3_a2_Husung0.42", "N3_a2_Husung0.21",
        "N3_a2_noCOS_noCtm_noCWTI",
        "N3_a2_a4WTI",    "N3_a2_a4tm",    "N3_a2_a4OS",
        "N3_a2_a4_noa4WTI","N3_a2_a4_noa4tm","N3_a2_a4_noa4OS",
        "N3_a2_HusungWTI", "N3_a2_Husungtm", "N3_a2_HusungOS"
    };
    // double (*funcArray_fpi[])(int, int, double*, int, double*) = { rhs_a2_N3, rhs_a2_N3, rhs_a2_N3_a4 , rhs_a2_N3_Husung<0.42>, rhs_a2_N3_Husung<0.21>,
    //  rhs_a2_N3,
    // rhs_a2_N3_a4WTI,
    // rhs_a2_N3_a4tm,
    // rhs_a2_N3_a4OS
    // };
    std::vector<double (*)(int, int, double*, int, double*)> funcArray_fpi = { rhs_a2_N3,    rhs_a2_N3,
    rhs_a2_N3_a4,
    rhs_a2_N3_Husung<0.42>,    rhs_a2_N3_Husung<0.21>,
    rhs_a2_N3,
    rhs_a2_N3_a4WTI,    rhs_a2_N3_a4tm,    rhs_a2_N3_a4OS,
    rhs_a2_N3_a4_noa4WTI, rhs_a2_N3_a4_noa4tm, rhs_a2_N3_a4_noa4OS,
    rhs_a2_N3_HusungWTI<0.42>, rhs_a2_N3_Husungtm<0.42>, rhs_a2_N3_HusungOS<0.42>
    };
    std::vector<int> Npars_fpi = { 4, 4,
        7,
        7 ,7,
        4,
        5,5,5,
        6,6,6,
        5,5,5
    };

    error(fits_fpi.size() != Npars_fpi.size(), 1, "main", "fits_fpi  %zu  and Npars_fpi  %zu  must have the same size", fits_fpi.size(), Npars_fpi.size());
    error(fits_fpi.size() != points.size(), 1, "main", "fits_fpi  %zu  and points  %zu  must have the same size", fits_fpi.size(), points.size());
    error(fits_fpi.size() != funcArray_fpi.size(), 1, "main", "fits_fpi  %zu  and funcArray_fpi  %zu  must have the same size", fits_fpi.size(), funcArray_fpi.size());

    for (int ic = 0;ic < coeffs.size();ic++) {
        int id_a_of_ic = 95 + ic * 8;
        for (auto [ifit, fit] : std::views::enumerate(fits_fpi)) {

            fit_type fit_info;
            fit_info.corr_id = {
                96 + ic * 8,
                sid_fpi_A0 + ic * 2,
                id_fpi_OS_dWTI + ic/* sid_fpi_A0 + 0 * 2 + 1  */ };
            fit_info.N = 3;
            fit_info.Nxen = std::vector<std::vector<int>>(fit_info.N);
            for (int n = 0; n < fit_info.N; n++) {
                // fit_info.Nxen[n].resize(ens_to_fit[ifit].size());
                // for (int b = 0; b < ens_to_fit[ifit].size(); b++) {
                //     fit_info.Nxen[n][b] = ens_to_fit[ifit][b];  // the charm is in myen[beta][0]
                // }
                fit_info.Nxen[n] = points[ifit][n];
            }
            fit_info.init_N_etot_form_Nxen();
            fit_info.function = funcArray_fpi[ifit];
            fit_info.linear_fit = true;
            fit_info.Npar = Npars_fpi[ifit];
            fit_info.Nvar = 1; // a2
            fit_info.Njack = jackall.en[0].Njack;
            fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);

            int count = 0;
            for (int n = 0; n < fit_info.N; n++) {
                for (int e : fit_info.Nxen[n]) {
                    for (int j = 0; j < Njack; j++) {
                        fit_info.x[0][count][j] = pow(jackall.en[e].jack[id_a_of_ic][j], 2); // a^2 fm^2
                    }
                    count++;
                }
            }
            // fit_info.linear_fit = false;
            fit_info.verbosity = 0;
            fit_info.covariancey = true;
            fit_info.compute_cov_fit(argv, jackall, lhs_fun_N3);
            // fit_info.make_covariance_block_diagonal_in_n();
            fit_info.make_covariance_block_diagonal_in_e();
            fit_info.compute_cov1_fit();

            std::string namefit = "fit_fpi_C" + std::to_string(coeffs[ic]) + "_" + fit;

            fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun_N3, fit_info, namefit.c_str());
            fit_info.band_range = { 0, 0.008145209846823482 };
            print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});
            der_fpi_const_full.clear();

            fit_info.restore_default();

        }
    }

    //////////////////////////////////////////////////////////////
    // interpolation w0
    //////////////////////////////////////////////////////////////
    double* fpi_MeV_wp25_Cm5_j = myres->create_fake(fpi_MeV_wp25_Cm5, fpi_MeV_err_wp25_Cm5, 1234);
    double* w0_interpolated_j = myres->create_zero();
    {
        fit_type fit_info;
        fit_info.corr_id = { 1 };
        fit_info.N = 1;
        fit_info.Nxen = { {0,1} };
        fit_info.init_N_etot_form_Nxen();
        fit_info.function = rhs_a2;
        fit_info.linear_fit = true;
        fit_info.Npar = 2;
        fit_info.Nvar = 1; // fpi
        fit_info.Njack = jack_interpol.en[0].Njack;
        fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
        int count = 0;
        for (int n = 0; n < fit_info.N; n++) {
            for (int e : fit_info.Nxen[n]) {
                for (int j = 0; j < Njack; j++) {
                    fit_info.x[0][count][j] = jack_interpol.en[e].jack[0][j]; // fpi MeV
                }
                count++;
            }
        }
        fit_info.verbosity = 0;
        std::string namefit = "fit_w0_interpolation_fpi";
        fit_result der_fpi_const_full = fit_all_data(argv, jack_interpol, lhs_identity, fit_info, namefit.c_str());
        fit_info.band_range = { jack_interpol.en[0].jack[0][Njack - 1] - 0.1,jack_interpol.en[1].jack[0][Njack - 1] + 0.1 };
        print_fit_band(argv, jack_interpol, fit_info, fit_info, namefit.c_str(), "fpi", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.1, {});

        der_fpi_const_full.fit_to_tif();
        std::vector<double> tmpx(1);
        for (int j = 0; j < Njack;j++) {
            tmpx[0] = fpi_MeV_wp25_Cm5_j[j]; // fpi MeV
            w0_interpolated_j[j] = fit_info.function(0, fit_info.Nvar, tmpx.data(), fit_info.Npar, der_fpi_const_full.tif[j]);
        }

        // shuffling test
        double* w0_tmp = myres->create_zero();
        printf("Nshuffle  w0   dw0   fp   dfpi\n");
        for (int j = 0; j < Njack;j++) {
            tmpx[0] = fpi_MeV_wp25_Cm5_j[j]; // fpi MeV
            w0_tmp[j] = fit_info.function(0, fit_info.Nvar, tmpx.data(), fit_info.Npar, der_fpi_const_full.tif[j]);
        }
        printf("%-5d  %-17.12g   %-17.12g   %-17.12g   %-17.12g\n", 0, myres->mean(w0_tmp), myres->comp_error(w0_tmp), myres->mean(fpi_MeV_wp25_Cm5_j), myres->comp_error(fpi_MeV_wp25_Cm5_j));
        for (int i = 0; i < 100;i++) {
            shuffle_jackknife(fpi_MeV_wp25_Cm5_j, Njack);
            for (int j = 0; j < Njack;j++) {
                tmpx[0] = fpi_MeV_wp25_Cm5_j[j]; // fpi MeV
                w0_tmp[j] = fit_info.function(0, fit_info.Nvar, tmpx.data(), fit_info.Npar, der_fpi_const_full.tif[j]);
            }
            printf("%-5d  %-17.12g   %-17.12g   %-17.12g   %-17.12g\n", i + 1, myres->mean(w0_tmp), myres->comp_error(w0_tmp), myres->mean(fpi_MeV_wp25_Cm5_j), myres->comp_error(fpi_MeV_wp25_Cm5_j));
        }
        printf("seed  w0   dw0   fp   dfpi\n");
        for (int i = 0; i < 100;i++) {
            free(fpi_MeV_wp25_Cm5_j);
            fpi_MeV_wp25_Cm5_j = myres->create_fake(fpi_MeV_wp25_Cm5, fpi_MeV_err_wp25_Cm5, i + 1);
            for (int j = 0; j < Njack;j++) {
                tmpx[0] = fpi_MeV_wp25_Cm5_j[j]; // fpi MeV
                w0_tmp[j] = fit_info.function(0, fit_info.Nvar, tmpx.data(), fit_info.Npar, der_fpi_const_full.tif[j]);
            }
            printf("%-5d  %-17.12g   %-17.12g   %-17.12g   %-17.12g\n", i + 1, myres->mean(w0_tmp), myres->comp_error(w0_tmp), myres->mean(fpi_MeV_wp25_Cm5_j), myres->comp_error(fpi_MeV_wp25_Cm5_j));
        }

        der_fpi_const_full.clear();

        // interpolate at fpi_MeV_interpol_j
    }

    printf("w0_fpi_interpol    %.12g    %.12g\n", myres->mean(w0_interpolated_j), myres->comp_error(w0_interpolated_j));

    FILE* fout = fopen("test_seed_and_Njack.txt", "w");
    fprintf(fout, "seed  Njack  mean   err\n");
    for (int n : {20, 100, 1000})
        for (int i = 0;i < 100;i++) {
            resampling_jack myresampling(n);
            double* tmp;
            tmp = myresampling.create_fake(1, 0.1, i);
            fprintf(fout, "%-5d  %-5d %.12g   %.12g\n", i, n, myresampling.mean(tmp), myresampling.comp_error(tmp));
            free(tmp);
        }
    fclose(fout);

    // auto it = std::find(coeffs.begin(), coeffs.end(), -5.0);
    auto it = std::ranges::find(coeffs, -5.0);
    printf("Finding -5.0 in coeffs: %s   %g\n", (it != coeffs.end()) ? "Found" : "Not Found", *it);
    std::size_t id_Cm5;
    if (it != coeffs.end()) {
        id_Cm5 = std::distance(coeffs.begin(), it);
    }
    else {
        printf("Error: -5.0 not found in coeffs\n");
        exit(1);
    }
    printf("Index of -5.0 in coeffs: %zu\n", id_Cm5);
    int id_sqrtt0_Cm5 = id_sqrtt0_from_fpi + 1 + id_Cm5;
    printf("Index of  strtt0 wp25 -5.0 in coeffs: %zu\n", id_Cm5);


    //////////////////////////////////////////////////////////////
    // fits sqrtt0
    //////////////////////////////////////////////////////////////
    //fits
    points = {
        {{0,1,2,3}},
        {{1,2,3}},
        {{0,2,3}},
        {{0,1,2,3}},
        {{0,1,2,3}},
        {{0,1,3}},
        {{0,1,2}}
    };
    std::vector<std::string> fits_sqrtt0 = { "a2", "a2_noB", "a2_noC", "a2_a4", "a2_Husung0.42" , "a2_noD", "a2_noE" };
    std::vector<double (*)(int, int, double*, int, double*)> funcArray_sqrtt0 = { rhs_a2,rhs_a2,rhs_a2, rhs_a2_a4 , rhs_a2_Husung<0.42>,rhs_a2,rhs_a2 };
    std::vector<int> Npars_sqrtt0 = { 2,2,2, 3, 3, 2, 2 };

    //obs to fit
    id_obs = { id_sqrtt0_from_fpi ,id_sqrtt0_Cm5 , id_w0_lin_deriv };
    obs = { "sqrtt0_FLAG", "sqrtt0_wp25_Cm5", "w0_lin_deriv_FLAG" };
    id_a = { 33 , 95 + static_cast<int>(id_Cm5) * 8, 33 };

    for (int i = 0; i < obs.size(); i++) {
        for (auto [ifit, fit] : std::views::enumerate(fits_sqrtt0)) {

            fit_type fit_info;

            fit_info.corr_id = { id_obs[i] };
            fit_info.Nxen = points[ifit];
            fit_info.init_N_etot_form_Nxen();
            fit_info.function = funcArray_sqrtt0[ifit];
            fit_info.linear_fit = true;
            fit_info.Npar = Npars_sqrtt0[ifit];
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
            std::string namefit = "fit_" + obs[i] + "_" + fit;
            fit_result res_sqrtt0 = fit_all_data(argv, jackall, lhs_fun_N3, fit_info, namefit.c_str());
            fit_info.band_range = { 0, 0.008145209846823482 };
            print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", res_sqrtt0, res_sqrtt0, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});
            res_sqrtt0.clear();

            fit_info.restore_default();


        }
    }
    //////////////////////////////////////////////////////////////
    // w0 sim * aiso
    //////////////////////////////////////////////////////////////
    {
        int id_a = 33;
        int id_w0_small_L = id_fpi_sim_big_L + 9;
        int id_w0_big_L = id_fpi_sim_big_L + 6;
        fit_type fit_info;

        fit_info.corr_id = { id_w0_small_L, id_w0_big_L , id_a}; // w0, fpi
        fit_info.Nxen = { {0,1,2,3}, {0,2} };
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
                    double a =  jackall.en[e].jack[id_a][j];
                    fit_info.x[0][count][j] = a*a; // a^2 fm^2
                }
                count++;
            }
        }
        std::string namefit = "fit_w0_sim_aiso_a2";
        fit_result res_sqrtt0 = fit_all_data(argv, jackall, lhs_w0_sim_aiso, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", res_sqrtt0, res_sqrtt0, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});
        res_sqrtt0.clear();

        fit_info.restore_default();

    }
    //////////////////////////////////////////////////////////////
    // w0 max twist
    //////////////////////////////////////////////////////////////
    {
        int id_fpi = 3;
        int id_w0 = 4;
        fit_type fit_info;

        fit_info.corr_id = { id_fpi_sim_big_L+1, id_fpi_sim_big_L+2 , id_fpi_sim_big_L+3}; // w0, fpi
        fit_info.Nxen = { {0,1,2,3}, {0,1} };
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
                    double fpi = (n == 0) ? jackall.en[e].jack[fit_info.corr_id[1]][j] : jackall.en[e].jack[fit_info.corr_id[2]][j];
                    fit_info.x[0][count][j] = pow(fpi / (fpi_MeV / hbarc), 2); // a^2 fm^2
                }
                count++;
            }
        }
        std::string namefit = "fit_w0_max_twist_a2";
        fit_result res_sqrtt0 = fit_all_data(argv, jackall, lhs_w0_sim, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", res_sqrtt0, res_sqrtt0, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});
        res_sqrtt0.clear();

        fit_info.restore_default();

    }
    //////////////////////////////////////////////////////////////
    // w0 Linf
    //////////////////////////////////////////////////////////////
    {
        int id_fpi = 3;
        int id_w0 = 4;
        fit_type fit_info;

        fit_info.corr_id = { id_fpi_sim_big_L+4, id_fpi_sim_big_L+5}; // w0, fpi
        fit_info.Nxen = { {0,1,2,3} };
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
                    double fpi = (n == 0) ? jackall.en[e].jack[fit_info.corr_id[1]][j] : jackall.en[e].jack[fit_info.corr_id[2]][j];
                    fit_info.x[0][count][j] = pow(fpi / (fpi_MeV / hbarc), 2); // a^2 fm^2
                }
                count++;
            }
        }
        std::string namefit = "fit_w0_Linf_a2";
        fit_result res_sqrtt0 = fit_all_data(argv, jackall, lhs_w0_sim, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", res_sqrtt0, res_sqrtt0, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});
        res_sqrtt0.clear();

        fit_info.restore_default();

    }


    // done with B96 and C112
    //////////////////////////////////////////////////////////////
    // w0 sim
    //////////////////////////////////////////////////////////////
    {
        int id_fpi_small_L = id_fpi_sim_big_L + 8;
        int id_w0_small_L = id_fpi_sim_big_L + 9;
        int id_w0_big_L = id_fpi_sim_big_L + 6;
        fit_type fit_info;

        fit_info.corr_id = { id_w0_small_L, id_fpi_small_L , id_fpi_sim_big_L, id_w0_big_L}; // w0, fpi
        fit_info.Nxen = { {0,1,2,3}, {0,1,2} };
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
                    double fpi = (n == 0) ? jackall.en[e].jack[id_fpi_small_L][j] : jackall.en[e].jack[id_fpi_sim_big_L][j];
                    fit_info.x[0][count][j] = pow(fpi / (fpi_MeV / hbarc), 2); // a^2 fm^2
                }
                count++;
            }
        }
        std::string namefit = "fit_w0_sim_L_a2";
        fit_result res_sqrtt0 = fit_all_data(argv, jackall, lhs_w0_sim1, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", res_sqrtt0, res_sqrtt0, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});
        res_sqrtt0.clear();

        fit_info.restore_default();

    }
    //////////////////////////////////////////////////////////////
    // w0 full
    //////////////////////////////////////////////////////////////
    {
        int id_a = 33;
        int id_w0 = id_w0_lin_deriv;
        fit_type fit_info;

        fit_info.corr_id = { id_w0, id_a }; // w0, fpi
        fit_info.Nxen = { {0,1,2,3} };
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
                    double a =  jackall.en[e].jack[id_a][j];
                    // fit_info.x[0][count][j] = pow(fpi / (fpi_MeV / hbarc), 2); // a^2 fm^2
                    fit_info.x[0][count][j] = a*a;
                }
                count++;
            }
        }
        std::string namefit = "fit_w0_full_a2";
        fit_result res_sqrtt0 = fit_all_data(argv, jackall, lhs_w0_full, fit_info, namefit.c_str());
        fit_info.band_range = { 0, 0.008145209846823482 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "a2", res_sqrtt0, res_sqrtt0, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.001, {});
        res_sqrtt0.clear();

        fit_info.restore_default();

    }
    
    printf("Njack: %d\n", Njack);
}
