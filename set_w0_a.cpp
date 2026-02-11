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
#include "fve.hpp"
#include "correlators_analysis.hpp"

#include "functions_w0.hpp"

#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

enum enum_ensembles {
    B72_64,
    B72_96,
    C06,
    D54,
    A53,
    A40,
    A30,
    E112,
    C112
};
constexpr double fpi_MeV = 130.5;
constexpr double fpi_MeV_err = 0.04;

constexpr double Mpi_MeV = 135;
constexpr double Mpi_MeV_err = 0.2;

constexpr double w0_fm = 0.17236;
constexpr double w0_fm_err = 0.0000002; // 0.00070

constexpr double w0_MeV = w0_fm / hbarc;
constexpr double w0_MeV_err = w0_fm_err / hbarc; // 0.00070

constexpr double MK_MeV = 494.6;
constexpr double MK_MeV_err = 0.3;

constexpr double MDs_MeV = 1967.0;
constexpr double MDs_MeV_err = 0.4;

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
    // if (j==fit_info.Njack-1) {
    //     printf("lhs_fun: n=%d e=%d j=%d fit_info.corr_id[0]=%d fit_info.corr_id[1]=%d fit_info.corr_id[2]=%d\n",
    //            n, e, j, fit_info.corr_id[0], fit_info.corr_id[1], fit_info.corr_id[2]);
    //     printf("lhs_fun: gjack.en[e].jack[fit_info.corr_id[0]][j]=%g\n", gjack.en[e].jack[fit_info.corr_id[0]][j]);
    //     printf("lhs_fun: gjack.en[e].jack[fit_info.corr_id[1]][j]=%g\n", gjack.en[e].jack[fit_info.corr_id[1]][j]);
    //     printf("lhs_fun: gjack.en[e].jack[fit_info.corr_id[2]][j]=%g\n", gjack.en[e].jack[fit_info.corr_id[2]][j]);
    // }
    double r = gjack.en[e].jack[fit_info.corr_id[0]][j]; // d(w0/a)/d(a*mu)
    r *= gjack.en[e].jack[fit_info.corr_id[1]][j];       // a*mu_l
    r *= gjack.en[e].jack[fit_info.corr_id[2]][j];       // a
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
struct file_out_name {
    char path[NAMESIZE];
    char basename[NAMESIZE];
    char namefile[NAMESIZE];
    char label[NAMESIZE];
    file_out_name(const char* path, const char* label) {
        mysprintf(this->path, NAMESIZE, "%s", path);
        mysprintf(this->basename, NAMESIZE, "%s", path);
        mysprintf(this->label, NAMESIZE, "%s", label);
        mysprintf(this->namefile, NAMESIZE, "%s/%s_fit_extra.txt", path, label);
    }
    file_out_name(const file_out_name& other) {
        mysprintf(this->path, NAMESIZE, "%s", other.path);
        mysprintf(this->basename, NAMESIZE, "%s", other.basename);
        mysprintf(this->label, NAMESIZE, "%s", other.label);
        mysprintf(this->namefile, NAMESIZE, "%s", other.namefile);
    }
    file_out_name& operator=(const file_out_name& other) {
        mysprintf(this->path, NAMESIZE, "%s", other.path);
        mysprintf(this->basename, NAMESIZE, "%s", other.basename);
        mysprintf(this->label, NAMESIZE, "%s", other.label);
        mysprintf(this->namefile, NAMESIZE, "%s", other.namefile);
        return *this;
    }

    ~file_out_name() {
        // free(this->namefile);
        // free(this->path);
        // free(this->basename);
    }
};

void read_file_debug(double* jack, const char* file) {
    printf("reading file %s\n", file);
    myres->read_jack_from_file(jack, file);
}
int id_deriv(int Meson, int deriv, int quark, int val_sea) {
    return Meson + deriv * (5 + 5 * (quark + 3 * val_sea));
}

int main(int argc, char** argv) {
    error(argc != 5, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir    file_inputs");
    char namefile[NAMESIZE];

    std::vector<std::string> files;
    std::ifstream file(argv[4]);
    if (!file.is_open()) {
        error(1, 1, "main", "Could not open file with input files: %s", argv[4]);
    }
    std::string line;
    std::string basename;
    int beta_count = 0;
    int file_count = 0;
    // std::vector<std::string> beta_names;
    std::vector<std::vector<int>> myen(1, std::vector<int>());

    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '#') // skip empty lines and comments
        {

            std::vector<std::string> word = split(line, ' ');

            if (word.size() == 1) {
                basename = word[0];
                mysprintf(namefile, NAMESIZE, "%s", word[0].c_str());
                // printf("adding file %s\n", namefile);
                files.emplace_back(namefile);
                myen[beta_count].emplace_back(file_count);
                file_count++;
            }
            else if (word.size() == 1) {
                if (strcmp(word[0].c_str(), "new_beta") == 0) {
                    myen.emplace_back(std::vector<int>());
                    beta_count++;
                }
                else {
                    printf("impossible to give a meaning to the line: %s\n", word[0].c_str());
                    exit(1);
                }
            }
            else if (word.size() == 0) {
            }
            else {
                error(1, 1, "main", "Invalid line in input file: %s", line.c_str());
            }
        }
    }
    error(files.size() != 37, 1, "main", "input files in  file %s must be 37 instead they are %d", argv[4], files.size());


    // mysprintf(namefile, NAMESIZE, "%s/%s_onlinemeas_B64.dat_reweight_charm_0.1_OS_B64.dat", argv[2], argv[1]);
    // files.emplace_back(namefile);
    // mysprintf(namefile, NAMESIZE, "%s/%s_onlinemeas_B64.dat_reweight_step_0.1232229005_to_0.1234969805_OS_B64.dat", argv[2], argv[1]);
    // files.emplace_back(namefile);
    // mysprintf(namefile, NAMESIZE, "%s/%s_onlinemeas_B64.dat_reweight_step_0.03125_to_0.03152408_OS_B64.dat", argv[2], argv[1]);
    // files.emplace_back(namefile);
    // mysprintf(namefile, NAMESIZE, "%s/%s_onlinemeas_B64.dat_reweight_strange_OS_B64.dat", argv[2], argv[1]);
    // files.emplace_back(namefile);
    // mysprintf(namefile, NAMESIZE, "%s/%s_onlinemeas_B64.dat_reweight_step_0.01_to_0.01027408_OS_B64.dat", argv[2], argv[1]);
    // files.emplace_back(namefile);

    // std::vector<int> myen(files.size());
    // for (int e = 0; e < files.size(); e++)
    // {
    //     myen[e] = e;
    // }

    // error(files.size() != 4, 1, "main", "No input files found in  file %s we need 4 lines:\n file for w0 \nml derivative \nms derivative  \nmc derivative", argv[4]);
    // data_all jackall = read_all_the_files(files, argv[1]);
    // printf("we read all\n");
    // jackall.create_generalised_resampling();
    // // data_all jackall;


    // jack setup
    int Njack = 21;

    if (strcmp(argv[1], "jack") == 0) {
        myres = new resampling_jack(Njack - 1);
    }
    else if (strcmp(argv[1], "boot") == 0) {
        myres = new resampling_boot(Njack - 1);
    }


    //////////////////////////////////////////////////////////////
    //  read m^iso
    //////////////////////////////////////////////////////////////
    char** option_read;
    option_read = (char**)malloc(sizeof(char*) * 7);
    option_read[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option_read[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option_read[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option_read[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option_read[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option_read[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option_read[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option_read[1], NAMESIZE, "read_plateaux");        // blind/see/read_plateaux
    mysprintf(option_read[2], NAMESIZE, "-p");                   // -p
    mysprintf(option_read[3], NAMESIZE, "../../data/");          // path
    mysprintf(option_read[4], NAMESIZE, argv[1]);                // resampling
    mysprintf(option_read[5], NAMESIZE, "no");                   // pdf
    mysprintf(option_read[6], NAMESIZE, "%s", files[35].c_str()); // infile

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");
    double mean, err;
    int seed;
    line_read_param(option_read, "seed", mean, err, seed, namefile_plateaux);
    double* tmp = myres->create_fake(mean, err, seed); // just to set the seed
    free(tmp);

    std::vector<double*> amuiso(3);
    line_read_param(option_read, "muliso", mean, err, seed, namefile_plateaux);
    amuiso[0] = myres->create_fake(mean, err, -1);
    line_read_param(option_read, "musiso", mean, err, seed, namefile_plateaux);
    amuiso[1] = myres->create_fake(mean, err, -1);
    line_read_param(option_read, "muciso", mean, err, seed, namefile_plateaux);
    amuiso[2] = myres->create_fake(mean, err, -1);

    std::vector<double*> amusim(3);
    line_read_param(option_read, "mulsim", mean, err, seed, namefile_plateaux);
    amusim[0] = myres->create_fake(mean, err, -1);
    line_read_param(option_read, "mussim", mean, err, seed, namefile_plateaux);
    amusim[1] = myres->create_fake(mean, err, -1);
    line_read_param(option_read, "mucsim", mean, err, seed, namefile_plateaux);
    amusim[2] = myres->create_fake(mean, err, -1);

    double* mpcac_over_mu;
    line_read_param(option_read, "mpcac_over_mu", mean, err, seed, namefile_plateaux);
    mpcac_over_mu = myres->create_fake(mean, err, -1);

    double* ZA;
    line_read_param(option_read, "ZA", mean, err, seed, namefile_plateaux);
    ZA = myres->create_fake(mean, err, -1);

    int L, Lerr;
    line_read_param(option_read, "L", L, Lerr, seed, namefile_plateaux);

    double* previous_a;
    line_read_param(option_read, "a", mean, err, seed, namefile_plateaux);
    previous_a = myres->create_fake(mean, err, -1);



    //////////////////////////////////////////////////////////////
    // read the jacks
    //////////////////////////////////////////////////////////////
    double** data = malloc_2<double>(5 * (1 + 3 * 2), Njack);

    // data_all jackall_chi = read_all_the_files(files, argv[1]);
    // jackall_chi.create_generalised_resampling();

    // sim values
    printf("////////////// readding sim jacks\n");
    read_file_debug(data[0], files[0].c_str());// Mpi
    read_file_debug(data[1], files[1].c_str());// MK
    read_file_debug(data[2], files[2].c_str());// MDs
    read_file_debug(data[3], files[3].c_str());// fpi
    read_file_debug(data[4], files[4].c_str());// w0

    // val deriv ml
    printf("////////////// readding val deriv ml jacks\n");
    read_file_debug(data[0 + 5], files[0 + 5].c_str());// dRpi2 / dmul
    read_file_debug(data[1 + 5], files[1 + 5].c_str());// dRK2  / dmul
    read_file_debug(data[2 + 5], files[2 + 5].c_str());// dRDs / dmul
    read_file_debug(data[3 + 5], files[3 + 5].c_str());// ddpi / dmul
    read_file_debug(data[4 + 5], files[4 + 5].c_str());// dw0 / dmul
    for (int j = 0; j < Njack;j++) {
        for (int i = 0; i < 3; i++) {
            double fpi = data[3][j];
            double M = data[i][j];
            double dR2 = data[i + 5][j];
            double df = data[3 + 5][j];
            if (i == 2) { // for the D2 we do not have the square but actually R
                data[i + 5][j] = (dR2 + M * df / (fpi * fpi)) * fpi;
            }
            else {
                data[i + 5][j] = dR2 + 2.0 * M * M * df / (fpi * fpi * fpi);
                data[i + 5][j] *= fpi * fpi / (2.0 * M);
            }
        }
    }
    // val deriv ms
    printf("////////////// readding val deriv ms jacks\n");
    read_file_debug(data[0 + 5 * 2], files[0 + 5 * 2].c_str());// dMpi / dmus
    read_file_debug(data[1 + 5 * 2], files[1 + 5 * 2].c_str());// dMK2  / dmus
    read_file_debug(data[2 + 5 * 2], files[2 + 5 * 2].c_str());// dMDs / dmus
    read_file_debug(data[3 + 5 * 2], files[3 + 5 * 2].c_str());// ddpi / dmus
    read_file_debug(data[4 + 5 * 2], files[4 + 5 * 2].c_str());// dw0 / dmus
    for (int j = 0; j < Njack;j++) {
        double MK = data[1][j];
        data[1 + 5 * 2][j] /= (2 * MK);
    }
    // val deriv mc
    printf("////////////// readding val deriv mc jacks\n");
    read_file_debug(data[0 + 5 * 3], files[0 + 5 * 3].c_str());// dMpi / dmuc
    read_file_debug(data[1 + 5 * 3], files[1 + 5 * 3].c_str());// dMK  / dmuc
    read_file_debug(data[2 + 5 * 3], files[2 + 5 * 3].c_str());// dMDs / dmuc
    read_file_debug(data[3 + 5 * 3], files[3 + 5 * 3].c_str());// ddpi / dmuc
    read_file_debug(data[4 + 5 * 3], files[4 + 5 * 3].c_str());// dw0 / dmuc
    // sea deriv ml
    printf("////////////// readding sea deriv ml jacks\n");
    read_file_debug(data[0 + 5 * 4], files[0 + 5 * 4].c_str());// dMpi / dmul(sea)
    read_file_debug(data[1 + 5 * 4], files[1 + 5 * 4].c_str());// dMK  / dmul(sea)
    read_file_debug(data[2 + 5 * 4], files[2 + 5 * 4].c_str());// dMDs / dmul(sea)
    read_file_debug(data[3 + 5 * 4], files[3 + 5 * 4].c_str());// dpi / dmul(sea)
    read_file_debug(data[4 + 5 * 4], files[4 + 5 * 4].c_str());// dw0 / dmul(sea)
    for (int j = 0; j < Njack;j++) { // we need a factor 2 in w0 because of the two quarks u and d
        data[4 + 5 * 4][j] *= 2;
    }
    // sea deriv ms
    printf("////////////// readding sea deriv ms jacks\n");
    read_file_debug(data[0 + 5 * 5], files[0 + 5 * 5].c_str());// dMpi / dmus(sea)
    read_file_debug(data[1 + 5 * 5], files[1 + 5 * 5].c_str());// dMK2  / dmus(sea)
    read_file_debug(data[2 + 5 * 5], files[2 + 5 * 5].c_str());// dMDs / dmus(sea)
    read_file_debug(data[3 + 5 * 5], files[3 + 5 * 5].c_str());// ddpi / dmus(sea)
    read_file_debug(data[4 + 5 * 5], files[4 + 5 * 5].c_str());// dw0 / dmus(sea)
    for (int j = 0; j < Njack;j++) {
        double mul = amuiso[0][j];
        double f = data[3][j];
        double w0 = data[4][j];
        double* dR = &data[3 + 5 * 5][j]; // fpi
        *dR = (*dR) * (f) / mul;
        dR = &data[4 + 5 * 5][j]; // w0
        *dR = (*dR) * (w0) / mul;
        for (int i = 0; i < 3; i++) {
            // first we remove the prefacor
            double* dR = &data[i + 5 * 5][j];
            double M = data[i][j];
            if (j == Njack - 1 && i == 2)    printf(" %g     %g   %g  %g \n", *dR, f, M, mul);
            double df = data[3 + 5 * 5][j];
            *dR = (*dR) * (M / f) / mul;
            if (j == Njack - 1 && i == 2)    printf(" %g    %g \n", *dR, df);
            (*dR) = ((*dR) + M * df / (f * f)) * f;
            if (j == Njack - 1 && i == 2)    printf(" %g     \n", *dR);
        }
    }
    // sea deriv mc
    printf("////////////// readding sea deriv mc jacks\n");
    read_file_debug(data[0 + 5 * 6], files[0 + 5 * 6].c_str());// dMpi / dmuc(sea)
    read_file_debug(data[1 + 5 * 6], files[1 + 5 * 6].c_str());// dMK2 / dmuc(sea)
    read_file_debug(data[2 + 5 * 6], files[2 + 5 * 6].c_str());// dMDs / dmuc(sea)
    read_file_debug(data[3 + 5 * 6], files[3 + 5 * 6].c_str());// ddpi / dmuc(sea)
    read_file_debug(data[4 + 5 * 6], files[4 + 5 * 6].c_str());// dw0  / dmuc(sea)
    for (int j = 0; j < Njack;j++) {
        double mul = amuiso[0][j];
        double f = data[3][j];
        double w0 = data[4][j];
        double* dR = &data[3 + 5 * 6][j]; // femove prefactor in dfpi/dmc(sea)
        *dR = (*dR) * (f) / mul;
        dR = &data[4 + 5 * 6][j]; // remove prefactor in dw0/dmc(sea)
        *dR = (*dR) * (w0) / mul;
        for (int i = 0; i < 3; i++) {
            // first we remove the prefacor
            double* dR = &data[i + 5 * 6][j];
            double M = data[i][j];
            double df = data[3 + 5 * 6][j];
            *dR = (*dR) * (M / f) / mul;
            (*dR) = ((*dR) + M * df / (f * f)) * f;
        }
    }

    mysprintf(namefile, NAMESIZE, "%s/%s_scale_setting_system_%s", argv[2], argv[1], files[36].c_str());
    FILE* jack_file = open_file(namefile, "w+");
    generic_header head;
    head.Njack = Njack;
    head.T = L * 2;
    head.L = L;
    head.ncorr = 1;
    head.beta = 0;
    head.kappa = 0;
    head.mus = std::vector<double>(0);
    head.thetas = std::vector<double>(0);
    head.rs = std::vector<double>(0);
    head.gammas = std::vector<std::string>(0);
    head.smearing = std::vector<std::string>(0);
    head.bananas = std::vector<int>(0);
    head.oranges = std::vector<double>(0);
    head.size =head.ncorr * 2 * head.T;
    // write_header_g2(jack_file, head);
    head.write_header(jack_file);

    corr_counter = -1;
    for (int i = 0;i < 30;i++) {
        write_jack(data[i], Njack, jack_file);
    }

    //////////////////////////////////////////////////////////////
    // max twist correction
    //////////////////////////////////////////////////////////////
    printf("////////////// applying max twist correction\n");
    printf("Mpi before max twist correction: %.12g  %.12g\n", data[0][Njack - 1], myres->comp_error(data[0]));
    printf("fpi before max twist correction: %.12g  %.12g\n", data[3][Njack - 1], myres->comp_error(data[3]));
    for (int j = 0; j < Njack;j++) {
        double mr = ZA[j] * mpcac_over_mu[j];
        double cl = sqrt(1 + mr * mr);
        double* Mpi = &data[0][j];
        double* fpi = &data[3][j];
        *Mpi = *Mpi / sqrt(cl);
        *fpi = *fpi * cl;
    }
    printf("Mpi after max twist correction: %.12g  %.12g\n", data[0][Njack - 1], myres->comp_error(data[0]));
    printf("fpi after max twist correction: %.12g  %.12g\n", data[3][Njack - 1], myres->comp_error(data[3]));
    //////////////////////////////////////////////////////////////
    // FVE
    //////////////////////////////////////////////////////////////
    printf("////////////// applying max twist correction\n");
    printf("Mpi before FVE correction: %.12g  %.12g\n", data[0][Njack - 1], myres->comp_error(data[0]));
    printf("fpi before FVE correction: %.12g  %.12g\n", data[3][Njack - 1], myres->comp_error(data[3]));
    for (int j = 0; j < Njack;j++) {
        double* Mpi = &data[0][j];
        double* fpi = &data[3][j];
        double xi = (*Mpi) * (*Mpi) / ((*fpi) * (*fpi));
        double delta_FVE = FVE_GL_Mpi(L, xi, (*fpi));
        *Mpi /= (1 - 0.25 * delta_FVE);
        *fpi /= (1 + delta_FVE);

    }
    printf("Mpi after FVE correction: %.12g  %.12g\n", data[0][Njack - 1], myres->comp_error(data[0]));
    printf("fpi after FVE correction: %.12g  %.12g\n", data[3][Njack - 1], myres->comp_error(data[3]));

    //////////////////////////////////////////////////////////////
    // sistemone fpi
    //////////////////////////////////////////////////////////////
    enum class Meson_enum {
        Pi = 0,
        K = 1,
        Ds = 2,
        fPi = 3,
        w0 = 4
    };
    enum class quark_enum {
        ml = 0,
        ms = 1,
        mc = 2
    };

    for (int iM = 0; iM < 5; iM++) {

        double* M = data[id_deriv(iM, 0, 0, 0)];
        printf("Obs %d val=%-8.3g (%-8.3g)\t", iM, M[Njack - 1], myres->comp_error(M));
        // valence
        for (int im = 0; im < 3; im++) {
            double* dM = data[id_deriv(iM, 1, im, 0)];
            printf(" dval%d = %-8.3g (%-8.3g)", im, dM[Njack - 1], myres->comp_error(dM));
        }
        for (int im = 0; im < 3; im++) {
            double* dM = data[id_deriv(iM, 1, im, 1)];
            printf(" dsea%d = %-8.3g (%-8.3g)", im, dM[Njack - 1], myres->comp_error(dM));
        }
        printf("\n");
    }
    double* Rds = myres->create_copy(data[id_deriv(2, 0, 0, 0)]);
    myres->div(Rds, Rds, data[id_deriv(3, 0, 0, 0)]);
    // printf("Rds = %.12g  %.12g\n", Rds[Njack-1], myres->comp_error(Rds));

    double** Mat = malloc_2<double>(3, 3);
    double* y = (double*)malloc(sizeof(double) * 3);
    double** miso = malloc_2<double>(3, Njack);
    for (int j = 0; j < Njack;j++) {

        y[0] = Mpi_MeV / fpi_MeV;
        y[1] = MK_MeV / fpi_MeV;
        y[2] = MDs_MeV / fpi_MeV;
        int ifpi = 3;
        double f = data[id_deriv(ifpi, 0, 0, 0)][j];

        for (int iM = 0; iM < 2; iM++) {

            double M = data[id_deriv(iM, 0, 0, 0)][j];

            for (int im = 0; im < 3; im++) {
                Mat[iM][im] = 0.0;
                // valence
                double dM = data[id_deriv(iM, 1, im, 0)][j];
                double df = data[id_deriv(ifpi, 1, im, 0)][j];
                Mat[iM][im] += 2 * M * dM / (f * f) - 2 * M * M * df / (f * f * f);
                // sea
                dM = data[id_deriv(iM, 1, im, 1)][j];
                df = data[id_deriv(ifpi, 1, im, 1)][j];
                Mat[iM][im] += 2 * M * dM / (f * f) - 2 * M * M * df / (f * f * f);
            }
            y[iM] *= y[iM];
            // if (j==Njack-1) printf(" %.12g   %g   %g   %g\n ",y[iM], M / f, f, M);
            y[iM] -= (M / f) * (M / f);

        }
        for (int iM = 2; iM < 3; iM++) {

            double M = data[id_deriv(iM, 0, 0, 0)][j];

            for (int im = 0; im < 3; im++) {
                Mat[iM][im] = 0.0;
                // valence
                double dM = data[id_deriv(iM, 1, im, 0)][j];
                double df = data[id_deriv(ifpi, 1, im, 0)][j];
                Mat[iM][im] += dM / f - M * df / (f * f);
                // sea
                dM = data[id_deriv(iM, 1, im, 1)][j];
                df = data[id_deriv(ifpi, 1, im, 1)][j];
                Mat[iM][im] += dM / f - M * df / (f * f);
            }
            // if (j==Njack-1) printf(" %.12g   %g   %g   %g\n ",y[iM], M / f, f, M);
            y[iM] -= M / f;

        }
        double* P = LU_decomposition_solver(3, Mat, y);
        miso[0][j] = (amusim[0][j] + P[0]);
        miso[1][j] = (amusim[1][j] + P[1]);
        miso[2][j] = (amusim[2][j] + P[2]);
        if (j == Njack - 1) {
            printf("Matrix for m^iso solution (jackknife %d):\n", j);
            for (int ii = 0; ii < 3; ii++) {
                for (int jj = 0; jj < 3; jj++) {
                    printf("%g ", Mat[ii][jj]);
                }
                printf("\n");
            }
            printf("RHS:\n");
            for (int ii = 0; ii < 3; ii++) {
                printf("%g\n", y[ii]);
            }
            printf("solution:\n");
            for (int ii = 0; ii < 3; ii++) {
                printf("%g\n", P[ii]);
            }
        }
        free(P);
    }
    printf("Results for m^iso (MeV):\n");
    for (int i = 0; i < 3; i++) {
        double mean = myres->mean(miso[i]);
        double err = myres->comp_error(miso[i]);
        printf("m^iso_%d = %-12g +/- %-12g   starting from %-12g +/- %-12g\n", i, mean, err, myres->mean(amuiso[i]), myres->comp_error(amuiso[i]));
    }
    printf("sim values:\n");
    for (int i = 0; i < 3; i++) {
        printf("m^sim_%d = %g \n", i, myres->mean(amusim[i]));
    }




    double* a_fm = myres->create_copy(data[4]);
    double* w0_from_fpi = myres->create_copy(data[4]);
    for (int j = 0; j < Njack;j++) {
        int ifpi = 3;
        double af = data[id_deriv(ifpi, 0, 0, 0)][j];
        int iw0 = 4;
        double w_a = data[id_deriv(iw0, 0, 0, 0)][j];
        for (int im = 0; im < 3; im++) {
            for (int val_sea = 0; val_sea < 2; val_sea++) {
                double df = data[id_deriv(ifpi, 1, im, val_sea)][j];
                double dm = (miso[im][j] - amusim[im][j]);
                af += dm * df;
                double dw = data[id_deriv(iw0, 1, im, val_sea)][j];
                w_a += dm * dw;
            }
        }
        a_fm[j] = af / (fpi_MeV / hbarc);
        w0_from_fpi[j] = w_a * a_fm[j];
    }
    double* diff_a = myres->create_copy(a_fm);
    myres->sub(diff_a, a_fm, previous_a);
    printf("lattice spacing (fm): %g +/- %g   vs previous value %g +/- %g  diff %g +/- %g  \n", myres->mean(a_fm), myres->comp_error(a_fm),
        myres->mean(previous_a), myres->comp_error(previous_a),
        myres->mean(diff_a), myres->comp_error(diff_a));
    printf("w0 (fm): %g +/- %g\n", myres->mean(w0_from_fpi), myres->comp_error(w0_from_fpi));
    char name_out[NAMESIZE];
    mysprintf(name_out, NAMESIZE, "scale_setting/%s_a_from_fpi.jack", files[36].c_str());
    myres->write_jack_in_file(a_fm, name_out);
    mysprintf(name_out, NAMESIZE, "scale_setting/%s_w0_from_fpi.jack", files[36].c_str());
    myres->write_jack_in_file(w0_from_fpi, name_out);

    write_jack(miso[0], Njack, jack_file);     check_correlatro_counter(30);
    write_jack(miso[1], Njack, jack_file);     check_correlatro_counter(31);
    write_jack(miso[2], Njack, jack_file);     check_correlatro_counter(32);
    write_jack(a_fm, Njack, jack_file);     check_correlatro_counter(33);
    write_jack(w0_from_fpi, Njack, jack_file);     check_correlatro_counter(34);

    printf("//////////////////////////////////////////////////////////////\n");
    printf("// w0 systemone\n");
    printf("//////////////////////////////////////////////////////////////\n");

    double** miso_w0 = malloc_2<double>(3, Njack);
    for (int j = 0; j < Njack;j++) {

        y[0] = Mpi_MeV * w0_MeV;
        y[1] = MK_MeV * w0_MeV;
        y[2] = MDs_MeV * w0_MeV;
        int iw0 = 4;
        double w0 = data[id_deriv(iw0, 0, 0, 0)][j];

        for (int iM = 0; iM < 3; iM++) {

            double M = data[id_deriv(iM, 0, 0, 0)][j];

            for (int im = 0; im < 3; im++) {
                Mat[iM][im] = 0.0;
                // valence
                double dM = data[id_deriv(iM, 1, im, 0)][j];
                double dw = data[id_deriv(iw0, 1, im, 0)][j];
                Mat[iM][im] += dM * w0 + M * dw;
                // sea
                dM = data[id_deriv(iM, 1, im, 1)][j];
                dw = data[id_deriv(iw0, 1, im, 1)][j];
                Mat[iM][im] += dM * w0 + M * dw;
            }
            y[iM] -= M * w0;

        }
        double* P = LU_decomposition_solver(3, Mat, y);
        miso_w0[0][j] = (amusim[0][j] + P[0]);
        miso_w0[1][j] = (amusim[1][j] + P[1]);
        miso_w0[2][j] = (amusim[2][j] + P[2]);
        if (j == Njack - 1) {
            printf("Matrix for m^iso solution (jackknife %d):\n", j);
            for (int ii = 0; ii < 3; ii++) {
                for (int jj = 0; jj < 3; jj++) {
                    printf("%g ", Mat[ii][jj]);
                }
                printf("\n");
            }
            printf("RHS:\n");
            for (int ii = 0; ii < 3; ii++) {
                printf("%g\n", y[ii]);
            }
            printf("solution:\n");
            for (int ii = 0; ii < 3; ii++) {
                printf("%g\n", P[ii]);
            }
        }
        free(P);
    }
    printf("Results for m^iso (MeV):\n");
    for (int i = 0; i < 3; i++) {
        double mean = myres->mean(miso_w0[i]);
        double err = myres->comp_error(miso_w0[i]);
        printf("m^iso_%d = %-12g +/- %-12g   starting from %-12g +/- %-12g\n", i, mean, err, myres->mean(amuiso[i]), myres->comp_error(amuiso[i]));
    }
    printf("sim values:\n");
    for (int i = 0; i < 3; i++) {
        printf("m^sim_%d = %g \n", i, myres->mean(amusim[i]));
    }

    double* a_from_w0 = myres->create_copy(data[4]);
    double* fpi_from_w0 = myres->create_copy(data[4]);
    for (int j = 0; j < Njack;j++) {
        int ifpi = 3;
        double af = data[id_deriv(ifpi, 0, 0, 0)][j];
        int iw0 = 4;
        double w_a = data[id_deriv(iw0, 0, 0, 0)][j];
        for (int im = 0; im < 3; im++) {
            for (int val_sea = 0; val_sea < 2; val_sea++) {
                double df = data[id_deriv(ifpi, 1, im, val_sea)][j];
                double dm = (miso_w0[im][j] - amusim[im][j]);
                af += dm * df;
                double dw = data[id_deriv(iw0, 1, im, val_sea)][j];
                w_a += dm * dw;
            }
        }
        a_from_w0[j] = w0_fm / w_a;
        fpi_from_w0[j] = af / (a_from_w0[j] / hbarc);
    }
    printf("lattice spacing (fm): %g +/- %g\n", myres->mean(a_from_w0), myres->comp_error(a_from_w0));
    printf("fpi (fm): %g +/- %g\n", myres->mean(fpi_from_w0), myres->comp_error(fpi_from_w0));

    mysprintf(name_out, NAMESIZE, "scale_setting/%s_a_from_w0.jack", files[36].c_str());
    myres->write_jack_in_file(a_from_w0, name_out);
    mysprintf(name_out, NAMESIZE, "scale_setting/%s_fpi_from_w0.jack", files[36].c_str());
    myres->write_jack_in_file(fpi_from_w0, name_out);

    write_jack(miso_w0[0], Njack, jack_file);     check_correlatro_counter(35);
    write_jack(miso_w0[1], Njack, jack_file);     check_correlatro_counter(36);
    write_jack(miso_w0[2], Njack, jack_file);     check_correlatro_counter(37);
    write_jack(a_from_w0, Njack, jack_file);     check_correlatro_counter(38);
    write_jack(fpi_from_w0, Njack, jack_file);     check_correlatro_counter(39);


    return 0;
}
