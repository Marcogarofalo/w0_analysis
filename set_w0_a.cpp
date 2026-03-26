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


int id_deriv_fpi(int reg, int deriv, int quark, int val_sea) {
    return 46 + reg * 7 + deriv * (1 + (quark + 3 * val_sea));
}

void weighted_average(double* M0, double* M1) {
    double dM0 = myres->comp_error(M0);
    double dM1 = myres->comp_error(M1);
    double wM0 = 1.0 / (dM0 * dM0);
    double wM1 = 1.0 / (dM1 * dM1);

    for (int j = 0; j < myres->Njack; j++) {
        M0[j] = (wM0 * M0[j] + wM1 * M1[j]) / (wM0 + wM1);
    }
}

double get_linear_deriv_w0c(double** data, std::vector<double*>amuiso, double* previous_a, int j) {
    double dw;
    double P0 = data[41][j];
    double P1 = data[42][j];
    // remove the prefactor in dw0/dmc(sea)
    double mul = amuiso[0][j];
    double w0 = data[4][j];
    double a = previous_a[j];
    P0 *= w0 / mul;
    P1 *= w0 / mul;
    dw = P0 + P1 * a * a;
    return dw;
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
    // error(files.size() != 41, 1, "main", "input files in  file %s must be 41 instead they are %d", argv[4], files.size());


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

    std::string ensemble = files[36];

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
    int file_to_read = 5 * (1 + 3 * 2) + 4 + 2 + 2 + 1 + 3 + 7 * 2 + 2;
    printf("%d\n", file_to_read);
    error(files.size() - 2 != file_to_read, 1, "main", "No input files found in  file %s we need %d lines but we have %d", argv[4], file_to_read, files.size() - 2);
    double** data = malloc_2<double>(file_to_read, Njack);

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
        data[0 + 5 * 4][j] *= 2;
        data[1 + 5 * 4][j] *= 2;
        data[2 + 5 * 4][j] *= 2; // this should be zero but we put the factor 2 for consistency with the other two
        data[3 + 5 * 4][j] *= 2;
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

    // decorrelate jacks
    // for (int iM = 0; iM < 5;iM++) {
    //     for (int im = 0;im < 3;im++) {
    //         for (int sv = 0;sv < 2;sv++) {
    //             int idw0 = iM;
    //             int id = id_deriv(idw0, 1, im, sv);
    //             double* dw0 = data[id];
    //             double m = dw0[Njack - 1];
    //             double err = myres->comp_error(dw0);
    //             printf("old jack = %g  %g\n", m, err);
    //             free(dw0);
    //             data[id] = myres->create_fake(m, err, -1);
    //             printf("new jack = %g  %g\n", data[id][Njack - 1], myres->comp_error(data[id]));
    //         }
    //     }
    // }

    // printing
    // std::vector<std::string>  obs = { "Mpi","MK","MDs","fpi","w0",
    // "Mpi_ml_val","dMK_ml_val","MDs_ml_val","fpi_ml_val","w0_ml_val",
    // "Mpi_ms_val","dMK_ms_val","MDs_ms_val","fpi_ms_val","w0_ms_val",
    // "Mpi_mc_val","dMK_mc_val","MDs_mc_val","fpi_mc_val","w0_mc_val",
    // "Mpi_ml_sea","dMK_ml_sea","MDs_ml_sea","fpi_ml_sea","w0_ml_sea",
    // "Mpi_ms_sea","dMK_ms_sea","MDs_ms_sea","fpi_ms_sea","w0_ms_sea",
    // "Mpi_mc_sea","dMK_mc_sea","MDs_mc_sea","fpi_mc_sea","w0_mc_sea"
    // };
    // double *md=(double*) malloc(sizeof(double)*35);
    // double** cov = myres->comp_cov(35, data);
    // for (int i = 0;i < 35;i++) {
    //     for (int j = i + 1;j < 35;j++) {
    //         double corr = cov[i][j] / sqrt(cov[i][i] * cov[j][j]);
    //         if (std::fabs(corr) > 0.4)
    //             printf("corr  %s   %s  %-7.3f\n", obs[i].c_str(), obs[j].c_str(), corr);
    //     }
    //     md[i] = myres->mean(data[i]); 
    // }

    // set off-diag to zero
    // for (int i = 0;i < 35;i++) {
    //     for (int j = 0 ;j < 35;j++) {
    //         if (i!=j) cov[i][j]=0;
    //         else cov[i][j]+=1e-14;
    //     }
    // }
    // printf("after\n");

    // regenerate fake jackk
    // data = myres->create_fake_covariance(md,35,cov,-1);
    // cov = myres->comp_cov(35, data);
    // for (int i = 0;i < 35;i++) {
    //     for (int j = i + 1;j < 35;j++) {
    //         double corr = cov[i][j] / sqrt(cov[i][i] * cov[j][j]);
    //         if (std::fabs(corr) > 0.4)
    //             printf("corr  %s   %s  %-7.3f\n", obs[i].c_str(), obs[j].c_str(), corr);
    //     }
    //     md[i] = myres->mean(data[i]); 

    // }

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
    head.size = head.ncorr * 2 * head.T;
    // write_header_g2(jack_file, head);
    head.write_header(jack_file);

    corr_counter = -1;
    for (int i = 0;i < 30;i++) { // why untill 30, now I can not change
        write_jack(data[i], Njack, jack_file);
    }

    // read delta m0 and m0 derivative
    read_file_debug(data[35], files[37].c_str());// dm0
    read_file_debug(data[36], files[38].c_str());// dw0/dm0l
    read_file_debug(data[37], files[39].c_str());// dw0/dm0s
    read_file_debug(data[38], files[40].c_str());// dw0/dm0c

    read_file_debug(data[39], files[41].c_str());// dw0/dms_sea
    read_file_debug(data[40], files[42].c_str());// dw0/dmc_sea


    static constexpr int iw0 = 4;

    read_file_debug(data[41], files[43].c_str());// dw0/dms_sea P0
    read_file_debug(data[42], files[44].c_str());// dw0/dmc_sea P1

    read_file_debug(data[43], files[45].c_str());// mpcac
    for (int j = 0; j < Njack;j++) {
        mpcac_over_mu[j] = data[43][j] / amusim[0][j];
    }

    std::vector<std::vector<int>> id_to_correct = { {0,3} };
    std::vector<int> Ls = { L };
    std::vector<bool> correct_M = { true };
    std::vector<bool> correct_M_twist = { true };
    std::vector<std::vector<int>> id_average;
    if (strcmp(files[46].c_str(), "skip") != 0 &&
        strcmp(files[47].c_str(), "skip") != 0 &&
        strcmp(files[48].c_str(), "skip") != 0) {
        // id_to_correct = { {0,3}, {44,45} };
        id_to_correct.push_back({ 44, 45 }); // 
        read_file_debug(data[44], files[46].c_str());// mpi 
        read_file_debug(data[45], files[47].c_str());// fpi
        Ls.push_back(std::stoi(files[48]));
        correct_M.push_back(true);
        correct_M_twist.push_back(true);
        id_average.push_back({ 0, 44 });
        id_average.push_back({ 3, 45 });
    }

    read_file_debug(data[46], files[49].c_str());// fpi A0 tm
    read_file_debug(data[47], files[50].c_str());// deriv val ml fpi A0 
    read_file_debug(data[48], files[51].c_str());// deriv val ms fpi A0 
    read_file_debug(data[49], files[52].c_str());// deriv val mc fpi A0 
    read_file_debug(data[50], files[53].c_str());// deriv sea ml fpi A0 
    read_file_debug(data[51], files[54].c_str());// deriv sea ms fpi A0 
    read_file_debug(data[52], files[55].c_str());// deriv sea mc fpi A0

    read_file_debug(data[53], files[56].c_str());// fpi A0 OS
    read_file_debug(data[54], files[57].c_str());// deriv val ml fpi A0 
    read_file_debug(data[55], files[58].c_str());// deriv val ms fpi A0 
    read_file_debug(data[56], files[59].c_str());// deriv val mc fpi A0 
    read_file_debug(data[57], files[60].c_str());// deriv sea ml fpi A0 
    read_file_debug(data[58], files[61].c_str());// deriv sea ms fpi A0 
    read_file_debug(data[59], files[62].c_str());// deriv sea mc fpi A0

    error(id_deriv_fpi(0, 0, 0, 0) != 46, 1, "main", "id_deriv_fpi(0, 0, 0, 0) should be 46 but it is %d", id_deriv_fpi(0, 0, 0, 0));
    error(id_deriv_fpi(0, 1, 0, 0) != 47, 1, "main", "id_deriv_fpi(0, 1, 0, 0) should be 47 but it is %d", id_deriv_fpi(0, 1, 0, 0));
    error(id_deriv_fpi(0, 1, 1, 0) != 48, 1, "main", "id_deriv_fpi(0, 1, 1, 0) should be 48 but it is %d", id_deriv_fpi(0, 1, 1, 0));
    error(id_deriv_fpi(0, 1, 2, 0) != 49, 1, "main", "id_deriv_fpi(0, 1, 2, 0) should be 49 but it is %d", id_deriv_fpi(0, 1, 2, 0));
    error(id_deriv_fpi(0, 1, 0, 1) != 50, 1, "main", "id_deriv_fpi(0, 1, 0, 1) should be 50 but it is %d", id_deriv_fpi(0, 1, 0, 1));
    error(id_deriv_fpi(0, 1, 1, 1) != 51, 1, "main", "id_deriv_fpi(0, 1, 1, 1) should be 51 but it is %d", id_deriv_fpi(0, 1, 1, 1));
    error(id_deriv_fpi(0, 1, 2, 1) != 52, 1, "main", "id_deriv_fpi(0, 1, 2, 1) should be 52 but it is %d", id_deriv_fpi(0, 1, 2, 1));

    error(id_deriv_fpi(1, 0, 0, 0) != 53, 1, "main", "id_deriv_fpi(1, 0, 0, 0) should be 53 but it is %d", id_deriv_fpi(1, 0, 0, 0));
    error(id_deriv_fpi(1, 1, 0, 0) != 54, 1, "main", "id_deriv_fpi(1, 1, 0, 0) should be 54 but it is %d", id_deriv_fpi(1, 1, 0, 0));
    error(id_deriv_fpi(1, 1, 1, 0) != 55, 1, "main", "id_deriv_fpi(1, 1, 1, 0) should be 55 but it is %d", id_deriv_fpi(1, 1, 1, 0));
    error(id_deriv_fpi(1, 1, 2, 0) != 56, 1, "main", "id_deriv_fpi(1, 1, 2, 0) should be 56 but it is %d", id_deriv_fpi(1, 1, 2, 0));
    error(id_deriv_fpi(1, 1, 0, 1) != 57, 1, "main", "id_deriv_fpi(1, 1, 0, 1) should be 57 but it is %d", id_deriv_fpi(1, 1, 0, 1));
    error(id_deriv_fpi(1, 1, 1, 1) != 58, 1, "main", "id_deriv_fpi(1, 1, 1, 1) should be 58 but it is %d", id_deriv_fpi(1, 1, 1, 1));
    error(id_deriv_fpi(1, 1, 2, 1) != 59, 1, "main", "id_deriv_fpi(1, 1, 2, 1) should be 59 but it is %d", id_deriv_fpi(1, 1, 2, 1));


    for (int j = 0; j < Njack;j++) {
        for (int r = 0; r < 2; r++) {

            int iddf = id_deriv_fpi(r, 1, 0, 1);// ml sea needs a factor 2 because of the two quarks u and d
            double* dR = &data[iddf][j];
            *dR *= 2.0;

            for (int sc = 1;sc < 3;sc++) {
                int iddf = id_deriv_fpi(r, 1, sc, 1);
                double* dR = &data[iddf][j];
                int idf = id_deriv_fpi(r, 0, 0, 0);
                double f = data[idf][j];
                *dR = (*dR) * (f) / amuiso[0][j];
            }
        }
    }

    id_to_correct.push_back({ id_deriv(0, 0, 0, 0)  ,id_deriv_fpi(0, 0, 0, 0) });
    correct_M.push_back(false);
    correct_M_twist.push_back(false);
    Ls.push_back(L);

    // OS does not need max twist correction
    std::vector<std::vector<int>> id_correct_max_twist = id_to_correct;

    id_to_correct.push_back({ id_deriv(0, 0, 0, 0), id_deriv_fpi(1, 0, 0, 0) });
    correct_M.push_back(false);
    correct_M_twist.push_back(false);
    Ls.push_back(L);

    // B96 fpi OS and tm
    if (strcmp(files[63].c_str(), "skip") != 0 &&
        strcmp(files[64].c_str(), "skip") != 0) {

        read_file_debug(data[60], files[63].c_str());// fpi tm B96 
        id_correct_max_twist.push_back({ 44, 60 });// mass , fpi
        id_to_correct.push_back({ 44, 60 });
        correct_M.push_back(false);
        correct_M_twist.push_back(true);
        Ls.push_back(std::stoi(files[48]));
        id_average.push_back({ id_deriv_fpi(0, 0, 0, 0), 60 });


        read_file_debug(data[61], files[64].c_str());// fpi OS B96
        id_to_correct.push_back({ 44, 61 });
        correct_M.push_back(false);
        correct_M_twist.push_back(false);
        Ls.push_back(std::stoi(files[48]));
        id_average.push_back({ id_deriv_fpi(1, 0, 0, 0), 61 });
    }



    //////////////////////////////////////////////////////////////
    // max twist correction
    //////////////////////////////////////////////////////////////
    printf("////////////// applying max twist correction\n");

    for (int i = 0; i < id_correct_max_twist.size(); i++) {
        int idM;
        int idf;

        idM = id_correct_max_twist[i][0];
        idf = id_correct_max_twist[i][1];

        printf("Mpi before max twist correction id=%d: %.12g  %.12g\n", idM, data[idM][Njack - 1], myres->comp_error(data[idM]));
        printf("fpi before max twist correction id=%d: %.12g  %.12g\n", idf, data[idf][Njack - 1], myres->comp_error(data[idf]));
        for (int j = 0; j < Njack;j++) {
            double mr = ZA[j] * mpcac_over_mu[j];
            double cl = sqrt(1 + mr * mr);
            double* Mpi = &data[idM][j];
            double* fpi = &data[idf][j];
            if (correct_M_twist[i]) *Mpi = *Mpi / sqrt(cl);
            *fpi = *fpi * cl;
        }
        printf("Mpi after max twist correction id=%d: %.12g  %.12g\n", idM, data[idM][Njack - 1], myres->comp_error(data[idM]));
        printf("fpi after max twist correction id=%d: %.12g  %.12g\n", idf, data[idf][Njack - 1], myres->comp_error(data[idf]));
    }

    for (int j = 0; j < Njack;j++) {
        double dm0 = data[35][j];
        double dw0_dm0 = data[36][j] + 0.5 * (data[37][j] + data[38][j]);
        double delta = dm0 * dw0_dm0;
        data[iw0][j] += delta;
    }

    //////////////////////////////////////////////////////////////
    // FVE
    //////////////////////////////////////////////////////////////
    printf("////////////// applying FVE correction\n");
    for (int i = 0; i < id_to_correct.size(); i++) {
        int idM;
        int idf;

        idM = id_to_correct[i][0];
        idf = id_to_correct[i][1];

        for (int j = 0; j < Njack;j++) {
            double* Mpi = &data[idM][j];
            double* fpi = &data[idf][j];
            double xi = (*Mpi) * (*Mpi) / ((4 * M_PI * *fpi) * (4 * M_PI * *fpi));
            double delta_FVE = FVE_GL_Mpi(Ls[i], xi, (*fpi));
            // *Mpi /= (1 - 0.25 * delta_FVE);
            // *fpi /= (1 + delta_FVE);

            delta_fve_NNLO_CDH  d = Delta_pi_NNLO_CDH(Ls[i], xi, (*fpi));
            // if (j == Njack - 1) printf("delta FVE M (L=%d) = %g   %g\n", Ls[i], -0.25 * delta_FVE, d.dM);
            // if (j == Njack - 1) printf("delta FVE f (L=%d) = %g   %g\n", Ls[i], delta_FVE, d.df);
            if (correct_M[i])*Mpi /= (1 + d.dM);
            *fpi /= (1 + d.df);

        }
        printf("Mpi after FVE correction id[%d] L=%d: %.12g  %.12g\n", idM, Ls[i], data[idM][Njack - 1], myres->comp_error(data[idM]));
        printf("fpi after FVE correction id[%d] L=%d: %.12g  %.12g\n", idf, Ls[i], data[idf][Njack - 1], myres->comp_error(data[idf]));
    }


    // if (strcmp(files[46].c_str(), "skip") != 0 && strcmp(files[47].c_str(), "skip") != 0) {
    //     printf("averaging volumes:\n");
    //     int idM0 = id_to_correct[0][0];
    //     int idf0 = id_to_correct[0][1];
    //     int idM1 = id_to_correct[1][0];
    //     int idf1 = id_to_correct[1][1];
    //     // double dM0 = myres->comp_error(data[idM0]);
    //     // double df0 = myres->comp_error(data[idf0]);
    //     // double dM1 = myres->comp_error(data[idM1]);
    //     // double df1 = myres->comp_error(data[idf1]);
    //     // double wM0 = 1.0 / (dM0 * dM0);
    //     // double wM1 = 1.0 / (dM1 * dM1);
    //     // double wf0 = 1.0 / (df0 * df0);
    //     // double wf1 = 1.0 / (df1 * df1);

    //     // for (int j = 0; j < Njack;j++) {
    //     //     double* M0 = &data[idM0][j];
    //     //     double* f0 = &data[idf0][j];
    //     //     double* M1 = &data[idM1][j];
    //     //     double* f1 = &data[idf1][j];
    //     //     *M0 = (wM0 * (*M0) + wM1 * (*M1)) / (wM0 + wM1);
    //     //     *f0 = (wf0 * (*f0) + wf1 * (*f1)) / (wf0 + wf1);
    //     // }
    //     weighted_average(data[idM0], data[idM1]);
    //     weighted_average(data[idf0], data[idf1]);
    // }
    if (id_average.size() > 0)
        printf("averaging volumes:\n");
    for (int i = 0; i < id_average.size(); i++) {
        int idO0 = id_average[i][0];
        int idO1 = id_average[i][1];
        weighted_average(data[idO0], data[idO1]);
        printf("obs %d after averaging with %d volumes: %.12g  %.12g\n", idO0, idO1, data[idO0][Njack - 1], myres->comp_error(data[idO0]));
    }

    printf("Mpi after averaging volumes: %.12g  %.12g\n", data[0][Njack - 1], myres->comp_error(data[0]));
    printf("fpi after averaging volumes: %.12g  %.12g\n", data[3][Njack - 1], myres->comp_error(data[3]));

    double* Rpi = myres->create_copy(data[0]);
    myres->div(Rpi, Rpi, data[3]);
    printf("Rpi after averaging volumes: %.12g  %.12g\n", Rpi[Njack - 1], myres->comp_error(Rpi));

    //////////////////////////////////////////////////////////////
    // sistemone fpi
    //////////////////////////////////////////////////////////////

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
    for (int ir = 0;ir < 2;ir++) {
        double* M = data[id_deriv_fpi(ir, 0, 0, 0)];
        printf("fpi_r%d val=%-8.3g (%-8.3g)\t", ir, M[Njack - 1], myres->comp_error(M));
        // valence
        for (int im = 0; im < 3; im++) {
            double* dM = data[id_deriv_fpi(ir, 1, im, 0)];
            printf(" dval%d = %-8.3g (%-8.3g)", im, dM[Njack - 1], myres->comp_error(dM));
        }
        for (int im = 0; im < 3; im++) {
            double* dM = data[id_deriv_fpi(ir, 1, im, 1)];
            printf(" dsea%d = %-8.3g (%-8.3g)", im, dM[Njack - 1], myres->comp_error(dM));
        }
        printf("\n");
    }

    double* w0_lin_deriv = myres->create_copy(data[4]);

    double* Rds = myres->create_copy(data[id_deriv(2, 0, 0, 0)]);
    myres->div(Rds, Rds, data[id_deriv(3, 0, 0, 0)]);
    // printf("Rds = %.12g  %.12g\n", Rds[Njack-1], myres->comp_error(Rds));

    double** Mat = malloc_2<double>(3, 3);
    double*** Matj = malloc_3<double>(3, 3, Njack);
    double* y = (double*)malloc(sizeof(double) * 3);
    double** yj = malloc_2<double>(3, Njack);
    double** miso = malloc_2<double>(3, Njack);
    double** dm_fpi = malloc_2<double>(3, Njack);
    double* w0_from_fpi_ensemble = myres->create_copy(data[iw0]);
    double* w0_from_fpi_hybrid = myres->create_copy(data[iw0]);
    {
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
                    Matj[iM][im][j] = Mat[iM][im];
                }
                y[iM] *= y[iM];
                // if (j==Njack-1) printf(" %.12g   %g   %g   %g\n ",y[iM], M / f, f, M);
                y[iM] -= (M / f) * (M / f);
                yj[iM][j] = y[iM];

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
                    Matj[iM][im][j] = Mat[iM][im];

                }
                // if (j == Njack - 1) printf(" %.12g   %g   %g   %g\n ", y[iM], M / f, f, M);
                y[iM] -= M / f;
                yj[iM][j] = y[iM];

            }
            double* P = LU_decomposition_solver(3, Mat, y);
            miso[0][j] = (amusim[0][j] + P[0]);
            miso[1][j] = (amusim[1][j] + P[1]);
            miso[2][j] = (amusim[2][j] + P[2]);
            for (int im = 0; im < 3; im++) {
                dm_fpi[im][j] = P[im];
            }
            if (j == Njack - 1) {
                printf("Matrix for m^iso solution (jackknife %d):\n", j);
                for (int ii = 0; ii < 3; ii++) {
                    for (int jj = 0; jj < 3; jj++) {
                        printf("%g (%g) ", Matj[ii][jj][Njack - 1], myres->comp_error(Matj[ii][jj]));
                    }
                    printf("\n");
                }
                printf("RHS:\n");
                for (int ii = 0; ii < 3; ii++) {
                    printf("%g (%g)\n", yj[ii][Njack - 1], myres->comp_error(yj[ii]));
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

        // decorrelate miso
        // for (int i = 0; i < 3; i++) {
        //     double mean = myres->mean(miso[i]);
        //     double err = myres->comp_error(miso[i]);
        //     free(miso[i]);
        //     miso[i] = myres->create_fake(mean, err, -1);
        // }


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
                    // if (val_sea != 1 && im != 2)
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

        // w0 from fpi with strange and charm deriv ensemble
        for (int j = 0; j < Njack;j++) {
            double dw = data[id_deriv(iw0, 1, 0, 1)][j];
            double dm = (miso[0][j] - amusim[0][j]);
            w0_from_fpi_ensemble[j] += dm * dw;
            dw = data[39][j];
            dm = (miso[1][j] - amusim[1][j]);
            w0_from_fpi_ensemble[j] += dm * dw;
            dw = data[40][j];
            dm = (miso[2][j] - amusim[2][j]);
            w0_from_fpi_ensemble[j] += dm * dw;

            w0_from_fpi_ensemble[j] *= a_fm[j];
        }
        // hybrid approach
        for (int j = 0; j < Njack;j++) {
            double dw = data[id_deriv(iw0, 1, 0, 1)][j];
            double dm = (miso[0][j] - amusim[0][j]);
            w0_from_fpi_hybrid[j] += dm * dw;

            if (ensemble.compare("C80") == 0 || ensemble.compare("B64") == 0) {
                if (j == Njack - 1) printf("HYBRID APPROACH: USINGE SMALL VOLUME DERIV\n");
                dw = data[id_deriv(iw0, 1, 1, 1)][j]; // small volume deriv strange
                dm = (miso[1][j] - amusim[1][j]);
                w0_from_fpi_hybrid[j] += dm * dw;
                dw = data[id_deriv(iw0, 1, 2, 1)][j]; // small volume deriv charm
                dm = (miso[2][j] - amusim[2][j]);
                w0_from_fpi_hybrid[j] += dm * dw;

                w0_from_fpi_hybrid[j] *= a_fm[j];
            }
            else { // add the ensemble deriv
                dw = data[39][j];
                dm = (miso[1][j] - amusim[1][j]);
                w0_from_fpi_hybrid[j] += dm * dw;
                dw = data[40][j];
                dm = (miso[2][j] - amusim[2][j]);
                w0_from_fpi_hybrid[j] += dm * dw;

                w0_from_fpi_hybrid[j] *= a_fm[j];
            }
        }
        // linear deriv mc
        for (int j = 0; j < Njack;j++) {
            for (int im = 0; im < 3; im++) {
                int val_sea = 1;
                double dm = (miso[im][j] - amusim[im][j]);
                double dw = data[id_deriv(iw0, 1, im, val_sea)][j];
                if (im == 2 && val_sea == 1) {
                    if (j == Njack - 1) printf("replacing dw0/dmc = %g\n", dw);
                    dw = get_linear_deriv_w0c(data, amuiso, previous_a, j);
                    if (j == Njack - 1) printf("with      dw0/dmc = %g\n", dw);
                }

                w0_lin_deriv[j] += dm * dw;
            }
            w0_lin_deriv[j] *= a_fm[j];
        }

        double** data_m_a = malloc_2<double>(4, Njack);
        for (int j = 0; j < Njack;j++) {
            data_m_a[0][j] = miso[0][j];
            data_m_a[1][j] = miso[1][j];
            data_m_a[2][j] = miso[2][j];
            data_m_a[3][j] = a_fm[j];
        }
        double** cov_m_a = myres->comp_cov(4, data_m_a);
        printf("covariance matrix for (in order) m^iso l,s,c and a, ens: %s\n", ensemble.c_str());
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                // double corr = cov_m_a[i][j] / sqrt(cov_m_a[i][i] * cov_m_a[j][j]);
                // if (std::fabs(corr) > 0.4)
                printf("%-22.12g", cov_m_a[i][j]);
            }
            printf("\n");
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        printf("//////////////////////////////////////////////////////////////\n");
        printf("// w0 systemone\n");
        printf("//////////////////////////////////////////////////////////////\n");

        double** miso_w0 = malloc_2<double>(3, Njack);
        double** dm_w0 = malloc_2<double>(3, Njack);

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
                    // if (im!=2) 
                    Mat[iM][im] += dM * w0 + M * dw;
                }
                y[iM] -= M * w0;

            }
            double* P = LU_decomposition_solver(3, Mat, y);
            miso_w0[0][j] = (amusim[0][j] + P[0]);
            miso_w0[1][j] = (amusim[1][j] + P[1]);
            miso_w0[2][j] = (amusim[2][j] + P[2]);
            for (int im = 0; im < 3; im++) {
                dm_w0[im][j] = P[im];
            }
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
                    // if (val_sea != 2 && im != 2)
                    w_a += dm * dw;
                }
            }
            a_from_w0[j] = w0_fm / w_a;
            fpi_from_w0[j] = af / (a_from_w0[j] / hbarc);
        }
        printf("lattice spacing (fm): %g +/- %g\n", myres->mean(a_from_w0), myres->comp_error(a_from_w0));
        printf("fpi (fm): %g +/- %g\n", myres->mean(fpi_from_w0), myres->comp_error(fpi_from_w0));

        char name_out[NAMESIZE];
        mysprintf(name_out, NAMESIZE, "scale_setting/%s_a_from_w0.jack", files[36].c_str());
        myres->write_jack_in_file(a_from_w0, name_out);
        mysprintf(name_out, NAMESIZE, "scale_setting/%s_fpi_from_w0.jack", files[36].c_str());
        myres->write_jack_in_file(fpi_from_w0, name_out);

        write_jack(miso_w0[0], Njack, jack_file);     check_correlatro_counter(35);
        write_jack(miso_w0[1], Njack, jack_file);     check_correlatro_counter(36);
        write_jack(miso_w0[2], Njack, jack_file);     check_correlatro_counter(37);
        write_jack(a_from_w0, Njack, jack_file);     check_correlatro_counter(38);
        write_jack(fpi_from_w0, Njack, jack_file);     check_correlatro_counter(39);

        write_jack(dm_fpi[0], Njack, jack_file);     check_correlatro_counter(40);
        write_jack(dm_fpi[1], Njack, jack_file);     check_correlatro_counter(41);
        write_jack(dm_fpi[2], Njack, jack_file);     check_correlatro_counter(42);

        write_jack(dm_w0[0], Njack, jack_file);     check_correlatro_counter(43);
        write_jack(dm_w0[1], Njack, jack_file);     check_correlatro_counter(44);
        write_jack(dm_w0[2], Njack, jack_file);     check_correlatro_counter(45);


        write_jack(w0_from_fpi_ensemble, Njack, jack_file);     check_correlatro_counter(46);
        write_jack(w0_from_fpi_hybrid, Njack, jack_file);     check_correlatro_counter(47);

    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        printf("//////////////////////////////////////////////////////////////\n");
        printf("// w0 systemone hybrid\n");
        printf("//////////////////////////////////////////////////////////////\n");

        double** miso_w0_h = malloc_2<double>(3, Njack);
        double** dm_w0_h = malloc_2<double>(3, Njack);

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
                    dw = data[id_deriv(iw0, 1, im, 1)][j]; // small volume deriv strange
                    if (!(ensemble.compare("C80") == 0 || ensemble.compare("B64") == 0)) {
                        if (im == 1) dw = data[39][j]; // ensemble deriv    
                        if (im == 2) dw = data[40][j]; // ensemble deriv    
                    }
                    Mat[iM][im] += dM * w0 + M * dw;
                    Matj[iM][im][j] = Mat[iM][im];

                }
                y[iM] -= M * w0;
                yj[iM][j] = y[iM];

            }
            double* P = LU_decomposition_solver(3, Mat, y);
            miso_w0_h[0][j] = (amusim[0][j] + P[0]);
            miso_w0_h[1][j] = (amusim[1][j] + P[1]);
            miso_w0_h[2][j] = (amusim[2][j] + P[2]);
            for (int im = 0; im < 3; im++) {
                dm_w0_h[im][j] = P[im];
            }
            if (j == Njack - 1) {
                printf("Matrix for m^iso solution (jackknife %d):\n", j);
                for (int ii = 0; ii < 3; ii++) {
                    for (int jj = 0; jj < 3; jj++) {
                        printf("%g (%g) ", Matj[ii][jj][Njack - 1], myres->comp_error(Matj[ii][jj]));
                    }
                    printf("\n");
                }
                printf("RHS:\n");
                for (int ii = 0; ii < 3; ii++) {
                    printf("%g (%g)\n", yj[ii][Njack - 1], myres->comp_error(yj[ii]));
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
            double mean = myres->mean(miso_w0_h[i]);
            double err = myres->comp_error(miso_w0_h[i]);
            printf("m^iso_%d = %-12g +/- %-12g   starting from %-12g +/- %-12g\n", i, mean, err, myres->mean(amuiso[i]), myres->comp_error(amuiso[i]));
        }
        printf("sim values:\n");
        for (int i = 0; i < 3; i++) {
            printf("m^sim_%d = %g \n", i, myres->mean(amusim[i]));
        }

        double* a_from_w0_hybrid = myres->create_copy(data[iw0]);
        for (int j = 0; j < Njack;j++) {
            double w0_hybrid = data[iw0][j];
            double dw = data[id_deriv(iw0, 1, 0, 1)][j];
            double dm = (miso_w0_h[0][j] - amusim[0][j]);
            w0_hybrid += dm * dw;

            if (ensemble.compare("C80") == 0 || ensemble.compare("B64") == 0) {
                if (j == Njack - 1) printf("HYBRID APPROACH: USINGE SMALL VOLUME DERIV\n");
                dw = data[id_deriv(iw0, 1, 1, 1)][j]; // small volume deriv strange
                dm = (miso_w0_h[1][j] - amusim[1][j]);
                w0_hybrid += dm * dw;
                dw = data[id_deriv(iw0, 1, 2, 1)][j]; // small volume deriv charm
                dm = (miso_w0_h[2][j] - amusim[2][j]);
                w0_hybrid += dm * dw;

            }
            else { // add the ensemble deriv
                dw = data[39][j];
                dm = (miso_w0_h[1][j] - amusim[1][j]);
                w0_hybrid += dm * dw;
                dw = data[40][j];
                dm = (miso_w0_h[2][j] - amusim[2][j]);
                w0_hybrid += dm * dw;

            }
            a_from_w0_hybrid[j] = w0_fm / w0_hybrid;
        }

        double* a_from_Mpi_wp25_hybrid = myres->create_copy(data[iw0]);
        for (int j = 0; j < Njack;j++) {
            int iMpi = 0;
            double aMpi = data[id_deriv(iMpi, 0, 0, 0)][j];
            if (j == Njack - 1) printf("a from Mpi wp25 hybrid: starting value from Mpi %g\n", aMpi);
            for (int im = 0; im < 3; im++) {
                for (int val_sea = 0; val_sea < 2; val_sea++) {
                    double dM = data[id_deriv(iMpi, 1, im, val_sea)][j];
                    double dm = (miso_w0_h[im][j] - amusim[im][j]);
                    aMpi += dm * dM;
                }
            }
            a_from_Mpi_wp25_hybrid[j] = aMpi / (Mpi_MeV / hbarc);
        }


        double* fpi_from_w0_h = myres->create_copy(data[4]);
        double* fpi_from_Mpi_wp25_hybrid = myres->create_copy(data[4]);
        for (int j = 0; j < Njack;j++) {
            int ifpi = 3;
            double af = data[id_deriv(ifpi, 0, 0, 0)][j];
            int iw0 = 4;
            double w_a = data[id_deriv(iw0, 0, 0, 0)][j];
            for (int im = 0; im < 3; im++) {
                for (int val_sea = 0; val_sea < 2; val_sea++) {
                    double df = data[id_deriv(ifpi, 1, im, val_sea)][j];
                    double dm = (miso_w0_h[im][j] - amusim[im][j]);
                    af += dm * df;
                    double dw = data[id_deriv(iw0, 1, im, val_sea)][j];
                    // if (val_sea != 2 && im != 2)
                    w_a += dm * dw;
                }
            }
            fpi_from_w0_h[j] = af / (a_from_w0_hybrid[j] / hbarc);
            fpi_from_Mpi_wp25_hybrid[j] = af / (a_from_Mpi_wp25_hybrid[j] / hbarc);
        }
        printf("lattice spacing (fm): %g +/- %g\n", myres->mean(a_from_w0_hybrid), myres->comp_error(a_from_w0_hybrid));
        printf("fpi (fm): %g +/- %g\n", myres->mean(fpi_from_w0_h), myres->comp_error(fpi_from_w0_h));

        char name_out[NAMESIZE];
        mysprintf(name_out, NAMESIZE, "scale_setting/%s_a_from_w0_hybrid.jack", files[36].c_str());
        myres->write_jack_in_file(a_from_w0_hybrid, name_out);
        mysprintf(name_out, NAMESIZE, "scale_setting/%s_fpi_from_w0_hybrid.jack", files[36].c_str());
        myres->write_jack_in_file(fpi_from_w0_h, name_out);

        write_jack(miso_w0_h[0], Njack, jack_file);     check_correlatro_counter(48);
        write_jack(miso_w0_h[1], Njack, jack_file);     check_correlatro_counter(49);
        write_jack(miso_w0_h[2], Njack, jack_file);     check_correlatro_counter(50);
        write_jack(a_from_w0_hybrid, Njack, jack_file);     check_correlatro_counter(51);
        write_jack(fpi_from_w0_h, Njack, jack_file);     check_correlatro_counter(52);

        write_jack(dm_w0_h[0], Njack, jack_file);     check_correlatro_counter(53);
        write_jack(dm_w0_h[1], Njack, jack_file);     check_correlatro_counter(54);
        write_jack(dm_w0_h[2], Njack, jack_file);     check_correlatro_counter(55);

        mysprintf(name_out, NAMESIZE, "scale_setting/%s_dmu_l_wp25_hybrid.jack", files[36].c_str());
        myres->write_jack_in_file(dm_w0_h[0], name_out);
        mysprintf(name_out, NAMESIZE, "scale_setting/%s_dmu_s_wp25_hybrid.jack", files[36].c_str());
        myres->write_jack_in_file(dm_w0_h[1], name_out);
        mysprintf(name_out, NAMESIZE, "scale_setting/%s_dmu_c_wp25_hybrid.jack", files[36].c_str());
        myres->write_jack_in_file(dm_w0_h[2], name_out);


        printf("lattice spacing  wpMpi (fm): %g +/- %g\n", myres->mean(a_from_Mpi_wp25_hybrid), myres->comp_error(a_from_Mpi_wp25_hybrid));
        write_jack(a_from_Mpi_wp25_hybrid, Njack, jack_file);     check_correlatro_counter(56);
        write_jack(fpi_from_Mpi_wp25_hybrid, Njack, jack_file);     check_correlatro_counter(57);
    }

    //////////////////////////////////////////////////////////////
    // w0 system linear
    //////////////////////////////////////////////////////////////
    {
        printf("//////////////////////////////////////////////////////////////\n");
        printf("// w0 systemone linear deriv mc\n");
        printf("//////////////////////////////////////////////////////////////\n");

        double** miso_w0 = malloc_2<double>(3, Njack);
        double** dm_w0 = malloc_2<double>(3, Njack);

        for (int j = 0; j < Njack;j++) {

            y[0] = Mpi_MeV * w0_MeV;
            y[1] = MK_MeV * w0_MeV;
            y[2] = MDs_MeV * w0_MeV;
            static constexpr int iw0 = 4;
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
                    if (im == 2) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    // if (im!=2) 
                    Mat[iM][im] += dM * w0 + M * dw;
                }
                y[iM] -= M * w0;

            }
            double* P = LU_decomposition_solver(3, Mat, y);
            miso_w0[0][j] = (amusim[0][j] + P[0]);
            miso_w0[1][j] = (amusim[1][j] + P[1]);
            miso_w0[2][j] = (amusim[2][j] + P[2]);
            for (int im = 0; im < 3; im++) {
                dm_w0[im][j] = P[im];
            }
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
                    // add special case for charm sea
                    if (im == 2 && val_sea == 1) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    w_a += dm * dw;
                }
            }
            a_from_w0[j] = w0_fm / w_a;
            fpi_from_w0[j] = af / (a_from_w0[j] / hbarc);
        }
        printf("lattice spacing (fm): %g +/- %g\n", myres->mean(a_from_w0), myres->comp_error(a_from_w0));
        printf("fpi (fm): %g +/- %g\n", myres->mean(fpi_from_w0), myres->comp_error(fpi_from_w0));

        char name_out[NAMESIZE];
        mysprintf(name_out, NAMESIZE, "scale_setting/%s_a_from_w0.jack", files[36].c_str());
        myres->write_jack_in_file(a_from_w0, name_out);
        mysprintf(name_out, NAMESIZE, "scale_setting/%s_fpi_from_w0.jack", files[36].c_str());
        myres->write_jack_in_file(fpi_from_w0, name_out);

        write_jack(miso_w0[0], Njack, jack_file);     check_correlatro_counter(58);
        write_jack(miso_w0[1], Njack, jack_file);     check_correlatro_counter(59);
        write_jack(miso_w0[2], Njack, jack_file);     check_correlatro_counter(60);
        write_jack(a_from_w0, Njack, jack_file);     check_correlatro_counter(61);
        write_jack(fpi_from_w0, Njack, jack_file);     check_correlatro_counter(62);

        write_jack(dm_w0[0], Njack, jack_file);     check_correlatro_counter(63);
        write_jack(dm_w0[1], Njack, jack_file);     check_correlatro_counter(64);
        write_jack(dm_w0[2], Njack, jack_file);     check_correlatro_counter(65);

    }

    //////////////////////////////////////////////////////////////
    // w0 system linear plus lattice artefact RDs
    //////////////////////////////////////////////////////////////
    {
        printf("//////////////////////////////////////////////////////////////\n");
        printf("// w0 systemone linear deriv mc la RDs\n");
        printf("//////////////////////////////////////////////////////////////\n");

        double** miso_w0 = malloc_2<double>(3, Njack);
        double** dm_w0 = malloc_2<double>(3, Njack);
        double coeff = -5.0;
        for (int j = 0; j < Njack;j++) {

            y[0] = Mpi_MeV * w0_MeV;
            y[1] = MK_MeV * w0_MeV;
            y[2] = MDs_MeV * w0_MeV;
            static constexpr int iw0 = 4;
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
                    if (im == 2) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    // if (im!=2) 
                    Mat[iM][im] += dM * w0 + M * dw;
                }
                y[iM] -= M * w0;

            }
            // add the lattice artefact RDs
            y[2] += coeff * previous_a[j] * previous_a[j];

            double* P = LU_decomposition_solver(3, Mat, y);
            miso_w0[0][j] = (amusim[0][j] + P[0]);
            miso_w0[1][j] = (amusim[1][j] + P[1]);
            miso_w0[2][j] = (amusim[2][j] + P[2]);
            for (int im = 0; im < 3; im++) {
                dm_w0[im][j] = P[im];
            }
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
                    // add special case for charm sea
                    if (im == 2 && val_sea == 1) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    w_a += dm * dw;
                }
            }
            a_from_w0[j] = w0_fm / w_a;
            fpi_from_w0[j] = af / (a_from_w0[j] / hbarc);
        }
        printf("lattice spacing (fm): %g +/- %g\n", myres->mean(a_from_w0), myres->comp_error(a_from_w0));
        printf("fpi (fm): %g +/- %g\n", myres->mean(fpi_from_w0), myres->comp_error(fpi_from_w0));

        // char name_out[NAMESIZE];
        // mysprintf(name_out, NAMESIZE, "scale_setting/%s_a_from_w0.jack", files[36].c_str());
        // myres->write_jack_in_file(a_from_w0, name_out);
        // mysprintf(name_out, NAMESIZE, "scale_setting/%s_fpi_from_w0.jack", files[36].c_str());
        // myres->write_jack_in_file(fpi_from_w0, name_out);

        write_jack(miso_w0[0], Njack, jack_file);     check_correlatro_counter(66);
        write_jack(miso_w0[1], Njack, jack_file);     check_correlatro_counter(67);
        write_jack(miso_w0[2], Njack, jack_file);     check_correlatro_counter(68);
        write_jack(a_from_w0, Njack, jack_file);     check_correlatro_counter(69);
        write_jack(fpi_from_w0, Njack, jack_file);     check_correlatro_counter(70);


        write_jack(dm_w0[0], Njack, jack_file);     check_correlatro_counter(71);
        write_jack(dm_w0[1], Njack, jack_file);     check_correlatro_counter(72);
        write_jack(dm_w0[2], Njack, jack_file);     check_correlatro_counter(73);

    }

    // //////////////////////////////////////////////////////////////
    // // w0 system linear plus lattice artefact RDs
    // //////////////////////////////////////////////////////////////
    // double* a_old = myres->create_copy(previous_a);
    // double* a_new = myres->create_fake(0, 1e-20, -1);
    // int iter = 0;
    // while (std::fabs(a_new[Njack - 1] - a_old[Njack - 1]) > 0.0000001 * myres->comp_error(a_new)) {
    // // while (iter <1000){
    //     printf("###############################  iter  %d\n", iter++);
    //     printf("//////////////////////////////////////////////////////////////\n");
    //     printf("// w0 systemone linear deriv mc la RDs\n");
    //     printf("//////////////////////////////////////////////////////////////\n");

    //     double** miso_w0 = malloc_2<double>(3, Njack);
    //     double** dm_w0 = malloc_2<double>(3, Njack);
    //     double coeff = -5.0;
    //     for (int j = 0; j < Njack;j++) {

    //         y[0] = Mpi_MeV * w0_MeV;
    //         y[1] = MK_MeV * w0_MeV;
    //         y[2] = MDs_MeV * w0_MeV;
    //         static constexpr int iw0 = 4;
    //         double w0 = data[id_deriv(iw0, 0, 0, 0)][j];

    //         for (int iM = 0; iM < 3; iM++) {

    //             double M = data[id_deriv(iM, 0, 0, 0)][j];

    //             for (int im = 0; im < 3; im++) {
    //                 Mat[iM][im] = 0.0;
    //                 // valence
    //                 double dM = data[id_deriv(iM, 1, im, 0)][j];
    //                 double dw = data[id_deriv(iw0, 1, im, 0)][j];
    //                 Mat[iM][im] += dM * w0 + M * dw;
    //                 // sea
    //                 dM = data[id_deriv(iM, 1, im, 1)][j];
    //                 dw = data[id_deriv(iw0, 1, im, 1)][j];
    //                 if (im == 2) {
    //                     double P0 = data[41][j];
    //                     double P1 = data[42][j];
    //                     // remove the prefactor in dw0/dmc(sea)
    //                     double mul = amuiso[0][j];
    //                     double w0 = data[iw0][j];
    //                     double a = previous_a[j];
    //                     P0 *= w0 / mul;
    //                     P1 *= w0 / mul;
    //                     dw = P0 + P1 * a * a;
    //                 }
    //                 // if (im!=2) 
    //                 Mat[iM][im] += dM * w0 + M * dw;
    //             }
    //             y[iM] -= M * w0;

    //         }
    //         // add the lattice artefact RDs
    //         y[2] += coeff * a_old[j] * a_old[j];

    //         double* P = LU_decomposition_solver(3, Mat, y);
    //         miso_w0[0][j] = (amusim[0][j] + P[0]);
    //         miso_w0[1][j] = (amusim[1][j] + P[1]);
    //         miso_w0[2][j] = (amusim[2][j] + P[2]);
    //         for (int im = 0; im < 3; im++) {
    //             dm_w0[im][j] = P[im];
    //         }

    //         free(P);
    //     }
    //     printf("Results for m^iso (MeV):\n");
    //     for (int i = 0; i < 3; i++) {
    //         double mean = myres->mean(miso_w0[i]);
    //         double err = myres->comp_error(miso_w0[i]);
    //         printf("m^iso_%d = %-12g +/- %-12g   starting from %-12g +/- %-12g\n", i, mean, err, myres->mean(amuiso[i]), myres->comp_error(amuiso[i]));
    //     }
    //     printf("sim values:\n");
    //     for (int i = 0; i < 3; i++) {
    //         printf("m^sim_%d = %g \n", i, myres->mean(amusim[i]));
    //     }

    //     double* a_from_w0 = myres->create_copy(data[4]);
    //     double* fpi_from_w0 = myres->create_copy(data[4]);
    //     for (int j = 0; j < Njack;j++) {
    //         int ifpi = 3;
    //         double af = data[id_deriv(ifpi, 0, 0, 0)][j];
    //         int iw0 = 4;
    //         double w_a = data[id_deriv(iw0, 0, 0, 0)][j];
    //         for (int im = 0; im < 3; im++) {
    //             for (int val_sea = 0; val_sea < 2; val_sea++) {
    //                 double df = data[id_deriv(ifpi, 1, im, val_sea)][j];
    //                 double dm = (miso_w0[im][j] - amusim[im][j]);
    //                 af += dm * df;
    //                 double dw = data[id_deriv(iw0, 1, im, val_sea)][j];
    //                 // add special case for charm sea
    //                 if (im == 2 && val_sea == 1) {
    //                     double P0 = data[41][j];
    //                     double P1 = data[42][j];
    //                     // remove the prefactor in dw0/dmc(sea)
    //                     double mul = amuiso[0][j];
    //                     double w0 = data[iw0][j];
    //                     double a = previous_a[j];
    //                     P0 *= w0 / mul;
    //                     P1 *= w0 / mul;
    //                     dw = P0 + P1 * a * a;
    //                 }
    //                 w_a += dm * dw;
    //             }
    //         }
    //         a_from_w0[j] = w0_fm / w_a;
    //         fpi_from_w0[j] = af / (a_from_w0[j] / hbarc);
    //     }
    //     printf("lattice spacing (fm): %g +/- %g\n", myres->mean(a_from_w0), myres->comp_error(a_from_w0));
    //     printf("fpi (fm): %g +/- %g\n", myres->mean(fpi_from_w0), myres->comp_error(fpi_from_w0));

    //     myres->copy(a_old, a_new);
    //     myres->copy(a_new, a_from_w0);

    // }


    //////////////////////////////////////////////////////////////
    // w0 system linear plus lattice artefact RDs exact 
    //////////////////////////////////////////////////////////////
    {
        printf("//////////////////////////////////////////////////////////////\n");
        printf("// w0 systemone linear deriv mc la RDs exact\n");
        printf("//////////////////////////////////////////////////////////////\n");

        double** miso_w0 = malloc_2<double>(3, Njack);
        double** dm_w0 = malloc_2<double>(3, Njack);
        double coeff = -5.0;
        for (int j = 0; j < Njack;j++) {

            y[0] = Mpi_MeV * w0_MeV;
            y[1] = MK_MeV * w0_MeV;
            y[2] = MDs_MeV * w0_MeV;
            static constexpr int iw0 = 4;
            double w0 = data[id_deriv(iw0, 0, 0, 0)][j];
            // double a2_sim = previous_a[j] * previous_a[j];
            double a2_sim = std::pow(w0_fm / w0, 2);

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
                    if (im == 2) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    // if (im!=2) 
                    Mat[iM][im] += dM * w0 + M * dw;
                    //a2 coeff in RDs
                    if (iM == 2)
                        Mat[iM][im] += coeff * a2_sim * 2.0 * dw / w0;
                }
                y[iM] -= M * w0;

            }
            // add the lattice artefact RDs
            y[2] += coeff * a2_sim;

            double* P = LU_decomposition_solver(3, Mat, y);
            miso_w0[0][j] = (amusim[0][j] + P[0]);
            miso_w0[1][j] = (amusim[1][j] + P[1]);
            miso_w0[2][j] = (amusim[2][j] + P[2]);
            for (int im = 0; im < 3; im++) {
                dm_w0[im][j] = P[im];
            }
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
        int iMDs = 2;
        double* MDs = myres->create_copy(data[iMDs]);
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
                    // add special case for charm sea
                    if (im == 2 && val_sea == 1) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    w_a += dm * dw;

                    // MDs
                    double dM = data[id_deriv(iMDs, 1, im, val_sea)][j];
                    MDs[j] += dm * dM;

                }
            }
            a_from_w0[j] = w0_fm / w_a;
            fpi_from_w0[j] = af / (a_from_w0[j] / hbarc);
            MDs[j] /= (a_from_w0[j] / hbarc);
        }
        printf("lattice spacing (fm): %g +/- %g\n", myres->mean(a_from_w0), myres->comp_error(a_from_w0));
        printf("fpi (fm): %g +/- %g\n", myres->mean(fpi_from_w0), myres->comp_error(fpi_from_w0));

        // char name_out[NAMESIZE];
        // mysprintf(name_out, NAMESIZE, "scale_setting/%s_a_from_w0.jack", files[36].c_str());
        // myres->write_jack_in_file(a_from_w0, name_out);
        // mysprintf(name_out, NAMESIZE, "scale_setting/%s_fpi_from_w0.jack", files[36].c_str());
        // myres->write_jack_in_file(fpi_from_w0, name_out);

        write_jack(miso_w0[0], Njack, jack_file);     check_correlatro_counter(74);
        write_jack(miso_w0[1], Njack, jack_file);     check_correlatro_counter(75);
        write_jack(miso_w0[2], Njack, jack_file);     check_correlatro_counter(76);
        write_jack(a_from_w0, Njack, jack_file);     check_correlatro_counter(77);
        write_jack(fpi_from_w0, Njack, jack_file);     check_correlatro_counter(78);


        write_jack(dm_w0[0], Njack, jack_file);     check_correlatro_counter(79);
        write_jack(dm_w0[1], Njack, jack_file);     check_correlatro_counter(80);
        write_jack(dm_w0[2], Njack, jack_file);     check_correlatro_counter(81);

        write_jack(MDs, Njack, jack_file);     check_correlatro_counter(82);

    }

    //////////////////////////////////////////////////////////////
    // w0 system linear plus lattice artefact mc
    //////////////////////////////////////////////////////////////
    {
        printf("//////////////////////////////////////////////////////////////\n");
        printf("// w0 systemone linear deriv mc la mc\n");
        printf("//////////////////////////////////////////////////////////////\n");

        double** miso_w0 = malloc_2<double>(3, Njack);
        double** dm_w0 = malloc_2<double>(3, Njack);
        // double coeff = +0.75;
        double coeff = -260.70;

        for (int j = 0; j < Njack;j++) {

            y[0] = Mpi_MeV * w0_MeV;
            y[1] = MK_MeV * w0_MeV;
            y[2] = MDs_MeV * w0_MeV;
            double la[3] = { 0,0,0 };
            double rhs[3] = { 0,0,0 };

            static constexpr int iw0 = 4;
            double w0 = data[id_deriv(iw0, 0, 0, 0)][j];
            // double a2_sim = previous_a[j] * previous_a[j];
            double a3_sim = std::pow(w0_fm / w0, 3);

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
                    if (im == 2) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    // if (im!=2) 
                    Mat[iM][im] += dM * w0 + M * dw;
                    //a2 coeff in RDs
                    if (iM == 2)
                        Mat[iM][im] += coeff * a3_sim * 3.0 * dw / w0;
                }
                y[iM] -= M * w0;

            }
            // add the lattice artefact RDs
            la[2] += coeff * a3_sim;
            for (int ii = 0;ii < 3;ii++) {
                for (int ij = 0;ij < 3;ij++) {
                    rhs[ii] += Mat[ii][ij] * la[ij];
                }
                rhs[ii] += y[ii];
            }


            double* P = LU_decomposition_solver(3, Mat, y);
            miso_w0[0][j] = (amusim[0][j] + P[0]);
            miso_w0[1][j] = (amusim[1][j] + P[1]);
            miso_w0[2][j] = (amusim[2][j] + P[2]);
            for (int im = 0; im < 3; im++) {
                dm_w0[im][j] = P[im];
            }
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
        int iMDs = 2;
        double* MDs = myres->create_copy(data[iMDs]);
        for (int j = 0; j < Njack;j++) {
            int ifpi = 3;
            double af = data[id_deriv(ifpi, 0, 0, 0)][j];
            int iw0 = 4;
            double w_a = data[id_deriv(iw0, 0, 0, 0)][j];
            for (int im = 0; im < 3; im++) {
                for (int val_sea = 0; val_sea < 2; val_sea++) {
                    double dm = (miso_w0[im][j] - amusim[im][j]);
                    double dw = data[id_deriv(iw0, 1, im, val_sea)][j];
                    // add special case for charm sea
                    if (im == 2 && val_sea == 1) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    w_a += dm * dw;

                }
            }
            a_from_w0[j] = w0_fm / w_a;
            // // add lattice artefact to mc
            // miso_w0[2][j] += coeff * a_from_w0[j] * a_from_w0[j] * a_from_w0[j];
            // dm_w0[2][j] += coeff * a_from_w0[j] * a_from_w0[j] * a_from_w0[j];

            for (int im = 0; im < 3; im++) {
                for (int val_sea = 0; val_sea < 2; val_sea++) {
                    double dm = (miso_w0[im][j] - amusim[im][j]);
                    double df = data[id_deriv(ifpi, 1, im, val_sea)][j];
                    af += dm * df;

                    // MDs
                    double dM = data[id_deriv(iMDs, 1, im, val_sea)][j];
                    MDs[j] += dm * dM;

                }
            }

            fpi_from_w0[j] = af / (a_from_w0[j] / hbarc);
            MDs[j] /= (a_from_w0[j] / hbarc);
        }
        printf("lattice spacing (fm): %g +/- %g\n", myres->mean(a_from_w0), myres->comp_error(a_from_w0));
        printf("fpi (fm): %g +/- %g\n", myres->mean(fpi_from_w0), myres->comp_error(fpi_from_w0));


        write_jack(miso_w0[0], Njack, jack_file);     check_correlatro_counter(83);
        write_jack(miso_w0[1], Njack, jack_file);     check_correlatro_counter(84);
        write_jack(miso_w0[2], Njack, jack_file);     check_correlatro_counter(85);
        write_jack(a_from_w0, Njack, jack_file);     check_correlatro_counter(86);
        write_jack(fpi_from_w0, Njack, jack_file);     check_correlatro_counter(87);


        write_jack(dm_w0[0], Njack, jack_file);     check_correlatro_counter(88);
        write_jack(dm_w0[1], Njack, jack_file);     check_correlatro_counter(89);
        write_jack(dm_w0[2], Njack, jack_file);     check_correlatro_counter(90);

        write_jack(MDs, Njack, jack_file);     check_correlatro_counter(91);

    }

    //////////////////////////////////////////////////////////////
    // loop over the different coeff in the lattice artefact and see how the results change
    //////////////////////////////////////////////////////////////
    std::vector<double*> fpi_tm(coeffs.size());
    std::vector<double*> fpi_OS(coeffs.size());
    std::vector<double*> afpi_OS(coeffs.size());
    // a name for fpi_OS where i take the derivative from the WTI identity:
    std::vector<double*> fpi_OS_dWTI(coeffs.size());
    double** fpi_tm_split = malloc_2<double>(7, Njack);
    double** fpi_OS_split = malloc_2<double>(7, Njack);
    double** afpi_OS_split = malloc_2<double>(7, Njack);
    double** fpi_from_w0_split = malloc_2<double>(7, Njack);

    for (std::size_t ic = 0; ic < coeffs.size(); ++ic) {
        const double coeff = coeffs[ic];
        printf("//////////////////////////////////////////////////////////////\n");
        printf("// w0 systemone linear deriv mc la RDs  coef=%g\n", coeff);
        printf("//////////////////////////////////////////////////////////////\n");

        double** miso_w0 = malloc_2<double>(3, Njack);
        double** dm_w0 = malloc_2<double>(3, Njack);

        for (int j = 0; j < Njack;j++) {

            y[0] = Mpi_MeV * w0_MeV;
            y[1] = MK_MeV * w0_MeV;
            y[2] = MDs_MeV * w0_MeV;
            static constexpr int iw0 = 4;
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
                    if (im == 2) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    // if (im!=2) 
                    Mat[iM][im] += dM * w0 + M * dw;
                }
                y[iM] -= M * w0;

            }
            // add the lattice artefact RDs
            y[2] += coeff * previous_a[j] * previous_a[j];

            double* P = LU_decomposition_solver(3, Mat, y);
            miso_w0[0][j] = (amusim[0][j] + P[0]);
            miso_w0[1][j] = (amusim[1][j] + P[1]);
            miso_w0[2][j] = (amusim[2][j] + P[2]);
            for (int im = 0; im < 3; im++) {
                dm_w0[im][j] = P[im];
            }

            free(P);
        }
        printf("Results for m^iso (MeV):\n");
        for (int i = 0; i < 3; i++) {
            double mean = myres->mean(miso_w0[i]);
            double err = myres->comp_error(miso_w0[i]);
            printf("m^iso_%d = %-12g +/- %-12g   starting from %-12g +/- %-12g\n", i, mean, err, myres->mean(amuiso[i]), myres->comp_error(amuiso[i]));
        }


        double* a_from_w0 = myres->create_copy(data[4]);
        double* fpi_from_w0 = myres->create_copy(data[4]);
        fpi_tm[ic] = myres->create_copy(data[id_deriv_fpi(0, 0, 0, 0)]);
        fpi_OS[ic] = myres->create_copy(data[id_deriv_fpi(1, 0, 0, 0)]);
        afpi_OS[ic] = myres->create_copy(fpi_OS[ic]);
        fpi_OS_dWTI[ic] = myres->create_copy(data[id_deriv_fpi(1, 0, 0, 0)]);

        myres->copy(fpi_tm_split[0], data[id_deriv_fpi(0, 0, 0, 0)]);
        myres->copy(fpi_OS_split[0], data[id_deriv_fpi(1, 0, 0, 0)]);
        myres->copy(afpi_OS_split[0], fpi_OS_split[0]);
        myres->copy(fpi_from_w0_split[0], data[id_deriv(3, 0, 0, 0)]);
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
                    fpi_from_w0_split[1 + im + val_sea * 3][j] = dm * df;
                    double dw = data[id_deriv(iw0, 1, im, val_sea)][j];
                    // add special case for charm sea
                    if (im == 2 && val_sea == 1) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    w_a += dm * dw;

                    fpi_tm[ic][j] += dm * data[id_deriv_fpi(0, 1, im, val_sea)][j];
                    fpi_OS[ic][j] += dm * data[id_deriv_fpi(1, 1, im, val_sea)][j];
                    if (im == 0)
                        fpi_OS_dWTI[ic][j] += dm * data[id_deriv_fpi(1, 1, im, val_sea)][j];
                    else
                        fpi_OS_dWTI[ic][j] += dm * data[id_deriv(ifpi, 1, im, val_sea)][j];

                    fpi_tm_split[1 + im + val_sea * 3][j] = dm * data[id_deriv_fpi(0, 1, im, val_sea)][j];
                    fpi_OS_split[1 + im + val_sea * 3][j] = dm * data[id_deriv_fpi(1, 1, im, val_sea)][j];
                    afpi_OS_split[1 + im + val_sea * 3][j] = fpi_OS_split[1 + im + val_sea * 3][j];
                }
            }

            a_from_w0[j] = w0_fm / w_a;
            fpi_from_w0[j] = af / (a_from_w0[j] / hbarc);
            fpi_tm[ic][j] /= (a_from_w0[j] / hbarc);
            afpi_OS[ic][j] = fpi_OS[ic][j];
            fpi_OS[ic][j] /= (a_from_w0[j] / hbarc);
            fpi_OS_dWTI[ic][j] /= (a_from_w0[j] / hbarc);

            for (int ii = 0; ii < 7; ii++) {
                fpi_tm_split[ii][j] /= (a_from_w0[j] / hbarc);
                fpi_OS_split[ii][j] /= (a_from_w0[j] / hbarc);
                fpi_from_w0_split[ii][j] /= (a_from_w0[j] / hbarc);
            }
        }
        printf("lattice spacing (fm): %g +/- %g\n", myres->mean(a_from_w0), myres->comp_error(a_from_w0));
        printf("fpi (fm): %.6g +/- %.3g =", myres->mean(fpi_from_w0), myres->comp_error(fpi_from_w0));
        for (int ii = 0; ii < 7; ii++) {
            printf("%+.3g  (%.3g) ", myres->mean(fpi_from_w0_split[ii]), myres->comp_error(fpi_from_w0_split[ii]));
        }
        printf("\n");
        printf("fpi_tm(fm): %.6g +/- %.3g =", myres->mean(fpi_tm[ic]), myres->comp_error(fpi_tm[ic]));
        for (int ii = 0; ii < 7; ii++) {
            printf("%+.3g  (%.3g) ", myres->mean(fpi_tm_split[ii]), myres->comp_error(fpi_tm_split[ii]));
        }
        printf("\n");
        printf("fpi_OS(fm): %.6g +/- %.3g =", myres->mean(fpi_OS[ic]), myres->comp_error(fpi_OS[ic]));
        for (int ii = 0; ii < 7; ii++) {
            printf("%+.3g  (%.3g) ", myres->mean(fpi_OS_split[ii]), myres->comp_error(fpi_OS_split[ii]));
        }
        printf("\n");
        printf("afpi_OS(fm): %.6g +/- %.3g =", myres->mean(afpi_OS[ic]), myres->comp_error(afpi_OS[ic]));
        for (int ii = 0; ii < 7; ii++) {
            printf("%+.3g  (%.3g) ", myres->mean(afpi_OS_split[ii]), myres->comp_error(afpi_OS_split[ii]));
        }
        printf("\n");
        printf("fpi_OS_dWTI(fm): %.6g +/- %.3g \n", myres->mean(fpi_OS_dWTI[ic]), myres->comp_error(fpi_OS_dWTI[ic]));

        write_jack(miso_w0[0], Njack, jack_file);     check_correlatro_counter(92 + ic * 8);
        write_jack(miso_w0[1], Njack, jack_file);     check_correlatro_counter(93 + ic * 8);
        write_jack(miso_w0[2], Njack, jack_file);     check_correlatro_counter(94 + ic * 8);

        write_jack(a_from_w0, Njack, jack_file);     check_correlatro_counter(95 + ic * 8);
        write_jack(fpi_from_w0, Njack, jack_file);     check_correlatro_counter(96 + ic * 8);

        write_jack(dm_w0[0], Njack, jack_file);     check_correlatro_counter(97 + ic * 8);
        write_jack(dm_w0[1], Njack, jack_file);     check_correlatro_counter(98 + ic * 8);
        write_jack(dm_w0[2], Njack, jack_file);     check_correlatro_counter(99 + ic * 8);

    }

    std::vector<double*> fpi_mc_tm(coeffs_mc.size());
    std::vector<double*> fpi_mc_OS(coeffs_mc.size());

    //////////////////////////////////////////////////////////////
    // w0 system linear plus lattice artefact mc
    //////////////////////////////////////////////////////////////
    for (std::size_t ic = 0; ic < coeffs_mc.size(); ++ic) {
        const double coeff = coeffs_mc[ic];
        printf("//////////////////////////////////////////////////////////////\n");
        printf("// w0 systemone linear deriv mc la mc  coef=%g\n", coeff);
        printf("//////////////////////////////////////////////////////////////\n");

        double** miso_w0 = malloc_2<double>(3, Njack);
        double** dm_w0 = malloc_2<double>(3, Njack);
        // double coeff = +0.75;

        for (int j = 0; j < Njack;j++) {

            y[0] = Mpi_MeV * w0_MeV;
            y[1] = MK_MeV * w0_MeV;
            y[2] = MDs_MeV * w0_MeV;
            double la[3] = { 0,0,0 };
            double rhs[3] = { 0,0,0 };

            static constexpr int iw0 = 4;
            double w0 = data[id_deriv(iw0, 0, 0, 0)][j];
            // double a2_sim = previous_a[j] * previous_a[j];
            double a3_sim = std::pow(w0_fm / w0, 3);

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
                    if (im == 2) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    // if (im!=2) 
                    Mat[iM][im] += dM * w0 + M * dw;
                    //a2 coeff in RDs
                    if (iM == 2)
                        Mat[iM][im] += coeff * a3_sim * 3.0 * dw / w0;
                }
                y[iM] -= M * w0;

            }
            // add the lattice artefact RDs
            la[2] += coeff * a3_sim;
            for (int ii = 0;ii < 3;ii++) {
                for (int ij = 0;ij < 3;ij++) {
                    rhs[ii] += Mat[ii][ij] * la[ij];
                }
                rhs[ii] += y[ii];
            }


            double* P = LU_decomposition_solver(3, Mat, y);
            miso_w0[0][j] = (amusim[0][j] + P[0]);
            miso_w0[1][j] = (amusim[1][j] + P[1]);
            miso_w0[2][j] = (amusim[2][j] + P[2]);
            for (int im = 0; im < 3; im++) {
                dm_w0[im][j] = P[im];
            }
            free(P);
        }
        printf("Results for m^iso (MeV):\n");
        for (int i = 0; i < 3; i++) {
            double mean = myres->mean(miso_w0[i]);
            double err = myres->comp_error(miso_w0[i]);
            printf("m^iso_%d = %-12g +/- %-12g   starting from %-12g +/- %-12g\n", i, mean, err, myres->mean(amuiso[i]), myres->comp_error(amuiso[i]));
        }


        double* a_from_w0 = myres->create_copy(data[4]);
        double* fpi_from_w0 = myres->create_copy(data[4]);
        int iMDs = 2;
        double* MDs = myres->create_copy(data[iMDs]);
        fpi_mc_tm[ic] = myres->create_copy(data[id_deriv_fpi(0, 0, 0, 0)]);
        fpi_mc_OS[ic] = myres->create_copy(data[id_deriv_fpi(1, 0, 0, 0)]);
        for (int j = 0; j < Njack;j++) {
            int ifpi = 3;
            double af = data[id_deriv(ifpi, 0, 0, 0)][j];
            int iw0 = 4;
            double w_a = data[id_deriv(iw0, 0, 0, 0)][j];
            for (int im = 0; im < 3; im++) {
                for (int val_sea = 0; val_sea < 2; val_sea++) {
                    double dm = (miso_w0[im][j] - amusim[im][j]);
                    double dw = data[id_deriv(iw0, 1, im, val_sea)][j];
                    // add special case for charm sea
                    if (im == 2 && val_sea == 1) {
                        double P0 = data[41][j];
                        double P1 = data[42][j];
                        // remove the prefactor in dw0/dmc(sea)
                        double mul = amuiso[0][j];
                        double w0 = data[iw0][j];
                        double a = previous_a[j];
                        P0 *= w0 / mul;
                        P1 *= w0 / mul;
                        dw = P0 + P1 * a * a;
                    }
                    w_a += dm * dw;

                }
            }
            a_from_w0[j] = w0_fm / w_a;

            // // add lattice artefact to mc
            // miso_w0[2][j] += coeff * a_from_w0[j] * a_from_w0[j] * a_from_w0[j];
            // dm_w0[2][j] += coeff * a_from_w0[j] * a_from_w0[j] * a_from_w0[j];

            for (int im = 0; im < 3; im++) {
                for (int val_sea = 0; val_sea < 2; val_sea++) {
                    double dm = (miso_w0[im][j] - amusim[im][j]);
                    double df = data[id_deriv(ifpi, 1, im, val_sea)][j];
                    af += dm * df;

                    // MDs
                    double dM = data[id_deriv(iMDs, 1, im, val_sea)][j];
                    MDs[j] += dm * dM;

                    fpi_mc_tm[ic][j] += dm * data[id_deriv_fpi(0, 1, im, val_sea)][j];
                    fpi_mc_OS[ic][j] += dm * data[id_deriv_fpi(1, 1, im, val_sea)][j];
                }
            }

            fpi_from_w0[j] = af / (a_from_w0[j] / hbarc);
            MDs[j] /= (a_from_w0[j] / hbarc);
            fpi_mc_tm[ic][j] /= (a_from_w0[j] / hbarc);
            fpi_mc_OS[ic][j] /= (a_from_w0[j] / hbarc);
        }
        printf("lattice spacing (fm): %g +/- %g\n", myres->mean(a_from_w0), myres->comp_error(a_from_w0));
        printf("fpi (fm): %g +/- %g\n", myres->mean(fpi_from_w0), myres->comp_error(fpi_from_w0));
        printf("fpi_tm(fm): %g +/- %g\n", myres->mean(fpi_mc_tm[ic]), myres->comp_error(fpi_mc_tm[ic]));
        printf("fpi_OS(fm): %g +/- %g\n", myres->mean(fpi_mc_OS[ic]), myres->comp_error(fpi_mc_OS[ic]));


        write_jack(miso_w0[0], Njack, jack_file);     check_correlatro_counter(92 + coeffs.size() * 8 + ic * 9);
        write_jack(miso_w0[1], Njack, jack_file);     check_correlatro_counter(93 + coeffs.size() * 8 + ic * 9);
        write_jack(miso_w0[2], Njack, jack_file);     check_correlatro_counter(94 + coeffs.size() * 8 + ic * 9);
        write_jack(a_from_w0, Njack, jack_file);     check_correlatro_counter(95 + coeffs.size() * 8 + ic * 9);
        write_jack(fpi_from_w0, Njack, jack_file);     check_correlatro_counter(96 + coeffs.size() * 8 + ic * 9);


        write_jack(dm_w0[0], Njack, jack_file);     check_correlatro_counter(97 + coeffs.size() * 8 + ic * 9);
        write_jack(dm_w0[1], Njack, jack_file);     check_correlatro_counter(98 + coeffs.size() * 8 + ic * 9);
        write_jack(dm_w0[2], Njack, jack_file);     check_correlatro_counter(99 + coeffs.size() * 8 + ic * 9);

        write_jack(MDs, Njack, jack_file);     check_correlatro_counter(100 + coeffs.size() * 8 + ic * 9);

    }
    for (int ic = 0; ic < coeffs.size(); ic++) {
        write_jack(fpi_tm[ic], Njack, jack_file);     check_correlatro_counter(sid_fpi_A0 + ic * 2);
        write_jack(fpi_OS[ic], Njack, jack_file);     check_correlatro_counter(sid_fpi_A0 + ic * 2 + 1);
    }
    for (int ic = 0; ic < coeffs_mc.size(); ic++) {
        write_jack(fpi_mc_tm[ic], Njack, jack_file);     check_correlatro_counter(sid_fpi_A0 + coeffs.size() * 2 + ic * 2);
        write_jack(fpi_mc_OS[ic], Njack, jack_file);     check_correlatro_counter(sid_fpi_A0 + coeffs.size() * 2 + ic * 2 + 1);
    }
    for (int ic = 0; ic < coeffs.size(); ic++) {
        write_jack(fpi_OS_dWTI[ic], Njack, jack_file);     check_correlatro_counter(id_fpi_OS_dWTI + ic);
    }

    write_jack(w0_lin_deriv, Njack, jack_file);     check_correlatro_counter(id_w0_lin_deriv);
    return 0;
}
