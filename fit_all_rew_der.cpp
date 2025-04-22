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

#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

enum enum_ensembles
{
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
constexpr double Mpi_MeV = 135;
constexpr double Mpi_MeV_err = 0.2;
// #include "do_analysis_charm.hpp"

generic_header read_header(FILE *stream)
{
    generic_header header;
    int ir = 0;
    ir += fread(&header.T, sizeof(int), 1, stream);
    ir += fread(&header.L, sizeof(int), 1, stream);
    int s;
    ir += fread(&s, sizeof(int), 1, stream);
    header.mus = std::vector<double>(s);
    for (int i = 0; i < s; i++)
    {
        ir += fread(&header.mus[i], sizeof(double), 1, stream);
    }
    ir += fread(&s, sizeof(int), 1, stream);
    header.thetas = std::vector<double>(s);
    for (int i = 0; i < s; i++)
    {
        ir += fread(&header.thetas[i], sizeof(double), 1, stream);
    }

    ir += fread(&header.Njack, sizeof(int), 1, stream);
    header.struct_size = ftell(stream);
    return header;
}

double read_single_Nobs(FILE *stream, int header_size, int Njack)
{
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

data_single read_single_dataj(FILE *stream)
{

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
    for (int obs = 0; obs < dj.Nobs; obs++)
    {
        i += fread(dj.jack[obs], sizeof(double), dj.Njack, stream);
    }
    return dj;
}

data_all read_all_the_files(std::vector<std::string> files, const char *resampling)
{
    data_all jackall;
    jackall.resampling = resampling;
    // jackall->en = (data_single*)malloc(sizeof(data_single) * files.size());
    jackall.en = new data_single[files.size()];
    jackall.ens = files.size();
    int count = 0;
    for (std::string s : files)
    {
        printf("reading %s\n", s.c_str());
        FILE *f = open_file(s.c_str(), "r");

        // read_single_dataj(f, params, &(jackall->en[count]));
        jackall.en[count] = read_single_dataj(f);
        jackall.en[count].resampling = resampling;
        count++;
        fclose(f);
    }
    return jackall;
}

double lhs_fun(int n, int e, int j, data_all gjack, struct fit_type fit_info)
{

    double r = gjack.en[e].jack[fit_info.corr_id[0]][j];
    return r;
}
double rhs_1overamu(int n, int Nvar, double *x, int Npar, double *P)
{
    double amu = x[0];
    double r = P[0] / amu;

    return r;
}

double rhs_1overamu_1overamu2(int n, int Nvar, double *x, int Npar, double *P)
{
    double amu = x[0];
    double r = P[0] / amu + P[1] / (amu * amu);

    return r;
}


int main(int argc, char **argv)
{
    error(argc != 4, 1, "main ",
          "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");
    char namefile[NAMESIZE];

    std::vector<std::string> files;
    mysprintf(namefile, NAMESIZE, "%s/%s_onlinemeas_B64.dat_reweight_charm_0.1_OS_B64.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_onlinemeas_B64.dat_reweight_step_0.01_to_0.01027408_OS_B64.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_onlinemeas_B64.dat_reweight_step_0.03125_to_0.03152408_OS_B64.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_onlinemeas_B64.dat_reweight_step_0.1232229005_to_0.1234969805_OS_B64.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_onlinemeas_B64.dat_reweight_strange_OS_B64.dat", argv[2], argv[1]);
    files.emplace_back(namefile);

    std::vector<int> myen(files.size());
    for (int e = 0; e < files.size(); e++)
    {
        myen[e] = e;
    }

    data_all jackall = read_all_the_files(files, argv[1]);
    printf("we read all\n");

    jackall.create_generalised_resampling();
    int Njack = jackall.en[0].Njack;
    if (strcmp(argv[1], "jack") == 0)
    {
        myres = new resampling_jack(Njack - 1);
    }
    else if (strcmp(argv[1], "boot") == 0)
    {
        myres = new resampling_boot(Njack - 1);
    }

    //
    int id_dfpi_dmu = 9;
    //////////////////////////////////////////////////////////////
    // fitting
    //////////////////////////////////////////////////////////////
    fit_type fit_info;

    fit_info.corr_id = {id_dfpi_dmu};
    fit_info.Nxen = { {0, 1, 2, 3, 4} };
     fit_info.N = fit_info.Nxen.size();
    fit_info.Npar = 1;
    fit_info.function = rhs_1overamu ;

    fit_info.Nvar = 1;
    fit_info.Njack = Njack;
    fit_info.init_N_etot_form_Nxen();
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);

    int count = 0;
    for (int n = 0; n < fit_info.N; n++)
    {
        for (int e : fit_info.Nxen[n])
        {
            for (int j = 0; j < Njack; j++)
            {
                printf(" %g  %d  %d\n", jackall.en[e].jack[11][j],e,j);
                fit_info.x[0][count][j] = jackall.en[e].jack[11][j];
            }
            count++;
        }
    }
    fit_info.linear_fit = true;
    fit_info.verbosity = 0;
    fit_info.covariancey = false;
    std::string namefit( "der_fpi_vs_mu");

    fit_result der_fpi = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
    fit_info.band_range = {0.005, 0.4};
    print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "amu", der_fpi, der_fpi, 0, fit_info.Nxen[0].size() - 1, 0.002, {});

    free_fit_result(fit_info, der_fpi);

    fit_info.Npar = 2;
    fit_info.function = rhs_1overamu_1overamu2;
    namefit= "der_fpi_vs_mu_mu2";

    fit_result der_fpi_mu_mu2 = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
    print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "amu", der_fpi_mu_mu2, der_fpi_mu_mu2, 0, fit_info.Nxen[0].size() - 1, 0.002, {});

    free_fit_result(fit_info, der_fpi_mu_mu2);


}
