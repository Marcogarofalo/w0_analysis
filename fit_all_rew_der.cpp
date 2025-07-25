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

    double r = gjack.en[e].jack[fit_info.corr_id[0]][j]; // fpi
    r *= gjack.en[e].jack[fit_info.corr_id[1]][j]; // amu_l
    r /= gjack.en[e].jack[fit_info.corr_id[2]][j]/hbarc; // a
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
double rhs_1overamu2(int n, int Nvar, double *x, int Npar, double *P)
{
    double amu = x[0];
    double r = P[0] / (amu * amu);
    return r;
}

double rhs_1overamu_1overamu2_1overamu3(int n, int Nvar, double *x, int Npar, double *P)
{
    double amu = x[0];
    double r = P[0] / amu + P[1] / (amu * amu) + P[2] / (amu * amu * amu);
    return r;
}
struct file_out_name
{
    char path[NAMESIZE];
    char basename[NAMESIZE];
    char namefile[NAMESIZE];
    file_out_name(const char *path, const char *label)
    {
        mysprintf(this->path, NAMESIZE, "%s", path);
        mysprintf(this->basename, NAMESIZE, "%s", path);
        mysprintf(this->namefile, NAMESIZE, "%s/%s_fit_extra.txt", path, label);
    }
    file_out_name(const file_out_name &other)
    {
        mysprintf(this->path, NAMESIZE, "%s", other.path);
        mysprintf(this->basename, NAMESIZE, "%s", other.basename);
        mysprintf(this->namefile, NAMESIZE, "%s", other.namefile);
    }
    file_out_name &operator=(const file_out_name &other)
    {
        mysprintf(this->path, NAMESIZE, "%s", other.path);
        mysprintf(this->basename, NAMESIZE, "%s", other.basename);
        mysprintf(this->namefile, NAMESIZE, "%s", other.namefile);
        return *this;
    }

    ~file_out_name()
    {
        // free(this->namefile);
        // free(this->path);
        // free(this->basename);
    }
};

void compute_fpi_at_mciso(data_all jackall, fit_type fit_info, fit_result fit_res, std::vector<double *> amusim, std::vector<double *> amuiso, int id_dfpi, file_out_name f_name)
{
    int Njack = jackall.en[0].Njack;
    error(fit_info.Njack != jackall.en[0].Njack, 1, "compute_fpi_at_mciso", "fit_info.Njack != jackall.en[0].Njack");
    error(fit_info.Njack != myres->Njack, 1, "compute_fpi_at_mciso", "fit_info.Njack != myres->Njack");
    double *fpi_mciso = (double *)malloc(sizeof(double) * Njack);
    double *relative_error = (double *)malloc(sizeof(double) * Njack);
    double *tmp_x = (double *)malloc(sizeof(double) * fit_info.Nvar);
    double *tif = (double *)malloc(sizeof(double) * fit_info.Npar);
    for (int j = 0; j < Njack; j++)
    {
        tmp_x[0] = amuiso[2][j]/amuiso[0][j]; // amuciso/amuliso
        for (int p = 0; p < fit_info.Npar; p++)
            tif[p] = fit_res.P[p][j];
        double der = fit_info.function(0, fit_info.Nvar, tmp_x, fit_info.Npar, tif)* (jackall.en[0].jack[16][j] /hbarc) /amuiso[0][j];
        fpi_mciso[j] = jackall.en[0].jack[id_dfpi][j] + (amuiso[2][j] - amusim[2][j]) * der;
        relative_error[j] = der * (amuiso[2][j] - amusim[2][j]) / jackall.en[0].jack[id_dfpi][j];
    }
    free(tmp_x);
    free(tif);
    printf("fpi(mc_sim) %g  %g\n", jackall.en[0].jack[id_dfpi][Njack - 1], myres->comp_error(jackall.en[0].jack[id_dfpi]));
    printf("fpi(mc_iso) %g  %g\n", fpi_mciso[Njack - 1], myres->comp_error(fpi_mciso));

    printf("writing in file %s\n", f_name.namefile);
    FILE *f = open_file(f_name.namefile, "w+");
    fprintf(f, "%-20.12g %-20.12g\n", jackall.en[0].jack[id_dfpi][Njack - 1], myres->comp_error(jackall.en[0].jack[id_dfpi]));
    fprintf(f, "%-20.12g %-20.12g\n", fpi_mciso[Njack - 1], myres->comp_error(fpi_mciso));
    fprintf(f, "%-20.12g %-20.12g\n", relative_error[Njack - 1], myres->comp_error(relative_error));
    fclose(f);
}

int main(int argc, char **argv)
{
    error(argc != 5, 1, "main ",
          "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir file_inputs");
    char namefile[NAMESIZE];

    std::vector<std::string> files;
    std::ifstream file(argv[4]);
    if (!file.is_open())
    {
        error(1, 1, "main", "Could not open file with input files: %s", argv[4]);
    }
    std::string line;
    std::string basename;
    while (std::getline(file, line))
    {
        if (!line.empty() && line[0] != '#') // skip empty lines and comments
        {

            std::vector<std::string> word = split(line, ' ');

            if (word.size() == 2)
            {
                basename = word[0];
                mysprintf(namefile, NAMESIZE, "%s/%s_%s_%s", argv[2], argv[1], word[0].c_str(), word[1].c_str());
                printf("adding file %s\n", namefile);
                files.emplace_back(namefile);
            }
            else
            {
                error(1, 1, "main", "Invalid line in input file: %s", line.c_str());
            }
        }
    }

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

    // read miso
    //////////////////////////////////////////////////////////////
    //  read m^iso
    //////////////////////////////////////////////////////////////
    char **option_read;
    option_read = (char **)malloc(sizeof(char *) * 7);
    option_read[0] = (char *)malloc(sizeof(char) * NAMESIZE);
    option_read[1] = (char *)malloc(sizeof(char) * NAMESIZE);
    option_read[2] = (char *)malloc(sizeof(char) * NAMESIZE);
    option_read[3] = (char *)malloc(sizeof(char) * NAMESIZE);
    option_read[4] = (char *)malloc(sizeof(char) * NAMESIZE);
    option_read[5] = (char *)malloc(sizeof(char) * NAMESIZE);
    option_read[6] = (char *)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option_read[1], NAMESIZE, "read_plateaux");        // blind/see/read_plateaux
    mysprintf(option_read[2], NAMESIZE, "-p");                   // -p
    mysprintf(option_read[3], NAMESIZE, "../../data/");          // path
    mysprintf(option_read[4], NAMESIZE, argv[1]);                // resampling
    mysprintf(option_read[5], NAMESIZE, "no");                   // pdf
    mysprintf(option_read[6], NAMESIZE, "%s", basename.c_str()); // infile

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");
    double mean, err;
    int seed;
    std::vector<double *> amuiso(3);
    line_read_param(option_read, "muliso", mean, err, seed, namefile_plateaux);
    amuiso[0] = myres->create_fake(mean, err, seed);
    line_read_param(option_read, "musiso", mean, err, seed, namefile_plateaux);
    amuiso[1] = myres->create_fake(mean, err, seed);
    line_read_param(option_read, "muciso", mean, err, seed, namefile_plateaux);
    amuiso[2] = myres->create_fake(mean, err, seed);

    std::vector<double *> amusim(3);
    line_read_param(option_read, "mulsim", mean, err, seed, namefile_plateaux);
    amusim[0] = myres->create_fake(mean, err, seed);
    line_read_param(option_read, "mussim", mean, err, seed, namefile_plateaux);
    amusim[1] = myres->create_fake(mean, err, seed);
    line_read_param(option_read, "mucsim", mean, err, seed, namefile_plateaux);
    amusim[2] = myres->create_fake(mean, err, seed);
    //
    int id_dfpi_dmu = 9;
    int id_dMpi_dmu = 4;
    int id_dfpi = 6;
    //////////////////////////////////////////////////////////////
    // fitting
    //////////////////////////////////////////////////////////////
    std::vector<std::string> der_name = {"fpi", "Mpi"};
    std::vector<int> id_der = {id_dfpi_dmu, id_dMpi_dmu};
    int id_amuliso=13;
    int id_a_fm=16;
    for (int i = 0; i < 2; i++)
    {
        fit_type fit_info;

        fit_info.corr_id = {id_der[i], id_amuliso, id_a_fm};

        fit_info.Nxen.resize(1, std::vector<int>(files.size()));
        for (int e = 0; e < files.size(); e++)
        {
            fit_info.Nxen[0][e] = e;
        }
        fit_info.N = fit_info.Nxen.size();
        fit_info.Npar = 1;
        fit_info.function = rhs_1overamu;

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
                    // printf(" %g  %d  %d\n", jackall.en[e].jack[11][j], e, j);
                    fit_info.x[0][count][j] = jackall.en[e].jack[11][j]/jackall.en[e].jack[id_amuliso][j];// mu/muliso
                }
                count++;
            }
        }
        fit_info.linear_fit = false;
        // fit_info.verbosity = 1;
        fit_info.covariancey = true;
        fit_info.compute_cov_fit(argv, jackall, lhs_fun);
        fit_info.compute_cov1_fit();
        std::string namefit("der_"+der_name[i]+"_vs_mu");

        fit_result der_fpi = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
        fit_info.band_range = {0.1, 350};
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "amu", der_fpi, der_fpi, 0, fit_info.Nxen[0].size() - 1, 1, {});

        // compute fpi at the charm point
        file_out_name f_name(argv[3], namefit.c_str());
        compute_fpi_at_mciso(jackall, fit_info, der_fpi, amusim, amuiso, id_dfpi, f_name);

        free_fit_result(fit_info, der_fpi);

        // add 1/amu^2
        fit_info.Npar = 2;
        fit_info.function = rhs_1overamu_1overamu2;
        namefit = "der_"+der_name[i]+"_vs_mu_mu2";

        fit_result der_fpi_mu_mu2 = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "amu", der_fpi_mu_mu2, der_fpi_mu_mu2, 0, fit_info.Nxen[0].size() - 1, 1, {});
        f_name = file_out_name(argv[3], namefit.c_str());
        compute_fpi_at_mciso(jackall, fit_info, der_fpi_mu_mu2, amusim, amuiso, id_dfpi, f_name);
        free_fit_result(fit_info, der_fpi_mu_mu2);

        // only 1/amu^2
        fit_info.Npar = 1;
        fit_info.function = rhs_1overamu2;
        namefit = "der_"+der_name[i]+"_vs_mu2";

        fit_result der_fpi_mu2 = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "amu", der_fpi_mu2, der_fpi_mu2, 0, fit_info.Nxen[0].size() - 1, 1, {});
        f_name = file_out_name(argv[3], namefit.c_str());
        compute_fpi_at_mciso(jackall, fit_info, der_fpi_mu2, amusim, amuiso, id_dfpi, f_name);
        free_fit_result(fit_info, der_fpi_mu2);

        // add 1/amu^2 + 1/amu^3
        fit_info.Npar = 3;
        fit_info.function = rhs_1overamu_1overamu2_1overamu3;
        namefit = "der_"+der_name[i]+"_vs_mu_mu2_mu3";

        fit_result der_fpi_mu_mu2_mu3 = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "amu", der_fpi_mu_mu2_mu3, der_fpi_mu_mu2_mu3, 0, fit_info.Nxen[0].size() - 1, 1, {});
        f_name = file_out_name(argv[3], namefit.c_str());
        compute_fpi_at_mciso(jackall, fit_info, der_fpi_mu_mu2_mu3, amusim, amuiso, id_dfpi, f_name);
        free_fit_result(fit_info, der_fpi_mu_mu2_mu3);

        // add 1/amu^2
        fit_info.Nxen.resize(1, std::vector<int>(files.size()));
        for (int e = 0; e < files.size() - 1; e++)
        {
            fit_info.Nxen[0][e] = e;
        }
        fit_info.Npar = 2;
        fit_info.function = rhs_1overamu_1overamu2;
        fit_info.compute_cov_fit(argv, jackall, lhs_fun);
        fit_info.compute_cov1_fit();
        namefit = "der_"+der_name[i]+"_0123_vs_mu_mu2";

        fit_result der_fpi_0123_mu_mu2 = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "amu", der_fpi_0123_mu_mu2, der_fpi_0123_mu_mu2, 0, fit_info.Nxen[0].size() - 1, 1, {});
        f_name = file_out_name(argv[3], namefit.c_str());
        compute_fpi_at_mciso(jackall, fit_info, der_fpi_0123_mu_mu2, amusim, amuiso, id_dfpi, f_name);
        free_fit_result(fit_info, der_fpi_0123_mu_mu2);

        // add 1/amu^2
        fit_info.Nxen.resize(1, std::vector<int>(files.size()));
        for (int e = 0; e < files.size() - 2; e++)
        {
            fit_info.Nxen[0][e] = e;
        }
        fit_info.Npar = 2;
        fit_info.function = rhs_1overamu_1overamu2;
        namefit = "der_"+der_name[i]+"_012_vs_mu_mu2";

        fit_result der_fpi_012_mu_mu2 = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "amu", der_fpi_012_mu_mu2, der_fpi_012_mu_mu2, 0, fit_info.Nxen[0].size() - 1, 1, {});
        f_name = file_out_name(argv[3], namefit.c_str());
        compute_fpi_at_mciso(jackall, fit_info, der_fpi_012_mu_mu2, amusim, amuiso, id_dfpi, f_name);
        free_fit_result(fit_info, der_fpi_012_mu_mu2);
    }
}
