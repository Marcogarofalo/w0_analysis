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
    // if (j==fit_info.Njack-1) {
    //     printf("lhs_fun: n=%d e=%d j=%d fit_info.corr_id[0]=%d fit_info.corr_id[1]=%d fit_info.corr_id[2]=%d\n",
    //            n, e, j, fit_info.corr_id[0], fit_info.corr_id[1], fit_info.corr_id[2]);
    //     printf("lhs_fun: gjack.en[e].jack[fit_info.corr_id[0]][j]=%g\n", gjack.en[e].jack[fit_info.corr_id[0]][j]);
    //     printf("lhs_fun: gjack.en[e].jack[fit_info.corr_id[1]][j]=%g\n", gjack.en[e].jack[fit_info.corr_id[1]][j]);
    //     printf("lhs_fun: gjack.en[e].jack[fit_info.corr_id[2]][j]=%g\n", gjack.en[e].jack[fit_info.corr_id[2]][j]);
    // }
    double r = gjack.en[e].jack[fit_info.corr_id[0]][j];   // fpi
    r *= gjack.en[e].jack[fit_info.corr_id[1]][j];         // amu_l
    r /= gjack.en[e].jack[fit_info.corr_id[2]][j] / hbarc; // a
    return r;
}
double lhs_fun_fpi(int n, int e, int j, data_all gjack, struct fit_type fit_info)
{

    double r = gjack.en[e].jack[fit_info.corr_id[0]][j]; // fpi
    return r;
}
double rhs_1overamu(int n, int Nvar, double *x, int Npar, double *P)
{
    double amu = x[0];
    double r = P[0] / amu;
    return r;
}

double rhs_1overamu_mu2_mu3(int n, int Nvar, double *x, int Npar, double *P)
{
    double amu = x[0];
    double r = P[0] / amu + P[1] / (amu * amu) + P[2] / (amu * amu * amu);
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
    char label[NAMESIZE];
    file_out_name(const char *path, const char *label)
    {
        mysprintf(this->path, NAMESIZE, "%s", path);
        mysprintf(this->basename, NAMESIZE, "%s", path);
        mysprintf(this->label, NAMESIZE, "%s", label);
        mysprintf(this->namefile, NAMESIZE, "%s/%s_fit_extra.txt", path, label);
    }
    file_out_name(const file_out_name &other)
    {
        mysprintf(this->path, NAMESIZE, "%s", other.path);
        mysprintf(this->basename, NAMESIZE, "%s", other.basename);
        mysprintf(this->label, NAMESIZE, "%s", other.label);
        mysprintf(this->namefile, NAMESIZE, "%s", other.namefile);
    }
    file_out_name &operator=(const file_out_name &other)
    {
        mysprintf(this->path, NAMESIZE, "%s", other.path);
        mysprintf(this->basename, NAMESIZE, "%s", other.basename);
        mysprintf(this->label, NAMESIZE, "%s", other.label);
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

void compute_fpi_at_mciso(data_all jackall, fit_type fit_info, fit_result fit_res, int id_dfpi, file_out_name f_name)
{
    int Njack = jackall.en[0].Njack;
    error(fit_info.Njack != jackall.en[0].Njack, 1, "compute_fpi_at_mciso", "fit_info.Njack != jackall.en[0].Njack");
    error(fit_info.Njack != myres->Njack, 1, "compute_fpi_at_mciso", "fit_info.Njack != myres->Njack");
    double *fpi_mciso = (double *)malloc(sizeof(double) * Njack);
    // double *fpi_mcisochi = (double *)malloc(sizeof(double) * Njack);
    double *relative_error = (double *)malloc(sizeof(double) * Njack);
    double *tmp_x = (double *)malloc(sizeof(double) * fit_info.Nvar);
    double *tif = (double *)malloc(sizeof(double) * fit_info.Npar);

    double *fpi_sim_mev = myres->create_copy(jackall.en[0].jack[id_dfpi]);

    double *fpi_mcisochidof = (double *)malloc(sizeof(double) * Njack);
    double **Pchidof = malloc_2<double>(fit_info.Npar, Njack);
    double *relative_error_chidof = (double *)malloc(sizeof(double) * Njack);
    for (int p = 0; p < fit_info.Npar; p++)
        myres->mult_error(Pchidof[p], fit_res.P[p], std::sqrt(fit_res.chi2[Njack - 1]));
    double **tifchidof = swap_indices(fit_info.Npar, Njack, Pchidof);

    double *fpi_mcisochi = (double *)malloc(sizeof(double) * Njack);
    double **Pchi = malloc_2<double>(fit_info.Npar, Njack);
    double *relative_error_chi = (double *)malloc(sizeof(double) * Njack);
    for (int p = 0; p < fit_info.Npar; p++)
        myres->mult_error(Pchi[p], fit_res.P[p], std::sqrt(fit_res.chi2[Njack - 1] * fit_res.dof));
    double **tifchi = swap_indices(fit_info.Npar, Njack, Pchi);

    printf("Pchidof:\n");
    for (int p = 0; p < fit_info.Npar; p++)
    {
        printf("P[%d]      : %g   %g\n", p, fit_res.P[p][Njack - 1], myres->comp_error(fit_res.P[p]));
        printf("Pchidof[%d]: %g   %g\n", p, Pchidof[p][Njack - 1], myres->comp_error(Pchidof[p]));
        printf("Pchi[%d]   : %g   %g\n", p, Pchi[p][Njack - 1], myres->comp_error(Pchi[p]));
    }
    for (int n = 0; n < fit_info.N; n++)
    {
        int id = fit_info.Nxen[n][0];
        myres->copy(fpi_sim_mev, jackall.en[id].jack[id_dfpi]);                         // w0
        myres->mult(fpi_sim_mev, fpi_sim_mev, jackall.en[id].jack[fit_info.corr_id[2]]);// convert to fm
        
        for (int j = 0; j < Njack; j++)
        {
            double amuliso = jackall.en[id].jack[fit_info.corr_id[1]][j];
            double amuciso = jackall.en[id].jack[15][j];
            double amucsim = jackall.en[id].jack[19][j];
            double a_fm = jackall.en[id].jack[fit_info.corr_id[2]][j]; // a

            tmp_x[0] = amuciso / amuliso; // amuciso/amuliso
            for (int p = 0; p < fit_info.Npar; p++)
                tif[p] = fit_res.P[p][j];
            double der = fit_info.function(n, fit_info.Nvar, tmp_x, fit_info.Npar, tif) * (a_fm / hbarc) / amuliso;
            fpi_mciso[j] = (jackall.en[id].jack[id_dfpi][j] + (amuciso - amucsim) * der) * (a_fm );
            relative_error[j] = der * (amuciso - amucsim) / jackall.en[id].jack[id_dfpi][j];

            double derchidof = fit_info.function(n, fit_info.Nvar, tmp_x, fit_info.Npar, tifchidof[j]) * (a_fm / hbarc) / amuliso;
            fpi_mcisochidof[j] = (jackall.en[id].jack[id_dfpi][j] + (amuciso - amucsim) * derchidof) * (a_fm );
            relative_error_chidof[j] = derchidof * (amuciso - amucsim) / jackall.en[id].jack[id_dfpi][j];

            double derchi = fit_info.function(n, fit_info.Nvar, tmp_x, fit_info.Npar, tifchi[j]) * (a_fm / hbarc) / amuliso;
            fpi_mcisochi[j] = (jackall.en[id].jack[id_dfpi][j] + (amuciso - amucsim) * derchi) * (a_fm );
            relative_error_chi[j] = derchi * (amuciso - amucsim) / jackall.en[id].jack[id_dfpi][j];

            // relative_error[j] = fit_res.P[0][j] * (amuciso - amucsim);
            // relative_error_chidof[j] = Pchidof[0][j] * (amuciso - amucsim);
        }

        printf("fpi(mc_sim) %g  %g\n", fpi_sim_mev[Njack - 1], myres->comp_error(fpi_sim_mev));
        printf("fpi(mc_iso) %g  %g\n", fpi_mciso[Njack - 1], myres->comp_error(fpi_mciso));

        char nameextra[NAMESIZE];
        mysprintf(nameextra, NAMESIZE, "%s/%s_fit_extra_n%d.txt", f_name.path, f_name.label, n);
        printf("writing in file %s\n", nameextra);

        FILE *f = open_file(nameextra, "w+");
        fprintf(f, "%-20.12g %-20.12g\n", fpi_sim_mev[Njack - 1], myres->comp_error(fpi_sim_mev));
        fprintf(f, "%-20.12g %-20.12g\n", fpi_mciso[Njack - 1], myres->comp_error(fpi_mciso));
        fprintf(f, "%-20.12g %-20.12g\n", relative_error[Njack - 1], myres->comp_error(relative_error));
        fprintf(f, "%-20.12g %-20.12g\n", fpi_mcisochidof[Njack - 1], myres->comp_error(fpi_mcisochidof));
        fprintf(f, "%-20.12g %-20.12g\n", relative_error_chidof[Njack - 1], myres->comp_error(relative_error_chidof));
        fprintf(f, "%-20.12g %-20.12g\n", fpi_mcisochi[Njack - 1], myres->comp_error(fpi_mcisochi));
        fprintf(f, "%-20.12g %-20.12g\n", relative_error_chi[Njack - 1], myres->comp_error(relative_error_chi));

        printf("chi2/dof=%g  dof=%d\n", fit_res.chi2[Njack - 1], fit_res.dof);
        printf("%-20.12g %-20.12g\n", fpi_sim_mev[Njack - 1], myres->comp_error(fpi_sim_mev));
        printf("%-20.12g %-20.12g\n", fpi_mciso[Njack - 1], myres->comp_error(fpi_mciso));
        printf("%-20.12g %-20.12g\n", relative_error[Njack - 1], myres->comp_error(relative_error));
        printf("%-20.12g %-20.12g\n", fpi_mcisochidof[Njack - 1], myres->comp_error(fpi_mcisochidof));
        printf("%-20.12g %-20.12g\n", relative_error_chidof[Njack - 1], myres->comp_error(relative_error_chidof));
        printf("%-20.12g %-20.12g\n", fpi_mcisochi[Njack - 1], myres->comp_error(fpi_mcisochi));
        printf("%-20.12g %-20.12g\n", relative_error_chi[Njack - 1], myres->comp_error(relative_error_chi));

        fclose(f);
    }
    free(tif);
    free(tmp_x);
    free(fpi_sim_mev);
    free(fpi_mciso);
    free(relative_error);
    free(fpi_mcisochidof);
    free_2(fit_info.Npar, Pchidof);
    free_2(Njack, tifchidof);
    free(fpi_mcisochi);
    free_2(fit_info.Npar, Pchi);
    free_2(Njack, tifchi);
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
    int beta_count = 0;
    int file_count = 0;
    // std::vector<std::string> beta_names;
    std::vector<std::vector<int>> myen(1, std::vector<int>());

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
                myen[beta_count].emplace_back(file_count);
                file_count++;
            }
            else if (word.size() == 1)
            {
                if (strcmp(word[0].c_str(), "new_beta") == 0)
                {
                    myen.emplace_back(std::vector<int>());
                    beta_count++;
                }
                else
                {
                    printf("impossible to give a meaning to the line: %s\n", word[0].c_str());
                    exit(1);
                }
            }
            else if (word.size() == 0)
            {
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

    // std::vector<int> myen(files.size());
    // for (int e = 0; e < files.size(); e++)
    // {
    //     myen[e] = e;
    // }

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
    // int seed;
    // std::vector<double *> amuiso(3);
    // line_read_param(option_read, "muliso", mean, err, seed, namefile_plateaux);
    // amuiso[0] = myres->create_fake(mean, err, seed);
    // line_read_param(option_read, "musiso", mean, err, seed, namefile_plateaux);
    // amuiso[1] = myres->create_fake(mean, err, seed);
    // line_read_param(option_read, "muciso", mean, err, seed, namefile_plateaux);
    // amuiso[2] = myres->create_fake(mean, err, seed);

    // std::vector<double *> amusim(3);
    // line_read_param(option_read, "mulsim", mean, err, seed, namefile_plateaux);
    // amusim[0] = myres->create_fake(mean, err, seed);
    // line_read_param(option_read, "mussim", mean, err, seed, namefile_plateaux);
    // amusim[1] = myres->create_fake(mean, err, seed);
    // line_read_param(option_read, "mucsim", mean, err, seed, namefile_plateaux);
    // amusim[2] = myres->create_fake(mean, err, seed);
    //
    int id_dw0_dmu = 4;

    //////////////////////////////////////////////////////////////
    // fitting
    //////////////////////////////////////////////////////////////
    std::vector<std::string> der_name = {"w0"};
    std::vector<int> id_der = {id_dw0_dmu};

    int id_amuliso = 13;
    int id_a_fm = 16;
    int id_w0 = 1;
    int id_mu = 11;
    for (int i = 0; i < 1; i++)
    {
        fit_type fit_info;

        fit_info.corr_id = {id_dw0_dmu, id_amuliso, id_a_fm};

        fit_info.Nxen = std::vector<std::vector<int>>(myen.size());
        for (int n = 0; n < myen.size(); n++)
        {
            fit_info.Nxen[n].resize(myen[n].size());
            for (int e = 0; e < myen[n].size(); e++)
            {
                fit_info.Nxen[n][e] = myen[n][e];
            }
        }
        fit_info.init_N_etot_form_Nxen();
        fit_info.function = rhs_1overamu_mu2_mu3;
        fit_info.linear_fit = true;
        fit_info.Npar = 3;
        fit_info.Nvar = 1; // amuciso/amuliso
        fit_info.Njack = jackall.en[0].Njack;
        fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);

        int count = 0;
        for (int n = 0; n < fit_info.N; n++)
        {
            for (int e : fit_info.Nxen[n])
            {
                for (int j = 0; j < Njack; j++)
                {
                    // printf(" %g  %d  %d\n", jackall.en[e].jack[11][j], e, j);
                    fit_info.x[0][count][j] = jackall.en[e].jack[id_mu][j] / jackall.en[e].jack[id_amuliso][j]; // mu/muliso
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

        std::string namefit = "der_" + der_name[i] + "_full_mu_mu2_mu3";

        fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
        fit_info.band_range = {1, 350};
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "amu", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 1, {});

        // // compute fpi at the charm point
        file_out_name f_name_c_full(argv[3], namefit.c_str());
        compute_fpi_at_mciso(jackall, fit_info, der_fpi_const_full, id_w0, f_name_c_full);
    }
}
