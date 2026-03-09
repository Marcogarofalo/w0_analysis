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

double lhs_fun(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r = gjack.en[e].jack[0][j]; // d(w0/a)/d(a*mu)
    return r;
}
double lhs_fun_m05(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r = gjack.en[e].jack[1][j]; // d(w0/a)/d(a*mu)
    return r;
}
double lhs_fun_m05_mx2(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r = gjack.en[e].jack[2][j]; // d(w0/a)/d(a*mu)
    return r;
}


double full_fun(double x) {
    return 0.1 + 0.2 * x + 10 * x * x;
}


int main() {
    data_all jackall;
    jackall.resampling = "jack";

    int Ne = 3;
    jackall.en = new data_single[Ne];
    jackall.ens = Ne;

    int Njack = 50;
    myres = new resampling_jack(Njack - 1);

    std::vector<double> myx = { 0.1,0.2,0.3 };

    for (int i = 0; i < Ne; i++) {
        jackall.en[i].Nobs = 3;
        jackall.en[i].Njack = Njack;
        jackall.en[i].jack = double_malloc_2(jackall.en[i].Nobs, Njack);
    }


    // y original
    jackall.en[0].jack[0] = myres->create_fake(full_fun(myx[0]), 0.01, 1);
    jackall.en[1].jack[0] = myres->create_fake(full_fun(myx[1]), 0.01, 1);
    jackall.en[2].jack[0] = myres->create_fake(full_fun(myx[2]), 0.01, 1);


    // y -0.05*x
    jackall.en[0].jack[1] = myres->create_fake(full_fun(myx[0]) - 2.5 * myx[0], 0.01, 1);
    jackall.en[1].jack[1] = myres->create_fake(full_fun(myx[1]) - 2.5 * myx[1], 0.01, 1);
    jackall.en[2].jack[1] = myres->create_fake(full_fun(myx[2]) - 2.5 * myx[2], 0.01, 1);


    // y -0.05*x + x^2
    jackall.en[0].jack[2] = myres->create_fake(full_fun(myx[0]) - 2.5 * myx[0] - 10.0 * myx[0] * myx[0], 0.01, 1);
    jackall.en[1].jack[2] = myres->create_fake(full_fun(myx[1]) - 2.5 * myx[1] - 10.0 * myx[1] * myx[1], 0.01, 1);
    jackall.en[2].jack[2] = myres->create_fake(full_fun(myx[2]) - 2.5 * myx[2] - 10.0 * myx[2] * myx[2], 0.01, 1);

    fit_type fit_info;


    fit_info.N = 1;
    fit_info.Nxen = std::vector<std::vector<int>>(fit_info.N);
    for (int n = 0; n < fit_info.N; n++) {
        fit_info.Nxen[n].resize(Ne);
        for (int b = 0; b < Ne; b++) {
            fit_info.Nxen[n][b] = b;  // the charm is in myen[beta][0]
        }
    }
    fit_info.init_N_etot_form_Nxen();

    fit_info.function = rhs_a2;
    fit_info.linear_fit = true;
    fit_info.Npar = 2;
    fit_info.Nvar = 1; // a2
    fit_info.Njack = Njack;
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);


    for (int n = 0; n < fit_info.N; n++) {
        for (int e : fit_info.Nxen[n]) {
            for (int j = 0; j < Njack; j++) {
                fit_info.x[0][0][j] = myx[0];
                fit_info.x[0][1][j] = myx[1];
                fit_info.x[0][2][j] = myx[2];
            }
        }
    }
    fit_info.verbosity = 0;
    char** argv = (char**)malloc(sizeof(char*) * 7);
    argv[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    argv[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    argv[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    mysprintf(argv[1], NAMESIZE, "jack");        // blind/see/read_plateaux
    mysprintf(argv[3], NAMESIZE, "test_fit");        // blind/see/read_plateaux
    fit_info.band_range = { 0, 0.31 };

    {
        std::string namefit = "fit_test_ori_data";
        fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "x", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.002, {});
    }

    {
        std::string namefit = "fit_test_m05_data";
        fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun_m05, fit_info, namefit.c_str());
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "x", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.002, {});
    }
    {
        std::string namefit = "fit_test_m05_mx2_data";
        fit_result der_fpi_const_full = fit_all_data(argv, jackall, lhs_fun_m05_mx2, fit_info, namefit.c_str());
        print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "x", der_fpi_const_full, der_fpi_const_full, 0, fit_info.Nxen[0][0] /* set the other variables to the first of the n*/, 0.002, {});
    }
    
    return 0;

}