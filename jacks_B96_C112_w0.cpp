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
#include "correlators_analysis.hpp"
#include "fit_all.hpp"
#include "fve.hpp"
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

int main() {

  for (int Nj : {20, 25, 30, 35, 40, 45}) {
    myres = new resampling_jack(Nj);
    // double *w0B96 = myres->create_fake_exact(2.1288775, 0.0007680, 67); //BK
    double *w0B96 = myres->create_fake_exact(2.12888,  0.000655983, 67); // MG
    // double *t0B96 = myres->create_fake_exact(1.8035479, 0.0004070, -1); //BK
    double *t0B96 = myres->create_fake_exact(1.80355,  0.000353077, -1); // MG
    
    std::string name =std::string("deriv/w0_flow_B96.dat_jack") + std::to_string(Nj)  + std::string(".dat");
    myres->write_jack_in_file(w0B96,  name.c_str());
    name =std::string("deriv/sqrtt0_flow_B96.dat_jack") + std::to_string(Nj)  + std::string(".dat");
    myres->write_jack_in_file(t0B96, name.c_str());

    double *w0C112 = myres->create_fake_exact(2.5099262, 0.0016857, -1); // BK
    // double *w0C112 = myres->create_fake_exact(2.50697,  0.00181551, -1);    // MG
    double *t0C112 = myres->create_fake_exact(2.1112191, 0.0009974, -1);    // BK
    // double *t0C112 = myres->create_fake_exact(2.10995,  0.000833663, -1);   // MG
    name =std::string("deriv/w0_flow_C112.dat_jack") + std::to_string(Nj)  + std::string(".dat");
    myres->write_jack_in_file(w0C112,  name.c_str());
    name =std::string("deriv/sqrtt0_flow_C112.dat_jack") + std::to_string(Nj)  + std::string(".dat");
    myres->write_jack_in_file(t0C112, name.c_str());

    double *w0D128 = myres->create_fake_exact(3.01568,  0.00101822, -1); 
    double *t0D128 = myres->create_fake_exact(2.51797,  0.000382474, -1);   
    double *MPSD128 = myres->create_fake_exact(0.0405759070698846,    1.844897606897e-05, -1); 
    double *fPSD128 = myres->create_fake_exact(0.0377153212176888,    2.25775254422091e-05, -1);   
    name =std::string("deriv/w0_flow_D128.dat_jack") + std::to_string(Nj)  + std::string(".dat");
    myres->write_jack_in_file(w0D128,  name.c_str());
    name =std::string("deriv/sqrtt0_flow_D128.dat_jack") + std::to_string(Nj)  + std::string(".dat");
    myres->write_jack_in_file(t0D128, name.c_str());
    name =std::string("deriv/M_PS_D128.dat_jack") + std::to_string(Nj)  + std::string(".dat");
    myres->write_jack_in_file(MPSD128,  name.c_str());
    name =std::string("deriv/f_PS_D128.dat_jack") + std::to_string(Nj)  + std::string(".dat");
    myres->write_jack_in_file(fPSD128,  name.c_str());
    
    delete myres;

  }
}
