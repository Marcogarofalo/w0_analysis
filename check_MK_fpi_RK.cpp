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


int main(){

    std::vector<std::string> files = {
        "../../data_l_35/E112/fpi_sim.dat",
        "../../data_l_35/E112/MK_sim.dat",
        "../../data_l_35/E112/derivatives/sea/ml/E112_dfPi_dml.dat",
        "../../data_l_35/E112/derivatives/sea/ml/E112_dMK_dml.dat",
        "../../data_l_35/E112/derivatives/sea/ml/E112_dRK_dml.dat"
    };
    std::ifstream filep(files[0].c_str());
    if (!filep.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        return 1;
    }

    int line_count = 0;
    std::string linep;

    // Read the file line by line
    while (std::getline(filep, linep)) {
        line_count++;
    }
    filep.close();
    std::cout << "Total lines: " << line_count << std::endl;
    // jack setup
    int Njack = line_count - 1;

    // if (strcmp(argv[1], "jack") == 0) {
        myres = new resampling_jack(Njack - 1);
    // }
    // else if (strcmp(argv[1], "boot") == 0) {
    //     myres = new resampling_boot(Njack - 1);
    // }

    double *fpi = new double[Njack];
    double *MK = new double[Njack];
    double *dfpi = new double[Njack];
    double *dMK = new double[Njack];
    double *dRK = new double[Njack];
    myres->read_jack_from_file(fpi, files[0].c_str());
    myres->read_jack_from_file(MK, files[1].c_str());
    myres->read_jack_from_file(dfpi, files[2].c_str());
    myres->read_jack_from_file(dMK, files[3].c_str());
    myres->read_jack_from_file(dRK, files[4].c_str());

    printf("fpi = %.12g   (%.12g)\n", myres->mean(fpi), myres->comp_error(fpi));
    printf("MK  = %.12g   (%.12g)\n", myres->mean(MK), myres->comp_error(MK));
    printf("dfpi = %.12g   (%.12g)\n", myres->mean(dfpi), myres->comp_error(dfpi));
    printf("dMK  = %.12g   (%.12g)\n", myres->mean(dMK), myres->comp_error(dMK));
    printf("dRK = %.12g   (%.12g)\n", myres->mean(dRK), myres->comp_error(dRK));

    double *tmp = new double[Njack];
    for(int j = 0; j < Njack ;j++){
        tmp[j] = dMK[j]/fpi[j] + MK[j]*dfpi[j]/(fpi[j]*fpi[j]) ;
     }
    printf("reconstructed dRK = %.12g   (%.12g)\n", myres->mean(tmp), myres->comp_error(tmp));

}