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
    mysprintf(option[6], NAMESIZE, argv[3]);         // infile
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
    error(argc != 9, 1, "main ",
        "usage:././w0_rew -p path file -bin $bin  jack/boot   file_mu  out\n separate "
        "path and file please   argc =%d", argc);

    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, argv[6]);
    printf("resampling: %s\n", resampling);

    char** option = argv_to_options(argv);


    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");




    //////////////////////////////////// read and setup header
    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], option[6]);
    FILE* infile_A0 = open_file(namefile, "r");
    generic_header head_A0;
    head_A0.read_header_debug(infile_A0);
    init_global_head(head_A0);

    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[7]);
    FILE* infile_A0_mu = open_file(namefile, "r");
    generic_header head_A0_mu;
    head_A0_mu.read_header_debug(infile_A0_mu);
    init_global_head(head_A0_mu);


    // out argv mu
    double mu0 = head_A0.mus[0];
    double mu1 = head_A0_mu.mus[0];
    head_A0.mus = { mu0, mu1 };
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], argv[8]);
    FILE* outfile = open_file(namefile, "w+");
    head_A0.write_header(outfile);


    //////////////////////////////////////////////////////////////
    // reading flow
    //////////////////////////////////////////////////////////////
    error(head_A0.Njack != head_A0_mu.Njack, 1, "main", "A0 file does not have the same confs of P5 file\n");
    int ncorr_new = head_A0.ncorr;        // current number of correlators
    int Max_corr = head_A0.ncorr; // max number of correlators

    double**** data = calloc_corr(head_A0.Njack, Max_corr, head_A0.T);

    for (int iconf = 0; iconf < head_A0.Njack; iconf++) {
        read_twopt(infile_A0, data[iconf], head_A0);

    }

    double**** data_mu = calloc_corr(head_A0_mu.Njack, 2, head_A0_mu.T);
    for (int iconf = 0; iconf < head_A0_mu.Njack; iconf++) {
        read_twopt(infile_A0_mu, data_mu[iconf], head_A0_mu);

    }
    double dmu = head_A0.mus[0] - head_A0_mu.mus[0];
    for (int iconf = 0; iconf < head_A0.Njack; iconf++) {
        for (int iv = 0; iv < head_A0.ncorr; iv++) {
            for (int t = 0; t < head_A0.T; t++) {
                data[iconf][iv][t][0] = (data[iconf][iv][t][0] - data_mu[iconf][0][t][0]) / dmu;
            }
        }

    }

    for (int iconf = 0; iconf < head_A0.Njack; iconf++) {
        fwrite(&iconf, sizeof(int), 1, outfile);
        for (int iv = 0; iv < head_A0.ncorr; iv++) {
            for (int t = 0; t < head_A0.T; t++) {
                fwrite(data[iconf][iv][t], sizeof(double), 2, outfile);
            }
        }
    }

    fclose(outfile);


}