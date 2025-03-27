#define functions_w0_C
#include "functions_w0.hpp"
#include "read.hpp"

void make_ratio_of_jacks(double ****final, int Njack, int i, int T, double ****num, int i1, double ****den,  int i2)
{

    // double sum_r = 0;
    for (int j = 0; j < Njack; j++)
    {
        for (int tf = 0; tf < T; tf++)
        {
            for (int j1 = 0; j1 < Njack-1; j1++)
            {
                if (j1 == j)
                    continue;
                double r = 0;
                for (int j2 = 0; j2 < Njack-1; j2++)
                {
                    if (j2 == j)
                        continue;
                    else
                    {
                        r += exp(den[j2][i2][0][0] - num[j1][i1][tf][0]);
                        // printf("r=%g   %g    %g \n", r, conf_jack_r[j2][0][0][0], conf_jack_rO[j1][0][tf][0]);
                    }
                }
                final[j][i][tf][0] += 1.0 / r;
            }
        }
    }
}

double ****bin_intoN_exp(double ****data, int ivar, int T, int Nconf_in, int Nb)
{

    double clustSize = ((double)Nconf_in) / ((double)Nb);
    // double clustSize = ((double)confs.confs_after_binning) / ((double)Nb);
    double ****to_write = calloc_corr(Nb, 1, T);
    for (size_t iClust = 0; iClust < Nb; iClust++)
    {

        /// Initial time of the bin
        const double binBegin = iClust * clustSize;
        /// Final time of the bin
        const double binEnd = binBegin + clustSize;
        double binPos = binBegin;
        const size_t iConf0 = floor(binPos + 1e-10);
        do
        {
            /// Index of the configuration related to the time
            const size_t iConf = floor(binPos + 1e-10);

            /// Rectangle left point
            const double beg = binPos;

            /// Rectangle right point
            const double end = std::min(binEnd, iConf + 1.0);

            /// Rectangle horizontal size
            const double weight = end - beg;

            // Perform the operation passing the info
            for (int t = 0; t < T; t++)
            {
                // printf("iConf %ld  iConf0 %ld ivar %d  t %d\n", iConf, iConf0, ivar, t);
                to_write[iClust][ivar][t][0] += weight * exp(data[iConf][ivar][t][0] - data[iConf0][ivar][t][0]);
                // to_write[iClust][ivar][t][1] += weight * exp(data[iConf][ivar][t][1]- data[0][ivar][t][1]); // imag part not implemented
            }
            // printf("Cluster=%ld  iConf=%ld  weight=%g  size=%g  end=%g  beg=%g\n",iClust,iConf,weight,clustSize, end,beg);
            // Updates the position
            binPos = end;
        } while (binEnd - binPos > 1e-10);
        for (int t = 0; t < T; t++)
        {
            to_write[iClust][ivar][t][0] /= ((double)clustSize);
            // to_write[iClust][ivar][t][1] /= ((double)clustSize);

            to_write[iClust][ivar][t][0] = log(to_write[iClust][ivar][t][0]);
            to_write[iClust][ivar][t][0] += data[iConf0][ivar][t][0];
        }
    }
    return to_write;
}

double ****bin_intoN_exp1(double ****data, double ****data_noexp, int ivar, int ivar_noexp, int T, int Nconf_in, int Nb)
{

    double clustSize = ((double)Nconf_in) / ((double)Nb);
    // double clustSize = ((double)confs.confs_after_binning) / ((double)Nb);
    double ****to_write = calloc_corr(Nb, 1, T);
    for (size_t iClust = 0; iClust < Nb; iClust++)
    {

        /// Initial time of the bin
        const double binBegin = iClust * clustSize;
        /// Final time of the bin
        const double binEnd = binBegin + clustSize;
        double binPos = binBegin;
        const size_t iConf0 = floor(binPos + 1e-10);
        do
        {
            /// Index of the configuration related to the time
            const size_t iConf = floor(binPos + 1e-10);

            /// Rectangle left point
            const double beg = binPos;

            /// Rectangle right point
            const double end = std::min(binEnd, iConf + 1.0);

            /// Rectangle horizontal size
            const double weight = end - beg;

            // Perform the operation passing the info
            for (int t = 0; t < T; t++)
            {
                to_write[iClust][ivar][t][0] += weight * data_noexp[iConf][ivar_noexp][t][0] * exp(data[iConf][ivar][0][0] - data[iConf0][ivar][0][0]);
                // to_write[iClust][ivar][t][1] += weight * exp(data[iConf][ivar][t][1]- data[0][ivar][t][1]); // imag part not implemented
            }
            // printf("Cluster=%ld  iConf=%ld  weight=%g  size=%g  end=%g  beg=%g\n",iClust,iConf,weight,clustSize, end,beg);
            // Updates the position
            binPos = end;
        } while (binEnd - binPos > 1e-10);
        for (int t = 0; t < T; t++)
        {
            to_write[iClust][ivar][t][0] /= ((double)clustSize);
            // to_write[iClust][ivar][t][1] /= ((double)clustSize);

            to_write[iClust][ivar][t][0] = log(to_write[iClust][ivar][t][0]);
            to_write[iClust][ivar][t][0] += data[iConf0][ivar][0][0];
        }
    }
    return to_write;
}

double lhs_function_w0_eg(int j, double ****in, int t, struct fit_type fit_info)
{
    int id = fit_info.corr_id[0];
    double r = in[j][id][t][0];
    // if(j==fit_info.Njack-1) printf("%d   %g\n",t,r);
    return r;
}

double lhs_function_Wt_p_dmcorr(int j, double ****in, int t, struct fit_type fit_info)
{
    int id = fit_info.corr_id[0];
    int id_cor = fit_info.corr_id[1];
    int reim = fit_info.myen[0];
    double dmu = fit_info.ext_P[0][j];

    double Wt = in[j][id][t][0];
    double loop = in[j][id_cor][0][reim];
    double loop_Wt = in[j][id_cor][t][reim];

    double r = Wt - dmu * (loop_Wt - loop * Wt);
    // if(j==fit_info.Njack-1) printf("%d   %g\n",t,r);
    return r;
}

double lhs_function_Wt_der_mu(int j, double ****in, int t, struct fit_type fit_info)
{
    int id = fit_info.corr_id[0];
    int id_cor = fit_info.corr_id[1];
    int reim = fit_info.myen[0];
    // double dmu = fit_info.ext_P[0][j];

    double Wt = in[j][id][t][0];
    double loop = in[j][id_cor][0][reim];
    double loop_Wt = in[j][id_cor][t][reim];

    double r = - (loop_Wt - loop * Wt);
    // if(j==fit_info.Njack-1) printf("%d   %g\n",t,r);
    return r;
}

double lhs_function_Wt_p_dm_all_corr(int j, double ****in, int t, struct fit_type fit_info)
{
    int id = fit_info.corr_id[0];
    int id_cor_l = fit_info.corr_id[1];
    int id_cor_s = fit_info.corr_id[2];
    int id_cor_c = fit_info.corr_id[3];
    int reim = fit_info.myen[0];
    double dmul = fit_info.ext_P[0][j];
    double dmus = fit_info.ext_P[1][j];
    double dmuc = fit_info.ext_P[2][j];

    double Wt = in[j][id][t][0];

    double loop_l = in[j][id_cor_l][0][reim];
    double loop_Wt_l = in[j][id_cor_l][t][reim];

    double loop_s = in[j][id_cor_s][0][reim];
    double loop_Wt_s = in[j][id_cor_s][t][reim];

    double loop_c = in[j][id_cor_l][0][reim];
    double loop_Wt_c = in[j][id_cor_l][t][reim];

    double r = Wt - dmul * (loop_Wt_l - loop_l * Wt) -
               dmus * (loop_Wt_s - loop_s * Wt) -
               dmuc * (loop_Wt_c - loop_c * Wt);
    // if(j==fit_info.Njack-1) printf("%d   %g\n",t,r);
    return r;
}
