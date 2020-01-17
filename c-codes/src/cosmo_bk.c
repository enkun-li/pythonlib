/* ===============================================
 * File Name: cosmo_bk.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 14:46:42
 * =============================================== 
 */

#include "common.h"
#include "utils.h"

/* @@@@@@@@@@@@@@@@@@@@ */
/*      LCDM model      */
/* @@@@@@@@@@@@@@@@@@@@ */

double omm, h0;

int initial_lcdm(double omegam, double hubble)
{
    omm = omegam;
    h0 = hubble;
    return _SUCCESS_;
}

double lcdm_hz(double z)
{
    return h0 * sqrt(omm * pow(1+z, 3.0) + 1-omm);
}

double lcdm_one_over_hz(double z, double *args)
{
    return pow(1+z, *(args))/pow(lcdm_hz(z), *(args+1));
}

double lcdm_comoving_distance(double z)
{
    double pars[] = {0.0, 1.0};

    return qgausleg(&lcdm_one_over_hz, 0.0, z, &pars[0]);
}

double lcdm_growth(double z, double sig8, double gam)
{
    double pars1[] = {3*gam, 2*gam};
    double pars2[] = {3*gam-1, 2*gam};
    double B0 = pow(omm*h0*h0, gam);
    double integ = qgausleg(&lcdm_one_over_hz, 0.0, z, &pars2[0]);
    
    return sig8 * B0 * lcdm_one_over_hz(z, &pars1[0]) * exp(-B0* integ);
}

double lcdm_mu_SN(double z)
{
    return 5*log10((1+z)*_C_/1000.0*lcdm_comoving_distance(z));
}
