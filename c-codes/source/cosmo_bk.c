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

double omm, h0, omk;
double fK;

int initial_lcdm(double omegam, double omegak, double hubble)
{
    omm = omegam;
    h0 = hubble;
    omk = omegak;
    
    if (fabs(omk) > 1.0e-8)
        fK = h0*sqrt(fabs(omk))/(_C_/1000);
    else
        fK = 1.0;

    return _SUCCESS_;
}

double lcdm_hz(double z)
{
    return h0 * sqrt(omm * pow(1+z, 3.0) + 1-omm + omk*pow(1+z,2) );
}

double lcdm_one_over_hz(double z, double *args)
{
    return pow(1+z, *(args))/pow(lcdm_hz(z), *(args+1));
}

double lcdm_comoving_distance(double z)
{
    double pars[] = {0.0, 1.0};

    return _C_/1000*qgausleg(&lcdm_one_over_hz, 0.0, z, &pars[0]);
}

double lcdm_sinn(double x)
{
    if (omk > 1.0e-8) {
        return sinh(x);
    } else if (omk < -1.0e-8) {
        return sin(x);
    } else {
        return x;
    }
}

double rofchi(double chi)
{
    return 1.0/fK*lcdm_sinn(fK*chi);
}

double lcdm_growth(double z, double sig8, double gam)
{
    double pars1[] = {3*gam, 2*gam};
    double pars2[] = {3*gam-1, 2*gam};
    double B0 = pow(omm*h0*h0, gam);
    double integ = qgausleg(&lcdm_one_over_hz, 0.0, z, &pars2[0]);
    
    return sig8 * B0 * lcdm_one_over_hz(z, &pars1[0]) * exp(-B0* integ);
}

double lcdm_como_ang_diam_distance(double z) // in Mpc
{
    return rofchi(lcdm_comoving_distance(z));
    //if (omk > 1.0e-8) {
    //    return _C_/1000.0/h0/sqrt(omk) * 
    //        sinh(h0*sqrt(omk)*lcdm_comoving_distance(z));
    //} else if (omk < -1.0e-8) {
    //    return _C_/1000.0/h0/sqrt(-omk) * 
    //        sin(h0*sqrt(-omk)*lcdm_comoving_distance(z));
    //} else {
    //    return _C_/1000.0*lcdm_comoving_distance(z);
    //}
}

double lcdm_luminosity_distance(double z)
{
    return (1+z)*rofchi(lcdm_comoving_distance(z));
}

double lcdm_angulardiameter_distance(double z)
{
    return 1/(1+z)*rofchi(lcdm_comoving_distance(z));
}

double lcdm_angulardiameter_distance2(double z1, double z2) // z2 > z1
{
    return 1/(1+z2)*rofchi(lcdm_comoving_distance(z2) - lcdm_comoving_distance(z1));
}

double lcdm_mu_SN(double z)
{
    return 5*log10(lcdm_luminosity_distance(z)) +25;
}

double lcdm_Dobs_SL(double zl, double zs)
{
    return lcdm_angulardiameter_distance2(zl, zs) /
        lcdm_angulardiameter_distance(zs);
}


