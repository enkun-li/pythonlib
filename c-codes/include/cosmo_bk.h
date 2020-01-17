/* ===============================================
 * File Name: cosmo_bk.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 14:54:08
 * =============================================== 
 */

#ifndef __COSMO_BK__
#define __COSMO_BK__

int initial_lcdm(double omegam, double omegak, double hubble);
double lcdm_hz(double z);
double lcdm_one_over_hz(double z, double *args);
double lcdm_comoving_distance(double z);
double lcdm_growth(double z, double sig8, double gam);

double lcdm_como_ang_diam_distance(double z); // in Mpc
double lcdm_luminosity_distance(double z);
double lcdm_angulardiameter_distance(double z);

double lcdm_mu_SN(double z);
double lcdm_Dobs_SL(double zl, double zs);

#endif /* __COSMO_BK__ */

