/* ===============================================
 * File Name: gapp_rec.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-28 22:32:03
 * =============================================== 
 */

#ifndef __GAPP_REC__
#define __GAPP_REC__
double rec_mu0_over_mu(double x, double * args);
double rec_covxy_over_muxy(double x, double y);
double rec_cov_mu0_over_mu(double x, double y);
double rec_covmu0mu_over_fmu0mu(double x, double y, double * args);

double int_mu0_over_mu(double x);
double cov_int_mu0_over_mu(double x, double y);

int return_Pantheon_dz(void);

#endif /* __GAPP_REC__ */

