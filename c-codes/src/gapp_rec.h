/* ===============================================
 * File Name: gapp_rec.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 14:03:42
 * =============================================== 
 */

#ifndef __GAPP_REC__
#define __GAPP_REC__

double func_one_over_mu(double x, double *args);

double func_covxy_over_muxy(double x, double y);
double func_fxy_covxy_over_muxy(double x, double y, double *args);
double func_fx_covxj_over_muxj(double x, double *args);

double func_covxy_over_mux_muy(double x, double y, double *args);
double func_covxj_over_mux_muj(double x, double *args);

double integrate_fx_over_fmu(double a, double b, double *args);
double integrate_fxcovxj_over_fmux(double a, double b, double *args);
double integrate_fxcovxy_over_fmuxy(double a, double b, double c, double d, double *args);

double integrate_one_over_mu(double a, double b);
double integrate_cov_one_over_mu(double a, double b, double c, double d);

double rec_mu_over_mu0(double x);
double rec_cov_mu_over_mu0(double x, double y);
double func_mu0_over_mu(double x, double *args);
double func_mu0_mu0_covxy_over_mux_muy(double x, double y, double *args);
double func_mu0_mu0_covxj_over_mux_muj(double x, double *args);

double integrate_fxmu0_over_fmu(double a, double b, double *args);
double integrate_fxmu0covxj_over_fmux(double a, double b, double *args);
double integrate_fxmu0covxy_over_fmuxy(double a, double b, double c, double d, double *args);

double rec_mu_int_one_over_mu(double x);
double rec_cov_mu_int_one_over_mu(double x, double y);

#endif /* __GAPP_REC__ */

