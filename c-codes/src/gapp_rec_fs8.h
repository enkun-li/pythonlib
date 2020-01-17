/* ===============================================
 * File Name: gapp_rec_fs8.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 22:10:33
 * =============================================== 
 */

#ifndef __GAPP_REC_FS8__
#define __GAPP_REC_FS8__

double rec_fs8_with_Hz(double z, double sig8, double ommh2, double gam);
double rec_cov_fs8_with_Hz(double zi, double zj, double sig8, double ommh2, double gam);

double rec_fs8_with_Ez(double z, double sig8, double omm, double gam);
double rec_cov_fs8_with_Ez(double zi, double zj, double sig8, double omm, double gam);

int initial_fs8_obs(int N, char *file_fz, char *file_cov);
double fs8_loglikelihood_Hz(double sig8, double ommh2, double gam);
double fs8_loglikelihood_Ez(double sig8, double omm, double gam);

double fs8_loglikelihood_Hz_t(double sig8, double ommh2, double gam);
double fs8_loglikelihood_Hz_not(double sig8, double ommh2, double gam);
double fs8_loglikelihood_Ez_not(double sig8, double omm, double gam);

double rec_fs8_Ii_Hz(double z, double gam);
double rec_fs8_Mij_Hz(double x, double y, double gam);
double rec_fs8_Nij_Hz(double zi, double zj, double gam);

double rec_fs8_Ii_Ez(double z, double gam);
double rec_fs8_Mij_Ez(double x, double y, double gam);
double rec_fs8_Nij_Ez(double zi, double zj, double gam);

int free_fs8obs();

#endif /* __GAPP_REC_FS8__ */

