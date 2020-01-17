/* ===============================================
 * File Name: gapp_sn_lens.h
 * Author: ekli
 * Mail: ekli_091@mail.dlut.edu.cn  
 * Created Time: 2019-06-22 21:26:37
 * ===============================================
 */

#ifndef __GAPP_SN_LENS__
#define __GAPP_SN_LENS__

double rec_distance_noM(double z);
double rec_cov_distance_noM(double zi, double zj);

int output_mean_cov_of_SL(int N, char * filename, char *out_mean, char *out_cov);
int output_mean_cov_of_TDSL(int N, char * filename);

double rec_distance(double z, double MB);
double rec_cov_distance(double zi, double zj, double MB);

double rec_distance_ls(double zl, double zs, double MB, double omk);
double rec_pd_dls_over_dl(double zl, double zs, double MB, double omk);
double rec_pd_dls_over_ds(double zl, double zs, double MB, double omk);
double rec_Einstein_rings(double zl, double zs, double MB, double omk);
double rec_pd_D_over_dl(double zl, double zs, double MB, double omk);
double rec_pd_D_over_ds(double zl, double zs, double MB, double omk);

int initial_SL(int N, char * filename);
int initial_TDSL(int N, char * filename);
int free_SL_TDSL(void);

double SL_loglikelihood(double MB, double omk, double fE);

#endif /* __GAPP_SN_LENS__ */

