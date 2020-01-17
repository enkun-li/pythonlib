/* ===============================================
 * File Name: gapp.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 13:08:24
 * =============================================== 
 */

#ifndef __GAPP__
#define __GAPP__

int initial_gapp(int N, char * filename);
int initial_gapp_cov(int N, char * filename);
int setup_gapp(double sigf, double lenf);
double rec_mu(double x);
double rec_covariance(double x, double y);
double loglikelihood();
int free_gapp();

#endif /* __GAPP__ */

