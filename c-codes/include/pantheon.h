/* ===============================================
 * File Name: pantheon.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-27 16:51:39
 * =============================================== 
 */

#ifndef __PANTHEON__
#define __PANTHEON__

int initial_Pantheon(int N, char * file_dat, char * file_cov);
int free_Pantheon(void);
double sn_loglikelihood(double (*func)(double));

double lcdm_pantheon(double h0, double omm, double omk);

#endif /* __PANTHEON__ */

