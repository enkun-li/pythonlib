/* ===============================================
 * File Name: stronglens.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-28 08:54:55
 * =============================================== 
 */

#ifndef __STRONGLENS__
#define __STRONGLENS__

int initial_SL_data(int N, char *filename);
double sl_loglikelihood(double (*func)(double, double), double fE);
int free_stronglens(void);

double lcdm_stronglens(double h0, double omm, double omk, double fE);
double cosmography_stronglens(double omk, double a1, double a2,  double fE);

#endif /* __STRONGLENS__ */

