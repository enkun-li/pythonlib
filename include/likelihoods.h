/* ===============================================
 * File Name: likelihoods.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-10-21 17:13:08
 * =============================================== 
 */

#ifndef __LIKELIHOOD__
#define __LIKELIHOOD__

#include <stdio.h>

// ==============================
// HST likelihood
// ------------------------------
int initial_HST(char * dat_file);
double Likelihood_HST(double * param, void * model);

#endif /* __LIKELIHOOD__ */

