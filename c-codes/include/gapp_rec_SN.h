/* ===============================================
 * File Name: gapp_rec_SN.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-28 16:02:30
 * =============================================== 
 */

#ifndef __GAPP_REC_SN__
#define __GAPP_REC_SN__

int initial_SL_usedat(int num);

//int initial_SL_loglike(int N, char * file_dat, char * file_cov);
int initial_SL_loglike(int N, char * file_cov);

int free_SL(void);
//double SL_loglike(double MB, double H0, double omk, double fE);
double SL_loglike(double MB, double omk, double fE);
int initial_TDSL_loglike(int N, char * file_dat, char * file_cov);
int free_TDSL(void);
double TDSL_loglike(double MB, double H0, double omk);

//double return_SL_loglike(double MB, double omk, double fE, double H0);
double return_SL_loglike(double MB, double omk, double fE);
double return_TDSL_loglike(double MB, double omk, double H0);

//double margin_of_MB_SL_loglike(double omk, double fE, double H0);
double margin_of_MB_SL_loglike(double omk, double fE);
double margin_of_MB_TDSL_loglike(double omk, double H0);

#endif /* __GAPP_REC_SN__ */

