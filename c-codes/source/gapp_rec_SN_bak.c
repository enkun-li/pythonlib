/* ===============================================
 * File Name: gapp_rec_SN.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-28 14:20:55
 * =============================================== 
 */

#include "common.h"
#include "gapp.h"
#include "gapp_sn_lens.h"
#include "utils.h"

typedef struct {
    int n;
    double *zl;
    double *zs;
    double *dl;
    double *ds;
    double *obs;
    double *sig;
    double *lam;
    double *cll;
    double *cls;
    double *csl;
    double *css;
} DATA;

DATA *sl;
DATA *tdsl;

int flag_SL = 0;
int flag_TDSL = 0;

/* @@@@@@@@@@@@@@@@@@@@ */
/*   SL_loglike   */
/* -------------------- */
int initial_SL_loglike(
        int N, 
        char * file_dat, 
        char * file_cov)
{
    int i;
    FILE *fp = NULL;
    FILE *fps = NULL;
    
    // initial zl, zs and dl, ds
    sl = (DATA *)malloc(sizeof(DATA));
    sl->n = N;
    sl->zl  = (double *)calloc(N, sizeof(double));
    sl->zs  = (double *)calloc(N, sizeof(double));
    sl->dl  = (double *)calloc(N, sizeof(double));
    sl->ds  = (double *)calloc(N, sizeof(double));
    sl->obs = (double *)calloc(N, sizeof(double));
    sl->sig = (double *)calloc(N, sizeof(double));
    sl->lam = (double *)calloc(N, sizeof(double));
    sl->cll = (double *)calloc(N*N, sizeof(double));
    sl->cls = (double *)calloc(N*N, sizeof(double));
    sl->csl = (double *)calloc(N*N, sizeof(double));
    sl->css = (double *)calloc(N*N, sizeof(double));

    // open file
    fp = fopen(file_dat, "r");
    fps = fopen(file_cov, "r");

    // read in means and obs data
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf\n",
                &sl->zl[i], &sl->zs[i], &sl->dl[i], &sl->ds[i], 
                &sl->obs[i], &sl->sig[i]);
    }

    // read in cov data
    for (i=0; i<N*N; i++) {
        fscanf(fps, "%lf %lf %lf %lf\n",
                &sl->cll[i], &sl->cls[i], &sl->csl[i], &sl->css[i]);
    }

    // close file
    fclose(fp);
    fclose(fps);

    return 0;

}

int free_SL(void)
{
    free(sl);
    flag_SL = 0;
    return _SUCCESS_;
}

double SL_loglike(
        double MB, 
        double H0, 
        double omk, 
        double fE)
{
    int N = sl->n;
    int i, j, k;
    double chisq, logDet;
    double mMB = MB - 5*log10(H0/(_C_/1000.0)) +25;
    double FMB = pow(10.0, 0.2*mMB);
    double FMBsq = pow(FMB,2.0);
    double dl, ds;//, dli, dlj, dsi, dsj, cll, cls, csl, css;
    double Dth[N], pdDdl[N], pdDds[N],
           VD[N], Cov[N*N];

    // give Dth - Dobs vector
//# pragma omp parallel for
    for (i=0; i<N; i++) {
        dl = sl->dl[i]/FMB;
        ds = sl->ds[i]/FMB;
        Dth[i] = sqrt(1+omk*pow(dl, 2.0)) - dl/ds*sqrt(1+omk*pow(ds,2.0));
        pdDdl[i] = omk * dl/sqrt(1+omk*pow(dl, 2.0)) - sqrt(1+omk*pow(ds, 2))/ds;
        pdDds[i] = dl/pow(ds, 2.0)/sqrt(1+omk*pow(ds, 2.0));
        VD[i] = Dth[i] * pow(fE, 2.0) - sl->obs[i];
    }

    // give total Cov
    for (i=0; i<N; i++) {
//# pragma omp parallel for
        for (j=0; j<=i; j++) {
            k = i*N+j;
            Cov[k] = pdDdl[i] * pdDdl[j] * sl->cll[k]/FMBsq
                + pdDdl[i] * pdDds[j] * sl->cls[k]/FMBsq
                + pdDds[i] * pdDdl[j] * sl->csl[k]/FMBsq
                + pdDds[i] * pdDds[j] * sl->css[k]/FMBsq;
            if(i==j) {
                Cov[k] += pow(sl->sig[i], 2.0);
            }
            Cov[j*N+i] = Cov[k];
        }
    }

    // get inv Cov and logDet
    gsl_matrix_view M = gsl_matrix_view_array(Cov, N, N);
    gsl_vector_view V = gsl_vector_view_array(VD, N);

    gsl_linalg_cholesky_decomp(&M.matrix);

    logDet = 0.0;
//# pragma omp parallel for
    for (i=0; i<N; i++) {
        logDet += log(Cov[i*N+i]);
    }
    logDet = 2*logDet;

    gsl_linalg_cholesky_invert(&M.matrix);

    // get dif M dif : chisq

    gsl_vector *C = gsl_vector_calloc(N);

    gsl_blas_dgemv(CblasNoTrans, 
            1.0, &M.matrix, &V.vector,
            0.0, C);
    
    gsl_blas_ddot(&V.vector, C, &chisq);

    chisq += logDet + N * log(2*_PI_);
    if(gsl_isnan(chisq)) 
        return -3.0e30;

    gsl_vector_free(C);
    return -0.5*chisq;
}

/* @@@@@@@@@@@@@@@@@@@@ */
/*   TDSL_loglike   */
/* -------------------- */
int initial_TDSL_loglike(
        int N, 
        char * file_dat, 
        char * file_cov)
{
    int i;
    FILE *fp = NULL;
    FILE *fps = NULL;
    
    // initial zl, zs and dl, ds
    tdsl = (DATA *)malloc(sizeof(DATA));
    tdsl->n = N;
    tdsl->zl  = (double *)calloc(N, sizeof(double));
    tdsl->zs  = (double *)calloc(N, sizeof(double));
    tdsl->dl  = (double *)calloc(N, sizeof(double));
    tdsl->ds  = (double *)calloc(N, sizeof(double));
    tdsl->obs = (double *)calloc(N, sizeof(double));
    tdsl->sig = (double *)calloc(N, sizeof(double));
    tdsl->lam = (double *)calloc(N, sizeof(double));
    tdsl->cll = (double *)calloc(N*N, sizeof(double));
    tdsl->cls = (double *)calloc(N*N, sizeof(double));
    tdsl->csl = (double *)calloc(N*N, sizeof(double));
    tdsl->css = (double *)calloc(N*N, sizeof(double));

    // open file
    fp = fopen(file_dat, "r");
    fps = fopen(file_cov, "r");

    // read in means and obs data
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf\n",
                &tdsl->zl[i], &tdsl->zs[i], &tdsl->dl[i], &tdsl->ds[i], 
                &tdsl->obs[i], &tdsl->sig[i], &tdsl->lam[i]);
    }

    // read in cov data
    for (i=0; i<N*N; i++) {
        fscanf(fps, "%lf %lf %lf %lf\n",
                &tdsl->cll[i], &tdsl->cls[i], &tdsl->csl[i], &tdsl->css[i]);
    }

    // close file
    fclose(fp);
    fclose(fps);

    return 0;
}

int free_TDSL(void)
{
    free(tdsl);
    flag_TDSL = 0;
    return _SUCCESS_;
}

double TDSL_loglike(
        double MB, 
        double H0, 
        double omk)
{
    int N = tdsl->n;
    int i, j, k;
    double chisq, logDet;
    double mMB = MB - 5*log10(H0/(_C_/1000.0)) +25;
    double FMB = pow(10.0, 0.2*mMB);
    double FMBsq = pow(FMB,2.0);
    double dl, ds, dls; //, dlj, dsi, dsj, cll, cls, csl, css;
    double Dth[N], pdDdl[N], pdDds[N],
           VD[N], Cov[N*N], lCov[N*N];

    // give Dth - Dobs vector
//# pragma omp parallel for
    for (i=0; i<N; i++) {
        dl = tdsl->dl[i]/FMB;
        ds = tdsl->ds[i]/FMB;
        dls =  ds*sqrt(1+omk*pow(dl, 2.0))
            - dl*sqrt(1+omk*pow(ds,2.0));
        Dth[i] = (_C_/1000.0)/H0 * dl*ds / dls;
        pdDdl[i] = (_C_/1000.0)/H0 * pow(ds, 2)/pow(dls, 2) /
            sqrt(1+omk*pow(dl, 2.0));
        pdDds[i] = (_C_/1000.0)/H0 * pow(dl, 2)/pow(dls, 2) /
            sqrt(1+omk*pow(ds, 2.0));
        VD[i] = log(Dth[i] - tdsl->lam[i])  - tdsl->obs[i];
    }

    // give total Cov
    for (i=0; i<N; i++) {
//# pragma omp parallel for
        for (j=0; j<=i; j++) {
            k = i*N+j;
            lCov[k] = pdDdl[i] * pdDdl[j] * tdsl->cll[k]/FMBsq
                + pdDdl[i] * pdDds[j] * tdsl->cls[k]/FMBsq
                + pdDds[i] * pdDdl[j] * tdsl->csl[k]/FMBsq
                + pdDds[i] * pdDds[j] * tdsl->css[k]/FMBsq;
            Cov[k] = lCov[k]/(Dth[i] - tdsl->lam[i])/(Dth[j] - tdsl->lam[j]);
            if(i==j) {
                Cov[k] += pow(tdsl->sig[i], 2.0);
                lCov[k] += pow(tdsl->sig[i], 2.0) *
                    pow((Dth[i] - tdsl->lam[i]), 2.0);
            }
            lCov[j*N+i] = lCov[k];
            Cov[j*N+i] = Cov[k];
        }
    }

    // get inv Cov and logDet
    gsl_matrix_view M = gsl_matrix_view_array(Cov, N, N);
    gsl_vector_view V = gsl_vector_view_array(VD, N);

    gsl_linalg_cholesky_decomp(&M.matrix);
    gsl_linalg_cholesky_invert(&M.matrix);

    logDet = 0.0;
    gsl_matrix_view MM = gsl_matrix_view_array(lCov, N, N);
    gsl_linalg_cholesky_decomp(&MM.matrix);
//# pragma omp parallel for
    for (i=0; i<N; i++) {
        logDet += log(lCov[i*N+i]);
    }
    logDet = 2*logDet;

    // get dif M dif : chisq

    gsl_vector *C = gsl_vector_calloc(N);

    gsl_blas_dgemv(CblasNoTrans, 
            1.0, &M.matrix, &V.vector,
            0.0, C);
    
    gsl_blas_ddot(&V.vector, C, &chisq);

    chisq += logDet + N * log(2*_PI_);
    if(gsl_isnan(chisq)) 
        return -3.0e30;

    gsl_vector_free(C);
    return -0.5*chisq;
}

/* @@@@@@@@@@@@@@@@@@@@@@ */
/*  output loglike  */
/* ---------------------- */
double return_SL_loglike(
        double MB, 
        double omk,
        double fE,
        double H0)
{
    //int N = 167;
    //char * file_dat = "SL_mean.txt";
    //char * file_cov = "SL_cov.txt";
    //int N = 205;
    //char * file_dat = "SL_mean_obs.txt";
    //char * file_cov = "SL_cov_obs.txt";
    //int N = 138; // z <= 2.26
    //char * file_dat = "SL_mean_real_snz.txt";
    //char * file_cov = "SL_cov_real_snz.txt";
    //int N = 111;
    //char * file_dat = "SL_mean_real_1.4.txt";
    //char * file_cov = "SL_cov_real_1.4.txt";
    //int N = 120;
    //char * file_dat = "SL_mean_real_1.55.txt";
    //char * file_cov = "SL_cov_real_1.55.txt";
    //int N = 126;
    //char * file_dat = "SL_mean_real_1.8.txt";
    //char * file_cov = "SL_cov_real_1.8.txt";
    //int N = 129;
    //char * file_dat = "SL_mean_real_2.0.txt";
    //char * file_cov = "SL_cov_real_2.0.txt";
    //int N = 54;
    //char * file_dat = "SL_mean_lowmass.txt";
    //char * file_cov = "SL_cov_lowmass.txt";
    //int N = 136;
    //char * file_dat = "SL_mean_intmass.txt";
    //char * file_cov = "SL_cov_intmass.txt";
    int N = 15;
    char * file_dat = "SL_mean_higmass.txt";
    char * file_cov = "SL_cov_higmass.txt";

    if(flag_SL == 0) {
        initial_SL_loglike(N,file_dat,file_cov);
        flag_SL = 1;
    }

    return SL_loglike(MB,H0,omk,fE);
}

// margin over MB
double int_SL_loglike(double MB, double * pars)
{
    return exp( 
            SL_loglike(MB,*pars, *(pars+1), *(pars+2)) 
            );
}

double margin_of_MB_SL_loglike(
        double omk,
        double fE,
        double H0)
{
    double pars[] = {H0,omk,fE};
    double chisq;
    //int N = 120;
    //char * file_dat = "./data/SL_mean_real_1.55.txt";
    //char * file_cov = "./data/SL_cov_real_1.55.txt";
    //int N = 126;
    //char * file_dat = "./data/SL_mean_real_1.8.txt";
    //char * file_cov = "./data/SL_cov_real_1.8.txt";
    int N = 138; // z <= 2.26
    char * file_dat = "./data/SL_mean_real_snz.txt";
    char * file_cov = "./data/SL_cov_real_snz.txt";
    //int N = 205;
    //char * file_dat = "./data/SL_mean_obs.txt";
    //char * file_cov = "./data/SL_cov_obs.txt";
    //int N = 54;
    //char * file_dat = "SL_mean_lowmass.txt";
    //char * file_cov = "SL_cov_lowmass.txt";
    //int N = 136;
    //char * file_dat = "./data/SL_mean_intmass.txt";
    //char * file_cov = "./data/SL_cov_intmass.txt";
    //int N = 15;
    //char * file_dat = "./data/SL_mean_higmass.txt";
    //char * file_cov = "./data/SL_cov_higmass.txt";

    if(flag_SL == 0) {
        initial_SL_loglike(N,file_dat,file_cov);
        flag_SL = 1;
    }

    chisq = log( qgausleg(&int_SL_loglike, -22.0, -17.0, &pars[0]) );
    if (gsl_isinf(chisq) == -1 || gsl_isinf(chisq) == 1)
        return -3e30;
    return chisq;
}

// @@@@@@@@@@@@@@@@@@@@@@
// TDSL loglike
double return_TDSL_loglike(
        double MB, 
        double omk,
        double H0)
{
    int N = 4;
    char * file_dat = "TDSL_mean.txt";
    char * file_cov = "TDSL_cov.txt";
    //int N = 2;
    //char * file_dat = "TDSL_mean_1.55.txt";
    //char * file_cov = "TDSL_cov_1.55.txt";

    //printf("%s\n", file_dat);

    if(flag_TDSL == 0) {
        initial_TDSL_loglike(N,file_dat,file_cov);
        flag_TDSL = 1;
    }
    return TDSL_loglike(MB,H0,omk);
}

// margin over MB
double int_TDSL_loglike(double MB, double * pars)
{
    return exp( 
            TDSL_loglike(MB,*pars, *(pars+1)) 
            );
}

double margin_of_MB_TDSL_loglike(
        double omk,
        double H0)
{
    double pars[] = {H0,omk};
    double chisq;
    int N = 4;
    char * file_dat = "./data/TDSL_mean.txt";
    char * file_cov = "./data/TDSL_cov.txt";
    //int N = 2;
    //char * file_dat = "./data/TDSL_mean_1.55.txt";
    //char * file_cov = "./data/TDSL_cov_1.55.txt";

    
    if(flag_TDSL == 0) {
        initial_TDSL_loglike(N,file_dat,file_cov);
        flag_TDSL = 1;
    }

    chisq = log( qgausleg(&int_TDSL_loglike, -22.0, -17.0, &pars[0]) );
    if (gsl_isinf(chisq) == -1 || gsl_isinf(chisq) == 1)
        return -3e30;
    return chisq;
}

