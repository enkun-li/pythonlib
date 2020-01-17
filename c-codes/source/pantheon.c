/* ===============================================
 * File Name: pantheon.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-27 16:09:42
 * =============================================== 
 */

#include "common.h"
#include "utils.h"
#include "cosmo_bk.h"

typedef struct {
    int n; // number
    double * x; // x
    double * y; // f(x)
    double * s; // sigma
    double * cov; // covariance
} DATA;

DATA * sn;

// initial Pantheon data
int initial_Pantheon(int N, char * file_dat, char * file_cov)
{
    int i, j;
    FILE * fp=NULL;

    // initial data array
    sn = (DATA *)malloc(sizeof(DATA));
    sn->n = N;
    sn->x = (double *)calloc(N, sizeof(double));
    sn->y = (double *)calloc(N, sizeof(double));
    sn->s = (double *)calloc(N, sizeof(double));
    sn->cov = (double *)calloc(N*N, sizeof(double));

    // open data file
    fp = fopen(file_dat, "r");

    // read in data
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf\n", &sn->x[i], &sn->y[i], &sn->s[i]);
    }

    // close data file
    fclose(fp);

    // open cov file
    fp = fopen(file_cov, "r");

    // read in cov
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            fscanf(fp, "%lf\n", &sn->cov[i*N+j]);
            if(i==j)
                sn->cov[i*N+j] += pow(sn->s[i], 2.0);
            //printf("%8d %16.8f\n", i*N+j, sn->cov[i*N+j]);
        }
    }

    // close cov file
    fclose(fp);

    // inverse the cov
    gsl_matrix_view M = gsl_matrix_view_array(sn->cov, sn->n, sn->n);

    gsl_linalg_cholesky_decomp(&M.matrix);
    gsl_linalg_cholesky_invert(&M.matrix);

    return _SUCCESS_;
}

// free Pantheon
int free_Pantheon(void)
{
    free(sn);
    return _SUCCESS_;
}

// sn_loglikelihood
double sn_loglikelihood(double (*func)(double))
{
    int i;
    double magA, magB, magC, chisq;
    double difmu[sn->n];
    //double logDet;

    for (i=0; i<sn->n; i++) {
        difmu[i] = (*func)(sn->x[i]) - sn->y[i];
    }

    // get inv Cov
    gsl_matrix_view M = gsl_matrix_view_array(sn->cov, sn->n, sn->n);
    gsl_vector_view V = gsl_vector_view_array(difmu, sn->n);

    magC = 0.0;
//# pragma omp parallel for
    for (i=0; i<sn->n*sn->n; i++) {
        magC += sn->cov[i];
    }

    gsl_vector *C = gsl_vector_calloc(sn->n);
    gsl_blas_dgemv(CblasNoTrans,
            1.0, &M.matrix, &V.vector,
            0.0, C);

    magB = 0.0;
//# pragma omp parallel for
    for (i=0; i<sn->n; i++) {
        magB += gsl_vector_get(C, i);
    }

    gsl_blas_ddot(&V.vector, C, &magA);

    gsl_vector_free(C);
    chisq = magA - pow(magB, 2.0)/magC + log(magC/_PI_);

    // check whether is nan
    if(gsl_isnan(chisq))
        return -3.0e30;
    return -0.5*chisq;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* output LCDM model loglikelihood with Pantheon data */
/* -------------------------------------------------- */
/* int N = 1048;
 * char * sn_dat = "/home/ekli/myworks/cosmodata/sn_full_dat.txt";
 * char * sn_cov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
 * initial_Pantheon(N, sn_dat, sn_cov);
 **/
double lcdm_pantheon(double h0, double omm, double omk)
{
   initial_lcdm(omm,omk,h0);
   return sn_loglikelihood(&lcdm_mu_SN);
}


