/* ===============================================
 * File Name: gapp.c
 * Author: ekli
 * Mail: ekli_091@mail.dlut.edu.cn  
 * Created Time: 2019-06-21 15:27:24
 * ===============================================
 */

#include "common.h"
#include "kernel.h"
#include "utils.h"

/* @@@@@@@@@@@@@@@@@@@ */
/* typedef data struct */
/* ~~~~~~~~~~~~~~~~~~~ */
typedef struct {
    int n;
    gsl_vector * x;
    gsl_vector * y;
    gsl_vector * sig;
    gsl_matrix * cov;
} DATA;

static DATA * obs ;

/* @@@@@@@@@@@@@@@@@@@@@@ */
/* typedef gp pars struct */
/* ~~~~~~~~~~~~~~~~~~~~~~ */
typedef struct {
    double sigf;
    double lenf;
    double logDet;
    gsl_matrix * IMcov;
} GAPP_pars;

static GAPP_pars *gp;

/* @@@@@@@@@@@@@@@@ */
/* initial obs data */
/* @@@@@@@@@@@@@@@@ */
int initial_gapp(int N, char * filename)
{
    int i;
    double x, y, sig;
    FILE *fp = NULL;

    // initial obs array
    obs = (DATA *)malloc(sizeof(DATA));
    obs->n = N;
    obs->x = gsl_vector_calloc(N);
    obs->y = gsl_vector_calloc(N);
    obs->sig = gsl_vector_calloc(N);
    obs->cov = gsl_matrix_calloc(N,N);

    // open data file
    fp = fopen(filename, "r");

    // read in obs data
    for (i=0; i<N; i++) {
        fscanf(fp, "%lf %lf %lf\n", &x, &y, &sig);
        gsl_vector_set(obs->x, i, x);
        gsl_vector_set(obs->y, i, y);
        gsl_vector_set(obs->sig, i, sig);
    }   

    // close data file
    fclose(fp);

    return _SUCCESS_;
}

/* @@@@@@@@@@@@@@@ */
/* initial obs cov */
/* ~~~~~~~~~~~~~~~ */
int initial_gapp_cov(int N, char *filename)
{
    int i, j;
    double covij;
    FILE *fp = NULL;

    // open data file
    fp = fopen(filename, "r");

    // read in obs data
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            fscanf(fp, "%lf\n", &covij);
            gsl_matrix_set(obs->cov, i, j, covij);
            // printf("%lf\n", obs->cov[i]);
        }
    }

    // close data file
    fclose(fp);

    return _SUCCESS_;
}

/* @@@@@@@@@@@@@@@@@@@@@@ */
/* set_up gapp covariance */
/* ~~~~~~~~~~~~~~~~~~~~~~ */
int setup_gapp(double sigf, double lenf)
{
    int i, j;
    double x, y, cij, mij, sig;

    // initial gp pars
    gp = (GAPP_pars *)malloc(sizeof(GAPP_pars));
    gp->sigf = sigf;
    gp->lenf = lenf;
    gp->IMcov = gsl_matrix_calloc(obs->n, obs->n);

    // initial M = Kernel + Cij
    for (i=0; i<obs->n; i++) {
        x = gsl_vector_get(obs->x, i);
        for (j=0; j<obs->n; j++) {
            y = gsl_vector_get(obs->x, j);
            cij = gsl_matrix_get(obs->cov, i, j);
            
            if (i==j)
                sig = gsl_vector_get(obs->sig, i);
            else
                sig = 0;

            mij = kernel(x, y, sigf, lenf) + cij + pow(sig, 2.0);
            gsl_matrix_set(gp->IMcov, i, j, mij);
            //printf("%lf", gp->IMcov->data[i*gp->IMcov->tda +j]);
        }
        //printf("\n");
    }

    //gp->logDet = gsl_matrix_inv(gp->IMcov);
    gsl_linalg_cholesky_decomp(gp->IMcov);
    gp->logDet = 0;
# pragma omp parallel for
    for (i=0; i<obs->n; i++) {
        gp->logDet += log( gp->IMcov->data[i*gp->IMcov->tda +i] );
    }
    gp->logDet = 2.0*gp->logDet;
    gsl_linalg_cholesky_invert(gp->IMcov);

    return _SUCCESS_;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* free malloc/calloc struct/array */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
int free_gapp(void)
{
    gsl_vector_free(obs->x);
    gsl_vector_free(obs->y);
    gsl_vector_free(obs->sig);
    gsl_matrix_free(obs->cov);
    free(obs);
    gsl_matrix_free(gp->IMcov);
    free(gp);

    return _SUCCESS_;
}

/* @@@@@@@@ */
/* rec meam */
/* ~~~~~~~~ */
double rec_mu(double x)
{
    int i;
    double result;
    gsl_vector *V = gsl_vector_calloc(obs->n);
    gsl_vector *C = gsl_vector_calloc(obs->n);

# pragma omp parallel for
    for (i=0; i<obs->n; i++) {
        gsl_vector_set(V, i, 
                kernel(x, obs->x->data[i], gp->sigf, gp->lenf));
    }

    gsl_blas_dgemv(CblasNoTrans, 
            1.0, gp->IMcov, obs->y,
            0.0, C);
    
    gsl_blas_ddot(V, C, &result);

    gsl_vector_free(V);
    gsl_vector_free(C);
    
    return result;
}

/* @@@@@@@@@@@@@@ */
/* rec covariance */
/* ~~~~~~~~~~~~~~ */
double rec_covariance(double x, double y)
{
    double cov;
    int i;
    gsl_vector *Vi = gsl_vector_calloc(obs->n);
    gsl_vector *Vj = gsl_vector_calloc(obs->n);
    gsl_vector *C = gsl_vector_calloc(obs->n);
    
# pragma omp parallel for
    for (i=0; i<obs->n; i++) {
        gsl_vector_set(Vi, i, 
                kernel(x, obs->x->data[i], gp->sigf, gp->lenf));
        gsl_vector_set(Vj, i, 
                kernel(obs->x->data[i], y, gp->sigf, gp->lenf));
    }
    
    gsl_blas_dgemv(CblasNoTrans, 
            1.0, gp->IMcov, Vj,
            0.0, C);
    
    gsl_blas_ddot(Vi, C, &cov);
    cov = kernel(x,y, gp->sigf, gp->lenf) - cov;
    
    gsl_vector_free(Vi);
    gsl_vector_free(Vj);
    gsl_vector_free(C);
    return cov;
}

/* @@@@@@@@@@@@@@@@@@ */
/* gapp loglikelihood */
/* ~~~~~~~~~~~~~~~~~~ */
double loglikelihood(void)
{
    double chisq;
    gsl_vector *C = gsl_vector_calloc(obs->n);
    
    gsl_blas_dgemv(CblasNoTrans, 
            1.0, gp->IMcov, obs->y,
            0.0, C);
    
    gsl_blas_ddot(obs->y, C, &chisq);

    chisq += gp->logDet + obs->n*log(2*_PI_);

    gsl_vector_free(C);
    return -0.5*chisq;
}

