/* ===============================================
 * File Name: gapp.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 11:08:19
 * =============================================== 
 */

#include "common.h"
#include "kernel.h"
#include <omp.h>

/* typedef data struct */
struct DATA {
    int n;
    gsl_vector *x;
    gsl_vector *y;
    gsl_vector *sig;
    gsl_matrix *cov;
};

static struct DATA obs;

/* typedef gapp pars struct */
struct GAPP_pars {
    double sigf;
    double lenf;
    double logDet;
    gsl_matrix * IMcov;
};

static struct GAPP_pars gp;

/* initial obs data */
int initial_gapp(int N, char * filename)
{
    int i;
    double x, y, sig;
    FILE * fp = NULL;

    obs.n = N;
    obs.x = gsl_vector_calloc(N);
    obs.y = gsl_vector_calloc(N);
    obs.sig = gsl_vector_calloc(N);
    obs.cov = gsl_matrix_calloc(N,N);
    
    gp.IMcov = gsl_matrix_calloc(N,N);

    // open file
    fp = fopen(filename, "r");

    // read data
    for (i = 0; i< N; i++){
        fscanf(fp, "%lf %lf %lf\n", &x, &y, &sig);
        gsl_vector_set(obs.x, i, x);
        gsl_vector_set(obs.y, i, y);
        gsl_vector_set(obs.sig, i, sig);
    }

    // close file
    fclose(fp);

    return _SUCCESS_;
}

int initial_gapp_cov(char * filename)
{
    FILE *fp = NULL;
    double cij;
    int i, j;

    // open file
    fp = fopen(filename, "r");

    // read data
    for (i=0; i<obs.n; i++) {
        for (j=0; j<obs.n; j++) {
            fscanf(fp, "%lf\n", &cij);
            gsl_matrix_set(obs.cov, i, j, cij);
        }
    }

    // close file
    fclose(fp);

    return _SUCCESS_;
}

// set_up covariance
int setup_gapp(double sigf, double lenf)
{
    int i, j;
    double x,y,sig, mij;
    
    gp.sigf = sigf;
    gp.lenf = lenf;

    for (i = 0; i< obs.n; i++){
        x = gsl_vector_get(obs.x, i);
        sig = gsl_vector_get(obs.sig, i);
        for (j = 0; j< obs.n; j++){
            y = gsl_vector_get(obs.x, j);
            mij = kernel(x,y, gp.sigf, gp.lenf) + gsl_matrix_get(obs.cov, i, j);
            if(i == j)
                gsl_matrix_set(gp.IMcov, i, i, sig*sig + mij);
            else
                gsl_matrix_set(gp.IMcov, i, j, mij);
        }
    }

    gsl_linalg_cholesky_decomp(gp.IMcov);
    gp.logDet = 0;
    for (i=0; i<obs.n; i++) {
        x = gsl_matrix_get(gp.IMcov, i, i);
        gp.logDet += 2*log(x);
    }
    gsl_linalg_cholesky_invert(gp.IMcov);

    return 0;
}

// rec mu
double rec_mu(double x)
{
    double mu, xi, yj, M;
    int i,j;

    mu = 0.0;

//#   pragma omp parallel for
    for (i=0; i<obs.n; i++){
        xi = gsl_vector_get(obs.x, i);
        for (j=0; j<obs.n; j++){
            yj = gsl_vector_get(obs.y, j);
            M = gsl_matrix_get(gp.IMcov, i, j);
            mu += kernel(x, xi, gp.sigf, gp.lenf) * M * yj;
        }
    }
    return mu;
}

// rec covariance
double rec_covariance(double x, double y)
{
    double cov, xi, xj, M;
    int i,j;

    cov = 0.0;

    for (i=0; i<obs.n; i++){
        xi = gsl_vector_get(obs.x, i);
        for (j=0; j<obs.n; j++){
            xj = gsl_vector_get(obs.x, j);
            M = gsl_matrix_get(gp.IMcov, i, j);
            cov += kernel(x, xi, gp.sigf, gp.lenf) * M * kernel(xj, y, gp.sigf, gp.lenf);
        }
    }
    cov = kernel(x,y, gp.sigf, gp.lenf) - cov;
    return cov;
}

// gapp loglikelihood
double loglikelihood()
{
    int i, j;
    double chisq, yi, yj, mij;

    chisq = 0;
    for (i=0; i<obs.n; i++) {
        yi = gsl_vector_get(obs.y, i);
        for (j=0; j<obs.n; j++) {
            yj = gsl_vector_get(obs.y, j);
            mij = gsl_matrix_get(gp.IMcov, i, j);
            chisq += yi*mij*yj;
        }
    }
    chisq += gp.logDet + obs.n*log(2*_PI_);
    return -0.5*chisq;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
// free gapp
int free_gapp()
{
    gsl_vector_free(obs.x);
    gsl_vector_free(obs.y);
    gsl_vector_free(obs.sig);
    gsl_matrix_free(gp.IMcov);
    return 0;
}


